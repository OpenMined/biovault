use crate::cli::commands::messages::init_message_system_quiet;
use crate::config::Config;
use crate::syftbox::app::SyftBoxApp;
use anyhow::anyhow;
use notify::{Config as NotifyConfig, RecommendedWatcher, RecursiveMode, Watcher};
use std::sync::{mpsc, Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

pub struct MessageRpcWatcherHandle {
    stop_tx: mpsc::Sender<()>,
    handle: Option<JoinHandle<()>>,
}

impl MessageRpcWatcherHandle {
    pub fn stop(&mut self) {
        let _ = self.stop_tx.send(());
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

impl Drop for MessageRpcWatcherHandle {
    fn drop(&mut self) {
        self.stop();
    }
}

/// Start watching the SyftBox RPC message endpoint for incoming/outgoing message activity.
/// Calls `on_sync` with any new message IDs after each filesystem event + an initial sync.
pub fn start_message_rpc_watcher<F>(
    config: Config,
    on_sync: F,
) -> crate::Result<MessageRpcWatcherHandle>
where
    F: FnMut(&[String]) + Send + 'static,
{
    let cfg = config.clone();
    let callback = Arc::new(Mutex::new(on_sync));
    let (stop_tx, stop_rx) = mpsc::channel();

    let handle = thread::spawn(move || {
        if let Err(err) = run_watcher(cfg, callback, stop_rx) {
            eprintln!("Message watcher error: {}", err);
        }
    });

    Ok(MessageRpcWatcherHandle {
        stop_tx,
        handle: Some(handle),
    })
}

fn run_watcher<F>(
    config: Config,
    callback: Arc<Mutex<F>>,
    stop_rx: mpsc::Receiver<()>,
) -> crate::Result<()>
where
    F: FnMut(&[String]) + Send + 'static,
{
    // Build SyftBox endpoint path for message RPC
    let data_dir = config.get_syftbox_data_dir()?;
    let app = SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
    let watch_path = app.register_endpoint("/message")?;

    // Initialize message system within watcher thread (DB + MessageSync)
    let (_db, sync) = init_message_system_quiet(&config)?;

    let (tx, rx) = mpsc::channel();
    let mut watcher = RecommendedWatcher::new(
        move |res| {
            let _ = tx.send(res);
        },
        NotifyConfig::default(),
    )
    .map_err(|e| anyhow!(e))?;

    watcher
        .watch(&watch_path, RecursiveMode::Recursive)
        .map_err(|e| anyhow!(e))?;

    // Initial sync to establish baseline and deliver any queued items
    if let Ok((ids, _)) = sync.sync_quiet() {
        if let Ok(mut cb) = callback.lock() {
            cb(&ids);
        }
    }

    // React to filesystem events by syncing and notifying listener
    loop {
        if stop_rx.try_recv().is_ok() {
            break;
        }

        match rx.recv_timeout(Duration::from_millis(250)) {
            Ok(Ok(_event)) => {
                if let Ok((ids, _)) = sync.sync_quiet() {
                    if let Ok(mut cb) = callback.lock() {
                        cb(&ids);
                    }
                }
            }
            Ok(Err(err)) => {
                eprintln!("Message watcher notify error: {}", err);
            }
            Err(mpsc::RecvTimeoutError::Timeout) => {}
            Err(mpsc::RecvTimeoutError::Disconnected) => break,
        }
    }

    Ok(())
}
