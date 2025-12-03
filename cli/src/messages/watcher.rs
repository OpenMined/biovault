use crate::cli::commands::messages::init_message_system_quiet;
use crate::config::Config;
use crate::syftbox::app::SyftBoxApp;
use anyhow::anyhow;
use notify::{Config as NotifyConfig, RecommendedWatcher, RecursiveMode, Watcher};
use std::sync::{mpsc, Arc, Mutex};
use std::thread::{self, JoinHandle};

/// Start watching the SyftBox RPC message endpoint for incoming/outgoing message activity.
/// Calls `on_sync` with any new message IDs after each filesystem event + an initial sync.
pub fn start_message_rpc_watcher<F>(config: Config, on_sync: F) -> crate::Result<JoinHandle<()>>
where
    F: FnMut(&[String]) + Send + 'static,
{
    let cfg = config.clone();
    let callback = Arc::new(Mutex::new(on_sync));

    let handle = thread::spawn(move || {
        if let Err(err) = run_watcher(cfg, callback) {
            eprintln!("Message watcher error: {}", err);
        }
    });

    Ok(handle)
}

fn run_watcher<F>(config: Config, callback: Arc<Mutex<F>>) -> crate::Result<()>
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
    for res in rx {
        match res {
            Ok(_event) => {
                if let Ok((ids, _)) = sync.sync_quiet() {
                    if let Ok(mut cb) = callback.lock() {
                        cb(&ids);
                    }
                }
            }
            Err(err) => {
                eprintln!("Message watcher notify error: {}", err);
            }
        }
    }

    Ok(())
}
