use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use notify::{Config as NotifyConfig, RecommendedWatcher, RecursiveMode, Watcher};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use tokio::signal;
use tracing::{debug, error, info, warn};

use crate::config::Config;
use crate::messages::sync::MessageSync;
use crate::syftbox::app::SyftBoxApp;

#[derive(Debug, Serialize, Deserialize)]
pub struct DaemonStatus {
    pid: u32,
    started_at: DateTime<Utc>,
    last_sync: Option<DateTime<Utc>>,
    message_count: usize,
    status: String,
}

impl DaemonStatus {
    fn new(pid: u32) -> Self {
        Self {
            pid,
            started_at: Utc::now(),
            last_sync: None,
            message_count: 0,
            status: "running".to_string(),
        }
    }
}

pub struct Daemon {
    config: Config,
    sync: MessageSync,
    log_writer: Arc<Mutex<std::fs::File>>,
    status: Arc<Mutex<DaemonStatus>>,
    status_path: PathBuf,
}

impl Daemon {
    pub fn new(config: &Config) -> Result<Self> {
        let biovault_dir = get_biovault_dir(config)?;
        let logs_dir = biovault_dir.join("logs");
        std::fs::create_dir_all(&logs_dir)
            .with_context(|| format!("Failed to create logs directory: {:?}", logs_dir))?;

        let log_path = logs_dir.join("daemon.log");
        let log_writer = Arc::new(Mutex::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&log_path)
                .with_context(|| format!("Failed to create log file: {:?}", log_path))?,
        ));

        let db_path = super::messages::get_message_db_path(config)?;
        let data_dir = config.get_syftbox_data_dir()?;
        let app = SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
        let sync = MessageSync::new(&db_path, app)?;

        let status_path = biovault_dir.join("daemon.status");
        let status = Arc::new(Mutex::new(DaemonStatus::new(std::process::id())));

        Ok(Self {
            config: config.clone(),
            sync,
            log_writer,
            status,
            status_path,
        })
    }

    fn log(&self, level: &str, message: &str) {
        let timestamp = Utc::now().format("%Y-%m-%d %H:%M:%S%.3f UTC");
        let log_line = format!("[{}] [{}] {}\n", timestamp, level, message);

        if let Ok(mut writer) = self.log_writer.lock() {
            let _ = writer.write_all(log_line.as_bytes());
            let _ = writer.flush();
        }

        match level {
            "ERROR" => error!("{}", message),
            "WARN" => warn!("{}", message),
            "INFO" => info!("{}", message),
            "DEBUG" => debug!("{}", message),
            _ => info!("{}", message),
        }
    }

    fn update_status(&self, last_sync: Option<DateTime<Utc>>, message_count: usize) {
        if let Ok(mut status) = self.status.lock() {
            if let Some(sync_time) = last_sync {
                status.last_sync = Some(sync_time);
            }
            status.message_count += message_count;

            let status_json = serde_json::to_string_pretty(&*status).unwrap_or_default();
            let _ = std::fs::write(&self.status_path, status_json);
        }
    }

    async fn sync_messages(&self) -> Result<()> {
        self.log("INFO", "Starting message sync");

        match self.sync.sync_quiet() {
            Ok((new_messages, count)) => {
                if count > 0 {
                    self.log("INFO", &format!("Processed {} new messages", count));
                    for msg_id in new_messages {
                        self.log("DEBUG", &format!("New message: {}", msg_id));
                    }
                } else {
                    self.log("DEBUG", "No new messages");
                }
                self.update_status(Some(Utc::now()), count);
                Ok(())
            }
            Err(e) => {
                self.log("ERROR", &format!("Message sync failed: {}", e));
                Err(e)
            }
        }
    }

    fn is_in_sbenv(&self) -> Result<bool> {
        // Check if we're in an sbenv by looking for .sbenv file
        let data_dir = self.config.get_syftbox_data_dir()?;
        let sbenv_file = data_dir.join(".sbenv");
        Ok(sbenv_file.exists())
    }

    fn check_syftbox_running(&self) -> Result<bool> {
        let data_dir = self.config.get_syftbox_data_dir()?;
        let config_path = self.config.get_syftbox_config_path()?;

        // Try to find syftbox process using the specific config file
        let output = Command::new("ps").args(["aux"]).output()?;

        let ps_output = String::from_utf8_lossy(&output.stdout);
        let config_str = config_path.to_string_lossy();

        // Check if there's a syftbox process using our config
        let is_running = ps_output.lines().any(|line| {
            line.contains("syftbox")
                && (line.contains(&*config_str)
                    || line.contains(data_dir.to_string_lossy().as_ref()))
        });

        Ok(is_running)
    }

    async fn start_syftbox(&self) -> Result<()> {
        let data_dir = self.config.get_syftbox_data_dir()?;
        let is_sbenv = self.is_in_sbenv()?;

        if is_sbenv {
            self.log("INFO", "Starting SyftBox via sbenv...");

            // Change to the data directory and run sbenv start
            let output = Command::new("sbenv")
                .arg("start")
                .arg("--skip-login-check")
                .current_dir(&data_dir)
                .output();

            match output {
                Ok(out) => {
                    if out.status.success() {
                        self.log("INFO", "Successfully started SyftBox via sbenv");
                        tokio::time::sleep(Duration::from_secs(2)).await;
                        Ok(())
                    } else {
                        let stderr = String::from_utf8_lossy(&out.stderr);
                        self.log("ERROR", &format!("sbenv start failed: {}", stderr));
                        Err(anyhow::anyhow!(
                            "Failed to start SyftBox via sbenv: {}",
                            stderr
                        ))
                    }
                }
                Err(e) => {
                    self.log("ERROR", &format!("Failed to run sbenv: {}", e));
                    Err(anyhow::anyhow!(
                        "Failed to run sbenv: {}. Is sbenv installed?",
                        e
                    ))
                }
            }
        } else {
            self.log("INFO", "Starting SyftBox directly...");

            let syftbox_config_path = self.config.get_syftbox_config_path()?;

            // Check if syftbox command is available
            let syftbox_check = Command::new("which").arg("syftbox").output();

            if syftbox_check.is_err() || !syftbox_check.as_ref().unwrap().status.success() {
                self.log("ERROR", "syftbox command not found in PATH");
                return Err(anyhow::anyhow!(
                    "SyftBox is not installed. Please install it: pip install syftbox"
                ));
            }

            // Start SyftBox with the config
            let mut start_cmd = Command::new("syftbox");
            start_cmd.arg("-c").arg(&syftbox_config_path);

            let start_output = start_cmd
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn();

            match start_output {
                Ok(mut child) => {
                    tokio::time::sleep(Duration::from_secs(3)).await;

                    match child.try_wait() {
                        Ok(None) => {
                            self.log("INFO", "Successfully started SyftBox");
                            Ok(())
                        }
                        Ok(Some(status)) => {
                            let msg = format!("SyftBox exited with status: {}", status);
                            self.log("ERROR", &msg);
                            Err(anyhow::anyhow!(msg))
                        }
                        Err(e) => {
                            let msg = format!("Failed to check SyftBox status: {}", e);
                            self.log("ERROR", &msg);
                            Err(anyhow::anyhow!(msg))
                        }
                    }
                }
                Err(e) => {
                    let msg = format!("Failed to start SyftBox: {}", e);
                    self.log("ERROR", &msg);
                    Err(anyhow::anyhow!(msg))
                }
            }
        }
    }

    async fn ensure_syftbox_running(&self) -> Result<()> {
        if self.check_syftbox_running()? {
            self.log("INFO", "SyftBox is already running");
            Ok(())
        } else {
            self.log("WARN", "SyftBox is not running, attempting to start it...");
            self.start_syftbox().await
        }
    }

    pub async fn run(&self) -> Result<()> {
        self.log("INFO", "BioVault daemon starting");

        // Log system resource info at startup
        let thread_count = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        self.log("INFO", &format!("System CPU cores: {}", thread_count));
        self.log("INFO", "Tokio worker threads: 2 (limited)");

        // Log resource limits
        let pid = std::process::id();
        if let Ok(limits) = std::fs::read_to_string(format!("/proc/{}/limits", pid)) {
            if let Some(proc_line) = limits.lines().find(|l| l.starts_with("Max processes")) {
                self.log("INFO", &format!("Resource limit: {}", proc_line.trim()));
            }
            if let Some(thread_line) = limits.lines().find(|l| l.contains("threads")) {
                self.log("INFO", &format!("Resource limit: {}", thread_line.trim()));
            }
        }

        self.log("INFO", &format!("Config email: {}", self.config.email));
        self.log(
            "INFO",
            &format!("Config syftbox_config: {:?}", self.config.syftbox_config),
        );

        // Debug: check environment variables
        if let Ok(syftbox_data_dir) = std::env::var("SYFTBOX_DATA_DIR") {
            self.log(
                "INFO",
                &format!("ENV SYFTBOX_DATA_DIR: {}", syftbox_data_dir),
            );
        } else {
            self.log("WARN", "ENV SYFTBOX_DATA_DIR not set");
        }

        let data_dir = self.config.get_syftbox_data_dir()?;
        self.log("INFO", &format!("SyftBox data_dir: {:?}", data_dir));

        let is_sbenv = self.is_in_sbenv()?;
        self.log("INFO", &format!("Is sbenv: {}", is_sbenv));

        // Ensure SyftBox is running first
        if let Err(e) = self.ensure_syftbox_running().await {
            self.log(
                "WARN",
                &format!(
                    "Could not ensure SyftBox is running: {}. Continuing anyway...",
                    e
                ),
            );
        }

        let app = SyftBoxApp::new(&data_dir, &self.config.email, "biovault")?;
        let watch_path = app
            .data_dir
            .join("datasites")
            .join(&app.email)
            .join("app_data")
            .join("biovault")
            .join("rpc")
            .join("message");

        if !watch_path.exists() {
            std::fs::create_dir_all(&watch_path)
                .with_context(|| format!("Failed to create watch directory: {:?}", watch_path))?;
        }

        self.log("INFO", &format!("Watching directory: {:?}", watch_path));

        let (tx, rx) = mpsc::channel();
        let rx = Arc::new(Mutex::new(rx));

        let mut watcher = RecommendedWatcher::new(
            move |res| {
                if let Err(e) = tx.send(res) {
                    eprintln!("Watch error: {}", e);
                }
            },
            NotifyConfig::default(),
        )?;

        watcher.watch(&watch_path, RecursiveMode::Recursive)?;

        let mut last_sync = Utc::now();
        let mut last_syftbox_check = Utc::now();
        let mut last_stats_log = Utc::now();
        let mut syftbox_restart_attempts = 0;
        const MAX_RESTART_ATTEMPTS: u32 = 3;
        const SYFTBOX_CHECK_INTERVAL_SECS: i64 = 10; // Health check every 10 seconds for testing
        const MESSAGE_SYNC_INTERVAL_SECS: i64 = 30;
        const STATS_LOG_INTERVAL_SECS: i64 = 300; // Log stats every 5 minutes

        // Use tokio interval for periodic checks (runs every 1 second)
        let mut check_interval = tokio::time::interval(Duration::from_secs(1));
        check_interval.set_missed_tick_behavior(tokio::time::MissedTickBehavior::Skip);

        self.log(
            "INFO",
            &format!(
                "Health check interval: {} seconds",
                SYFTBOX_CHECK_INTERVAL_SECS
            ),
        );
        self.log(
            "INFO",
            &format!(
                "Message sync interval: {} seconds",
                MESSAGE_SYNC_INTERVAL_SECS
            ),
        );

        loop {
            tokio::select! {
                _ = signal::ctrl_c() => {
                    self.log("INFO", "Received shutdown signal");
                    break;
                }
                _ = check_interval.tick() => {
                    let now = Utc::now();

                    // Check for file system events (checked every second for fast response)
                    let rx_clone = Arc::clone(&rx);
                    let has_event = tokio::task::spawn_blocking(move || {
                        if let Ok(rx) = rx_clone.lock() {
                            // Drain all pending events
                            let mut found_event = false;
                            while rx.try_recv().is_ok() {
                                found_event = true;
                            }
                            found_event
                        } else {
                            false
                        }
                    }).await.unwrap_or(false);

                    if has_event {
                        self.log("DEBUG", "File system event detected");
                        match self.sync_messages().await {
                            Ok(_) => {}
                            Err(e) => {
                                self.log("ERROR", &format!("Event-triggered sync failed: {}", e));
                                // Check if it's a resource error
                                if e.to_string().contains("Resource") || e.to_string().contains("thread") {
                                    self.log("FATAL", "Resource exhaustion detected - check thread/process limits");
                                }
                            }
                        }
                        last_sync = now;
                    }

                    // Periodic SyftBox health check
                    let seconds_since_check = (now - last_syftbox_check).num_seconds();
                    if seconds_since_check >= SYFTBOX_CHECK_INTERVAL_SECS {
                        self.log("INFO", &format!("Running SyftBox health check (last check {} seconds ago)...", seconds_since_check));

                        match self.ensure_syftbox_running().await {
                            Ok(_) => {
                                // Reset restart attempts on successful check
                                syftbox_restart_attempts = 0;
                                self.log("INFO", "SyftBox health check passed");
                            }
                            Err(e) => {
                                syftbox_restart_attempts += 1;
                                self.log("ERROR", &format!(
                                    "SyftBox health check failed (attempt {}/{}): {}",
                                    syftbox_restart_attempts, MAX_RESTART_ATTEMPTS, e
                                ));

                                if syftbox_restart_attempts >= MAX_RESTART_ATTEMPTS {
                                    self.log("FATAL", &format!(
                                        "Failed to restart SyftBox after {} attempts. Daemon will exit.",
                                        MAX_RESTART_ATTEMPTS
                                    ));
                                    return Err(anyhow::anyhow!(
                                        "Unable to keep SyftBox running after {} attempts",
                                        MAX_RESTART_ATTEMPTS
                                    ));
                                }
                            }
                        }
                        last_syftbox_check = now;
                    }

                    // Regular message sync
                    let seconds_since_sync = (now - last_sync).num_seconds();
                    if seconds_since_sync >= MESSAGE_SYNC_INTERVAL_SECS {
                        self.log("DEBUG", &format!("Running message sync (last sync {} seconds ago)...", seconds_since_sync));
                        if let Err(e) = self.sync_messages().await {
                            self.log("ERROR", &format!("Scheduled sync failed: {}", e));
                        }
                        last_sync = now;
                    }

                    // Periodic stats logging (every 5 minutes)
                    let seconds_since_stats = (now - last_stats_log).num_seconds();
                    if seconds_since_stats >= STATS_LOG_INTERVAL_SECS {
                        // Log process stats
                        let pid = std::process::id();

                        // Try to get thread count from /proc
                        let thread_count = std::fs::read_to_string(format!("/proc/{}/status", pid))
                            .ok()
                            .and_then(|status| {
                                status.lines()
                                    .find(|line| line.starts_with("Threads:"))
                                    .and_then(|line| line.split_whitespace().nth(1))
                                    .and_then(|s| s.parse::<usize>().ok())
                            })
                            .unwrap_or(0);

                        self.log("INFO", &format!(
                            "Stats: PID={}, Threads={}, Uptime={}min",
                            pid,
                            thread_count,
                            (now - self.status.lock().unwrap().started_at).num_minutes()
                        ));
                        last_stats_log = now;
                    }
                }
            }
        }

        self.log("INFO", "BioVault daemon shutting down");
        let _ = std::fs::remove_file(&self.status_path);
        Ok(())
    }
}

fn get_biovault_dir(config: &Config) -> Result<PathBuf> {
    // Use the syftbox data dir as the base, then add .biovault
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(data_dir.join(".biovault"))
}

fn get_pid_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    Ok(biovault_dir.join("daemon.pid"))
}

fn get_status_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    Ok(biovault_dir.join("daemon.status"))
}

fn get_log_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    let logs_dir = biovault_dir.join("logs");
    Ok(logs_dir.join("daemon.log"))
}

pub fn is_daemon_running(config: &Config) -> Result<bool> {
    let pid_path = get_pid_file_path(config)?;

    if !pid_path.exists() {
        return Ok(false);
    }

    let pid_str = std::fs::read_to_string(&pid_path)?;
    let pid: u32 = match pid_str.trim().parse() {
        Ok(p) => p,
        Err(_) => {
            // Invalid PID file, clean it up
            let _ = std::fs::remove_file(&pid_path);
            return Ok(false);
        }
    };

    // Check if process exists and is actually a bv daemon process
    let is_running = check_process_running(pid)?;

    if !is_running {
        // Clean up stale PID file
        let _ = std::fs::remove_file(&pid_path);
        let _ = std::fs::remove_file(get_status_file_path(config)?);
    }

    Ok(is_running)
}

fn check_process_running(pid: u32) -> Result<bool> {
    #[cfg(unix)]
    {
        // First check if process exists
        let exists = unsafe {
            let result = libc::kill(pid as i32, 0);
            result == 0
        };

        if !exists {
            return Ok(false);
        }

        // Verify it's actually a bv daemon process
        let output = Command::new("ps")
            .args(["-p", &pid.to_string(), "-o", "command="])
            .output()?;

        if !output.status.success() {
            return Ok(false);
        }

        let cmd = String::from_utf8_lossy(&output.stdout);
        Ok(cmd.contains("bv daemon") || cmd.contains("biovault"))
    }

    #[cfg(windows)]
    {
        let output = Command::new("tasklist")
            .args(["/FI", &format!("PID eq {}", pid), "/NH", "/FO", "CSV"])
            .output()?;

        if !output.status.success() {
            return Ok(false);
        }

        let output_str = String::from_utf8_lossy(&output.stdout);
        Ok(output_str.contains(&pid.to_string()) && output_str.contains("bv"))
    }

    #[cfg(not(any(unix, windows)))]
    {
        Ok(false)
    }
}

pub async fn start(config: &Config, foreground: bool) -> Result<()> {
    // Clean up stale PID files first
    cleanup_stale_pid_files(config)?;

    // Check if daemon is already running
    if is_daemon_running(config)? {
        let pid_path = get_pid_file_path(config)?;
        if let Ok(pid_str) = std::fs::read_to_string(&pid_path) {
            if let Ok(pid) = pid_str.trim().parse::<u32>() {
                println!("❌ Daemon is already running (PID: {})", pid);
                return Ok(());
            }
        }
        println!("❌ Daemon is already running");
        return Ok(());
    }

    let biovault_dir = get_biovault_dir(config)?;
    std::fs::create_dir_all(&biovault_dir)
        .with_context(|| format!("Failed to create biovault directory: {:?}", biovault_dir))?;

    if foreground {
        println!("🚀 Starting BioVault daemon in foreground mode");

        let pid = std::process::id();

        // Check if we're spawned from background start
        let pid_file_path = if let Ok(path_str) = std::env::var("BV_DAEMON_PID_FILE") {
            // We're a spawned daemon, use the provided PID file path
            Some(PathBuf::from(path_str))
        } else {
            // We're a manual foreground daemon, use the default PID path
            let path = get_pid_file_path(config)?;
            Some(path)
        };

        // Write PID file
        if let Some(ref pid_path) = pid_file_path {
            std::fs::write(pid_path, pid.to_string())
                .with_context(|| format!("Failed to write PID file: {:?}", pid_path))?;
        }

        let daemon = Daemon::new(config)?;

        // Wrap daemon.run() to catch and log any errors before exiting
        let result = match daemon.run().await {
            Ok(()) => {
                daemon.log("INFO", "Daemon stopped normally");
                Ok(())
            }
            Err(e) => {
                // Log the error in detail
                daemon.log("FATAL", &format!("Daemon crashed with error: {}", e));
                daemon.log("FATAL", &format!("Error chain: {:?}", e));

                // Log the backtrace
                let backtrace = e.backtrace();
                daemon.log("FATAL", &format!("Backtrace:\n{}", backtrace));

                Err(e)
            }
        };

        // Clean up PID file on exit
        if let Some(ref pid_path) = pid_file_path {
            let _ = std::fs::remove_file(pid_path);
        }

        result?;
    } else {
        println!("🚀 Starting BioVault daemon in background");

        let config_json = serde_json::to_string(config).context("Failed to serialize config")?;

        // Get the syftbox data dir now (before spawning) to ensure it's available to child
        let syftbox_data_dir = config.get_syftbox_data_dir()?;

        let current_exe =
            std::env::current_exe().context("Failed to get current executable path")?;

        // Redirect output to log file instead of /dev/null for debugging
        let log_path = get_log_file_path(config)?;

        // Ensure logs directory exists
        if let Some(log_dir) = log_path.parent() {
            std::fs::create_dir_all(log_dir)
                .with_context(|| format!("Failed to create logs directory: {:?}", log_dir))?;
        }

        let log_file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&log_path)
            .context("Failed to open log file for daemon spawn")?;

        let log_file_stderr = log_file.try_clone()?;

        // Pass the PID file path as an environment variable so the child can write it
        let pid_path = get_pid_file_path(config)?;

        let mut child = Command::new(current_exe)
            .args(["daemon", "start", "--foreground"])
            .env("BV_DAEMON_CONFIG", config_json)
            .env("BV_DAEMON_PID_FILE", pid_path.to_string_lossy().to_string())
            .env(
                "SYFTBOX_DATA_DIR",
                syftbox_data_dir.to_string_lossy().to_string(),
            )
            .stdin(Stdio::null())
            .stdout(Stdio::from(log_file))
            .stderr(Stdio::from(log_file_stderr))
            .spawn()
            .context("Failed to spawn daemon process")?;

        let pid = child.id();

        // Wait a moment for the child to start and write its PID file
        tokio::time::sleep(Duration::from_millis(1500)).await;

        // Check if child exited
        match child.try_wait() {
            Ok(Some(status)) => {
                return Err(anyhow::anyhow!(
                    "Daemon process exited with status: {}",
                    status
                ));
            }
            Ok(None) => {
                // Child is still running, verify PID file was written
                if pid_path.exists() {
                    println!("✅ Daemon started successfully (PID: {})", pid);
                    println!("📝 Use 'bv daemon logs' to view daemon logs");
                } else {
                    return Err(anyhow::anyhow!(
                        "Daemon started but PID file was not created"
                    ));
                }
            }
            Err(e) => {
                return Err(anyhow::anyhow!("Failed to check daemon status: {}", e));
            }
        }
    }

    Ok(())
}

fn cleanup_stale_pid_files(config: &Config) -> Result<()> {
    let pid_path = get_pid_file_path(config)?;
    let status_path = get_status_file_path(config)?;

    if !pid_path.exists() {
        return Ok(());
    }

    // Try to read and validate PID
    if let Ok(pid_str) = std::fs::read_to_string(&pid_path) {
        if let Ok(pid) = pid_str.trim().parse::<u32>() {
            // Check if process is actually running
            if let Ok(false) = check_process_running(pid) {
                // Process not running, clean up stale files
                let _ = std::fs::remove_file(&pid_path);
                let _ = std::fs::remove_file(&status_path);
            }
        } else {
            // Invalid PID file, clean it up
            let _ = std::fs::remove_file(&pid_path);
            let _ = std::fs::remove_file(&status_path);
        }
    }

    Ok(())
}

pub async fn logs(config: &Config, follow: bool, lines: Option<usize>) -> Result<()> {
    let log_path = get_log_file_path(config)?;

    if !log_path.exists() {
        println!("📝 No log file found. Start the daemon with 'bv daemon start' to generate logs.");
        return Ok(());
    }

    if follow {
        println!("📖 Following daemon logs (Ctrl+C to stop):");
        println!("════════════════════════════════════════");

        let mut child = Command::new("tail")
            .args(["-f", &log_path.to_string_lossy()])
            .stdout(Stdio::piped())
            .spawn()
            .context("Failed to start tail command")?;

        if let Some(stdout) = child.stdout.take() {
            let reader = BufReader::new(stdout);

            for line in reader.lines() {
                match line {
                    Ok(content) => println!("{}", content),
                    Err(_) => break,
                }
            }
        }

        let _ = child.wait();
    } else {
        let tail_lines = lines.unwrap_or(50);
        println!("📖 Last {} lines of daemon logs:", tail_lines);
        println!("═══════════════════════════════════");

        let output = Command::new("tail")
            .args(["-n", &tail_lines.to_string(), &log_path.to_string_lossy()])
            .output()
            .context("Failed to read log file")?;

        print!("{}", String::from_utf8_lossy(&output.stdout));
    }

    Ok(())
}

pub fn get_daemon_status(config: &Config) -> Result<Option<DaemonStatus>> {
    let status_path = get_status_file_path(config)?;

    if !status_path.exists() {
        return Ok(None);
    }

    let status_json = std::fs::read_to_string(&status_path)?;
    let status: DaemonStatus =
        serde_json::from_str(&status_json).context("Failed to parse daemon status")?;

    Ok(Some(status))
}

pub async fn stop(config: &Config) -> Result<()> {
    if !is_daemon_running(config)? {
        println!("❌ Daemon is not running");
        return Ok(());
    }

    let pid_path = get_pid_file_path(config)?;
    let pid_str = std::fs::read_to_string(&pid_path)?;
    let pid: u32 = pid_str.trim().parse().context("Invalid PID file")?;

    // Platform-specific process termination
    #[cfg(unix)]
    {
        unsafe {
            let result = libc::kill(pid as i32, libc::SIGTERM);
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to stop daemon"));
            }
        }
    }

    #[cfg(windows)]
    {
        let output = Command::new("taskkill")
            .args(["/F", "/PID", &pid.to_string()])
            .output()?;
        if !output.status.success() {
            return Err(anyhow::anyhow!("Failed to stop daemon"));
        }
    }

    #[cfg(not(any(unix, windows)))]
    {
        return Err(anyhow::anyhow!(
            "Stopping daemon not supported on this platform"
        ));
    }

    // Clean up PID and status files
    let _ = std::fs::remove_file(&pid_path);
    let _ = std::fs::remove_file(get_status_file_path(config)?);

    println!("✅ Daemon stopped successfully");
    Ok(())
}

fn check_systemd_available() -> Result<()> {
    // Runtime OS check instead of compile-time
    if !cfg!(target_os = "linux") {
        return Err(anyhow::anyhow!(
            "Service installation is only supported on Linux systems"
        ));
    }

    // Check if systemd is available
    let output = Command::new("systemctl")
        .arg("--version")
        .output()
        .context("systemctl not found. This system doesn't appear to use systemd")?;

    if !output.status.success() {
        return Err(anyhow::anyhow!(
            "systemd is not available on this system. Service installation requires systemd."
        ));
    }
    Ok(())
}

fn get_service_name(config: &Config) -> String {
    // Create unique service name per email/sbenv
    let safe_email = config.email.replace('@', "-at-").replace('.', "-");
    format!("biovault-daemon-{}.service", safe_email)
}

fn generate_systemd_service_content(config: &Config) -> Result<String> {
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;

    let user = std::env::var("USER").unwrap_or_else(|_| "nobody".to_string());

    // Get the data dir for this sbenv
    let data_dir = config.get_syftbox_data_dir()?;
    let data_dir_str = data_dir.to_string_lossy();

    // Check if we're in an sbenv
    let sbenv_file = data_dir.join(".sbenv");
    let is_sbenv = sbenv_file.exists();

    // Find the full path to bv (current executable)
    let bv_path = std::env::current_exe()
        .context("Failed to get current executable path")?
        .to_string_lossy()
        .to_string();

    // Generate the ExecStart command
    let exec_start = if is_sbenv {
        // Find the full path to sbenv
        let sbenv_path = Command::new("which")
            .arg("sbenv")
            .output()
            .ok()
            .and_then(|out| {
                if out.status.success() {
                    String::from_utf8(out.stdout)
                        .ok()
                        .map(|s| s.trim().to_string())
                } else {
                    None
                }
            })
            .unwrap_or_else(|| {
                // Default fallback to common installation paths
                let cargo_bin = format!("{}/.cargo/bin/sbenv", home_dir.display());
                let local_bin = format!("{}/.local/bin/sbenv", home_dir.display());
                if std::path::Path::new(&cargo_bin).exists() {
                    cargo_bin
                } else {
                    local_bin
                }
            });

        // Use sbenv exec to run bv within the sbenv context
        format!(
            "{} exec {} {} daemon start --foreground",
            sbenv_path, config.email, bv_path
        )
    } else {
        // Use the bv path directly
        format!("{} daemon start --foreground", bv_path)
    };

    // Set working directory to the data dir for sbenv context
    let working_dir = if is_sbenv {
        data_dir_str.to_string()
    } else {
        home_dir.display().to_string()
    };

    let service_content = format!(
        r#"[Unit]
Description=BioVault Daemon ({email})
After=network.target
Wants=network-online.target

[Service]
Type=simple
User={user}
Group={user}
WorkingDirectory={working_dir}
ExecStart={exec_start}
Restart=on-failure
RestartSec=10
StandardOutput=journal
StandardError=journal
SyslogIdentifier=biovault-{safe_email}
Environment="HOME={home_dir}"
Environment="PATH=/usr/local/bin:/usr/bin:/bin:{home_dir}/.local/bin:{home_dir}/.cargo/bin"

# Security settings
PrivateTmp=true
NoNewPrivileges=true
ProtectSystem=full

[Install]
WantedBy=multi-user.target
"#,
        user = user,
        home_dir = home_dir.display(),
        email = config.email,
        safe_email = config.email.replace('@', "-at-").replace('.', "-"),
        exec_start = exec_start,
        working_dir = working_dir,
    );

    Ok(service_content)
}

pub async fn install_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    // Check if service is already installed
    let service_name = get_service_name(config);
    let check_output = Command::new("systemctl")
        .args(["status", &service_name])
        .output()?;

    if check_output.status.success() || check_output.status.code() == Some(3) {
        // Status code 3 means service exists but is inactive
        println!("⚠️  Service '{}' is already installed", service_name);
        println!("   Use 'bv daemon uninstall' first if you want to reinstall");
        return Ok(());
    }

    println!("📦 Installing BioVault daemon as systemd service...");

    // Generate service file content
    let service_content = generate_systemd_service_content(config)?;

    // Write service file to temporary location
    let temp_service_path = format!("/tmp/{}", service_name);
    std::fs::write(&temp_service_path, service_content)
        .context("Failed to write temporary service file")?;

    // Move service file to systemd directory (requires sudo)
    println!("🔐 Installing service (requires sudo)...");
    let install_output = Command::new("sudo")
        .args([
            "mv",
            &temp_service_path,
            &format!("/etc/systemd/system/{}", service_name),
        ])
        .output()
        .context("Failed to install service file. Make sure you have sudo privileges")?;

    if !install_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to install service: {}",
            String::from_utf8_lossy(&install_output.stderr)
        ));
    }

    // Reload systemd daemon
    println!("🔄 Reloading systemd daemon...");
    let reload_output = Command::new("sudo")
        .args(["systemctl", "daemon-reload"])
        .output()
        .context("Failed to reload systemd daemon")?;

    if !reload_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to reload systemd: {}",
            String::from_utf8_lossy(&reload_output.stderr)
        ));
    }

    // Enable service to start on boot
    println!("🚀 Enabling service to start on boot...");
    let enable_output = Command::new("sudo")
        .args(["systemctl", "enable", &service_name])
        .output()
        .context("Failed to enable service")?;

    if !enable_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to enable service: {}",
            String::from_utf8_lossy(&enable_output.stderr)
        ));
    }

    // Start the service
    println!("▶️  Starting service...");
    let start_output = Command::new("sudo")
        .args(["systemctl", "start", &service_name])
        .output()
        .context("Failed to start service")?;

    if !start_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to start service: {}",
            String::from_utf8_lossy(&start_output.stderr)
        ));
    }

    println!("✅ BioVault daemon installed and started successfully!");
    println!("\n📊 Service Management Commands:");
    println!("   • Status:  sudo systemctl status {}", service_name);
    println!("   • Stop:    sudo systemctl stop {}", service_name);
    println!("   • Start:   sudo systemctl start {}", service_name);
    println!("   • Restart: sudo systemctl restart {}", service_name);
    println!("   • Logs:    sudo journalctl -u {} -f", service_name);
    println!("\n   Or use 'bv daemon status' for a quick check");

    Ok(())
}

pub async fn uninstall_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    let service_name = get_service_name(config);

    println!("🗑️  Uninstalling BioVault daemon service...");

    // Stop the service if running
    println!("⏹️  Stopping service...");
    let stop_output = Command::new("sudo")
        .args(["systemctl", "stop", &service_name])
        .output()?;

    if !stop_output.status.success() {
        // Service might not be running, continue anyway
        println!("   (Service was not running)");
    }

    // Disable the service
    println!("🚫 Disabling service...");
    let disable_output = Command::new("sudo")
        .args(["systemctl", "disable", &service_name])
        .output()?;

    if !disable_output.status.success() {
        println!("   (Service was not enabled)");
    }

    // Remove service file
    println!("🗑️  Removing service file...");
    let remove_output = Command::new("sudo")
        .args(["rm", "-f", &format!("/etc/systemd/system/{}", service_name)])
        .output()?;

    if !remove_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to remove service file: {}",
            String::from_utf8_lossy(&remove_output.stderr)
        ));
    }

    // Reload systemd daemon
    println!("🔄 Reloading systemd daemon...");
    let reload_output = Command::new("sudo")
        .args(["systemctl", "daemon-reload"])
        .output()?;

    if !reload_output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to reload systemd: {}",
            String::from_utf8_lossy(&reload_output.stderr)
        ));
    }

    println!("✅ BioVault daemon service uninstalled successfully!");
    Ok(())
}

pub async fn list_services() -> Result<()> {
    // Check if systemd is available (runtime check)
    if !cfg!(target_os = "linux") {
        println!("⚠️  Service listing is only supported on Linux systems");
        return Ok(());
    }

    // Check if systemctl is available
    let check = Command::new("systemctl").arg("--version").output();

    if check.is_err() || !check.as_ref().unwrap().status.success() {
        println!("⚠️  systemd is not available on this system");
        return Ok(());
    }

    println!("🔍 Searching for BioVault daemon services...\n");

    // List all biovault-daemon services
    let output = Command::new("systemctl")
        .args([
            "list-units",
            "--all",
            "--no-pager",
            "biovault-daemon-*.service",
        ])
        .output()?;

    if !output.status.success() {
        println!("⚠️  Failed to list services");
        return Ok(());
    }

    let output_str = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = output_str.lines().collect();

    // Parse the systemctl output to extract service information
    let mut services = Vec::new();
    for line in lines.iter().skip(1) {
        // Skip header
        if line.contains("biovault-daemon-") && line.contains(".service") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if let Some(name) = parts.first() {
                if name.starts_with("biovault-daemon-") {
                    // Extract email from service name
                    let email_part = name
                        .trim_start_matches("biovault-daemon-")
                        .trim_end_matches(".service")
                        .replace("-at-", "@")
                        .replace("-", ".");

                    // Get status
                    let status = if line.contains("active") && line.contains("running") {
                        "RUNNING"
                    } else if line.contains("failed") {
                        "FAILED"
                    } else if line.contains("inactive") || line.contains("dead") {
                        "STOPPED"
                    } else {
                        "UNKNOWN"
                    };

                    services.push((name.to_string(), email_part, status));
                }
            }
        }
    }

    if services.is_empty() {
        println!("📝 No BioVault daemon services found");
        println!("   Use 'bv daemon install' to install a service");
    } else {
        println!("📋 Found {} BioVault daemon service(s):\n", services.len());
        for (service_name, email, status) in services {
            let status_icon = match status {
                "RUNNING" => "✅",
                "FAILED" => "❌",
                "STOPPED" => "⏹️",
                _ => "❓",
            };
            println!("   {} {} ({})", status_icon, email, status);
            println!("      Service: {}", service_name);
            println!();
        }
        println!("💡 Commands:");
        println!("   • Status:  sudo systemctl status <service-name>");
        println!("   • Stop:    sudo systemctl stop <service-name>");
        println!("   • Start:   sudo systemctl start <service-name>");
        println!("   • Logs:    sudo journalctl -u <service-name> -f");
    }

    Ok(())
}

pub async fn service_status(config: &Config) -> Result<()> {
    // First check if daemon is running the old way (manual start)
    let manual_running = is_daemon_running(config).unwrap_or(false);

    if manual_running {
        println!("🤖 Daemon Status: RUNNING (manual mode)");
        if let Some(status) = get_daemon_status(config)? {
            println!("   • PID: {}", status.pid);
            println!(
                "   • Started: {}",
                status.started_at.format("%Y-%m-%d %H:%M:%S UTC")
            );
            if let Some(last_sync) = status.last_sync {
                println!(
                    "   • Last sync: {}",
                    last_sync.format("%Y-%m-%d %H:%M:%S UTC")
                );
            }
            println!("   • Messages processed: {}", status.message_count);
        }
        println!("\n   ℹ️  Note: Daemon was started manually with 'bv daemon start'");
        return Ok(());
    }

    // Check systemd service status (runtime check)
    if cfg!(target_os = "linux") {
        let service_name = get_service_name(config);
        let output = Command::new("systemctl")
            .args(["status", &service_name, "--no-pager"])
            .output()?;

        if output.status.success() || output.status.code() == Some(3) {
            // Parse the output to get status
            let output_str = String::from_utf8_lossy(&output.stdout);

            if output_str.contains("Active: active (running)") {
                println!("🤖 Daemon Status: RUNNING (systemd service)");
            } else if output_str.contains("Active: inactive") || output_str.contains("Active: dead")
            {
                println!("⚠️  Daemon Status: STOPPED");
            } else if output_str.contains("Active: failed") {
                println!("❌ Daemon Status: FAILED");
            } else {
                println!("❓ Daemon Status: UNKNOWN");
            }

            // Print the full systemctl status output
            println!("\n{}", output_str);

            if output_str.contains("could not be found") {
                println!(
                    "\n   ℹ️  Service is not installed. Use 'bv daemon install' to install it."
                );
            }
        } else if output.status.code() == Some(4) {
            println!("❌ Service '{}' not found", service_name);
            println!("   Use 'bv daemon install' to install the service");
        } else {
            println!("⚠️  Could not determine service status");
            println!("   Error: {}", String::from_utf8_lossy(&output.stderr));
        }
    } else {
        println!("⚠️  Daemon Status: NOT RUNNING");
        println!("   Service installation is only supported on Linux");
        println!("   Use 'bv daemon start' to run the daemon manually");
    }

    Ok(())
}

pub async fn show_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    let service_name = get_service_name(config);
    let service_path = format!("/etc/systemd/system/{}", service_name);

    if !std::path::Path::new(&service_path).exists() {
        println!("❌ Service file not found: {}", service_path);
        println!("   Use 'bv daemon install' to install the service");
        return Ok(());
    }

    println!("📄 Service file: {}", service_path);
    println!("════════════════════════════════════════");

    let content = std::fs::read_to_string(&service_path)
        .with_context(|| format!("Failed to read service file: {}", service_path))?;

    println!("{}", content);

    Ok(())
}

pub async fn reinstall_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    println!("🔄 Reinstalling BioVault daemon service...\n");

    // Check if service exists
    let service_name = get_service_name(config);
    let check_output = Command::new("systemctl")
        .args(["status", &service_name])
        .output()?;

    if check_output.status.success() || check_output.status.code() == Some(3) {
        // Service exists, uninstall it first
        println!("📤 Uninstalling existing service...");
        uninstall_service(config).await?;
        println!();
    }

    // Install the service
    println!("📥 Installing service...");
    install_service(config).await?;

    Ok(())
}
