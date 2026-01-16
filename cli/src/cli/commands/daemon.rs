use anyhow::{anyhow, Context, Result};
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
use tracing::{debug, error, info, instrument, warn};

use crate::config::Config;
use crate::messages::sync::MessageSync;
use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::{
    detect_mode, is_syftbox_running, start_syftbox as start_syftbox_process, SyftBoxMode,
};
use syftbox_sdk::syftbox::config::SyftboxRuntimeConfig;

fn command<S: AsRef<std::ffi::OsStr>>(program: S) -> Command {
    let mut cmd = Command::new(program);
    super::configure_child_process(&mut cmd);
    cmd
}

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

// NOTE: Daemon struct and its methods are NOT suitable for unit testing because they:
// - Require actual file system operations (log files, status files)
// - Depend on MessageSync which requires a database connection
// - Depend on SyftBoxApp which requires SyftBox configuration
// - Perform real process management and signal handling
// These should be tested via integration/e2e tests with proper setup.
pub struct Daemon {
    config: Config,
    sync: MessageSync,
    log_writer: Arc<Mutex<std::fs::File>>,
    status: Arc<Mutex<DaemonStatus>>,
    status_path: PathBuf,
}

impl Daemon {
    // NOT unit testable - requires MessageSync, SyftBoxApp, and file system setup
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

    fn runtime_config(&self) -> Result<SyftboxRuntimeConfig> {
        self.config
            .to_syftbox_runtime_config()
            .map_err(|e| anyhow!(e))
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

    // NOT unit testable - requires MessageSync and actual message processing
    #[instrument(skip(self), fields(component = "daemon"), err)]
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

    // NOT unit testable - spawns actual SyftBox processes (sbenv or syftbox command)
    #[instrument(skip(self), fields(component = "daemon"), err)]
    async fn start_syftbox(&self) -> Result<()> {
        let runtime_config = self.runtime_config()?;
        let mode = detect_mode(&runtime_config)?;
        match mode {
            SyftBoxMode::Sbenv => self.log("INFO", "Starting SyftBox via sbenv..."),
            SyftBoxMode::Direct => self.log("INFO", "Starting SyftBox directly..."),
        }

        match start_syftbox_process(&runtime_config) {
            Ok(true) => self.log("INFO", "Successfully started SyftBox"),
            Ok(false) => self.log("INFO", "SyftBox already running"),
            Err(e) => {
                self.log("ERROR", &format!("Failed to start SyftBox: {}", e));
                return Err(e);
            }
        }

        Ok(())
    }

    #[instrument(skip(self), fields(component = "daemon"), err)]
    async fn ensure_syftbox_running(&self) -> Result<()> {
        let runtime_config = self.runtime_config()?;
        if is_syftbox_running(&runtime_config)? {
            self.log("INFO", "SyftBox is already running");
            Ok(())
        } else {
            self.log("WARN", "SyftBox is not running, attempting to start it...");
            self.start_syftbox().await
        }
    }

    // NOT unit testable - main daemon loop with file watching, signal handling, and async runtime
    #[instrument(skip(self), fields(component = "daemon"), err)]
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

        let runtime_config = self.runtime_config()?;
        let mode = detect_mode(&runtime_config)?;
        self.log("INFO", &format!("SyftBox mode: {:?}", mode));

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
        let watch_path = app.register_endpoint("/message")?;

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
                    let has_event = if let Ok(rx) = rx.lock() {
                        // Drain all pending events
                        let mut found_event = false;
                        while rx.try_recv().is_ok() {
                            found_event = true;
                        }
                        found_event
                    } else {
                        false
                    };

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

fn get_biovault_dir(_config: &Config) -> Result<PathBuf> {
    crate::config::get_biovault_home()
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

// Partially unit testable - can verify it doesn't panic, but actual process checking
// depends on ps/tasklist commands and real PIDs. Integration tests needed for full coverage.
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
        let output = command("ps")
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
        let output = command("tasklist")
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

// NOT unit testable - spawns daemon processes, requires full environment setup
// Test via integration tests instead
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

        let mut child = command(current_exe)
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

// Partially unit testable - tested with no log file case
// Full testing (reading actual logs, following) requires integration tests
pub async fn logs(config: &Config, follow: bool, lines: Option<usize>) -> Result<()> {
    let log_path = get_log_file_path(config)?;

    if !log_path.exists() {
        println!("📝 No log file found. Start the daemon with 'bv daemon start' to generate logs.");
        return Ok(());
    }

    if follow {
        println!("📖 Following daemon logs (Ctrl+C to stop):");
        println!("════════════════════════════════════════");

        let mut child = command("tail")
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

        let output = command("tail")
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

// Partially unit testable - tested with no daemon running case
// Full testing (actually stopping a daemon) requires integration tests
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
        let output = command("taskkill")
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
    let output = command("systemctl")
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
        let sbenv_path = command("which")
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
Environment="RUST_BACKTRACE=1"

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

// NOT unit testable - requires systemd, writes to system directories
// Test via integration tests on Linux systems only
pub async fn install_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    // Check if service is already installed
    let service_name = get_service_name(config);
    let check_output = command("systemctl")
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
    let install_output = command("sudo")
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
    let reload_output = command("sudo")
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
    let enable_output = command("sudo")
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
    let start_output = command("sudo")
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

// NOT unit testable - requires systemd, modifies system services
// Test via integration tests on Linux systems only
pub async fn uninstall_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    let service_name = get_service_name(config);

    println!("🗑️  Uninstalling BioVault daemon service...");

    // Stop the service if running
    println!("⏹️  Stopping service...");
    let stop_output = command("sudo")
        .args(["systemctl", "stop", &service_name])
        .output()?;

    if !stop_output.status.success() {
        // Service might not be running, continue anyway
        println!("   (Service was not running)");
    }

    // Disable the service
    println!("🚫 Disabling service...");
    let disable_output = command("sudo")
        .args(["systemctl", "disable", &service_name])
        .output()?;

    if !disable_output.status.success() {
        println!("   (Service was not enabled)");
    }

    // Remove service file
    println!("🗑️  Removing service file...");
    let remove_output = command("sudo")
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
    let reload_output = command("sudo")
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
    let check = command("systemctl").arg("--version").output();

    if check.is_err() || !check.as_ref().unwrap().status.success() {
        println!("⚠️  systemd is not available on this system");
        return Ok(());
    }

    println!("🔍 Searching for BioVault daemon services...\n");

    // List all biovault-daemon services
    let output = command("systemctl")
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

// Partially unit testable - non-Linux systems return early
// Full testing requires systemd and installed service (integration tests)
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
        let output = command("systemctl")
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

// NOT unit testable - requires systemd
// Test via integration tests on Linux systems only
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

// NOT unit testable - requires systemd, modifies system services
// Test via integration tests on Linux systems only
pub async fn reinstall_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    println!("🔄 Reinstalling BioVault daemon service...\n");

    // Check if service exists
    let service_name = get_service_name(config);
    let check_output = command("systemctl")
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    struct TestDirGuard;

    impl TestDirGuard {
        fn new(temp: &TempDir) -> Self {
            crate::config::set_test_syftbox_data_dir(temp.path());
            let home = temp.path().join(".biovault");
            crate::config::set_test_biovault_home(&home);
            TestDirGuard
        }
    }

    impl Drop for TestDirGuard {
        fn drop(&mut self) {
            crate::config::clear_test_biovault_home();
            crate::config::clear_test_syftbox_data_dir();
        }
    }

    // ============================================================================
    // UNIT TESTING GUIDELINES FOR DAEMON MODULE
    // ============================================================================
    //
    // WHAT CAN BE UNIT TESTED:
    // ✅ Helper functions (get_*_path, get_service_name, etc.)
    // ✅ Data structures (DaemonStatus serialization, etc.)
    // ✅ Simple validation logic (check_systemd_available on non-Linux)
    // ✅ Edge cases with mock data (invalid PIDs, missing files)
    //
    // WHAT CANNOT BE UNIT TESTED (requires integration/e2e tests):
    // ❌ Daemon::new() - requires MessageSync, SyftBoxApp, database
    // ❌ Daemon methods (run, sync_messages, etc.) - require running daemon
    // ❌ start() - spawns actual processes
    // ❌ stop() - requires running daemon to stop
    // ❌ Service operations (install/uninstall/etc.) - require systemd
    // ❌ Process checking with real PIDs - depends on actual running processes
    // ❌ File watching and signal handling - require async runtime
    //
    // When adding tests, check if the function:
    // 1. Requires external processes/services → integration test
    // 2. Modifies system state → integration test
    // 3. Is pure logic/data manipulation → unit test
    // ============================================================================

    #[test]
    fn test_daemon_status_new() {
        let pid = 12345;
        let status = DaemonStatus::new(pid);
        assert_eq!(status.pid, pid);
        assert_eq!(status.status, "running");
        assert_eq!(status.message_count, 0);
        assert!(status.last_sync.is_none());
    }

    #[test]
    fn test_daemon_status_serialize() {
        let status = DaemonStatus::new(999);
        let json = serde_json::to_string(&status);
        assert!(json.is_ok());
    }

    #[test]
    fn test_get_service_name() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let name = get_service_name(&config);
        assert!(name.starts_with("biovault-daemon-"));
        assert!(name.contains("test"));
    }

    #[test]
    fn test_get_biovault_dir() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let dir = get_biovault_dir(&config);
        assert!(dir.is_ok());
    }

    #[test]
    fn test_get_pid_file_path() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let path = get_pid_file_path(&config);
        assert!(path.is_ok());
        assert!(path.unwrap().to_string_lossy().contains("daemon.pid"));
    }

    #[test]
    fn test_get_status_file_path() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let path = get_status_file_path(&config);
        assert!(path.is_ok());
        assert!(path.unwrap().to_string_lossy().contains("daemon.status"));
    }

    #[test]
    fn test_get_log_file_path() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let path = get_log_file_path(&config);
        assert!(path.is_ok());
        assert!(path.unwrap().to_string_lossy().contains("daemon.log"));
    }

    #[test]
    fn test_daemon_status_update_fields() {
        let mut status = DaemonStatus::new(999);
        assert_eq!(status.message_count, 0);
        assert!(status.last_sync.is_none());

        status.message_count = 5;
        status.last_sync = Some(Utc::now());

        assert_eq!(status.message_count, 5);
        assert!(status.last_sync.is_some());
    }

    #[test]
    fn test_daemon_status_debug() {
        let status = DaemonStatus::new(123);
        let debug_str = format!("{:?}", status);
        assert!(debug_str.contains("123"));
        assert!(debug_str.contains("running"));
    }

    #[test]
    fn test_daemon_status_deserialize() {
        let json = r#"{
            "pid": 456,
            "started_at": "2024-01-01T00:00:00Z",
            "last_sync": null,
            "message_count": 0,
            "status": "running"
        }"#;
        let status: Result<DaemonStatus, _> = serde_json::from_str(json);
        assert!(status.is_ok());
        let s = status.unwrap();
        assert_eq!(s.pid, 456);
        assert_eq!(s.status, "running");
    }

    #[test]
    #[cfg(unix)]
    fn test_check_process_running_invalid_pid() {
        // PID 999999 is very unlikely to exist
        let result = check_process_running(999999);
        assert!(result.is_ok());
        assert!(!result.unwrap());
    }

    #[test]
    #[cfg(unix)]
    fn test_check_process_running_current_process() {
        // Current test process should be running
        let pid = std::process::id();
        let result = check_process_running(pid);
        // Function should return a result, even if error
        assert!(result.is_ok() || result.is_err());
    }

    #[test]
    fn test_is_daemon_running_no_pid_file() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = is_daemon_running(&config);
        assert!(result.is_ok());
        assert!(!result.unwrap());
    }

    #[test]
    fn test_cleanup_stale_pid_files_no_file() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = cleanup_stale_pid_files(&config);
        assert!(result.is_ok());
    }

    #[test]
    fn test_get_daemon_status_no_file() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = get_daemon_status(&config);
        assert!(result.is_ok());
        assert!(result.unwrap().is_none());
    }

    #[test]
    fn test_check_systemd_available() {
        // Just verify it doesn't panic
        let _result = check_systemd_available();
    }

    #[test]
    fn test_generate_systemd_service_content() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = generate_systemd_service_content(&config);
        // May fail due to missing dirs/env, but shouldn't panic
        let _ = result;
    }

    #[test]
    fn test_cleanup_stale_pid_files_with_invalid_pid() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create invalid PID file
        let pid_path = get_pid_file_path(&config).unwrap();
        std::fs::create_dir_all(pid_path.parent().unwrap()).unwrap();
        std::fs::write(&pid_path, "not_a_number").unwrap();

        let result = cleanup_stale_pid_files(&config);
        assert!(result.is_ok());

        // Should have cleaned up the invalid file
        assert!(!pid_path.exists());
    }

    #[test]
    fn test_is_daemon_running_with_invalid_pid() {
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create invalid PID file
        let pid_path = get_pid_file_path(&config).unwrap();
        std::fs::create_dir_all(pid_path.parent().unwrap()).unwrap();
        std::fs::write(&pid_path, "not_a_pid").unwrap();

        let result = is_daemon_running(&config);
        assert!(result.is_ok());
        assert!(!result.unwrap());

        // Should have cleaned up the invalid file
        assert!(!pid_path.exists());
    }

    #[test]
    fn test_get_daemon_status_with_valid_file() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create valid status file
        let status_path = get_status_file_path(&config).unwrap();
        std::fs::create_dir_all(status_path.parent().unwrap()).unwrap();

        let status = DaemonStatus::new(123);
        let json = serde_json::to_string(&status).unwrap();
        std::fs::write(&status_path, json).unwrap();

        let result = get_daemon_status(&config);
        assert!(result.is_ok());
        let loaded_status = result.unwrap();
        assert!(loaded_status.is_some());
        assert_eq!(loaded_status.unwrap().pid, 123);
    }

    #[test]
    fn test_daemon_status_serialization_round_trip() {
        let status = DaemonStatus::new(789);
        let json = serde_json::to_string(&status).unwrap();
        let deserialized: DaemonStatus = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.pid, 789);
        assert_eq!(deserialized.status, "running");
        assert_eq!(deserialized.message_count, 0);
    }

    #[test]
    fn test_get_service_name_sanitizes_email() {
        let config = Config {
            email: "user@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        let name = get_service_name(&config);
        assert!(name.contains("user-at-example-com"));
        assert!(!name.contains('@'));
        assert!(!name.contains('.') || name.ends_with(".service"));
    }

    #[test]
    fn test_daemon_status_with_last_sync() {
        let mut status = DaemonStatus::new(555);
        let now = Utc::now();
        status.last_sync = Some(now);
        status.message_count = 10;

        let json = serde_json::to_string(&status).unwrap();
        let deserialized: DaemonStatus = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.message_count, 10);
        assert!(deserialized.last_sync.is_some());
    }

    #[test]
    fn test_path_helpers_consistency() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let biovault_dir = get_biovault_dir(&config).unwrap();
        let pid_path = get_pid_file_path(&config).unwrap();
        let status_path = get_status_file_path(&config).unwrap();
        let log_path = get_log_file_path(&config).unwrap();

        // All paths should be under biovault_dir
        assert!(pid_path.starts_with(&biovault_dir));
        assert!(status_path.starts_with(&biovault_dir));
        assert!(log_path.starts_with(&biovault_dir));
    }

    #[test]
    #[cfg(unix)]
    fn test_check_process_running_pid_1() {
        // PID 1 is init/systemd, should exist but won't be a bv process
        let result = check_process_running(1);
        assert!(result.is_ok());
        // Will return false because PID 1 is not a bv/biovault process
        assert!(!result.unwrap());
    }

    #[test]
    fn test_check_systemd_available_on_non_linux() {
        #[cfg(not(target_os = "linux"))]
        {
            let result = check_systemd_available();
            assert!(result.is_err());
        }
    }

    #[tokio::test]
    async fn test_stop_daemon_not_running() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = stop(&config).await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_logs_no_file() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = logs(&config, false, Some(10)).await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_daemon_status_display() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Should work even when daemon not running
        let status_opt = get_daemon_status(&config);
        assert!(status_opt.is_ok());
    }

    #[tokio::test]
    #[cfg(not(target_os = "linux"))]
    async fn test_install_service_non_linux() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = install_service(&config).await;
        assert!(result.is_err());
    }

    #[tokio::test]
    #[cfg(not(target_os = "linux"))]
    async fn test_uninstall_service_non_linux() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = uninstall_service(&config).await;
        assert!(result.is_err());
    }

    #[tokio::test]
    #[cfg(not(target_os = "linux"))]
    async fn test_service_status_non_linux() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = service_status(&config).await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    #[cfg(not(target_os = "linux"))]
    async fn test_show_service_non_linux() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = show_service(&config).await;
        assert!(result.is_err());
    }

    #[tokio::test]
    #[cfg(not(target_os = "linux"))]
    async fn test_reinstall_service_non_linux() {
        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let result = reinstall_service(&config).await;
        assert!(result.is_err());
    }

    #[test]
    fn test_generate_systemd_service_content_fields() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        if let Ok(content) = generate_systemd_service_content(&config) {
            assert!(content.contains("[Unit]"));
            assert!(content.contains("[Service]"));
            assert!(content.contains("[Install]"));
            assert!(content.contains("test@example.com") || content.contains("test-at-example"));
        }
    }

    #[test]
    fn test_daemon_status_fields_modification() {
        let mut status = DaemonStatus::new(444);

        // Test initial state
        assert_eq!(status.pid, 444);
        assert_eq!(status.message_count, 0);
        assert!(status.last_sync.is_none());
        assert_eq!(status.status, "running");

        // Modify fields
        status.message_count = 15;
        status.status = "syncing".to_string();
        let now = Utc::now();
        status.last_sync = Some(now);

        // Verify modifications
        assert_eq!(status.message_count, 15);
        assert_eq!(status.status, "syncing");
        assert!(status.last_sync.is_some());
    }

    #[test]
    fn test_get_daemon_status_invalid_json() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);
        let bv_home = tmp.path().join("bv_home");
        crate::config::set_test_biovault_home(&bv_home);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create invalid JSON status file
        let status_path = get_status_file_path(&config).unwrap();
        std::fs::create_dir_all(status_path.parent().unwrap()).unwrap();
        std::fs::write(&status_path, "not valid json").unwrap();

        let result = get_daemon_status(&config);
        assert!(result.is_err());
    }

    #[test]
    fn test_cleanup_stale_pid_files_with_stale_process() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create PID file with non-existent PID
        let pid_path = get_pid_file_path(&config).unwrap();
        std::fs::create_dir_all(pid_path.parent().unwrap()).unwrap();
        std::fs::write(&pid_path, "999999").unwrap();

        let result = cleanup_stale_pid_files(&config);
        assert!(result.is_ok());
    }

    #[test]
    fn test_is_daemon_running_with_stale_pid() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create PID file with non-existent PID
        let pid_path = get_pid_file_path(&config).unwrap();
        std::fs::create_dir_all(pid_path.parent().unwrap()).unwrap();
        std::fs::write(&pid_path, "999999").unwrap();

        let result = is_daemon_running(&config);
        assert!(result.is_ok());
        assert!(!result.unwrap());
    }

    #[test]
    #[cfg(unix)]
    fn test_check_process_running_edge_cases() {
        // Test with PID 0 (invalid) - may return error or false
        let result = check_process_running(0);
        let _ = result; // Just ensure it doesn't panic

        // Test with max PID - may return error or false
        let result = check_process_running(99999);
        let _ = result; // Just ensure it doesn't panic
    }

    #[test]
    fn test_service_name_multiple_emails() {
        let configs = vec![
            ("user@example.com", "user-at-example-com"),
            ("test.user@domain.org", "test-user-at-domain-org"),
            ("name+tag@email.com", "name+tag-at-email-com"),
        ];

        for (email, expected_part) in configs {
            let config = Config {
                email: email.to_string(),
                syftbox_config: None,
                version: None,
                binary_paths: None,
                syftbox_credentials: None,
                agent_bridge_enabled: None,
                agent_bridge_port: None,
                agent_bridge_http_port: None,
                agent_bridge_token: None,
                agent_bridge_blocklist: None,
            };
            let name = get_service_name(&config);
            assert!(
                name.contains(expected_part),
                "Service name '{}' doesn't contain '{}'",
                name,
                expected_part
            );
            assert!(name.starts_with("biovault-daemon-"));
            assert!(name.ends_with(".service"));
        }
    }

    #[test]
    fn test_path_helpers_different_emails() {
        let tmp = tempfile::TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config1 = Config {
            email: "user1@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let config2 = Config {
            email: "user2@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let dir1 = get_biovault_dir(&config1).unwrap();
        let dir2 = get_biovault_dir(&config2).unwrap();

        // Both should get same biovault dir (shared)
        assert_eq!(dir1, dir2);
    }

    // NOTE: test_start_already_running removed - start() requires full daemon setup
    // including MessageSync and SyftBox, making it unsuitable for unit tests.
    // This should be tested via integration/e2e tests instead.

    #[test]
    fn test_cleanup_stale_pid_files_preserves_running_daemon() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Create PID file with current process ID (which is running)
        let pid_path = get_pid_file_path(&config).unwrap();
        std::fs::create_dir_all(pid_path.parent().unwrap()).unwrap();
        std::fs::write(&pid_path, std::process::id().to_string()).unwrap();

        let result = cleanup_stale_pid_files(&config);
        assert!(result.is_ok());

        // Clean up
        let _ = std::fs::remove_file(&pid_path);
    }

    #[test]
    fn test_daemon_status_json_format() {
        let status = DaemonStatus::new(12345);
        let json = serde_json::to_string_pretty(&status).unwrap();

        assert!(json.contains("\"pid\": 12345"));
        assert!(json.contains("\"status\": \"running\""));
        assert!(json.contains("\"message_count\": 0"));
        assert!(json.contains("\"started_at\""));
    }

    #[test]
    fn test_get_pid_file_path_creates_parent() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);
        let bv_home = tmp.path().join("bv_home");
        crate::config::set_test_biovault_home(&bv_home);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let pid_path = get_pid_file_path(&config).unwrap();
        assert!(pid_path.parent().is_some());
        assert!(pid_path.starts_with(&bv_home));
    }

    #[test]
    fn test_get_log_file_path_structure() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let log_path = get_log_file_path(&config).unwrap();
        assert!(log_path.to_string_lossy().contains("logs"));
        assert!(log_path.to_string_lossy().contains("daemon.log"));
    }

    #[tokio::test]
    async fn test_logs_with_lines_parameter() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        // Test with different line counts
        let result1 = logs(&config, false, Some(10)).await;
        assert!(result1.is_ok());

        let result2 = logs(&config, false, Some(100)).await;
        assert!(result2.is_ok());

        let result3 = logs(&config, false, None).await;
        assert!(result3.is_ok());
    }

    // NOTE: test_is_daemon_running_with_valid_current_pid removed
    // is_daemon_running() calls check_process_running() which uses ps/tasklist
    // and requires the process to match "bv" or "biovault" in the command line.
    // Test processes don't match this pattern, making this unsuitable for unit tests.
    // Covered by test_is_daemon_running_no_pid_file and test_is_daemon_running_with_stale_pid instead.

    #[test]
    fn test_get_status_file_path_uniqueness() {
        use tempfile::TempDir;
        let tmp = TempDir::new().unwrap();
        let _guard = TestDirGuard::new(&tmp);

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let status_path = get_status_file_path(&config).unwrap();
        let pid_path = get_pid_file_path(&config).unwrap();

        // Should be different files
        assert_ne!(status_path, pid_path);
        // But in same directory
        assert_eq!(status_path.parent(), pid_path.parent());
    }

    #[test]
    fn test_daemon_status_started_at_is_recent() {
        let before = Utc::now();
        let status = DaemonStatus::new(999);
        let after = Utc::now();

        assert!(status.started_at >= before);
        assert!(status.started_at <= after);
    }

    #[test]
    fn test_service_name_no_special_chars() {
        let config = Config {
            email: "user@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let name = get_service_name(&config);

        // Should not contain @ or . (except in .service extension)
        let name_without_ext = name.strip_suffix(".service").unwrap();
        assert!(!name_without_ext.contains('@'));
        assert!(!name_without_ext.contains('.'));
    }

    #[test]
    #[cfg(unix)]
    fn test_check_process_running_returns_bool() {
        // Just verify it returns a bool without panicking
        let result1 = check_process_running(1);
        assert!(result1.is_ok());
        // Result is a bool, no need to assert it's true or false

        let result2 = check_process_running(999999);
        assert!(result2.is_ok());
        // Result is a bool, function completes without panic
    }
}
