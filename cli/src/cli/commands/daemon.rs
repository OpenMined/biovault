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
use tokio::time::sleep;
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

    async fn ensure_syftbox_running(&self) -> Result<()> {
        // Check if SyftBox is running
        let output = Command::new("pgrep").arg("-f").arg("syftbox").output();

        let is_running = match output {
            Ok(out) => out.status.success() && !out.stdout.is_empty(),
            Err(_) => {
                // pgrep might not be available, try alternative check
                Command::new("ps")
                    .args(["aux"])
                    .output()
                    .map(|o| String::from_utf8_lossy(&o.stdout).contains("syftbox"))
                    .unwrap_or(false)
            }
        };

        if !is_running {
            self.log("INFO", "SyftBox is not running, attempting to start it...");

            // Try to start SyftBox
            let syftbox_config_path = self.config.get_syftbox_config_path()?;

            // First, check if syftbox command is available
            let syftbox_check = Command::new("which").arg("syftbox").output();

            if syftbox_check.is_err() || !syftbox_check.unwrap().status.success() {
                self.log(
                    "WARN",
                    "syftbox command not found in PATH. Please install SyftBox first.",
                );
                return Err(anyhow::anyhow!(
                    "SyftBox is not installed. Please install it first: pip install syftbox"
                ));
            }

            // Start SyftBox with the config from BioVault's configuration
            let mut start_cmd = Command::new("syftbox");

            // If we have a specific config path, use it
            if syftbox_config_path.exists() {
                start_cmd.arg("--config").arg(&syftbox_config_path);
            }

            let start_output = start_cmd
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn();

            match start_output {
                Ok(mut child) => {
                    // Give it a moment to start
                    tokio::time::sleep(Duration::from_secs(3)).await;

                    // Check if it's still running
                    match child.try_wait() {
                        Ok(None) => {
                            self.log("INFO", "Successfully started SyftBox");
                            Ok(())
                        }
                        Ok(Some(status)) => {
                            self.log("ERROR", &format!("SyftBox exited with status: {}", status));
                            Err(anyhow::anyhow!("SyftBox failed to start: {}", status))
                        }
                        Err(e) => {
                            self.log("ERROR", &format!("Failed to check SyftBox status: {}", e));
                            Err(anyhow::anyhow!("Failed to check SyftBox status: {}", e))
                        }
                    }
                }
                Err(e) => {
                    self.log("ERROR", &format!("Failed to start SyftBox: {}", e));
                    Err(anyhow::anyhow!("Failed to start SyftBox: {}", e))
                }
            }
        } else {
            self.log("INFO", "SyftBox is already running");
            Ok(())
        }
    }

    pub async fn run(&self) -> Result<()> {
        self.log("INFO", "BioVault daemon starting");

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

        let data_dir = self.config.get_syftbox_data_dir()?;
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

        let sync_interval = Duration::from_secs(30);
        let mut last_sync = Utc::now();
        let mut last_syftbox_check = Utc::now();

        loop {
            tokio::select! {
                _ = signal::ctrl_c() => {
                    self.log("INFO", "Received shutdown signal");
                    break;
                }
                _ = sleep(sync_interval) => {
                    let now = Utc::now();

                    // Periodic SyftBox health check
                    if (now - last_syftbox_check).num_seconds() >= 300 {
                        if let Err(e) = self.ensure_syftbox_running().await {
                            self.log("WARN", &format!("SyftBox health check failed: {}", e));
                        }
                        last_syftbox_check = now;
                    }

                    // Regular message sync
                    if (now - last_sync).num_seconds() >= 30 {
                        if let Err(e) = self.sync_messages().await {
                            self.log("ERROR", &format!("Scheduled sync failed: {}", e));
                        }
                        last_sync = now;
                    }
                }
                _ = tokio::time::sleep(Duration::from_millis(100)) => {
                    let rx_clone = Arc::clone(&rx);
                    let has_event = tokio::task::spawn_blocking(move || {
                        if let Ok(rx) = rx_clone.lock() {
                            rx.try_recv().is_ok()
                        } else {
                            false
                        }
                    }).await.unwrap_or(false);

                    if has_event {
                        self.log("DEBUG", "File system event detected");
                        if let Err(e) = self.sync_messages().await {
                            self.log("ERROR", &format!("Event-triggered sync failed: {}", e));
                        }
                        last_sync = Utc::now();
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
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
    Ok(home_dir.join(".biovault"))
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
    let status_path = get_status_file_path(config)?;

    if !pid_path.exists() || !status_path.exists() {
        return Ok(false);
    }

    let pid_str = std::fs::read_to_string(&pid_path)?;
    let pid: u32 = pid_str.trim().parse().context("Invalid PID file")?;

    // Runtime OS check with compile-time guards for platform-specific code
    #[cfg(unix)]
    {
        unsafe {
            let result = libc::kill(pid as i32, 0);
            Ok(result == 0)
        }
    }

    #[cfg(windows)]
    {
        let output = Command::new("tasklist")
            .args(["/FI", &format!("PID eq {}", pid)])
            .output()?;
        Ok(String::from_utf8_lossy(&output.stdout).contains(&pid.to_string()))
    }

    #[cfg(not(any(unix, windows)))]
    {
        Ok(false)
    }
}

pub async fn start(config: &Config, foreground: bool) -> Result<()> {
    if is_daemon_running(config)? {
        println!("❌ Daemon is already running");
        return Ok(());
    }

    let biovault_dir = get_biovault_dir(config)?;
    std::fs::create_dir_all(&biovault_dir)
        .with_context(|| format!("Failed to create biovault directory: {:?}", biovault_dir))?;

    if foreground {
        println!("🚀 Starting BioVault daemon in foreground mode");
        let daemon = Daemon::new(config)?;
        daemon.run().await?;
    } else {
        println!("🚀 Starting BioVault daemon in background");

        let config_json = serde_json::to_string(config).context("Failed to serialize config")?;

        let current_exe =
            std::env::current_exe().context("Failed to get current executable path")?;

        let mut child = Command::new(current_exe)
            .args(["start", "--foreground"])
            .env("BV_DAEMON_CONFIG", config_json)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .context("Failed to spawn daemon process")?;

        let pid = child.id();
        let pid_path = get_pid_file_path(config)?;
        std::fs::write(&pid_path, pid.to_string())
            .with_context(|| format!("Failed to write PID file: {:?}", pid_path))?;

        tokio::time::sleep(Duration::from_millis(500)).await;

        match child.try_wait() {
            Ok(Some(status)) => {
                let _ = std::fs::remove_file(&pid_path);
                return Err(anyhow::anyhow!(
                    "Daemon process exited with status: {}",
                    status
                ));
            }
            Ok(None) => {
                println!("✅ Daemon started successfully (PID: {})", pid);
                println!("📝 Use 'bv logs' to view daemon logs");
            }
            Err(e) => {
                let _ = std::fs::remove_file(&pid_path);
                return Err(anyhow::anyhow!("Failed to check daemon status: {}", e));
            }
        }
    }

    Ok(())
}

pub async fn logs(config: &Config, follow: bool, lines: Option<usize>) -> Result<()> {
    let log_path = get_log_file_path(config)?;

    if !log_path.exists() {
        println!("📝 No log file found. Start the daemon with 'bv start' to generate logs.");
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

fn get_service_name() -> String {
    "biovault-daemon.service".to_string()
}

fn generate_systemd_service_content(config: &Config) -> Result<String> {
    let exe_path = std::env::current_exe().context("Failed to get current executable path")?;

    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;

    let user = std::env::var("USER").unwrap_or_else(|_| "nobody".to_string());

    let service_content = format!(
        r#"[Unit]
Description=BioVault Daemon - Automatic message processing for BioVault
After=network.target
Wants=network-online.target

[Service]
Type=simple
User={user}
Group={user}
WorkingDirectory={home_dir}
ExecStart={exe_path} start --foreground
Restart=on-failure
RestartSec=10
StandardOutput=journal
StandardError=journal
SyslogIdentifier=biovault-daemon
Environment="HOME={home_dir}"
Environment="BV_DAEMON_EMAIL={email}"
Environment="PATH=/usr/local/bin:/usr/bin:/bin:{home_dir}/.local/bin"

# Security settings
PrivateTmp=true
NoNewPrivileges=true
ProtectSystem=strict
ProtectHome=read-only
ReadWritePaths={home_dir}/.biovault {home_dir}/.syftbox {home_dir}/dev

[Install]
WantedBy=multi-user.target
"#,
        user = user,
        home_dir = home_dir.display(),
        exe_path = exe_path.display(),
        email = config.email,
    );

    Ok(service_content)
}

pub async fn install_service(config: &Config) -> Result<()> {
    check_systemd_available()?;

    // Check if service is already installed
    let service_name = get_service_name();
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

pub async fn uninstall_service(_config: &Config) -> Result<()> {
    check_systemd_available()?;

    let service_name = get_service_name();

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
        println!("\n   ℹ️  Note: Daemon was started manually with 'bv start'");
        return Ok(());
    }

    // Check systemd service status (runtime check)
    if cfg!(target_os = "linux") {
        let service_name = get_service_name();
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
        println!("   Use 'bv start' to run the daemon manually");
    }

    Ok(())
}
