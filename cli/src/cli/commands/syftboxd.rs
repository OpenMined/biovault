use anyhow::{anyhow, Context, Result};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::Duration;
use tokio::signal;

use crate::config::Config;
use syftbox_sdk::syftbox::config::SyftboxRuntimeConfig;
use syftbox_sdk::syftbox::control::{is_syftbox_running, start_syftbox, stop_syftbox};

#[derive(Debug, Serialize, Deserialize)]
pub struct SyftboxdStatus {
    pid: u32,
    started_at: DateTime<Utc>,
    status: String,
}

impl SyftboxdStatus {
    fn new(pid: u32) -> Self {
        Self {
            pid,
            started_at: Utc::now(),
            status: "running".to_string(),
        }
    }
}

fn get_biovault_dir(_config: &Config) -> Result<PathBuf> {
    Ok(crate::config::get_biovault_home()?)
}

fn get_pid_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    Ok(biovault_dir.join("syftboxd.pid"))
}

fn get_status_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    Ok(biovault_dir.join("syftboxd.status"))
}

fn get_log_file_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = get_biovault_dir(config)?;
    let logs_dir = biovault_dir.join("logs");
    Ok(logs_dir.join("syftboxd.log"))
}

fn write_status(config: &Config, status: &SyftboxdStatus) -> Result<()> {
    let status_path = get_status_file_path(config)?;
    let json =
        serde_json::to_string_pretty(status).context("Failed to serialize syftboxd status")?;
    let _ = std::fs::write(status_path, json);
    Ok(())
}

pub fn is_syftboxd_running(config: &Config) -> Result<bool> {
    let pid_path = get_pid_file_path(config)?;
    if !pid_path.exists() {
        return Ok(false);
    }

    let pid_str = std::fs::read_to_string(&pid_path)?;
    let pid: u32 = match pid_str.trim().parse() {
        Ok(p) => p,
        Err(_) => {
            let _ = std::fs::remove_file(&pid_path);
            let _ = std::fs::remove_file(get_status_file_path(config)?);
            return Ok(false);
        }
    };

    let is_running = check_process_running(pid)?;
    if !is_running {
        let _ = std::fs::remove_file(&pid_path);
        let _ = std::fs::remove_file(get_status_file_path(config)?);
    }
    Ok(is_running)
}

fn check_process_running(pid: u32) -> Result<bool> {
    #[cfg(unix)]
    {
        let exists = unsafe { libc::kill(pid as i32, 0) == 0 };
        if !exists {
            return Ok(false);
        }

        let output = Command::new("ps")
            .args(["-p", &pid.to_string(), "-o", "command="])
            .output()
            .context("Failed to run ps")?;
        if !output.status.success() {
            return Ok(false);
        }
        let cmd = String::from_utf8_lossy(&output.stdout);
        Ok(cmd.contains("bv syftboxd") || cmd.contains("biovault"))
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

fn runtime_config(config: &Config) -> Result<SyftboxRuntimeConfig> {
    config.to_syftbox_runtime_config().map_err(|e| anyhow!(e))
}

pub async fn start(config: &Config, foreground: bool) -> Result<()> {
    // syftboxd is explicitly an embedded-host command; never fall back to spawning an external `syftbox` binary.
    std::env::set_var("BV_SYFTBOX_BACKEND", "embedded");

    if is_syftboxd_running(config)? {
        println!("❌ syftboxd is already running");
        return Ok(());
    }

    let biovault_dir = get_biovault_dir(config)?;
    std::fs::create_dir_all(&biovault_dir)
        .with_context(|| format!("Failed to create biovault directory: {:?}", biovault_dir))?;

    if foreground {
        let pid = std::process::id();

        let pid_file_path = if let Ok(path_str) = std::env::var("BV_SYFTBOXD_PID_FILE") {
            Some(PathBuf::from(path_str))
        } else {
            Some(get_pid_file_path(config)?)
        };

        if let Some(ref pid_path) = pid_file_path {
            std::fs::write(pid_path, pid.to_string())
                .with_context(|| format!("Failed to write PID file: {:?}", pid_path))?;
        }

        let status = SyftboxdStatus::new(pid);
        write_status(config, &status)?;

        let runtime = runtime_config(config)?;
        match start_syftbox(&runtime) {
            Ok(true) => {}
            Ok(false) => {}
            Err(e) => {
                return Err(e);
            }
        }

        // Keep running until we get a shutdown signal.
        #[cfg(unix)]
        {
            let mut term =
                tokio::signal::unix::signal(tokio::signal::unix::SignalKind::terminate())
                    .context("Failed to listen for SIGTERM")?;
            loop {
                tokio::select! {
                    _ = signal::ctrl_c() => break,
                    _ = term.recv() => break,
                    _ = tokio::time::sleep(Duration::from_secs(2)) => {
                        // Best-effort: if syftbox is not running anymore, try to restart it.
                        if !is_syftbox_running(&runtime).unwrap_or(false) {
                            let _ = start_syftbox(&runtime);
                        }
                    }
                }
            }
        }

        #[cfg(not(unix))]
        {
            loop {
                tokio::select! {
                    _ = signal::ctrl_c() => break,
                    _ = tokio::time::sleep(Duration::from_secs(2)) => {
                        if !is_syftbox_running(&runtime).unwrap_or(false) {
                            let _ = start_syftbox(&runtime);
                        }
                    }
                }
            }
        }

        let _ = stop_syftbox(&runtime);
        let _ = std::fs::remove_file(get_status_file_path(config)?);
        if let Some(ref pid_path) = pid_file_path {
            let _ = std::fs::remove_file(pid_path);
        }

        return Ok(());
    }

    // Background: spawn a foreground syftboxd process and detach.
    let config_json = serde_json::to_string(config).context("Failed to serialize config")?;
    let syftbox_data_dir = config.get_syftbox_data_dir()?;
    let current_exe = std::env::current_exe().context("Failed to get current executable path")?;

    let log_path = get_log_file_path(config)?;
    if let Some(log_dir) = log_path.parent() {
        std::fs::create_dir_all(log_dir)
            .with_context(|| format!("Failed to create logs directory: {:?}", log_dir))?;
    }
    let log_file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&log_path)
        .context("Failed to open log file for syftboxd spawn")?;
    let log_file_stderr = log_file.try_clone()?;

    let pid_path = get_pid_file_path(config)?;

    let mut child = Command::new(current_exe)
        .args(["syftboxd", "start", "--foreground"])
        .env("BV_SYFTBOXD_CONFIG", config_json)
        .env(
            "BV_SYFTBOXD_PID_FILE",
            pid_path.to_string_lossy().to_string(),
        )
        .env(
            "SYFTBOX_DATA_DIR",
            syftbox_data_dir.to_string_lossy().to_string(),
        )
        .stdin(Stdio::null())
        .stdout(Stdio::from(log_file))
        .stderr(Stdio::from(log_file_stderr))
        .spawn()
        .context("Failed to spawn syftboxd process")?;

    let pid = child.id();
    tokio::time::sleep(Duration::from_millis(750)).await;

    match child.try_wait() {
        Ok(Some(status)) => {
            return Err(anyhow!(
                "syftboxd exited immediately with status: {}",
                status
            ));
        }
        Ok(None) => {
            if pid_path.exists() {
                println!("✅ syftboxd started successfully (PID: {})", pid);
            } else {
                return Err(anyhow!("syftboxd started but PID file was not created"));
            }
        }
        Err(e) => {
            return Err(anyhow!("Failed to check syftboxd status: {}", e));
        }
    }

    Ok(())
}

pub async fn stop(config: &Config) -> Result<()> {
    if !is_syftboxd_running(config)? {
        println!("❌ syftboxd is not running");
        return Ok(());
    }

    let pid_path = get_pid_file_path(config)?;
    let pid_str = std::fs::read_to_string(&pid_path)?;
    let pid: u32 = pid_str.trim().parse().context("Invalid PID file")?;

    #[cfg(unix)]
    unsafe {
        let result = libc::kill(pid as i32, libc::SIGTERM);
        if result != 0 {
            return Err(anyhow!("Failed to stop syftboxd"));
        }
    }

    #[cfg(windows)]
    {
        let output = Command::new("taskkill")
            .args(["/F", "/PID", &pid.to_string()])
            .output()?;
        if !output.status.success() {
            return Err(anyhow!("Failed to stop syftboxd"));
        }
    }

    let _ = std::fs::remove_file(&pid_path);
    let _ = std::fs::remove_file(get_status_file_path(config)?);

    println!("✅ syftboxd stopped");
    Ok(())
}

pub async fn status(config: &Config) -> Result<()> {
    if is_syftboxd_running(config)? {
        println!("✓ syftboxd running");
    } else {
        println!("✗ syftboxd not running");
    }
    Ok(())
}

pub async fn logs(config: &Config, follow: bool, lines: Option<usize>) -> Result<()> {
    let log_path = get_log_file_path(config)?;
    if !log_path.exists() {
        println!("📝 No syftboxd log file found.");
        return Ok(());
    }

    if follow {
        let mut child = Command::new("tail")
            .args(["-f", &log_path.to_string_lossy()])
            .stdout(Stdio::piped())
            .spawn()
            .context("Failed to start tail")?;

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
        let output = Command::new("tail")
            .args(["-n", &tail_lines.to_string(), &log_path.to_string_lossy()])
            .output()
            .context("Failed to read log file")?;
        print!("{}", String::from_utf8_lossy(&output.stdout));
        let _ = std::io::stdout().flush();
    }

    Ok(())
}
