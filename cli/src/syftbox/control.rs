use crate::config::Config;
use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum SyftBoxMode {
    Sbenv,
    Direct,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftBoxState {
    pub running: bool,
    pub mode: SyftBoxMode,
}

pub fn detect_mode(config: &Config) -> Result<SyftBoxMode> {
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(if data_dir.join(".sbenv").exists() {
        SyftBoxMode::Sbenv
    } else {
        SyftBoxMode::Direct
    })
}

pub fn state(config: &Config) -> Result<SyftBoxState> {
    let mode = detect_mode(config)?;
    let running = is_running_with_mode(config, mode)?;
    Ok(SyftBoxState { running, mode })
}

pub fn is_syftbox_running(config: &Config) -> Result<bool> {
    let mode = detect_mode(config)?;
    is_running_with_mode(config, mode)
}

pub fn start_syftbox(config: &Config) -> Result<bool> {
    let mode = detect_mode(config)?;
    if is_running_with_mode(config, mode)? {
        return Ok(false);
    }

    match mode {
        SyftBoxMode::Sbenv => start_with_sbenv(config)?,
        SyftBoxMode::Direct => start_direct(config)?,
    }

    if !wait_for(
        || is_running_with_mode(config, mode),
        true,
        Duration::from_secs(5),
    )? {
        return Err(anyhow!("SyftBox did not start in time"));
    }

    Ok(true)
}

pub fn stop_syftbox(config: &Config) -> Result<bool> {
    let mode = detect_mode(config)?;
    let pids = running_pids(config, mode)?;
    if pids.is_empty() {
        return Ok(false);
    }

    match mode {
        SyftBoxMode::Sbenv => stop_with_sbenv(config)?,
        SyftBoxMode::Direct => stop_direct(&pids)?,
    }

    if !wait_for(
        || is_running_with_mode(config, mode),
        false,
        Duration::from_secs(5),
    )? {
        return Err(anyhow!("SyftBox did not stop in time"));
    }

    Ok(true)
}

fn wait_for<F>(mut check: F, expected: bool, timeout: Duration) -> Result<bool>
where
    F: FnMut() -> Result<bool>,
{
    let deadline = Instant::now() + timeout;
    loop {
        let current = check()?;
        if current == expected {
            return Ok(true);
        }
        if Instant::now() >= deadline {
            return Ok(false);
        }
        thread::sleep(Duration::from_millis(250));
    }
}

fn start_with_sbenv(config: &Config) -> Result<()> {
    let data_dir = config.get_syftbox_data_dir()?;
    let status = Command::new("sbenv")
        .arg("start")
        .arg("--skip-login-check")
        .current_dir(&data_dir)
        .status()
        .context("Failed to execute sbenv start")?;

    if !status.success() {
        return Err(anyhow!("sbenv start exited with status {}", status));
    }

    Ok(())
}

fn stop_with_sbenv(config: &Config) -> Result<()> {
    let data_dir = config.get_syftbox_data_dir()?;
    let status = Command::new("sbenv")
        .arg("stop")
        .current_dir(&data_dir)
        .status()
        .context("Failed to execute sbenv stop")?;

    if !status.success() {
        return Err(anyhow!("sbenv stop exited with status {}", status));
    }

    Ok(())
}

fn start_direct(config: &Config) -> Result<()> {
    let config_path = config.get_syftbox_config_path()?;
    let binary_path = resolve_syftbox_binary(config)?;
    eprintln!("🔧 Requested SyftBox binary: {}", binary_path.display());

    if !config_path.exists() {
        return Err(anyhow!(
            "SyftBox config file does not exist: {}",
            config_path.display()
        ));
    }

    eprintln!("📄 Using SyftBox config: {}", config_path.display());

    let mut child = Command::new(&binary_path)
        .arg("-c")
        .arg(&config_path)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .with_context(|| {
            format!(
                "Failed to spawn syftbox process using '{}'",
                binary_path.display()
            )
        })?;

    thread::sleep(Duration::from_secs(2));

    if let Some(status) = child
        .try_wait()
        .context("Failed to check syftbox child status")?
    {
        if status.success() {
            return Ok(());
        }
        return Err(anyhow!("SyftBox exited immediately with status {}", status));
    }

    std::mem::forget(child);
    Ok(())
}

fn stop_direct(pids: &[u32]) -> Result<()> {
    for pid in pids {
        let mut cmd = Command::new("kill");
        cmd.arg("-TERM").arg(pid.to_string());
        let status = cmd
            .status()
            .with_context(|| format!("Failed to send TERM to process {}", pid))?;
        if !status.success() {
            return Err(anyhow!("Failed to terminate syftbox process {}", pid));
        }
    }
    Ok(())
}

fn is_running_with_mode(config: &Config, mode: SyftBoxMode) -> Result<bool> {
    Ok(!running_pids(config, mode)?.is_empty())
}

fn running_pids(config: &Config, mode: SyftBoxMode) -> Result<Vec<u32>> {
    #[cfg(unix)]
    {
        let output = Command::new("ps")
            .arg("aux")
            .output()
            .context("Failed to execute ps command")?;

        if !output.status.success() {
            return Err(anyhow!("ps command failed"));
        }

        let ps_output = String::from_utf8_lossy(&output.stdout);

        let config_path = config.get_syftbox_config_path()?;
        let data_dir = config.get_syftbox_data_dir()?;
        let config_str = config_path.to_string_lossy();
        let data_dir_str = data_dir.to_string_lossy();

        let mut pids = Vec::new();
        for line in ps_output.lines() {
            if !line.contains("syftbox") {
                continue;
            }

            let matches_mode = match mode {
                SyftBoxMode::Sbenv => line.contains(data_dir_str.as_ref()),
                SyftBoxMode::Direct => {
                    line.contains(config_str.as_ref()) || line.contains(data_dir_str.as_ref())
                }
            };

            if !matches_mode {
                continue;
            }

            if let Some(pid) = parse_pid(line) {
                pids.push(pid);
            }
        }

        Ok(pids)
    }

    #[cfg(not(unix))]
    {
        let _ = config;
        let _ = mode;
        Err(anyhow!(
            "SyftBox process inspection is only supported on Unix-like platforms"
        ))
    }
}

fn resolve_syftbox_binary(config: &Config) -> Result<PathBuf> {
    if let Some(path) = config.get_binary_path("syftbox").and_then(|p| {
        let trimmed = p.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(PathBuf::from(trimmed))
        }
    }) {
        if path.is_absolute() && !path.exists() {
            return Err(anyhow!(
                "Configured SyftBox binary not found at {}",
                path.display()
            ));
        }
        eprintln!("ℹ️  Using configured SyftBox binary from config");
        return Ok(path);
    }

    if let Ok(env_path) = env::var("SYFTBOX_BINARY") {
        let path = PathBuf::from(env_path.trim());
        if path.is_absolute() && !path.exists() {
            return Err(anyhow!(
                "SYFTBOX_BINARY points to missing path: {}",
                path.display()
            ));
        }
        eprintln!("ℹ️  Using SyftBox binary from SYFTBOX_BINARY env var");
        return Ok(path);
    }

    if let Some(path) = find_syftbox_in_sbenv() {
        eprintln!("ℹ️  Detected SyftBox in ~/.sbenv: {}", path.display());
        return Ok(path);
    }

    eprintln!("ℹ️  No custom SyftBox path found; falling back to 'syftbox' in PATH");
    Ok(PathBuf::from("syftbox"))
}

fn find_syftbox_in_sbenv() -> Option<PathBuf> {
    let home = dirs::home_dir()?;
    let binaries_dir = home.join(".sbenv").join("binaries");

    if !binaries_dir.exists() {
        return None;
    }

    let mut candidates = Vec::new();

    if let Ok(entries) = fs::read_dir(&binaries_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                let syftbox_path = path.join("syftbox");
                if syftbox_path.is_file() {
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::PermissionsExt;
                        if let Ok(metadata) = syftbox_path.metadata() {
                            if metadata.permissions().mode() & 0o111 != 0 {
                                candidates.push(syftbox_path);
                            }
                        }
                    }
                    #[cfg(not(unix))]
                    {
                        candidates.push(syftbox_path);
                    }
                }
            }
        }
    }

    if candidates.is_empty() {
        return None;
    }

    candidates.sort_by(|a, b| {
        let a_parent = a
            .parent()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().into_owned());
        let b_parent = b
            .parent()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().into_owned());
        b_parent.cmp(&a_parent)
    });

    candidates.into_iter().next()
}

#[cfg(unix)]
fn parse_pid(line: &str) -> Option<u32> {
    line.split_whitespace()
        .nth(1)
        .and_then(|pid| pid.parse::<u32>().ok())
}

#[cfg(not(unix))]
fn parse_pid(_line: &str) -> Option<u32> {
    None
}

pub fn syftbox_paths(config: &Config) -> Result<(PathBuf, PathBuf)> {
    let config_path = config.get_syftbox_config_path()?;
    let data_dir = config.get_syftbox_data_dir()?;
    Ok((config_path, data_dir))
}
