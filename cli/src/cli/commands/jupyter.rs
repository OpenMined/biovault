use crate::data::BioVaultDb;
use crate::error::Result;
use anyhow::anyhow;
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::net::TcpListener;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;
use std::env;
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::time::{Duration, SystemTime};
use tokio::net::TcpStream;
use tracing::{info, warn};

const PYFORY_MACOS_INTEL_WHEEL: &str = "https://files.pythonhosted.org/packages/35/c5/b2de2a2dc0d2b74002924cdd46a6e6d3bccc5380181ca0dc850855608bfe/pyfory-0.13.2-cp312-cp312-macosx_10_13_x86_64.whl";

fn hide_console_window(_cmd: &mut Command) {
    #[cfg(target_os = "windows")]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x08000000;
        _cmd.creation_flags(CREATE_NO_WINDOW);
    }
}

fn is_pid_alive(pid: i32) -> bool {
    if pid <= 0 {
        return false;
    }

    if cfg!(windows) {
        // Use CSV format for more reliable parsing across locales.
        let filter = format!("PID eq {}", pid);
        let mut cmd = Command::new("tasklist");
        cmd.args(["/FI", &filter, "/FO", "CSV", "/NH"]);
        hide_console_window(&mut cmd);
        let output = cmd.output();

        let Ok(output) = output else {
            return false;
        };

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stdout = stdout.trim();
        if stdout.is_empty() {
            return false;
        }

        // When no task matches, tasklist prints an INFO line (often with exit code 0).
        if stdout.starts_with("INFO:") {
            return false;
        }

        let pid_needle = format!(",\"{}\",", pid);
        stdout.lines().any(|line| line.contains(&pid_needle))
    } else {
        Command::new("kill")
            .args(["-0", &pid.to_string()])
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
    }
}

fn terminate_pid(pid: i32) -> Result<()> {
    if pid <= 0 || !is_pid_alive(pid) {
        return Ok(());
    }

    if cfg!(windows) {
        let mut cmd = Command::new("taskkill");
        cmd.args(["/PID", &pid.to_string(), "/T", "/F"]);
        hide_console_window(&mut cmd);
        let output = cmd.output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            return Err(anyhow!(
                "Failed to stop Jupyter process (PID: {}) via taskkill (status: {}). stdout='{}' stderr='{}'",
                pid,
                output.status,
                stdout.trim(),
                stderr.trim()
            )
            .into());
        }

        Ok(())
    } else {
        let kill_status = Command::new("kill").arg(pid.to_string()).status();
        match kill_status {
            Ok(status) if status.success() => Ok(()),
            _ => {
                let status = Command::new("kill")
                    .args(["-9", &pid.to_string()])
                    .status()?;
                if status.success() {
                    Ok(())
                } else {
                    Err(anyhow!("Failed to stop Jupyter process (PID: {})", pid).into())
                }
            }
        }
    }
}

async fn wait_for_pid_exit(pid: i32, timeout: Duration) -> bool {
    let deadline = SystemTime::now() + timeout;
    while SystemTime::now() < deadline {
        if !is_pid_alive(pid) {
            return true;
        }
        tokio::time::sleep(Duration::from_millis(200)).await;
    }
    !is_pid_alive(pid)
}

async fn remove_dir_all_with_retry(path: &Path, max_wait: Duration) -> io::Result<()> {
    let deadline = SystemTime::now() + max_wait;
    let mut last_err: Option<io::Error> = None;

    while SystemTime::now() < deadline {
        match fs::remove_dir_all(path) {
            Ok(()) => return Ok(()),
            Err(err) => {
                let raw = err.raw_os_error();
                let retryable = matches!(raw, Some(5) | Some(32) | Some(145))
                    || err.kind() == io::ErrorKind::PermissionDenied;
                if !retryable {
                    return Err(err);
                }
                last_err = Some(err);
                tokio::time::sleep(Duration::from_millis(300)).await;
            }
        }
    }

    Err(last_err.unwrap_or_else(|| io::Error::other("Timed out while removing directory")))
}

fn format_process_output(output: &std::process::Output) -> String {
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    let stdout_trim = stdout.trim();
    let stderr_trim = stderr.trim();

    let mut msg = String::new();
    if !stdout_trim.is_empty() {
        msg.push_str("stdout:\n");
        msg.push_str(stdout_trim);
        msg.push('\n');
    }
    if !stderr_trim.is_empty() {
        msg.push_str("stderr:\n");
        msg.push_str(stderr_trim);
        msg.push('\n');
    }
    msg
}

fn output_indicates_windows_file_lock(output: &std::process::Output) -> bool {
    let stdout = String::from_utf8_lossy(&output.stdout).to_lowercase();
    let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
    let combined = format!("{}\n{}", stdout, stderr);

    combined.contains("os error 32")
        || combined.contains("os error 5")
        || combined.contains("the process cannot access the file")
        || combined.contains("access is denied")
        || combined.contains("being used by another process")
}

fn best_effort_stop_jupyter_for_module(module_dir: &Path) {
    // Try DB PID first (may be a wrapper PID), then runtime jpserver PID(s).
    if let Ok(db) = BioVaultDb::new() {
        if let Ok(canonical) = module_dir.canonicalize() {
            if let Ok(Some(env)) = db.get_dev_env(canonical.to_string_lossy().as_ref()) {
                if let Some(pid) = env.jupyter_pid {
                    let _ = terminate_pid(pid);
                }
            }
        }
    }

    let runtime_dir = module_dir.join(".jupyter-runtime");
    for pid in find_module_runtime_pids(&runtime_dir) {
        let _ = terminate_pid(pid);
    }
}

fn resolve_uv_path() -> Result<String> {
    if let Ok(env_path) = env::var("BIOVAULT_BUNDLED_UV") {
        if !env_path.trim().is_empty() {
            let p = PathBuf::from(env_path.trim());
            if p.exists() {
                println!(
                    "🔧 Using bundled uv from BIOVAULT_BUNDLED_UV={}",
                    p.display()
                );
                info!("Using bundled uv: {}", p.display());
                return Ok(p.display().to_string());
            } else {
                warn!("BIOVAULT_BUNDLED_UV is set but missing: {}", p.display());
                return Err(anyhow!(
                    "BIOVAULT_BUNDLED_UV is set to '{}' but the file does not exist",
                    p.display()
                )
                .into());
            }
        }
    }

    if let Ok(p) = which::which("uv") {
        let s = p.display().to_string();
        println!("🔧 Using uv from PATH: {}", s);
        info!("Using uv from PATH: {}", s);
        return Ok(s);
    }

    warn!("uv not found in BIOVAULT_BUNDLED_UV or PATH");
    Err(anyhow!(
        "uv not found. Set BIOVAULT_BUNDLED_UV to the bundled binary or install uv on PATH."
    )
    .into())
}

fn is_macos_intel() -> bool {
    cfg!(target_os = "macos") && cfg!(target_arch = "x86_64")
}

fn ensure_virtualenv(module_dir: &Path, python_version: &str, uv_bin: &str) -> Result<()> {
    let venv_path = module_dir.join(".venv");
    let capture_output = std::env::var_os("BIOVAULT_DESKTOP_LOG_FILE").is_some();

    // Version of biovault-beaver from PyPI - auto-detected from submodule at compile time
    // Falls back to hardcoded version if env var not set (e.g., when building biovault CLI standalone)
    let beaver_version = option_env!("BEAVER_VERSION").unwrap_or("0.1.30");

    // Marker file to track if dependencies are already installed
    // Includes version to trigger reinstall on beaver version changes
    let deps_marker = venv_path.join(format!(".deps-installed-{}", beaver_version));
    let pyfory_marker = venv_path.join(format!(".deps-installed-{}-pyfory-x86_64", beaver_version));
    let needs_base_install = !venv_path.exists() || !deps_marker.exists();
    let needs_pyfory_fix = is_macos_intel() && (needs_base_install || !pyfory_marker.exists());

    // Check if venv exists AND deps are already installed
    if !needs_base_install && !needs_pyfory_fix {
        println!("✅ Using existing virtualenv with dependencies");
        // Ensure uv is available inside the virtualenv PATH
        link_uv_into_venv(&venv_path, uv_bin);
        return Ok(());
    }

    if !venv_path.exists() {
        println!("📦 Creating virtualenv with Python {}...", python_version);

        if capture_output {
            let mut cmd = Command::new(uv_bin);
            cmd.args(["venv", "--python", python_version, ".venv"]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let output = cmd.output()?;
            if !output.status.success() {
                return Err(anyhow!(
                    "Failed to create virtualenv. Try: bv python install {}\n{}",
                    python_version,
                    format_process_output(&output)
                )
                .into());
            }
        } else {
            let mut cmd = Command::new(uv_bin);
            cmd.args(["venv", "--python", python_version, ".venv"]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let status = cmd.status()?;
            if !status.success() {
                return Err(anyhow!(
                    "Failed to create virtualenv. Try: bv python install {}",
                    python_version
                )
                .into());
            }
        }
    } else if needs_base_install {
        println!("📦 Virtualenv exists but dependencies need install/update...");
    }

    // Ensure uv is available inside the virtualenv PATH (symlink/copy bundled uv)
    link_uv_into_venv(&venv_path, uv_bin);

    if needs_base_install {
        println!(
            "📦 Installing packages: jupyterlab cleon biovault-beaver[lib-support]=={}",
            beaver_version
        );

        // Install base packages from PyPI including pinned biovault-beaver
        // Note: syftbox-sdk is a direct dependency of biovault-beaver, not an extra
        let beaver_pkg = format!("biovault-beaver[lib-support]=={}", beaver_version);
        if capture_output {
            let mut cmd = Command::new(uv_bin);
            cmd.args([
                "pip",
                "install",
                "--python",
                ".venv",
                "jupyterlab",
                "cleon",
                &beaver_pkg,
            ]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let output = cmd.output()?;

            if !output.status.success() {
                if cfg!(windows) && output_indicates_windows_file_lock(&output) {
                    // Windows often locks entrypoint exes (like .venv/Scripts/jupyter.exe) while Jupyter is running.
                    // Try a best-effort stop and retry once.
                    best_effort_stop_jupyter_for_module(module_dir);
                    std::thread::sleep(Duration::from_millis(500));

                    let mut cmd = Command::new(uv_bin);
                    cmd.args([
                        "pip",
                        "install",
                        "--python",
                        ".venv",
                        "jupyterlab",
                        "cleon",
                        &beaver_pkg,
                    ]);
                    cmd.current_dir(module_dir);
                    hide_console_window(&mut cmd);
                    let retry = cmd.output()?;

                    if retry.status.success() {
                        // Proceed with the rest of environment setup (DEV overlays, marker file, etc.)
                        // now that package installation succeeded.
                    } else {
                        return Err(anyhow!(
                            "Failed to install required Python packages (jupyterlab/cleon/biovault-beaver) after retrying. {}\n{}",
                            "Ensure all Jupyter/Python processes for this session are stopped, then retry.",
                            format_process_output(&retry)
                        )
                        .into());
                    }
                }

                return Err(anyhow!(
                    "Failed to install required Python packages (jupyterlab/cleon/biovault-beaver). If you see file access errors on Windows, stop all running Jupyter/Python processes and retry.\n{}",
                    format_process_output(&output)
                )
                .into());
            }
        } else {
            let mut cmd = Command::new(uv_bin);
            cmd.args([
                "pip",
                "install",
                "--python",
                ".venv",
                "jupyterlab",
                "cleon",
                &beaver_pkg,
            ]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let status = cmd.status()?;

            if !status.success() {
                return Err(anyhow!(
                    "Failed to install required Python packages (jupyterlab/cleon/biovault-beaver)"
                )
                .into());
            }
        }
    }

    // Path structure: module_dir is BIOVAULT_HOME/sessions/<session_id>
    // BIOVAULT_HOME is like: workspace3/biovault/sandbox/client1@sandbox.local
    // biovault submodule root is at BIOVAULT_HOME/../.. (e.g., workspace3/biovault)
    // So from module_dir (BIOVAULT_HOME/sessions/<id>): 4 levels up
    let biovault_root = module_dir.join("..").join("..").join("..").join("..");

    // DEV MODE: If local source exists, install editable versions on top of PyPI packages
    // This allows developers to test local changes while maintaining compatibility
    // Note: biovault_root IS the biovault submodule, so paths are directly under it
    let syftbox_path = biovault_root.join("syftbox-sdk").join("python");
    let beaver_path = biovault_root.join("biovault-beaver").join("python");

    if needs_base_install && (syftbox_path.exists() || beaver_path.exists()) {
        println!("🔧 DEV MODE: Local source detected, installing editable packages...");

        // Install syftbox-sdk first (beaver depends on it)
        // Prefer pre-built wheel to avoid slow maturin rebuilds
        if syftbox_path.exists() {
            let wheels_dir = syftbox_path.join("target").join("wheels");
            let wheel_exists = wheels_dir.exists()
                && std::fs::read_dir(&wheels_dir)
                    .map(|mut d| {
                        d.any(|e| {
                            e.map(|e| {
                                e.path()
                                    .extension()
                                    .map(|ext| ext == "whl")
                                    .unwrap_or(false)
                            })
                            .unwrap_or(false)
                        })
                    })
                    .unwrap_or(false);

            if wheel_exists {
                // Use pre-built wheel (much faster than editable install)
                println!("📦 Installing syftbox-sdk from pre-built wheel...");
                let wheels_canonical = wheels_dir.canonicalize().unwrap_or(wheels_dir);
                let mut cmd = Command::new(uv_bin);
                cmd.args([
                    "pip",
                    "install",
                    "--python",
                    ".venv",
                    "--find-links",
                    wheels_canonical.to_str().unwrap_or("."),
                    "--reinstall-package",
                    "syftbox-sdk",
                    "syftbox-sdk",
                ]);
                cmd.current_dir(module_dir);
                hide_console_window(&mut cmd);
                let status = cmd.status()?;

                if status.success() {
                    println!(
                        "✅ syftbox-sdk installed from wheel: {}",
                        wheels_canonical.display()
                    );
                } else {
                    println!(
                        "⚠️ Failed to install syftbox-sdk from wheel, falling back to editable..."
                    );
                    // Fall back to editable install
                    let syftbox_canonical =
                        syftbox_path.canonicalize().unwrap_or(syftbox_path.clone());
                    let mut cmd = Command::new(uv_bin);
                    cmd.args([
                        "pip",
                        "install",
                        "--python",
                        ".venv",
                        "-e",
                        syftbox_canonical.to_str().unwrap_or("."),
                    ]);
                    cmd.current_dir(module_dir);
                    hide_console_window(&mut cmd);
                    let _ = cmd.status();
                }
            } else {
                // No wheel found, use editable install (slower, triggers maturin)
                println!("📦 Installing syftbox-sdk from source (no pre-built wheel found)...");
                let syftbox_canonical = syftbox_path.canonicalize().unwrap_or(syftbox_path.clone());
                let mut cmd = Command::new(uv_bin);
                cmd.args([
                    "pip",
                    "install",
                    "--python",
                    ".venv",
                    "-e",
                    syftbox_canonical.to_str().unwrap_or("."),
                ]);
                cmd.current_dir(module_dir);
                hide_console_window(&mut cmd);
                let status = cmd.status()?;

                if status.success() {
                    println!(
                        "✅ syftbox-sdk installed from: {}",
                        syftbox_canonical.display()
                    );
                } else {
                    println!("⚠️ Failed to install syftbox-sdk from local path");
                }
            }
        }

        // Install beaver from local editable path (overwrites PyPI version)
        if beaver_path.exists() {
            println!("🦫 Installing beaver from local editable path (overwriting PyPI version)...");
            let beaver_canonical = beaver_path.canonicalize().unwrap_or(beaver_path);
            let beaver_with_extras =
                format!("{}[lib-support]", beaver_canonical.to_str().unwrap_or("."));
            let mut cmd = Command::new(uv_bin);
            cmd.args([
                "pip",
                "install",
                "--python",
                ".venv",
                "-e",
                &beaver_with_extras,
            ]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let status = cmd.status()?;

            if status.success() {
                println!("✅ beaver installed from: {}", beaver_canonical.display());
            } else {
                println!("⚠️ Failed to install beaver from local path");
            }
        }

        println!("✅ Virtualenv ready with jupyterlab, cleon, and DEV beaver/syftbox-sdk");
    } else if needs_base_install {
        println!(
            "✅ Virtualenv ready with jupyterlab, cleon, biovault-beaver=={}, and syftbox-sdk",
            beaver_version
        );
    }

    if needs_pyfory_fix {
        println!("🛠️ macOS Intel detected; forcing pyfory reinstall...");
        if capture_output {
            let mut cmd = Command::new(uv_bin);
            cmd.args([
                "pip",
                "install",
                "--python",
                ".venv",
                "--force-reinstall",
                PYFORY_MACOS_INTEL_WHEEL,
            ]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let output = cmd.output()?;
            if !output.status.success() {
                return Err(anyhow!(
                    "Failed to force-reinstall pyfory (macOS Intel). Try reinstalling Jupyter dependencies.\n{}",
                    format_process_output(&output)
                )
                .into());
            }
        } else {
            let mut cmd = Command::new(uv_bin);
            cmd.args([
                "pip",
                "install",
                "--python",
                ".venv",
                "--force-reinstall",
                PYFORY_MACOS_INTEL_WHEEL,
            ]);
            cmd.current_dir(module_dir);
            hide_console_window(&mut cmd);
            let status = cmd.status()?;
            if !status.success() {
                return Err(anyhow!(
                    "Failed to force-reinstall pyfory (macOS Intel). Try reinstalling Jupyter dependencies."
                )
                .into());
            }
        }

        if let Err(e) = fs::write(&pyfory_marker, "pyfory") {
            warn!("Failed to create pyfory marker file: {}", e);
        }
    }

    // Create marker file to skip reinstall on next launch
    if needs_base_install {
        if let Err(e) = fs::write(&deps_marker, beaver_version) {
            warn!("Failed to create deps marker file: {}", e);
        }
    }

    Ok(())
}

/// Place the resolved uv binary into the virtualenv's bin directory for notebook use.
fn link_uv_into_venv(venv_path: &Path, uv_bin: &str) {
    let uv_src = PathBuf::from(uv_bin);
    if !uv_src.exists() {
        println!(
            "⚠️  uv binary not found at {}; skipping venv link",
            uv_src.display()
        );
        return;
    }

    let target = if cfg!(windows) {
        venv_path.join("Scripts").join("uv.exe")
    } else {
        venv_path.join("bin").join("uv")
    };
    if target.exists() {
        return;
    }

    #[cfg(unix)]
    {
        if let Err(e) = std::os::unix::fs::symlink(&uv_src, &target) {
            println!(
                "⚠️  Failed to symlink uv into venv ({} -> {}): {}",
                uv_src.display(),
                target.display(),
                e
            );
        } else {
            println!(
                "✅ Linked uv into venv: {} -> {}",
                target.display(),
                uv_src.display()
            );
        }
    }

    #[cfg(windows)]
    {
        if let Err(e) = std::fs::copy(&uv_src, &target) {
            println!(
                "⚠️  Failed to copy uv into venv ({} -> {}): {}",
                uv_src.display(),
                target.display(),
                e
            );
        } else {
            println!(
                "✅ Copied uv into venv: {} -> {}",
                target.display(),
                uv_src.display()
            );
        }
    }
}

#[derive(Clone, Debug)]
struct JupyterRuntimeInfo {
    port: Option<i32>,
    url: Option<String>,
    token: Option<String>,
}

fn candidate_runtime_dirs() -> Vec<PathBuf> {
    let mut dirs = Vec::new();

    if let Ok(dir) = env::var("JUPYTER_RUNTIME_DIR") {
        if !dir.is_empty() {
            dirs.push(PathBuf::from(dir));
        }
    }

    if let Some(home) = dirs::home_dir() {
        dirs.push(home.join(".local/share/jupyter/runtime"));
        dirs.push(home.join("Library/Jupyter/runtime"));
    }

    if cfg!(target_os = "windows") {
        if let Ok(appdata) = env::var("APPDATA") {
            dirs.push(PathBuf::from(appdata).join("jupyter/runtime"));
        }
    }

    dirs.push(std::env::temp_dir());
    dirs.push(PathBuf::from("/tmp"));

    dirs.sort();
    dirs.dedup();
    dirs
}

fn parse_runtime_file(path: &Path) -> Option<(JupyterRuntimeInfo, Option<i32>)> {
    let content = fs::read_to_string(path).ok()?;
    let value: Value = serde_json::from_str(&content).ok()?;
    let port = value.get("port").and_then(|v| v.as_i64()).map(|v| v as i32);
    let url = value
        .get("url")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());
    let token = value
        .get("token")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());
    let pid = value.get("pid").and_then(|v| v.as_i64()).map(|v| v as i32);
    Some((JupyterRuntimeInfo { port, url, token }, pid))
}

fn find_module_runtime_pids(runtime_dir: &Path) -> Vec<i32> {
    if !runtime_dir.exists() {
        return Vec::new();
    }

    let mut pids = BTreeSet::<i32>::new();
    let Ok(entries) = fs::read_dir(runtime_dir) else {
        return Vec::new();
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
            continue;
        }
        let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        if !file_name.starts_with("jpserver-") {
            continue;
        }

        let Some((_info, pid_opt)) = parse_runtime_file(&path) else {
            continue;
        };

        if let Some(pid) = pid_opt {
            pids.insert(pid);
        }
    }

    pids.into_iter().collect()
}

fn find_runtime_info(pid: u32) -> Option<JupyterRuntimeInfo> {
    let mut latest: Option<((bool, SystemTime), JupyterRuntimeInfo)> = None;

    for dir in candidate_runtime_dirs() {
        if !dir.exists() {
            continue;
        }

        if let Ok(entries) = fs::read_dir(&dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
                    continue;
                }

                if let Some(file_name) = path.file_name().and_then(|name| name.to_str()) {
                    if !file_name.starts_with("jpserver-") {
                        continue;
                    }
                } else {
                    continue;
                }

                let (info, info_pid_opt) = match parse_runtime_file(&path) {
                    Some(result) => result,
                    None => continue,
                };

                // Instead of matching exact PID (which fails with uv run wrapper),
                // find the most recent jpserver file within the last 2 minutes.
                //
                // If the runtime file includes a PID, use it to prefer an exact match,
                // but don't hard-require it (uv may wrap/spawn a different PID).
                if let Ok(metadata) = entry.metadata() {
                    if let Ok(modified) = metadata.modified() {
                        if let Ok(elapsed) = modified.elapsed() {
                            if elapsed.as_secs() > 120 {
                                continue;
                            }
                        }

                        let prefers_exact_pid =
                            info_pid_opt.map(|p| p as u32 == pid).unwrap_or(false);
                        let score = (prefers_exact_pid, modified);

                        match &latest {
                            Some((best_score, _)) if score <= *best_score => {}
                            _ => latest = Some((score, info.clone())),
                        }
                    }
                }
            }
        }
    }

    latest.map(|(_, info)| info)
}

/// Parse Jupyter URL from output line like:
/// "[I 2025-10-22 22:57:28.681 ServerApp] http://localhost:8888/lab?token=abc123"
fn parse_jupyter_url_from_line(line: &str) -> Option<(i32, String, String)> {
    // Look for lines containing "http://localhost:" or "http://127.0.0.1:"
    if let Some(start) = line.find("http://") {
        let url_part = &line[start..];
        // Extract URL until whitespace or end of line
        let url = url_part.split_whitespace().next()?;

        // Extract port from "http://localhost:8888" or "http://127.0.0.1:8888"
        let port = if let Some(port_start) = url.find("://localhost:").or(url.find("://127.0.0.1:"))
        {
            let after_host = &url[port_start + 13..]; // Skip "://localhost:" or "://127.0.0.1:"
            after_host
                .split(&['/', '?'][..])
                .next()?
                .parse::<i32>()
                .ok()?
        } else {
            8888 // Default port
        };

        // Extract token from "?token=abc123"
        let token = if let Some(token_start) = url.find("?token=") {
            let after_token = &url[token_start + 7..]; // Skip "?token="
            after_token.split('&').next()?.to_string()
        } else {
            return None;
        };

        return Some((port, url.to_string(), token));
    }
    None
}

fn find_available_port() -> Option<i32> {
    // Bind to port 0 to let the OS pick an ephemeral port, then release it.
    // Retry a few times in case of transient failures.
    for _ in 0..10 {
        if let Ok(listener) = TcpListener::bind("127.0.0.1:0") {
            if let Ok(addr) = listener.local_addr() {
                return Some(addr.port() as i32);
            }
        }
    }
    None
}

async fn wait_for_server_ready(port: i32) -> Result<()> {
    let port_u16 = u16::try_from(port).map_err(|_| anyhow!("Invalid Jupyter port: {}", port))?;

    for _ in 0..60 {
        if TcpStream::connect(("127.0.0.1", port_u16)).await.is_ok() {
            return Ok(());
        }

        tokio::time::sleep(Duration::from_millis(500)).await;
    }

    Err(anyhow!("Timed out waiting for Jupyter to start on port {}", port).into())
}

pub async fn start(module_path: &str, python_version: &str) -> Result<()> {
    let uv_bin = resolve_uv_path()?;

    // Check if module_path is a number (list index)
    let module_dir = if let Ok(index) = module_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Module index must be >= 1").into());
        }

        let modules = get_modules_with_venvs()?;

        if modules.is_empty() {
            return Err(anyhow!("No modules with virtualenvs found. Run 'bv jupyter list' to see available modules.").into());
        }

        if index > modules.len() {
            return Err(anyhow!(
                "Module index {} out of range. Only {} module(s) available.",
                index,
                modules.len()
            )
            .into());
        }

        modules[index - 1].clone()
    } else {
        PathBuf::from(module_path)
    };

    if !module_dir.exists() {
        return Err(anyhow!("Module directory does not exist: {}", module_dir.display()).into());
    }

    let venv_path = module_dir.join(".venv");

    let db = BioVaultDb::new()?;

    // Check if Jupyter is already running for this module
    let canonical_path = module_dir.canonicalize()?;
    if let Some(env) = db.get_dev_env(canonical_path.to_str().unwrap())? {
        if let Some(pid) = env.jupyter_pid {
            // Check if process is still alive
            let is_alive = is_pid_alive(pid);

            if is_alive {
                println!("✅ Jupyter Lab already running (PID: {})", pid);
                if let Some(url) = env.jupyter_url.as_ref() {
                    println!("   Access at: {}", url);
                } else if let Some(port) = env.jupyter_port {
                    println!("   Access at: http://localhost:{}", port);
                } else {
                    println!("   Access at: <unknown>");
                }
                println!("   Use 'bv jupyter stop' with module path or index to stop it first");
                return Ok(());
            } else {
                // Clear stale session info
                db.update_jupyter_session(&module_dir, None, None, None, None)?;
            }
        }
    }

    ensure_virtualenv(&module_dir, python_version, &uv_bin)?;

    // Register/update in database
    db.register_dev_env(&module_dir, python_version, "jupyter", true)?;

    // Launch Jupyter Lab
    let chosen_port = find_available_port();

    // Isolate Jupyter runtime files per session to avoid cross-instance detection
    let runtime_dir = module_dir.join(".jupyter-runtime");
    let _ = fs::create_dir_all(&runtime_dir);
    std::env::set_var("JUPYTER_RUNTIME_DIR", &runtime_dir);
    // Also point XDG_RUNTIME_DIR to keep jpserver files local (macOS/Linux)
    std::env::set_var("XDG_RUNTIME_DIR", &runtime_dir);

    println!("🚀 Launching Jupyter Lab with: uv run --python .venv jupyter lab");
    if let Some(port) = chosen_port {
        println!("🎯 Requested port: {} (random to reduce conflicts)", port);
    }

    if !venv_path
        .join(if cfg!(windows) {
            "Scripts/jupyter.exe"
        } else {
            "bin/jupyter"
        })
        .exists()
    {
        return Err(anyhow!(
            "Jupyter not found in virtualenv. Try: bv jupyter reset {}",
            module_path
        )
        .into());
    }

    use std::process::Stdio;

    let venv_bin_dir = if cfg!(windows) {
        venv_path.join("Scripts")
    } else {
        venv_path.join("bin")
    };
    let existing_path = std::env::var("PATH").unwrap_or_default();
    let path_sep = if cfg!(windows) { ";" } else { ":" };
    let combined_path = if existing_path.is_empty() {
        venv_bin_dir.to_string_lossy().to_string()
    } else {
        format!(
            "{}{}{}",
            venv_bin_dir.to_string_lossy(),
            path_sep,
            existing_path
        )
    };

    let mut args: Vec<String> = vec![
        "run",
        "--python",
        ".venv",
        "jupyter",
        "lab",
        "--no-browser",
        "--ServerApp.token=",
        "--ServerApp.password=",
        "--ServerApp.disable_check_xsrf=true",
        "--ServerApp.allow_origin=*",
    ]
    .into_iter()
    .map(String::from)
    .collect();

    if let Some(port) = chosen_port {
        args.push("--port".into());
        args.push(port.to_string());
        // Avoid retries so we don't silently hop to a different port under race conditions.
        args.push("--ServerApp.port_retries=0".into());
    }

    let mut cmd = Command::new(&uv_bin);
    cmd.args(&args)
        .current_dir(&module_dir)
        .env("JUPYTER_RUNTIME_DIR", &runtime_dir)
        .env("XDG_RUNTIME_DIR", &runtime_dir)
        .env("VIRTUAL_ENV", &venv_path)
        .env("PATH", &combined_path)
        // Ensure the runner's Python environment doesn't bleed into uv/venv processes.
        .env_remove("PYTHONHOME")
        .env_remove("PYTHONPATH")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    hide_console_window(&mut cmd);
    let mut child = cmd.spawn()?;
    info!(
        "Launching Jupyter with uv at {}: {:?}",
        uv_bin,
        Command::new(&uv_bin).args(&args)
    );

    let pid = child.id();
    println!("✅ Jupyter Lab started (PID: {})", pid);
    println!("   Waiting for server to start...");

    // Take stderr to read in background
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| anyhow!("Failed to capture stderr"))?;

    // Spawn background thread to parse stderr for token
    let (tx, mut rx) = tokio::sync::mpsc::channel::<JupyterRuntimeInfo>(1);
    std::thread::spawn(move || {
        let reader = BufReader::new(stderr);
        for line in reader.lines().map_while(|r| r.ok()) {
            eprintln!("{}", line); // Still print to stderr
            if let Some((port, url, token)) = parse_jupyter_url_from_line(&line) {
                let _ = tx.blocking_send(JupyterRuntimeInfo {
                    port: Some(port),
                    url: Some(url),
                    token: Some(token),
                });
                break;
            }
        }
    });

    let mut runtime_info: Option<JupyterRuntimeInfo> = None;
    for _ in 0..60 {
        // Try to get runtime info from stderr parsing
        if let Ok(info) = rx.try_recv() {
            runtime_info = Some(info);
            break;
        }

        // Also try jpserver files as fallback
        if let Some(info) = find_runtime_info(pid) {
            runtime_info = Some(info);
            break;
        }

        if let Some(status) = child.try_wait()? {
            return Err(anyhow!("Jupyter Lab exited immediately with status: {}", status).into());
        }

        tokio::time::sleep(Duration::from_millis(500)).await;
    }

    if runtime_info.is_none() {
        if let Some(port) = chosen_port {
            println!(
                "⚠️  Jupyter runtime info not detected; falling back to chosen port {}",
                port
            );
            runtime_info = Some(JupyterRuntimeInfo {
                port: Some(port),
                url: Some(format!("http://localhost:{}/lab", port)),
                token: None,
            });
        }
    }

    let mut runtime_info =
        runtime_info.ok_or_else(|| anyhow!("Timed out waiting for Jupyter runtime information"))?;

    // If Jupyter did not report a port, fall back to our chosen one
    if runtime_info.port.is_none() {
        runtime_info.port = chosen_port;
    }

    if let Some(port) = runtime_info.port {
        wait_for_server_ready(port).await?;
        println!("✅ Jupyter server is ready on port {}", port);
    }

    // Ensure the stored URL is usable: if token is empty/None, drop it; otherwise append if missing
    let token_opt = runtime_info
        .token
        .as_ref()
        .and_then(|t| if t.is_empty() { None } else { Some(t) });

    let url_with_token = match (runtime_info.url.as_ref(), token_opt) {
        (Some(url), Some(token)) => {
            if url.contains("token=") {
                Some(url.clone())
            } else {
                let mut new_url = url.clone();
                if url.contains('?') {
                    let _ = write!(new_url, "&token={}", token);
                } else {
                    let _ = write!(new_url, "?token={}", token);
                }
                Some(new_url)
            }
        }
        (Some(url), None) => Some(url.clone()),
        _ => runtime_info.url.clone(),
    };

    let store_token = runtime_info.token.clone().filter(|t| !t.is_empty());

    db.update_jupyter_session(
        &module_dir,
        runtime_info.port,
        Some(pid as i32),
        url_with_token.as_deref(),
        store_token.as_deref(),
    )?;

    // Determine the final URL to use
    let final_url = if let Some(url) = url_with_token.as_ref().or(runtime_info.url.as_ref()) {
        println!("   Access at: {}", url);
        Some(url.clone())
    } else if let Some(port) = runtime_info.port {
        let url = format!("http://localhost:{}", port);
        println!("   Access at: {}", url);
        Some(url)
    } else {
        println!("   Access at: <unknown>");
        None
    };

    // Open browser automatically (unless JUPYTER_SKIP_BROWSER is set for test mode)
    let skip_browser = env::var("JUPYTER_SKIP_BROWSER")
        .map(|v| !v.is_empty() && v != "0")
        .unwrap_or(false);

    if let Some(url) = final_url {
        if skip_browser {
            println!("🌐 Browser open skipped (JUPYTER_SKIP_BROWSER set)");
        } else {
            println!("🌐 Opening browser...");
            #[cfg(target_os = "macos")]
            let _ = Command::new("open").arg(&url).spawn();
            #[cfg(target_os = "linux")]
            let _ = Command::new("xdg-open").arg(&url).spawn();
            #[cfg(target_os = "windows")]
            {
                let mut cmd = Command::new("cmd");
                cmd.args(["/C", "start", "", &url]);
                hide_console_window(&mut cmd);
                let _ = cmd.spawn();
            }
        }
    }

    println!("   Press Ctrl+C in the terminal running Jupyter to stop");
    println!("\n💡 Tip: Jupyter Lab is running in the background");

    Ok(())
}

pub async fn stop(module_path: &str) -> Result<()> {
    // Check if module_path is a number (list index)
    let module_dir = if let Ok(index) = module_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Module index must be >= 1").into());
        }

        let modules = get_modules_with_venvs()?;

        if modules.is_empty() {
            return Err(anyhow!("No modules with virtualenvs found. Run 'bv jupyter list' to see available modules.").into());
        }

        if index > modules.len() {
            return Err(anyhow!(
                "Module index {} out of range. Only {} module(s) available.",
                index,
                modules.len()
            )
            .into());
        }

        modules[index - 1].clone()
    } else {
        PathBuf::from(module_path)
    };

    let venv_path = module_dir.join(".venv");

    if !venv_path.exists() {
        println!(
            "⚠️  Virtualenv not found for {}. Nothing to stop.",
            module_dir.display()
        );
    }

    // Get PID from database and kill the process directly
    let db = BioVaultDb::new()?;
    let canonical_path = module_dir
        .canonicalize()
        .unwrap_or_else(|_| module_dir.clone());
    let env_info = db.get_dev_env(canonical_path.to_str().unwrap_or(""))?;

    let mut stopped = false;
    let mut extra_pids: Vec<i32> = Vec::new();
    if let Some(env) = &env_info {
        if let Some(pid) = env.jupyter_pid {
            println!("🛑 Stopping Jupyter server (PID: {})...", pid);
            terminate_pid(pid)?;
            let _ = wait_for_pid_exit(pid, Duration::from_secs(5)).await;
            println!("✅ Jupyter server stopped");
            stopped = true;
        }
    }

    // Also stop the actual Jupyter ServerApp PID from jpserver runtime files (Windows frequently locks jupyter.exe).
    let runtime_dir = module_dir.join(".jupyter-runtime");
    extra_pids.extend(find_module_runtime_pids(&runtime_dir));
    for pid in extra_pids {
        if is_pid_alive(pid) {
            println!("🛑 Stopping Jupyter runtime PID: {}...", pid);
            let _ = terminate_pid(pid);
            let _ = wait_for_pid_exit(pid, Duration::from_secs(5)).await;
            stopped = true;
        }
    }

    if !stopped {
        println!("⚠️  No Jupyter process found to stop");
    }

    // Clear session info from database
    db.update_jupyter_session(&module_dir, None, None, None, None)?;

    Ok(())
}

pub async fn reset(module_path: &str, python_version: &str) -> Result<()> {
    // Check if module_path is a number (list index)
    let module_dir = if let Ok(index) = module_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Module index must be >= 1").into());
        }

        let modules = get_modules_with_venvs()?;

        if modules.is_empty() {
            return Err(anyhow!("No modules with virtualenvs found. Run 'bv jupyter list' to see available modules.").into());
        }

        if index > modules.len() {
            return Err(anyhow!(
                "Module index {} out of range. Only {} module(s) available.",
                index,
                modules.len()
            )
            .into());
        }

        modules[index - 1].clone()
    } else {
        PathBuf::from(module_path)
    };

    let venv_path = module_dir.join(".venv");

    // Stop Jupyter if running
    if let Some(path_str) = module_dir.to_str() {
        stop(path_str).await?;
    }

    // Remove old venv
    if venv_path.exists() {
        println!("🗑️  Removing old virtualenv...");
        remove_dir_all_with_retry(&venv_path, Duration::from_secs(10)).await?;
        println!("✅ Old virtualenv removed");
    }

    // Delete from database
    let db = BioVaultDb::new()?;
    if let Ok(canonical_path) = module_dir.canonicalize() {
        let _ = db.delete_dev_env(canonical_path.to_str().unwrap());
    }

    // Create fresh venv without launching Jupyter
    println!("🔄 Creating fresh virtualenv...");
    let uv_bin = resolve_uv_path()?;
    ensure_virtualenv(&module_dir, python_version, &uv_bin)?;

    db.register_dev_env(&module_dir, python_version, "jupyter", true)?;
    db.update_jupyter_session(&module_dir, None, None, None, None)?;

    println!("✅ Virtualenv rebuilt. Jupyter server is stopped.");
    Ok(())
}

pub async fn status() -> Result<()> {
    println!("📊 Checking Jupyter Lab status...");

    // Try to find running Jupyter processes
    let output = if cfg!(windows) {
        let mut cmd = Command::new("tasklist");
        cmd.args(["/FI", "IMAGENAME eq jupyter.exe"]);
        hide_console_window(&mut cmd);
        match cmd.output() {
            Ok(out) => out,
            Err(err) => {
                if err.kind() == io::ErrorKind::NotFound {
                    println!("⚪ 'tasklist' not available; skipping process check");
                    return Ok(());
                }
                return Err(err.into());
            }
        }
    } else {
        match Command::new("pgrep").args(["-f", "jupyter-lab"]).output() {
            Ok(out) => out,
            Err(err) => {
                if err.kind() == io::ErrorKind::NotFound {
                    println!("⚪ 'pgrep' not available; skipping process check");
                    return Ok(());
                }
                return Err(err.into());
            }
        }
    };

    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        if stdout.trim().is_empty() {
            println!("⚪ No running Jupyter Lab sessions found");
        } else {
            println!("🟢 Running Jupyter Lab sessions:");
            println!("{}", stdout);
        }
    } else {
        println!("⚪ No running Jupyter Lab sessions found");
    }

    Ok(())
}

fn get_modules_with_venvs() -> Result<Vec<PathBuf>> {
    let db = BioVaultDb::new()?;
    let envs = db.list_dev_envs()?;

    let modules: Vec<PathBuf> = envs
        .iter()
        .filter_map(|env| {
            let path = PathBuf::from(&env.module_path);
            if path.exists() {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    Ok(modules)
}

pub async fn list() -> Result<()> {
    println!("📁 Modules with Jupyter virtualenvs:");

    let modules = get_modules_with_venvs()?;

    if modules.is_empty() {
        println!("   No modules with virtualenvs found");
        println!("\n💡 Tip: Run 'bv jupyter start <module-path>' to create one");
    } else {
        for (i, module) in modules.iter().enumerate() {
            println!("   {}. {}", i + 1, module.display());
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    struct TestEnv {
        _tmp: TempDir,
    }

    impl TestEnv {
        fn new() -> Self {
            let tmp = TempDir::new().unwrap();
            crate::config::set_test_biovault_home(tmp.path());
            Self { _tmp: tmp }
        }
    }

    impl Drop for TestEnv {
        fn drop(&mut self) {
            crate::config::clear_test_biovault_home();
        }
    }

    #[tokio::test]
    async fn test_status_does_not_fail() {
        let _env = TestEnv::new();
        let result = status().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_list_does_not_fail() {
        let _env = TestEnv::new();
        let result = list().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_start_with_nonexistent_dir() {
        let _env = TestEnv::new();
        let result = start("/nonexistent/path", "3.12").await;
        assert!(result.is_err());
    }

    #[tokio::test]
    #[ignore = "requires UV and creates actual virtualenv"]
    async fn test_start_and_reset() {
        let _env = TestEnv::new();
        let tmp = TempDir::new().unwrap();
        let module_path = tmp.path().to_str().unwrap();

        let result = start(module_path, "3.12").await;
        assert!(result.is_ok());

        let venv_path = tmp.path().join(".venv");
        assert!(venv_path.exists());

        let result = reset(module_path, "3.12").await;
        assert!(result.is_ok());
    }
}
