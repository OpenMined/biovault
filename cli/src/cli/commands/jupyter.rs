use crate::data::BioVaultDb;
use crate::error::Result;
use anyhow::anyhow;
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

fn ensure_virtualenv(project_dir: &Path, python_version: &str) -> Result<()> {
    let venv_path = project_dir.join(".venv");

    if !venv_path.exists() {
        println!("📦 Creating virtualenv with Python {}...", python_version);

        let status = Command::new("uv")
            .args(["venv", "--python", python_version, ".venv"])
            .current_dir(project_dir)
            .status()?;

        if !status.success() {
            return Err(anyhow!(
                "Failed to create virtualenv. Try: bv python install {}",
                python_version
            )
            .into());
        }
    } else {
        println!("✅ Using existing virtualenv");
    }

    // Version of biovault-beaver from PyPI - auto-detected from submodule at compile time
    // Falls back to hardcoded version if env var not set (e.g., when building biovault CLI standalone)
    let beaver_version = option_env!("BEAVER_VERSION").unwrap_or("0.1.26");

    println!(
        "📦 Installing/Updating packages via: uv pip install -U jupyterlab cleon biovault-beaver[syftbox]=={}",
        beaver_version
    );

    // Install base packages from PyPI including pinned biovault-beaver with syftbox-sdk
    let beaver_pkg = format!("biovault-beaver[syftbox]=={}", beaver_version);
    let status = Command::new("uv")
        .args(["pip", "install", "-U", "jupyterlab", "cleon", &beaver_pkg])
        .current_dir(project_dir)
        .status()?;

    if !status.success() {
        return Err(anyhow!(
            "Failed to install required Python packages (jupyterlab/cleon/biovault-beaver)"
        )
        .into());
    }

    // Path structure: project_dir is BIOVAULT_HOME/sessions/<session_id>
    // biovault root is at BIOVAULT_HOME/../.. (e.g., workspace3/biovault)
    // So from project_dir: ../../../..
    let biovault_root = project_dir.join("..").join("..").join("..").join("..");

    // DEV MODE: If local source exists, install editable versions on top of PyPI packages
    // This allows developers to test local changes while maintaining compatibility
    let syftbox_path = biovault_root.join("syftbox-sdk").join("python");
    let beaver_path = biovault_root.join("biovault-beaver").join("python");

    if syftbox_path.exists() || beaver_path.exists() {
        println!("🔧 DEV MODE: Local source detected, installing editable packages...");

        // Install syftbox-sdk first (beaver depends on it)
        if syftbox_path.exists() {
            println!("📦 Installing syftbox-sdk from local editable path...");
            let syftbox_canonical = syftbox_path.canonicalize().unwrap_or(syftbox_path.clone());
            let status = Command::new("uv")
                .args([
                    "pip",
                    "install",
                    "-e",
                    syftbox_canonical.to_str().unwrap_or("."),
                ])
                .current_dir(project_dir)
                .status()?;

            if status.success() {
                println!(
                    "✅ syftbox-sdk installed from: {}",
                    syftbox_canonical.display()
                );
            } else {
                println!("⚠️ Failed to install syftbox-sdk from local path");
            }
        }

        // Install beaver from local editable path (overwrites PyPI version)
        if beaver_path.exists() {
            println!("🦫 Installing beaver from local editable path (overwriting PyPI version)...");
            let beaver_canonical = beaver_path.canonicalize().unwrap_or(beaver_path);
            let status = Command::new("uv")
                .args([
                    "pip",
                    "install",
                    "-e",
                    beaver_canonical.to_str().unwrap_or("."),
                ])
                .current_dir(project_dir)
                .status()?;

            if status.success() {
                println!("✅ beaver installed from: {}", beaver_canonical.display());
            } else {
                println!("⚠️ Failed to install beaver from local path");
            }
        }

        println!("✅ Virtualenv ready with jupyterlab, cleon, and DEV beaver/syftbox-sdk");
    } else {
        println!(
            "✅ Virtualenv ready with jupyterlab, cleon, biovault-beaver=={}, and syftbox-sdk",
            beaver_version
        );
    }

    Ok(())
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

fn find_runtime_info(pid: u32) -> Option<JupyterRuntimeInfo> {
    let mut latest: Option<(SystemTime, JupyterRuntimeInfo)> = None;

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

                if let Some(info_pid) = info_pid_opt {
                    if info_pid as u32 != pid {
                        // Skip runtime files for other processes to avoid cross-app collisions
                        continue;
                    }
                }

                // Instead of matching exact PID (which fails with uv run wrapper),
                // find the most recent jpserver file within the last 2 minutes
                if let Ok(metadata) = entry.metadata() {
                    if let Ok(modified) = metadata.modified() {
                        if let Ok(elapsed) = modified.elapsed() {
                            if elapsed.as_secs() > 120 {
                                continue;
                            }
                        }

                        match &latest {
                            Some((best_time, _)) if modified <= *best_time => {}
                            _ => latest = Some((modified, info.clone())),
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

pub async fn start(project_path: &str, python_version: &str) -> Result<()> {
    // Check if project_path is a number (list index)
    let project_dir = if let Ok(index) = project_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Project index must be >= 1").into());
        }

        let projects = get_projects_with_venvs()?;

        if projects.is_empty() {
            return Err(anyhow!("No projects with virtualenvs found. Run 'bv jupyter list' to see available projects.").into());
        }

        if index > projects.len() {
            return Err(anyhow!(
                "Project index {} out of range. Only {} project(s) available.",
                index,
                projects.len()
            )
            .into());
        }

        projects[index - 1].clone()
    } else {
        PathBuf::from(project_path)
    };

    if !project_dir.exists() {
        return Err(anyhow!(
            "Project directory does not exist: {}",
            project_dir.display()
        )
        .into());
    }

    let venv_path = project_dir.join(".venv");

    ensure_virtualenv(&project_dir, python_version)?;

    // Register/update in database
    let db = BioVaultDb::new()?;
    db.register_dev_env(&project_dir, python_version, "jupyter", true)?;

    // Check if Jupyter is already running for this project
    let canonical_path = project_dir.canonicalize()?;
    if let Some(env) = db.get_dev_env(canonical_path.to_str().unwrap())? {
        if let Some(pid) = env.jupyter_pid {
            // Check if process is still alive
            let is_alive = std::process::Command::new("kill")
                .args(["-0", &pid.to_string()])
                .status()
                .map(|s| s.success())
                .unwrap_or(false);

            if is_alive {
                println!("✅ Jupyter Lab already running (PID: {})", pid);
                if let Some(url) = env.jupyter_url.as_ref() {
                    println!("   Access at: {}", url);
                } else if let Some(port) = env.jupyter_port {
                    println!("   Access at: http://localhost:{}", port);
                } else {
                    println!("   Access at: <unknown>");
                }
                println!("   Use 'bv jupyter stop' with project path or index to stop it first");
                return Ok(());
            } else {
                // Clear stale session info
                db.update_jupyter_session(&project_dir, None, None, None, None)?;
            }
        }
    }

    // Launch Jupyter Lab
    let chosen_port = find_available_port();

    // Isolate Jupyter runtime files per session to avoid cross-instance detection
    let runtime_dir = project_dir.join(".jupyter-runtime");
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
            project_path
        )
        .into());
    }

    use std::process::Stdio;

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

    let mut child = Command::new("uv")
        .args(&args)
        .current_dir(&project_dir)
        .env("JUPYTER_RUNTIME_DIR", &runtime_dir)
        .env("XDG_RUNTIME_DIR", &runtime_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

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
        &project_dir,
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

    // Open browser automatically
    if let Some(url) = final_url {
        println!("🌐 Opening browser...");
        #[cfg(target_os = "macos")]
        let _ = Command::new("open").arg(&url).spawn();
        #[cfg(target_os = "linux")]
        let _ = Command::new("xdg-open").arg(&url).spawn();
        #[cfg(target_os = "windows")]
        let _ = Command::new("cmd").args(["/C", "start", &url]).spawn();
    }

    println!("   Press Ctrl+C in the terminal running Jupyter to stop");
    println!("\n💡 Tip: Jupyter Lab is running in the background");

    Ok(())
}

pub async fn stop(project_path: &str) -> Result<()> {
    // Check if project_path is a number (list index)
    let project_dir = if let Ok(index) = project_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Project index must be >= 1").into());
        }

        let projects = get_projects_with_venvs()?;

        if projects.is_empty() {
            return Err(anyhow!("No projects with virtualenvs found. Run 'bv jupyter list' to see available projects.").into());
        }

        if index > projects.len() {
            return Err(anyhow!(
                "Project index {} out of range. Only {} project(s) available.",
                index,
                projects.len()
            )
            .into());
        }

        projects[index - 1].clone()
    } else {
        PathBuf::from(project_path)
    };

    let venv_path = project_dir.join(".venv");

    if !venv_path.exists() {
        println!(
            "⚠️  Virtualenv not found for {}. Nothing to stop.",
            project_dir.display()
        );
        return Ok(());
    }

    println!("🛑 Stopping Jupyter Lab with: uv run --python .venv jupyter lab stop...");

    let status = Command::new("uv")
        .args(["run", "--python", ".venv", "jupyter", "lab", "stop"])
        .current_dir(&project_dir)
        .status()?;

    if status.success() {
        println!("✅ Jupyter Lab stopped");
    } else {
        println!("⚠️  Could not stop Jupyter Lab (may not be running)");
    }

    // Clear session info from database
    let db = BioVaultDb::new()?;
    db.update_jupyter_session(&project_dir, None, None, None, None)?;

    Ok(())
}

pub async fn reset(project_path: &str, python_version: &str) -> Result<()> {
    // Check if project_path is a number (list index)
    let project_dir = if let Ok(index) = project_path.parse::<usize>() {
        if index == 0 {
            return Err(anyhow!("Project index must be >= 1").into());
        }

        let projects = get_projects_with_venvs()?;

        if projects.is_empty() {
            return Err(anyhow!("No projects with virtualenvs found. Run 'bv jupyter list' to see available projects.").into());
        }

        if index > projects.len() {
            return Err(anyhow!(
                "Project index {} out of range. Only {} project(s) available.",
                index,
                projects.len()
            )
            .into());
        }

        projects[index - 1].clone()
    } else {
        PathBuf::from(project_path)
    };

    let venv_path = project_dir.join(".venv");

    // Stop Jupyter if running
    if let Some(path_str) = project_dir.to_str() {
        let _ = stop(path_str).await;
    }

    // Remove old venv
    if venv_path.exists() {
        println!("🗑️  Removing old virtualenv...");
        std::fs::remove_dir_all(&venv_path)?;
        println!("✅ Old virtualenv removed");
    }

    // Delete from database
    let db = BioVaultDb::new()?;
    if let Ok(canonical_path) = project_dir.canonicalize() {
        let _ = db.delete_dev_env(canonical_path.to_str().unwrap());
    }

    // Create fresh venv without launching Jupyter
    println!("🔄 Creating fresh virtualenv...");
    ensure_virtualenv(&project_dir, python_version)?;

    db.register_dev_env(&project_dir, python_version, "jupyter", true)?;
    db.update_jupyter_session(&project_dir, None, None, None, None)?;

    println!("✅ Virtualenv rebuilt. Jupyter server is stopped.");
    Ok(())
}

pub async fn status() -> Result<()> {
    println!("📊 Checking Jupyter Lab status...");

    // Try to find running Jupyter processes
    let output = if cfg!(windows) {
        match Command::new("tasklist")
            .args(["/FI", "IMAGENAME eq jupyter.exe"])
            .output()
        {
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

fn get_projects_with_venvs() -> Result<Vec<PathBuf>> {
    let db = BioVaultDb::new()?;
    let envs = db.list_dev_envs()?;

    let projects: Vec<PathBuf> = envs
        .iter()
        .filter_map(|env| {
            let path = PathBuf::from(&env.project_path);
            if path.exists() {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    Ok(projects)
}

pub async fn list() -> Result<()> {
    println!("📁 Projects with Jupyter virtualenvs:");

    let projects = get_projects_with_venvs()?;

    if projects.is_empty() {
        println!("   No projects with virtualenvs found");
        println!("\n💡 Tip: Run 'bv jupyter start <project-path>' to create one");
    } else {
        for (i, project) in projects.iter().enumerate() {
            println!("   {}. {}", i + 1, project.display());
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
        let project_path = tmp.path().to_str().unwrap();

        let result = start(project_path, "3.12").await;
        assert!(result.is_ok());

        let venv_path = tmp.path().join(".venv");
        assert!(venv_path.exists());

        let result = reset(project_path, "3.12").await;
        assert!(result.is_ok());
    }
}
