use crate::data::BioVaultDb;
use crate::error::Result;
use anyhow::anyhow;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;
use std::env;
use std::fs;
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

    println!(
        "📦 Installing/Updating packages via: uv pip install -U --python .venv jupyterlab bioscript"
    );

    let status = Command::new("uv")
        .args([
            "pip",
            "install",
            "-U",
            "--python",
            ".venv",
            "jupyterlab",
            "bioscript",
        ])
        .current_dir(project_dir)
        .status()?;

    if !status.success() {
        return Err(
            anyhow!("Failed to install required Python packages (jupyterlab/bioscript)").into(),
        );
    }

    println!("✅ Virtualenv ready with jupyterlab and bioscript");
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

                let (info, info_pid) = match parse_runtime_file(&path) {
                    Some(result) => result,
                    None => continue,
                };

                if let Some(info_pid) = info_pid {
                    if info_pid as u32 == pid {
                        return Some(info);
                    }
                }

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
    println!("🚀 Launching Jupyter Lab with: uv run --python .venv jupyter lab");

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

    let mut child = Command::new("uv")
        .args(["run", "--python", ".venv", "jupyter", "lab", "--no-browser"])
        .current_dir(&project_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    let pid = child.id();
    println!("✅ Jupyter Lab started (PID: {})", pid);
    println!("   Waiting for server to start...");

    let mut runtime_info: Option<JupyterRuntimeInfo> = None;
    for _ in 0..60 {
        if let Some(info) = find_runtime_info(pid) {
            runtime_info = Some(info);
            break;
        }

        if let Some(status) = child.try_wait()? {
            return Err(anyhow!("Jupyter Lab exited immediately with status: {}", status).into());
        }

        tokio::time::sleep(Duration::from_millis(500)).await;
    }

    let runtime_info =
        runtime_info.ok_or_else(|| anyhow!("Timed out waiting for Jupyter runtime information"))?;

    if let Some(port) = runtime_info.port {
        wait_for_server_ready(port).await?;
    }

    db.update_jupyter_session(
        &project_dir,
        runtime_info.port,
        Some(pid as i32),
        runtime_info.url.as_deref(),
        runtime_info.token.as_deref(),
    )?;

    if let Some(url) = runtime_info.url.as_ref() {
        println!("   Access at: {}", url);
    } else if let Some(port) = runtime_info.port {
        println!("   Access at: http://localhost:{}", port);
    } else {
        println!("   Access at: <unknown>");
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
        Command::new("tasklist")
            .args(["/FI", "IMAGENAME eq jupyter.exe"])
            .output()?
    } else {
        Command::new("pgrep").args(["-f", "jupyter-lab"]).output()?
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
