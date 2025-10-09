use crate::data::BioVaultDb;
use crate::error::Result;
use anyhow::anyhow;
use std::path::PathBuf;
use std::process::Command;

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

    // Check if venv exists
    if !venv_path.exists() {
        println!("📦 Creating virtualenv with Python {}...", python_version);

        // Create venv with UV-managed Python
        let status = Command::new("uv")
            .args(["venv", "--python", python_version, ".venv"])
            .current_dir(&project_dir)
            .status()?;

        if !status.success() {
            return Err(anyhow!(
                "Failed to create virtualenv. Try: bv python install {}",
                python_version
            )
            .into());
        }

        println!("📦 Installing Jupyter Lab...");

        // Install Jupyter using uv pip into the venv
        let status = Command::new("uv")
            .args(["pip", "install", "-U", "--python", ".venv", "jupyterlab"])
            .current_dir(&project_dir)
            .status()?;

        if !status.success() {
            return Err(anyhow!("Failed to install Jupyter Lab").into());
        }

        println!("✅ Virtualenv created and Jupyter installed");
    } else {
        println!("✅ Using existing virtualenv");
    }

    // Register/update in database
    let db = BioVaultDb::new()?;
    db.register_dev_env(&project_dir, python_version, "jupyter", true)?;

    // Check if Jupyter is already running for this project
    let canonical_path = project_dir.canonicalize()?;
    if let Some(env) = db.get_dev_env(canonical_path.to_str().unwrap())? {
        if let (Some(pid), Some(port)) = (env.jupyter_pid, env.jupyter_port) {
            // Check if process is still alive
            let is_alive = std::process::Command::new("kill")
                .args(["-0", &pid.to_string()])
                .status()
                .map(|s| s.success())
                .unwrap_or(false);

            if is_alive {
                println!("✅ Jupyter Lab already running (PID: {})", pid);
                println!("   Access at: http://localhost:{}", port);
                println!("   Use 'bv jupyter stop' with project path or index to stop it first");
                return Ok(());
            } else {
                // Clear stale session info
                db.update_jupyter_session(&project_dir, None, None)?;
            }
        }
    }

    // Launch Jupyter Lab
    println!("🚀 Launching Jupyter Lab...");

    let jupyter_bin = if cfg!(windows) {
        venv_path.join("Scripts").join("jupyter.exe")
    } else {
        venv_path.join("bin").join("jupyter")
    };

    if !jupyter_bin.exists() {
        return Err(anyhow!(
            "Jupyter not found in virtualenv. Try: bv jupyter reset {}",
            project_path
        )
        .into());
    }

    // Convert to absolute path for spawning
    let jupyter_bin_abs = std::fs::canonicalize(&jupyter_bin)?;

    use std::process::Stdio;

    let mut child = Command::new(jupyter_bin_abs)
        .args(["lab"])
        .current_dir(&project_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    let pid = child.id();
    println!("✅ Jupyter Lab started (PID: {})", pid);
    println!("   Waiting for server to start...");

    // Wait for Jupyter to start and capture port
    tokio::time::sleep(std::time::Duration::from_secs(2)).await;

    // Check if process is still running
    match child.try_wait()? {
        Some(status) => {
            return Err(anyhow!("Jupyter Lab exited immediately with status: {}", status).into());
        }
        None => {
            // Try to find the port from Jupyter runtime files
            let runtime_dir = std::env::var("HOME")
                .map(|h| std::path::PathBuf::from(h).join(".local/share/jupyter/runtime"))
                .unwrap_or_else(|_| std::path::PathBuf::from("/tmp"));

            let mut port = 8888; // Default port

            // Look for runtime files matching our PID
            if runtime_dir.exists() {
                if let Ok(entries) = std::fs::read_dir(&runtime_dir) {
                    for entry in entries.flatten() {
                        let filename = entry.file_name();
                        let filename_str = filename.to_string_lossy();

                        if filename_str.starts_with("jpserver-") && filename_str.ends_with(".json")
                        {
                            if let Ok(content) = std::fs::read_to_string(entry.path()) {
                                // Simple regex-free parsing - look for "port": <number>
                                if let Some(port_start) = content.find("\"port\":") {
                                    let after_colon = &content[port_start + 7..];
                                    if let Some(port_end) = after_colon.find([',', '}']) {
                                        if let Ok(parsed_port) =
                                            after_colon[..port_end].trim().parse::<i32>()
                                        {
                                            // Check if this runtime file was recently modified (within last 5 seconds)
                                            if let Ok(metadata) = entry.metadata() {
                                                if let Ok(modified) = metadata.modified() {
                                                    if let Ok(elapsed) = modified.elapsed() {
                                                        if elapsed.as_secs() < 5 {
                                                            port = parsed_port;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Store session info in database
            db.update_jupyter_session(&project_dir, Some(port), Some(pid as i32))?;

            println!("   Access at: http://localhost:{}", port);
            println!("   Press Ctrl+C in the terminal running Jupyter to stop");
            println!("\n💡 Tip: Jupyter Lab is running in the background");
        }
    }

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

    let jupyter_bin = if cfg!(windows) {
        venv_path.join("Scripts").join("jupyter.exe")
    } else {
        venv_path.join("bin").join("jupyter")
    };

    if !jupyter_bin.exists() {
        println!(
            "⚠️  No Jupyter installation found in {}",
            project_dir.display()
        );
        return Ok(());
    }

    println!("🛑 Stopping Jupyter Lab...");

    // Convert to absolute path for spawning
    let jupyter_bin_abs = std::fs::canonicalize(&jupyter_bin)?;

    let status = Command::new(jupyter_bin_abs)
        .args(["lab", "stop"])
        .current_dir(&project_dir)
        .status()?;

    if status.success() {
        println!("✅ Jupyter Lab stopped");
    } else {
        println!("⚠️  Could not stop Jupyter Lab (may not be running)");
    }

    // Clear session info from database
    let db = BioVaultDb::new()?;
    db.update_jupyter_session(&project_dir, None, None)?;

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

    // Create fresh venv
    println!("🔄 Creating fresh virtualenv...");
    if let Some(path_str) = project_dir.to_str() {
        start(path_str, python_version).await
    } else {
        Err(anyhow!("Invalid project path").into())
    }
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

    #[tokio::test]
    async fn test_status_does_not_fail() {
        let result = status().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_list_does_not_fail() {
        let result = list().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_start_with_nonexistent_dir() {
        let result = start("/nonexistent/path", "3.12").await;
        assert!(result.is_err());
    }

    #[tokio::test]
    #[ignore = "requires UV and creates actual virtualenv"]
    async fn test_start_and_reset() {
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
