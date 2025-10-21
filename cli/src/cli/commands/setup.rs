use super::check::DependencyConfig;
use crate::Result;
use anyhow::anyhow;
use std::collections::HashSet;
use std::env;
use std::io::{self, Write};
use std::process::{Command, ExitStatus, Stdio};

#[cfg(target_os = "macos")]
use std::fs;

struct InstallCommandOutput {
    status: ExitStatus,
    stderr: Vec<u8>,
}

impl InstallCommandOutput {
    fn from_output(output: std::process::Output) -> Self {
        Self {
            status: output.status,
            stderr: output.stderr,
        }
    }
}

fn skip_install_commands() -> bool {
    env::var("BIOVAULT_SKIP_INSTALLS")
        .map(|v| v != "0")
        .unwrap_or(false)
}

#[derive(Debug)]
enum SystemType {
    GoogleColab,
    MacOs,
    Ubuntu,
    ArchLinux,
    Windows,
    Unknown,
}

pub async fn execute(dependencies: Vec<String>, force: bool) -> Result<()> {
    eprintln!(
        "🔧 BioVault Setup: execute() called with dependencies: {:?}, force: {}",
        dependencies, force
    );
    println!("BioVault Environment Setup");
    println!("==========================\n");

    // Validate dependency names
    let valid_deps = ["java", "docker", "nextflow", "syftbox", "uv"];
    if !dependencies.is_empty() {
        for dep in &dependencies {
            if !valid_deps.contains(&dep.as_str()) {
                return Err(anyhow!(
                    "Unknown dependency '{}'. Valid options are: {}",
                    dep,
                    valid_deps.join(", ")
                )
                .into());
            }
        }
        println!(
            "Installing specific dependencies: {}\n",
            dependencies.join(", ")
        );
    } else {
        eprintln!("🔧 No specific dependencies requested - will install all missing");
    }

    if force {
        println!("🔄 Force mode enabled - reinstalling even if already present\n");
    }

    let system_type = detect_system();
    eprintln!("🔧 Detected system type: {:?}", system_type);

    match system_type {
        SystemType::GoogleColab => {
            println!("✓ Detected Google Colab environment");
            eprintln!("🔧 Calling setup_google_colab()");
            setup_google_colab(dependencies, force).await?;
            eprintln!("✅ setup_google_colab() completed");
        }
        SystemType::MacOs => {
            println!("✓ Detected macOS environment");
            eprintln!("🔧 Calling setup_macos()");
            setup_macos(dependencies, force).await?;
            eprintln!("✅ setup_macos() completed");
        }
        SystemType::Ubuntu => {
            println!("✓ Detected Ubuntu/Debian environment");
            eprintln!("🔧 Calling setup_ubuntu()");
            setup_ubuntu(dependencies, force).await?;
            eprintln!("✅ setup_ubuntu() completed");
        }
        SystemType::ArchLinux => {
            println!("✓ Detected Arch Linux environment");
            eprintln!("🔧 Calling setup_arch()");
            setup_arch(dependencies, force).await?;
            eprintln!("✅ setup_arch() completed");
        }
        SystemType::Windows => {
            println!("✓ Detected Windows environment");
            eprintln!("🔧 Calling setup_windows()");
            setup_windows(dependencies, force).await?;
            eprintln!("✅ setup_windows() completed");
        }
        SystemType::Unknown => {
            eprintln!("⚠️  System type is Unknown - cannot proceed with installation");
            println!("ℹ️  System type not detected or not supported for automated setup");
            println!("   This command currently supports:");
            println!("   - Google Colab");
            println!("   - macOS (Homebrew)");
            println!("   - Ubuntu/Debian (apt)");
            println!("   - Arch Linux (pacman)");
            println!("   - Windows (WinGet)");
            println!("\n   For manual setup, please run: bv check");
        }
    }

    eprintln!("✅ execute() completed successfully");
    Ok(())
}

/// Install a single dependency programmatically (for use by desktop app).
/// Returns the path to the installed binary if it can be detected.
pub async fn install_single_dependency(name: &str) -> Result<Option<String>> {
    eprintln!("🔧 install_single_dependency('{}') called", name);

    // Install the dependency
    execute(vec![name.to_string()], false).await?;

    // Try to detect the path using 'which' command
    let path = Command::new("which")
        .arg(name)
        .output()
        .ok()
        .and_then(|output| {
            if output.status.success() {
                String::from_utf8(output.stdout)
                    .ok()
                    .map(|s| s.trim().to_string())
            } else {
                None
            }
        });

    eprintln!(
        "✅ install_single_dependency('{}') completed, path: {:?}",
        name, path
    );
    Ok(path)
}

/// Install multiple dependencies programmatically (for use by desktop app).
pub async fn install_dependencies(names: &[String]) -> Result<()> {
    eprintln!("🔧 install_dependencies({:?}) called", names);
    let result = execute(names.to_vec(), false).await;
    if result.is_ok() {
        eprintln!(
            "✅ install_dependencies({:?}) completed successfully",
            names
        );
    } else {
        eprintln!("❌ install_dependencies({:?}) failed: {:?}", names, result);
    }
    result
}

fn detect_system() -> SystemType {
    let target_os = std::env::consts::OS;

    // Short-circuit for Windows before inspecting other environment hints
    if target_os == "windows" {
        return SystemType::Windows;
    }

    // Check for Google Colab environment variables (Colab runs on Linux)
    if is_google_colab() {
        return SystemType::GoogleColab;
    }

    // Detect macOS
    if target_os == "macos" {
        return SystemType::MacOs;
    }

    // Detect Linux distributions
    if target_os == "linux" {
        // Check for apt (Ubuntu/Debian)
        let has_apt = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v apt-get >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if has_apt {
            return SystemType::Ubuntu;
        }

        // Check for pacman (Arch Linux)
        let has_pacman = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v pacman >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if has_pacman {
            return SystemType::ArchLinux;
        }
    }

    SystemType::Unknown
}

fn is_google_colab() -> bool {
    // Check for COLAB_RELEASE_TAG which is specific to Colab
    if env::var("COLAB_RELEASE_TAG").is_ok() {
        return true;
    }

    // Fallback: check for any COLAB_ prefixed environment variable
    for (key, _) in env::vars() {
        if key.starts_with("COLAB_") {
            return true;
        }
    }

    false
}

async fn setup_google_colab(dependencies: Vec<String>, _force: bool) -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Google Colab setup commands");
        return Ok(());
    }

    println!("\nSetting up Google Colab environment...\n");

    // Load the deps.yaml file to get environment-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Filter and reorder dependencies
    let deps_to_install: Vec<_> = if dependencies.is_empty() {
        // Install all dependencies, but put docker last
        let mut all_deps: Vec<_> = config.dependencies.iter().collect();
        if let Some(docker_idx) = all_deps.iter().position(|d| d.name == "docker") {
            let docker = all_deps.remove(docker_idx);
            all_deps.push(docker);
        }
        all_deps
    } else {
        // Only install specified dependencies
        let dep_set: HashSet<_> = dependencies.iter().map(|s| s.as_str()).collect();
        config
            .dependencies
            .iter()
            .filter(|d| dep_set.contains(d.name.as_str()))
            .collect()
    };

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in deps_to_install {
        // Check if this dependency has google_colab environment config
        if let Some(environments) = &dep.environments {
            if let Some(env_config) = environments.get("google_colab") {
                if env_config.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_config
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_config.install_commands {
                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_succeeded = true;

                    for cmd in install_commands {
                        println!("   Running: {}", cmd);

                        // For Colab, we need to run these commands with sh -c
                        let output = Command::new("sh").arg("-c").arg(cmd).output();

                        match output {
                            Ok(output) => {
                                if output.status.success() {
                                    println!("   ✓ Command succeeded");
                                } else {
                                    println!("   ❌ Command failed");
                                    if !output.stderr.is_empty() {
                                        println!(
                                            "   Error: {}",
                                            String::from_utf8_lossy(&output.stderr)
                                        );
                                    }
                                    all_succeeded = false;
                                    break;
                                }
                            }
                            Err(e) => {
                                println!("   ❌ Failed to execute: {}", e);
                                all_succeeded = false;
                                break;
                            }
                        }
                    }

                    // Verify installation if verification command is provided
                    if all_succeeded {
                        if let Some(verify_cmd) = &env_config.verify_command {
                            print!("   Verifying installation... ");
                            let output = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::piped())
                                .stderr(Stdio::piped())
                                .output();

                            if let Ok(output) = output {
                                if output.status.success() {
                                    println!("✓");
                                    success_count += 1;
                                } else {
                                    println!("❌ Verification failed");
                                    fail_count += 1;
                                }
                            } else {
                                println!("❌ Could not verify");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    // Add PATH export instructions for Colab
    println!("📝 Final setup steps for Google Colab:\n");
    println!("   Add these lines to your notebook for persistence:");
    println!("   ```python");
    println!("   import os");
    println!("   os.environ['PATH'] = f\"/usr/local/bin:{{os.environ['PATH']}}\"");
    println!("   ```");
    println!();
    println!("   Or in a shell cell:");
    println!("   ```bash");
    println!("   !export PATH=\"/usr/local/bin:$PATH\"");
    println!("   ```");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
    } else {
        println!("\n✅ Setup completed successfully!");
        println!("   Run 'bv check' to verify all dependencies.");
    }

    Ok(())
}

async fn setup_macos(dependencies: Vec<String>, force: bool) -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping macOS setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up macOS environment...\n");

    // Check if we're in CI mode (non-interactive)
    let is_ci = env::var("CI").is_ok() || env::var("GITHUB_ACTIONS").is_ok();
    let assume_yes = env::var("BIOVAULT_SETUP_ASSUME_YES")
        .map(|value| {
            let normalized = value.to_ascii_lowercase();
            normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "y"
        })
        .unwrap_or(false);

    if assume_yes {
        env::set_var("NONINTERACTIVE", "1");
        env::set_var("HOMEBREW_NO_AUTO_UPDATE", "1");
        env::set_var("HOMEBREW_NO_ENV_HINTS", "1");
    }

    if let Ok(home_dir) = env::var("HOME") {
        env::set_var(
            "HOMEBREW_CASK_OPTS",
            format!("--appdir={}/Applications", home_dir),
        );
    }

    // Check for Homebrew
    let brew_in_path = which::which("brew").is_ok();
    let mut brew_path = None;

    if !brew_in_path {
        // Check common Homebrew locations even if not in PATH
        let common_brew_paths = vec![
            "/opt/homebrew/bin/brew", // Apple Silicon
            "/usr/local/bin/brew",    // Intel Mac
        ];

        for path in &common_brew_paths {
            if std::path::Path::new(path).exists() {
                brew_path = Some(path.to_string());
                println!("📦 Found Homebrew at {} (not in PATH)", path);
                break;
            }
        }
    }

    // Install Homebrew if not found
    if !brew_in_path && brew_path.is_none() {
        println!("📦 Homebrew not found. Would you like to install it? [Y/n]: ");

        if is_ci && !assume_yes {
            println!("   CI mode: Skipping Homebrew installation.");
            println!("   Please ensure Homebrew is pre-installed in CI environment.");
            return Ok(());
        }

        let should_install = if assume_yes {
            println!("   Auto-confirming Homebrew installation (BIOVAULT_SETUP_ASSUME_YES=1).");
            true
        } else {
            io::stdout().flush()?;
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let answer = input.trim().to_lowercase();
            answer.is_empty() || answer == "y" || answer == "yes"
        };

        if should_install {
            println!("Installing Homebrew...");
            let install_cmd = "/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"";
            let status = Command::new("sh").arg("-c").arg(install_cmd).status()?;

            if status.success() {
                println!("✓ Homebrew installed successfully!");
                if std::path::Path::new("/opt/homebrew/bin/brew").exists() {
                    brew_path = Some("/opt/homebrew/bin/brew".to_string());
                } else if std::path::Path::new("/usr/local/bin/brew").exists() {
                    brew_path = Some("/usr/local/bin/brew".to_string());
                }
            } else {
                println!("❌ Homebrew installation failed.");
                println!("Please install manually from: https://brew.sh");
                return Ok(());
            }
        } else {
            println!("Skipping Homebrew installation.");
            println!("Please install Homebrew manually from: https://brew.sh");
            println!("Then re-run: bv setup");
            return Ok(());
        }
    }

    // Use the brew command (either from PATH or specific path)
    let brew_cmd = if brew_in_path {
        "brew".to_string()
    } else if let Some(ref bp) = brew_path {
        bp.clone()
    } else {
        "brew".to_string() // Fallback
    };

    // Load deps.yaml and execute macOS-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Filter and reorder dependencies
    let deps_to_install: Vec<_> = if dependencies.is_empty() {
        // Install all dependencies, but put docker last
        let mut all_deps: Vec<_> = config.dependencies.iter().collect();
        // Find docker and move it to the end
        if let Some(docker_idx) = all_deps.iter().position(|d| d.name == "docker") {
            let docker = all_deps.remove(docker_idx);
            all_deps.push(docker);
        }
        all_deps
    } else {
        // Only install specified dependencies
        let dep_set: HashSet<_> = dependencies.iter().map(|s| s.as_str()).collect();
        config
            .dependencies
            .iter()
            .filter(|d| dep_set.contains(d.name.as_str()))
            .collect()
    };

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    // Use the brew command for installations
    for dep in deps_to_install {
        if let Some(environments) = &dep.environments {
            if let Some(env_config) = environments.get("macos") {
                if env_config.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_config
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on macOS".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_config.install_commands {
                    // Decide if install is necessary (skip check if force mode)
                    let need_install = if force {
                        true
                    } else {
                        let mut should_install = true;

                        // If verify_command is available, try it first
                        if let Some(verify_cmd) = &env_config.verify_command {
                            let verified = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .env("PATH", ensure_user_space_docker_path())
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if verified {
                                should_install = false;
                            }
                        } else {
                            // Fallback to which for simple presence
                            if which::which(&dep.name).is_ok() {
                                should_install = false;
                            }
                        }

                        // For Java, also enforce min_version if specified
                        if dep.name == "java" {
                            if let Some(min_v) = dep.min_version {
                                if let Some(current) = java_major_version() {
                                    if current >= min_v {
                                        println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                        should_install = false;
                                    }
                                }
                            }
                        }

                        should_install
                    };

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_succeeded = true;

                    for cmd in install_commands {
                        // Special handling for Docker - download and install from .dmg directly
                        // Check this BEFORE modifying the command
                        if dep.name == "docker" && cmd.contains("brew install --cask") {
                            println!("   Installing Docker Desktop from official .dmg");
                            match install_docker_desktop_from_dmg() {
                                Ok(_) => {
                                    // Verify installation by checking if Docker.app exists
                                    if std::path::Path::new("/Applications/Docker.app").exists() {
                                        println!("   ✓ Docker Desktop installed successfully");
                                        println!("   ℹ️  Docker Desktop has been installed to /Applications");
                                        println!("   ℹ️  Launch it manually to complete setup (Docker.app may not start in VMs)");
                                    } else {
                                        eprintln!("⚠️  Docker.app not found at /Applications/Docker.app after installation");
                                        println!("   ❌ Docker Desktop installation could not be verified");
                                        all_succeeded = false;
                                        break;
                                    }
                                }
                                Err(e) => {
                                    println!("   ❌ Failed to install Docker Desktop: {}", e);
                                    println!("   💡 Fallback: You can manually download Docker Desktop from:");
                                    println!(
                                        "      https://www.docker.com/products/docker-desktop/"
                                    );
                                    all_succeeded = false;
                                    break;
                                }
                            }
                            continue;
                        }

                        // Replace 'brew' with the actual brew path if needed
                        let mut adjusted_cmd = if !brew_in_path && cmd.starts_with("brew ") {
                            cmd.replace("brew ", &format!("{} ", brew_cmd))
                        } else {
                            cmd.clone()
                        };

                        // Add --force to brew install commands when force mode is enabled
                        if force
                            && adjusted_cmd.starts_with("brew install")
                            && !adjusted_cmd.contains("--force")
                        {
                            // Insert --force before the package name
                            if let Some(pkg_start) = adjusted_cmd.find("brew install") {
                                let prefix = &adjusted_cmd[..pkg_start + "brew install".len()];
                                let suffix = &adjusted_cmd[pkg_start + "brew install".len()..];
                                adjusted_cmd = format!("{} --force{}", prefix, suffix);
                            }
                        }

                        println!("   Running: {}", adjusted_cmd);
                        let output = Command::new("sh")
                            .arg("-c")
                            .arg(&adjusted_cmd)
                            .stdout(Stdio::piped())
                            .stderr(Stdio::piped())
                            .output()
                            .map(InstallCommandOutput::from_output);

                        match output {
                            Ok(output) => {
                                if output.status.success() {
                                    println!("   ✓ Command succeeded");
                                } else {
                                    println!("   ❌ Command failed");
                                    if !output.stderr.is_empty() {
                                        println!(
                                            "   Error: {}",
                                            String::from_utf8_lossy(&output.stderr)
                                        );
                                    }
                                    all_succeeded = false;
                                    break;
                                }
                            }
                            Err(e) => {
                                println!("   ❌ Failed to execute: {}", e);
                                all_succeeded = false;
                                break;
                            }
                        }
                    }

                    if all_succeeded {
                        // Special handling for Java: add brew path before verification
                        if dep.name == "java" {
                            if let Some(java_path) = check_java_in_brew_not_in_path() {
                                // Add the brew Java path to current environment for verification
                                let current_path = env::var("PATH").unwrap_or_default();
                                env::set_var("PATH", format!("{}:{}", java_path, current_path));
                            }
                        }

                        if dep.name == "docker" {
                            if !is_ci {
                                println!("   ⚠️  Skipping automated Docker Desktop privileged installer.");
                                println!(
                                    "      Please open Docker Desktop manually once to finish setup."
                                );
                            } else {
                                println!("   CI mode: Skipping Docker Desktop privileged setup.");
                            }

                            let updated_path = ensure_user_space_docker_path();
                            env::set_var("PATH", &updated_path);
                        }

                        if let Some(verify_cmd) = &env_config.verify_command {
                            print!("   Verifying installation... ");
                            let verify_result = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .env("PATH", ensure_user_space_docker_path())
                                .output();

                            match verify_result {
                                Ok(output) if output.status.success() => {
                                    println!("✓");
                                    success_count += 1;
                                }
                                Ok(output) => {
                                    // Special handling for Docker: check if Docker.app exists
                                    if dep.name == "docker"
                                        && std::path::Path::new("/Applications/Docker.app").exists()
                                    {
                                        println!("✓ (Docker.app installed)");
                                        println!("   ℹ️  Docker Desktop is installed but CLI tools aren't ready yet.");
                                        println!("      Launch Docker Desktop and complete setup to activate CLI tools.");
                                        success_count += 1;
                                    } else {
                                        // Installation commands succeeded but verification failed
                                        // This often happens because PATH hasn't been updated yet
                                        println!("⚠️  Verification failed (may need new shell to pick up PATH)");
                                        if !output.stdout.is_empty() {
                                            eprintln!(
                                                "   stdout: {}",
                                                String::from_utf8_lossy(&output.stdout).trim_end()
                                            );
                                        }
                                        if !output.stderr.is_empty() {
                                            eprintln!(
                                                "   stderr: {}",
                                                String::from_utf8_lossy(&output.stderr).trim_end()
                                            );
                                        }
                                        println!("   ℹ️  Installation commands succeeded - try 'bv check' to verify");
                                        // Count as success since install commands worked
                                        success_count += 1;
                                    }
                                }
                                Err(err) => {
                                    // Special handling for Docker: check if Docker.app exists
                                    if dep.name == "docker"
                                        && std::path::Path::new("/Applications/Docker.app").exists()
                                    {
                                        println!("✓ (Docker.app installed)");
                                        println!("   ℹ️  Docker Desktop is installed but CLI tools aren't ready yet.");
                                        println!("      Launch Docker Desktop and complete setup to activate CLI tools.");
                                        success_count += 1;
                                    } else {
                                        println!(
                                            "⚠️  Could not verify (failed to run '{}'): {}",
                                            verify_cmd, err
                                        );
                                        println!("   ℹ️  Installation commands succeeded - try 'bv check' to verify");
                                        // Count as success since install commands worked
                                        success_count += 1;
                                    }
                                }
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    // After all installations, check if Java needs PATH configuration
    // This handles the case where Java was already installed but not in PATH
    check_and_configure_java_path(is_ci).await?;

    println!("\nNotes:");
    println!("- If this is your first time installing Docker Desktop, open it once to finish setup and grant permissions.");

    // SyftBox info for manual setup later (we installed in setup-only mode)
    println!("\nSyftBox:");
    print_syftbox_instructions();

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

fn print_syftbox_instructions() {
    // Best-effort arch hint for user
    let arch = match std::env::consts::ARCH {
        "aarch64" => "arm64 (Apple Silicon)",
        "x86_64" => "x86_64 (Intel)",
        other => other,
    };
    println!(
        "Get the latest SyftBox for macOS ({}):\n  https://github.com/OpenMined/syftbox/releases/latest",
        arch
    );
    println!("After downloading, ensure the 'syftbox' binary is on your PATH (e.g., move to /usr/local/bin and chmod +x).");
}

#[cfg(target_os = "macos")]
fn install_docker_desktop_from_dmg() -> Result<()> {
    use std::env;

    eprintln!("🔧 install_docker_desktop_from_dmg() starting");

    // Detect architecture
    let arch = env::consts::ARCH;
    let download_url = match arch {
        "aarch64" => "https://desktop.docker.com/mac/main/arm64/Docker.dmg",
        "x86_64" => "https://desktop.docker.com/mac/main/amd64/Docker.dmg",
        _ => {
            eprintln!("❌ Unsupported architecture: {}", arch);
            return Err(anyhow!("Unsupported architecture: {}", arch).into());
        }
    };

    eprintln!("🔧 Downloading Docker Desktop for architecture: {}", arch);
    println!("   Downloading Docker Desktop for {}...", arch);

    // Create temp directory
    let temp_dir = env::temp_dir();
    let dmg_path = temp_dir.join("Docker.dmg");
    eprintln!("🔧 Download path: {:?}", dmg_path);

    // Download the .dmg
    let status = Command::new("curl")
        .arg("-L")
        .arg("-o")
        .arg(&dmg_path)
        .arg(download_url)
        .arg("--progress-bar")
        .status()
        .map_err(|e| {
            eprintln!("❌ Failed to execute curl: {}", e);
            anyhow!("Failed to download Docker Desktop: {}", e)
        })?;

    if !status.success() {
        eprintln!("❌ curl command failed with status: {:?}", status.code());
        return Err(anyhow!(
            "Failed to download Docker Desktop (curl exit code: {:?})",
            status.code()
        )
        .into());
    }

    eprintln!("✅ Download completed");
    println!("   Mounting Docker.dmg...");

    // Mount the .dmg
    let output = Command::new("hdiutil")
        .arg("attach")
        .arg(&dmg_path)
        .arg("-nobrowse")
        .output()
        .map_err(|e| {
            eprintln!("❌ Failed to execute hdiutil: {}", e);
            anyhow!("Failed to mount Docker.dmg: {}", e)
        })?;

    if !output.status.success() {
        eprintln!(
            "❌ hdiutil failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Err(anyhow!(
            "Failed to mount Docker.dmg: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    // Parse the mount point from hdiutil output
    let mount_output = String::from_utf8_lossy(&output.stdout);
    eprintln!("🔧 hdiutil output: {}", mount_output);

    let mount_point = mount_output
        .lines()
        .filter(|line| line.contains("/Volumes/"))
        .next_back()
        .and_then(|line| line.split_whitespace().last())
        .ok_or_else(|| {
            eprintln!("❌ Could not find mount point in hdiutil output");
            anyhow!("Could not find mount point")
        })?;

    eprintln!("✅ Mounted at: {}", mount_point);
    println!("   Copying Docker.app to /Applications...");
    println!("   (macOS will prompt for your password)");

    // Copy Docker.app to /Applications using osascript for GUI authentication
    let docker_app_source = format!("{}/Docker.app", mount_point);
    let docker_app_dest = "/Applications/Docker.app";

    eprintln!("🔧 Source: {}", docker_app_source);
    eprintln!("🔧 Destination: {}", docker_app_dest);

    // Build AppleScript command to remove old installation and copy new one in a single auth prompt
    // Combine both operations into one shell script so user only authenticates once
    let install_script = format!(
        "do shell script \"rm -rf '{}' ; cp -R '{}' '{}'\" with administrator privileges",
        docker_app_dest,
        docker_app_source.replace("'", "'\\''"),
        docker_app_dest
    );

    eprintln!("🔧 Installing Docker.app with administrator privileges...");
    let output = Command::new("osascript")
        .arg("-e")
        .arg(&install_script)
        .output()
        .map_err(|e| {
            eprintln!("❌ Failed to execute osascript: {}", e);
            anyhow!("Failed to install Docker.app: {}", e)
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("❌ osascript failed: {}", stderr);

        // Check if user cancelled
        if stderr.contains("User canceled") || stderr.contains("-128") {
            return Err(anyhow!("Installation cancelled by user").into());
        }

        return Err(anyhow!("Failed to copy Docker.app to /Applications: {}", stderr).into());
    }

    eprintln!("✅ Docker.app copied successfully");

    // Unmount the .dmg
    println!("   Cleaning up...");
    eprintln!("🔧 Unmounting DMG...");
    let _ = Command::new("hdiutil")
        .arg("detach")
        .arg(mount_point)
        .status();

    // Clean up downloaded .dmg
    eprintln!("🔧 Removing temporary DMG file...");
    let _ = fs::remove_file(&dmg_path);

    // Launch Docker Desktop
    println!("   Launching Docker Desktop...");
    eprintln!("🔧 Launching Docker Desktop...");
    let result = Command::new("open").arg(docker_app_dest).spawn();

    match result {
        Ok(_) => eprintln!("✅ Docker Desktop launched"),
        Err(e) => eprintln!("⚠️  Could not launch Docker Desktop: {}", e),
    }

    eprintln!("✅ install_docker_desktop_from_dmg() completed successfully");
    Ok(())
}

#[cfg(not(target_os = "macos"))]
fn install_docker_desktop_from_dmg() -> Result<()> {
    Err(anyhow!("Docker Desktop .dmg installation is only supported on macOS").into())
}

// Minimal java version detection to respect min_version in deps.yaml
fn java_major_version() -> Option<u32> {
    let out = Command::new("java").arg("-version").output().ok()?;
    let text = String::from_utf8_lossy(&out.stderr);
    parse_java_version(&text)
}

fn parse_java_version(output: &str) -> Option<u32> {
    for line in output.lines() {
        if line.contains("version") {
            if let Some(start) = line.find('"') {
                if let Some(end) = line[start + 1..].find('"') {
                    let version_str = &line[start + 1..start + 1 + end];
                    if let Some(stripped) = version_str.strip_prefix("1.") {
                        if let Some(dot_pos) = stripped.find('.') {
                            if let Ok(v) = stripped[..dot_pos].parse::<u32>() {
                                return Some(v);
                            }
                        }
                    } else {
                        let major_part = version_str.split('.').next().unwrap_or(version_str);
                        if let Ok(v) = major_part.parse::<u32>() {
                            return Some(v);
                        }
                    }
                }
            }
        }
    }
    None
}

async fn setup_ubuntu(dependencies: Vec<String>, _force: bool) -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Ubuntu/Debian setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up Ubuntu/Debian environment...\n");

    // Ensure apt-get exists
    let apt_exists = Command::new("sh")
        .arg("-c")
        .arg("command -v apt-get >/dev/null 2>&1")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !apt_exists {
        println!("❌ apt-get not found. This setup targets Ubuntu/Debian-based systems.");
        println!("Please ensure you're on an apt-based distribution.");
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Filter and reorder dependencies
    let deps_to_install: Vec<_> = if dependencies.is_empty() {
        // Install all dependencies, but put docker last
        let mut all_deps: Vec<_> = config.dependencies.iter().collect();
        if let Some(docker_idx) = all_deps.iter().position(|d| d.name == "docker") {
            let docker = all_deps.remove(docker_idx);
            all_deps.push(docker);
        }
        all_deps
    } else {
        // Only install specified dependencies
        let dep_set: HashSet<_> = dependencies.iter().map(|s| s.as_str()).collect();
        config
            .dependencies
            .iter()
            .filter(|d| dep_set.contains(d.name.as_str()))
            .collect()
    };

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in deps_to_install {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("ubuntu") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Ubuntu".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("sh")
                            .arg("-c")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    } else if which::which(&dep.name).is_ok() {
                        need_install = false;
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        // For apt commands on CI, we may need to run with sudo
                        let cmd_to_run = if cmd.starts_with("apt-get") && !cmd.starts_with("sudo") {
                            format!("sudo {}", cmd)
                        } else {
                            cmd.clone()
                        };

                        let status = Command::new("sh").arg("-c").arg(&cmd_to_run).status();
                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!("- For Docker on Ubuntu, you may need to add your user to the docker group: 'sudo usermod -aG docker $USER' and re-login.");
    println!("- If Docker service is not running: 'sudo systemctl start docker'");

    println!("\nSyftBox:");
    println!("The installer has been invoked in setup-only mode if needed.");
    println!("If you want to set up later: syftbox login; syftbox");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

async fn setup_arch(dependencies: Vec<String>, _force: bool) -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Arch Linux setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up Arch Linux environment...\n");

    // Ensure pacman exists
    let pacman_exists = Command::new("sh")
        .arg("-c")
        .arg("command -v pacman >/dev/null 2>&1")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !pacman_exists {
        println!("❌ pacman not found. This setup targets Arch Linux.");
        println!("Please ensure you're on Arch/Manjaro with pacman available.");
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Filter and reorder dependencies
    let deps_to_install: Vec<_> = if dependencies.is_empty() {
        // Install all dependencies, but put docker last
        let mut all_deps: Vec<_> = config.dependencies.iter().collect();
        if let Some(docker_idx) = all_deps.iter().position(|d| d.name == "docker") {
            let docker = all_deps.remove(docker_idx);
            all_deps.push(docker);
        }
        all_deps
    } else {
        // Only install specified dependencies
        let dep_set: HashSet<_> = dependencies.iter().map(|s| s.as_str()).collect();
        config
            .dependencies
            .iter()
            .filter(|d| dep_set.contains(d.name.as_str()))
            .collect()
    };

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in deps_to_install {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("arch") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Arch".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("sh")
                            .arg("-c")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    } else if which::which(&dep.name).is_ok() {
                        need_install = false;
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        let status = Command::new("sh").arg("-c").arg(cmd).status();
                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!("- For Docker on Arch, you may need to enable and start the daemon: 'sudo systemctl enable --now docker' and add your user to the docker group.");

    println!("\nSyftBox:");
    println!("The installer has been invoked in setup-only mode if needed.");
    println!("If you want to set up later: syftbox login; syftbox");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

async fn setup_windows(dependencies: Vec<String>, _force: bool) -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Windows setup commands");
        return Ok(());
    }

    println!("\nSetting up Windows environment...\n");

    // Check for WinGet availability
    let winget_exists = Command::new("winget")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    // Check for Chocolatey availability (fallback)
    let choco_exists = Command::new("choco")
        .arg("-v")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if winget_exists {
        println!("✓ WinGet found");
    } else if choco_exists {
        println!("❌ WinGet not found. Using Chocolatey fallback.");
    } else {
        println!("❌ Neither WinGet nor Chocolatey found.");
        println!("Automated installation is unavailable on this system.");
        println!("\nTo install WinGet:");
        println!("1. Update Windows to the latest version (WinGet comes with modern Windows)");
        println!("2. Or install from Microsoft Store: 'App Installer'");
        println!("3. Or download from: https://github.com/microsoft/winget-cli/releases");
        println!("\nAlternatively, install Chocolatey from https://chocolatey.org/install");
        print_windows_manual_instructions();
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Filter and reorder dependencies
    let deps_to_install: Vec<_> = if dependencies.is_empty() {
        // Install all dependencies, but put docker last
        let mut all_deps: Vec<_> = config.dependencies.iter().collect();
        if let Some(docker_idx) = all_deps.iter().position(|d| d.name == "docker") {
            let docker = all_deps.remove(docker_idx);
            all_deps.push(docker);
        }
        all_deps
    } else {
        // Only install specified dependencies
        let dep_set: HashSet<_> = dependencies.iter().map(|s| s.as_str()).collect();
        config
            .dependencies
            .iter()
            .filter(|d| dep_set.contains(d.name.as_str()))
            .collect()
    };

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in deps_to_install {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("windows") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Windows".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("powershell")
                            .arg("-Command")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        let status = if cmd.starts_with("winget") {
                            if winget_exists {
                                Command::new("winget")
                                    .args(cmd.split_whitespace().skip(1))
                                    .status()
                            } else {
                                // Chocolatey fallback for common packages
                                let mut parts = cmd.split_whitespace();
                                let _ = parts.next(); // winget
                                let _ = parts.next(); // install
                                let pkg = parts.next().unwrap_or("");
                                let choco_pkg = map_winget_pkg_to_choco(pkg);
                                if choco_pkg.is_empty() {
                                    Err(std::io::Error::other("No Chocolatey mapping for package"))
                                } else {
                                    Command::new("choco")
                                        .arg("install")
                                        .arg(choco_pkg)
                                        .arg("-y")
                                        .status()
                                }
                            }
                        } else {
                            Command::new("powershell").arg("-Command").arg(cmd).status()
                        };

                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("powershell")
                                .arg("-Command")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!(
        "- You may need to restart your terminal/PowerShell after installation to update PATH"
    );
    println!("- For Docker on Windows, Docker Desktop is required and may need manual setup");

    print_windows_manual_instructions();

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

fn print_windows_manual_instructions() {
    println!("\nManual Installation Options:");
    println!(
        "Java 17+: Download from https://openjdk.org/ or use 'winget install Microsoft.OpenJDK'"
    );
    println!(
        "Docker: Download Docker Desktop from https://www.docker.com/products/docker-desktop/"
    );
    println!("Nextflow: Use WSL or Docker - no native Windows support");
    println!("SyftBox: Download from https://github.com/OpenMined/syftbox/releases/latest");
    println!("UV: Run 'winget install --id=astral-sh.uv -e' or use PowerShell installer from https://docs.astral.sh/uv/");
}

// Map common WinGet package IDs to Chocolatey package names for fallback
fn map_winget_pkg_to_choco(pkg: &str) -> &'static str {
    match pkg.to_ascii_lowercase().as_str() {
        // Java/OpenJDK
        // WinGet: Microsoft.OpenJDK => Chocolatey: openjdk (generic)
        "microsoft.openjdk" => "openjdk",
        // UV - Python package installer
        // WinGet: astral-sh.uv => Chocolatey doesn't have UV yet, so we return empty
        "astral-sh.uv" => "",
        // Add other mappings here as needed
        _ => "",
    }
}

async fn check_and_configure_java_path(is_ci: bool) -> Result<()> {
    // Check if Java is already in PATH
    if which::which("java").is_ok() {
        return Ok(());
    }

    let assume_yes = env::var("BIOVAULT_SETUP_ASSUME_YES")
        .map(|value| {
            let normalized = value.to_ascii_lowercase();
            normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "y"
        })
        .unwrap_or(false);

    // Check if Java is installed via brew but not in PATH
    let java_brew_path = check_java_in_brew_not_in_path();

    if let Some(brew_path) = java_brew_path {
        println!("\n⚠️  Java is installed via Homebrew but not in your PATH.");
        println!("   Location: {}", brew_path);

        let shell = env::var("SHELL").unwrap_or_else(|_| "/bin/zsh".to_string());
        let shell_config = if shell.contains("zsh") {
            format!(
                "{}/.zshrc",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        } else if shell.contains("bash") {
            format!(
                "{}/.bash_profile",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        } else {
            format!(
                "{}/.profile",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        };

        let export_line = format!("export PATH=\"{}:$PATH\"", brew_path);

        if is_ci || assume_yes {
            if is_ci {
                println!("   CI mode: Automatically configuring PATH...");
            } else {
                println!(
                    "   Auto-confirming Java PATH configuration (BIOVAULT_SETUP_ASSUME_YES=1)."
                );
            }

            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&shell_config)?;
            writeln!(file, "\n# Added by BioVault setup")?;
            writeln!(file, "{}", export_line)?;

            println!("   ✓ Added to {}", shell_config);
            println!("   Note: You'll need to restart your shell or run 'source {}' for changes to take effect.", shell_config);
        } else {
            println!("\n   Would you like to automatically add Java to your PATH? [Y/n]: ");
            io::stdout().flush()?;

            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let answer = input.trim().to_lowercase();

            if answer.is_empty() || answer == "y" || answer == "yes" {
                let mut file = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&shell_config)?;
                writeln!(file, "\n# Added by BioVault setup")?;
                writeln!(file, "{}", export_line)?;

                println!("   ✓ Added to {}", shell_config);
                println!("   Note: You'll need to restart your shell or run 'source {}' for changes to take effect.", shell_config);
            } else {
                println!("   Skipped PATH configuration.");
                println!("   To manually add Java to your PATH, run:");
                println!("     echo '{}' >> {}", export_line, shell_config);
                println!("     source {}", shell_config);
            }
        }
    }

    Ok(())
}
#[cfg(target_os = "macos")]
fn check_java_in_brew_not_in_path() -> Option<String> {
    // Find brew command (in PATH or common locations)
    let brew_cmd = find_brew_command();
    brew_cmd.as_ref()?;
    let brew_cmd = brew_cmd.unwrap();

    // Check if Java/OpenJDK is installed via brew
    let output = Command::new(&brew_cmd)
        .args(["list", "--formula"])
        .output()
        .ok()?;

    let installed_packages = String::from_utf8_lossy(&output.stdout);

    // Look for any OpenJDK version
    let mut found_java_package = None;
    for line in installed_packages.lines() {
        if line.starts_with("openjdk") {
            found_java_package = Some(line.to_string());
            break;
        }
    }

    found_java_package.as_ref()?;

    // Get the actual path where brew installed Java
    let pkg = found_java_package.unwrap();
    let prefix_output = Command::new(&brew_cmd)
        .args(["--prefix", &pkg])
        .output()
        .ok()?;

    if !prefix_output.status.success() {
        return None;
    }

    let brew_prefix = String::from_utf8_lossy(&prefix_output.stdout)
        .trim()
        .to_string();
    let java_bin_path = format!("{}/bin", brew_prefix);

    // Check if this path contains java binary
    if std::path::Path::new(&format!("{}/java", java_bin_path)).exists() {
        Some(java_bin_path)
    } else {
        None
    }
}

#[cfg(not(target_os = "macos"))]
fn check_java_in_brew_not_in_path() -> Option<String> {
    None
}
#[cfg(target_os = "macos")]
fn find_brew_command() -> Option<String> {
    // First check if brew is in PATH
    if which::which("brew").is_ok() {
        return Some("brew".to_string());
    }

    // Check common locations
    let common_brew_paths = vec![
        "/opt/homebrew/bin/brew", // Apple Silicon
        "/usr/local/bin/brew",    // Intel Mac
    ];

    for path in &common_brew_paths {
        if std::path::Path::new(path).exists() {
            return Some(path.to_string());
        }
    }

    None
}

#[cfg(not(target_os = "macos"))]
#[allow(dead_code)]
fn find_brew_command() -> Option<String> {
    None
}

fn ensure_user_space_docker_path() -> String {
    let current_path = env::var("PATH").unwrap_or_default();

    if let Ok(home) = env::var("HOME") {
        let entry = format!("{}/.docker/bin", home);
        if current_path.split(':').any(|segment| segment == entry) {
            current_path
        } else if current_path.is_empty() {
            entry
        } else {
            format!("{}:{}", entry, current_path)
        }
    } else {
        current_path
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;

    struct SkipInstallGuard(Option<String>);

    impl SkipInstallGuard {
        fn new() -> Self {
            let previous = env::var("BIOVAULT_SKIP_INSTALLS").ok();
            env::set_var("BIOVAULT_SKIP_INSTALLS", "1");
            Self(previous)
        }
    }

    impl Drop for SkipInstallGuard {
        fn drop(&mut self) {
            if let Some(ref value) = self.0 {
                env::set_var("BIOVAULT_SKIP_INSTALLS", value);
            } else {
                env::remove_var("BIOVAULT_SKIP_INSTALLS");
            }
        }
    }

    #[test]
    #[serial_test::serial]
    fn skip_install_commands_env_behavior() {
        // Save original state to restore later
        let original = env::var("BIOVAULT_SKIP_INSTALLS").ok();

        // Test 1: When not set, should return false
        env::remove_var("BIOVAULT_SKIP_INSTALLS");
        assert!(
            !super::skip_install_commands(),
            "Expected false when env var not set"
        );

        // Test 2: When set to "1", should return true
        env::set_var("BIOVAULT_SKIP_INSTALLS", "1");
        assert!(
            super::skip_install_commands(),
            "Expected true when env var set to '1'"
        );

        // Test 3: When set to "0", should return false
        env::set_var("BIOVAULT_SKIP_INSTALLS", "0");
        let result = super::skip_install_commands();
        let actual_value =
            env::var("BIOVAULT_SKIP_INSTALLS").unwrap_or_else(|_| "NOT_SET".to_string());
        assert!(
            !result,
            "Expected false when env var set to '0', but got true. Actual env value: '{}'",
            actual_value
        );

        // Restore original state
        if let Some(val) = original {
            env::set_var("BIOVAULT_SKIP_INSTALLS", val);
        } else {
            env::remove_var("BIOVAULT_SKIP_INSTALLS");
        }
    }

    #[test]
    fn java_parse_various_formats() {
        let cases = [
            ("openjdk version \"17.0.2\" 2022-01-18", Some(17)),
            ("java version \"1.8.0_321\"", Some(8)),
            ("openjdk version \"11.0.14\" 2022-01-18", Some(11)),
            ("java version \"21\"", Some(21)),
            ("garbage", None),
        ];
        for (s, want) in cases {
            assert_eq!(parse_java_version(s), want);
        }
    }

    #[test]
    #[serial_test::serial]
    fn google_colab_detection_via_env() {
        // Ensure variable not set
        std::env::remove_var("COLAB_RELEASE_TAG");
        for (k, _) in std::env::vars() {
            if k.starts_with("COLAB_") {
                std::env::remove_var(k);
            }
        }
        assert!(!is_google_colab());
        // Set specific var and detect
        std::env::set_var("COLAB_RELEASE_TAG", "test");
        assert!(is_google_colab());
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    #[serial_test::serial]
    #[cfg_attr(target_os = "windows", ignore = "Colab detection is Linux-specific")]
    fn detect_system_prefers_colab_env() {
        // Clear any existing COLAB_* variables first
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in &keys {
            std::env::remove_var(k);
        }

        // Force Colab-like environment
        std::env::set_var("COLAB_RELEASE_TAG", "1");
        match detect_system() {
            SystemType::GoogleColab => {}
            other => panic!("expected GoogleColab, got {:?}", other),
        }

        // Clean up
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    fn print_syftbox_instructions_runs() {
        // Just ensure it doesn't panic; covers simple printing logic
        super::print_syftbox_instructions();
    }

    #[test]
    #[serial_test::serial]
    fn is_google_colab_detects_prefix_env() {
        std::env::remove_var("COLAB_RELEASE_TAG");
        std::env::set_var("COLAB_FOO", "1");
        assert!(super::is_google_colab());
        std::env::remove_var("COLAB_FOO");
    }

    #[test]
    #[serial_test::serial]
    #[cfg(target_os = "macos")]
    fn detect_system_reports_macos_on_macos() {
        // Ensure no COLAB_* noise affects detection
        std::env::remove_var("COLAB_RELEASE_TAG");
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in keys {
            std::env::remove_var(k);
        }
        match super::detect_system() {
            super::SystemType::MacOs => {}
            other => panic!("expected MacOs, got {:?}", other),
        }
    }

    #[tokio::test]
    #[cfg(target_os = "macos")]
    async fn setup_ubuntu_returns_ok_when_apt_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_ubuntu(vec![], false).await.unwrap();
    }

    #[tokio::test]
    #[cfg(target_os = "macos")]
    async fn setup_arch_returns_ok_when_pacman_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_arch(vec![], false).await.unwrap();
    }

    #[test]
    fn winget_to_choco_mapping() {
        assert_eq!(
            super::map_winget_pkg_to_choco("Microsoft.OpenJDK"),
            "openjdk"
        );
        // Unknown returns empty mapping
        assert_eq!(super::map_winget_pkg_to_choco("Unknown.Package"), "");
    }

    #[tokio::test]
    async fn setup_google_colab_runs_without_panic() {
        let _guard = SkipInstallGuard::new();
        super::setup_google_colab(vec![], false).await.unwrap();
    }

    #[tokio::test]
    async fn setup_macos_returns_ok_without_brew() {
        let _guard = SkipInstallGuard::new();
        // Only run when brew is not available; otherwise skip to avoid invoking installs
        let brew_exists = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v brew >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !brew_exists {
            super::setup_macos(vec![], false).await.unwrap();
        }
    }

    #[tokio::test]
    #[cfg(target_os = "windows")]
    #[cfg_attr(
        not(feature = "e2e-tests"),
        ignore = "runs installer commands; e2e-only"
    )]
    async fn setup_windows_returns_ok_when_tools_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_windows(vec![], false).await.unwrap();
    }

    #[tokio::test]
    #[serial_test::serial]
    #[cfg_attr(
        target_os = "windows",
        ignore = "Colab execution path installs Linux tools"
    )]
    async fn setup_execute_colab_branch() {
        let _guard = SkipInstallGuard::new();
        std::env::set_var("COLAB_RELEASE_TAG", "1");
        super::execute(vec![], false).await.unwrap();
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    fn print_windows_manual_instructions_runs() {
        super::print_windows_manual_instructions();
    }

    #[test]
    fn test_system_type_debug() {
        let s = format!("{:?}", SystemType::MacOs);
        assert_eq!(s, "MacOs");
        let s2 = format!("{:?}", SystemType::GoogleColab);
        assert_eq!(s2, "GoogleColab");
    }

    #[test]
    #[serial_test::serial]
    fn test_detect_system_on_windows() {
        if cfg!(target_os = "windows") {
            match detect_system() {
                SystemType::Windows => {}
                _ => panic!("Expected Windows on windows platform"),
            }
        }
    }

    #[test]
    #[serial_test::serial]
    fn test_skip_install_commands_not_set() {
        std::env::remove_var("BIOVAULT_SKIP_INSTALLS");
        // Just verify it returns a bool without panicking
        let _result = skip_install_commands();
        // Function completes without panic - test passes
    }

    #[test]
    fn test_parse_java_version_edge_cases() {
        assert_eq!(parse_java_version(""), None);
        assert_eq!(parse_java_version("no version here"), None);
        // Just test that it doesn't panic on weird input
        let _ = parse_java_version("version 999");
    }

    #[test]
    fn test_map_winget_pkg_to_choco_all_mappings() {
        // Only Microsoft.OpenJDK is mapped
        assert_eq!(map_winget_pkg_to_choco("Microsoft.OpenJDK"), "openjdk");
        assert_eq!(map_winget_pkg_to_choco("microsoft.openjdk"), "openjdk");
        // Others return empty
        assert_eq!(map_winget_pkg_to_choco("Git.Git"), "");
        assert_eq!(map_winget_pkg_to_choco("RandomPackage"), "");
    }

    #[test]
    fn test_is_google_colab_without_env() {
        std::env::remove_var("COLAB_RELEASE_TAG");
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in keys {
            std::env::remove_var(&k);
        }
        assert!(!is_google_colab());
    }

    // NOTE: Do NOT add unit tests for execute() - it runs actual installation commands
    // and should only be tested via e2e/integration tests.
    // The execute() function performs real system operations (downloads, installs, etc.)
    // which are not suitable for unit tests.
}
