use super::check::DependencyConfig;
use crate::Result;
use anyhow::anyhow;
use std::env;
use std::process::{Command, Stdio};

#[derive(Debug)]
enum SystemType {
    GoogleColab,
    MacOs,
    Ubuntu,
    ArchLinux,
    Windows,
    Unknown,
}

pub async fn execute() -> Result<()> {
    println!("BioVault Environment Setup");
    println!("==========================\n");

    let system_type = detect_system();

    match system_type {
        SystemType::GoogleColab => {
            println!("✓ Detected Google Colab environment");
            setup_google_colab().await?;
        }
        SystemType::MacOs => {
            println!("✓ Detected macOS environment");
            setup_macos().await?;
        }
        SystemType::Ubuntu => {
            println!("✓ Detected Ubuntu/Debian environment");
            setup_ubuntu().await?;
        }
        SystemType::ArchLinux => {
            println!("✓ Detected Arch Linux environment");
            setup_arch().await?;
        }
        SystemType::Windows => {
            println!("✓ Detected Windows environment");
            setup_windows().await?;
        }
        SystemType::Unknown => {
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

    Ok(())
}

fn detect_system() -> SystemType {
    // Check for Google Colab environment variables
    if is_google_colab() {
        return SystemType::GoogleColab;
    }

    // Detect macOS
    if std::env::consts::OS == "macos" {
        return SystemType::MacOs;
    }

    // Detect Windows
    if std::env::consts::OS == "windows" {
        return SystemType::Windows;
    }

    // Detect Linux distributions
    if std::env::consts::OS == "linux" {
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

async fn setup_google_colab() -> Result<()> {
    println!("\nSetting up Google Colab environment...\n");

    // Load the deps.yaml file to get environment-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
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

async fn setup_macos() -> Result<()> {
    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up macOS environment...\n");

    // Check for Homebrew (most macOS installs use brew per deps.yaml)
    let brew_exists = Command::new("sh")
        .arg("-c")
        .arg("command -v brew >/dev/null 2>&1")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if !brew_exists {
        println!("❌ Homebrew not found.");
        println!("Please install Homebrew first: https://brew.sh");
        println!("Install command:");
        println!("  /bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"");
        println!("\nAfter installing Homebrew, re-run: bv setup\n");
        // Still provide SyftBox link
        print_syftbox_instructions();
        return Ok(());
    }

    // Load deps.yaml and execute macOS-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
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
                    // Decide if install is necessary
                    let mut need_install = true;

                    // If verify_command is available, try it first
                    if let Some(verify_cmd) = &env_config.verify_command {
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
                    } else {
                        // Fallback to which for simple presence
                        if which::which(&dep.name).is_ok() {
                            need_install = false;
                        }
                    }

                    // For Java, also enforce min_version if specified
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

                    let mut all_succeeded = true;

                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        let output = Command::new("sh")
                            .arg("-c")
                            .arg(cmd)
                            .stdout(Stdio::piped())
                            .stderr(Stdio::piped())
                            .output();
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
                        if let Some(verify_cmd) = &env_config.verify_command {
                            print!("   Verifying installation... ");
                            let output = Command::new("sh").arg("-c").arg(verify_cmd).output();
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

    println!("\nNotes:");
    println!("- If this is your first time installing Docker Desktop, open it once to finish setup and grant permissions.");
    println!("- You may need to ensure the OpenJDK 17 binaries are on PATH. Brew usually prints a caveat like adding a PATH export.");

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

async fn setup_ubuntu() -> Result<()> {
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

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
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

async fn setup_arch() -> Result<()> {
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

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
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

async fn setup_windows() -> Result<()> {
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

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
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
    println!("Nextflow: Download from https://www.nextflow.io/ or use PowerShell script");
    println!("SyftBox: Download from https://github.com/OpenMined/syftbox/releases/latest");
}

// Map common WinGet package IDs to Chocolatey package names for fallback
fn map_winget_pkg_to_choco(pkg: &str) -> &'static str {
    match pkg.to_ascii_lowercase().as_str() {
        // Java/OpenJDK
        // WinGet: Microsoft.OpenJDK => Chocolatey: openjdk (generic)
        "microsoft.openjdk" => "openjdk",
        // Add other mappings here as needed
        _ => "",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
