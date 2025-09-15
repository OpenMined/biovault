use super::check::DependencyConfig;
use crate::Result;
use std::env;
use std::process::Command;

#[derive(Debug)]
enum SystemType {
    GoogleColab,
    MacOs,
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
        SystemType::Unknown => {
            println!("ℹ️  System type not detected or not supported for automated setup");
            println!("   This command currently supports:");
            println!("   - Google Colab");
            println!("   - macOS (Homebrew)");
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
