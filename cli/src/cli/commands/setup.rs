use super::check::DependencyConfig;
use crate::Result;
use std::env;
use std::process::Command;

#[derive(Debug)]
enum SystemType {
    GoogleColab,
    Unknown,
}

impl SystemType {
    fn as_str(&self) -> &str {
        match self {
            SystemType::GoogleColab => "google_colab",
            SystemType::Unknown => "unknown",
        }
    }
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
        SystemType::Unknown => {
            println!("ℹ️  System type not detected or not supported for automated setup");
            println!("   This command currently supports:");
            println!("   - Google Colab");
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
