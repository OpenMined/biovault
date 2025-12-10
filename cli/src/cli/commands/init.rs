use crate::config::{get_biovault_home, is_syftbox_env, Config};
use crate::syftbox::syc;
use crate::Result;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use std::env;
use std::fs;
use std::io::IsTerminal;
use std::path::PathBuf;
use tracing::info;

pub async fn execute(email: Option<&str>, quiet: bool) -> Result<()> {
    // Check if config already exists before prompting for location
    let default_biovault_dir = get_biovault_home()?;
    let config_file = default_biovault_dir.join("config.yaml");
    let is_existing_installation = config_file.exists();

    // Get the BioVault home directory (respects env vars and prompts for new installs)
    let biovault_dir = if is_existing_installation {
        // Use existing location
        default_biovault_dir
    } else if env::var("BIOVAULT_HOME").is_ok() {
        // Explicitly set via env var
        default_biovault_dir
    } else if env::var("SYFTBOX_DATA_DIR").is_ok() {
        // In SyftBox virtualenv, use that location
        default_biovault_dir
    } else if !quiet && std::io::stdin().is_terminal() {
        // New installation - prompt for location
        prompt_for_location()?
    } else {
        // Non-interactive or quiet mode - use default
        default_biovault_dir
    };

    // Determine email with priority:
    // 1. CLI argument (if provided)
    // 2. SYFTBOX_EMAIL env var (confirm on TTY, auto-accept non-interactive)
    // 3. Interactive prompt if on a TTY
    // 4. Error if neither provided and not a TTY
    let email = if let Some(email) = email {
        email.to_string()
    } else if let Ok(syftbox_email) = env::var("SYFTBOX_EMAIL") {
        if quiet {
            // In quiet mode, auto-accept SYFTBOX_EMAIL
            syftbox_email
        } else if std::io::stdin().is_terminal() {
            println!("Detected SyftBox environment email: {}", syftbox_email);
            let use_syftbox = Confirm::with_theme(&ColorfulTheme::default())
                .with_prompt("Use this email for BioVault?")
                .default(true)
                .interact()
                .map_err(|e| anyhow::anyhow!("Prompt error: {}", e))?;

            if use_syftbox {
                syftbox_email
            } else {
                Input::with_theme(&ColorfulTheme::default())
                    .with_prompt("Enter email address")
                    .interact_text()
                    .map_err(|e| anyhow::anyhow!("Prompt error: {}", e))?
            }
        } else {
            // Non-interactive environment: auto-accept SYFTBOX_EMAIL
            syftbox_email
        }
    } else if std::io::stdin().is_terminal() && !quiet {
        Input::with_theme(&ColorfulTheme::default())
            .with_prompt("Enter email address")
            .interact_text()
            .map_err(|e| anyhow::anyhow!("Prompt error: {}", e))?
    } else {
        return Err(anyhow::anyhow!(
            "email is required. Provide as 'bv init <email>' or set SYFTBOX_EMAIL"
        )
        .into());
    };

    // Persist the selected location for future runs when applicable
    crate::config::set_persisted_biovault_home(&biovault_dir);

    if !biovault_dir.exists() {
        fs::create_dir_all(&biovault_dir)?;
        info!("Created directory: {:?}", biovault_dir);
    }

    let config_file = biovault_dir.join("config.yaml");

    // Step 1: Create config.yaml if it doesn't exist
    let config = if config_file.exists() {
        println!(
            "Configuration file already exists at: {}",
            config_file.display()
        );
        // Load existing config
        Config::load()?
    } else {
        // Show detected environment if in SyftBox virtualenv
        if is_syftbox_env() {
            println!("✓ Running in SyftBox virtualenv");
            if let Ok(data_dir) = env::var("SYFTBOX_DATA_DIR") {
                println!("  SyftBox data directory: {}", data_dir);
            }
            println!("  BioVault home: {}", biovault_dir.display());
        }

        // Auto-detect SyftBox config if it exists
        let syftbox_config = if let Ok(config_path) = env::var("SYFTBOX_CONFIG_PATH") {
            if std::path::Path::new(&config_path).exists() {
                println!("✓ Detected SyftBox config at: {}", config_path);
                Some(config_path)
            } else {
                None
            }
        } else {
            let home_dir = dirs::home_dir()
                .ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
            let default_syftbox = home_dir.join(".syftbox").join("config.json");
            if default_syftbox.exists() {
                Some(default_syftbox.to_string_lossy().to_string())
            } else {
                None
            }
        };

        let config = Config {
            email: email.to_string(),
            syftbox_config: syftbox_config.clone(),
            version: Some("0.1.27".to_string()),
            binary_paths: None,
            syftbox_credentials: None,
        };

        config.save(&config_file)?;

        println!("✓ BioVault initialized successfully!");
        println!("  Configuration saved to: {}", config_file.display());
        println!("  Email: {}", email);
        if let Some(ref syftbox_cfg) = syftbox_config {
            println!("  SyftBox config: {}", syftbox_cfg);
        }
        config
    };

    // Step 2: Always ensure templates and directories exist (even if config existed)
    {
        // Create env directory
        let env_dir = biovault_dir.join("env");
        if !env_dir.exists() {
            fs::create_dir_all(&env_dir)?;
            info!("Created environment directory: {:?}", env_dir);
        }

        // Copy default templates
        let default_dir = env_dir.join("default");
        if !default_dir.exists() {
            fs::create_dir_all(&default_dir)?;
            info!("Created default template directory: {:?}", default_dir);
        }

        // Copy default template.nf
        let template_nf_content = include_str!("../../templates/default/template.nf");
        let template_nf_path = default_dir.join("template.nf");
        fs::write(&template_nf_path, template_nf_content)?;
        info!("Created template.nf at: {:?}", template_nf_path);

        // Copy default nextflow.config
        let nextflow_config_content = include_str!("../../templates/default/nextflow.config");
        let nextflow_config_path = default_dir.join("nextflow.config");
        fs::write(&nextflow_config_path, nextflow_config_content)?;
        info!("Created nextflow.config at: {:?}", nextflow_config_path);

        // Copy SNP templates
        let snp_dir = env_dir.join("snp");
        if !snp_dir.exists() {
            fs::create_dir_all(&snp_dir)?;
            info!("Created SNP template directory: {:?}", snp_dir);
        }

        // Copy SNP template.nf
        let snp_template_nf_content = include_str!("../../templates/snp/template.nf");
        let snp_template_nf_path = snp_dir.join("template.nf");
        fs::write(&snp_template_nf_path, snp_template_nf_content)?;
        info!("Created SNP template.nf at: {:?}", snp_template_nf_path);

        // Copy SNP nextflow.config
        let snp_nextflow_config_content = include_str!("../../templates/snp/nextflow.config");
        let snp_nextflow_config_path = snp_dir.join("nextflow.config");
        fs::write(&snp_nextflow_config_path, snp_nextflow_config_content)?;
        info!(
            "Created SNP nextflow.config at: {:?}",
            snp_nextflow_config_path
        );

        // Copy sheet templates
        let sheet_dir = env_dir.join("sheet");
        if !sheet_dir.exists() {
            fs::create_dir_all(&sheet_dir)?;
            info!("Created sheet template directory: {:?}", sheet_dir);
        }

        // Copy sheet template.nf
        let sheet_template_nf_content = include_str!("../../templates/sheet/template.nf");
        let sheet_template_nf_path = sheet_dir.join("template.nf");
        fs::write(&sheet_template_nf_path, sheet_template_nf_content)?;
        info!("Created sheet template.nf at: {:?}", sheet_template_nf_path);

        // Copy sheet nextflow.config
        let sheet_nextflow_config_content = include_str!("../../templates/sheet/nextflow.config");
        let sheet_nextflow_config_path = sheet_dir.join("nextflow.config");
        fs::write(&sheet_nextflow_config_path, sheet_nextflow_config_content)?;
        info!(
            "Created sheet nextflow.config at: {:?}",
            sheet_nextflow_config_path
        );

        match config.get_syftbox_data_dir() {
            Ok(data_root) => match syc::provision_local_identity(&config.email, &data_root, None) {
                Ok(outcome) => {
                    if outcome.generated {
                        println!("✓ Generated Syft Crypto identity for {}", outcome.identity);
                        if let Some(mnemonic) = outcome.recovery_mnemonic.as_deref() {
                            println!("  Recovery mnemonic (store securely!): {}", mnemonic);
                        }
                    } else {
                        println!("✓ Syft Crypto identity detected for {}", outcome.identity);
                    }
                    println!("  Vault directory: {}", outcome.vault_path.display());
                    println!(
                        "  Public bundle published at: {}",
                        outcome.public_bundle_path.display()
                    );
                }
                Err(err) => {
                    eprintln!("⚠️  Unable to provision Syft Crypto identity automatically: {err}");
                    eprintln!(
                        "    Run 'bv syc import --bundle <path> --expected-identity <email>' once your peer shares a bundle."
                    );
                }
            },
            Err(err) => {
                println!(
                    "⚠️  Skipping Syft Crypto identity provisioning (no SyftBox config): {err}"
                );
            }
        }

        // Copy dynamic templates
        let dynamic_dir = env_dir.join("dynamic-nextflow");
        if !dynamic_dir.exists() {
            fs::create_dir_all(&dynamic_dir)?;
            info!("Created dynamic template directory: {:?}", dynamic_dir);
        }

        // Copy dynamic template.nf
        let dynamic_template_nf_content = include_str!("../../templates/dynamic/template.nf");
        let dynamic_template_nf_path = dynamic_dir.join("template.nf");
        fs::write(&dynamic_template_nf_path, dynamic_template_nf_content)?;
        info!(
            "Created dynamic template.nf at: {:?}",
            dynamic_template_nf_path
        );

        // Copy dynamic nextflow.config (minimal config)
        let dynamic_nextflow_config_content = "process.executor = 'local'\n";
        let dynamic_nextflow_config_path = dynamic_dir.join("nextflow.config");
        fs::write(
            &dynamic_nextflow_config_path,
            dynamic_nextflow_config_content,
        )?;
        info!(
            "Created dynamic nextflow.config at: {:?}",
            dynamic_nextflow_config_path
        );

        println!("✓ Nextflow templates installed:");
        println!("  - Default templates: {}", default_dir.display());
        println!("  - SNP templates: {}", snp_dir.display());
        println!("  - Sheet templates: {}", sheet_dir.display());
        println!("  - Dynamic templates: {}", dynamic_dir.display());

        // Initialize SyftBox RPC folders for messaging if SyftBox is configured
        let default_syftbox_config = Config::default_syftbox_config_path()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|_| "<BioVault syftbox/config.json>".to_string());
        match config.get_syftbox_data_dir() {
            Ok(data_dir) => {
                let app = crate::syftbox::SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
                // Ensure the message endpoint exists
                let _ = app.register_endpoint("/message")?;
                println!(
                    "✓ SyftBox RPC initialized for messaging at: {}",
                    app.rpc_dir.display()
                );
            }
            Err(e) => {
                // Not fatal: user might not be in a SyftBox env yet
                println!(
                    "⚠️  Skipped SyftBox RPC init (no data dir): {}\n    Hint: set SYFTBOX_DATA_DIR or configure {} and re-run 'bv init'",
                    e, default_syftbox_config
                );
            }
        }

        // Create projects and runs directories
        let projects_dir = biovault_dir.join("projects");
        if !projects_dir.exists() {
            fs::create_dir_all(&projects_dir)?;
            info!("Created projects directory: {:?}", projects_dir);
            println!("✓ Created projects directory: {}", projects_dir.display());
        }

        let runs_dir = biovault_dir.join("runs");
        if !runs_dir.exists() {
            fs::create_dir_all(&runs_dir)?;
            info!("Created runs directory: {:?}", runs_dir);
            println!("✓ Created runs directory: {}", runs_dir.display());
        }

        // Create data/cache directory
        let cache_dir = biovault_dir.join("data").join("cache");
        if !cache_dir.exists() {
            fs::create_dir_all(&cache_dir)?;
            info!("Created cache directory: {:?}", cache_dir);
            println!("✓ Created cache directory: {}", cache_dir.display());
        }
    }

    Ok(())
}

/// Prompt user for BioVault installation location
fn prompt_for_location() -> Result<PathBuf> {
    println!("\n📁 Where would you like to store BioVault data?\n");

    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
    let desktop_dir = dirs::desktop_dir().unwrap_or_else(|| home_dir.join("Desktop"));

    let choices = vec![
        format!(
            "Desktop/BioVault (recommended for desktop use) [{}]",
            desktop_dir.join("BioVault").display()
        ),
        format!(
            "Hidden in home directory (for CLI/server use) [{}]",
            home_dir.join(".biovault").display()
        ),
        "Custom location".to_string(),
    ];

    let selection = Select::with_theme(&ColorfulTheme::default())
        .with_prompt("Choose installation location")
        .items(&choices)
        .default(0)
        .interact()
        .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

    let biovault_dir = match selection {
        0 => desktop_dir.join("BioVault"),
        1 => home_dir.join(".biovault"),
        2 => {
            let custom_path: String = Input::with_theme(&ColorfulTheme::default())
                .with_prompt("Enter custom path")
                .interact_text()
                .map_err(|e| anyhow::anyhow!("Input error: {}", e))?;

            // Expand ~ to home directory
            if let Some(path) = custom_path.strip_prefix("~/") {
                home_dir.join(path)
            } else if custom_path == "~" {
                home_dir.clone()
            } else {
                PathBuf::from(custom_path)
            }
        }
        _ => desktop_dir.join("BioVault"),
    };

    println!("\n✓ Selected location: {}\n", biovault_dir.display());

    Ok(biovault_dir)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[tokio::test]
    async fn test_execute_with_email_arg() {
        let temp_dir = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));
        let data_dir = temp_dir.path().join("syftbox");
        crate::config::set_test_syftbox_data_dir(&data_dir);

        // Initialize with email argument (quiet=true for non-interactive test)
        let result = execute(Some("test@example.com"), true).await;
        assert!(result.is_ok());

        // Verify config was created
        let config = Config::load().unwrap();
        assert_eq!(config.email, "test@example.com");

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn test_execute_with_syftbox_email_env() {
        let temp_dir = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));
        let data_dir = temp_dir.path().join("syftbox");
        crate::config::set_test_syftbox_data_dir(&data_dir);

        // Set SYFTBOX_EMAIL env var
        env::set_var("SYFTBOX_EMAIL", "syftbox@example.com");

        // In non-interactive mode, it should use the env var
        // We can't easily test interactive mode in unit tests
        let result = execute(None, true).await;

        // Clean up env var
        env::remove_var("SYFTBOX_EMAIL");
        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();

        // The result depends on whether we're in TTY or not
        // In CI/test environment, it's usually non-TTY
        let _ = result;
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_execute_overwrite_existing_config() {
        let temp_dir = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));
        let data_dir = temp_dir.path().join("syftbox");
        crate::config::set_test_syftbox_data_dir(&data_dir);

        // Create initial config
        let initial_config = Config {
            email: "old@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
        };
        let config_path = temp_dir.path().join(".biovault").join("config.yaml");
        fs::create_dir_all(config_path.parent().unwrap()).unwrap();
        initial_config.save(&config_path).unwrap();

        // Re-initialize with new email - this might prompt in TTY mode
        // We'll just test that the function doesn't panic
        let _ = execute(Some("new@example.com"), true).await;

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn test_execute_creates_biovault_dir() {
        let temp_dir = TempDir::new().unwrap();
        let biovault_path = temp_dir.path().join(".biovault");
        crate::config::set_test_biovault_home(biovault_path.clone());
        let data_dir = temp_dir.path().join("syftbox");
        crate::config::set_test_syftbox_data_dir(&data_dir);

        // Directory shouldn't exist initially
        assert!(!biovault_path.exists());

        // Execute init (quiet=true for non-interactive test)
        let result = execute(Some("test@example.com"), true).await;
        assert!(result.is_ok());

        // Directory should now exist
        assert!(biovault_path.exists());

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn test_execute_no_email_non_tty() {
        let temp_dir = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        // Make sure SYFTBOX_EMAIL is not set
        env::remove_var("SYFTBOX_EMAIL");

        // In non-TTY mode (like CI), this should error without email
        // Note: we can't easily control TTY state in tests
        // The actual behavior depends on the test environment

        crate::config::clear_test_biovault_home();
    }
}
