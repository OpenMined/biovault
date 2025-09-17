use crate::config::{get_biovault_home, is_syftbox_env, Config};
use crate::Result;
use dialoguer::{theme::ColorfulTheme, Confirm, Input};
use std::env;
use std::fs;
use std::io::IsTerminal;
use tracing::info;

pub async fn execute(email: Option<&str>) -> Result<()> {
    // Get the BioVault home directory (respects env vars)
    let biovault_dir = get_biovault_home()?;

    // Determine email with priority:
    // 1. CLI argument (if provided)
    // 2. SYFTBOX_EMAIL env var (confirm on TTY, auto-accept non-interactive)
    // 3. Interactive prompt if on a TTY
    // 4. Error if neither provided and not a TTY
    let email = if let Some(email) = email {
        email.to_string()
    } else if let Ok(syftbox_email) = env::var("SYFTBOX_EMAIL") {
        if std::io::stdin().is_terminal() {
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
    } else if std::io::stdin().is_terminal() {
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

    if !biovault_dir.exists() {
        fs::create_dir_all(&biovault_dir)?;
        info!("Created directory: {:?}", biovault_dir);
    }

    let config_file = biovault_dir.join("config.yaml");

    if config_file.exists() {
        println!(
            "Configuration file already exists at: {}",
            config_file.display()
        );
        println!("Skipping initialization.");
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
            version: Some("0.1.16".to_string()),
        };

        config.save(&config_file)?;

        println!("✓ BioVault initialized successfully!");
        println!("  Configuration saved to: {}", config_file.display());
        println!("  Email: {}", email);
        if let Some(ref syftbox_cfg) = syftbox_config {
            println!("  SyftBox config: {}", syftbox_cfg);
        }

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

        println!("✓ Nextflow templates installed:");
        println!("  - Default templates: {}", default_dir.display());
        println!("  - SNP templates: {}", snp_dir.display());

        // Initialize SyftBox RPC folders for messaging if SyftBox is configured
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
                    "⚠️  Skipped SyftBox RPC init (no data dir): {}\n    Hint: set SYFTBOX_DATA_DIR or configure ~/.syftbox/config.json and re-run 'bv init'",
                    e
                );
            }
        }
    }

    Ok(())
}
