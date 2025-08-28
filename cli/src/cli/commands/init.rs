use crate::config::Config;
use crate::Result;
use std::fs;
use tracing::info;

pub async fn execute(email: &str) -> Result<()> {
    // For testing: allow overriding home directory via environment variable
    let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
        std::path::PathBuf::from(test_home)
    } else {
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
    };

    let biovault_dir = home_dir.join(".biovault");

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
        let syftbox_config = {
            let default_syftbox = home_dir.join(".syftbox").join("config.json");
            if default_syftbox.exists() {
                Some(default_syftbox.to_string_lossy().to_string())
            } else {
                None
            }
        };

        let config = Config {
            email: email.to_string(),
            syftbox_config,
        };

        config.save(&config_file)?;

        println!("✓ BioVault initialized successfully!");
        println!("  Configuration saved to: {}", config_file.display());
        println!("  Email: {}", email);

        // Create env/default directory and copy templates
        let env_dir = biovault_dir.join("env").join("default");
        if !env_dir.exists() {
            fs::create_dir_all(&env_dir)?;
            info!("Created environment directory: {:?}", env_dir);
        }

        // Copy template.nf
        let template_nf_content = include_str!("../../templates/template.nf");
        let template_nf_path = env_dir.join("template.nf");
        fs::write(&template_nf_path, template_nf_content)?;
        info!("Created template.nf at: {:?}", template_nf_path);

        // Copy nextflow.config
        let nextflow_config_content = include_str!("../../templates/nextflow.config");
        let nextflow_config_path = env_dir.join("nextflow.config");
        fs::write(&nextflow_config_path, nextflow_config_content)?;
        info!("Created nextflow.config at: {:?}", nextflow_config_path);

        println!("✓ Nextflow templates installed to: {}", env_dir.display());
    }

    Ok(())
}
