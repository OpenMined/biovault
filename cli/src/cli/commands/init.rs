use crate::Result;
use serde::{Deserialize, Serialize};
use std::fs;
use tracing::info;

#[derive(Debug, Serialize, Deserialize)]
struct BioVaultConfig {
    email: String,
}

pub async fn execute(email: &str) -> Result<()> {
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;

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
        let config = BioVaultConfig {
            email: email.to_string(),
        };

        let yaml_content = serde_yaml::to_string(&config)?;
        fs::write(&config_file, yaml_content)?;

        println!("✓ BioVault initialized successfully!");
        println!("  Configuration saved to: {}", config_file.display());
        println!("  Email: {}", email);
    }

    Ok(())
}
