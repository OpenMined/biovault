use crate::config::{get_biovault_home, Config};
use std::fs;
use tracing::info;

const CURRENT_VERSION: &str = "0.1.16";

pub fn check_and_upgrade() -> anyhow::Result<()> {
    let biovault_home = get_biovault_home()?;
    let config_file = biovault_home.join("config.yaml");

    if !config_file.exists() {
        // Fresh install, no upgrade needed
        return Ok(());
    }

    let mut config = Config::from_file(&config_file)?;
    let stored_version = config
        .version
        .clone()
        .unwrap_or_else(|| "0.0.0".to_string());

    if stored_version != CURRENT_VERSION {
        info!(
            "Upgrading BioVault from {} to {}",
            stored_version, CURRENT_VERSION
        );

        // Perform upgrade tasks
        upgrade_templates(&biovault_home)?;

        // Update version in config
        config.version = Some(CURRENT_VERSION.to_string());
        config.save(&config_file)?;

        println!("✓ BioVault upgraded to version {}", CURRENT_VERSION);
    }

    Ok(())
}

fn upgrade_templates(biovault_home: &std::path::Path) -> anyhow::Result<()> {
    let env_dir = biovault_home.join("env");

    // Ensure env directory exists
    if !env_dir.exists() {
        fs::create_dir_all(&env_dir)?;
    }

    // Copy/update default templates
    let default_dir = env_dir.join("default");
    if !default_dir.exists() {
        fs::create_dir_all(&default_dir)?;
    }

    // Update default template.nf
    let template_nf_content = include_str!("../templates/default/template.nf");
    let template_nf_path = default_dir.join("template.nf");
    fs::write(&template_nf_path, template_nf_content)?;
    info!("Updated default template.nf");

    // Update default nextflow.config
    let nextflow_config_content = include_str!("../templates/default/nextflow.config");
    let nextflow_config_path = default_dir.join("nextflow.config");
    fs::write(&nextflow_config_path, nextflow_config_content)?;
    info!("Updated default nextflow.config");

    // Copy/update SNP templates
    let snp_dir = env_dir.join("snp");
    if !snp_dir.exists() {
        fs::create_dir_all(&snp_dir)?;
    }

    // Update SNP template.nf
    let snp_template_nf_content = include_str!("../templates/snp/template.nf");
    let snp_template_nf_path = snp_dir.join("template.nf");
    fs::write(&snp_template_nf_path, snp_template_nf_content)?;
    info!("Updated SNP template.nf");

    // Update SNP nextflow.config
    let snp_nextflow_config_content = include_str!("../templates/snp/nextflow.config");
    let snp_nextflow_config_path = snp_dir.join("nextflow.config");
    fs::write(&snp_nextflow_config_path, snp_nextflow_config_content)?;
    info!("Updated SNP nextflow.config");

    println!("✓ Templates updated:");
    println!("  - Default templates: {}", default_dir.display());
    println!("  - SNP templates: {}", snp_dir.display());

    Ok(())
}
