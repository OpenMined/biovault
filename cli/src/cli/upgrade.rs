use crate::config::{get_biovault_home, Config};
use std::fs;
use std::path::Path;
use tracing::info;

const CURRENT_VERSION: &str = "0.1.27";

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

        // Fix permission file YAML indentation (added in 0.1.27)
        if version_less_than(&stored_version, "0.1.27") {
            fix_permission_yaml_indentation(&config)?;
        }

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

fn version_less_than(v1: &str, v2: &str) -> bool {
    let parse_version =
        |v: &str| -> Vec<u32> { v.split('.').filter_map(|s| s.parse::<u32>().ok()).collect() };

    let v1_parts = parse_version(v1);
    let v2_parts = parse_version(v2);

    for i in 0..3 {
        let v1_part = v1_parts.get(i).unwrap_or(&0);
        let v2_part = v2_parts.get(i).unwrap_or(&0);
        if v1_part < v2_part {
            return true;
        } else if v1_part > v2_part {
            return false;
        }
    }
    false
}

fn fix_permission_yaml_indentation(config: &Config) -> anyhow::Result<()> {
    println!("🔧 Fixing YAML permission file indentation...");

    // Try to get the SyftBox data directory
    let data_dir = match config.get_syftbox_data_dir() {
        Ok(dir) => dir,
        Err(e) => {
            info!("Skipping permission file fix (no SyftBox data dir): {}", e);
            return Ok(());
        }
    };

    // Fix RPC permission file
    let rpc_permission_file = data_dir
        .join("datasites")
        .join(&config.email)
        .join("app_data")
        .join("biovault")
        .join("rpc")
        .join("syft.pub.yaml");

    if rpc_permission_file.exists() {
        fix_yaml_file(&rpc_permission_file)?;
        println!(
            "  ✓ Fixed RPC permission file: {}",
            rpc_permission_file.display()
        );
    }

    // Fix app-level permission file
    let app_permission_file = data_dir
        .join("datasites")
        .join(&config.email)
        .join("app_data")
        .join("biovault")
        .join("syft.pub.yaml");

    if app_permission_file.exists() {
        fix_yaml_file(&app_permission_file)?;
        println!(
            "  ✓ Fixed app permission file: {}",
            app_permission_file.display()
        );
    }

    Ok(())
}

fn fix_yaml_file(path: &Path) -> anyhow::Result<()> {
    // Just replace the entire file with the correct format
    // This is safer than trying to fix various indentation issues

    // Determine which template to use based on the path
    let new_content = if path.to_string_lossy().contains("/rpc/") {
        // RPC permission file - allow read/write for requests
        r#"rules:
  - pattern: '**/*.request'
    access:
      admin: []
      read:
        - '*'
      write:
        - '*'
"#
    } else {
        // App-level permission file - allow read for user's files
        r#"rules:
  - pattern: '{{.UserEmail}}/*'
    access:
      admin: []
      read:
        - '*'
      write: []
"#
    };

    // Always write the correct content
    fs::write(path, new_content)?;

    Ok(())
}
