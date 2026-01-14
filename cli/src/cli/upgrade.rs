use crate::config::{get_biovault_home, Config};
use crate::syftbox::storage::SyftBoxStorage;
use std::fs;
use std::path::Path;
use tracing::info;

const CURRENT_VERSION: &str = "0.1.28";

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

    // Copy/update sheet templates
    let sheet_dir = env_dir.join("sheet");
    if !sheet_dir.exists() {
        fs::create_dir_all(&sheet_dir)?;
    }

    // Update sheet template.nf
    let sheet_template_nf_content = include_str!("../templates/sheet/template.nf");
    let sheet_template_nf_path = sheet_dir.join("template.nf");
    fs::write(&sheet_template_nf_path, sheet_template_nf_content)?;
    info!("Updated sheet template.nf");

    // Update sheet nextflow.config
    let sheet_nextflow_config_content = include_str!("../templates/sheet/nextflow.config");
    let sheet_nextflow_config_path = sheet_dir.join("nextflow.config");
    fs::write(&sheet_nextflow_config_path, sheet_nextflow_config_content)?;
    info!("Updated sheet nextflow.config");

    println!("✓ Templates updated:");
    println!("  - Default templates: {}", default_dir.display());
    println!("  - SNP templates: {}", snp_dir.display());
    println!("  - Sheet templates: {}", sheet_dir.display());

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

    let storage = SyftBoxStorage::new(&data_dir);

    // Fix RPC permission file
    let rpc_permission_file = data_dir
        .join("datasites")
        .join(&config.email)
        .join("app_data")
        .join("biovault")
        .join("rpc")
        .join("syft.pub.yaml");

    if rpc_permission_file.exists() {
        fix_yaml_file(&storage, &rpc_permission_file)?;
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
        fix_yaml_file(&storage, &app_permission_file)?;
        println!(
            "  ✓ Fixed app permission file: {}",
            app_permission_file.display()
        );
    }

    Ok(())
}

fn fix_yaml_file(storage: &SyftBoxStorage, path: &Path) -> anyhow::Result<()> {
    // Just replace the entire file with the correct format
    // This is safer than trying to fix various indentation issues

    // Determine which template to use based on the path
    let s = path.to_string_lossy();
    let is_rpc =
        s.contains("/rpc/") || s.contains("\\rpc\\") || s.ends_with("/rpc") || s.ends_with("\\rpc");
    let new_content = if is_rpc {
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
    storage.write_plaintext_file(path, new_content.as_bytes(), true)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use tempfile::TempDir;

    #[test]
    fn version_compare_works() {
        assert!(version_less_than("0.1.0", "0.1.27"));
        assert!(!version_less_than("0.2.0", "0.1.27"));
        assert!(!version_less_than("0.1.27", "0.1.27"));
        assert!(version_less_than("0.1", "0.1.1"));
    }

    #[test]
    fn upgrade_writes_templates_and_fixes_permissions() {
        let tmp = TempDir::new().unwrap();
        // Point BIOVAULT home to temp
        config::set_test_biovault_home(tmp.path().join(".bv"));
        std::fs::create_dir_all(tmp.path().join(".bv")).unwrap();

        // Seed config with old version and syftbox config path
        let email = "user@example.com";
        let cfg_path = tmp.path().join(".bv/config.yaml");
        let syft_dir = tmp.path().join("syft");
        std::fs::create_dir_all(&syft_dir).unwrap();
        let syft_cfg = syft_dir.join("config.json");
        let data_root = tmp.path().join("data");
        std::fs::create_dir_all(&data_root).unwrap();
        std::fs::write(
            &syft_cfg,
            serde_json::json!({"data_dir": data_root.to_string_lossy()}).to_string(),
        )
        .unwrap();

        let cfg = Config {
            email: email.into(),
            syftbox_config: Some(syft_cfg.to_string_lossy().to_string()),
            version: Some("0.1.0".into()),
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        cfg.save(&cfg_path).unwrap();

        // Seed permission files with dummy content to be fixed
        let rpc_perm = data_root
            .join("datasites")
            .join(email)
            .join("app_data/biovault/rpc/syft.pub.yaml");
        let app_perm = data_root
            .join("datasites")
            .join(email)
            .join("app_data/biovault/syft.pub.yaml");
        std::fs::create_dir_all(rpc_perm.parent().unwrap()).unwrap();
        std::fs::create_dir_all(app_perm.parent().unwrap()).unwrap();
        std::fs::write(&rpc_perm, "bad").unwrap();
        std::fs::write(&app_perm, "bad").unwrap();

        // Run upgrade
        check_and_upgrade().unwrap();

        // Config version updated
        let updated = Config::from_file(&cfg_path).unwrap();
        assert_eq!(updated.version.as_deref(), Some(CURRENT_VERSION));

        // Templates written
        let env_dir = tmp.path().join(".bv/env");
        assert!(env_dir.join("default/template.nf").exists());
        assert!(env_dir.join("default/nextflow.config").exists());
        assert!(env_dir.join("snp/template.nf").exists());
        assert!(env_dir.join("snp/nextflow.config").exists());

        // Permission files fixed with expected content prefixes
        let rpc_fixed = std::fs::read_to_string(&rpc_perm).unwrap();
        assert!(rpc_fixed.contains("**/*.request"));
        let app_fixed = std::fs::read_to_string(&app_perm).unwrap();
        assert!(app_fixed.contains("rules:"));

        // Cleanup thread-local override
        config::clear_test_biovault_home();
    }
}
