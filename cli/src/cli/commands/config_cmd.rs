use crate::config::Config;
use crate::error::Error;
use anyhow::Context;
use serde::{Deserialize, Serialize};

fn default_syftbox_path_str() -> String {
    Config::default_syftbox_config_path()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "<BioVault syftbox/config.json>".to_string())
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ConfigDisplay {
    pub config_file: String,
    pub email: String,
    pub syftbox_config: String,
    pub syftbox_data_dir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binary_paths: Option<crate::config::BinaryPaths>,
}

pub async fn show(json: bool) -> crate::error::Result<()> {
    let config = Config::load()?;
    let config_path = Config::get_config_path()?;

    let syftbox_config_str = if let Some(ref syftbox_path) = config.syftbox_config {
        syftbox_path.clone()
    } else {
        format!("{} (default)", default_syftbox_path_str())
    };

    let syftbox_data_dir_str = match config.get_syftbox_data_dir() {
        Ok(data_dir) => Some(data_dir.display().to_string()),
        Err(_) => None,
    };

    let display = ConfigDisplay {
        config_file: config_path.display().to_string(),
        email: config.email.clone(),
        syftbox_config: syftbox_config_str.clone(),
        syftbox_data_dir: syftbox_data_dir_str.clone(),
        binary_paths: config.binary_paths.clone(),
    };

    if json {
        println!("{}", serde_json::to_string_pretty(&display)?);
    } else {
        println!("BioVault Configuration");
        println!("----------------------");
        println!("Config file: {}", config_path.display());
        println!("Email: {}", config.email);
        println!("SyftBox config: {}", syftbox_config_str);

        match config.get_syftbox_data_dir() {
            Ok(data_dir) => println!("SyftBox data dir: {}", data_dir.display()),
            Err(e) => println!("SyftBox data dir: <error: {}>", e),
        }

        if let Some(ref bp) = config.binary_paths {
            println!("\nCustom Binary Paths:");
            if let Some(ref java) = bp.java {
                println!("  Java: {}", java);
            }
            if let Some(ref docker) = bp.docker {
                println!("  Docker: {}", docker);
            }
            if let Some(ref nextflow) = bp.nextflow {
                println!("  Nextflow: {}", nextflow);
            }
            if let Some(ref syftbox) = bp.syftbox {
                println!("  SyftBox: {}", syftbox);
            }
            if let Some(ref uv) = bp.uv {
                println!("  UV: {}", uv);
            }
        }
    }

    Ok(())
}

pub async fn set_email(email: String) -> crate::error::Result<()> {
    let mut config = Config::load()?;
    let config_path = Config::get_config_path()?;

    config.email = email.clone();
    config.save(&config_path)?;

    println!("✓ Email updated to: {}", email);
    Ok(())
}

pub async fn set_syftbox(path: Option<String>) -> crate::error::Result<()> {
    let mut config = Config::load()?;
    let config_path = Config::get_config_path()?;

    if let Some(path) = path {
        let syftbox_path = std::path::PathBuf::from(&path);
        if !syftbox_path.exists() {
            return Err(Error::SyftBoxConfigMissing(
                syftbox_path.display().to_string(),
            ));
        }

        let content = std::fs::read_to_string(&syftbox_path).with_context(|| {
            format!("Failed to read SyftBox config: {}", syftbox_path.display())
        })?;

        let _: serde_json::Value = serde_json::from_str(&content).with_context(|| {
            format!("Invalid JSON in SyftBox config: {}", syftbox_path.display())
        })?;

        config.syftbox_config = Some(path.clone());
        config.save(&config_path)?;

        println!("✓ SyftBox config path updated to: {}", path);
    } else {
        config.syftbox_config = None;
        config.save(&config_path)?;

        println!(
            "✓ SyftBox config path reset to default: {}",
            default_syftbox_path_str()
        );
    }

    match config.get_syftbox_data_dir() {
        Ok(data_dir) => println!("  Data directory: {}", data_dir.display()),
        Err(e) => println!("  Warning: {}", e),
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[tokio::test]
    #[serial_test::serial]
    async fn test_config_operations() -> crate::error::Result<()> {
        let temp_dir = TempDir::new()?;
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        let config_dir = temp_dir.path().join(".biovault");
        fs::create_dir_all(&config_dir)?;

        let initial_config = Config {
            email: "initial@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        initial_config.save(config_dir.join("config.yaml"))?;

        set_email("new@example.com".to_string()).await?;

        let updated_config = Config::load()?;
        assert_eq!(updated_config.email, "new@example.com");

        let syftbox_dir = temp_dir.path().join(".syftbox");
        fs::create_dir_all(&syftbox_dir)?;
        let syftbox_config_path = syftbox_dir.join("config.json");

        let syftbox_data = serde_json::json!({
            "data_dir": "/test/data"
        });
        fs::write(
            &syftbox_config_path,
            serde_json::to_string_pretty(&syftbox_data)?,
        )?;

        set_syftbox(Some(syftbox_config_path.to_string_lossy().to_string())).await?;

        let final_config = Config::load()?;
        assert!(final_config.syftbox_config.is_some());

        crate::config::clear_test_biovault_home();
        Ok(())
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_syftbox_reset_and_errors() -> crate::error::Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        fs::create_dir_all(tmp.path().join(".bv")).unwrap();
        // Seed config
        let cfg = Config {
            email: "e@example".into(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        cfg.save(tmp.path().join(".bv/config.yaml"))?;

        // set_syftbox(None) resets to default
        set_syftbox(None).await?;
        let cfg2 = Config::load()?;
        assert!(cfg2.syftbox_config.is_none());

        // Missing file path returns error
        let err = set_syftbox(Some(
            tmp.path()
                .join("no/such.json")
                .to_string_lossy()
                .to_string(),
        ))
        .await;
        assert!(err.is_err());

        // Invalid JSON returns error
        let bad = tmp.path().join("bad.json");
        fs::write(&bad, "not json").unwrap();
        let err2 = set_syftbox(Some(bad.to_string_lossy().to_string())).await;
        assert!(err2.is_err());
        crate::config::clear_test_biovault_home();
        Ok(())
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_show_config() -> crate::error::Result<()> {
        let temp_dir = TempDir::new()?;
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        let config_dir = temp_dir.path().join(".biovault");
        fs::create_dir_all(&config_dir)?;

        let config = Config {
            email: "show@example.com".to_string(),
            syftbox_config: Some("/custom/syftbox.json".to_string()),
            version: Some("1.0.0".to_string()),
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        config.save(config_dir.join("config.yaml"))?;

        // This function prints to stdout, so we just verify it doesn't panic
        let result = show(false).await;
        assert!(result.is_ok());

        crate::config::clear_test_biovault_home();
        Ok(())
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_email_validation() -> crate::error::Result<()> {
        let temp_dir = TempDir::new()?;
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        let config_dir = temp_dir.path().join(".biovault");
        fs::create_dir_all(&config_dir)?;

        // Create initial config
        let initial_config = Config {
            email: "initial@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        initial_config.save(config_dir.join("config.yaml"))?;

        // Set various email formats
        set_email("user@example.com".to_string()).await?;
        let config = Config::load()?;
        assert_eq!(config.email, "user@example.com");

        set_email("user.name+tag@example.co.uk".to_string()).await?;
        let config = Config::load()?;
        assert_eq!(config.email, "user.name+tag@example.co.uk");

        crate::config::clear_test_biovault_home();
        Ok(())
    }

    #[test]
    fn test_show_without_config() {
        // Test that show handles missing config gracefully
        let temp_dir = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        // Don't create any config file
        let runtime = tokio::runtime::Runtime::new().unwrap();
        let _result = runtime.block_on(show(false));

        // Should either succeed with defaults or fail gracefully
        // We're just making sure it doesn't panic

        crate::config::clear_test_biovault_home();
    }
}
