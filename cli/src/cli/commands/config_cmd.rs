use crate::config::Config;
use crate::error::Error;
use anyhow::Context;

pub async fn show() -> crate::error::Result<()> {
    let config = Config::load()?;
    let config_path = Config::get_config_path()?;

    println!("BioVault Configuration");
    println!("----------------------");
    println!("Config file: {}", config_path.display());
    println!("Email: {}", config.email);

    if let Some(ref syftbox_path) = config.syftbox_config {
        println!("SyftBox config: {}", syftbox_path);
    } else {
        println!("SyftBox config: ~/.syftbox/config.json (default)");
    }

    match config.get_syftbox_data_dir() {
        Ok(data_dir) => println!("SyftBox data dir: {}", data_dir.display()),
        Err(e) => println!("SyftBox data dir: <error: {}>", e),
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

        println!("✓ SyftBox config path reset to default: ~/.syftbox/config.json");
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
    async fn test_config_operations() -> crate::error::Result<()> {
        let temp_dir = TempDir::new()?;
        std::env::set_var("BIOVAULT_TEST_HOME", temp_dir.path());

        let config_dir = temp_dir.path().join(".biovault");
        fs::create_dir_all(&config_dir)?;

        let initial_config = Config {
            email: "initial@example.com".to_string(),
            syftbox_config: None,
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

        std::env::remove_var("BIOVAULT_TEST_HOME");
        Ok(())
    }
}
