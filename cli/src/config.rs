use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub email: String,
}

impl Config {
    pub fn new(email: String) -> Self {
        Self { email }
    }

    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;

        let config: Config = serde_yaml::from_str(&content)
            .with_context(|| format!("Failed to parse config file: {}", path.display()))?;

        Ok(config)
    }

    pub fn save_to<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();

        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!("Failed to create config directory: {}", parent.display())
            })?;
        }

        let yaml_content =
            serde_yaml::to_string(self).context("Failed to serialize configuration")?;

        fs::write(path, yaml_content)
            .with_context(|| format!("Failed to write config file: {}", path.display()))?;

        Ok(())
    }

    pub fn merge(&mut self, other: Config) {
        if !other.email.is_empty() {
            self.email = other.email;
        }
    }

    pub fn get(&self, key: &str) -> Result<String> {
        let parts: Vec<&str> = key.split('.').collect();

        match parts[0] {
            "email" => Ok(self.email.clone()),
            _ => Err(anyhow::anyhow!("Unknown config key: {}", key).into()),
        }
    }

    pub fn set(&mut self, key: &str, value: &str) -> Result<()> {
        let parts: Vec<&str> = key.split('.').collect();

        match parts[0] {
            "email" => {
                self.email = value.to_string();
                Ok(())
            }
            _ => Err(anyhow::anyhow!("Unknown config key: {}", key).into()),
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.email.is_empty() {
            return Err(anyhow::anyhow!("Email is required").into());
        }

        if !self.email.contains('@') {
            return Err(anyhow::anyhow!("Invalid email format: {}", self.email).into());
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_config_new() {
        let config = Config::new("test@example.com".to_string());
        assert_eq!(config.email, "test@example.com");
    }

    #[test]
    fn test_config_validate_valid_email() {
        let config = Config::new("test@example.com".to_string());
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_config_validate_empty_email() {
        let config = Config::new("".to_string());
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_config_validate_invalid_email() {
        let config = Config::new("invalid-email".to_string());
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_config_save_and_load() {
        let temp_dir = TempDir::new().unwrap();
        let config_file = temp_dir.path().join("config.yaml");

        let config = Config::new("test@example.com".to_string());
        config.save_to(&config_file).unwrap();

        assert!(config_file.exists());

        let loaded_config = Config::from_file(&config_file).unwrap();
        assert_eq!(loaded_config.email, "test@example.com");
    }

    #[test]
    fn test_config_get() {
        let config = Config::new("test@example.com".to_string());
        assert_eq!(config.get("email").unwrap(), "test@example.com");
        assert!(config.get("nonexistent").is_err());
    }

    #[test]
    fn test_config_set() {
        let mut config = Config::new("test@example.com".to_string());
        config.set("email", "new@example.com").unwrap();
        assert_eq!(config.email, "new@example.com");

        assert!(config.set("invalid_key", "value").is_err());
    }

    #[test]
    fn test_config_merge() {
        let mut config1 = Config::new("first@example.com".to_string());
        let config2 = Config::new("second@example.com".to_string());

        config1.merge(config2);
        assert_eq!(config1.email, "second@example.com");
    }

    #[test]
    fn test_config_merge_empty_email() {
        let mut config1 = Config::new("first@example.com".to_string());
        let config2 = Config::new("".to_string());

        config1.merge(config2);
        assert_eq!(config1.email, "first@example.com");
    }
}
