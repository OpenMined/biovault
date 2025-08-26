use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub email: String,
}

impl Config {
    pub fn new(email: String) -> Self {
        Self {
            email
        }
    }

    pub fn config_dir() -> Result<PathBuf> {
        let home = dirs::home_dir()
            .ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
        Ok(home.join(".biovault"))
    }

    pub fn config_path() -> Result<PathBuf> {
        Ok(Self::config_dir()?.join("config.yaml"))
    }

    pub fn load() -> Result<Self> {
        let config_path = Self::config_path()?;
        
        if !config_path.exists() {
            return Err(anyhow::anyhow!(
                "Configuration file not found at {}. Run 'bv init <email>' to create one.",
                config_path.display()
            ).into());
        }
        
        Self::from_file(&config_path)
    }

    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;
        
        let config: Config = serde_yaml::from_str(&content)
            .with_context(|| format!("Failed to parse config file: {}", path.display()))?;
        
        Ok(config)
    }

    pub fn save(&self) -> Result<()> {
        let config_path = Self::config_path()?;
        self.save_to(&config_path)
    }

    pub fn save_to<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create config directory: {}", parent.display()))?;
        }
        
        let yaml_content = serde_yaml::to_string(self)
            .context("Failed to serialize configuration")?;
        
        fs::write(path, yaml_content)
            .with_context(|| format!("Failed to write config file: {}", path.display()))?;
        
        Ok(())
    }

    pub fn merge_from_file<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let override_config = Self::from_file(path)?;
        self.merge(override_config);
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
            },
            _ => Err(anyhow::anyhow!("Unknown config key: {}", key).into()),
        }
    }

    pub fn display(&self) {
        println!("Email: {}", self.email);
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