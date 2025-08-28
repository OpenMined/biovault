use crate::Result;
use anyhow::{anyhow, Context};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub email: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox_config: Option<String>,
}

impl Config {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;

        let config: Config = serde_yaml::from_str(&content)
            .with_context(|| format!("Failed to parse config file: {}", path.display()))?;

        Ok(config)
    }
}

pub fn get_config_path() -> Result<PathBuf> {
    let home = dirs::home_dir().ok_or_else(|| anyhow!("Could not determine home directory"))?;
    Ok(home.join(".biovault").join("config.yaml"))
}

pub fn get_config() -> anyhow::Result<Config> {
    let config_path = get_config_path()?;
    if !config_path.exists() {
        return Err(anyhow!(
            "BioVault config not found. Please run 'bv init' first."
        ));
    }
    let config = Config::from_file(&config_path)?;
    Ok(config)
}
