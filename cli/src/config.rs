use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub email: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox_config: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftBoxConfig {
    pub data_dir: String,
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

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let yaml_content = serde_yaml::to_string(&self)?;
        fs::write(path, yaml_content)
            .with_context(|| format!("Failed to write config file: {}", path.display()))?;
        Ok(())
    }

    pub fn get_syftbox_data_dir(&self) -> Result<PathBuf> {
        let syftbox_config_path = if let Some(ref path) = self.syftbox_config {
            PathBuf::from(path)
        } else {
            let home_dir = dirs::home_dir()
                .ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
            home_dir.join(".syftbox").join("config.json")
        };

        if !syftbox_config_path.exists() {
            return Err(anyhow::anyhow!(
                "SyftBox config file not found at: {}",
                syftbox_config_path.display()
            )
            .into());
        }

        let content = fs::read_to_string(&syftbox_config_path).with_context(|| {
            format!(
                "Failed to read SyftBox config: {}",
                syftbox_config_path.display()
            )
        })?;

        let syftbox_config: SyftBoxConfig = serde_json::from_str(&content).with_context(|| {
            format!(
                "Failed to parse SyftBox config: {}",
                syftbox_config_path.display()
            )
        })?;

        Ok(PathBuf::from(syftbox_config.data_dir))
    }

    pub fn get_config_path() -> Result<PathBuf> {
        let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
            PathBuf::from(test_home)
        } else {
            dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
        };
        Ok(home_dir.join(".biovault").join("config.yaml"))
    }

    pub fn load() -> Result<Self> {
        let config_path = Self::get_config_path()?;
        if !config_path.exists() {
            return Err(
                anyhow::anyhow!("BioVault not initialized. Run 'bv init <email>' first.").into(),
            );
        }
        Self::from_file(config_path)
    }
}
