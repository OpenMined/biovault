use crate::error::Error;
use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::env;
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

    pub fn save<P: AsRef<Path>>(&self, path: P) -> crate::error::Result<()> {
        let path = path.as_ref();
        let yaml_content = serde_yaml::to_string(&self)?;
        fs::write(path, yaml_content)
            .with_context(|| format!("Failed to write config file: {}", path.display()))?;
        Ok(())
    }

    pub fn get_syftbox_data_dir(&self) -> crate::error::Result<PathBuf> {
        // Check for SYFTBOX_DATA_DIR environment variable first
        if let Ok(data_dir) = env::var("SYFTBOX_DATA_DIR") {
            return Ok(PathBuf::from(data_dir));
        }

        let syftbox_config_path = self.get_syftbox_config_path()?;

        if !syftbox_config_path.exists() {
            return Err(Error::SyftBoxConfigMissing(
                syftbox_config_path.display().to_string(),
            ));
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

    pub fn get_config_path() -> crate::error::Result<PathBuf> {
        Ok(get_biovault_home()?.join("config.yaml"))
    }

    pub fn get_syftbox_config_path(&self) -> crate::error::Result<PathBuf> {
        // Priority order:
        // 1. SYFTBOX_CONFIG_PATH env var
        // 2. Config specified in BioVault config
        // 3. ~/.syftbox/config.json (default)

        if let Ok(config_path) = env::var("SYFTBOX_CONFIG_PATH") {
            return Ok(PathBuf::from(config_path));
        }

        if let Some(ref path) = self.syftbox_config {
            return Ok(PathBuf::from(path));
        }

        let home_dir = dirs::home_dir()
            .ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
        Ok(home_dir.join(".syftbox").join("config.json"))
    }

    pub fn load() -> crate::error::Result<Self> {
        let config_path = Self::get_config_path()?;
        if !config_path.exists() {
            return Err(Error::NotInitialized);
        }
        let config = Self::from_file(config_path)?;
        Ok(config)
    }

    pub fn get_datasite_path(&self) -> crate::error::Result<PathBuf> {
        let data_dir = self.get_syftbox_data_dir()?;
        Ok(data_dir.join("datasites").join(&self.email))
    }

    pub fn get_shared_submissions_path(&self) -> crate::error::Result<PathBuf> {
        let datasite_path = self.get_datasite_path()?;
        Ok(datasite_path
            .join("shared")
            .join("biovault")
            .join("submissions"))
    }
}

pub fn get_config() -> anyhow::Result<Config> {
    // Use the Config method which respects BIOVAULT_TEST_HOME
    match Config::load() {
        Ok(config) => Ok(config),
        Err(e) => Err(anyhow::anyhow!("{}", e)),
    }
}

/// Get the BioVault home directory
/// Priority order:
/// 1. BIOVAULT_HOME env var
/// 2. SYFTBOX_DATA_DIR/.biovault (if in virtualenv)
/// 3. BIOVAULT_TEST_HOME (for tests)
/// 4. ~/.biovault (default)
pub fn get_biovault_home() -> anyhow::Result<PathBuf> {
    // Check for explicit BIOVAULT_HOME
    if let Ok(biovault_home) = env::var("BIOVAULT_HOME") {
        return Ok(PathBuf::from(biovault_home));
    }

    // Check for SyftBox virtualenv
    if let Ok(syftbox_data_dir) = env::var("SYFTBOX_DATA_DIR") {
        return Ok(PathBuf::from(syftbox_data_dir).join(".biovault"));
    }

    // Check for test environment
    if let Ok(test_home) = env::var("BIOVAULT_TEST_HOME") {
        return Ok(PathBuf::from(test_home).join(".biovault"));
    }

    // Default to home directory
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
    Ok(home_dir.join(".biovault"))
}

/// Get the shared cache directory
/// Priority order:
/// 1. BIOVAULT_CACHE_DIR env var
/// 2. ~/.biovault/data/cache (default shared)
pub fn get_cache_dir() -> anyhow::Result<PathBuf> {
    // Check for explicit cache directory
    if let Ok(cache_dir) = env::var("BIOVAULT_CACHE_DIR") {
        return Ok(PathBuf::from(cache_dir));
    }

    // Always use the shared cache in user's home directory
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
    Ok(home_dir.join(".biovault").join("data").join("cache"))
}

/// Check if running in a SyftBox virtualenv
pub fn is_syftbox_env() -> bool {
    env::var("SYFTBOX_DATA_DIR").is_ok() || env::var("SYFTBOX_EMAIL").is_ok()
}
