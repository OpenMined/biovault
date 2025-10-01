use crate::error::Error;
use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

thread_local! {
    static TEST_CONFIG: RefCell<Option<Config>> = const { RefCell::new(None) };
    static TEST_SYFTBOX_DATA_DIR: RefCell<Option<PathBuf>> = const { RefCell::new(None) };
    static TEST_BIOVAULT_HOME: RefCell<Option<PathBuf>> = const { RefCell::new(None) };
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub email: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox_config: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binary_paths: Option<BinaryPaths>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinaryPaths {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub java: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub docker: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nextflow: Option<String>,
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
        // Thread-local test override takes highest priority
        if let Some(dir) = TEST_SYFTBOX_DATA_DIR.with(|p| p.borrow().clone()) {
            return Ok(dir);
        }

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

    pub fn get_biovault_dir(&self) -> crate::error::Result<PathBuf> {
        Ok(get_biovault_home()?)
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
    // Thread-local test override first
    if let Some(cfg) = TEST_CONFIG.with(|c| c.borrow().clone()) {
        return Ok(cfg);
    }

    // Use the Config method which respects BIOVAULT_TEST_HOME
    match Config::load() {
        Ok(config) => Ok(config),
        Err(e) => Err(anyhow::anyhow!("{}", e)),
    }
}

/// Get the BioVault home directory
/// Priority order:
/// 1. Thread-local test override (internal)
/// 2. BIOVAULT_TEST_HOME (for tests)
/// 3. BIOVAULT_HOME env var
/// 4. SYFTBOX_DATA_DIR/.biovault (if in virtualenv)
/// 5. Walk up from cwd looking for .biovault/config.yaml (for sbenv)
/// 6. ~/.biovault (default)
pub fn get_biovault_home() -> anyhow::Result<PathBuf> {
    // Thread-local test override first
    if let Some(home) = TEST_BIOVAULT_HOME.with(|h| h.borrow().clone()) {
        return Ok(home);
    }

    // Check for test environment (must be before directory walk to prevent finding global config)
    if let Ok(test_home) = env::var("BIOVAULT_TEST_HOME") {
        return Ok(PathBuf::from(test_home).join(".biovault"));
    }

    // Check for explicit BIOVAULT_HOME
    if let Ok(biovault_home) = env::var("BIOVAULT_HOME") {
        return Ok(PathBuf::from(biovault_home));
    }

    // Check for SyftBox virtualenv
    if let Ok(syftbox_data_dir) = env::var("SYFTBOX_DATA_DIR") {
        return Ok(PathBuf::from(syftbox_data_dir).join(".biovault"));
    }

    // Walk up from current directory looking for .biovault/config.yaml (for sbenv)
    if let Ok(current_dir) = env::current_dir() {
        let mut dir = current_dir.as_path();
        loop {
            let biovault_dir = dir.join(".biovault");
            let config_file = biovault_dir.join("config.yaml");
            if config_file.exists() {
                return Ok(biovault_dir);
            }

            match dir.parent() {
                Some(parent) => dir = parent,
                None => break,
            }
        }
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

/// Get the SyftBox binary path from environment
pub fn get_syftbox_binary() -> Option<String> {
    env::var("SYFTBOX_BINARY").ok()
}

/// Get the SyftBox version from environment
pub fn get_syftbox_version() -> Option<String> {
    env::var("SYFTBOX_VERSION").ok()
}

// Test utilities to isolate config and paths per test thread.
// Safe to call from tests; no-ops in normal usage if unset.
pub fn set_test_config(config: Config) {
    TEST_CONFIG.with(|c| {
        *c.borrow_mut() = Some(config);
    });
}

pub fn clear_test_config() {
    TEST_CONFIG.with(|c| {
        *c.borrow_mut() = None;
    });
}

pub fn set_test_syftbox_data_dir<P: Into<PathBuf>>(path: P) {
    TEST_SYFTBOX_DATA_DIR.with(|p| {
        *p.borrow_mut() = Some(path.into());
    });
}

pub fn clear_test_syftbox_data_dir() {
    TEST_SYFTBOX_DATA_DIR.with(|p| {
        *p.borrow_mut() = None;
    });
}

pub fn set_test_biovault_home<P: Into<PathBuf>>(path: P) {
    TEST_BIOVAULT_HOME.with(|h| {
        *h.borrow_mut() = Some(path.into());
    });
}

pub fn clear_test_biovault_home() {
    TEST_BIOVAULT_HOME.with(|h| {
        *h.borrow_mut() = None;
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn config_save_and_load_round_trip() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("config.yaml");
        let cfg = Config {
            email: "user@example.com".into(),
            syftbox_config: None,
            version: Some("1.0.0".into()),

            binary_paths: None,
        };
        cfg.save(&path).unwrap();
        let loaded = Config::from_file(&path).unwrap();
        assert_eq!(loaded.email, "user@example.com");
        assert_eq!(loaded.version.as_deref(), Some("1.0.0"));
    }

    #[test]
    fn config_syftbox_dir_from_json_and_paths() {
        let tmp = TempDir::new().unwrap();
        let data_dir = tmp.path().join("syftbox_data");
        fs::create_dir_all(&data_dir).unwrap();

        let syft_cfg = SyftBoxConfig {
            data_dir: data_dir.to_string_lossy().to_string(),
        };
        let syft_cfg_path = tmp.path().join("syftbox_config.json");
        fs::write(&syft_cfg_path, serde_json::to_string(&syft_cfg).unwrap()).unwrap();

        let cfg = Config {
            email: "user@example.com".into(),
            syftbox_config: Some(syft_cfg_path.to_string_lossy().to_string()),
            version: None,

            binary_paths: None,
        };

        let dir = cfg.get_syftbox_data_dir().unwrap();
        assert!(dir.ends_with("syftbox_data"));

        let fake = tmp.path().join("override");
        fs::create_dir_all(&fake).unwrap();
        set_test_syftbox_data_dir(&fake);
        let dir2 = cfg.get_syftbox_data_dir().unwrap();
        assert_eq!(dir2, fake);
        clear_test_syftbox_data_dir();

        let datasite = cfg.get_datasite_path().unwrap();
        assert!(datasite.ends_with(PathBuf::from("datasites/user@example.com")));
        let shared = cfg.get_shared_submissions_path().unwrap();
        assert!(shared.ends_with(PathBuf::from(
            "datasites/user@example.com/shared/biovault/submissions"
        )));
    }

    #[test]
    fn config_biovault_home_threadlocal_override() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join("home");
        fs::create_dir_all(&home).unwrap();
        set_test_biovault_home(&home);
        let got = get_biovault_home().unwrap();
        assert_eq!(got, home);
        clear_test_biovault_home();
    }

    #[test]
    fn config_get_config_path_uses_biovault_home() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join("bv_home");
        fs::create_dir_all(&home).unwrap();
        set_test_biovault_home(&home);
        let p = Config::get_config_path().unwrap();
        assert!(p.ends_with("config.yaml"));
        assert!(p.starts_with(&home));
        clear_test_biovault_home();
    }
}
