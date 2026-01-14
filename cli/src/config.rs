use crate::error::Error;
use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::cell::RefCell;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[cfg(test)]
use std::collections::HashMap;
use syftbox_sdk::syftbox::config::SyftboxRuntimeConfig;
use syftbox_sdk::syftbox::syc;

thread_local! {
    static TEST_CONFIG: RefCell<Option<Config>> = const { RefCell::new(None) };
    static TEST_SYFTBOX_DATA_DIR: RefCell<Option<PathBuf>> = const { RefCell::new(None) };
    static TEST_BIOVAULT_HOME: RefCell<Option<PathBuf>> = const { RefCell::new(None) };
}

#[cfg(test)]
thread_local! {
    static TEST_ENV_OVERRIDES: RefCell<HashMap<String, Vec<Option<String>>>> =
        RefCell::new(HashMap::new());
}

const PROFILES_STORE_DIR: &str = ".bvprofiles";
const PROFILES_STORE_FILE: &str = "profiles.json";

fn get_env_var(key: &str) -> Option<String> {
    #[cfg(test)]
    {
        if let Some(value) = TEST_ENV_OVERRIDES.with(|map| {
            map.borrow()
                .get(key)
                .and_then(|stack| stack.last().cloned())
        }) {
            return value;
        }
    }

    env::var(key).ok()
}

fn profiles_store_path() -> Option<PathBuf> {
    if let Ok(path) = env::var("BIOVAULT_PROFILES_PATH") {
        let trimmed = path.trim();
        if !trimmed.is_empty() {
            return Some(PathBuf::from(trimmed));
        }
    }
    if let Ok(dir) = env::var("BIOVAULT_PROFILES_DIR") {
        let trimmed = dir.trim();
        if !trimmed.is_empty() {
            let expanded = if trimmed.starts_with("~/") {
                dirs::home_dir()
                    .map(|h| h.join(&trimmed[2..]))
                    .unwrap_or_else(|| PathBuf::from(trimmed))
            } else {
                PathBuf::from(trimmed)
            };
            return Some(expanded.join(PROFILES_STORE_FILE));
        }
    }
    dirs::home_dir().map(|h| h.join(PROFILES_STORE_DIR).join(PROFILES_STORE_FILE))
}

fn read_current_profile_home() -> Option<PathBuf> {
    let store_path = profiles_store_path()?;
    let contents = fs::read_to_string(&store_path).ok()?;
    let json: serde_json::Value = serde_json::from_str(&contents).ok()?;

    let current_id = json.get("current_profile_id")?.as_str()?;
    let profiles = json.get("profiles")?.as_array()?;

    for profile in profiles {
        if profile.get("id")?.as_str()? == current_id {
            let home = profile.get("biovault_home")?.as_str()?;
            if !home.is_empty() {
                return Some(PathBuf::from(home));
            }
        }
    }
    None
}

fn detect_existing_home(home_dir: &Path) -> Option<PathBuf> {
    let desktop_root = resolve_desktop_dir(home_dir);
    let desktop_candidate = desktop_root.join("BioVault");

    // Only check Desktop/BioVault, not .biovault legacy folder
    if desktop_candidate.join("config.yaml").exists()
        || desktop_candidate.join("env").exists()
        || desktop_candidate.join("data").exists()
    {
        Some(desktop_candidate)
    } else {
        None
    }
}

fn resolve_desktop_dir(home_dir: &Path) -> PathBuf {
    if let Some(desktop_dir) = dirs::desktop_dir() {
        if desktop_dir.starts_with(home_dir) {
            return desktop_dir;
        }
    }
    home_dir.join("Desktop")
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox_credentials: Option<SyftboxCredentials>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub agent_bridge_enabled: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub agent_bridge_port: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub agent_bridge_http_port: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub agent_bridge_token: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub agent_bridge_blocklist: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct BinaryPaths {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub java: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub docker: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nextflow: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uv: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftBoxConfig {
    pub data_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SyftboxCredentials {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub refresh_token: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub access_token: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub server_url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub client_url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub data_dir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub email: Option<String>,
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
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!("Failed to create config directory: {}", parent.display())
            })?;
        }
        let yaml_content = serde_yaml::to_string(&self)?;

        if let Some(dir) = path.parent() {
            fs::create_dir_all(dir)
                .with_context(|| format!("Failed to create config directory: {}", dir.display()))?;
        }

        fs::write(path, yaml_content)
            .with_context(|| format!("Failed to write config file: {}", path.display()))?;
        Ok(())
    }

    pub fn new(email: String) -> Self {
        Self {
            email,
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        }
    }

    pub fn load_or_new(default_email: &str) -> crate::error::Result<Self> {
        match Self::load() {
            Ok(config) => Ok(config),
            Err(Error::NotInitialized) => Ok(Self::new(default_email.to_string())),
            Err(err) => Err(err),
        }
    }

    pub fn set_binary_path(
        &mut self,
        name: &str,
        path: Option<String>,
    ) -> crate::error::Result<()> {
        let paths = self.binary_paths.get_or_insert_with(BinaryPaths::default);

        match name {
            "java" => paths.java = path,
            "docker" => paths.docker = path,
            "nextflow" => paths.nextflow = path,
            "syftbox" => paths.syftbox = path,
            "uv" => paths.uv = path,
            _ => return Err(anyhow::anyhow!("Unknown dependency: {}", name).into()),
        }

        Ok(())
    }

    pub fn get_binary_path(&self, name: &str) -> Option<String> {
        self.binary_paths.as_ref().and_then(|paths| match name {
            "java" => paths.java.clone(),
            "docker" => paths.docker.clone(),
            "nextflow" => paths.nextflow.clone(),
            "syftbox" => paths.syftbox.clone(),
            "uv" => paths.uv.clone(),
            _ => None,
        })
    }

    pub fn save_binary_path(name: &str, path: Option<String>) -> crate::error::Result<Self> {
        // Try to load existing config, but preserve the email if it exists
        let mut config = match Self::load() {
            Ok(existing) => existing,
            Err(_) => Self::new(String::new()),
        };

        // Avoid persisting placeholder onboarding values when dependency paths are saved early.
        if config.email.trim() == "setup@pending" {
            config.email.clear();
        }

        let normalized = path.and_then(|p| {
            let trimmed = p.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        });

        config.set_binary_path(name, normalized)?;
        let config_path = Self::get_config_path()?;
        config.save(&config_path)?;
        Ok(config)
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

        let mut json: Value = serde_json::from_str(&content).with_context(|| {
            format!(
                "Failed to parse SyftBox config: {}",
                syftbox_config_path.display()
            )
        })?;
        let raw_data_dir = json
            .get("data_dir")
            .and_then(|v| v.as_str())
            .unwrap_or_default()
            .trim();
        let default_data_dir = Self::default_syftbox_data_dir()?;
        let legacy_default = dirs::home_dir().map(|h| h.join("SyftBox"));
        let mut dir = if raw_data_dir.is_empty() {
            default_data_dir.clone()
        } else {
            PathBuf::from(raw_data_dir)
        };
        if dir.as_os_str().is_empty() {
            dir = default_data_dir.clone();
        } else if let Some(legacy) = &legacy_default {
            if &dir == legacy {
                dir = default_data_dir.clone();
            }
        }

        // If we normalized away the legacy path, persist the new default to the config file
        if dir != Path::new(raw_data_dir) {
            if let Some(map) = json.as_object_mut() {
                map.insert(
                    "data_dir".to_string(),
                    Value::String(dir.to_string_lossy().into_owned()),
                );
                if let Ok(serialized) = serde_json::to_string_pretty(&json) {
                    if let Some(parent) = syftbox_config_path.parent() {
                        let _ = fs::create_dir_all(parent);
                    }
                    let _ = fs::write(&syftbox_config_path, serialized);
                }
            }
        }

        Ok(dir)
    }

    pub fn get_config_path() -> crate::error::Result<PathBuf> {
        Ok(get_biovault_home()?.join("config.yaml"))
    }

    pub fn get_syftbox_config_path(&self) -> crate::error::Result<PathBuf> {
        // Priority order:
        // 1. SYFTBOX_CONFIG_PATH env var
        // 2. Config specified in BioVault config
        // 3. {BIOVAULT_HOME}/syftbox/config.json (default)

        if let Ok(config_path) = env::var("SYFTBOX_CONFIG_PATH") {
            return Ok(PathBuf::from(config_path));
        }

        if let Some(ref path) = self.syftbox_config {
            return Ok(PathBuf::from(path));
        }

        Self::default_syftbox_config_path()
    }

    pub fn default_syftbox_config_path() -> crate::error::Result<PathBuf> {
        Ok(get_biovault_home()?.join("syftbox").join("config.json"))
    }

    pub fn default_syftbox_data_dir() -> crate::error::Result<PathBuf> {
        Ok(get_biovault_home()?)
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

    pub fn to_syftbox_runtime_config(&self) -> crate::error::Result<SyftboxRuntimeConfig> {
        let config_path = self.get_syftbox_config_path()?;
        let data_dir = self.get_syftbox_data_dir()?;
        let mut runtime = SyftboxRuntimeConfig::new(self.email.clone(), config_path, data_dir);
        // Prefer explicit env override for the syftbox binary, then config, else default
        if let Ok(env_bin) = std::env::var("SYFTBOX_BINARY") {
            let trimmed = env_bin.trim();
            if !trimmed.is_empty() {
                runtime = runtime.with_binary_path(Some(PathBuf::from(trimmed)));
            }
        } else if let Some(path) = self.get_binary_path("syftbox") {
            let trimmed = path.trim();
            if !trimmed.is_empty() {
                runtime = runtime.with_binary_path(Some(PathBuf::from(trimmed)));
            }
        }

        // Allow opt-out of crypto for testing/local overrides
        let disable_crypto = std::env::var("SYFTBOX_DISABLE_SYC")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        runtime = runtime.with_disable_crypto(disable_crypto);

        Ok(runtime)
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
/// 6. ~/Desktop/BioVault (default for desktop users)
///
/// Creates the directory if it doesn't exist.
pub fn get_biovault_home() -> anyhow::Result<PathBuf> {
    // Thread-local test override first
    if let Some(home) = TEST_BIOVAULT_HOME.with(|h| h.borrow().clone()) {
        // Ensure test directory exists
        fs::create_dir_all(&home).with_context(|| {
            format!(
                "Failed to create test biovault directory: {}",
                home.display()
            )
        })?;
        return Ok(home);
    }

    // Check for test environment (must be before directory walk to prevent finding global config)
    if let Some(test_home) = get_env_var("BIOVAULT_TEST_HOME") {
        // Use the path directly without appending .biovault
        let path = PathBuf::from(test_home);
        fs::create_dir_all(&path).with_context(|| {
            format!(
                "Failed to create test biovault directory: {}",
                path.display()
            )
        })?;
        return Ok(path);
    }

    // Check for explicit BIOVAULT_HOME
    if let Some(biovault_home) = get_env_var("BIOVAULT_HOME") {
        let path = PathBuf::from(biovault_home);
        fs::create_dir_all(&path)
            .with_context(|| format!("Failed to create biovault directory: {}", path.display()))?;
        return Ok(path);
    }

    // If the user has already picked a home via profiles system, respect it even
    // inside a SyftBox virtualenv to avoid creating a second .biovault.
    if let Some(profile_home) = read_current_profile_home() {
        fs::create_dir_all(&profile_home).with_context(|| {
            format!(
                "Failed to create profile biovault directory: {}",
                profile_home.display()
            )
        })?;
        return Ok(profile_home);
    }

    // Check for SyftBox virtualenv
    if let Some(syftbox_data_dir) = get_env_var("SYFTBOX_DATA_DIR") {
        let path = PathBuf::from(syftbox_data_dir).join(".biovault");
        fs::create_dir_all(&path)
            .with_context(|| format!("Failed to create biovault directory: {}", path.display()))?;
        return Ok(path);
    }

    // Check if config.yaml exists in current directory (convenient for CLI usage in dev environments)
    // This is checked AFTER SYFTBOX_DATA_DIR so sandbox tests work correctly
    if let Ok(current_dir) = env::current_dir() {
        if current_dir.join("config.yaml").exists() {
            return Ok(current_dir);
        }
    }

    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;

    // Check Desktop/BioVault BEFORE walking up directory tree
    // This ensures desktop app uses Desktop/BioVault instead of finding legacy .biovault
    if let Some(existing) = detect_existing_home(&home_dir) {
        fs::create_dir_all(&existing).with_context(|| {
            format!(
                "Failed to create biovault directory: {}",
                existing.display()
            )
        })?;
        return Ok(existing);
    }

    let desktop_dir = resolve_desktop_dir(&home_dir);
    let desktop_path = desktop_dir.join("BioVault");
    if desktop_path.exists() {
        // Desktop/BioVault folder exists, use it
        fs::create_dir_all(&desktop_path).with_context(|| {
            format!(
                "Failed to create biovault directory: {}",
                desktop_path.display()
            )
        })?;
        return Ok(desktop_path);
    }

    // Walk up from current directory looking for .biovault/config.yaml (for sbenv)
    // Only do this if Desktop/BioVault doesn't exist, to avoid finding legacy .biovault
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

    // Create Desktop/BioVault as final default
    fs::create_dir_all(&desktop_path).with_context(|| {
        format!(
            "Failed to create biovault directory: {}",
            desktop_path.display()
        )
    })?;

    Ok(desktop_path)
}

/// Get the shared cache directory
/// Priority order:
/// 1. BIOVAULT_CACHE_DIR env var
/// 2. {BioVault home}/data/cache (default shared, e.g., ~/Desktop/BioVault/data/cache)
pub fn get_cache_dir() -> anyhow::Result<PathBuf> {
    // Check for explicit cache directory
    if let Ok(cache_dir) = env::var("BIOVAULT_CACHE_DIR") {
        return Ok(PathBuf::from(cache_dir));
    }

    // Always use the shared cache within the BioVault home directory
    Ok(get_biovault_home()?.join("data").join("cache"))
}

/// Check if running in a SyftBox virtualenv
pub fn is_syftbox_env() -> bool {
    env::var("SYFTBOX_DATA_DIR").is_ok() || env::var("SYFTBOX_EMAIL").is_ok()
}

/// Resolve the default Syft Crypto vault path, preferring the BioVault home location
/// but falling back to the legacy global ~/.syc when present.
pub fn resolve_default_syc_vault_path() -> anyhow::Result<PathBuf> {
    let home = get_biovault_home()?;
    let colocated = syc::vault_path_for_home(&home);
    if colocated.exists() {
        return Ok(colocated);
    }

    let legacy = dirs::home_dir()
        .map(|h| h.join(".syc"))
        .unwrap_or_else(|| PathBuf::from(".syc"));
    if legacy.exists() {
        return Ok(legacy);
    }

    Ok(colocated)
}

/// Ensure SYC_VAULT is set so SyftBox uses the same vault as BioVault.
pub fn ensure_syc_vault_env() -> anyhow::Result<PathBuf> {
    if let Ok(current) = env::var("SYC_VAULT") {
        let path = PathBuf::from(current);
        if !path.as_os_str().is_empty() {
            return Ok(path);
        }
    }

    let vault_path = resolve_default_syc_vault_path()?;
    env::set_var("SYC_VAULT", &vault_path);
    Ok(vault_path)
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
pub fn set_test_env_override(key: &str, value: Option<&str>) {
    TEST_ENV_OVERRIDES.with(|overrides| {
        overrides
            .borrow_mut()
            .entry(key.to_string())
            .or_default()
            .push(value.map(|v| v.to_string()));
    });
}

#[cfg(test)]
pub fn clear_test_env_override(key: &str) {
    TEST_ENV_OVERRIDES.with(|overrides| {
        let mut map = overrides.borrow_mut();
        if let Some(stack) = map.get_mut(key) {
            stack.pop();
            if stack.is_empty() {
                map.remove(key);
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    struct EnvGuard {
        key: String,
        previous: Option<String>,
    }

    impl EnvGuard {
        fn set(key: &str, value: &str) -> Self {
            let previous = env::var(key).ok();
            env::set_var(key, value);
            set_test_env_override(key, Some(value));
            Self {
                key: key.to_string(),
                previous,
            }
        }

        fn unset(key: &str) -> Self {
            let previous = env::var(key).ok();
            env::remove_var(key);
            set_test_env_override(key, None);
            Self {
                key: key.to_string(),
                previous,
            }
        }
    }

    impl Drop for EnvGuard {
        fn drop(&mut self) {
            if let Some(ref value) = self.previous {
                env::set_var(&self.key, value);
            } else {
                env::remove_var(&self.key);
            }
            clear_test_env_override(&self.key);
        }
    }

    #[test]
    fn config_save_and_load_round_trip() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("config.yaml");
        let cfg = Config {
            email: "user@example.com".into(),
            syftbox_config: None,
            version: Some("1.0.0".into()),
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
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
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
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

    #[test]
    #[serial_test::serial]
    fn profiles_store_home_is_respected() {
        let tmp = TempDir::new().unwrap();
        let home_dir = tmp.path().join("home");
        let profiles_dir = tmp.path().join("profiles");
        fs::create_dir_all(&home_dir).unwrap();
        fs::create_dir_all(&profiles_dir).unwrap();
        fs::create_dir_all(home_dir.join("Desktop")).unwrap();

        let profiles_path = profiles_dir.join("profiles.json");
        let custom_home = home_dir.join("custom_location");
        fs::create_dir_all(&custom_home).unwrap();

        // Create a profiles.json with a current profile pointing to custom_home
        let profiles_json = serde_json::json!({
            "version": 1,
            "current_profile_id": "test-profile-id",
            "profiles": [{
                "id": "test-profile-id",
                "biovault_home": custom_home.to_string_lossy(),
                "created_at": "2025-01-01T00:00:00Z"
            }]
        });
        fs::write(
            &profiles_path,
            serde_json::to_string_pretty(&profiles_json).unwrap(),
        )
        .unwrap();

        let _profiles_guard = EnvGuard::set(
            "BIOVAULT_PROFILES_PATH",
            profiles_path.to_string_lossy().as_ref(),
        );
        let _unset_biovault_home = EnvGuard::unset("BIOVAULT_HOME");
        let _unset_syftbox_dir = EnvGuard::unset("SYFTBOX_DATA_DIR");
        let _unset_test_home = EnvGuard::unset("BIOVAULT_TEST_HOME");
        clear_test_biovault_home();

        let resolved_home = get_biovault_home().unwrap();
        assert_eq!(resolved_home, custom_home);
        assert!(resolved_home.exists());
    }

    #[test]
    #[serial_test::serial]
    #[cfg(not(windows))]
    fn default_home_creates_desktop_biovault() {
        let tmp = TempDir::new().unwrap();
        let home_dir = tmp.path().join("home");
        let profiles_dir = tmp.path().join("profiles");
        fs::create_dir_all(home_dir.join("Desktop")).unwrap();
        fs::create_dir_all(&profiles_dir).unwrap();

        let home_guard = EnvGuard::set("HOME", home_dir.to_string_lossy().as_ref());
        // Point to a non-existent profiles file so it falls through to default
        let profiles_path = profiles_dir.join("profiles.json");
        let _profiles_guard = EnvGuard::set(
            "BIOVAULT_PROFILES_PATH",
            profiles_path.to_string_lossy().as_ref(),
        );
        let _unset_biovault_home = EnvGuard::unset("BIOVAULT_HOME");
        let _unset_syftbox_dir = EnvGuard::unset("SYFTBOX_DATA_DIR");
        let _unset_test_home = EnvGuard::unset("BIOVAULT_TEST_HOME");
        clear_test_biovault_home();
        clear_test_syftbox_data_dir();

        let expected_home = home_dir.join("Desktop").join("BioVault");
        let _ = fs::remove_dir_all(&expected_home);

        let resolved = get_biovault_home().unwrap();
        assert_eq!(resolved, expected_home);
        assert!(resolved.exists());

        drop(home_guard);
    }
}
