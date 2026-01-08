use crate::config::Config;
#[cfg(target_os = "macos")]
use crate::error::Error as CliError;
use crate::Result;
use anyhow::anyhow;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
#[cfg(target_os = "macos")]
use std::io::Write;
#[cfg(target_os = "macos")]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;
#[cfg(target_os = "macos")]
use std::sync::{Arc, OnceLock};
#[cfg(target_os = "macos")]
use tempfile::NamedTempFile;

#[cfg(target_os = "windows")]
fn configure_child_process(cmd: &mut Command) {
    use std::os::windows::process::CommandExt;
    const CREATE_NO_WINDOW: u32 = 0x08000000;
    cmd.creation_flags(CREATE_NO_WINDOW);
}

#[cfg(not(target_os = "windows"))]
fn configure_child_process(_cmd: &mut Command) {}

fn syftbox_backend_is_embedded() -> bool {
    env::var("BV_SYFTBOX_BACKEND")
        .ok()
        .map(|v| v.eq_ignore_ascii_case("embedded"))
        .unwrap_or(false)
}

#[cfg(target_os = "macos")]
type HomebrewLogger = Arc<dyn Fn(&str) + Send + Sync>;

#[cfg(target_os = "macos")]
static HOMEBREW_INSTALL_LOGGER: OnceLock<HomebrewLogger> = OnceLock::new();

#[cfg(target_os = "macos")]
const BREW_COMMON_PATHS: [&str; 3] = [
    "/opt/homebrew/bin/brew",
    "/usr/local/bin/brew",
    "/home/linuxbrew/.linuxbrew/bin/brew",
];

#[cfg(target_os = "macos")]
const OPENJDK_OPT_BASES: [&str; 3] = [
    "/opt/homebrew/opt",
    "/usr/local/opt",
    "/home/linuxbrew/.linuxbrew/opt",
];

#[cfg(target_os = "macos")]
fn resolve_brew_command() -> Option<String> {
    if let Ok(path) = which::which("brew") {
        return Some(path.display().to_string());
    }

    BREW_COMMON_PATHS
        .iter()
        .map(Path::new)
        .find(|candidate| candidate.exists())
        .map(|path| path.display().to_string())
}

#[cfg(target_os = "macos")]
fn java_bin_from_prefix(prefix: &str) -> Option<String> {
    let prefix_path = Path::new(prefix);
    let direct = prefix_path.join("bin/java");
    if direct.exists() {
        return Some(prefix_path.join("bin").display().to_string());
    }

    let alt = prefix_path
        .join("libexec")
        .join("openjdk.jdk")
        .join("Contents")
        .join("Home")
        .join("bin/java");
    if alt.exists() {
        let bin_dir = alt.parent()?;
        return Some(bin_dir.display().to_string());
    }

    None
}

fn is_valid_java_binary(path: &str) -> bool {
    let mut cmd = Command::new(path);
    cmd.arg("-version");
    configure_child_process(&mut cmd);

    cmd.output()
        .map(|output| {
            if !output.status.success() {
                return false;
            }
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            let combined = format!("{}\n{}", stderr, stdout);
            combined.contains("version")
                && (combined.contains("java") || combined.contains("openjdk"))
                && !combined.contains("Unable to locate a Java Runtime")
        })
        .unwrap_or(false)
}

fn bundled_env_key(dep: &str) -> Option<&'static str> {
    match dep {
        "java" => Some("BIOVAULT_BUNDLED_JAVA"),
        "nextflow" => Some("BIOVAULT_BUNDLED_NEXTFLOW"),
        "syftbox" => Some("SYFTBOX_BINARY"),
        "uv" => Some("BIOVAULT_BUNDLED_UV"),
        _ => None,
    }
}

fn resolve_bundled_env_path(dep: &str) -> Option<String> {
    let env_key = bundled_env_key(dep)?;
    let env_bin = std::env::var(env_key).ok()?;
    let trimmed = env_bin.trim();
    if trimmed.is_empty() {
        return None;
    }
    let p = Path::new(trimmed);
    if p.exists() {
        Some(trimmed.to_string())
    } else {
        None
    }
}

fn adjust_java_binary(mut current: Option<String>) -> Option<String> {
    if let Some(ref path) = current {
        if !is_valid_java_binary(path) {
            current = None;
        }
    }

    #[cfg(target_os = "macos")]
    {
        if current.is_none() {
            if let Some(brew_bin) = check_java_in_brew_not_in_path() {
                let candidate = format!("{}/java", brew_bin);
                if is_valid_java_binary(&candidate) {
                    return Some(candidate);
                }
            }

            if let Ok(output) = Command::new("/usr/libexec/java_home").output() {
                if output.status.success() {
                    let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
                    if !path.is_empty() {
                        let candidate = format!("{}/bin/java", path);
                        if is_valid_java_binary(&candidate) {
                            return Some(candidate);
                        }
                    }
                }
            }

            let fallback_paths = [
                "/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java",
                "/Library/Java/JavaVirtualMachines/openjdk.jdk/Contents/Home/bin/java",
                "/Library/Java/JavaVirtualMachines/temurin-21.jdk/Contents/Home/bin/java",
                "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home/bin/java",
            ];

            for path in fallback_paths.iter() {
                if Path::new(path).exists() && is_valid_java_binary(path) {
                    return Some(path.to_string());
                }
            }
        }
    }

    #[cfg(target_os = "windows")]
    {
        if current.is_none() {
            if let Some(java_path) = check_java_in_windows_not_in_path() {
                if is_valid_java_binary(&java_path) {
                    return Some(java_path);
                }
            }
        }
    }

    current
}

fn version_path_from_sources<'a>(
    primary: &'a Option<String>,
    fallback: &'a Option<String>,
) -> Option<&'a str> {
    primary
        .as_ref()
        .map(|s| s.as_str())
        .or_else(|| fallback.as_ref().map(|s| s.as_str()))
}

fn get_java_version_string(path: Option<&str>) -> Option<String> {
    let command = path.unwrap_or("java");
    let mut cmd = Command::new(command);
    cmd.arg("-version");
    configure_child_process(&mut cmd);
    let output = cmd.output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stderr);
    parse_java_version(&version_str).map(|v| v.to_string())
}

fn check_java_version_with_info_at(path: Option<&str>, min_version: u32) -> (bool, Option<String>) {
    let command = path.unwrap_or("java");
    let mut cmd = Command::new(command);
    cmd.arg("-version");
    configure_child_process(&mut cmd);
    let output = cmd.output();

    match output {
        Ok(output) => {
            let version_str = String::from_utf8_lossy(&output.stderr);
            if let Some(version) = parse_java_version(&version_str) {
                (version >= min_version, Some(version.to_string()))
            } else {
                (false, None)
            }
        }
        Err(_) => (false, None),
    }
}

pub fn set_homebrew_install_logger<F>(logger: F)
where
    F: Fn(&str) + Send + Sync + 'static,
{
    #[cfg(target_os = "macos")]
    {
        let _ = HOMEBREW_INSTALL_LOGGER.set(Arc::new(logger));
    }

    #[cfg(not(target_os = "macos"))]
    {
        let _ = logger;
    }
}

#[cfg(target_os = "macos")]
fn log_homebrew_install(message: &str) {
    if let Some(callback) = HOMEBREW_INSTALL_LOGGER.get() {
        callback(message);
    } else {
        eprintln!("🍺 {message}");
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DependencyConfig {
    pub dependencies: Vec<Dependency>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Dependency {
    pub name: String,
    pub check_running: bool,
    #[serde(default)]
    pub min_version: Option<u32>,
    pub install_instructions: String,
    pub description: String,
    #[serde(default)]
    pub website: Option<String>,
    #[serde(default)]
    pub environments: Option<HashMap<String, EnvironmentConfig>>,
    #[serde(default)]
    pub package_manager: bool,
    #[serde(default)]
    pub os_filter: Option<Vec<String>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct EnvironmentConfig {
    #[serde(default)]
    pub install_commands: Option<Vec<String>>,
    #[serde(default)]
    pub verify_command: Option<String>,
    #[serde(default)]
    pub skip: bool,
    #[serde(default)]
    pub skip_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DependencyCheckResult {
    pub dependencies: Vec<DependencyResult>,
    pub all_satisfied: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DependencyResult {
    pub name: String,
    pub found: bool,
    pub path: Option<String>,
    pub version: Option<String>,
    pub running: Option<bool>,
    pub skipped: bool,
    pub skip_reason: Option<String>,
    pub description: Option<String>,
    pub website: Option<String>,
    pub install_instructions: Option<String>,
}

/// Check a single dependency with optional custom path (for library use)
pub fn check_single_dependency(
    name: &str,
    custom_path: Option<String>,
) -> Result<DependencyResult> {
    // Load the deps.yaml file embedded in the binary
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Find the dependency by name
    let dep = config
        .dependencies
        .iter()
        .find(|d| d.name == name)
        .ok_or_else(|| anyhow!("Dependency '{}' not found", name))?;

    // Detect if we're in Google Colab
    let is_colab = is_google_colab();

    // Check if this dependency should be skipped in Colab
    if is_colab {
        if let Some(environments) = &dep.environments {
            if let Some(colab_config) = environments.get("google_colab") {
                if colab_config.skip {
                    return Ok(DependencyResult {
                        name: dep.name.clone(),
                        found: false,
                        path: None,
                        version: None,
                        running: None,
                        skipped: true,
                        skip_reason: colab_config.skip_reason.clone(),
                        description: Some(dep.description.clone()),
                        website: dep.website.clone(),
                        install_instructions: Some(dep.install_instructions.clone()),
                    });
                }
            }
        }
    }

    let forced_missing: Option<HashSet<String>> = std::env::var("BIOVAULT_FORCE_MISSING_DEPS")
        .ok()
        .map(|value| {
            value
                .split(',')
                .map(|item| item.trim().to_lowercase())
                .filter(|item| !item.is_empty())
                .collect()
        });

    let should_force_missing = forced_missing.as_ref().is_some_and(|set| {
        set.contains("all") || set.contains("*") || set.contains(&dep.name.to_lowercase())
    });

    if should_force_missing {
        eprintln!(
            "🔧 Forcing dependency '{}' to appear missing (BIOVAULT_FORCE_MISSING_DEPS)",
            dep.name
        );

        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: false,
            path: custom_path.clone(),
            version: None,
            running: None,
            skipped: false,
            skip_reason: None,
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    }

    if dep.name == "syftbox" && syftbox_backend_is_embedded() {
        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: true,
            path: None,
            version: None,
            running: None,
            skipped: true,
            skip_reason: Some("Provided by BioVault (embedded)".to_string()),
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    }

    // On Windows, Nextflow and Java are run via Docker, so we treat them as satisfied (not required)
    // to avoid blocking onboarding/UI and to avoid prompting for installation.
    if std::env::consts::OS == "windows" && (dep.name == "java" || dep.name == "nextflow") {
        let reason = dep
            .environments
            .as_ref()
            .and_then(|envs| envs.get("windows"))
            .and_then(|cfg| cfg.skip_reason.clone())
            .unwrap_or_else(|| "Not required on Windows (runs via Docker)".to_string());

        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: true,
            path: None,
            version: None,
            running: None,
            skipped: true,
            skip_reason: Some(reason),
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    }

    // Check if the binary exists at custom path or in PATH
    let mut fallback_path: Option<String> = custom_path.clone();

    let direct_path = if let Some(custom) = custom_path.as_ref() {
        if Path::new(custom).exists() {
            // Verify this is actually the right tool by checking its version
            let is_valid = match dep.name.as_str() {
                "java" => is_valid_java_binary(custom),
                "docker" => {
                    let mut cmd = Command::new(custom);
                    cmd.arg("--version");
                    configure_child_process(&mut cmd);
                    cmd.output()
                        .map(|output| {
                            let version_str = String::from_utf8_lossy(&output.stdout);
                            version_str.contains("Docker")
                        })
                        .unwrap_or(false)
                }
                "nextflow" => {
                    let mut cmd = Command::new(custom);
                    cmd.arg("-version");
                    configure_child_process(&mut cmd);
                    cmd.output()
                        .map(|output| {
                            let version_str = String::from_utf8_lossy(&output.stdout);
                            version_str.contains("nextflow") || version_str.contains("version")
                        })
                        .unwrap_or(false)
                }
                "syftbox" => {
                    let mut cmd = Command::new(custom);
                    cmd.arg("--version");
                    configure_child_process(&mut cmd);
                    cmd.output()
                        .map(|output| {
                            let version_str = String::from_utf8_lossy(&output.stdout);
                            version_str.contains("syftbox") || version_str.contains("version")
                        })
                        .unwrap_or(false)
                }
                _ => true,
            };

            if is_valid {
                Some(custom.clone())
            } else {
                None
            }
        } else {
            None
        }
    } else {
        // Check config for custom path
        let bv_config = Config::load().ok();
        if let Some(config_path) = bv_config
            .as_ref()
            .and_then(|cfg| cfg.get_binary_path(&dep.name))
        {
            fallback_path = Some(config_path.clone());
            if Path::new(&config_path).exists() {
                Some(config_path)
            } else {
                None
            }
        } else {
            None
        }
    };

    let mut binary_path = direct_path
        .or_else(|| resolve_bundled_env_path(&dep.name))
        .or_else(|| {
            which::which(&dep.name)
                .ok()
                .map(|p| p.display().to_string())
        })
        .or_else(|| find_in_well_known_locations(&dep.name));

    if dep.name == "java" {
        binary_path = adjust_java_binary(binary_path);
        if binary_path.is_some() {
            fallback_path = binary_path.clone();
        }
    }

    // Special handling for Docker on macOS - check if Docker Desktop is installed
    #[cfg(target_os = "macos")]
    let docker_desktop_installed = dep.name == "docker" && is_docker_desktop_installed();
    #[cfg(not(target_os = "macos"))]
    let docker_desktop_installed = false;

    // Special handling for syftbox - check in ~/.sbenv/binaries
    let mut syftbox_sbenv_path: Option<String> = None;
    if binary_path.is_none() && dep.name == "syftbox" {
        syftbox_sbenv_path = check_syftbox_in_sbenv();
    }

    if binary_path.is_none() && !docker_desktop_installed && syftbox_sbenv_path.is_none() {
        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: false,
            path: fallback_path,
            version: None,
            running: None,
            skipped: false,
            skip_reason: None,
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    } else if docker_desktop_installed && binary_path.is_none() {
        // Docker Desktop is installed but docker CLI is not available in PATH
        // This typically happens when Docker Desktop is installed but not running
        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: true,
            path: Some("/Applications/Docker.app".to_string()),
            version: None,
            running: Some(false),
            skipped: false,
            skip_reason: None,
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some("Docker Desktop is installed but not running. Please start Docker Desktop from your Applications folder.".to_string()),
        });
    } else if let Some(sbenv_path) = syftbox_sbenv_path {
        // Syftbox found in ~/.sbenv/binaries
        // Get version from the syftbox binary
        let mut cmd = Command::new(&sbenv_path);
        cmd.arg("--version");
        configure_child_process(&mut cmd);
        let version = cmd
            .output()
            .ok()
            .filter(|output| output.status.success())
            .and_then(|output| {
                let version_str = String::from_utf8_lossy(&output.stdout);
                version_str.split_whitespace().nth(2).map(|v| v.to_string())
            });

        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: true,
            path: Some(sbenv_path),
            version,
            running: None,
            skipped: false,
            skip_reason: None,
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    } else if let Some(ref path) = binary_path {
        // Binary found - check version
        let version_source = version_path_from_sources(&binary_path, &fallback_path);
        let (version_ok, version_str) = if let Some(min_version) = dep.min_version {
            if dep.name == "java" {
                check_java_version_with_info_at(version_source, min_version)
            } else {
                check_version_with_info(&dep.name, min_version, version_source)
            }
        } else if dep.name == "java" {
            (true, get_java_version_string(version_source))
        } else {
            (true, get_version_string(&dep.name, version_source))
        };

        // Check if it needs to be running
        let is_running = if dep.check_running {
            Some(check_if_running(&dep.name))
        } else {
            None
        };

        return Ok(DependencyResult {
            name: dep.name.clone(),
            found: version_ok,
            path: Some(path.clone()),
            version: version_str,
            running: is_running,
            skipped: false,
            skip_reason: None,
            description: Some(dep.description.clone()),
            website: dep.website.clone(),
            install_instructions: Some(dep.install_instructions.clone()),
        });
    }

    // Fallback
    Ok(DependencyResult {
        name: dep.name.clone(),
        found: false,
        path: None,
        version: None,
        running: None,
        skipped: false,
        skip_reason: None,
        description: Some(dep.description.clone()),
        website: dep.website.clone(),
        install_instructions: Some(dep.install_instructions.clone()),
    })
}

/// Check dependencies and return the result (for library use)
pub fn check_dependencies_result() -> Result<DependencyCheckResult> {
    // Load the deps.yaml file embedded in the binary
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Detect if we're in Google Colab
    let is_colab = is_google_colab();

    // Load BioVault config to check for custom binary paths
    let bv_config = Config::load().ok();

    let mut all_found = true;
    let mut all_running = true;
    let mut results: Vec<DependencyResult> = Vec::new();

    // Get current OS for filtering
    let current_os = std::env::consts::OS;

    for dep in &config.dependencies {
        // Skip package managers (they're checked separately during onboarding)
        if dep.package_manager {
            continue;
        }

        // Check OS filter - skip if this dependency is not for the current OS
        if let Some(os_filter) = &dep.os_filter {
            if !os_filter.contains(&current_os.to_string()) {
                continue;
            }
        }
        // Check if this dependency should be skipped in Colab
        if is_colab {
            if let Some(environments) = &dep.environments {
                if let Some(colab_config) = environments.get("google_colab") {
                    if colab_config.skip {
                        results.push(DependencyResult {
                            name: dep.name.clone(),
                            found: false,
                            path: None,
                            version: None,
                            running: None,
                            skipped: true,
                            skip_reason: colab_config.skip_reason.clone(),
                            description: Some(dep.description.clone()),
                            website: dep.website.clone(),
                            install_instructions: Some(dep.install_instructions.clone()),
                        });
                        continue;
                    }
                }
            }
        }

        // Check for custom binary path in config
        let custom_path = bv_config
            .as_ref()
            .and_then(|cfg| cfg.get_binary_path(&dep.name));

        // On Windows, Nextflow and Java are run via Docker, so we treat them as satisfied (not required).
        if current_os == "windows" && (dep.name == "java" || dep.name == "nextflow") {
            let reason = dep
                .environments
                .as_ref()
                .and_then(|envs| envs.get("windows"))
                .and_then(|cfg| cfg.skip_reason.clone())
                .unwrap_or_else(|| "Not required on Windows (runs via Docker)".to_string());

            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: None,
                version: None,
                running: None,
                skipped: true,
                skip_reason: Some(reason),
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
            continue;
        }

        if dep.name == "syftbox" && syftbox_backend_is_embedded() {
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: None,
                version: None,
                running: None,
                skipped: true,
                skip_reason: Some("Provided by BioVault (embedded)".to_string()),
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
            continue;
        }

        // Priority for syftbox: SYFTBOX_BINARY env -> custom path -> PATH -> well-known locations -> sbenv
        let mut binary_path = None;

        if dep.name == "syftbox" {
            if let Ok(env_bin) = std::env::var("SYFTBOX_BINARY") {
                if !env_bin.trim().is_empty() && Path::new(&env_bin).exists() {
                    binary_path = Some(env_bin);
                }
            }
        }

        if binary_path.is_none() {
            if let Some(env_key) = bundled_env_key(&dep.name) {
                if let Ok(env_bin) = env::var(env_key) {
                    if !env_bin.trim().is_empty() && Path::new(&env_bin).exists() {
                        binary_path = Some(env_bin);
                    }
                }
            }
        }

        if binary_path.is_none() {
            binary_path = if custom_path
                .as_ref()
                .map(|p| Path::new(p).exists())
                .unwrap_or(false)
            {
                custom_path.clone()
            } else {
                which::which(&dep.name)
                    .ok()
                    .map(|p| p.display().to_string())
            }
            .or_else(|| find_in_well_known_locations(&dep.name));
        }
        if dep.name == "java" {
            binary_path = adjust_java_binary(binary_path);
        }

        // Special handling for Docker on macOS - check if Docker Desktop is installed
        #[cfg(target_os = "macos")]
        let docker_desktop_installed = dep.name == "docker" && is_docker_desktop_installed();
        #[cfg(not(target_os = "macos"))]
        let docker_desktop_installed = false;

        // Special handling for syftbox - check in ~/.sbenv/binaries
        let mut syftbox_sbenv_path: Option<String> = None;
        if binary_path.is_none() && dep.name == "syftbox" {
            syftbox_sbenv_path = check_syftbox_in_sbenv();
        }

        if binary_path.is_none() && !docker_desktop_installed && syftbox_sbenv_path.is_none() {
            all_found = false;
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: false,
                path: custom_path.clone(),
                version: None,
                running: None,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        } else if docker_desktop_installed && binary_path.is_none() {
            // Docker Desktop is installed but docker CLI is not available
            all_running = false; // Mark as not running since docker CLI is not available
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some("/Applications/Docker.app".to_string()),
                version: None,
                running: Some(false),
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some("Docker Desktop is installed but not running. Please start Docker Desktop from your Applications folder.".to_string()),
            });
        } else if let Some(sbenv_path) = syftbox_sbenv_path {
            // Syftbox found in ~/.sbenv/binaries
            let mut cmd = Command::new(&sbenv_path);
            cmd.arg("--version");
            configure_child_process(&mut cmd);
            let version = cmd
                .output()
                .ok()
                .filter(|output| output.status.success())
                .and_then(|output| {
                    let version_str = String::from_utf8_lossy(&output.stdout);
                    version_str.split_whitespace().nth(2).map(|v| v.to_string())
                });

            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some(sbenv_path),
                version,
                running: None,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        } else if let Some(ref path) = binary_path {
            // Binary found - check version
            let version_source = Some(path.as_str());
            let (version_ok, version_str) = if let Some(min_version) = dep.min_version {
                if dep.name == "java" {
                    check_java_version_with_info_at(version_source, min_version)
                } else {
                    check_version_with_info(&dep.name, min_version, version_source)
                }
            } else if dep.name == "java" {
                (true, get_java_version_string(version_source))
            } else {
                (true, get_version_string(&dep.name, version_source))
            };

            if !version_ok {
                all_found = false;
            }

            // Check if it needs to be running
            let is_running = if dep.check_running {
                Some(check_if_running(&dep.name))
            } else {
                None
            };

            if let Some(running) = is_running {
                if !running {
                    all_running = false;
                }
            }

            results.push(DependencyResult {
                name: dep.name.clone(),
                found: version_ok,
                path: Some(path.clone()),
                version: version_str,
                running: is_running,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        }
    }

    Ok(DependencyCheckResult {
        dependencies: results,
        all_satisfied: all_found && all_running,
    })
}

fn find_in_well_known_locations(executable: &str) -> Option<String> {
    let mut search_dirs: Vec<PathBuf> = Vec::new();

    #[cfg(target_os = "macos")]
    {
        search_dirs.extend([
            PathBuf::from("/opt/homebrew/bin"),
            PathBuf::from("/opt/homebrew/sbin"),
            PathBuf::from("/usr/local/bin"),
            PathBuf::from("/usr/local/sbin"),
            PathBuf::from("/usr/bin"),
            PathBuf::from("/usr/sbin"),
            PathBuf::from("/bin"),
            PathBuf::from("/sbin"),
        ]);

        if executable == "docker" {
            search_dirs.extend([
                PathBuf::from("/Applications/Docker.app/Contents/Resources/bin"),
                PathBuf::from("/usr/local/bin"), // Docker CLI installed via brew
                PathBuf::from("/opt/homebrew/bin"), // Docker CLI on Apple Silicon
            ]);
        }
    }

    #[cfg(target_os = "linux")]
    {
        search_dirs.extend([
            PathBuf::from("/usr/local/bin"),
            PathBuf::from("/usr/bin"),
            PathBuf::from("/bin"),
            PathBuf::from("/usr/local/sbin"),
            PathBuf::from("/usr/sbin"),
            PathBuf::from("/sbin"),
            PathBuf::from("/snap/bin"),
            PathBuf::from("/var/lib/snapd/snap/bin"),
        ]);
    }

    #[cfg(target_os = "windows")]
    {
        search_dirs.extend([
            PathBuf::from("C:/Program Files"),
            PathBuf::from("C:/Program Files (x86)"),
        ]);

        if executable == "docker" {
            search_dirs.push(PathBuf::from(
                "C:/Program Files/Docker/Docker/resources/bin",
            ));
        }
    }

    if let Some(home) = dirs::home_dir() {
        search_dirs.push(home.join(".local/bin"));
        search_dirs.push(home.join(".cargo/bin"));
        search_dirs.push(home.join("bin"));
    }

    let binary_name = if cfg!(target_os = "windows") {
        format!("{}.exe", executable)
    } else {
        executable.to_string()
    };

    for dir in search_dirs {
        let candidate = dir.join(&binary_name);
        if candidate.is_file() {
            return Some(candidate.display().to_string());
        }
    }

    None
}

pub async fn execute(json: bool) -> Result<()> {
    // Load the deps.yaml file embedded in the binary
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    if !json {
        println!("BioVault Dependency Check");
        println!("=========================\n");
    }

    // Detect if we're in Google Colab
    let is_colab = is_google_colab();
    if is_colab && !json {
        println!("ℹ️  Google Colab environment detected\n");
    }

    // Check if we're in CI mode (non-interactive)
    let is_ci = env::var("CI").is_ok() || env::var("GITHUB_ACTIONS").is_ok();

    // Load BioVault config to check for custom binary paths
    let bv_config = Config::load().ok();

    let mut all_found = true;
    let mut all_running = true;
    let mut results: Vec<DependencyResult> = Vec::new();

    // Get current OS for filtering
    let current_os = std::env::consts::OS;

    for dep in &config.dependencies {
        // Skip package managers (they're checked separately during onboarding)
        if dep.package_manager {
            continue;
        }

        // Check OS filter - skip if this dependency is not for the current OS
        if let Some(os_filter) = &dep.os_filter {
            if !os_filter.contains(&current_os.to_string()) {
                continue;
            }
        }
        // Check if this dependency should be skipped in Colab
        if is_colab {
            if let Some(environments) = &dep.environments {
                if let Some(colab_config) = environments.get("google_colab") {
                    if colab_config.skip {
                        if !json {
                            println!("Checking {}... ⏭️  SKIPPED", dep.name);
                            println!(
                                "  Reason: {}",
                                colab_config
                                    .skip_reason
                                    .as_ref()
                                    .unwrap_or(&"Not available in Colab".to_string())
                            );
                            println!();
                        }
                        results.push(DependencyResult {
                            name: dep.name.clone(),
                            found: false,
                            path: None,
                            version: None,
                            running: None,
                            skipped: true,
                            skip_reason: colab_config.skip_reason.clone(),
                            description: Some(dep.description.clone()),
                            website: dep.website.clone(),
                            install_instructions: Some(dep.install_instructions.clone()),
                        });
                        continue;
                    }
                }
            }
        }

        if !json {
            print!("Checking {}... ", dep.name);
        }

        // Check for custom binary path in config
        let custom_path = bv_config
            .as_ref()
            .and_then(|cfg| cfg.binary_paths.as_ref())
            .and_then(|bp| match dep.name.as_str() {
                "java" => bp.java.clone(),
                "docker" => bp.docker.clone(),
                "nextflow" => bp.nextflow.clone(),
                _ => None,
            });

        // Check if the binary exists in PATH or use custom path
        let mut binary_path = if let Some(custom) = custom_path {
            if std::path::Path::new(&custom).exists() {
                Some(custom)
            } else {
                None
            }
        } else {
            which::which(&dep.name)
                .ok()
                .map(|p| p.display().to_string())
        };

        if dep.name == "java" {
            binary_path = adjust_java_binary(binary_path);
        }

        // Special handling for UV on Windows (WinGet installs but PATH needs refresh)
        let mut uv_windows_path: Option<String> = None;
        if binary_path.is_none() && dep.name == "uv" && std::env::consts::OS == "windows" {
            uv_windows_path = check_uv_in_windows_not_in_path();
        }

        // Special handling for Docker on macOS - check if Docker Desktop is installed
        #[cfg(target_os = "macos")]
        let docker_desktop_installed = dep.name == "docker" && is_docker_desktop_installed();
        #[cfg(not(target_os = "macos"))]
        let docker_desktop_installed = false;

        // Special handling for syftbox - check in ~/.sbenv/binaries
        let mut syftbox_sbenv_path: Option<String> = None;
        if binary_path.is_none() && dep.name == "syftbox" {
            syftbox_sbenv_path = check_syftbox_in_sbenv();
        }

        if binary_path.is_none()
            && uv_windows_path.is_none()
            && !docker_desktop_installed
            && syftbox_sbenv_path.is_none()
        {
            all_found = false;
            if !json {
                println!("❌ NOT FOUND");
                println!("  Description: {}", dep.description);
                println!("  Installation instructions:");
                for line in dep.install_instructions.lines() {
                    if !line.trim().is_empty() {
                        println!("    {}", line);
                    }
                }
                println!();
            }
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: false,
                path: None,
                version: None,
                running: None,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        } else if docker_desktop_installed && binary_path.is_none() {
            // Docker Desktop is installed but not running
            all_running = false;
            if !json {
                println!("⚠️  Found (Docker Desktop installed but not running)");
                println!("  Docker Desktop is installed on your system.");
                println!("  To start Docker, open Docker Desktop from your Applications folder.");
                println!("  Location: /Applications/Docker.app");
                println!();
            }
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some("/Applications/Docker.app".to_string()),
                version: None,
                running: Some(false),
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some("Docker Desktop is installed but not running. Please start Docker Desktop from your Applications folder.".to_string()),
            });
        } else if let Some(uv_path) = uv_windows_path {
            // UV is installed on Windows but not in PATH
            if !json {
                println!("⚠️  Found (not in PATH)");
                println!("  UV is installed at: {}", uv_path);
                println!("  But it's not available in your PATH.");
                println!("  To make it available:");
                println!("  Option 1: Open a new Command Prompt or PowerShell window");
                println!("  Option 2: Manually refresh your PATH:");
                println!("    PowerShell: $env:Path = [System.Environment]::GetEnvironmentVariable(\"Path\",\"User\")");
                if !is_ci {
                    println!("  Option 3: Run 'bv setup' again after restarting your terminal");
                }
                println!();
            }
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some(uv_path.clone()),
                version: get_uv_version(&uv_path),
                running: None,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        } else if let Some(sbenv_path) = syftbox_sbenv_path {
            // Syftbox found in ~/.sbenv/binaries
            let mut cmd = Command::new(&sbenv_path);
            cmd.arg("--version");
            configure_child_process(&mut cmd);
            let version = cmd
                .output()
                .ok()
                .filter(|output| output.status.success())
                .and_then(|output| {
                    let version_str = String::from_utf8_lossy(&output.stdout);
                    version_str.split_whitespace().nth(2).map(|v| v.to_string())
                });

            if !json {
                if let Some(ref ver) = version {
                    println!("✓ Found (version {})", ver);
                } else {
                    println!("✓ Found");
                }
                println!("  Path: {}", sbenv_path);
                println!("  Note: Found in ~/.sbenv/binaries");
            }

            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some(sbenv_path),
                version,
                running: None,
                skipped: false,
                skip_reason: None,
                description: Some(dep.description.clone()),
                website: dep.website.clone(),
                install_instructions: Some(dep.install_instructions.clone()),
            });
        } else if let Some(ref path) = binary_path {
            // Binary found - check version
            let version_source = Some(path.as_str());
            let (version_ok, version_str) = if let Some(min_version) = dep.min_version {
                if dep.name == "java" {
                    check_java_version_with_info_at(version_source, min_version)
                } else {
                    check_version_with_info(&dep.name, min_version, version_source)
                }
            } else if dep.name == "java" {
                (true, get_java_version_string(version_source))
            } else {
                (true, get_version_string(&dep.name, version_source))
            };

            if !version_ok {
                all_found = false;
                if !json {
                    let min_ver = dep.min_version.unwrap_or(0);
                    println!("❌ Version too old (requires {} or higher)", min_ver);
                    println!("  Description: {}", dep.description);
                    println!("  Installation instructions:");
                    for line in dep.install_instructions.lines() {
                        if !line.trim().is_empty() {
                            println!("    {}", line);
                        }
                    }
                    println!();
                }
                results.push(DependencyResult {
                    name: dep.name.clone(),
                    found: true,
                    path: Some(path.clone()),
                    version: version_str,
                    running: None,
                    skipped: false,
                    skip_reason: None,
                    description: Some(dep.description.clone()),
                    website: dep.website.clone(),
                    install_instructions: Some(dep.install_instructions.clone()),
                });
            } else {
                // Check if it needs to be running
                let is_running = if dep.check_running {
                    Some(check_if_running(&dep.name))
                } else {
                    None
                };

                if let Some(running) = is_running {
                    if !running {
                        all_running = false;
                    }
                }

                if !json {
                    if let Some(ref ver) = version_str {
                        print!("✓ Found (version {})", ver);
                    } else {
                        print!("✓ Found");
                    }

                    if let Some(running) = is_running {
                        if running {
                            println!(" (running)");
                        } else {
                            println!(" (NOT RUNNING)");
                            println!(
                                "  To start {}, run: {}",
                                dep.name,
                                get_start_command(&dep.name)
                            );
                        }
                    } else {
                        println!();
                    }

                    println!("  Path: {}", path);
                }

                results.push(DependencyResult {
                    name: dep.name.clone(),
                    found: true,
                    path: Some(path.clone()),
                    version: version_str,
                    running: is_running,
                    skipped: false,
                    skip_reason: None,
                    description: Some(dep.description.clone()),
                    website: dep.website.clone(),
                    install_instructions: Some(dep.install_instructions.clone()),
                });
            }
        }
    }

    if json {
        let result = DependencyCheckResult {
            dependencies: results,
            all_satisfied: all_found && all_running,
        };
        println!("{}", serde_json::to_string_pretty(&result)?);
        if all_found && all_running {
            Ok(())
        } else if !all_found {
            Err(anyhow!("Dependencies missing").into())
        } else {
            Err(anyhow!("Services not running").into())
        }
    } else {
        println!("\n=========================");
        if all_found && all_running {
            println!("✓ All dependencies satisfied!");
            Ok(())
        } else if !all_found {
            println!(
                "⚠️  Some dependencies are missing. Please install them using the instructions above."
            );
            Err(anyhow!("Dependencies missing").into())
        } else {
            // !all_running
            println!(
                "⚠️  Some services are not running. Please start them using the commands above."
            );
            Err(anyhow!("Services not running").into())
        }
    }
}

fn check_if_running(service: &str) -> bool {
    match service {
        "docker" => {
            // Check if Docker daemon is running
            let mut cmd = Command::new("docker");
            cmd.arg("info");
            configure_child_process(&mut cmd);
            cmd.output()
                .map(|output| output.status.success())
                .unwrap_or(false)
        }
        _ => false,
    }
}

#[cfg(target_os = "macos")]
fn is_docker_desktop_installed() -> bool {
    // Check if Docker Desktop app exists on macOS
    let docker_app = Path::new("/Applications/Docker.app");
    if docker_app.exists() {
        return true;
    }

    // Also check user-specific Applications folder
    if let Some(home) = dirs::home_dir() {
        let user_docker_app = home.join("Applications/Docker.app");
        if user_docker_app.exists() {
            return true;
        }
    }

    false
}

#[cfg(not(target_os = "macos"))]
#[allow(dead_code)]
fn is_docker_desktop_installed() -> bool {
    // On non-macOS systems, we rely on the docker binary check
    false
}

fn check_syftbox_in_sbenv() -> Option<String> {
    // Check for syftbox binaries in ~/.sbenv/binaries
    let home = dirs::home_dir()?;
    let sbenv_binaries = home.join(".sbenv").join("binaries");

    if !sbenv_binaries.exists() {
        return None;
    }

    // Find all syftbox binaries and pick the latest one
    let mut syftbox_binaries = Vec::new();

    // Look in version subdirectories
    if let Ok(entries) = std::fs::read_dir(&sbenv_binaries) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                // Check for syftbox binary in each version directory
                let syftbox_path = path.join("syftbox");
                if syftbox_path.is_file() {
                    // Check if it's executable on Unix systems
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::PermissionsExt;
                        if let Ok(metadata) = syftbox_path.metadata() {
                            let permissions = metadata.permissions();
                            if permissions.mode() & 0o111 != 0 {
                                syftbox_binaries.push(syftbox_path);
                            }
                        }
                    }
                    #[cfg(not(unix))]
                    {
                        syftbox_binaries.push(syftbox_path);
                    }
                }
            }
        }
    }

    if syftbox_binaries.is_empty() {
        return None;
    }

    // Sort by path (which includes version) to get the latest version
    // Versions like 0.8.6 > 0.8.5 > 0.8.3
    syftbox_binaries.sort_by(|a, b| {
        // Extract version from path
        let a_parent = a
            .parent()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy());
        let b_parent = b
            .parent()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy());

        // Compare versions - this will work for simple semantic versions
        b_parent.cmp(&a_parent)
    });

    syftbox_binaries.first().map(|p| p.display().to_string())
}

fn get_start_command(service: &str) -> String {
    match service {
        "docker" => {
            // Provide OS-specific instructions
            if std::env::consts::OS == "macos" {
                "Open Docker Desktop from your Applications folder (/Applications/Docker.app)"
                    .to_string()
            } else if std::env::consts::OS == "windows" {
                "Open Docker Desktop from your Start menu".to_string()
            } else {
                "Start Docker daemon with 'sudo systemctl start docker' or open Docker Desktop"
                    .to_string()
            }
        }
        _ => format!("Start {}", service),
    }
}

#[allow(dead_code)]
fn check_version(tool: &str, min_version: u32) -> bool {
    match tool {
        "java" => check_java_version(min_version),
        _ => true, // For tools without version checking, assume OK
    }
}

#[allow(dead_code)]
fn check_java_version(min_version: u32) -> bool {
    // Try to execute java -version
    let mut command = Command::new("java");
    command.arg("-version");
    configure_child_process(&mut command);
    let output = command.output();

    match output {
        Ok(output) => {
            // Java version info is typically printed to stderr
            let version_str = String::from_utf8_lossy(&output.stderr);

            // Parse the version from the output
            if let Some(version) = parse_java_version(&version_str) {
                if version >= min_version {
                    print!(" (version {})", version);
                    true
                } else {
                    false
                }
            } else {
                false
            }
        }
        Err(_) => false,
    }
}

fn parse_java_version(output: &str) -> Option<u32> {
    // Java version output can be in different formats:
    // - openjdk version "17.0.2" 2022-01-18
    // - java version "1.8.0_321"
    // - openjdk version "11.0.14" 2022-01-18
    // - java version "17" 2021-09-14

    for line in output.lines() {
        if line.contains("version") {
            // Extract the version string in quotes
            if let Some(start) = line.find('"') {
                if let Some(end) = line[start + 1..].find('"') {
                    let version_str = &line[start + 1..start + 1 + end];

                    // Handle "1.x" format (Java 8 and earlier)
                    if let Some(stripped) = version_str.strip_prefix("1.") {
                        // Extract the minor version (e.g., "1.8.0_321" -> 8)
                        if let Some(dot_pos) = stripped.find('.') {
                            if let Ok(version) = stripped[..dot_pos].parse::<u32>() {
                                return Some(version);
                            }
                        }
                    } else {
                        // Modern format: extract major version
                        // Handle both "17" and "17.0.2" formats
                        let major_part = version_str.split('.').next().unwrap_or(version_str);
                        if let Ok(version) = major_part.parse::<u32>() {
                            return Some(version);
                        }
                    }
                }
            }
        }
    }

    None
}

fn check_version_with_info(
    tool: &str,
    min_version: u32,
    version_source: Option<&str>,
) -> (bool, Option<String>) {
    match tool {
        "java" => check_java_version_with_info(min_version),
        _ => (true, get_version_string(tool, version_source)),
    }
}

fn get_version_string(tool: &str, version_source: Option<&str>) -> Option<String> {
    match tool {
        "docker" => get_docker_version(version_source),
        "syftbox" => get_syftbox_version(version_source),
        "nextflow" => get_nextflow_version(version_source),
        "uv" => get_uv_version(version_source.unwrap_or(tool)),
        _ => None,
    }
}

fn command_for_version(path: Option<&str>, fallback: &str) -> Option<Command> {
    Some(Command::new(path.unwrap_or(fallback)))
}

fn get_docker_version(version_source: Option<&str>) -> Option<String> {
    let mut command = command_for_version(version_source, "docker")?;
    let output = command.arg("--version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    version_str
        .split_whitespace()
        .nth(2)
        .map(|v| v.trim_end_matches(',').to_string())
}

fn get_syftbox_version(version_source: Option<&str>) -> Option<String> {
    let mut command = command_for_version(version_source, "syftbox")?;

    // Keep syftbox probes sandboxed to the BioVault home so we don't litter the user's real HOME with ~/.syftbox logs/config.
    if let Ok(home) = crate::config::get_biovault_home() {
        let syftbox_dir = home.join(".syftbox");
        let _ = fs::create_dir_all(syftbox_dir.join("logs"));
        command.env("HOME", &home);
        command.env("SYFTBOX_DATA_DIR", &home);
        command.env("SYFTBOX_CONFIG_PATH", syftbox_dir.join("config.json"));
    }

    let output = command.arg("--version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    version_str.split_whitespace().nth(2).map(|v| v.to_string())
}

fn get_nextflow_version(version_source: Option<&str>) -> Option<String> {
    let mut command = command_for_version(version_source, "nextflow")?;
    configure_child_process(&mut command);
    let output = command.arg("-version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    for line in version_str.lines() {
        if line.trim().starts_with("version ") {
            return line.split_whitespace().nth(1).map(|v| v.to_string());
        }
    }
    None
}

fn check_java_version_with_info(min_version: u32) -> (bool, Option<String>) {
    let mut command = Command::new("java");
    command.arg("-version");
    configure_child_process(&mut command);
    let output = command.output();

    match output {
        Ok(output) => {
            let version_str = String::from_utf8_lossy(&output.stderr);
            if let Some(version) = parse_java_version(&version_str) {
                (version >= min_version, Some(version.to_string()))
            } else {
                (false, None)
            }
        }
        Err(_) => (false, None),
    }
}

fn is_google_colab() -> bool {
    // Check for COLAB_RELEASE_TAG which is specific to Colab
    if env::var("COLAB_RELEASE_TAG").is_ok() {
        return true;
    }

    // Fallback: check for any COLAB_ prefixed environment variable
    for (key, _) in env::vars() {
        if key.starts_with("COLAB_") {
            return true;
        }
    }

    false
}

#[cfg(target_os = "macos")]
fn check_java_in_brew_not_in_path() -> Option<String> {
    // Check known Homebrew opt prefixes first (covers keg-only installs)
    for base in OPENJDK_OPT_BASES.iter() {
        let base_path = Path::new(base);
        if !base_path.exists() {
            continue;
        }

        if let Ok(entries) = std::fs::read_dir(base_path) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                    if name.starts_with("openjdk") {
                        if let Some(bin_dir) = java_bin_from_prefix(path.to_string_lossy().as_ref())
                        {
                            let java_path = format!("{}/java", bin_dir);
                            if is_valid_java_binary(&java_path) {
                                return Some(bin_dir);
                            }
                        }
                    }
                }
            }
        }
    }

    if let Some(brew_cmd) = resolve_brew_command() {
        let formulas = [
            "openjdk",
            "openjdk@25",
            "openjdk@24",
            "openjdk@23",
            "openjdk@21",
            "openjdk@20",
            "openjdk@19",
            "openjdk@18",
            "openjdk@17",
            "openjdk@11",
            "openjdk@8",
        ];

        for formula in formulas.iter() {
            let output = Command::new(&brew_cmd)
                .args(["--prefix", formula])
                .output()
                .ok();

            let output = match output {
                Some(out) if out.status.success() => out,
                _ => continue,
            };

            let prefix = String::from_utf8_lossy(&output.stdout).trim().to_string();
            if prefix.is_empty() {
                continue;
            }

            if let Some(bin_dir) = java_bin_from_prefix(&prefix) {
                let java_path = format!("{}/java", bin_dir);
                if is_valid_java_binary(&java_path) {
                    return Some(bin_dir);
                }
            }
        }
    }

    None
}

fn check_uv_in_windows_not_in_path() -> Option<String> {
    // Check common UV installation locations on Windows
    let possible_locations = [
        // WinGet Links directory (primary location for winget-installed apps)
        env::var("LOCALAPPDATA")
            .ok()
            .map(|p| format!("{}\\Microsoft\\WinGet\\Links\\uv.exe", p)),
        // WinGet typically installs here
        env::var("LOCALAPPDATA")
            .ok()
            .map(|p| format!("{}\\Programs\\uv\\uv.exe", p)),
        // Alternative WinGet location
        env::var("LOCALAPPDATA")
            .ok()
            .map(|p| format!("{}\\uv\\bin\\uv.exe", p)),
        // WinGet may also install to ProgramFiles
        env::var("PROGRAMFILES")
            .ok()
            .map(|p| format!("{}\\uv\\uv.exe", p)),
        // Cargo install location
        env::var("USERPROFILE")
            .ok()
            .map(|p| format!("{}\\.cargo\\bin\\uv.exe", p)),
        // Direct installer location
        env::var("USERPROFILE")
            .ok()
            .map(|p| format!("{}\\.local\\bin\\uv.exe", p)),
        // Another common WinGet location (user-specific)
        env::var("LOCALAPPDATA")
            .ok()
            .map(|p| format!("{}\\Microsoft\\WinGet\\Packages\\astral-sh.uv_Microsoft.Winget.Source_8wekyb3d8bbwe\\uv.exe", p)),
    ];

    for location in possible_locations.iter().flatten() {
        if std::path::Path::new(location).exists() {
            return Some(location.clone());
        }
    }

    None
}

#[cfg(target_os = "windows")]
fn check_java_in_windows_not_in_path() -> Option<String> {
    // Check Microsoft OpenJDK installation locations on Windows
    let program_files =
        env::var("PROGRAMFILES").unwrap_or_else(|_| "C:\\Program Files".to_string());
    let microsoft_dir = std::path::Path::new(&program_files).join("Microsoft");

    if microsoft_dir.exists() {
        // Look for jdk-* directories (e.g., jdk-17.0.17.10-hotspot)
        if let Ok(entries) = std::fs::read_dir(&microsoft_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                    if name.starts_with("jdk-") {
                        let java_exe = path.join("bin").join("java.exe");
                        if java_exe.exists() {
                            return Some(java_exe.display().to_string());
                        }
                    }
                }
            }
        }
    }

    // Also check Oracle JDK locations
    let java_dir = std::path::Path::new(&program_files).join("Java");
    if java_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&java_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                    if name.starts_with("jdk") || name.starts_with("jre") {
                        let java_exe = path.join("bin").join("java.exe");
                        if java_exe.exists() {
                            return Some(java_exe.display().to_string());
                        }
                    }
                }
            }
        }
    }

    // Check AdoptOpenJDK/Eclipse Temurin locations
    let eclipse_dir = std::path::Path::new(&program_files).join("Eclipse Adoptium");
    if eclipse_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&eclipse_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let java_exe = path.join("bin").join("java.exe");
                if java_exe.exists() {
                    return Some(java_exe.display().to_string());
                }
            }
        }
    }

    None
}

fn get_uv_version(uv_path: &str) -> Option<String> {
    let mut command = Command::new(uv_path);
    command.arg("--version");
    configure_child_process(&mut command);
    let output = command.output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    // UV version output format: "uv 0.4.29"
    version_str.split_whitespace().nth(1).map(|v| v.to_string())
}

/// Check if Homebrew is installed on macOS (for library use)
pub fn check_brew_installed() -> Result<bool> {
    // Only relevant on macOS
    #[cfg(not(target_os = "macos"))]
    {
        Ok(true) // Not needed on other platforms
    }

    #[cfg(target_os = "macos")]
    {
        Ok(resolve_brew_command().is_some())
    }
}

/// Install Homebrew on macOS (for library use)
pub fn install_brew() -> Result<String> {
    #[cfg(not(target_os = "macos"))]
    {
        Err(anyhow!("Homebrew installation is only supported on macOS").into())
    }

    #[cfg(target_os = "macos")]
    {
        // If Homebrew is already available, return early with its path.
        if let Ok(existing) = which::which("brew") {
            return Ok(existing.display().to_string());
        }

        log_homebrew_install("Homebrew installer: brew not detected, preparing installation");

        fn command_output(program: &str, args: &[&str]) -> Result<String> {
            let output = Command::new(program)
                .args(args)
                .output()
                .map_err(|e| CliError::Anyhow(anyhow!("Failed to execute '{}': {}", program, e)))?;

            if !output.status.success() {
                return Err(CliError::Anyhow(anyhow!(
                    "'{}' exited with status {:?}",
                    program,
                    output.status.code()
                )));
            }

            let value = String::from_utf8(output.stdout).map_err(|e| {
                CliError::Anyhow(anyhow!(
                    "Output from '{}' was not valid UTF-8: {}",
                    program,
                    e
                ))
            })?;

            Ok(value.trim().to_string())
        }

        fn escape_for_shell_double_quotes(input: &str) -> String {
            input.replace("\\", "\\\\").replace('"', "\\\"")
        }

        fn escape_for_applescript(input: &str) -> String {
            input.replace('"', "\\\"")
        }

        let prompt_message =
            "BioVault needs Homebrew to manage required dependencies. macOS will ask for your administrator password.";

        let target_user = match std::env::var("SUDO_USER") {
            Ok(value) if !value.is_empty() => value,
            _ => match std::env::var("USER") {
                Ok(value) if !value.is_empty() => value,
                _ => {
                    log_homebrew_install(
                        "Homebrew installer: SUDO_USER/USER missing, falling back to 'id -un'",
                    );
                    command_output("id", &["-un"])?
                }
            },
        };

        let user_uid = command_output("id", &["-u", &target_user])?;
        let user_group = command_output("id", &["-gn", &target_user])?;
        log_homebrew_install(&format!(
            "Homebrew installer: target user '{}', uid {}, group {}",
            target_user, user_uid, user_group
        ));

        let user_groups_output = command_output("id", &["-Gn", &target_user])?;
        log_homebrew_install(&format!(
            "Homebrew installer: '{}' groups => {}",
            target_user, user_groups_output
        ));
        let is_admin = user_groups_output
            .split_whitespace()
            .any(|group| group == "admin");
        if !is_admin {
            log_homebrew_install(&format!(
                "Homebrew installer: user '{}' is not a member of the 'admin' group — cannot proceed",
                target_user
            ));
            return Err(CliError::Anyhow(anyhow!(
                "User '{target_user}' does not have administrator privileges required to install Homebrew automatically."
            )));
        }

        let clt_ready = Command::new("xcode-select")
            .arg("-p")
            .status()
            .map(|status| status.success())
            .unwrap_or(false);

        if !clt_ready {
            log_homebrew_install(
                "Homebrew installer: Command Line Tools not detected. Triggering installer",
            );

            let status = Command::new("xcode-select")
                .arg("--install")
                .status()
                .map_err(|e| {
                    CliError::Anyhow(anyhow!(
                        "Failed to launch Command Line Tools installer: {}",
                        e
                    ))
                })?;

            log_homebrew_install(&format!(
                "Homebrew installer: xcode-select --install exited with status {:?}",
                status.code()
            ));

            return Err(CliError::Anyhow(anyhow!(
                "macOS Command Line Tools must be installed. The installer has been launched—complete it, then click 'Check Again'."
            )));
        }

        let arch = command_output("uname", &["-m"])?;
        log_homebrew_install(&format!(
            "Homebrew installer: detected architecture {}",
            arch
        ));

        let brew_prefix = if arch.trim() == "arm64" {
            "/opt/homebrew"
        } else {
            "/usr/local"
        };
        log_homebrew_install(&format!("Homebrew installer: using prefix {}", brew_prefix));

        let home_dir = dirs::home_dir().ok_or_else(|| {
            CliError::Anyhow(anyhow!(
                "Could not determine home directory for '{}'. Please ensure BIOVAULT_HOME is set.",
                target_user
            ))
        })?;
        log_homebrew_install(&format!(
            "Homebrew installer: resolved home directory {}",
            home_dir.to_string_lossy()
        ));
        let home_dir = home_dir.to_string_lossy();

        let escaped_user = escape_for_shell_double_quotes(&target_user);
        let escaped_group = escape_for_shell_double_quotes(&user_group);
        let escaped_uid = escape_for_shell_double_quotes(&user_uid);
        let escaped_home = escape_for_shell_double_quotes(&home_dir);
        let escaped_prefix = escape_for_shell_double_quotes(brew_prefix);

        let user_script_contents = r#"#!/bin/bash
set -euo pipefail

/usr/bin/env git --version >/dev/null 2>&1 || {
  echo "git is required to install Homebrew. Please install Xcode Command Line Tools." >&2
  exit 2
}

if [ -n "${HOMEBREW_PREFIX:-}" ]; then
  BREW_PREFIX="${HOMEBREW_PREFIX}"
else
  ARCH="$(uname -m)"
  if [ "$ARCH" = "arm64" ]; then
    BREW_PREFIX="/opt/homebrew"
  else
    BREW_PREFIX="/usr/local"
  fi
fi

BREW_REPOSITORY="${BREW_PREFIX}/Homebrew"
BREW_BIN_DIR="${BREW_PREFIX}/bin"
BREW_BIN="${BREW_BIN_DIR}/brew"

if [ ! -d "${BREW_REPOSITORY}/.git" ]; then
  rm -rf "${BREW_REPOSITORY}"
  /usr/bin/env git clone https://github.com/Homebrew/brew "${BREW_REPOSITORY}"
else
  /usr/bin/env git -C "${BREW_REPOSITORY}" fetch origin --force
  /usr/bin/env git -C "${BREW_REPOSITORY}" reset --hard origin/HEAD
fi

/bin/mkdir -p "${BREW_BIN_DIR}"
/bin/ln -sf "${BREW_REPOSITORY}/bin/brew" "${BREW_BIN}"

PROFILE_LINE="eval \"$(${BREW_BIN} shellenv)\""

append_if_missing() {
  local file="$1"
  local line="$2"
  if [ -f "$file" ]; then
    if /usr/bin/grep -F "$line" "$file" >/dev/null 2>&1; then
      return
    fi
    printf '\n%s\n' "$line" >>"$file"
  else
    printf '%s\n' "$line" >"$file"
  fi
}

/bin/mkdir -p "$HOME/.config"
append_if_missing "$HOME/.zprofile" "$PROFILE_LINE"
append_if_missing "$HOME/.bash_profile" "$PROFILE_LINE"

"${BREW_BIN}" update --force --quiet
"${BREW_BIN}" analytics off >/dev/null 2>&1 || true
"#;

        let mut user_script = NamedTempFile::new().map_err(|e| {
            CliError::Anyhow(anyhow!(
                "Failed to create Homebrew user installer script: {}",
                e
            ))
        })?;
        user_script
            .write_all(user_script_contents.as_bytes())
            .map_err(|e| {
                CliError::Anyhow(anyhow!(
                    "Failed to write Homebrew user installer script: {}",
                    e
                ))
            })?;
        user_script.flush().map_err(|e| {
            CliError::Anyhow(anyhow!(
                "Failed to flush Homebrew user installer script: {}",
                e
            ))
        })?;
        std::fs::set_permissions(user_script.path(), PermissionsExt::from_mode(0o755)).map_err(
            |e| {
                CliError::Anyhow(anyhow!(
                    "Failed to set permissions on Homebrew user installer script: {}",
                    e
                ))
            },
        )?;

        let user_script_path = user_script.path().to_string_lossy().to_string();
        log_homebrew_install(&format!(
            "Homebrew installer: created user script at {}",
            user_script_path
        ));
        let escaped_user_script_path = escape_for_shell_double_quotes(&user_script_path);

        log_homebrew_install("Homebrew installer: preparing privileged bootstrap steps");

        log_homebrew_install("Homebrew installer: bootstrap will invoke user script via sudo -u");

        let script_contents = format!(
            r#"#!/bin/bash
set -euo pipefail

USER_NAME="{user}"
USER_GROUP="{group}"
USER_UID="{uid}"
USER_HOME="{home}"
USER_INSTALL_SCRIPT="{user_script}"

trap 'rm -f "${{USER_INSTALL_SCRIPT}}"' EXIT

BREW_PREFIX="{prefix}"

paths=(
  "${{BREW_PREFIX}}"
  "${{BREW_PREFIX}}/bin"
  "${{BREW_PREFIX}}/etc"
  "${{BREW_PREFIX}}/include"
  "${{BREW_PREFIX}}/lib"
  "${{BREW_PREFIX}}/sbin"
  "${{BREW_PREFIX}}/share"
  "${{BREW_PREFIX}}/share/man"
  "${{BREW_PREFIX}}/share/man/man1"
  "${{BREW_PREFIX}}/share/zsh"
  "${{BREW_PREFIX}}/share/zsh/site-functions"
  "${{BREW_PREFIX}}/var"
  "${{BREW_PREFIX}}/Homebrew"
  "${{BREW_PREFIX}}/Cellar"
  "${{BREW_PREFIX}}/Caskroom"
)

for path in "${{paths[@]}}"; do
  if [ ! -d "$path" ]; then
    mkdir -p "$path"
  fi
  chown "${{USER_NAME}}:${{USER_GROUP}}" "$path" 2>/dev/null || true
done

/bin/chmod 0755 "${{USER_INSTALL_SCRIPT}}"

sudo -u "${{USER_NAME}}" -H /usr/bin/env HOME="${{USER_HOME}}" USER="${{USER_NAME}}" LOGNAME="${{USER_NAME}}" HOMEBREW_PREFIX="${{BREW_PREFIX}}" NONINTERACTIVE=1 /bin/bash "${{USER_INSTALL_SCRIPT}}"
"#,
            user = escaped_user,
            group = escaped_group,
            uid = escaped_uid,
            home = escaped_home,
            user_script = escaped_user_script_path,
            prefix = escaped_prefix
        );

        let mut temp_script = NamedTempFile::new().map_err(|e| {
            CliError::Anyhow(anyhow!(
                "Failed to create temporary Homebrew installer script: {}",
                e
            ))
        })?;
        temp_script
            .write_all(script_contents.as_bytes())
            .map_err(|e| {
                CliError::Anyhow(anyhow!("Failed to write Homebrew installer script: {}", e))
            })?;
        temp_script.flush().map_err(|e| {
            CliError::Anyhow(anyhow!("Failed to flush Homebrew installer script: {}", e))
        })?;

        let script_path = temp_script.path().to_string_lossy().to_string();
        log_homebrew_install(&format!(
            "Homebrew installer: created privileged bootstrap script at {}",
            script_path
        ));

        let applescript = format!(
            r#"set shellCommand to "/bin/bash " & quoted form of "{}"
do shell script shellCommand with prompt "{}" with administrator privileges"#,
            escape_for_applescript(&script_path),
            escape_for_applescript(prompt_message)
        );

        log_homebrew_install("Homebrew installer: invoking osascript to run bootstrap script");

        let output = Command::new("osascript")
            .arg("-e")
            .arg(&applescript)
            .output()
            .map_err(|e| {
                CliError::Anyhow(anyhow!(
                    "Failed to request administrator privileges for Homebrew: {}",
                    e
                ))
            })?;

        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();

        if !stdout.trim().is_empty() {
            eprintln!("Homebrew installer stdout: {}", stdout.trim());
            log_homebrew_install(&format!("Homebrew installer stdout: {}", stdout.trim()));
        }
        if !stderr.trim().is_empty() {
            eprintln!("Homebrew installer stderr: {}", stderr.trim());
            log_homebrew_install(&format!("Homebrew installer stderr: {}", stderr.trim()));
        }

        if output.status.success() {
            let brew_path = which::which("brew")
                .ok()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "/opt/homebrew/bin/brew".to_string());

            log_homebrew_install(&format!(
                "Homebrew installer: brew detected at {}",
                brew_path
            ));

            Ok(brew_path)
        } else {
            let combined = format!("{}\n{}", stdout, stderr).to_lowercase();
            let was_cancelled = output.status.code() == Some(1)
                || combined.contains("user canceled")
                || combined.contains("user cancelled")
                || combined.contains("canceled")
                || combined.contains("cancelled");

            let reason = if was_cancelled {
                "Homebrew installation was cancelled."
            } else {
                "Homebrew installation failed."
            };

            let mut detail = String::new();
            if !stderr.trim().is_empty() {
                detail = format!(" Details: {}", stderr.trim());
            } else if !stdout.trim().is_empty() {
                detail = format!(" Details: {}", stdout.trim());
            }

            log_homebrew_install(&format!(
                "Homebrew installer: osascript exit code {:?}, cancelled: {}",
                output.status.code(),
                was_cancelled
            ));

            Err(anyhow!(
                "{} BioVault needs Homebrew to manage dependencies on macOS. Please install Homebrew manually from https://brew.sh and try again.{}",
                reason,
                detail
            )
            .into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn java_version_parsing_various_formats() {
        let samples = [
            ("openjdk version \"17.0.2\" 2022-01-18", Some(17)),
            ("java version \"1.8.0_321\"", Some(8)),
            ("openjdk version \"11.0.14\" 2022-01-18", Some(11)),
            ("java version \"21\" 2023-09-19", Some(21)),
            ("no version here", None),
        ];
        for (out, expected) in samples {
            assert_eq!(parse_java_version(out), expected);
        }
    }

    #[test]
    fn start_command_and_running_checks() {
        // Test Docker command - it should return OS-specific instructions
        let docker_cmd = get_start_command("docker");

        // Check that it returns something meaningful for Docker
        // The exact message depends on the OS running the test
        if std::env::consts::OS == "macos" {
            assert_eq!(
                docker_cmd,
                "Open Docker Desktop from your Applications folder (/Applications/Docker.app)"
            );
        } else if std::env::consts::OS == "windows" {
            assert_eq!(docker_cmd, "Open Docker Desktop from your Start menu");
        } else {
            assert_eq!(
                docker_cmd,
                "Start Docker daemon with 'sudo systemctl start docker' or open Docker Desktop"
            );
        }

        // Unknown service -> generic
        assert_eq!(get_start_command("xyz"), "Start xyz");

        // Unknown service not considered running by default
        assert!(!check_if_running("xyz"));
    }

    #[test]
    fn check_version_non_java_defaults_true() {
        assert!(check_version("not-java", 9999));
    }

    #[test]
    fn test_parse_java_version_edge_cases() {
        // Test various edge cases
        assert_eq!(parse_java_version(""), None);
        assert_eq!(parse_java_version("version"), None);
        assert_eq!(parse_java_version("version \"\""), None);
        assert_eq!(parse_java_version("version \"not a number\""), None);

        // Test versions without quotes
        assert_eq!(parse_java_version("java 17"), None);

        // Test malformed versions
        assert_eq!(parse_java_version("version \"1.\""), None);
        assert_eq!(parse_java_version("version \"1.x.0\""), None);
    }

    #[test]
    fn test_parse_java_version_modern_formats() {
        // Test various modern Java version formats
        assert_eq!(
            parse_java_version("openjdk version \"17\" 2021-09-14"),
            Some(17)
        );
        assert_eq!(
            parse_java_version("openjdk version \"21.0.1\" 2023-10-17 LTS"),
            Some(21)
        );
        assert_eq!(
            parse_java_version("java version \"19.0.2\" 2023-01-17"),
            Some(19)
        );
    }

    #[test]
    fn test_parse_java_version_legacy_formats() {
        // Test legacy Java version formats (1.x style)
        assert_eq!(parse_java_version("java version \"1.7.0_80\""), Some(7));
        assert_eq!(parse_java_version("java version \"1.6.0_45\""), Some(6));
        assert_eq!(parse_java_version("java version \"1.8.0_351\""), Some(8));
    }

    #[test]
    #[serial_test::serial]
    fn test_google_colab_detection() {
        // Save current env state
        let was_set = env::var("COLAB_RELEASE_TAG").is_ok();
        let old_value = env::var("COLAB_RELEASE_TAG").ok();

        // Test detection when COLAB_RELEASE_TAG is set
        env::set_var("COLAB_RELEASE_TAG", "release-123");
        assert!(is_google_colab());

        // Clean up
        if was_set {
            if let Some(val) = old_value {
                env::set_var("COLAB_RELEASE_TAG", val);
            }
        } else {
            env::remove_var("COLAB_RELEASE_TAG");
        }
    }

    #[test]
    #[serial_test::serial]
    fn test_google_colab_detection_prefix() {
        // Test detection with other COLAB_ prefixed vars
        let was_set = env::var("COLAB_TEST_VAR").is_ok();

        env::set_var("COLAB_TEST_VAR", "test");
        assert!(is_google_colab());

        // Clean up
        if !was_set {
            env::remove_var("COLAB_TEST_VAR");
        }
    }

    #[test]
    #[serial_test::serial]
    fn execute_reports_missing_in_clean_env() {
        // Force an empty PATH so no dependencies are found
        let old_path = std::env::var("PATH").unwrap_or_default();
        std::env::set_var("PATH", "");

        let rt = tokio::runtime::Runtime::new().unwrap();
        let res = rt.block_on(execute(false));
        assert!(res.is_err());

        // Restore PATH
        std::env::set_var("PATH", old_path);
    }

    #[test]
    #[serial_test::serial]
    #[cfg_attr(
        not(feature = "slow-tests"),
        ignore = "env-dependent; covered in slow/integration"
    )]
    #[cfg_attr(windows, ignore = "Windows shell semantics; covered in integration CI")]
    fn execute_with_all_tools_available_returns_ok() {
        let dir = TempDir::new().unwrap();
        // Create fake tools in a temp dir and prepend to PATH
        let make_exec = |name: &str, body: &str| {
            let p = dir.path().join(name);
            fs::write(&p, format!("#!/bin/sh\n{}\n", body)).unwrap();
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut perm = fs::metadata(&p).unwrap().permissions();
                perm.set_mode(0o755);
                fs::set_permissions(&p, perm).unwrap();
            }
            p
        };

        // java prints version to stderr
        make_exec("java", "echo 'openjdk version \"21\" 2024-01-01' 1>&2");
        // docker returns success on 'info'
        make_exec("docker", "[ \"$1\" = \"info\" ] && exit 0; exit 0");
        // nextflow and syftbox stubs
        make_exec("nextflow", "exit 0");
        make_exec("syftbox", "exit 0");
        // uv stub
        make_exec("uv", "exit 0");

        // Prepend temp dir to PATH
        let old_path = std::env::var("PATH").unwrap_or_default();
        let new_path = format!("{}:{}", dir.path().display(), old_path);
        std::env::set_var("PATH", &new_path);

        // Run
        let rt = tokio::runtime::Runtime::new().unwrap();
        let res = rt.block_on(execute(false));
        assert!(res.is_ok());

        // Restore PATH
        std::env::set_var("PATH", old_path);
    }

    #[test]
    #[serial_test::serial]
    fn execute_with_docker_not_running_reports_warning() {
        let dir = TempDir::new().unwrap();
        let make_exec = |name: &str, body: &str| {
            let p = dir.path().join(name);
            fs::write(&p, format!("#!/bin/sh\n{}\n", body)).unwrap();
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut perm = fs::metadata(&p).unwrap().permissions();
                perm.set_mode(0o755);
                fs::set_permissions(&p, perm).unwrap();
            }
            p
        };

        // Provide java and nextflow and syftbox and uv present
        make_exec("java", "echo 'openjdk version \"21\"' 1>&2");
        make_exec("nextflow", "exit 0");
        make_exec("syftbox", "exit 0");
        make_exec("uv", "exit 0");
        // docker exists but 'info' fails
        make_exec("docker", "[ \"$1\" = \"info\" ] && exit 1; exit 0");

        // Prepend
        let old_path = std::env::var("PATH").unwrap_or_default();
        let new_path = format!("{}:{}", dir.path().display(), old_path);
        std::env::set_var("PATH", &new_path);

        // Run: should return Err because a service is not running
        let rt = tokio::runtime::Runtime::new().unwrap();
        let res = rt.block_on(execute(false));
        assert!(res.is_err());

        std::env::set_var("PATH", old_path);
    }
}
