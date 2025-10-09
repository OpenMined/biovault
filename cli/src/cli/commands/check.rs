use crate::config::Config;
use crate::Result;
use anyhow::anyhow;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::process::Command;

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
    pub environments: Option<HashMap<String, EnvironmentConfig>>,
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

#[derive(Debug, Serialize, Deserialize)]
pub struct DependencyCheckResult {
    pub dependencies: Vec<DependencyResult>,
    pub all_satisfied: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DependencyResult {
    pub name: String,
    pub found: bool,
    pub path: Option<String>,
    pub version: Option<String>,
    pub running: Option<bool>,
    pub skipped: bool,
    pub skip_reason: Option<String>,
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

    for dep in &config.dependencies {
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
        let binary_path = if let Some(custom) = custom_path {
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

        // Special handling for Java on macOS
        let mut java_brew_path: Option<String> = None;
        if binary_path.is_none() && dep.name == "java" && std::env::consts::OS == "macos" {
            // Check if Java is installed via Homebrew but not in PATH
            java_brew_path = check_java_in_brew_not_in_path();
        }

        if binary_path.is_none() && java_brew_path.is_none() {
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
            });
        } else if let Some(brew_path) = java_brew_path {
            // Java is installed via brew but not in PATH
            if !json {
                println!("⚠️  Found (not in PATH)");
                println!("  Java is installed via Homebrew at: {}", brew_path);
                println!("  But it's not available in your PATH.");
                println!("  To fix this, add the following to your shell config:");
                let shell = env::var("SHELL").unwrap_or_else(|_| "/bin/zsh".to_string());
                let shell_config = if shell.contains("zsh") {
                    "~/.zshrc"
                } else if shell.contains("bash") {
                    "~/.bash_profile"
                } else {
                    "your shell config file"
                };
                println!(
                    "    echo 'export PATH=\"{}:$PATH\"' >> {}",
                    brew_path, shell_config
                );
                println!("    source {}", shell_config);
                if !is_ci {
                    println!("  Or run 'bv setup' to configure this automatically.");
                }
                println!();
            }
            results.push(DependencyResult {
                name: dep.name.clone(),
                found: true,
                path: Some(format!("{}/java", brew_path)),
                version: None,
                running: None,
                skipped: false,
                skip_reason: None,
            });
        } else if let Some(path) = binary_path {
            // Binary found - check version
            let (version_ok, version_str) = if let Some(min_version) = dep.min_version {
                check_version_with_info(&dep.name, min_version)
            } else {
                (true, get_version_string(&dep.name))
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
                    path: Some(path),
                    version: version_str,
                    running: None,
                    skipped: false,
                    skip_reason: None,
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
                    path: Some(path),
                    version: version_str,
                    running: is_running,
                    skipped: false,
                    skip_reason: None,
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
            Command::new("docker")
                .arg("info")
                .output()
                .map(|output| output.status.success())
                .unwrap_or(false)
        }
        _ => false,
    }
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
    let output = Command::new("java").arg("-version").output();

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

fn check_version_with_info(tool: &str, min_version: u32) -> (bool, Option<String>) {
    match tool {
        "java" => check_java_version_with_info(min_version),
        _ => (true, get_version_string(tool)),
    }
}

fn get_version_string(tool: &str) -> Option<String> {
    match tool {
        "docker" => get_docker_version(),
        "syftbox" => get_syftbox_version(),
        "nextflow" => get_nextflow_version(),
        _ => None,
    }
}

fn get_docker_version() -> Option<String> {
    let output = Command::new("docker").arg("--version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    version_str
        .split_whitespace()
        .nth(2)
        .map(|v| v.trim_end_matches(',').to_string())
}

fn get_syftbox_version() -> Option<String> {
    let output = Command::new("syftbox").arg("--version").output().ok()?;
    if !output.status.success() {
        return None;
    }
    let version_str = String::from_utf8_lossy(&output.stdout);
    version_str.split_whitespace().nth(2).map(|v| v.to_string())
}

fn get_nextflow_version() -> Option<String> {
    let output = Command::new("nextflow").arg("-version").output().ok()?;
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
    let output = Command::new("java").arg("-version").output();

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

fn check_java_in_brew_not_in_path() -> Option<String> {
    // Check if brew command exists
    if which::which("brew").is_err() {
        return None;
    }

    // Check if Java/OpenJDK is installed via brew
    let output = Command::new("brew")
        .args(["list", "--formula"])
        .output()
        .ok()?;

    let installed_packages = String::from_utf8_lossy(&output.stdout);

    // Look for any OpenJDK version
    let mut found_java_package = None;
    for line in installed_packages.lines() {
        if line.starts_with("openjdk") {
            found_java_package = Some(line.to_string());
            break;
        }
    }

    found_java_package.as_ref()?;

    // Get the actual path where brew installed Java
    let pkg = found_java_package.unwrap();
    let prefix_output = Command::new("brew")
        .args(["--prefix", &pkg])
        .output()
        .ok()?;

    if !prefix_output.status.success() {
        return None;
    }

    let brew_prefix = String::from_utf8_lossy(&prefix_output.stdout)
        .trim()
        .to_string();
    let java_bin_path = format!("{}/bin", brew_prefix);

    // Check if this path contains java binary
    if std::path::Path::new(&format!("{}/java", java_bin_path)).exists() {
        Some(java_bin_path)
    } else {
        None
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
