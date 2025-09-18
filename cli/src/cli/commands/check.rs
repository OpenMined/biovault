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

pub async fn execute() -> Result<()> {
    // Load the deps.yaml file embedded in the binary
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    println!("BioVault Dependency Check");
    println!("=========================\n");

    // Detect if we're in Google Colab
    let is_colab = is_google_colab();
    if is_colab {
        println!("ℹ️  Google Colab environment detected\n");
    }

    let mut all_found = true;
    let mut all_running = true;

    for dep in &config.dependencies {
        // Check if this dependency should be skipped in Colab
        if is_colab {
            if let Some(environments) = &dep.environments {
                if let Some(colab_config) = environments.get("google_colab") {
                    if colab_config.skip {
                        println!("Checking {}... ⏭️  SKIPPED", dep.name);
                        println!(
                            "  Reason: {}",
                            colab_config
                                .skip_reason
                                .as_ref()
                                .unwrap_or(&"Not available in Colab".to_string())
                        );
                        println!();
                        continue;
                    }
                }
            }
        }

        print!("Checking {}... ", dep.name);

        // Check if the binary exists in PATH
        let exists = which::which(&dep.name).is_ok();

        if !exists {
            all_found = false;
            println!("❌ NOT FOUND");
            println!("  Description: {}", dep.description);
            println!("  Installation instructions:");
            for line in dep.install_instructions.lines() {
                if !line.trim().is_empty() {
                    println!("    {}", line);
                }
            }
            println!();
        } else {
            // Check version requirement if specified
            if let Some(min_version) = dep.min_version {
                let version_ok = check_version(&dep.name, min_version);
                if !version_ok {
                    all_found = false;
                    println!("❌ Version too old (requires {} or higher)", min_version);
                    println!("  Description: {}", dep.description);
                    println!("  Installation instructions:");
                    for line in dep.install_instructions.lines() {
                        if !line.trim().is_empty() {
                            println!("    {}", line);
                        }
                    }
                    println!();
                } else {
                    print!("✓ Found");
                }
            } else {
                print!("✓ Found");
            }

            // Check if it needs to be running and if it is
            if dep.check_running {
                let is_running = check_if_running(&dep.name);
                if is_running {
                    println!(" (running)");
                } else {
                    all_running = false;
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
        }
    }

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
        println!("⚠️  Some services are not running. Please start them using the commands above.");
        Err(anyhow!("Services not running").into())
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
        "docker" => "Open Docker Desktop or run 'sudo dockerd' (Linux)".to_string(),
        _ => format!("Start {}", service),
    }
}

fn check_version(tool: &str, min_version: u32) -> bool {
    match tool {
        "java" => check_java_version(min_version),
        _ => true, // For tools without version checking, assume OK
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(
            get_start_command("docker"),
            "Open Docker Desktop or run 'sudo dockerd' (Linux)"
        );
        // Unknown service -> generic
        assert_eq!(get_start_command("xyz"), "Start xyz");

        // Unknown service not considered running by default
        assert!(!check_if_running("xyz"));
    }

    #[test]
    fn check_version_non_java_defaults_true() {
        assert!(check_version("not-java", 9999));
    }
}
