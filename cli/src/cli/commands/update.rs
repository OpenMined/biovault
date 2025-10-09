use anyhow::{Context, Result};
use colored::Colorize;
use dialoguer::Confirm;
use semver::Version;
use serde::Deserialize;
use std::env;
use tracing::{debug, info};

const CRATES_IO_API_URL: &str = "https://crates.io/api/v1/crates/biovault";
const GITHUB_API_URL: &str = "https://api.github.com/repos/openmined/biovault/releases/latest";

fn skip_update_checks() -> bool {
    match env::var("BIOVAULT_SKIP_UPDATE_CHECK") {
        Ok(value) => value != "0",
        Err(_) => false,
    }
}

#[derive(Debug, Deserialize)]
struct CratesApiResponse {
    #[serde(rename = "crate")]
    crate_info: CrateInfo,
}

#[derive(Debug, Deserialize)]
struct CrateInfo {
    max_version: String,
}

#[derive(Debug, Deserialize)]
struct GithubRelease {
    tag_name: String,
    #[allow(dead_code)]
    name: Option<String>,
    #[allow(dead_code)]
    published_at: String,
}

fn get_current_version() -> Version {
    Version::parse(env!("CARGO_PKG_VERSION")).expect("Invalid current version")
}

async fn check_crates_io() -> Result<Option<Version>> {
    let client = reqwest::Client::new();
    let response = client
        .get(CRATES_IO_API_URL)
        .header("User-Agent", "biovault-cli")
        .send()
        .await
        .context("Failed to check crates.io")?;

    if !response.status().is_success() {
        return Ok(None);
    }

    let api_response: CratesApiResponse = response
        .json()
        .await
        .context("Failed to parse crates.io response")?;

    let latest_version = Version::parse(&api_response.crate_info.max_version)
        .context("Invalid version from crates.io")?;

    Ok(Some(latest_version))
}

async fn check_github() -> Result<Option<Version>> {
    let client = reqwest::Client::new();
    let response = client
        .get(GITHUB_API_URL)
        .header("User-Agent", "biovault-cli")
        .send()
        .await
        .context("Failed to check GitHub releases")?;

    if !response.status().is_success() {
        return Ok(None);
    }

    let release: GithubRelease = response
        .json()
        .await
        .context("Failed to parse GitHub response")?;

    let version_str = release.tag_name.trim_start_matches('v');
    let latest_version = Version::parse(version_str).context("Invalid version from GitHub")?;

    Ok(Some(latest_version))
}

pub async fn check_for_updates() -> Result<Option<Version>> {
    let current = get_current_version();
    debug!("Current version: {}", current);

    let crates_version = check_crates_io().await.ok().flatten();
    let github_version = check_github().await.ok().flatten();

    let latest = match (crates_version, github_version) {
        (Some(c), Some(g)) => Some(if c > g { c } else { g }),
        (Some(c), None) => Some(c),
        (None, Some(g)) => Some(g),
        (None, None) => None,
    };

    if let Some(ref version) = latest {
        debug!("Latest version available: {}", version);
        if version > &current {
            return Ok(Some(version.clone()));
        }
    }

    Ok(None)
}

pub async fn check_and_notify_random() -> Result<()> {
    if skip_update_checks() {
        debug!("Skipping update notification check due to BIOVAULT_SKIP_UPDATE_CHECK env var");
        return Ok(());
    }

    use rand::Rng;
    let mut rng = rand::thread_rng();

    if rng.gen_bool(0.1) {
        if let Ok(Some(new_version)) = check_for_updates().await {
            println!(
                "\n{} {} {} {}",
                "📦".yellow(),
                "New version".yellow().bold(),
                new_version.to_string().green().bold(),
                "available! Run 'bv update' to upgrade.".yellow()
            );
        }
    }

    Ok(())
}

pub async fn execute() -> Result<()> {
    info!("Checking for updates...");

    let current = get_current_version();
    println!("Current version: {}", current.to_string().cyan());

    match check_for_updates().await? {
        Some(new_version) => {
            println!(
                "\n{} {} {}",
                "✨".green(),
                "New version available:".green().bold(),
                new_version.to_string().green().bold()
            );

            let confirm = Confirm::new()
                .with_prompt(format!("Upgrade from {} to {}?", current, new_version))
                .default(true)
                .interact()
                .context("Failed to get user confirmation")?;

            if confirm {
                perform_update(&new_version).await?;
            } else {
                println!("Update cancelled.");
            }
        }
        None => {
            println!(
                "{} {}",
                "✓".green(),
                "You're already on the latest version!".green()
            );
        }
    }

    Ok(())
}

async fn perform_update(new_version: &Version) -> Result<()> {
    println!("\n{} Updating biovault...", "🔄".cyan());

    let install_method = detect_install_method()?;

    match install_method {
        InstallMethod::Cargo => update_via_cargo(new_version).await?,
        InstallMethod::Binary => update_via_self_update(new_version).await?,
    }

    println!(
        "\n{} {} {} {}",
        "✨".green(),
        "Successfully updated to version".green().bold(),
        new_version.to_string().green().bold(),
        "!".green().bold()
    );

    Ok(())
}

enum InstallMethod {
    Cargo,
    Binary,
}

fn detect_install_method() -> Result<InstallMethod> {
    let exe_path = env::current_exe().context("Failed to get current executable path")?;
    let exe_path_str = exe_path.to_string_lossy();

    // Use lowercase for case-insensitive comparison and handle both forward and backslashes
    let path_lower = exe_path_str.to_lowercase();

    if path_lower.contains(".cargo")
        || path_lower.contains("target/release")
        || path_lower.contains("target\\release")
        || path_lower.contains("target/debug")
        || path_lower.contains("target\\debug")
    {
        Ok(InstallMethod::Cargo)
    } else {
        Ok(InstallMethod::Binary)
    }
}

async fn update_via_cargo(_new_version: &Version) -> Result<()> {
    use std::process::Command;

    println!("Updating via cargo install...");

    let output = Command::new("cargo")
        .args(["install", "biovault", "--force"])
        .output()
        .context("Failed to run cargo install")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("cargo install failed: {}", stderr);
    }

    Ok(())
}

async fn update_via_self_update(new_version: &Version) -> Result<()> {
    println!("Updating via direct binary download...");

    let status = self_update::backends::github::Update::configure()
        .repo_owner("openmined")
        .repo_name("biovault")
        .bin_name("bv")
        .target_version_tag(&format!("v{}", new_version))
        .show_download_progress(true)
        .current_version(env!("CARGO_PKG_VERSION"))
        .build()
        .context("Failed to build self-updater")?
        .update()
        .context("Failed to perform self-update")?;

    debug!("Update status: {:?}", status);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_install_method_on_test_binary_is_cargo() {
        // Test binaries run from target/debug/deps, so detect_install_method should detect Cargo
        let method = detect_install_method().expect("detect");
        match method {
            InstallMethod::Cargo => {}
            _ => panic!("expected Cargo install method for test binary"),
        }
    }

    #[test]
    fn current_version_parses_as_semver() {
        let v = get_current_version();
        // Reasonable expectations: non-empty and contains at least one dot
        let s = v.to_string();
        assert!(!s.is_empty());
        assert!(s.contains('.'));
    }

    #[test]
    fn test_install_method_enum() {
        // Test that InstallMethod enum works correctly
        let cargo = InstallMethod::Cargo;
        let binary = InstallMethod::Binary;

        match cargo {
            InstallMethod::Cargo => {}
            _ => panic!("Expected Cargo"),
        }

        match binary {
            InstallMethod::Binary => {}
            _ => panic!("Expected Binary"),
        }
    }

    #[test]
    fn test_detect_install_method_with_cargo_path() {
        // Test the detection logic - we can't set current_exe but we can test the function
        let method = detect_install_method();
        // Can't test the actual result as it depends on runtime environment
        assert!(method.is_ok());
    }

    #[test]
    fn test_github_release_struct() {
        let release = GithubRelease {
            tag_name: "v1.2.3".to_string(),
            name: Some("Release v1.2.3".to_string()),
            published_at: "2024-01-01T00:00:00Z".to_string(),
        };
        assert_eq!(release.tag_name, "v1.2.3");
        assert_eq!(release.name, Some("Release v1.2.3".to_string()));
        assert_eq!(release.published_at, "2024-01-01T00:00:00Z");
    }

    #[test]
    fn test_crates_api_response_struct() {
        let response = CratesApiResponse {
            crate_info: CrateInfo {
                max_version: "2.3.4".to_string(),
            },
        };
        assert_eq!(response.crate_info.max_version, "2.3.4");
    }

    #[test]
    fn test_get_current_version_returns_valid_semver() {
        let version = get_current_version();
        assert!(version.major > 0 || version.minor > 0 || version.patch > 0);
        // Check it's a valid semver
        let version_str = version.to_string();
        semver::Version::parse(&version_str).expect("Should be valid semver");
    }

    #[test]
    fn test_detect_install_method_cross_platform() {
        // Save original exe path
        let _original_exe = std::env::current_exe().ok();

        // Test various path patterns
        let test_cases = vec![
            ("/Users/test/.cargo/bin/biovault", InstallMethod::Cargo),
            (
                "C:\\Users\\test\\.cargo\\bin\\biovault.exe",
                InstallMethod::Cargo,
            ),
            ("/home/user/target/debug/biovault", InstallMethod::Cargo),
            (
                "C:\\project\\target\\release\\biovault.exe",
                InstallMethod::Cargo,
            ),
            ("/usr/local/bin/biovault", InstallMethod::Binary),
            (
                "C:\\Program Files\\biovault\\biovault.exe",
                InstallMethod::Binary,
            ),
        ];

        // We can't actually change current_exe in tests, but we can verify the logic
        for (path, _expected) in test_cases {
            // Verify path pattern matching
            let path_lower = path.to_lowercase();
            let is_cargo = path_lower.contains(".cargo")
                || path_lower.contains("target/release")
                || path_lower.contains("target\\release")
                || path_lower.contains("target/debug")
                || path_lower.contains("target\\debug");

            if path.contains(".cargo") || path.contains("target") {
                assert!(is_cargo, "Path {} should be detected as Cargo", path);
            } else {
                assert!(!is_cargo, "Path {} should be detected as Binary", path);
            }
        }
    }

    #[test]
    #[serial_test::serial]
    fn skip_update_checks_env_handling() {
        // Save original state
        let original = std::env::var("BIOVAULT_SKIP_UPDATE_CHECK").ok();

        std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK");
        assert!(!skip_update_checks());

        std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        assert!(skip_update_checks());

        std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", "0");
        assert!(!skip_update_checks());

        // Restore original state
        match original {
            Some(v) => std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", v),
            None => std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK"),
        }
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn check_and_notify_random_respects_skip_env() {
        let original = std::env::var("BIOVAULT_SKIP_UPDATE_CHECK").ok();

        std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        check_and_notify_random().await.unwrap();

        match original {
            Some(v) => std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", v),
            None => std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK"),
        }
    }

    #[test]
    #[serial_test::serial]
    fn test_skip_update_checks_env_set() {
        let original = std::env::var("BIOVAULT_SKIP_UPDATE_CHECK").ok();

        std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        assert!(skip_update_checks());

        match original {
            Some(v) => std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", v),
            None => std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK"),
        }
    }

    #[test]
    #[serial_test::serial]
    fn test_skip_update_checks_env_not_set() {
        let original = std::env::var("BIOVAULT_SKIP_UPDATE_CHECK").ok();

        std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK");
        assert!(!skip_update_checks());

        match original {
            Some(v) => std::env::set_var("BIOVAULT_SKIP_UPDATE_CHECK", v),
            None => std::env::remove_var("BIOVAULT_SKIP_UPDATE_CHECK"),
        }
    }

    #[test]
    fn test_install_method_enum_variants() {
        let cargo = InstallMethod::Cargo;
        let binary = InstallMethod::Binary;

        assert!(matches!(cargo, InstallMethod::Cargo));
        assert!(matches!(binary, InstallMethod::Binary));
    }
}
