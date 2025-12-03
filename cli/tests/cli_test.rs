use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use std::sync::Mutex;
use tempfile::TempDir;

// Use a mutex to ensure HOME modification tests don't run in parallel
static HOME_MUTEX: Mutex<()> = Mutex::new(());

#[test]
fn test_cli_version() {
    // Read version from Cargo.toml to avoid hardcoding
    let cargo_toml = include_str!("../Cargo.toml");
    let version = cargo_toml
        .lines()
        .find(|line| line.starts_with("version = "))
        .and_then(|line| line.split('"').nth(1))
        .expect("Version not found in Cargo.toml");

    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(version));
}

#[test]
fn test_cli_help() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("BioVault"))
        .stdout(predicate::str::contains("init"))
        .stdout(predicate::str::contains("info"))
        .stdout(predicate::str::contains("check"));
}

#[test]
fn test_init_command() {
    let _guard = HOME_MUTEX
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());

    let temp_dir = TempDir::new().unwrap();
    // BIOVAULT_TEST_HOME uses the path directly without .biovault suffix
    let config_file = temp_dir.path().join("config.yaml");

    // Save original HOME/USERPROFILE and set temporary one
    let original_home = if cfg!(windows) {
        // On Windows, dirs::home_dir() may use multiple env vars
        let original_profile = std::env::var("USERPROFILE").ok();
        std::env::set_var("USERPROFILE", temp_dir.path());
        original_profile
    } else {
        let original_home = std::env::var("HOME").ok();
        std::env::set_var("HOME", temp_dir.path());
        original_home
    };

    let mut cmd = Command::cargo_bin("bv").unwrap();
    // Use our test-specific env var that the init command will respect
    cmd.env("BIOVAULT_TEST_HOME", temp_dir.path());
    // Also set platform-specific home for other commands that might use dirs::home_dir()
    if cfg!(windows) {
        cmd.env("USERPROFILE", temp_dir.path());
        // Clear other Windows home-related env vars to prevent interference
        cmd.env_remove("HOMEDRIVE");
        cmd.env_remove("HOMEPATH");
    } else {
        cmd.env("HOME", temp_dir.path());
    }
    cmd.arg("init")
        .arg("test@example.com")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "BioVault initialized successfully",
        ));

    // Check that config file was created
    assert!(
        config_file.exists(),
        "Config file should exist at: {}",
        config_file.display()
    );

    // Check config file contents
    let contents = fs::read_to_string(&config_file).unwrap();
    assert!(contents.contains("email: test@example.com"));

    // Restore original HOME/USERPROFILE
    if cfg!(windows) {
        if let Some(home) = original_home {
            std::env::set_var("USERPROFILE", home);
        } else {
            std::env::remove_var("USERPROFILE");
        }
    } else if let Some(home) = original_home {
        std::env::set_var("HOME", home);
    } else {
        std::env::remove_var("HOME");
    }
}

#[test]
fn test_init_command_existing_config() {
    let _guard = HOME_MUTEX
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());

    let temp_dir = TempDir::new().unwrap();
    // BIOVAULT_TEST_HOME uses the path directly without .biovault suffix
    let config_file = temp_dir.path().join("config.yaml");

    // Create existing config
    fs::write(&config_file, "email: existing@example.com\n").unwrap();

    // Save original HOME/USERPROFILE and set temporary one
    let original_home = if cfg!(windows) {
        // On Windows, dirs::home_dir() may use multiple env vars
        let original_profile = std::env::var("USERPROFILE").ok();
        std::env::set_var("USERPROFILE", temp_dir.path());
        original_profile
    } else {
        let original_home = std::env::var("HOME").ok();
        std::env::set_var("HOME", temp_dir.path());
        original_home
    };

    let mut cmd = Command::cargo_bin("bv").unwrap();
    // Use our test-specific env var that the init command will respect
    cmd.env("BIOVAULT_TEST_HOME", temp_dir.path());
    // Also set platform-specific home for other commands that might use dirs::home_dir()
    if cfg!(windows) {
        cmd.env("USERPROFILE", temp_dir.path());
        // Clear other Windows home-related env vars to prevent interference
        cmd.env_remove("HOMEDRIVE");
        cmd.env_remove("HOMEPATH");
    } else {
        cmd.env("HOME", temp_dir.path());
    }
    cmd.arg("init")
        .arg("new@example.com")
        .assert()
        .success()
        .stdout(predicate::str::contains("already exists"));

    // Check that original config is unchanged
    let contents = fs::read_to_string(&config_file).unwrap();
    assert!(contents.contains("email: existing@example.com"));
    assert!(!contents.contains("email: new@example.com"));

    // Restore original HOME/USERPROFILE
    if cfg!(windows) {
        if let Some(home) = original_home {
            std::env::set_var("USERPROFILE", home);
        } else {
            std::env::remove_var("USERPROFILE");
        }
    } else if let Some(home) = original_home {
        std::env::set_var("HOME", home);
    } else {
        std::env::remove_var("HOME");
    }
}

#[test]
fn test_info_command() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("info")
        .assert()
        .success()
        .stdout(predicate::str::contains("System Information"))
        .stdout(predicate::str::contains("OS:"))
        .stdout(predicate::str::contains("CPU Arch:"))
        .stdout(predicate::str::contains("CPUs:"))
        .stdout(predicate::str::contains("RAM:"))
        .stdout(predicate::str::contains("DISK FREE:"));
}

#[test]
fn test_project_examples_cli() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("project")
        .arg("examples")
        .assert()
        .success()
        .stdout(predicate::str::contains("Available example templates"));
}

#[test]
fn test_sample_data_list_cli() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("sample-data")
        .arg("list")
        .assert()
        .success()
        .stdout(predicate::str::contains("Available sample data"));
}

#[test]
fn test_run_dry_run_cli() {
    // Skip test if Docker is not available (CI environments may not have it)
    if std::process::Command::new("docker")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("Skipping test_run_dry_run_cli: Docker not available");
        return;
    }

    let tmp = TempDir::new().unwrap();
    // Set explicit BIOVAULT_HOME for templates
    let bv_home = tmp.path().join(".bvhome");
    fs::create_dir_all(bv_home.join("env/dynamic-nextflow")).unwrap();
    // Copy embedded template
    let template_content = include_str!("../src/templates/dynamic/template.nf");
    fs::write(
        bv_home.join("env/dynamic-nextflow/template.nf"),
        template_content,
    )
    .unwrap();
    fs::write(
        bv_home.join("env/dynamic-nextflow/nextflow.config"),
        "process.executor = 'local'\n",
    )
    .unwrap();

    // Prepare minimal dynamic project with optional input
    let proj = tmp.path().join("proj");
    fs::create_dir_all(&proj).unwrap();
    fs::write(
        proj.join("project.yaml"),
        "name: p\nauthor: a\nworkflow: workflow.nf\ntemplate: dynamic-nextflow\ninputs:\n  - name: data\n    type: File?\n",
    )
    .unwrap();
    fs::write(proj.join("workflow.nf"), "workflow USER { }").unwrap();

    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.env("BIOVAULT_HOME", &bv_home)
        .arg("run")
        .arg("--dry-run")
        .arg(proj.to_string_lossy().to_string())
        .assert()
        .success()
        .stdout(predicate::str::contains("Dry run - would execute:"));
}

// Removed test_check_command: behavior is covered by OS/nightly installation tests.

#[test]
fn test_invalid_command() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("invalid-command")
        .assert()
        .failure()
        .stderr(predicate::str::contains("unrecognized subcommand"));
}

#[test]
fn test_init_missing_email() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("init")
        .assert()
        .failure()
        .stderr(predicate::str::contains("required"));
}
