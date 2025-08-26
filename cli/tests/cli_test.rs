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
    let _guard = HOME_MUTEX.lock().unwrap();

    let temp_dir = TempDir::new().unwrap();
    let config_dir = temp_dir.path().join(".biovault");

    // Save original HOME and set temporary HOME
    let original_home = std::env::var("HOME").ok();
    std::env::set_var("HOME", temp_dir.path());

    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.env("HOME", temp_dir.path())
        .arg("init")
        .arg("test@example.com")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "BioVault initialized successfully",
        ));

    // Check that config file was created
    let config_file = config_dir.join("config.yaml");
    assert!(config_file.exists(), "Config file should exist");

    // Check config file contents
    let contents = fs::read_to_string(&config_file).unwrap();
    assert!(contents.contains("email: test@example.com"));

    // Restore original HOME
    if let Some(home) = original_home {
        std::env::set_var("HOME", home);
    }
}

#[test]
fn test_init_command_existing_config() {
    let _guard = HOME_MUTEX.lock().unwrap();

    let temp_dir = TempDir::new().unwrap();
    let config_dir = temp_dir.path().join(".biovault");
    fs::create_dir_all(&config_dir).unwrap();

    // Create existing config
    let config_file = config_dir.join("config.yaml");
    fs::write(&config_file, "email: existing@example.com\n").unwrap();

    // Save original HOME and set temporary HOME
    let original_home = std::env::var("HOME").ok();
    std::env::set_var("HOME", temp_dir.path());

    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.env("HOME", temp_dir.path())
        .arg("init")
        .arg("new@example.com")
        .assert()
        .success()
        .stdout(predicate::str::contains("already exists"));

    // Check that original config is unchanged
    let contents = fs::read_to_string(&config_file).unwrap();
    assert!(contents.contains("email: existing@example.com"));
    assert!(!contents.contains("email: new@example.com"));

    // Restore original HOME
    if let Some(home) = original_home {
        std::env::set_var("HOME", home);
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
fn test_check_command() {
    let mut cmd = Command::cargo_bin("bv").unwrap();
    cmd.arg("check")
        .assert()
        .success()
        .stdout(predicate::str::contains("BioVault Dependency Check"))
        .stdout(predicate::str::contains("docker"))
        .stdout(predicate::str::contains("nextflow"));
    // Note: not checking for 'git' since deps.yaml was updated
}

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
