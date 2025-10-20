use super::check::DependencyConfig;
use crate::Result;
use anyhow::anyhow;
use std::collections::HashSet;
use std::env;
use std::io::{self, Write};
use std::process::{Command, ExitStatus, Stdio};
use std::thread;
use std::time::{Duration, Instant};

#[cfg(target_os = "macos")]
use std::fs;
#[cfg(target_os = "macos")]
use std::os::unix::fs::PermissionsExt;
#[cfg(target_os = "macos")]
use std::path::Path;
#[cfg(target_os = "macos")]
use std::sync::atomic::{AtomicBool, Ordering};
#[cfg(target_os = "macos")]
use std::sync::OnceLock;

use tracing::error;

const MAX_CAPTURED_OUTPUT: usize = 4000;

#[cfg(target_os = "macos")]
static SUDO_PRIMED: AtomicBool = AtomicBool::new(false);
#[cfg(target_os = "macos")]
static COMMAND_LINE_TOOLS_READY: AtomicBool = AtomicBool::new(false);

struct InstallCommandOutput {
    status: ExitStatus,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
}

impl InstallCommandOutput {
    fn from_output(output: std::process::Output) -> Self {
        Self {
            status: output.status,
            stdout: output.stdout,
            stderr: output.stderr,
        }
    }
}

fn truncate_output(text: &str) -> String {
    if text.len() > MAX_CAPTURED_OUTPUT {
        let mut truncated = text[..MAX_CAPTURED_OUTPUT].to_string();
        truncated.push_str("\n… (truncated)");
        truncated
    } else {
        text.to_string()
    }
}

fn shell_escape(segment: &str) -> String {
    if segment.is_empty() {
        "''".to_string()
    } else if segment.as_bytes().iter().all(|b| {
        matches!(b,
            b'0'..=b'9'
                | b'a'..=b'z'
                | b'A'..=b'Z'
                | b'.'
                | b'-'
                | b'_'
                | b'/'
        )
    }) {
        segment.to_string()
    } else {
        format!("'{}'", segment.replace("'", "'\\''"))
    }
}

fn skip_install_commands() -> bool {
    env::var("BIOVAULT_SKIP_INSTALLS")
        .map(|v| v != "0")
        .unwrap_or(false)
}

#[cfg(target_os = "macos")]
fn ensure_macos_command_line_tools_installed() -> Result<()> {
    if COMMAND_LINE_TOOLS_READY.load(Ordering::SeqCst) {
        if Command::new("xcode-select")
            .arg("-p")
            .status()
            .map(|status| status.success())
            .unwrap_or(false)
        {
            return Ok(());
        }

        COMMAND_LINE_TOOLS_READY.store(false, Ordering::SeqCst);
    }

    match Command::new("xcode-select").arg("-p").status() {
        Ok(status) if status.success() => {
            COMMAND_LINE_TOOLS_READY.store(true, Ordering::SeqCst);
            return Ok(());
        }
        Ok(_) => {}
        Err(err) => {
            return Err(anyhow!("Failed to check for macOS Command Line Tools: {}", err).into());
        }
    }

    eprintln!("🍎 macOS Command Line Tools missing – launching installer before proceeding",);

    let launch_result = Command::new("xcode-select").arg("--install").status();

    match launch_result {
        Ok(status) if status.success() => {
            eprintln!(
                "🍎 Command Line Tools installer launched. Complete the Apple dialog to continue."
            );
        }
        Ok(status) => {
            if status.code() == Some(1) {
                eprintln!(
                    "🍎 Command Line Tools installer already running. Complete the Apple dialog to continue."
                );
            } else {
                return Err(anyhow!(
                    "macOS Command Line Tools installer exited with status {:?}. Run 'xcode-select --install' manually, complete it, then rerun 'Install Missing'.",
                    status.code()
                )
                .into());
            }
        }
        Err(err) => {
            return Err(anyhow!(
                "Failed to launch the macOS Command Line Tools installer: {}. Run 'xcode-select --install' manually, complete it, then rerun 'Install Missing'.",
                err
            )
            .into());
        }
    }

    eprintln!(
        "🍎 Waiting for Command Line Tools installation to finish (this can take several minutes)..."
    );

    let start = Instant::now();
    let max_wait = Duration::from_secs(60 * 30); // 30 minutes maximum
    let poll_interval = Duration::from_secs(5);

    loop {
        if Command::new("xcode-select")
            .arg("-p")
            .status()
            .map(|status| status.success())
            .unwrap_or(false)
        {
            COMMAND_LINE_TOOLS_READY.store(true, Ordering::SeqCst);
            eprintln!("🍎 Command Line Tools detected. Continuing installation.");
            return Ok(());
        }

        if start.elapsed() >= max_wait {
            return Err(anyhow!(
                "macOS Command Line Tools installation did not complete within 30 minutes. Finish the Apple installer, then rerun 'Install Missing'."
            )
            .into());
        }

        thread::sleep(poll_interval);
    }
}

#[cfg(not(target_os = "macos"))]
fn ensure_macos_command_line_tools_installed() -> Result<()> {
    Ok(())
}

#[cfg(target_os = "macos")]
fn prime_macos_sudo_credentials(askpass_path: &str) -> io::Result<()> {
    if SUDO_PRIMED.load(Ordering::SeqCst) {
        if Command::new("/usr/bin/sudo")
            .arg("-n")
            .arg("true")
            .status()
            .map(|status| status.success())
            .unwrap_or(false)
        {
            return Ok(());
        }

        SUDO_PRIMED.store(false, Ordering::SeqCst);
    }

    let mut cmd = Command::new("/usr/bin/sudo");
    cmd.arg("-Av");
    cmd.env("SUDO_ASKPASS", askpass_path);
    cmd.env("HOMEBREW_SUDO_ASKPASS", askpass_path);

    match cmd.status() {
        Ok(status) if status.success() => {
            SUDO_PRIMED.store(true, Ordering::SeqCst);
            Ok(())
        }
        Ok(status) if status.code() == Some(1) => Err(io::Error::new(
            io::ErrorKind::Other,
            "Administrator authentication was cancelled.",
        )),
        Ok(status) => Err(io::Error::new(
            io::ErrorKind::Other,
            format!("sudo validation failed with status {:?}", status.code()),
        )),
        Err(err) => Err(err),
    }
}

#[cfg(not(target_os = "macos"))]
fn prime_macos_sudo_credentials(_: &str) -> io::Result<()> {
    Ok(())
}

fn run_install_command(
    system_type: &SystemType,
    command: &str,
) -> io::Result<InstallCommandOutput> {
    match system_type {
        SystemType::Windows => {
            if command.starts_with("winget") {
                let mut parts = command.split_whitespace();
                parts.next();
                Command::new("winget")
                    .args(parts)
                    .output()
                    .map(InstallCommandOutput::from_output)
            } else {
                Command::new("powershell")
                    .arg("-Command")
                    .arg(command)
                    .output()
                    .map(InstallCommandOutput::from_output)
            }
        }
        SystemType::MacOs => run_macos_install_command(command),
        _ => Command::new("sh")
            .arg("-c")
            .arg(command)
            .output()
            .map(InstallCommandOutput::from_output),
    }
}

#[cfg(target_os = "macos")]
fn run_macos_install_command(command: &str) -> io::Result<InstallCommandOutput> {
    if is_brew_command(command) {
        match ensure_macos_askpass_script()? {
            Some(askpass_path) => {
                let mut cmd = Command::new("sh");
                cmd.arg("-c").arg(command);
                prime_macos_sudo_credentials(&askpass_path)?;
                cmd.env("SUDO_ASKPASS", &askpass_path);
                cmd.env("HOMEBREW_SUDO_ASKPASS", &askpass_path);
                cmd.env("HOMEBREW_NO_AUTO_UPDATE", "1");

                if let Some(wrapper_path) = ensure_macos_sudo_wrapper()? {
                    let existing_path = env::var("PATH").unwrap_or_default();
                    let wrapper_dir = std::path::Path::new(&wrapper_path)
                        .parent()
                        .map(|p| p.display().to_string())
                        .unwrap_or_default();
                    let combined_path = if existing_path.is_empty() {
                        wrapper_dir
                    } else {
                        format!("{}:{}", wrapper_dir, existing_path)
                    };
                    cmd.env("PATH", combined_path);
                }

                return cmd.output().map(InstallCommandOutput::from_output);
            }
            None => {
                let prefixed = format!(
                    "PATH=\"/usr/local/bin:/opt/homebrew/bin:$PATH\"; export PATH; HOMEBREW_NO_AUTO_UPDATE=1 {}",
                    command.trim_start()
                );
                return run_macos_authorized_command(&prefixed);
            }
        }
    }

    Command::new("sh")
        .arg("-c")
        .arg(command)
        .output()
        .map(InstallCommandOutput::from_output)
}

#[cfg(target_os = "macos")]
fn escape_for_applescript(input: &str) -> String {
    input.replace("\\", "\\\\").replace("\"", "\\\"")
}

#[cfg(target_os = "macos")]
fn run_macos_authorized_command(command: &str) -> io::Result<InstallCommandOutput> {
    let applescript = format!(
        "do shell script \"{}\" with administrator privileges",
        escape_for_applescript(command.trim_start())
    );
    let output = Command::new("osascript")
        .arg("-e")
        .arg(applescript)
        .output()?;
    Ok(InstallCommandOutput::from_output(output))
}

#[cfg(not(target_os = "macos"))]
fn run_macos_install_command(command: &str) -> io::Result<InstallCommandOutput> {
    Command::new("sh")
        .arg("-c")
        .arg(command)
        .output()
        .map(InstallCommandOutput::from_output)
}

#[cfg(target_os = "macos")]
fn ensure_macos_askpass_script() -> io::Result<Option<String>> {
    static ASKPASS_PATH: OnceLock<Option<String>> = OnceLock::new();

    if let Some(existing) = ASKPASS_PATH.get() {
        return Ok(existing.clone());
    }

    if Path::new("/usr/libexec/ssh-askpass").exists() {
        let result = Some("/usr/libexec/ssh-askpass".to_string());
        let _ = ASKPASS_PATH.set(result.clone());
        return Ok(result);
    }

    let mut path = std::env::temp_dir();
    path.push("biovault_brew_askpass.sh");

    let script = r#"#!/bin/bash
result=$(/usr/bin/osascript <<'APPLESCRIPT'
on promptForPassword()
    tell application "System Events"
        activate
        set dialogText to "BioVault needs your administrator password to continue installing required tools. macOS is requesting this password on behalf of BioVault."
        repeat
            try
                set dialogResult to display dialog dialogText default answer "" with title "BioVault Installer" buttons {"Cancel", "Continue"} default button "Continue" with icon caution with hidden answer
                set userPassword to text returned of dialogResult
                if userPassword is not "" then
                    return userPassword
                end if
            on error number -128
                error number -128
            end try
        end repeat
    end tell
end promptForPassword

with timeout of 600 seconds
    promptForPassword()
end timeout
APPLESCRIPT
)
status=$?
if [ "$status" -ne 0 ]; then
    exit "$status"
fi
printf '%s\n' "$result"
"#;

    let result = match fs::write(&path, script) {
        Ok(()) => match fs::metadata(&path) {
            Ok(metadata) => {
                let mut perms = metadata.permissions();
                perms.set_mode(0o700);
                if let Err(err) = fs::set_permissions(&path, perms) {
                    eprintln!(
                        "⚠️  Failed to set askpass helper permissions at {}: {}",
                        path.display(),
                        err
                    );
                    None
                } else {
                    Some(path.to_string_lossy().to_string())
                }
            }
            Err(err) => {
                eprintln!(
                    "⚠️  Failed to read askpass helper metadata at {}: {}",
                    path.display(),
                    err
                );
                None
            }
        },
        Err(err) => {
            eprintln!(
                "⚠️  Failed to create askpass helper at {}: {}",
                path.display(),
                err
            );
            None
        }
    };

    let _ = ASKPASS_PATH.set(result.clone());
    Ok(result)
}

#[cfg(target_os = "macos")]
fn ensure_macos_sudo_wrapper() -> io::Result<Option<String>> {
    static WRAPPER_PATH: OnceLock<Option<String>> = OnceLock::new();

    if let Some(existing) = WRAPPER_PATH.get() {
        return Ok(existing.clone());
    }

    let mut path = std::env::temp_dir();
    path.push("biovault_sudo_wrapper.sh");

    let script = "#!/bin/bash\nexec /usr/bin/sudo -A \"$@\"\n";

    let result = match fs::write(&path, script) {
        Ok(()) => match fs::metadata(&path) {
            Ok(metadata) => {
                let mut perms = metadata.permissions();
                perms.set_mode(0o700);
                if let Err(err) = fs::set_permissions(&path, perms) {
                    eprintln!(
                        "⚠️  Failed to set sudo wrapper permissions at {}: {}",
                        path.display(),
                        err
                    );
                    None
                } else {
                    Some(path.to_string_lossy().to_string())
                }
            }
            Err(err) => {
                eprintln!(
                    "⚠️  Failed to read sudo wrapper metadata at {}: {}",
                    path.display(),
                    err
                );
                None
            }
        },
        Err(err) => {
            eprintln!(
                "⚠️  Failed to create sudo wrapper at {}: {}",
                path.display(),
                err
            );
            None
        }
    };

    let _ = WRAPPER_PATH.set(result.clone());
    Ok(result)
}

#[cfg(target_os = "macos")]
fn run_macos_sudo_command(args: &[&str]) -> io::Result<()> {
    let askpass = ensure_macos_askpass_script()?;
    if let Some(askpass_path) = askpass {
        let mut cmd = Command::new("/usr/bin/sudo");
        cmd.arg("-A").args(args);
        cmd.env("SUDO_ASKPASS", &askpass_path);
        cmd.env("HOMEBREW_SUDO_ASKPASS", &askpass_path);

        let status = cmd.status()?;
        if status.success() {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::Other,
                format!("sudo command failed: {:?}", args),
            ))
        }
    } else {
        let escaped = args
            .iter()
            .map(|arg| shell_escape(arg))
            .collect::<Vec<_>>()
            .join(" ");
        let command = format!("{}", escaped);
        let output = run_macos_authorized_command(&command)?;
        if output.status.success() {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::Other,
                format!("command failed: {:?}", args),
            ))
        }
    }
}

#[cfg(target_os = "macos")]
fn ensure_docker_cli_plugin_dir() -> io::Result<()> {
    let plugin_dir = Path::new("/usr/local/cli-plugins");

    if plugin_dir.exists() {
        return Ok(());
    }

    run_macos_sudo_command(&["/bin/mkdir", "-p", "/usr/local/cli-plugins"])?;

    if let Ok(username) = env::var("USER") {
        if !username.is_empty() {
            let ownership = format!("{}:staff", username);
            let _ =
                run_macos_sudo_command(&["/usr/sbin/chown", &ownership, "/usr/local/cli-plugins"]);
        }
    }

    let _ = run_macos_sudo_command(&["/bin/chmod", "755", "/usr/local/cli-plugins"]);

    Ok(())
}

fn is_brew_command(command: &str) -> bool {
    command
        .trim_start()
        .split_whitespace()
        .next()
        .map(|token| token == "brew" || token.ends_with("/brew"))
        .unwrap_or(false)
}

/// Install a single dependency and return the path to the installed binary (for library use)
pub async fn install_single_dependency(name: &str) -> Result<Option<String>> {
    // Load the deps.yaml file to get dependency information
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    // Find the dependency by name
    let dep = config
        .dependencies
        .iter()
        .find(|d| d.name == name)
        .ok_or_else(|| anyhow!("Dependency '{}' not found", name))?;

    // Detect system type
    let system_type = detect_system();
    let env_name = match system_type {
        SystemType::GoogleColab => "google_colab",
        SystemType::MacOs => "macos",
        SystemType::Ubuntu => "ubuntu",
        SystemType::ArchLinux => "arch",
        SystemType::Windows => "windows",
        SystemType::Unknown => {
            return Err(anyhow!(
                "System type not detected or not supported for automated installation"
            )
            .into());
        }
    };

    // Get environment-specific configuration
    let env_config = dep
        .environments
        .as_ref()
        .and_then(|envs| envs.get(env_name))
        .ok_or_else(|| anyhow!("No installation configuration for {} on this system", name))?;

    if env_config.skip {
        return Err(anyhow!(
            "Installation of {} is not supported on this system: {}",
            name,
            env_config
                .skip_reason
                .as_ref()
                .unwrap_or(&"Not available".to_string())
        )
        .into());
    }

    let install_commands = env_config
        .install_commands
        .as_ref()
        .ok_or_else(|| anyhow!("No installation commands available for {}", name))?;

    #[cfg(target_os = "macos")]
    if matches!(system_type, SystemType::MacOs)
        && install_commands
            .iter()
            .any(|cmd| cmd.trim_start().starts_with("brew "))
    {
        ensure_macos_command_line_tools_installed()?;
    }

    // Special handling for Homebrew on macOS
    let brew_cmd = if matches!(system_type, SystemType::MacOs) {
        find_brew_command()
    } else {
        None
    };

    #[cfg(target_os = "macos")]
    if matches!(system_type, SystemType::MacOs) {
        if let Some(ref brew) = brew_cmd {
            match Command::new(brew).arg("--version").status() {
                Ok(status) if status.success() => {}
                Ok(status) => {
                    return Err(anyhow!(
                        "Homebrew at '{}' is not working (brew --version exited with status {:?}). Reinstall Homebrew, then retry installing {}.",
                        brew,
                        status.code(),
                        name
                    )
                    .into());
                }
                Err(err) => {
                    return Err(anyhow!(
                        "Failed to execute '{} --version': {}. Reinstall Homebrew, then retry installing {}.",
                        brew,
                        err,
                        name
                    )
                    .into());
                }
            }
        }
    }

    #[cfg(target_os = "macos")]
    if matches!(system_type, SystemType::MacOs) && name == "docker" {
        ensure_docker_cli_plugin_dir()?;
    }

    // Execute installation commands
    for cmd in install_commands {
        let adjusted_cmd = if matches!(system_type, SystemType::MacOs) && cmd.starts_with("brew ") {
            if let Some(ref brew) = brew_cmd {
                cmd.replace("brew ", &format!("{} ", brew))
            } else {
                return Err(anyhow!("Homebrew not found. Please install Homebrew first.").into());
            }
        } else {
            cmd.clone()
        };

        let command_result = run_install_command(&system_type, &adjusted_cmd);

        match command_result {
            Ok(output) if output.status.success() => {
                if cfg!(debug_assertions) {
                    // In debug builds emit captured stdout/stderr to aid local testing
                    if !output.stdout.is_empty() {
                        println!("{}", String::from_utf8_lossy(&output.stdout).trim_end());
                    }
                    if !output.stderr.is_empty() {
                        eprintln!("{}", String::from_utf8_lossy(&output.stderr).trim_end());
                    }
                }
            }
            Ok(output) => {
                let stdout = String::from_utf8_lossy(&output.stdout);
                let stderr = String::from_utf8_lossy(&output.stderr);

                let stdout_trimmed = stdout.trim_end();
                let stderr_trimmed = stderr.trim_end();

                let mut details = Vec::new();
                if !stdout_trimmed.is_empty() {
                    details.push(format!("stdout:\n{}", truncate_output(stdout_trimmed)));
                }
                if !stderr_trimmed.is_empty() {
                    details.push(format!("stderr:\n{}", truncate_output(stderr_trimmed)));
                }
                if details.is_empty() {
                    details.push("(no output captured)".to_string());
                }

                let detail_text = details.join("\n\n");
                error!(
                    command = %adjusted_cmd,
                    stdout = %truncate_output(stdout_trimmed),
                    stderr = %truncate_output(stderr_trimmed),
                    "Installation command failed"
                );
                return Err(anyhow!(
                    "Failed to execute installation command: {}\n{}",
                    adjusted_cmd,
                    detail_text
                )
                .into());
            }
            Err(error) => {
                error!(
                    command = %adjusted_cmd,
                    ?error,
                    "Failed to spawn installation command"
                );
                return Err(anyhow!(
                    "Failed to spawn installation command: {} (error: {})",
                    adjusted_cmd,
                    error
                )
                .into());
            }
        }
    }

    // After successful installation, find the installed binary path
    let installed_path = match name {
        "java" => {
            // Special handling for Java on macOS (might be in brew path)
            if matches!(system_type, SystemType::MacOs) {
                if let Some(java_path) = check_java_in_brew_not_in_path() {
                    Some(format!("{}/java", java_path))
                } else {
                    which::which("java").ok().map(|p| p.display().to_string())
                }
            } else {
                which::which("java").ok().map(|p| p.display().to_string())
            }
        }
        _ => {
            // For other tools, try to find them in PATH
            which::which(name).ok().map(|p| p.display().to_string())
        }
    };

    Ok(installed_path)
}

pub async fn install_dependencies(names: &[String]) -> Result<Vec<(String, Option<String>)>> {
    #[cfg(target_os = "macos")]
    {
        if matches!(detect_system(), SystemType::MacOs) {
            ensure_macos_command_line_tools_installed()?;
        }
    }

    let mut results = Vec::new();
    let mut seen = HashSet::new();

    for name in names {
        if !seen.insert(name.clone()) {
            continue;
        }

        let installed = install_single_dependency(name).await?;
        results.push((name.clone(), installed));
    }

    Ok(results)
}

#[derive(Debug)]
enum SystemType {
    GoogleColab,
    MacOs,
    Ubuntu,
    ArchLinux,
    Windows,
    Unknown,
}

pub async fn execute() -> Result<()> {
    println!("BioVault Environment Setup");
    println!("==========================\n");

    let system_type = detect_system();

    match system_type {
        SystemType::GoogleColab => {
            println!("✓ Detected Google Colab environment");
            setup_google_colab().await?;
        }
        SystemType::MacOs => {
            println!("✓ Detected macOS environment");
            setup_macos().await?;
        }
        SystemType::Ubuntu => {
            println!("✓ Detected Ubuntu/Debian environment");
            setup_ubuntu().await?;
        }
        SystemType::ArchLinux => {
            println!("✓ Detected Arch Linux environment");
            setup_arch().await?;
        }
        SystemType::Windows => {
            println!("✓ Detected Windows environment");
            setup_windows().await?;
        }
        SystemType::Unknown => {
            println!("ℹ️  System type not detected or not supported for automated setup");
            println!("   This command currently supports:");
            println!("   - Google Colab");
            println!("   - macOS (Homebrew)");
            println!("   - Ubuntu/Debian (apt)");
            println!("   - Arch Linux (pacman)");
            println!("   - Windows (WinGet)");
            println!("\n   For manual setup, please run: bv check");
        }
    }

    Ok(())
}

fn detect_system() -> SystemType {
    let target_os = std::env::consts::OS;

    // Short-circuit for Windows before inspecting other environment hints
    if target_os == "windows" {
        return SystemType::Windows;
    }

    // Check for Google Colab environment variables (Colab runs on Linux)
    if is_google_colab() {
        return SystemType::GoogleColab;
    }

    // Detect macOS
    if target_os == "macos" {
        return SystemType::MacOs;
    }

    // Detect Linux distributions
    if target_os == "linux" {
        // Check for apt (Ubuntu/Debian)
        let has_apt = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v apt-get >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if has_apt {
            return SystemType::Ubuntu;
        }

        // Check for pacman (Arch Linux)
        let has_pacman = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v pacman >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if has_pacman {
            return SystemType::ArchLinux;
        }
    }

    SystemType::Unknown
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

async fn setup_google_colab() -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Google Colab setup commands");
        return Ok(());
    }

    println!("\nSetting up Google Colab environment...\n");

    // Load the deps.yaml file to get environment-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
        // Check if this dependency has google_colab environment config
        if let Some(environments) = &dep.environments {
            if let Some(env_config) = environments.get("google_colab") {
                if env_config.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_config
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_config.install_commands {
                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_succeeded = true;

                    for cmd in install_commands {
                        println!("   Running: {}", cmd);

                        // For Colab, we need to run these commands with sh -c
                        let output = Command::new("sh").arg("-c").arg(cmd).output();

                        match output {
                            Ok(output) => {
                                if output.status.success() {
                                    println!("   ✓ Command succeeded");
                                } else {
                                    println!("   ❌ Command failed");
                                    if !output.stderr.is_empty() {
                                        println!(
                                            "   Error: {}",
                                            String::from_utf8_lossy(&output.stderr)
                                        );
                                    }
                                    all_succeeded = false;
                                    break;
                                }
                            }
                            Err(e) => {
                                println!("   ❌ Failed to execute: {}", e);
                                all_succeeded = false;
                                break;
                            }
                        }
                    }

                    // Verify installation if verification command is provided
                    if all_succeeded {
                        if let Some(verify_cmd) = &env_config.verify_command {
                            print!("   Verifying installation... ");
                            let output = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::piped())
                                .stderr(Stdio::piped())
                                .output();

                            if let Ok(output) = output {
                                if output.status.success() {
                                    println!("✓");
                                    success_count += 1;
                                } else {
                                    println!("❌ Verification failed");
                                    fail_count += 1;
                                }
                            } else {
                                println!("❌ Could not verify");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    // Add PATH export instructions for Colab
    println!("📝 Final setup steps for Google Colab:\n");
    println!("   Add these lines to your notebook for persistence:");
    println!("   ```python");
    println!("   import os");
    println!("   os.environ['PATH'] = f\"/usr/local/bin:{{os.environ['PATH']}}\"");
    println!("   ```");
    println!();
    println!("   Or in a shell cell:");
    println!("   ```bash");
    println!("   !export PATH=\"/usr/local/bin:$PATH\"");
    println!("   ```");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
    } else {
        println!("\n✅ Setup completed successfully!");
        println!("   Run 'bv check' to verify all dependencies.");
    }

    Ok(())
}

async fn setup_macos() -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping macOS setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up macOS environment...\n");

    // Check if we're in CI mode (non-interactive)
    let is_ci = env::var("CI").is_ok() || env::var("GITHUB_ACTIONS").is_ok();

    // Check for Homebrew
    let brew_in_path = which::which("brew").is_ok();
    let mut brew_path = None;

    if !brew_in_path {
        // Check common Homebrew locations even if not in PATH
        let common_brew_paths = vec![
            "/opt/homebrew/bin/brew", // Apple Silicon
            "/usr/local/bin/brew",    // Intel Mac
        ];

        for path in &common_brew_paths {
            if std::path::Path::new(path).exists() {
                brew_path = Some(path.to_string());
                println!("📦 Found Homebrew at {} (not in PATH)", path);
                break;
            }
        }
    }

    // Install Homebrew if not found
    if !brew_in_path && brew_path.is_none() {
        println!("📦 Homebrew not found. Would you like to install it? [Y/n]: ");

        if is_ci {
            println!("   CI mode: Skipping Homebrew installation.");
            println!("   Please ensure Homebrew is pre-installed in CI environment.");
            return Ok(());
        } else {
            io::stdout().flush()?;
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let answer = input.trim().to_lowercase();

            if answer.is_empty() || answer == "y" || answer == "yes" {
                println!("Installing Homebrew...");
                let install_cmd = "/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"";
                let status = Command::new("sh").arg("-c").arg(install_cmd).status()?;

                if status.success() {
                    println!("✓ Homebrew installed successfully!");
                    // Detect where Homebrew was installed
                    if std::path::Path::new("/opt/homebrew/bin/brew").exists() {
                        brew_path = Some("/opt/homebrew/bin/brew".to_string());
                    } else if std::path::Path::new("/usr/local/bin/brew").exists() {
                        brew_path = Some("/usr/local/bin/brew".to_string());
                    }
                } else {
                    println!("❌ Homebrew installation failed.");
                    println!("Please install manually from: https://brew.sh");
                    return Ok(());
                }
            } else {
                println!("Skipping Homebrew installation.");
                println!("Please install Homebrew manually from: https://brew.sh");
                println!("Then re-run: bv setup");
                return Ok(());
            }
        }
    }

    // Use the brew command (either from PATH or specific path)
    let brew_cmd = if brew_in_path {
        "brew".to_string()
    } else if let Some(ref bp) = brew_path {
        bp.clone()
    } else {
        "brew".to_string() // Fallback
    };

    // Load deps.yaml and execute macOS-specific commands
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    // Use the brew command for installations
    for dep in &config.dependencies {
        if let Some(environments) = &dep.environments {
            if let Some(env_config) = environments.get("macos") {
                if env_config.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_config
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on macOS".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_config.install_commands {
                    // Decide if install is necessary
                    let mut need_install = true;

                    // If verify_command is available, try it first
                    if let Some(verify_cmd) = &env_config.verify_command {
                        let verified = Command::new("sh")
                            .arg("-c")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    } else {
                        // Fallback to which for simple presence
                        if which::which(&dep.name).is_ok() {
                            need_install = false;
                        }
                    }

                    // For Java, also enforce min_version if specified
                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_succeeded = true;

                    for cmd in install_commands {
                        // Replace 'brew' with the actual brew path if needed
                        let adjusted_cmd = if !brew_in_path && cmd.starts_with("brew ") {
                            cmd.replace("brew ", &format!("{} ", brew_cmd))
                        } else {
                            cmd.clone()
                        };

                        println!("   Running: {}", adjusted_cmd);
                        let output = Command::new("sh")
                            .arg("-c")
                            .arg(&adjusted_cmd)
                            .stdout(Stdio::piped())
                            .stderr(Stdio::piped())
                            .output();
                        match output {
                            Ok(output) => {
                                if output.status.success() {
                                    println!("   ✓ Command succeeded");
                                } else {
                                    println!("   ❌ Command failed");
                                    if !output.stderr.is_empty() {
                                        println!(
                                            "   Error: {}",
                                            String::from_utf8_lossy(&output.stderr)
                                        );
                                    }
                                    all_succeeded = false;
                                    break;
                                }
                            }
                            Err(e) => {
                                println!("   ❌ Failed to execute: {}", e);
                                all_succeeded = false;
                                break;
                            }
                        }
                    }

                    if all_succeeded {
                        // Special handling for Java: add brew path before verification
                        if dep.name == "java" {
                            if let Some(java_path) = check_java_in_brew_not_in_path() {
                                // Add the brew Java path to current environment for verification
                                let current_path = env::var("PATH").unwrap_or_default();
                                env::set_var("PATH", format!("{}:{}", java_path, current_path));
                            }
                        }

                        if let Some(verify_cmd) = &env_config.verify_command {
                            print!("   Verifying installation... ");
                            let verify_result =
                                Command::new("sh").arg("-c").arg(verify_cmd).output();

                            match verify_result {
                                Ok(output) if output.status.success() => {
                                    println!("✓");
                                    success_count += 1;
                                }
                                Ok(_)
                                    if dep.name == "docker" && docker_desktop_installed_macos() =>
                                {
                                    println!(
                                        "⚠️  Docker Desktop installed – launch it once to finish CLI setup"
                                    );
                                    success_count += 1;
                                }
                                Ok(output) => {
                                    println!("❌ Verification failed");
                                    if !output.stdout.is_empty() {
                                        println!(
                                            "   stdout: {}",
                                            String::from_utf8_lossy(&output.stdout).trim_end()
                                        );
                                    }
                                    if !output.stderr.is_empty() {
                                        println!(
                                            "   stderr: {}",
                                            String::from_utf8_lossy(&output.stderr).trim_end()
                                        );
                                    }
                                    fail_count += 1;
                                }
                                Err(err)
                                    if dep.name == "docker" && docker_desktop_installed_macos() =>
                                {
                                    println!(
                                        "⚠️  Docker Desktop installed – launch it once to finish CLI setup"
                                    );
                                    println!("   (could not run verification command: {})", err);
                                    success_count += 1;
                                }
                                Err(err) => {
                                    println!(
                                        "❌ Could not verify (failed to run '{}'): {}",
                                        verify_cmd, err
                                    );
                                    fail_count += 1;
                                }
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    // After all installations, check if Java needs PATH configuration
    // This handles the case where Java was already installed but not in PATH
    check_and_configure_java_path(is_ci).await?;

    println!("\nNotes:");
    println!("- If this is your first time installing Docker Desktop, open it once to finish setup and grant permissions.");

    // SyftBox info for manual setup later (we installed in setup-only mode)
    println!("\nSyftBox:");
    print_syftbox_instructions();

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

fn print_syftbox_instructions() {
    // Best-effort arch hint for user
    let arch = match std::env::consts::ARCH {
        "aarch64" => "arm64 (Apple Silicon)",
        "x86_64" => "x86_64 (Intel)",
        other => other,
    };
    println!(
        "Get the latest SyftBox for macOS ({}):\n  https://github.com/OpenMined/syftbox/releases/latest",
        arch
    );
    println!("After downloading, ensure the 'syftbox' binary is on your PATH (e.g., move to /usr/local/bin and chmod +x).");
}

#[cfg(target_os = "macos")]
fn docker_desktop_installed_macos() -> bool {
    if std::path::Path::new("/Applications/Docker.app").exists() {
        return true;
    }

    if let Some(home_dir) = dirs::home_dir() {
        let user_app = home_dir.join("Applications").join("Docker.app");
        return user_app.exists();
    }

    false
}

#[cfg(not(target_os = "macos"))]
fn docker_desktop_installed_macos() -> bool {
    false
}

// Minimal java version detection to respect min_version in deps.yaml
fn java_major_version() -> Option<u32> {
    let out = Command::new("java").arg("-version").output().ok()?;
    let text = String::from_utf8_lossy(&out.stderr);
    parse_java_version(&text)
}

fn parse_java_version(output: &str) -> Option<u32> {
    for line in output.lines() {
        if line.contains("version") {
            if let Some(start) = line.find('"') {
                if let Some(end) = line[start + 1..].find('"') {
                    let version_str = &line[start + 1..start + 1 + end];
                    if let Some(stripped) = version_str.strip_prefix("1.") {
                        if let Some(dot_pos) = stripped.find('.') {
                            if let Ok(v) = stripped[..dot_pos].parse::<u32>() {
                                return Some(v);
                            }
                        }
                    } else {
                        let major_part = version_str.split('.').next().unwrap_or(version_str);
                        if let Ok(v) = major_part.parse::<u32>() {
                            return Some(v);
                        }
                    }
                }
            }
        }
    }
    None
}

async fn setup_ubuntu() -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Ubuntu/Debian setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up Ubuntu/Debian environment...\n");

    // Ensure apt-get exists
    let apt_exists = Command::new("sh")
        .arg("-c")
        .arg("command -v apt-get >/dev/null 2>&1")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !apt_exists {
        println!("❌ apt-get not found. This setup targets Ubuntu/Debian-based systems.");
        println!("Please ensure you're on an apt-based distribution.");
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("ubuntu") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Ubuntu".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("sh")
                            .arg("-c")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    } else if which::which(&dep.name).is_ok() {
                        need_install = false;
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        // For apt commands on CI, we may need to run with sudo
                        let cmd_to_run = if cmd.starts_with("apt-get") && !cmd.starts_with("sudo") {
                            format!("sudo {}", cmd)
                        } else {
                            cmd.clone()
                        };

                        let status = Command::new("sh").arg("-c").arg(&cmd_to_run).status();
                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!("- For Docker on Ubuntu, you may need to add your user to the docker group: 'sudo usermod -aG docker $USER' and re-login.");
    println!("- If Docker service is not running: 'sudo systemctl start docker'");

    println!("\nSyftBox:");
    println!("The installer has been invoked in setup-only mode if needed.");
    println!("If you want to set up later: syftbox login; syftbox");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

async fn setup_arch() -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Arch Linux setup commands");
        return Ok(());
    }

    use super::check::DependencyConfig;
    use std::process::Command;

    println!("\nSetting up Arch Linux environment...\n");

    // Ensure pacman exists
    let pacman_exists = Command::new("sh")
        .arg("-c")
        .arg("command -v pacman >/dev/null 2>&1")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !pacman_exists {
        println!("❌ pacman not found. This setup targets Arch Linux.");
        println!("Please ensure you're on Arch/Manjaro with pacman available.");
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("arch") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Arch".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("sh")
                            .arg("-c")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    } else if which::which(&dep.name).is_ok() {
                        need_install = false;
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        let status = Command::new("sh").arg("-c").arg(cmd).status();
                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("sh")
                                .arg("-c")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!("- For Docker on Arch, you may need to enable and start the daemon: 'sudo systemctl enable --now docker' and add your user to the docker group.");

    println!("\nSyftBox:");
    println!("The installer has been invoked in setup-only mode if needed.");
    println!("If you want to set up later: syftbox login; syftbox");

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

async fn setup_windows() -> Result<()> {
    if skip_install_commands() {
        println!("(test mode) Skipping Windows setup commands");
        return Ok(());
    }

    println!("\nSetting up Windows environment...\n");

    // Check for WinGet availability
    let winget_exists = Command::new("winget")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    // Check for Chocolatey availability (fallback)
    let choco_exists = Command::new("choco")
        .arg("-v")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if winget_exists {
        println!("✓ WinGet found");
    } else if choco_exists {
        println!("❌ WinGet not found. Using Chocolatey fallback.");
    } else {
        println!("❌ Neither WinGet nor Chocolatey found.");
        println!("Automated installation is unavailable on this system.");
        println!("\nTo install WinGet:");
        println!("1. Update Windows to the latest version (WinGet comes with modern Windows)");
        println!("2. Or install from Microsoft Store: 'App Installer'");
        println!("3. Or download from: https://github.com/microsoft/winget-cli/releases");
        println!("\nAlternatively, install Chocolatey from https://chocolatey.org/install");
        print_windows_manual_instructions();
        return Ok(());
    }

    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;

    let mut success_count = 0;
    let mut skip_count = 0;
    let mut fail_count = 0;

    for dep in &config.dependencies {
        if let Some(envs) = &dep.environments {
            if let Some(env_cfg) = envs.get("windows") {
                if env_cfg.skip {
                    println!(
                        "⏭️  Skipping {}: {}",
                        dep.name,
                        env_cfg
                            .skip_reason
                            .as_ref()
                            .unwrap_or(&"Not needed on Windows".to_string())
                    );
                    skip_count += 1;
                    continue;
                }

                if let Some(install_commands) = &env_cfg.install_commands {
                    // Determine if installation is required
                    let mut need_install = true;
                    if let Some(verify_cmd) = &env_cfg.verify_command {
                        let verified = Command::new("powershell")
                            .arg("-Command")
                            .arg(verify_cmd)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .status()
                            .map(|s| s.success())
                            .unwrap_or(false);
                        if verified {
                            need_install = false;
                        }
                    }

                    if dep.name == "java" {
                        if let Some(min_v) = dep.min_version {
                            if let Some(current) = java_major_version() {
                                if current >= min_v {
                                    println!("   Java version {} already meets minimum requirement of {}", current, min_v);
                                    need_install = false;
                                }
                            }
                        }
                    }

                    if !need_install {
                        println!("✓ {} already meets requirements. Skipping.", dep.name);
                        skip_count += 1;
                        println!();
                        continue;
                    }

                    println!("📦 Installing {}...", dep.name);
                    println!("   {}", dep.description);

                    let mut all_ok = true;
                    for cmd in install_commands {
                        println!("   Running: {}", cmd);
                        let status = if cmd.starts_with("winget") {
                            if winget_exists {
                                Command::new("winget")
                                    .args(cmd.split_whitespace().skip(1))
                                    .status()
                            } else {
                                // Chocolatey fallback for common packages
                                let mut parts = cmd.split_whitespace();
                                let _ = parts.next(); // winget
                                let _ = parts.next(); // install
                                let pkg = parts.next().unwrap_or("");
                                let choco_pkg = map_winget_pkg_to_choco(pkg);
                                if choco_pkg.is_empty() {
                                    Err(std::io::Error::other("No Chocolatey mapping for package"))
                                } else {
                                    Command::new("choco")
                                        .arg("install")
                                        .arg(choco_pkg)
                                        .arg("-y")
                                        .status()
                                }
                            }
                        } else {
                            Command::new("powershell").arg("-Command").arg(cmd).status()
                        };

                        match status {
                            Ok(s) if s.success() => println!("   ✓ Command succeeded"),
                            Ok(_) | Err(_) => {
                                println!("   ❌ Command failed");
                                all_ok = false;
                                break;
                            }
                        }
                    }

                    if all_ok {
                        if let Some(verify_cmd) = &env_cfg.verify_command {
                            print!("   Verifying installation... ");
                            let ok = Command::new("powershell")
                                .arg("-Command")
                                .arg(verify_cmd)
                                .stdout(Stdio::null())
                                .stderr(Stdio::null())
                                .status()
                                .map(|s| s.success())
                                .unwrap_or(false);
                            if ok {
                                println!("✓");
                                success_count += 1;
                            } else {
                                println!("❌ Verification failed");
                                fail_count += 1;
                            }
                        } else {
                            success_count += 1;
                        }
                    } else {
                        fail_count += 1;
                    }

                    println!();
                }
            }
        }
    }

    println!("\nNotes:");
    println!(
        "- You may need to restart your terminal/PowerShell after installation to update PATH"
    );
    println!("- For Docker on Windows, Docker Desktop is required and may need manual setup");

    print_windows_manual_instructions();

    println!("\n==========================");
    println!("Setup Summary:");
    println!("  ✓ Installed: {}", success_count);
    println!("  ⏭️  Skipped: {}", skip_count);
    if fail_count > 0 {
        println!("  ❌ Failed: {}", fail_count);
        println!("\n⚠️  Some installations failed. Please check the errors above.");
        return Err(anyhow!("Some installations failed").into());
    } else {
        println!(
            "\n✅ Setup completed successfully!\n   Run 'bv check' to verify all dependencies."
        );
    }

    Ok(())
}

fn print_windows_manual_instructions() {
    println!("\nManual Installation Options:");
    println!(
        "Java 17+: Download from https://openjdk.org/ or use 'winget install Microsoft.OpenJDK'"
    );
    println!(
        "Docker: Download Docker Desktop from https://www.docker.com/products/docker-desktop/"
    );
    println!("Nextflow: Use WSL or Docker - no native Windows support");
    println!("SyftBox: Download from https://github.com/OpenMined/syftbox/releases/latest");
    println!("UV: Run 'winget install --id=astral-sh.uv -e' or use PowerShell installer from https://docs.astral.sh/uv/");
}

// Map common WinGet package IDs to Chocolatey package names for fallback
fn map_winget_pkg_to_choco(pkg: &str) -> &'static str {
    match pkg.to_ascii_lowercase().as_str() {
        // Java/OpenJDK
        // WinGet: Microsoft.OpenJDK => Chocolatey: openjdk (generic)
        "microsoft.openjdk" => "openjdk",
        // UV - Python package installer
        // WinGet: astral-sh.uv => Chocolatey doesn't have UV yet, so we return empty
        "astral-sh.uv" => "",
        // Add other mappings here as needed
        _ => "",
    }
}

async fn check_and_configure_java_path(is_ci: bool) -> Result<()> {
    // Check if Java is already in PATH
    if which::which("java").is_ok() {
        return Ok(());
    }

    // Check if Java is installed via brew but not in PATH
    let java_brew_path = check_java_in_brew_not_in_path();

    if let Some(brew_path) = java_brew_path {
        println!("\n⚠️  Java is installed via Homebrew but not in your PATH.");
        println!("   Location: {}", brew_path);

        let shell = env::var("SHELL").unwrap_or_else(|_| "/bin/zsh".to_string());
        let shell_config = if shell.contains("zsh") {
            format!(
                "{}/.zshrc",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        } else if shell.contains("bash") {
            format!(
                "{}/.bash_profile",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        } else {
            format!(
                "{}/.profile",
                env::var("HOME").unwrap_or_else(|_| "~".to_string())
            )
        };

        let export_line = format!("export PATH=\"{}:$PATH\"", brew_path);

        if is_ci {
            // In CI mode, automatically add to PATH configuration
            println!("   CI mode: Automatically configuring PATH...");

            // Add to the shell config file
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&shell_config)?;
            writeln!(file, "\n# Added by BioVault setup")?;
            writeln!(file, "{}", export_line)?;

            println!("   ✓ Added to {}", shell_config);
            println!("   Note: You'll need to restart your shell or run 'source {}' for changes to take effect.", shell_config);
        } else {
            // Interactive mode - prompt the user
            println!("\n   Would you like to automatically add Java to your PATH? [Y/n]: ");
            io::stdout().flush()?;

            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let answer = input.trim().to_lowercase();

            if answer.is_empty() || answer == "y" || answer == "yes" {
                // Add to the shell config file
                let mut file = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&shell_config)?;
                writeln!(file, "\n# Added by BioVault setup")?;
                writeln!(file, "{}", export_line)?;

                println!("   ✓ Added to {}", shell_config);
                println!("   Note: You'll need to restart your shell or run 'source {}' for changes to take effect.", shell_config);
            } else {
                println!("   Skipped PATH configuration.");
                println!("   To manually add Java to your PATH, run:");
                println!("     echo '{}' >> {}", export_line, shell_config);
                println!("     source {}", shell_config);
            }
        }
    }

    Ok(())
}
#[cfg(target_os = "macos")]
fn check_java_in_brew_not_in_path() -> Option<String> {
    // Find brew command (in PATH or common locations)
    let brew_cmd = find_brew_command();
    brew_cmd.as_ref()?;
    let brew_cmd = brew_cmd.unwrap();

    // Check if Java/OpenJDK is installed via brew
    let output = Command::new(&brew_cmd)
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
    let prefix_output = Command::new(&brew_cmd)
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

#[cfg(not(target_os = "macos"))]
fn check_java_in_brew_not_in_path() -> Option<String> {
    None
}
#[cfg(target_os = "macos")]
fn find_brew_command() -> Option<String> {
    // First check if brew is in PATH
    if which::which("brew").is_ok() {
        return Some("brew".to_string());
    }

    // Check common locations
    let common_brew_paths = vec![
        "/opt/homebrew/bin/brew", // Apple Silicon
        "/usr/local/bin/brew",    // Intel Mac
    ];

    for path in &common_brew_paths {
        if std::path::Path::new(path).exists() {
            return Some(path.to_string());
        }
    }

    None
}

#[cfg(not(target_os = "macos"))]
fn find_brew_command() -> Option<String> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;

    struct SkipInstallGuard(Option<String>);

    impl SkipInstallGuard {
        fn new() -> Self {
            let previous = env::var("BIOVAULT_SKIP_INSTALLS").ok();
            env::set_var("BIOVAULT_SKIP_INSTALLS", "1");
            Self(previous)
        }
    }

    impl Drop for SkipInstallGuard {
        fn drop(&mut self) {
            if let Some(ref value) = self.0 {
                env::set_var("BIOVAULT_SKIP_INSTALLS", value);
            } else {
                env::remove_var("BIOVAULT_SKIP_INSTALLS");
            }
        }
    }

    #[test]
    #[serial_test::serial]
    fn skip_install_commands_env_behavior() {
        // Save original state to restore later
        let original = env::var("BIOVAULT_SKIP_INSTALLS").ok();

        // Test 1: When not set, should return false
        env::remove_var("BIOVAULT_SKIP_INSTALLS");
        assert!(
            !super::skip_install_commands(),
            "Expected false when env var not set"
        );

        // Test 2: When set to "1", should return true
        env::set_var("BIOVAULT_SKIP_INSTALLS", "1");
        assert!(
            super::skip_install_commands(),
            "Expected true when env var set to '1'"
        );

        // Test 3: When set to "0", should return false
        env::set_var("BIOVAULT_SKIP_INSTALLS", "0");
        let result = super::skip_install_commands();
        let actual_value =
            env::var("BIOVAULT_SKIP_INSTALLS").unwrap_or_else(|_| "NOT_SET".to_string());
        assert!(
            !result,
            "Expected false when env var set to '0', but got true. Actual env value: '{}'",
            actual_value
        );

        // Restore original state
        if let Some(val) = original {
            env::set_var("BIOVAULT_SKIP_INSTALLS", val);
        } else {
            env::remove_var("BIOVAULT_SKIP_INSTALLS");
        }
    }

    #[test]
    fn java_parse_various_formats() {
        let cases = [
            ("openjdk version \"17.0.2\" 2022-01-18", Some(17)),
            ("java version \"1.8.0_321\"", Some(8)),
            ("openjdk version \"11.0.14\" 2022-01-18", Some(11)),
            ("java version \"21\"", Some(21)),
            ("garbage", None),
        ];
        for (s, want) in cases {
            assert_eq!(parse_java_version(s), want);
        }
    }

    #[test]
    #[serial_test::serial]
    fn google_colab_detection_via_env() {
        // Ensure variable not set
        std::env::remove_var("COLAB_RELEASE_TAG");
        for (k, _) in std::env::vars() {
            if k.starts_with("COLAB_") {
                std::env::remove_var(k);
            }
        }
        assert!(!is_google_colab());
        // Set specific var and detect
        std::env::set_var("COLAB_RELEASE_TAG", "test");
        assert!(is_google_colab());
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    #[serial_test::serial]
    #[cfg_attr(target_os = "windows", ignore = "Colab detection is Linux-specific")]
    fn detect_system_prefers_colab_env() {
        // Clear any existing COLAB_* variables first
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in &keys {
            std::env::remove_var(k);
        }

        // Force Colab-like environment
        std::env::set_var("COLAB_RELEASE_TAG", "1");
        match detect_system() {
            SystemType::GoogleColab => {}
            other => panic!("expected GoogleColab, got {:?}", other),
        }

        // Clean up
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    fn print_syftbox_instructions_runs() {
        // Just ensure it doesn't panic; covers simple printing logic
        super::print_syftbox_instructions();
    }

    #[test]
    #[serial_test::serial]
    fn is_google_colab_detects_prefix_env() {
        std::env::remove_var("COLAB_RELEASE_TAG");
        std::env::set_var("COLAB_FOO", "1");
        assert!(super::is_google_colab());
        std::env::remove_var("COLAB_FOO");
    }

    #[test]
    #[serial_test::serial]
    #[cfg(target_os = "macos")]
    fn detect_system_reports_macos_on_macos() {
        // Ensure no COLAB_* noise affects detection
        std::env::remove_var("COLAB_RELEASE_TAG");
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in keys {
            std::env::remove_var(k);
        }
        match super::detect_system() {
            super::SystemType::MacOs => {}
            other => panic!("expected MacOs, got {:?}", other),
        }
    }

    #[tokio::test]
    #[cfg(target_os = "macos")]
    async fn setup_ubuntu_returns_ok_when_apt_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_ubuntu().await.unwrap();
    }

    #[tokio::test]
    #[cfg(target_os = "macos")]
    async fn setup_arch_returns_ok_when_pacman_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_arch().await.unwrap();
    }

    #[test]
    fn winget_to_choco_mapping() {
        assert_eq!(
            super::map_winget_pkg_to_choco("Microsoft.OpenJDK"),
            "openjdk"
        );
        // Unknown returns empty mapping
        assert_eq!(super::map_winget_pkg_to_choco("Unknown.Package"), "");
    }

    #[tokio::test]
    async fn setup_google_colab_runs_without_panic() {
        let _guard = SkipInstallGuard::new();
        super::setup_google_colab().await.unwrap();
    }

    #[tokio::test]
    async fn setup_macos_returns_ok_without_brew() {
        let _guard = SkipInstallGuard::new();
        // Only run when brew is not available; otherwise skip to avoid invoking installs
        let brew_exists = std::process::Command::new("sh")
            .arg("-c")
            .arg("command -v brew >/dev/null 2>&1")
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !brew_exists {
            super::setup_macos().await.unwrap();
        }
    }

    #[tokio::test]
    #[cfg(target_os = "windows")]
    #[cfg_attr(
        not(feature = "e2e-tests"),
        ignore = "runs installer commands; e2e-only"
    )]
    async fn setup_windows_returns_ok_when_tools_missing() {
        let _guard = SkipInstallGuard::new();
        super::setup_windows().await.unwrap();
    }

    #[tokio::test]
    #[serial_test::serial]
    #[cfg_attr(
        target_os = "windows",
        ignore = "Colab execution path installs Linux tools"
    )]
    async fn setup_execute_colab_branch() {
        let _guard = SkipInstallGuard::new();
        std::env::set_var("COLAB_RELEASE_TAG", "1");
        super::execute().await.unwrap();
        std::env::remove_var("COLAB_RELEASE_TAG");
    }

    #[test]
    fn print_windows_manual_instructions_runs() {
        super::print_windows_manual_instructions();
    }

    #[test]
    fn test_system_type_debug() {
        let s = format!("{:?}", SystemType::MacOs);
        assert_eq!(s, "MacOs");
        let s2 = format!("{:?}", SystemType::GoogleColab);
        assert_eq!(s2, "GoogleColab");
    }

    #[test]
    #[serial_test::serial]
    fn test_detect_system_on_windows() {
        if cfg!(target_os = "windows") {
            match detect_system() {
                SystemType::Windows => {}
                _ => panic!("Expected Windows on windows platform"),
            }
        }
    }

    #[test]
    #[serial_test::serial]
    fn test_skip_install_commands_not_set() {
        std::env::remove_var("BIOVAULT_SKIP_INSTALLS");
        // Just verify it returns a bool without panicking
        let _result = skip_install_commands();
        // Function completes without panic - test passes
    }

    #[test]
    fn test_parse_java_version_edge_cases() {
        assert_eq!(parse_java_version(""), None);
        assert_eq!(parse_java_version("no version here"), None);
        // Just test that it doesn't panic on weird input
        let _ = parse_java_version("version 999");
    }

    #[test]
    fn test_map_winget_pkg_to_choco_all_mappings() {
        // Only Microsoft.OpenJDK is mapped
        assert_eq!(map_winget_pkg_to_choco("Microsoft.OpenJDK"), "openjdk");
        assert_eq!(map_winget_pkg_to_choco("microsoft.openjdk"), "openjdk");
        // Others return empty
        assert_eq!(map_winget_pkg_to_choco("Git.Git"), "");
        assert_eq!(map_winget_pkg_to_choco("RandomPackage"), "");
    }

    #[test]
    fn test_is_google_colab_without_env() {
        std::env::remove_var("COLAB_RELEASE_TAG");
        let keys: Vec<String> = std::env::vars()
            .filter(|(k, _)| k.starts_with("COLAB_"))
            .map(|(k, _)| k)
            .collect();
        for k in keys {
            std::env::remove_var(&k);
        }
        assert!(!is_google_colab());
    }

    // NOTE: Do NOT add unit tests for execute() - it runs actual installation commands
    // and should only be tested via e2e/integration tests.
    // The execute() function performs real system operations (downloads, installs, etc.)
    // which are not suitable for unit tests.
}
