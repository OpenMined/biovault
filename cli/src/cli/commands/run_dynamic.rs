use super::run::execute_with_logging;
use crate::error::Result;
use crate::project_spec::ProjectSpec;
use anyhow::Context;
use chrono::Local;
use colored::Colorize;
use serde_json::{json, Value as JsonValue};
use std::collections::{BTreeSet, HashMap};
use std::ffi::OsStr;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn append_desktop_log(message: &str) {
    if let Ok(path) = std::env::var("BIOVAULT_DESKTOP_LOG_FILE") {
        if path.is_empty() {
            return;
        }
        let path = PathBuf::from(path);
        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(&path) {
            let timestamp = Local::now().format("%Y-%m-%dT%H:%M:%S%:z");
            let _ = writeln!(file, "[{}][INFO] {}", timestamp, message);
        }
    }
}

fn shell_quote(value: &OsStr) -> String {
    let s = value.to_string_lossy();
    if s.is_empty() {
        return "''".to_string();
    }
    if s.chars()
        .all(|c| c.is_ascii_alphanumeric() || "-_./:@".contains(c))
    {
        return s.into_owned();
    }
    let escaped = s.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

fn format_command(cmd: &Command) -> String {
    let mut parts = Vec::new();
    parts.push(shell_quote(cmd.get_program()));
    for arg in cmd.get_args() {
        parts.push(shell_quote(arg));
    }
    parts.join(" ")
}

fn bundled_env_var(name: &str) -> Option<&'static str> {
    match name {
        "java" => Some("BIOVAULT_BUNDLED_JAVA"),
        "nextflow" => Some("BIOVAULT_BUNDLED_NEXTFLOW"),
        "uv" => Some("BIOVAULT_BUNDLED_UV"),
        "syftbox" => Some("SYFTBOX_BINARY"),
        _ => None,
    }
}

fn resolve_binary_path(cfg: Option<&crate::config::Config>, name: &str) -> Option<String> {
    if let Some(cfg) = cfg {
        if let Some(path) = cfg.get_binary_path(name) {
            if !path.is_empty() {
                return Some(path);
            }
        }
    }

    if let Some(env_key) = bundled_env_var(name) {
        if let Ok(env_path) = std::env::var(env_key) {
            let trimmed = env_path.trim();
            if !trimmed.is_empty() {
                return Some(trimmed.to_string());
            }
        }
    }

    None
}

fn build_augmented_path(cfg: Option<&crate::config::Config>) -> Option<String> {
    let mut entries = BTreeSet::new();
    for key in ["nextflow", "java", "docker"] {
        if let Some(bin_path) = resolve_binary_path(cfg, key) {
            if bin_path.is_empty() {
                continue;
            }
            if let Some(parent) = Path::new(&bin_path).parent() {
                entries.insert(parent.to_path_buf());
            }
        }
    }

    if entries.is_empty() {
        return None;
    }

    let mut paths: Vec<PathBuf> = entries.into_iter().collect();
    if let Some(existing) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&existing));
    }

    std::env::join_paths(paths)
        .ok()
        .and_then(|joined| joined.into_string().ok())
}

#[derive(Debug, Clone, Copy)]
pub struct RunSettings {
    /// Whether Docker must be running before launching Nextflow.
    pub require_docker: bool,
}

impl Default for RunSettings {
    fn default() -> Self {
        Self {
            // Current templates assume Docker; keep default strict.
            require_docker: true,
        }
    }
}

/// Check if we should use Docker to run Nextflow (Windows only)
fn should_use_docker_for_nextflow() -> bool {
    cfg!(target_os = "windows")
}

/// Convert a Windows path to a Docker-compatible path for volume mounts
/// e.g., C:\Users\foo -> /c/Users/foo
#[cfg(target_os = "windows")]
fn windows_path_to_docker(path: &Path) -> String {
    let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    let path_str = canonical.to_string_lossy();
    // Convert backslashes to forward slashes
    let unix_path = path_str.replace('\\', "/");
    // Convert drive letter: C:/... -> /c/...
    if unix_path.len() >= 2 && unix_path.chars().nth(1) == Some(':') {
        let drive = unix_path
            .chars()
            .next()
            .unwrap()
            .to_lowercase()
            .next()
            .unwrap();
        format!("/{}{}", drive, &unix_path[2..])
    } else if unix_path.starts_with("\\\\?\\") || unix_path.starts_with("//?/") {
        // Handle extended path prefix \\?\C:\... or //?/C:/...
        let stripped = unix_path
            .trim_start_matches("\\\\?\\")
            .trim_start_matches("//?/");
        let stripped = stripped.replace('\\', "/");
        if stripped.len() >= 2 && stripped.chars().nth(1) == Some(':') {
            let drive = stripped
                .chars()
                .next()
                .unwrap()
                .to_lowercase()
                .next()
                .unwrap();
            format!("/{}{}", drive, &stripped[2..])
        } else {
            format!("/{}", stripped)
        }
    } else {
        unix_path
    }
}

#[cfg(not(target_os = "windows"))]
fn windows_path_to_docker(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

/// Recursively remap Windows paths in JSON values to Docker-compatible paths
fn remap_json_paths_for_docker(value: &JsonValue) -> JsonValue {
    match value {
        JsonValue::String(s) => {
            // Check if it looks like a Windows path (contains backslash or drive letter)
            if s.contains('\\') || (s.len() >= 2 && s.chars().nth(1) == Some(':')) {
                JsonValue::String(windows_path_to_docker(Path::new(s)))
            } else {
                value.clone()
            }
        }
        JsonValue::Object(map) => {
            let mut new_map = serde_json::Map::new();
            for (k, v) in map {
                new_map.insert(k.clone(), remap_json_paths_for_docker(v));
            }
            JsonValue::Object(new_map)
        }
        JsonValue::Array(arr) => {
            JsonValue::Array(arr.iter().map(remap_json_paths_for_docker).collect())
        }
        _ => value.clone(),
    }
}

/// Normalize a Windows path string (strip extended prefix, convert to backslashes)
#[cfg(target_os = "windows")]
fn normalize_windows_path_str(s: &str) -> String {
    // Strip extended-length path prefix if present
    let stripped = if s.starts_with("\\\\?\\") {
        &s[4..]
    } else if s.starts_with("//?/") {
        &s[4..]
    } else {
        s
    };
    // Convert forward slashes to backslashes
    stripped.replace('/', "\\")
}

/// Extract all file/directory paths from a JSON value (for Docker volume mounting)
#[cfg(target_os = "windows")]
fn extract_paths_from_json(value: &JsonValue, paths: &mut Vec<PathBuf>) {
    match value {
        JsonValue::String(s) => {
            // Check if it looks like a Windows path (with drive letter)
            if looks_like_windows_absolute_path(s) {
                // Normalize: strip extended prefix and convert forward slashes to backslashes
                let normalized = normalize_windows_path_str(s);
                let path = Path::new(&normalized);
                append_desktop_log(&format!(
                    "[JSON Extract] Found path: {} (normalized: {}, exists: {})",
                    s,
                    normalized,
                    path.exists()
                ));
                if path.exists() {
                    // Get the parent directory for files, or the path itself for directories
                    if path.is_file() {
                        if let Some(parent) = path.parent() {
                            append_desktop_log(&format!(
                                "[JSON Extract] Adding parent dir: {}",
                                parent.display()
                            ));
                            paths.push(parent.to_path_buf());
                        }
                        // If it's a CSV file, also extract paths from inside it
                        if s.to_lowercase().ends_with(".csv") {
                            append_desktop_log(&format!(
                                "[JSON Extract] Reading CSV for embedded paths: {}",
                                s
                            ));
                            if let Ok(content) = fs::read_to_string(&path) {
                                extract_paths_from_csv(&content, paths);
                            } else {
                                append_desktop_log(&format!(
                                    "[JSON Extract] ERROR: Failed to read CSV: {}",
                                    s
                                ));
                            }
                        }
                    } else {
                        append_desktop_log(&format!(
                            "[JSON Extract] Adding directory: {}",
                            path.display()
                        ));
                        paths.push(path.to_path_buf());
                    }
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                extract_paths_from_json(v, paths);
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                extract_paths_from_json(v, paths);
            }
        }
        _ => {}
    }
}

/// Check if a string looks like a Windows absolute path
#[cfg(target_os = "windows")]
fn looks_like_windows_absolute_path(s: &str) -> bool {
    // Handle extended-length path prefix: \\?\C:\... or //?/C:/...
    let stripped = if s.starts_with("\\\\?\\") {
        &s[4..]
    } else if s.starts_with("//?/") {
        &s[4..]
    } else {
        s
    };

    if stripped.len() < 3 {
        return false;
    }

    let chars: Vec<char> = stripped.chars().take(3).collect();
    // Check for drive letter pattern: X:\ or X:/
    if chars.len() >= 3 && chars[1] == ':' {
        let first = chars[0];
        let third = chars[2];
        return first.is_ascii_alphabetic() && (third == '\\' || third == '/');
    }
    false
}

/// Extract Windows paths from CSV content
#[cfg(target_os = "windows")]
fn extract_paths_from_csv(content: &str, paths: &mut Vec<PathBuf>) {
    append_desktop_log(&format!(
        "[CSV Extract] Processing CSV content ({} bytes, {} lines)",
        content.len(),
        content.lines().count()
    ));

    let mut found_count = 0;
    let mut added_count = 0;

    for line in content.lines() {
        for field in line.split(',') {
            let field = field.trim().trim_matches('"');
            // Accept both backslash (C:\) and forward slash (C:/) Windows paths
            if looks_like_windows_absolute_path(field) {
                found_count += 1;
                // Normalize: strip extended prefix and convert forward slashes to backslashes
                let normalized = normalize_windows_path_str(field);
                let path = Path::new(&normalized);
                append_desktop_log(&format!(
                    "[CSV Extract] Found path: {} (normalized: {}, exists: {})",
                    field,
                    normalized,
                    path.exists()
                ));
                if path.exists() {
                    if path.is_file() {
                        if let Some(parent) = path.parent() {
                            append_desktop_log(&format!(
                                "[CSV Extract] Adding parent dir: {}",
                                parent.display()
                            ));
                            paths.push(parent.to_path_buf());
                            added_count += 1;
                        }
                    } else {
                        append_desktop_log(&format!(
                            "[CSV Extract] Adding directory: {}",
                            path.display()
                        ));
                        paths.push(path.to_path_buf());
                        added_count += 1;
                    }
                }
            }
        }
    }

    append_desktop_log(&format!(
        "[CSV Extract] Found {} paths, added {} mount candidates",
        found_count, added_count
    ));
}

/// Rewrite a CSV file converting Windows paths to Docker-compatible paths
#[cfg(target_os = "windows")]
fn rewrite_csv_with_docker_paths(csv_path: &Path) -> Result<()> {
    let content = fs::read_to_string(csv_path)
        .with_context(|| format!("Failed to read CSV: {}", csv_path.display()))?;

    append_desktop_log(&format!("[CSV Rewrite] Processing: {}", csv_path.display()));
    let mut converted_count = 0;

    let mut new_lines = Vec::new();
    for line in content.lines() {
        let mut new_fields = Vec::new();
        for field in line.split(',') {
            let trimmed = field.trim();
            let (was_quoted, inner) = if trimmed.starts_with('"') && trimmed.ends_with('"') {
                (true, &trimmed[1..trimmed.len() - 1])
            } else {
                (false, trimmed)
            };

            // Check if it's a Windows path (with backslash or forward slash)
            let new_value = if looks_like_windows_absolute_path(inner) {
                converted_count += 1;
                windows_path_to_docker(Path::new(inner))
            } else {
                inner.to_string()
            };

            if was_quoted {
                new_fields.push(format!("\"{}\"", new_value));
            } else {
                new_fields.push(new_value);
            }
        }
        new_lines.push(new_fields.join(","));
    }

    append_desktop_log(&format!(
        "[CSV Rewrite] Converted {} paths in {}",
        converted_count,
        csv_path.display()
    ));

    let new_content = new_lines.join("\n");
    fs::write(csv_path, new_content)
        .with_context(|| format!("Failed to write converted CSV: {}", csv_path.display()))?;

    Ok(())
}

/// Check if a string looks like a Windows path
#[cfg(target_os = "windows")]
fn is_windows_path(s: &str) -> bool {
    // Regular path: C:\...
    if s.len() >= 2 && s.chars().nth(1) == Some(':') {
        return true;
    }
    // Extended path: \\?\C:\...
    if s.starts_with("\\\\?\\") || s.starts_with("//?/") {
        return true;
    }
    false
}

/// Find and rewrite all CSV files referenced in inputs JSON
#[cfg(target_os = "windows")]
fn rewrite_input_csvs_for_docker(value: &JsonValue) -> Result<()> {
    match value {
        JsonValue::String(s) => {
            if is_windows_path(s) && s.to_lowercase().ends_with(".csv") {
                let path = Path::new(s);
                if path.exists() && path.is_file() {
                    append_desktop_log(&format!("[Pipeline] Rewriting CSV: {}", s));
                    rewrite_csv_with_docker_paths(path)?;
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                rewrite_input_csvs_for_docker(v)?;
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                rewrite_input_csvs_for_docker(v)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Get unique root directories that need to be mounted (deduplicates nested paths)
#[cfg(target_os = "windows")]
fn get_unique_mount_roots(paths: Vec<PathBuf>) -> Vec<PathBuf> {
    use std::collections::HashSet;

    let mut roots: Vec<PathBuf> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();

    for path in paths {
        // Canonicalize to resolve symlinks and normalize
        let canonical = path.canonicalize().unwrap_or(path);

        // Check if this path is already covered by an existing root
        let mut is_covered = false;
        for root in &roots {
            if canonical.starts_with(root) {
                is_covered = true;
                break;
            }
        }

        if !is_covered && !seen.contains(&canonical) {
            // Remove any existing roots that are children of this new path
            roots.retain(|r| !r.starts_with(&canonical));
            roots.push(canonical.clone());
            seen.insert(canonical);
        }
    }

    roots
}

/// Build a PATH that includes Docker's bin directory (needed for credential helpers on Windows)
#[cfg(target_os = "windows")]
fn build_docker_path(docker_bin: &str) -> Option<String> {
    let docker_path = Path::new(docker_bin);
    let docker_dir = docker_path.parent()?;

    let mut paths = vec![docker_dir.to_path_buf()];
    if let Some(existing) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&existing));
    }

    std::env::join_paths(paths)
        .ok()
        .and_then(|joined| joined.into_string().ok())
}

#[cfg(not(target_os = "windows"))]
fn build_docker_path(_docker_bin: &str) -> Option<String> {
    None
}

fn resolve_docker_config_path() -> Option<String> {
    if let Ok(config) = std::env::var("BIOVAULT_DOCKER_CONFIG") {
        let trimmed = config.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    if let Ok(config) = std::env::var("DOCKER_CONFIG") {
        let trimmed = config.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    None
}

fn apply_docker_config_arg(cmd: &mut Command) {
    if let Some(config) = resolve_docker_config_path() {
        append_desktop_log(&format!("[Pipeline] Using Docker config: {}", config));
        cmd.arg("--config").arg(config);
    }
}

/// Pull a Docker image if not already present (needed on Windows for credential helper PATH issues)
#[cfg(target_os = "windows")]
fn pull_docker_image_if_needed(docker_bin: &str, image: &str) -> Result<()> {
    append_desktop_log(&format!(
        "[Pipeline] Checking if image {} is available...",
        image
    ));

    // Check if image exists locally
    let mut check_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut check_cmd);
    apply_docker_config_arg(&mut check_cmd);
    check_cmd
        .arg("image")
        .arg("inspect")
        .arg(image)
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        check_cmd.env("PATH", &docker_path);
    }

    if check_cmd.status().map(|s| s.success()).unwrap_or(false) {
        append_desktop_log(&format!(
            "[Pipeline] Image {} is already available locally",
            image
        ));
        return Ok(());
    }

    // Pull the image
    append_desktop_log(&format!("[Pipeline] Pulling image {}...", image));
    println!("📦 Pulling Docker image: {}", image);

    let mut pull_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut pull_cmd);
    apply_docker_config_arg(&mut pull_cmd);
    pull_cmd.arg("pull").arg(image);

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        pull_cmd.env("PATH", &docker_path);
    }

    let status = pull_cmd.status().context("Failed to execute docker pull")?;

    if !status.success() {
        append_desktop_log(&format!("[Pipeline] Failed to pull image {}", image));
        return Err(anyhow::anyhow!("Failed to pull Docker image: {}", image).into());
    }

    append_desktop_log(&format!("[Pipeline] Successfully pulled image {}", image));
    Ok(())
}

#[cfg(target_os = "windows")]
fn ensure_nextflow_runner_image(docker_bin: &str) -> Result<&'static str> {
    // The stock `nextflow/nextflow` image includes an outdated Docker CLI (17.x), which cannot
    // talk to modern Docker daemons. When the workflow uses `container` directives, each task
    // fails with exit code 125 and messages like:
    //   "client version 1.32 is too old. Minimum supported API version is 1.44"
    //
    // We use a pre-built image from GitHub Container Registry that includes a modern docker client.
    const NEXTFLOW_RUNNER_IMAGE: &str = "ghcr.io/openmined/nextflow-runner:25.10.2";

    // Check if runner image exists locally
    let mut check_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut check_cmd);
    apply_docker_config_arg(&mut check_cmd);
    check_cmd
        .arg("image")
        .arg("inspect")
        .arg(NEXTFLOW_RUNNER_IMAGE)
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    if let Some(docker_path) = build_docker_path(docker_bin) {
        check_cmd.env("PATH", &docker_path);
    }

    if check_cmd.status().map(|s| s.success()).unwrap_or(false) {
        append_desktop_log(&format!(
            "[Pipeline] Using Nextflow runner image: {}",
            NEXTFLOW_RUNNER_IMAGE
        ));
        return Ok(NEXTFLOW_RUNNER_IMAGE);
    }

    // Pull the pre-built image from GitHub Container Registry
    append_desktop_log(&format!(
        "[Pipeline] Pulling Nextflow runner image from {}",
        NEXTFLOW_RUNNER_IMAGE
    ));

    pull_docker_image_if_needed(docker_bin, NEXTFLOW_RUNNER_IMAGE)?;

    append_desktop_log(&format!(
        "[Pipeline] Nextflow runner image {} ready",
        NEXTFLOW_RUNNER_IMAGE
    ));
    Ok(NEXTFLOW_RUNNER_IMAGE)
}

#[cfg(not(target_os = "windows"))]
fn ensure_nextflow_runner_image(_docker_bin: &str) -> Result<&'static str> {
    Ok("nextflow/nextflow:25.10.2")
}

#[cfg(target_os = "linux")]
fn check_docker_running(docker_bin: &str) -> Result<()> {
    let mut cmd = Command::new(docker_bin);
    apply_docker_config_arg(&mut cmd);
    cmd.arg("info")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        cmd.env("PATH", &docker_path);
    }

    let output = cmd
        .output()
        .with_context(|| format!("Failed to execute '{}'", docker_bin))?;

    if output.status.success() {
        return Ok(());
    }

    let stderr = String::from_utf8_lossy(&output.stderr);
    let mut base_msg = format!(
        "Docker daemon is not running ({} exited with {:?}). Please start Docker Desktop or the Docker service and retry.",
        docker_bin,
        output.status.code()
    );

    // Add a more helpful hint when the socket is reachable but permission is denied.
    let lowered = stderr.to_ascii_lowercase();
    if lowered.contains("permission denied") || lowered.contains("connect: permission denied") {
        if std::env::var_os("APPIMAGE").is_some() {
            base_msg.push_str(
                " It looks like this session cannot access /var/run/docker.sock. \
If you recently added your user to the docker group, log out and back in (do not rely on 'newgrp docker' inside the AppImage), then relaunch the app.",
            );
        } else {
            base_msg.push_str(
                " It looks like this session cannot access /var/run/docker.sock. \
If you recently added your user to the docker group, log out and back in, then retry.",
            );
        }
    }

    if !stderr.trim().is_empty() {
        base_msg.push_str(&format!(" Stderr: {}", stderr.trim()));
    }

    Err(anyhow::anyhow!(base_msg).into())
}

#[cfg(not(target_os = "linux"))]
fn check_docker_running(docker_bin: &str) -> Result<()> {
    let mut cmd = Command::new(docker_bin);
    super::configure_child_process(&mut cmd);
    apply_docker_config_arg(&mut cmd);
    cmd.arg("info").stdout(Stdio::null()).stderr(Stdio::null());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        cmd.env("PATH", &docker_path);
    }

    let status = cmd
        .status()
        .with_context(|| format!("Failed to execute '{}'", docker_bin))?;

    if status.success() {
        return Ok(());
    }

    Err(anyhow::anyhow!(
        "Docker daemon is not running ({} exited with {:?}). Please start Docker Desktop or the Docker service and retry.",
        docker_bin,
        status.code()
    )
    .into())
}

pub async fn execute_dynamic(
    project_folder: &str,
    args: Vec<String>,
    dry_run: bool,
    resume: bool,
    results_dir: Option<String>,
    run_settings: RunSettings,
) -> Result<()> {
    let project_path = Path::new(project_folder);
    if !project_path.exists() {
        return Err(anyhow::anyhow!("Project folder does not exist: {}", project_folder).into());
    }
    let project_abs = project_path
        .canonicalize()
        .context("Failed to resolve project path")?;

    let nextflow_log_path = project_path.join(".nextflow.log");
    fs::remove_file(&nextflow_log_path).ok();

    let spec_path = project_path.join("project.yaml");
    if !spec_path.exists() {
        return Err(anyhow::anyhow!(
            "project.yaml not found in {}. Use 'bv project create' first.",
            project_folder
        )
        .into());
    }

    let spec = ProjectSpec::load(&spec_path)?;

    if spec.template.as_deref() != Some("dynamic-nextflow") {
        return Err(anyhow::anyhow!(
            "This project uses template '{}'. Only 'dynamic-nextflow' is supported by the new run system.",
            spec.template.as_deref().unwrap_or("(none)")
        ).into());
    }

    println!("🚀 Running project: {}", spec.name.bold());

    let parsed_args = parse_cli_args(&args)?;
    let nextflow_args = parsed_args.passthrough.clone();

    validate_no_clashes(&spec, &parsed_args)?;

    let inputs_json = build_inputs_json(&spec, &parsed_args, project_path)?;
    let mut params_json = build_params_json(&spec, &parsed_args)?;

    let assets_dir_path = project_path.join("assets");
    let assets_dir_abs = assets_dir_path
        .canonicalize()
        .unwrap_or_else(|_| assets_dir_path.clone());

    params_json
        .entry("assets_dir".to_string())
        .or_insert_with(|| json!(assets_dir_abs.to_string_lossy().to_string()));

    let results_path = results_dir.as_deref().unwrap_or("results");
    let results_path_buf = if Path::new(results_path).is_absolute() {
        PathBuf::from(results_path)
    } else {
        project_path.join(results_path)
    };
    if !dry_run {
        fs::create_dir_all(&results_path_buf).context("Failed to create results directory")?;
    }
    let results_path_str = results_path_buf.to_string_lossy().to_string();

    // Check user workflow exists
    let workflow_path = project_path.join(&spec.workflow);
    if !workflow_path.exists() {
        return Err(anyhow::anyhow!(
            "Workflow file not found: {}. Expected at: {}",
            spec.workflow,
            workflow_path.display()
        )
        .into());
    }

    // Load template from .biovault/env/{template_name}/ (security boundary)
    let biovault_home = crate::config::get_biovault_home()?;
    let biovault_home_abs = biovault_home
        .canonicalize()
        .unwrap_or_else(|_| biovault_home.clone());
    let template_name = spec.template.as_deref().unwrap_or("dynamic-nextflow");
    let env_dir = biovault_home_abs.join("env").join(template_name);
    let template_path = env_dir.join("template.nf");

    if template_name == "dynamic-nextflow" {
        install_dynamic_template(&biovault_home_abs)?;
    }

    if !template_path.exists() {
        return Err(anyhow::anyhow!(
            "Template not found: {}. Run 'bv init' to install templates.",
            template_path.display()
        )
        .into());
    }

    // Canonicalize paths for Nextflow
    let template_abs = template_path
        .canonicalize()
        .context("Failed to resolve template path")?;

    let workflow_abs = workflow_path
        .canonicalize()
        .context("Failed to resolve workflow path")?;

    let project_spec_abs = spec_path
        .canonicalize()
        .context("Failed to resolve project spec path")?;

    let inputs_json_str =
        serde_json::to_string(&inputs_json).context("Failed to encode inputs metadata to JSON")?;
    let params_json_str = serde_json::to_string(&params_json)
        .context("Failed to encode parameters metadata to JSON")?;

    let config = crate::config::get_config().ok();
    let nextflow_bin =
        resolve_binary_path(config.as_ref(), "nextflow").unwrap_or_else(|| "nextflow".to_string());
    let docker_bin =
        resolve_binary_path(config.as_ref(), "docker").unwrap_or_else(|| "docker".to_string());

    // Check Docker availability (required for both Windows Docker execution and workflow containers)
    // Skip Docker checks in dry-run mode
    if !dry_run && (run_settings.require_docker || should_use_docker_for_nextflow()) {
        append_desktop_log("[Pipeline] Checking Docker availability...");
        if let Err(err) = check_docker_running(&docker_bin) {
            append_desktop_log(&format!("[Pipeline] Docker check failed: {}", err));
            return Err(err);
        }
        append_desktop_log("[Pipeline] Docker is running (docker info succeeded)");
    }

    // Build command - use Docker on Windows, native Nextflow elsewhere
    let mut cmd = if should_use_docker_for_nextflow() {
        append_desktop_log("[Pipeline] Using Docker to run Nextflow (Windows mode)");

        // Build/get Nextflow runner image with modern Docker CLI
        // Skip image pulling in dry-run mode - use placeholder
        let nextflow_image = if dry_run {
            "ghcr.io/openmined/nextflow-runner:latest"
        } else {
            ensure_nextflow_runner_image(&docker_bin)?
        };

        // Convert all paths to Docker-compatible format
        let docker_biovault_home = windows_path_to_docker(&biovault_home_abs);
        let docker_project_path = windows_path_to_docker(&project_abs);
        let docker_template = windows_path_to_docker(&template_abs);
        let docker_workflow = windows_path_to_docker(&workflow_abs);
        let docker_project_spec = windows_path_to_docker(&project_spec_abs);
        let docker_log_path = windows_path_to_docker(&nextflow_log_path);
        let results_abs = results_path_buf
            .canonicalize()
            .unwrap_or_else(|_| results_path_buf.clone());
        let docker_results = windows_path_to_docker(&results_abs);

        // Extract paths from inputs that need to be mounted (must do before rewriting CSVs)
        // This is Windows-specific: extract paths from CSV files and rewrite them for Docker
        let inputs_json_value: JsonValue =
            serde_json::to_value(&inputs_json).context("Failed to convert inputs to JSON value")?;

        #[cfg(target_os = "windows")]
        let mount_roots = {
            let mut data_paths: Vec<PathBuf> = Vec::new();
            extract_paths_from_json(&inputs_json_value, &mut data_paths);
            let roots = get_unique_mount_roots(data_paths);

            // Rewrite CSV files to convert Windows paths to Docker-compatible paths
            // Skip in dry-run mode since this modifies files
            if !dry_run {
                append_desktop_log(
                    "[Pipeline] Rewriting CSV files with Docker-compatible paths...",
                );
                append_desktop_log(&format!(
                    "[Pipeline] inputs_json: {}",
                    serde_json::to_string(&inputs_json_value).unwrap_or_default()
                ));
                rewrite_input_csvs_for_docker(&inputs_json_value)?;
            }

            roots
        };

        #[cfg(not(target_os = "windows"))]
        let mount_roots: Vec<PathBuf> = Vec::new();

        append_desktop_log("[Pipeline] Docker path mappings:");
        append_desktop_log(&format!(
            "  biovault_home: {} -> {}",
            biovault_home_abs.display(),
            docker_biovault_home
        ));
        append_desktop_log(&format!(
            "  project_path: {} -> {}",
            project_abs.display(),
            docker_project_path
        ));
        append_desktop_log(&format!(
            "  template: {} -> {}",
            template_abs.display(),
            docker_template
        ));
        append_desktop_log(&format!(
            "[Pipeline] Additional data mounts: {:?}",
            mount_roots
        ));

        let mut docker_cmd = Command::new(&docker_bin);
        super::configure_child_process(&mut docker_cmd);
        apply_docker_config_arg(&mut docker_cmd);

        // Add Docker PATH for credential helpers on Windows
        if let Some(docker_path) = build_docker_path(&docker_bin) {
            docker_cmd.env("PATH", docker_path);
        }

        docker_cmd
            .arg("run")
            .arg("--rm")
            // Mount Docker Desktop engine socket so Nextflow-in-Docker can launch workflow containers.
            //
            // With Docker Desktop (Linux containers), the daemon runs in a Linux VM and provides a Unix socket.
            // Binding it into the Nextflow container allows Nextflow to launch workflow containers.
            .arg("-v")
            .arg("/var/run/docker.sock:/var/run/docker.sock")
            // Mount the BioVault home directory
            .arg("-v")
            .arg(format!(
                "{}:{}",
                windows_path_to_docker(&biovault_home),
                docker_biovault_home
            ))
            // Mount the project path (may be same as above, Docker handles duplicates)
            .arg("-v")
            .arg(format!(
                "{}:{}",
                windows_path_to_docker(project_path),
                docker_project_path
            ));

        // Mount additional data directories discovered from inputs
        for mount_path in &mount_roots {
            let docker_mount = windows_path_to_docker(mount_path);
            append_desktop_log(&format!(
                "[Pipeline] Adding mount: {} -> {}",
                mount_path.display(),
                docker_mount
            ));
            docker_cmd
                .arg("-v")
                .arg(format!("{}:{}", docker_mount, docker_mount));
        }

        docker_cmd
            // Set working directory
            .arg("-w")
            .arg(&docker_project_path)
            // Use Nextflow runner image with modern Docker CLI
            .arg(nextflow_image)
            // Nextflow command (container entrypoint is bash, not nextflow)
            .arg("nextflow")
            // Nextflow arguments
            .arg("-log")
            .arg(&docker_log_path)
            .arg("run")
            .arg(&docker_template);

        if resume {
            docker_cmd.arg("-resume");
        }

        for extra in &nextflow_args {
            docker_cmd.arg(extra);
        }

        // Re-encode JSON with Docker paths (inputs_json_value already created above)
        let params_json_value: JsonValue =
            serde_json::to_value(&params_json).context("Failed to convert params to JSON value")?;
        let docker_inputs_json = remap_json_paths_for_docker(&inputs_json_value);
        let docker_params_json = remap_json_paths_for_docker(&params_json_value);
        let docker_inputs_json_str = serde_json::to_string(&docker_inputs_json)
            .context("Failed to encode Docker inputs metadata to JSON")?;
        let docker_params_json_str = serde_json::to_string(&docker_params_json)
            .context("Failed to encode Docker parameters metadata to JSON")?;

        docker_cmd
            .arg("--work_flow_file")
            .arg(&docker_workflow)
            .arg("--project_spec")
            .arg(&docker_project_spec)
            .arg("--inputs_json")
            .arg(docker_inputs_json_str)
            .arg("--params_json")
            .arg(docker_params_json_str)
            .arg("--results_dir")
            .arg(&docker_results);

        docker_cmd
    } else {
        // Native Nextflow execution (macOS, Linux)
        append_desktop_log(&format!(
            "[Pipeline] Using native nextflow binary: {}",
            nextflow_bin
        ));

        // Log original PATH from environment
        if let Some(original_path) = std::env::var_os("PATH") {
            append_desktop_log(&format!(
                "[Pipeline] Original PATH from environment: {}",
                original_path.to_string_lossy()
            ));
        } else {
            append_desktop_log("[Pipeline] WARNING: No PATH environment variable found!");
        }

        let mut native_cmd = Command::new(&nextflow_bin);
        super::configure_child_process(&mut native_cmd);

        append_desktop_log("[Pipeline] Preferred binary paths:");
        for binary in ["nextflow", "java", "docker"] {
            if let Some(path) = resolve_binary_path(config.as_ref(), binary) {
                append_desktop_log(&format!("  {} = {}", binary, path));
            } else {
                append_desktop_log(&format!("  {} = <not configured>", binary));
            }
        }

        if let Some(path_env) = build_augmented_path(config.as_ref()) {
            append_desktop_log(&format!(
                "[Pipeline] Final augmented PATH for nextflow: {}",
                path_env
            ));
            native_cmd.env("PATH", path_env);
        } else {
            append_desktop_log(
                "[Pipeline] WARNING: Could not build augmented PATH, using system PATH",
            );
        }

        native_cmd.arg("-log").arg(&nextflow_log_path);

        native_cmd.arg("run").arg(&template_abs);

        if resume {
            native_cmd.arg("-resume");
        }

        for extra in &nextflow_args {
            native_cmd.arg(extra);
        }

        native_cmd
            .arg("--work_flow_file")
            .arg(&workflow_abs)
            .arg("--project_spec")
            .arg(&project_spec_abs)
            .arg("--inputs_json")
            .arg(inputs_json_str)
            .arg("--params_json")
            .arg(params_json_str)
            .arg("--results_dir")
            .arg(&results_path_str);

        native_cmd
    };

    let display_cmd = format_command(&cmd);

    if dry_run {
        println!("\n🔍 Dry run - would execute:");
        println!("  {}", display_cmd.dimmed());
        append_desktop_log(&format!(
            "[Pipeline] (dry-run) Nextflow command: {}",
            display_cmd
        ));
        return Ok(());
    }

    println!("\n▶️  Executing Nextflow...\n");
    println!("  {}", display_cmd.dimmed());
    append_desktop_log(&format!("[Pipeline] Nextflow command: {}", display_cmd));

    cmd.current_dir(project_path);
    let work_dir = project_path.join("work");
    let status = execute_with_logging(
        cmd,
        Some(nextflow_log_path),
        Some(work_dir),
        Some(project_path.to_path_buf()),
    )
    .context("Failed to execute nextflow")?;

    if !status.success() {
        append_desktop_log(&format!(
            "[Pipeline] Nextflow exited with status: {:?}",
            status.code()
        ));
        return Err(
            anyhow::anyhow!("Nextflow execution failed with code: {:?}", status.code()).into(),
        );
    }

    println!("\n✅ Workflow completed successfully!");
    append_desktop_log("[Pipeline] Workflow completed successfully!");
    Ok(())
}

#[derive(Debug)]
struct ParsedArgs {
    inputs: HashMap<String, InputArg>,
    params: HashMap<String, String>,
    passthrough: Vec<String>,
}

#[derive(Debug)]
struct InputArg {
    value: String,
    format_override: Option<String>,
}

fn parse_cli_args(args: &[String]) -> Result<ParsedArgs> {
    let mut inputs = HashMap::new();
    let mut params = HashMap::new();
    let mut format_overrides = HashMap::new();
    let mut passthrough = Vec::new();

    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];

        if arg == "--" {
            passthrough.extend(args[i + 1..].iter().cloned());
            break;
        }

        if !arg.starts_with("--") {
            passthrough.push(arg.clone());
            i += 1;
            continue;
        }

        let key = arg.strip_prefix("--").unwrap();

        if key == "set" {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];

            let (target, val) = value.split_once('=').ok_or_else(|| {
                anyhow::anyhow!(
                    "Invalid --set assignment '{}'. Use inputs.name=value or params.name=value.",
                    value
                )
            })?;

            if let Some(input_name) = target.strip_prefix("inputs.") {
                inputs.insert(
                    input_name.to_string(),
                    InputArg {
                        value: val.to_string(),
                        format_override: None,
                    },
                );
            } else if let Some(param_name) = target.strip_prefix("params.") {
                params.insert(param_name.to_string(), val.to_string());
            } else if let Some(param_name) = target.strip_prefix("param.") {
                params.insert(param_name.to_string(), val.to_string());
            } else {
                return Err(anyhow::anyhow!(
                    "Unsupported --set target '{}'. Expected inputs.<name> or params.<name>.",
                    target
                )
                .into());
            }
            i += 2;
            continue;
        }

        if key.starts_with("param.") {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];
            let param_name = key.strip_prefix("param.").unwrap();
            params.insert(param_name.to_string(), value.clone());
            i += 2;
            continue;
        }

        if key.contains(".format") {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];
            let input_name = key.strip_suffix(".format").unwrap();
            format_overrides.insert(input_name.to_string(), value.clone());
            i += 2;
            continue;
        }

        if key.contains(".mapping.") {
            // Future: support inline mapping overrides
            return Err(
                anyhow::anyhow!("Inline mapping overrides not yet supported: {}", key).into(),
            );
        }

        match key {
            "results-dir" | "results_dir" => {
                i += 2;
                continue;
            }
            _ => {}
        }

        if i + 1 >= args.len() {
            return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
        }

        let value = &args[i + 1];
        inputs.insert(
            key.to_string(),
            InputArg {
                value: value.clone(),
                format_override: None,
            },
        );

        i += 2;
    }

    for (input_name, format) in &format_overrides {
        if let Some(input) = inputs.get_mut(input_name) {
            input.format_override = Some(format.clone());
        }
    }

    Ok(ParsedArgs {
        inputs,
        params,
        passthrough,
    })
}

fn validate_no_clashes(spec: &ProjectSpec, parsed: &ParsedArgs) -> Result<()> {
    let input_names: Vec<&str> = spec.inputs.iter().map(|i| i.name.as_str()).collect();
    let output_names: Vec<&str> = spec.outputs.iter().map(|o| o.name.as_str()).collect();

    for param_name in parsed.params.keys() {
        if input_names.contains(&param_name.as_str()) {
            return Err(anyhow::anyhow!(
                "Parameter '{}' clashes with input name. Use --param.{} instead.",
                param_name,
                param_name
            )
            .into());
        }
        if output_names.contains(&param_name.as_str()) {
            return Err(anyhow::anyhow!(
                "Parameter '{}' clashes with output name. Use --param.{} instead.",
                param_name,
                param_name
            )
            .into());
        }
    }

    for input_name in parsed.inputs.keys() {
        if !input_names.contains(&input_name.as_str()) {
            println!(
                "⚠️  Warning: Unknown input '{}'. Expected inputs: {}",
                input_name.yellow(),
                input_names.join(", ")
            );
        }
    }

    Ok(())
}

fn build_inputs_json(
    spec: &ProjectSpec,
    parsed: &ParsedArgs,
    _project_path: &Path,
) -> Result<HashMap<String, JsonValue>> {
    let mut inputs_json = HashMap::new();

    for input_spec in &spec.inputs {
        if let Some(input_arg) = parsed.inputs.get(&input_spec.name) {
            let path_str = &input_arg.value;
            let path = Path::new(path_str);

            if !path.exists() {
                return Err(anyhow::anyhow!("Input file not found: {}", path_str).into());
            }

            let format = input_arg
                .format_override
                .as_deref()
                .or(input_spec.format.as_deref())
                .or_else(|| detect_format(path))
                .unwrap_or("unknown");

            inputs_json.insert(
                input_spec.name.clone(),
                json!({
                    "path": path.canonicalize()?.to_string_lossy().to_string(),
                    "type": input_spec.raw_type,
                    "format": format,
                    "mapping": input_spec.mapping,
                }),
            );
        } else if !input_spec.raw_type.ends_with('?') {
            return Err(
                anyhow::anyhow!("Required input '{}' not provided", input_spec.name).into(),
            );
        }
    }

    Ok(inputs_json)
}

fn build_params_json(
    spec: &ProjectSpec,
    parsed: &ParsedArgs,
) -> Result<HashMap<String, JsonValue>> {
    let mut params_json = HashMap::new();

    for param_spec in &spec.parameters {
        let value = if let Some(v) = parsed.params.get(&param_spec.name) {
            match param_spec.raw_type.as_str() {
                "Bool" => {
                    let bool_val = v.parse::<bool>().context(format!(
                        "Parameter '{}' expects Bool, got '{}'",
                        param_spec.name, v
                    ))?;
                    json!(bool_val)
                }
                "String" => json!(v),
                ty if ty.starts_with("Enum") => json!(v),
                unsupported => {
                    return Err(
                        anyhow::anyhow!("Unsupported parameter type: {}", unsupported).into(),
                    );
                }
            }
        } else if let Some(default) = &param_spec.default {
            serde_json::to_value(default)
                .context("Failed to convert default param value to JSON")?
        } else {
            continue;
        };

        params_json.insert(param_spec.name.clone(), value);
    }

    Ok(params_json)
}

fn detect_format(path: &Path) -> Option<&'static str> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .and_then(|ext| match ext.to_lowercase().as_str() {
            "json" => Some("json"),
            "yaml" | "yml" => Some("yaml"),
            "csv" => Some("csv"),
            "tsv" => Some("tsv"),
            "vcf" | "vcf.gz" => Some("vcf"),
            _ => None,
        })
}

fn install_dynamic_template(biovault_home: &Path) -> Result<()> {
    let env_dir = biovault_home.join("env").join("dynamic-nextflow");
    if !env_dir.exists() {
        fs::create_dir_all(&env_dir).context("Failed to create dynamic template directory")?;
    }

    let template_path = env_dir.join("template.nf");
    let template_contents = include_str!("../../templates/dynamic/template.nf");
    fs::write(&template_path, template_contents)
        .context("Failed to install dynamic template.nf")?;

    println!("📦 Dynamic template ready at {}", template_path.display());

    let config_path = env_dir.join("nextflow.config");
    let config_contents = r#"process.executor = 'local'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
"#;
    fs::write(&config_path, config_contents)
        .context("Failed to install dynamic nextflow.config")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project_spec::{InputSpec, ProjectSpec};
    use tempfile::TempDir;

    fn sample_project_spec() -> ProjectSpec {
        ProjectSpec {
            name: "test".to_string(),
            author: "author".to_string(),
            workflow: "workflow.nf".to_string(),
            template: Some("dynamic-nextflow".to_string()),
            version: None,
            assets: vec![],
            parameters: vec![],
            inputs: vec![
                InputSpec {
                    name: "samplesheet".to_string(),
                    raw_type: "File".to_string(),
                    description: None,
                    format: Some("csv".to_string()),
                    path: None,
                    mapping: None,
                },
                InputSpec {
                    name: "data_dir".to_string(),
                    raw_type: "Directory".to_string(),
                    description: None,
                    format: None,
                    path: None,
                    mapping: None,
                },
            ],
            outputs: vec![],
        }
    }

    #[test]
    fn build_inputs_json_handles_file_and_directory() {
        let tmp = TempDir::new().unwrap();
        let file_path = tmp.path().join("participants.csv");
        std::fs::write(&file_path, "id,path\n1,a.txt\n").unwrap();
        let dir_path = tmp.path().join("data");
        std::fs::create_dir_all(&dir_path).unwrap();

        let parsed = ParsedArgs {
            inputs: HashMap::from([
                (
                    "samplesheet".to_string(),
                    InputArg {
                        value: file_path.to_string_lossy().to_string(),
                        format_override: None,
                    },
                ),
                (
                    "data_dir".to_string(),
                    InputArg {
                        value: dir_path.to_string_lossy().to_string(),
                        format_override: None,
                    },
                ),
            ]),
            params: HashMap::new(),
            passthrough: Vec::new(),
        };

        let project_spec = sample_project_spec();
        let inputs = build_inputs_json(&project_spec, &parsed, tmp.path()).unwrap();

        let sheet_entry = inputs.get("samplesheet").expect("samplesheet entry");
        assert_eq!(sheet_entry["type"], json!("File"));
        assert_eq!(sheet_entry["format"], json!("csv"));
        let sheet_path = sheet_entry["path"].as_str().unwrap();
        assert_eq!(
            sheet_path,
            file_path.canonicalize().unwrap().to_string_lossy()
        );

        let dir_entry = inputs.get("data_dir").expect("data_dir entry");
        assert_eq!(dir_entry["type"], json!("Directory"));
        let dir_json_path = dir_entry["path"].as_str().unwrap();
        assert_eq!(
            dir_json_path,
            dir_path.canonicalize().unwrap().to_string_lossy()
        );
    }

    #[test]
    fn parse_cli_args_supports_set_inputs_and_params() {
        let args = vec![
            "--set".to_string(),
            "inputs.samplesheet=/tmp/sheet.csv".to_string(),
            "--set".to_string(),
            "params.threshold=0.5".to_string(),
        ];

        let parsed = parse_cli_args(&args).expect("parse --set inputs");

        let sheet = parsed
            .inputs
            .get("samplesheet")
            .expect("samplesheet input parsed");
        assert_eq!(sheet.value, "/tmp/sheet.csv");

        let threshold = parsed
            .params
            .get("threshold")
            .expect("param threshold parsed");
        assert_eq!(threshold, "0.5");
        assert!(parsed.passthrough.is_empty());
    }

    #[test]
    fn parse_cli_args_ignores_results_dir() {
        let args = vec![
            "--results-dir".to_string(),
            "custom_results".to_string(),
            "--samplesheet".to_string(),
            "/tmp/sheet.csv".to_string(),
        ];

        let parsed = parse_cli_args(&args).unwrap();
        assert!(parsed.inputs.contains_key("samplesheet"));
        assert!(!parsed.inputs.contains_key("results-dir"));
        assert!(parsed.passthrough.is_empty());
    }

    #[test]
    fn parse_cli_args_captures_nextflow_flags() {
        let args = vec![
            "--samplesheet".to_string(),
            "/tmp/sheet.csv".to_string(),
            "-with-singularity".to_string(),
            "-profile".to_string(),
            "docker".to_string(),
        ];

        let parsed = parse_cli_args(&args).expect("parse passthrough flags");
        assert_eq!(
            parsed.passthrough,
            vec![
                "-with-singularity".to_string(),
                "-profile".to_string(),
                "docker".to_string()
            ]
        );
    }
}
