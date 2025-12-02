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
use std::process::Command;

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

pub async fn execute_dynamic(
    project_folder: &str,
    args: Vec<String>,
    dry_run: bool,
    resume: bool,
    results_dir: Option<String>,
) -> Result<()> {
    let project_path = Path::new(project_folder);
    if !project_path.exists() {
        return Err(anyhow::anyhow!("Project folder does not exist: {}", project_folder).into());
    }

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
    let template_name = spec.template.as_deref().unwrap_or("dynamic-nextflow");
    let env_dir = biovault_home.join("env").join(template_name);
    let template_path = env_dir.join("template.nf");

    if template_name == "dynamic-nextflow" {
        install_dynamic_template(&biovault_home)?;
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

    // Log environment details for debugging
    append_desktop_log(&format!(
        "[Pipeline] Using nextflow binary: {}",
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

    let mut cmd = Command::new(&nextflow_bin);

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
        cmd.env("PATH", path_env);
    } else {
        append_desktop_log("[Pipeline] WARNING: Could not build augmented PATH, using system PATH");
    }

    cmd.arg("-log").arg(&nextflow_log_path);

    cmd.arg("run").arg(&template_abs);

    if resume {
        cmd.arg("-resume");
    }

    for extra in &nextflow_args {
        cmd.arg(extra);
    }

    cmd.arg("--work_flow_file")
        .arg(&workflow_abs)
        .arg("--project_spec")
        .arg(&project_spec_abs)
        .arg("--inputs_json")
        .arg(inputs_json_str)
        .arg("--params_json")
        .arg(params_json_str)
        .arg("--results_dir")
        .arg(results_path);

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
    let status =
        execute_with_logging(cmd, Some(nextflow_log_path)).context("Failed to execute nextflow")?;

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
