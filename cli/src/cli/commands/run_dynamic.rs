use crate::error::Result;
use crate::project_spec::ProjectSpec;
use anyhow::Context;
use colored::Colorize;
use serde_json::{json, Value as JsonValue};
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Command;

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

    validate_no_clashes(&spec, &parsed_args)?;

    let inputs_json = build_inputs_json(&spec, &parsed_args, project_path)?;
    let params_json = build_params_json(&spec, &parsed_args)?;

    let runtime_spec = json!({
        "inputs": inputs_json,
        "parameters": params_json,
    });

    let runtime_json_path = project_path.join(".bv_runtime.json");
    fs::write(
        &runtime_json_path,
        serde_json::to_string_pretty(&runtime_spec)?,
    )
    .context("Failed to write runtime spec JSON")?;

    // Canonicalize runtime JSON path for Nextflow
    let runtime_json_abs = runtime_json_path
        .canonicalize()
        .context("Failed to resolve runtime JSON path")?;

    println!("📋 Runtime spec written to: {}", runtime_json_abs.display());

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

    if !template_path.exists() {
        return Err(anyhow::anyhow!(
            "Template not found: {}. Run 'bv init' to install templates.",
            template_path.display()
        )
        .into());
    }

    // Canonicalize template path for Nextflow
    let template_abs = template_path
        .canonicalize()
        .context("Failed to resolve template path")?;

    // Also canonicalize workflow path for --work_flow_file parameter
    let workflow_abs = workflow_path
        .canonicalize()
        .context("Failed to resolve workflow path")?;

    let mut cmd = Command::new("nextflow");
    cmd.arg("run")
        .arg(&template_abs)
        .arg("--dynamic_spec_json")
        .arg(&runtime_json_abs)
        .arg("--work_flow_file")
        .arg(&workflow_abs)
        .arg("--results_dir")
        .arg(results_path);

    if resume {
        cmd.arg("-resume");
    }

    if dry_run {
        println!("\n🔍 Dry run - would execute:");
        println!("  {}", format!("{:?}", cmd).dimmed());
        return Ok(());
    }

    println!("\n▶️  Executing Nextflow...\n");

    let status = cmd
        .current_dir(project_path)
        .status()
        .context("Failed to execute nextflow")?;

    if !status.success() {
        return Err(
            anyhow::anyhow!("Nextflow execution failed with code: {:?}", status.code()).into(),
        );
    }

    println!("\n✅ Workflow completed successfully!");
    Ok(())
}

#[derive(Debug)]
struct ParsedArgs {
    inputs: HashMap<String, InputArg>,
    params: HashMap<String, String>,
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

    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];

        if !arg.starts_with("--") {
            i += 1;
            continue;
        }

        let key = arg.strip_prefix("--").unwrap();

        if i + 1 >= args.len() {
            return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
        }

        let value = &args[i + 1];

        if key.starts_with("param.") {
            let param_name = key.strip_prefix("param.").unwrap();
            params.insert(param_name.to_string(), value.clone());
        } else if key.contains(".format") {
            let input_name = key.strip_suffix(".format").unwrap();
            format_overrides.insert(input_name.to_string(), value.clone());
        } else if key.contains(".mapping.") {
            // Future: support inline mapping overrides
            return Err(
                anyhow::anyhow!("Inline mapping overrides not yet supported: {}", key).into(),
            );
        } else {
            inputs.insert(
                key.to_string(),
                InputArg {
                    value: value.clone(),
                    format_override: None,
                },
            );
        }

        i += 2;
    }

    for (input_name, format) in &format_overrides {
        if let Some(input) = inputs.get_mut(input_name) {
            input.format_override = Some(format.clone());
        }
    }

    Ok(ParsedArgs { inputs, params })
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
