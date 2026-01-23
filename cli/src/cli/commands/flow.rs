fn append_desktop_log(message: &str) {
    if let Ok(path) = std::env::var("BIOVAULT_DESKTOP_LOG_FILE") {
        if path.is_empty() {
            return;
        }
        if let Some(parent) = std::path::Path::new(&path).parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        let timestamp = chrono::Local::now().format("%Y-%m-%dT%H:%M:%S%:z");
        let line = format!(
            "[{}][INFO] {}
",
            timestamp, message
        );
        let _ = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .and_then(|mut file| std::io::Write::write_all(&mut file, line.as_bytes()));
    }
}

use crate::cli::syft_url::SyftURL;
use crate::data::BioVaultDb;
use crate::error::Result;
use crate::flow_spec::{
    value_to_string, ConditionalInput, FlowInputSpec, FlowPermissionAccessSpec,
    FlowPermissionRuleSpec, FlowPermissionSpec, FlowShareSpec, FlowSpec, FlowSqlStoreSpec,
    FlowStepSpec, FlowStoreSpec,
};
use crate::module_spec::{InputSpec, ModuleSpec, OutputSpec};
use crate::types::{AccessControl, PermissionRule, SyftPermissions};
use anyhow::{anyhow, Context};
use chrono::Utc;
use colored::Colorize;
use csv::ReaderBuilder;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use rusqlite::params_from_iter;
use serde::{Deserialize, Serialize};
use serde_json;
use serde_yaml::Value as YamlValue;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::path::{Path, PathBuf};
use tokio::fs;
use tokio::time::{sleep, Duration};
use tracing::instrument;

use super::run_dynamic;

type StepOverrides = HashMap<(String, String), String>;
type FlowOverrides = HashMap<String, String>;
type ParseOverridesResult = (StepOverrides, FlowOverrides, Option<String>, Vec<String>);

async fn wait_for_path(path: &Path, timeout_secs: u64) -> Result<()> {
    let start = std::time::Instant::now();
    let poll_interval = Duration::from_millis(500);
    let timeout = Duration::from_secs(timeout_secs);

    while start.elapsed() < timeout {
        if path.exists() {
            return Ok(());
        }
        sleep(poll_interval).await;
    }

    Err(anyhow!(
        "Timeout waiting for path to exist: {} (waited {}s)",
        path.display(),
        timeout_secs
    )
    .into())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowStepProgress {
    pub status: String, // "pending", "in_progress", "completed", "failed"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub started_at: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub completed_at: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowProgress {
    pub datasite: String,
    pub flow_name: String,
    pub run_id: String,
    pub started_at: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub completed_at: Option<String>,
    pub steps: HashMap<String, FlowStepProgress>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_step: Option<String>,
}

impl FlowProgress {
    fn new(datasite: &str, flow_name: &str, run_id: &str) -> Self {
        Self {
            datasite: datasite.to_string(),
            flow_name: flow_name.to_string(),
            run_id: run_id.to_string(),
            started_at: Utc::now().to_rfc3339(),
            completed_at: None,
            steps: HashMap::new(),
            current_step: None,
        }
    }

    fn mark_step_started(&mut self, step_id: &str) {
        self.current_step = Some(step_id.to_string());
        self.steps.insert(
            step_id.to_string(),
            FlowStepProgress {
                status: "in_progress".to_string(),
                started_at: Some(Utc::now().to_rfc3339()),
                completed_at: None,
                error: None,
            },
        );
    }

    fn mark_step_completed(&mut self, step_id: &str) {
        if let Some(progress) = self.steps.get_mut(step_id) {
            progress.status = "completed".to_string();
            progress.completed_at = Some(Utc::now().to_rfc3339());
        }
        if self.current_step.as_deref() == Some(step_id) {
            self.current_step = None;
        }
    }

    fn mark_step_failed(&mut self, step_id: &str, error: &str) {
        if let Some(progress) = self.steps.get_mut(step_id) {
            progress.status = "failed".to_string();
            progress.completed_at = Some(Utc::now().to_rfc3339());
            progress.error = Some(error.to_string());
        }
    }

    fn mark_flow_completed(&mut self) {
        self.completed_at = Some(Utc::now().to_rfc3339());
        self.current_step = None;
    }
}

/// Get the flow run directory
/// Structure: shared/flows/{flow_name}/{run_id}/
fn get_flow_run_dir(shared_base: &Path, flow_name: &str, run_id: &str) -> PathBuf {
    shared_base.join("flows").join(flow_name).join(run_id)
}

/// Get the progress directory for a flow run
/// Structure: shared/flows/{flow_name}/{run_id}/_progress/
fn get_progress_dir(shared_base: &Path, flow_name: &str, run_id: &str) -> PathBuf {
    get_flow_run_dir(shared_base, flow_name, run_id).join("_progress")
}

/// Get the state.json file path (current snapshot of flow/step status)
/// Structure: shared/flows/{flow_name}/{run_id}/_progress/state.json
fn get_state_file(shared_base: &Path, flow_name: &str, run_id: &str) -> PathBuf {
    get_progress_dir(shared_base, flow_name, run_id).join("state.json")
}

fn write_state(
    shared_base: &Path,
    flow_name: &str,
    run_id: &str,
    progress: &FlowProgress,
) -> Result<()> {
    let progress_dir = get_progress_dir(shared_base, flow_name, run_id);
    std::fs::create_dir_all(&progress_dir)?;

    // Write structured state.json
    let state_file = progress_dir.join("state.json");
    let json = serde_json::to_string_pretty(progress)?;
    std::fs::write(&state_file, json)?;
    Ok(())
}

/// Append a log entry to progress.json (JSONL format for event streaming)
fn append_progress(
    shared_base: &Path,
    flow_name: &str,
    run_id: &str,
    event: &str,
    step_id: Option<&str>,
    message: &str,
) -> Result<()> {
    let progress_dir = get_progress_dir(shared_base, flow_name, run_id);
    std::fs::create_dir_all(&progress_dir)?;

    let log_file = progress_dir.join("progress.json");
    let timestamp = Utc::now().to_rfc3339();
    let log_entry = serde_json::json!({
        "timestamp": timestamp,
        "event": event,
        "step": step_id,
        "message": message,
    });

    use std::io::Write;
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&log_file)?;
    writeln!(file, "{}", log_entry.to_string())?;
    Ok(())
}

fn read_state(progress_file: &Path) -> Result<FlowProgress> {
    let content = std::fs::read_to_string(progress_file)?;
    let progress: FlowProgress = serde_json::from_str(&content)?;
    Ok(progress)
}

async fn wait_for_step_completion(
    shared_base: &Path,
    flow_name: &str,
    run_id: &str,
    source_datasite: &str,
    step_id: &str,
    timeout_secs: u64,
    verbose: bool,
) -> Result<()> {
    let progress_file = get_state_file(shared_base, flow_name, run_id);
    let start = std::time::Instant::now();
    let poll_interval = Duration::from_millis(1000);
    let timeout = Duration::from_secs(timeout_secs);

    if verbose {
        println!(
            "⏳ Waiting for {} to complete step '{}'...",
            source_datasite, step_id
        );
        println!("   Looking for: {}", progress_file.display());
    }

    while start.elapsed() < timeout {
        if progress_file.exists() {
            if let Ok(progress) = read_state(&progress_file) {
                if let Some(step_progress) = progress.steps.get(step_id) {
                    match step_progress.status.as_str() {
                        "completed" => {
                            if verbose {
                                println!(
                                    "✓  {} completed step '{}'",
                                    source_datasite, step_id
                                );
                            }
                            return Ok(());
                        }
                        "failed" => {
                            return Err(anyhow!(
                                "{} failed step '{}': {}",
                                source_datasite,
                                step_id,
                                step_progress.error.as_deref().unwrap_or("unknown error")
                            )
                            .into());
                        }
                        _ => {}
                    }
                }
            }
        }
        sleep(poll_interval).await;
    }

    Err(anyhow!(
        "Timeout waiting for {} to complete step '{}' (waited {}s)",
        source_datasite,
        step_id,
        timeout_secs
    )
    .into())
}

const RESULTS_TABLE_PREFIX: &str = "z_results_";
const DEFAULT_COLUMN_PREFIX: &str = "col";
const DEFAULT_TABLE_FALLBACK: &str = "t";

trait DialoguerResultExt<T> {
    fn cli_result(self) -> Result<T>;
}

impl<T> DialoguerResultExt<T> for std::result::Result<T, dialoguer::Error> {
    fn cli_result(self) -> Result<T> {
        self.map_err(|err| crate::error::Error::Anyhow(anyhow::anyhow!(err)))
    }
}

#[instrument(skip_all, fields(component = "flow", flow_name = ?name), err)]
pub async fn create(
    file: Option<String>,
    name: Option<String>,
    uses: Option<String>,
    step_id: Option<String>,
) -> Result<()> {
    let flow_path = resolve_flow_path(file);
    let db = BioVaultDb::new().ok();

    if flow_path.exists()
        && !Confirm::with_theme(&ColorfulTheme::default())
            .with_prompt(format!(
                "Flow file {} already exists. Overwrite?",
                flow_path.display()
            ))
            .default(false)
            .interact()
            .cli_result()?
    {
        println!("✋ Aborting flow creation.");
        return Ok(());
    }

    if let Some(uses_value) = uses {
        let flow_dir = flow_path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));

        let (reference, root) = normalize_module_reference(&uses_value, flow_dir)?;
        let choice = ModuleChoice::Path {
            reference: reference.clone(),
            root: root.clone(),
        };
        let module = load_module_spec(&choice)?;

        let mut spec = FlowSpec::default();
        spec.name = name.unwrap_or_else(|| module.spec.name.clone());

        let generated_id = generate_default_step_id(&module, &spec);
        let step_id_value = step_id.unwrap_or(generated_id);

        let mut with_map = BTreeMap::new();
        for input in &module.spec.inputs {
            let key = ensure_flow_input(&mut spec, &step_id_value, input);
            with_map.insert(
                input.name.clone(),
                YamlValue::String(format!("inputs.{}", key)),
            );
        }

        let publish_map = default_publish_map(&module.spec.outputs);

        spec.steps.push(FlowStepSpec {
            id: step_id_value,
            uses: Some(reference),
            where_exec: None,
            foreach: None,
            runs_on: None,
            order: None,
            coordination: None,
            with: with_map,
            publish: publish_map,
            share: BTreeMap::new(),
            store: BTreeMap::new(),
            permissions: Vec::new(),
            barrier: None,
        });

        spec.save(&flow_path)?;

        println!(
            "\n✅ Saved flow to {}",
            flow_path.display().to_string().bold()
        );

        let validation = validate_internal(&flow_path, &spec, db.as_ref())?;
        print_validation(&validation, true);

        return Ok(());
    }

    let wizard = FlowWizard::new(&flow_path, db.as_ref());
    let spec = wizard.run()?;
    spec.save(&flow_path)?;

    println!(
        "\n✅ Saved flow to {}",
        flow_path.display().to_string().bold()
    );

    let validation = validate_internal(&flow_path, &spec, db.as_ref())?;
    print_validation(&validation, true);

    Ok(())
}

pub async fn add_step(file: Option<String>) -> Result<()> {
    let flow_path = resolve_flow_path(file);
    if !flow_path.exists() {
        return Err(anyhow!(
            "Flow file not found: {}. Run 'bv flow create' first.",
            flow_path.display()
        )
        .into());
    }

    let mut spec = FlowSpec::load(&flow_path)?;
    let db = BioVaultDb::new().ok();

    let wizard = FlowWizard::new(&flow_path, db.as_ref());
    spec = wizard.add_step_to(spec)?;
    spec.save(&flow_path)?;

    println!(
        "\n✅ Updated flow at {}",
        flow_path.display().to_string().bold()
    );

    let validation = validate_internal(&flow_path, &spec, db.as_ref())?;
    print_validation(&validation, true);

    Ok(())
}

#[instrument(skip(extra_args), fields(component = "flow", flow = %flow_path, dry_run = %dry_run, resume = %resume), err)]
pub async fn run_flow(
    flow_path: &str,
    extra_args: Vec<String>,
    dry_run: bool,
    resume: bool,
    results_dir_override: Option<String>,
) -> Result<()> {
    let path = Path::new(flow_path);
    if !path.exists() {
        return Err(anyhow!("Flow file not found: {}", flow_path).into());
    }

    let (step_overrides, flow_overrides, explicit_results_dir, nextflow_passthrough) =
        parse_overrides(&extra_args, results_dir_override.clone())?;

    let spec = FlowSpec::load(path)?;
    let mut db = BioVaultDb::new().ok();
    let validation = validate_internal(path, &spec, db.as_ref())?;

    if validation.has_errors() {
        print_validation(&validation, true);
        return Err(anyhow!("Flow has validation errors. Fix them before running.").into());
    }

    if !validation.unresolved_steps.is_empty() {
        return Err(anyhow!(
            "Flow contains steps without a module. Set 'uses' or register the module before running."
        )
        .into());
    }

    let flow_dir = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    // Use shared run_id from env var (for distributed coordination) or generate new one
    let run_id = env::var("BIOVAULT_FLOW_RUN_ID")
        .ok()
        .filter(|s| !s.trim().is_empty())
        .unwrap_or_else(|| Utc::now().format("%Y%m%d%H%M%S").to_string());
    let run_msg = format!("🆔 Flow run {}", run_id);
    println!("{}", run_msg);
    append_desktop_log(&run_msg);

    for key in flow_overrides.keys() {
        if !spec.inputs.contains_key(key) {
            return Err(anyhow!("Unknown flow input '{}' in --set", key).into());
        }
    }

    let mut resolved_inputs: HashMap<String, String> = HashMap::new();
    for (name, input_spec) in &spec.inputs {
        if let Some(value) = flow_overrides.get(name) {
            let resolved = literal_to_value(value, input_spec.raw_type())?;
            resolved_inputs.insert(name.clone(), resolved);
        } else if let Some(default_literal) = input_spec.default_literal() {
            let resolved = literal_to_value(default_literal, input_spec.raw_type())?;
            resolved_inputs.insert(name.clone(), resolved);
        }
    }

    let mut step_outputs: HashMap<String, HashMap<String, HashMap<String, String>>> =
        HashMap::new();
    let mut step_output_order: HashMap<String, Vec<String>> = HashMap::new();

    let requested_results_dir = explicit_results_dir.or(results_dir_override);

    let base_results_dir = match requested_results_dir {
        Some(dir) => PathBuf::from(dir),
        None => {
            let mut base = PathBuf::from("results/flows");
            base.push(&spec.name);
            base
        }
    };
    let base_results_dir_abs = if base_results_dir.is_absolute() {
        base_results_dir.clone()
    } else {
        env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(&base_results_dir)
    };

    if !dry_run {
        fs::create_dir_all(&base_results_dir_abs).await?;
    }

    let current_datasite = resolve_current_datasite();
    let run_all_targets = env::var("BIOVAULT_FLOW_RUN_ALL")
        .map(|value| !value.trim().is_empty())
        .unwrap_or(false);
    let syftbox_data_dir = crate::config::Config::load()
        .ok()
        .and_then(|cfg| cfg.get_syftbox_data_dir().ok());
    let sandbox_root = env::var("BIOVAULT_SANDBOX_ROOT")
        .ok()
        .map(PathBuf::from)
        .and_then(|p| p.canonicalize().ok());

    // Progress tracking for distributed coordination
    let verbose_progress = env::var("BIOVAULT_FLOW_VERBOSE")
        .map(|v| !v.trim().is_empty())
        .unwrap_or(false);
    let await_timeout_secs: u64 = env::var("BIOVAULT_FLOW_AWAIT_TIMEOUT")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(120);

    // Determine shared base for progress logs
    // Each party writes to THEIR OWN shared folder
    // When reading, we look in the specific party's folder
    let progress_write_base: Option<PathBuf> = if let Some(ref root) = sandbox_root {
        if run_all_targets {
            // In sandbox mode with run_all, use first datasite's folder
            let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);
            if let Some(first) = all_datasites.first() {
                Some(
                    root.join(first)
                        .join("datasites")
                        .join(first)
                        .join("shared"),
                )
            } else {
                None
            }
        } else if let Some(ref ds) = current_datasite {
            // In sandbox mode with specific datasite, write to OUR OWN folder
            Some(
                root.join(ds)
                    .join("datasites")
                    .join(ds)
                    .join("shared"),
            )
        } else {
            None
        }
    } else if let Some(ref ds) = current_datasite {
        // In distributed mode, write to OUR OWN shared folder
        if let Some(ref data_dir) = syftbox_data_dir {
            Some(data_dir.join("datasites").join(ds).join("shared"))
        } else {
            None
        }
    } else {
        None
    };

    // For reading other parties' progress, we need to know where to look
    // In sandbox mode: read directly from source_datasite's own folder (no sync needed for tests)
    // In normal mode: use syftbox_data_dir/datasites/source_datasite/shared (synced via SyftBox)
    let get_progress_read_base = |_reader_datasite: &str, source_datasite: &str| -> Option<PathBuf> {
        if let Some(ref root) = sandbox_root {
            // In sandbox mode, read directly from the source's own shared folder
            // This bypasses the need for SyftBox sync in test environments
            let path = root.join(source_datasite).join("datasites").join(source_datasite).join("shared");
            eprintln!("DEBUG get_progress_read_base: sandbox_root={:?}, source={}, path={}", root, source_datasite, path.display());
            Some(path)
        } else if let Some(ref data_dir) = syftbox_data_dir {
            // In normal mode, syftbox syncs files to our local view
            let path = data_dir.join("datasites").join(source_datasite).join("shared");
            eprintln!("DEBUG get_progress_read_base: syftbox_data_dir={:?}, source={}, path={}", data_dir, source_datasite, path.display());
            Some(path)
        } else {
            eprintln!("DEBUG get_progress_read_base: no sandbox_root or syftbox_data_dir");
            None
        }
    };

    // Initialize progress tracking for this datasite (if we have a current datasite)
    let mut flow_progress: Option<FlowProgress> = if !run_all_targets {
        current_datasite
            .as_ref()
            .map(|ds| FlowProgress::new(ds, &spec.name, &run_id))
    } else {
        None
    };

    // Write initial progress
    if let (Some(ref progress), Some(ref shared_base)) = (&flow_progress, &progress_write_base) {
        if let Err(e) = write_state(shared_base, &spec.name, &run_id, progress) {
            if verbose_progress {
                println!("⚠️  Could not write initial progress: {}", e);
            }
        } else if verbose_progress {
            let progress_dir = get_progress_dir(shared_base, &spec.name, &run_id);
            println!("📝 Writing progress to: {}", progress_dir.display());
        }
    }

    // Handle top-level coordination setup (before any steps run)
    if let Some(ref coordination) = spec.coordination {
        let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);

        // Determine who can read the coordination files
        let readers: Vec<String> = match &coordination.share_with {
            crate::flow_spec::FlowCoordinationShareWith::All(s) if s == "all" => {
                vec!["{datasites[*]}".to_string()]
            }
            crate::flow_spec::FlowCoordinationShareWith::All(s) => vec![s.clone()],
            crate::flow_spec::FlowCoordinationShareWith::List(list) => list.clone(),
        };

        // Coordination runs on all datasites by default
        let coord_targets: Vec<String> = if let Some(ref targets) = coordination.targets {
            match targets {
                crate::flow_spec::FlowRunTargets::One(s) if s == "all" => all_datasites.clone(),
                crate::flow_spec::FlowRunTargets::One(s) => vec![s.clone()],
                crate::flow_spec::FlowRunTargets::Many(v) => v.clone(),
                crate::flow_spec::FlowRunTargets::Selector(sel) => sel.include.clone(),
            }
        } else {
            all_datasites.clone()
        };

        let coord_run_targets =
            resolve_run_targets(&coord_targets, current_datasite.as_deref(), run_all_targets);

        for target in coord_run_targets {
            let datasite_key = target.clone().unwrap_or_else(|| "local".to_string());
            let target_data_dir = if let (Some(ref root), Some(ref site)) = (&sandbox_root, &target)
            {
                root.join(site)
            } else {
                syftbox_data_dir.clone().ok_or_else(|| {
                    anyhow!("Coordination requires SYFTBOX_DATA_DIR or BIOVAULT_SANDBOX_ROOT")
                })?
            };

            println!(
                "\n🤝 Setting up coordination for {}",
                datasite_key.bold()
            );

            // Create coordination folder and permissions
            let url_template = coordination.url_template();
            let progress_url = render_flow_template(
                url_template,
                &datasite_key,
                &all_datasites,
                &run_id,
                &spec.name,
                None,
                None,
                &spec.vars,
            );

            let permission_spec = FlowPermissionSpec {
                url: progress_url,
                rules: vec![FlowPermissionRuleSpec {
                    pattern: "**".to_string(),
                    access: FlowPermissionAccessSpec {
                        read: readers.clone(),
                        write: vec!["{datasite.current}".to_string()],
                        admin: vec!["{datasite.current}".to_string()],
                    },
                }],
            };

            execute_permissions_step(
                &[permission_spec],
                &target_data_dir,
                &all_datasites,
                target.as_deref(),
                &spec.name,
                &run_id,
            )?;
        }
    }

    // Handle MPC setup (create communication channels between parties)
    if let Some(ref mpc) = spec.mpc {
        let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);
        let party_count = all_datasites.len();

        println!(
            "\n🔐 Setting up MPC channels (topology: {}, {} parties)",
            mpc.topology,
            party_count
        );

        // Generate channels based on topology
        let channels = mpc.generate_channels(party_count);
        let channel_pattern = mpc.channel_pattern();

        // MPC runs on all datasites
        let mpc_run_targets =
            resolve_run_targets(&all_datasites, current_datasite.as_deref(), run_all_targets);

        for target in mpc_run_targets {
            let datasite_key = target.clone().unwrap_or_else(|| "local".to_string());
            let target_data_dir = if let (Some(ref root), Some(ref site)) = (&sandbox_root, &target)
            {
                root.join(site)
            } else {
                syftbox_data_dir.clone().ok_or_else(|| {
                    anyhow!("MPC setup requires SYFTBOX_DATA_DIR or BIOVAULT_SANDBOX_ROOT")
                })?
            };

            // Find this datasite's party index
            let my_party_index = all_datasites
                .iter()
                .position(|ds| ds == &datasite_key)
                .unwrap_or(0);

            // Create outgoing channels (where this party is the sender)
            let my_outgoing: Vec<_> = channels
                .iter()
                .filter(|(from, _)| *from == my_party_index)
                .collect();

            for (from_idx, to_idx) in my_outgoing {
                let channel_name = channel_pattern
                    .replace("{from}", &from_idx.to_string())
                    .replace("{to}", &to_idx.to_string());

                let receiver_datasite = &all_datasites[*to_idx];

                // Create channel folder with permissions
                let url_template = mpc.url_template();
                let base_url = render_flow_template(
                    url_template,
                    &datasite_key,
                    &all_datasites,
                    &run_id,
                    &spec.name,
                    None,
                    None,
                    &spec.vars,
                );
                let channel_url = format!("{}/{}", base_url, channel_name);

                println!(
                    "  📡 {}@{} → {} (channel: {})",
                    my_party_index,
                    datasite_key,
                    receiver_datasite,
                    channel_name
                );

                let permission_spec = FlowPermissionSpec {
                    url: channel_url,
                    rules: vec![FlowPermissionRuleSpec {
                        pattern: "**".to_string(),
                        access: FlowPermissionAccessSpec {
                            read: vec![datasite_key.clone(), receiver_datasite.clone()],
                            write: vec![datasite_key.clone(), receiver_datasite.clone()],
                            admin: vec![datasite_key.clone()],
                        },
                    }],
                };

                execute_permissions_step(
                    &[permission_spec],
                    &target_data_dir,
                    &all_datasites,
                    target.as_deref(),
                    &spec.name,
                    &run_id,
                )?;
            }
        }

        println!("  ✓ MPC channels configured");
    }

    for (step_index, step) in spec.steps.iter().enumerate() {
        let step_number = step_index + 1; // 1-indexed for user display
        let module = resolve_module(step, flow_dir, db.as_ref())?;
        let is_permissions_only = module.is_none() && !step.permissions.is_empty();
        let is_coordination = step.coordination.is_some();
        let is_barrier = step.barrier.is_some();

        if module.is_none() && step.permissions.is_empty() && !is_coordination && !is_barrier {
            return Err(anyhow!(
                "Step '{}' has no module to run, no permissions to create, no coordination, and no barrier",
                step.id
            )
            .into());
        }

        let step_targets = resolve_step_targets(step);
        if !step_targets.is_empty() {
            step_output_order.insert(step.id.clone(), step_targets.clone());
        }
        if !step.share.is_empty() && step_targets.is_empty() {
            return Err(anyhow!(
                "Step '{}' declares shares but has no runs_on/foreach targets",
                step.id
            )
            .into());
        }

        let run_targets =
            resolve_run_targets(&step_targets, current_datasite.as_deref(), run_all_targets);

        let mut share_locations_by_target: HashMap<String, HashMap<String, ShareLocation>> =
            HashMap::new();
        if !step.share.is_empty() {
            for (share_name, share_spec) in &step.share {
                for target in &step_targets {
                    let data_dir = if let Some(ref root) = sandbox_root {
                        root.join(target)
                    } else {
                        syftbox_data_dir.clone().ok_or_else(|| {
                            anyhow!(
                                "Step '{}' requires SYFTBOX_DATA_DIR or BIOVAULT_SANDBOX_ROOT to resolve shared paths",
                                step.id
                            )
                        })?
                    };
                    let share_location = resolve_share_location(
                        &share_spec.url,
                        target,
                        &step_targets,
                        &run_id,
                        &spec.name,
                        &data_dir,
                        step_number,
                        &step.id,
                        &spec.vars,
                    )?;
                    share_locations_by_target
                        .entry(target.clone())
                        .or_default()
                        .insert(share_name.clone(), share_location);
                }
            }
        }

        if matches!(step.order.as_deref(), Some("parallel")) && run_targets.len() > 1 {
            println!(
                "⚠️  Step '{}' requested parallel execution; running sequentially.",
                step.id
            );
        }

        if run_targets.is_empty() {
            if !share_locations_by_target.is_empty() {
                let step_entry = step_outputs.entry(step.id.clone()).or_default();
                for (target, share_outputs) in &share_locations_by_target {
                    let outputs = step_entry.entry(target.clone()).or_default();
                    for (name, location) in share_outputs {
                        outputs
                            .entry(name.clone())
                            .or_insert_with(|| location.syft_url.clone());
                    }
                }
            }
            continue;
        }

        // Handle permissions-only steps separately
        if is_permissions_only {
            let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);
            for target in run_targets {
                let datasite_key = target.clone().unwrap_or_else(|| "local".to_string());
                let target_data_dir = if let (Some(ref root), Some(ref site)) = (&sandbox_root, &target) {
                    root.join(site)
                } else {
                    syftbox_data_dir.clone().ok_or_else(|| {
                        anyhow!("Permissions step requires SYFTBOX_DATA_DIR or BIOVAULT_SANDBOX_ROOT")
                    })?
                };

                let step_label = format!("{}@{}", step.id, datasite_key);
                println!("\n🔐 Creating permissions for step {}", step_label.bold());

                execute_permissions_step(
                    &step.permissions,
                    &target_data_dir,
                    &all_datasites,
                    target.as_deref(),
                    &spec.name,
                    &run_id,
                )?;
            }
            continue;
        }

        // Handle coordination steps - sets up _progress folder with permissions
        if is_coordination {
            let coordination = step.coordination.as_ref().unwrap();
            let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);

            // Determine targets from coordination spec or fall back to step runs_on targets
            let coord_targets: Vec<String> = if let Some(ref targets) = coordination.targets {
                match targets {
                    crate::flow_spec::FlowRunTargets::One(s) if s == "all" => {
                        all_datasites.clone()
                    }
                    crate::flow_spec::FlowRunTargets::One(s) => vec![s.clone()],
                    crate::flow_spec::FlowRunTargets::Many(v) => v.clone(),
                    crate::flow_spec::FlowRunTargets::Selector(sel) => sel.include.clone(),
                }
            } else if let Some(ref runs_on) = step.runs_on {
                runs_on.as_vec()
            } else {
                all_datasites.clone()
            };

            // Determine who can read the coordination files
            let readers: Vec<String> = match &coordination.share_with {
                crate::flow_spec::FlowCoordinationShareWith::All(s) if s == "all" => {
                    vec!["{datasites[*]}".to_string()]
                }
                crate::flow_spec::FlowCoordinationShareWith::All(s) => vec![s.clone()],
                crate::flow_spec::FlowCoordinationShareWith::List(list) => list.clone(),
            };

            let coord_run_targets =
                resolve_run_targets(&coord_targets, current_datasite.as_deref(), run_all_targets);

            for target in coord_run_targets {
                let datasite_key = target.clone().unwrap_or_else(|| "local".to_string());
                let target_data_dir = if let (Some(ref root), Some(ref site)) = (&sandbox_root, &target) {
                    root.join(site)
                } else {
                    syftbox_data_dir.clone().ok_or_else(|| {
                        anyhow!("Coordination step requires SYFTBOX_DATA_DIR or BIOVAULT_SANDBOX_ROOT")
                    })?
                };

                let step_label = format!("{}@{}", step.id, datasite_key);
                println!("\n🤝 Setting up coordination for {}", step_label.bold());

                // Create coordination folder and permissions
                let url_template = coordination.url_template();
                let progress_url = render_flow_template(
                    url_template,
                    &datasite_key,
                    &all_datasites,
                    &run_id,
                    &spec.name,
                    None,
                    None,
                    &spec.vars,
                );

                let permission_spec = FlowPermissionSpec {
                    url: progress_url,
                    rules: vec![FlowPermissionRuleSpec {
                        pattern: "**".to_string(),
                        access: FlowPermissionAccessSpec {
                            read: readers.clone(),
                            write: vec!["{datasite.current}".to_string()],
                            admin: vec!["{datasite.current}".to_string()],
                        },
                    }],
                };

                execute_permissions_step(
                    &[permission_spec],
                    &target_data_dir,
                    &all_datasites,
                    target.as_deref(),
                    &spec.name,
                    &run_id,
                )?;
            }
            continue;
        }

        // Handle barrier steps - wait for specified step to complete on specified targets
        if is_barrier {
            let barrier = step.barrier.as_ref().unwrap();
            let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);

            // Build groups map for target resolution
            let mut groups: std::collections::BTreeMap<String, Vec<String>> = std::collections::BTreeMap::new();
            groups.insert("all".to_string(), all_datasites.clone());
            if let Some(first) = all_datasites.first() {
                groups.insert("aggregator".to_string(), vec![first.clone()]);
            }
            if all_datasites.len() > 1 {
                groups.insert("clients".to_string(), all_datasites[1..].to_vec());
            }

            // Determine which targets to wait for (resolve groups)
            let barrier_targets: Vec<String> = if let Some(ref targets) = barrier.targets {
                match targets {
                    crate::flow_spec::FlowRunTargets::One(s) if s == "all" => all_datasites.clone(),
                    crate::flow_spec::FlowRunTargets::One(s) => {
                        // Check if it's a group name
                        groups.get(s).cloned().unwrap_or_else(|| vec![s.clone()])
                    },
                    crate::flow_spec::FlowRunTargets::Many(v) => v.clone(),
                    crate::flow_spec::FlowRunTargets::Selector(sel) => sel.include.clone()
                }
            } else {
                // Default to all datasites if not specified
                all_datasites.clone()
            };

            let timeout_secs = barrier.timeout.unwrap_or(300); // Default 5 minutes
            let wait_for_step = &barrier.wait_for;

            println!(
                "\n⏳ Barrier '{}': waiting for step '{}' to complete on {} target(s)...",
                step.id.bold(),
                wait_for_step.bold(),
                barrier_targets.len()
            );

            // Skip waiting for ourselves
            let targets_to_wait: Vec<&String> = barrier_targets
                .iter()
                .filter(|t| Some(t.as_str()) != current_datasite.as_deref())
                .collect();

            if targets_to_wait.is_empty() {
                println!("   ✓ No other targets to wait for (only self)");
            } else {
                for target in &targets_to_wait {
                    println!("   • Waiting for {} to complete '{}'...", target, wait_for_step);
                }

                // Wait for each target to complete the specified step
                for target in targets_to_wait {
                    if let Some(reader) = current_datasite.as_deref() {
                        if let Some(read_base) = get_progress_read_base(reader, target) {
                            wait_for_step_completion(
                                &read_base,
                                &spec.name,
                                &run_id,
                                target,
                                wait_for_step,
                                timeout_secs,
                                true, // verbose
                            )
                            .await?;
                        } else {
                            return Err(anyhow!(
                                "Barrier '{}': Cannot determine progress path for {}",
                                step.id,
                                target
                            ).into());
                        }
                    }
                }

                println!("   ✓ All {} target(s) completed step '{}'", barrier_targets.len(), wait_for_step);
            }

            // Update progress for this barrier step
            if let Some(ref mut progress) = flow_progress {
                progress.mark_step_completed(&step.id);
                if let Some(ref shared_base) = progress_write_base {
                    if let Err(e) = write_state(shared_base, &spec.name, &run_id, progress) {
                        if verbose_progress {
                            println!("⚠️  Could not write barrier progress: {}", e);
                        }
                    }
                    if let Err(e) = append_progress(
                        shared_base,
                        &spec.name,
                        &run_id,
                        "barrier_completed",
                        Some(&step.id),
                        &format!("Barrier {} completed (waited for {})", step.id, wait_for_step),
                    ) {
                        if verbose_progress {
                            println!("⚠️  Could not append barrier progress: {}", e);
                        }
                    }
                }
            }

            continue;
        }

        // At this point we have a module (permissions-only, coordination, and barrier steps already continued)
        let module = module.unwrap();
        let module_root = module.root.to_string_lossy().to_string();
        let module_spec = &module.spec;

        for target in run_targets {
            let current_datasite = target.as_deref();
            let target_data_dir = if let (Some(ref root), Some(site)) = (&sandbox_root, &target) {
                Some(root.join(site))
            } else {
                syftbox_data_dir.clone()
            };

            // Wait for upstream steps from other datasites if not running all targets
            if !run_all_targets {
                for input in &module_spec.inputs {
                    if let Some(binding) = step.with.get(&input.name).and_then(value_to_string) {
                        if binding.starts_with("step.") {
                            let parts: Vec<&str> = binding.split('.').collect();
                            // Handle both outputs and share namespaces
                            let is_valid_namespace = parts.len() >= 4
                                && (parts[2] == "outputs" || parts[2] == "share");
                            // Check if this is a url_list binding (needs all sites)
                            let is_url_list = parts.len() >= 5
                                && matches!(parts[4], "url_list" | "manifest" | "all" | "paths");

                            if is_valid_namespace {
                                let source_step = parts[1];
                                let output_name = parts[3];

                                if let Some(outputs_by_site) = step_outputs.get(source_step) {
                                    // Check if current datasite has its own output from this step
                                    let current_has_output = current_datasite
                                        .and_then(|site| outputs_by_site.get(site))
                                        .map(|outputs| outputs.contains_key(output_name))
                                        .unwrap_or(false);

                                    // Determine which sites to wait for:
                                    // - For url_list bindings: wait for ALL other sites
                                    // - For direct bindings: only wait if current doesn't have output
                                    let sites_to_wait: Vec<&String> = if is_url_list {
                                        // Need all sites for manifest
                                        outputs_by_site
                                            .keys()
                                            .filter(|site| Some(site.as_str()) != current_datasite)
                                            .collect()
                                    } else if !current_has_output {
                                        // Current site doesn't have output, wait for source sites
                                        outputs_by_site.keys().collect()
                                    } else {
                                        // Current site has its own output, no waiting needed
                                        Vec::new()
                                    };

                                    // Wait for step completion via progress logs
                                    for site in &sites_to_wait {
                                        if Some(site.as_str()) != current_datasite {
                                            // Look in the current datasite's synced view of the source datasite's folder
                                            if let Some(reader) = current_datasite {
                                                if let Some(read_base) = get_progress_read_base(reader, site) {
                                                    // Wait for this datasite to complete the source step
                                                    wait_for_step_completion(
                                                        &read_base,
                                                        &spec.name,
                                                        &run_id,
                                                        site,
                                                        source_step,
                                                        await_timeout_secs,
                                                        verbose_progress,
                                                    )
                                                    .await?;
                                                }
                                            }
                                        }
                                    }

                                    // Also wait for the actual files (belt and suspenders)
                                    let mut paths_to_wait: Vec<PathBuf> = Vec::new();
                                    for site in &sites_to_wait {
                                        if let Some(outputs) = outputs_by_site.get(*site) {
                                            if let Some(value) = outputs.get(output_name) {
                                                if value.starts_with("syft://") {
                                                    if let Ok(resolved) = resolve_syft_path(
                                                        value,
                                                        target_data_dir.as_ref(),
                                                        sandbox_root.as_ref(),
                                                    ) {
                                                        paths_to_wait.push(PathBuf::from(resolved));
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    for path in &paths_to_wait {
                                        if !path.exists() {
                                            if verbose_progress {
                                                println!(
                                                    "⏳ Waiting for file: {}",
                                                    path.display()
                                                );
                                            }
                                            wait_for_path(path, await_timeout_secs).await?;
                                            if verbose_progress {
                                                println!("✓  File available: {}", path.display());
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            let mut step_args = Vec::new();

            // Build groups map for conditional input resolution
            // Convention: "aggregator" = first datasite, "clients" = rest
            let datasite_groups: std::collections::BTreeMap<String, Vec<String>> = {
                let mut groups = std::collections::BTreeMap::new();
                if !step_targets.is_empty() {
                    groups.insert("all".to_string(), step_targets.clone());
                    if let Some(first) = step_targets.first() {
                        groups.insert("aggregator".to_string(), vec![first.clone()]);
                    }
                    if step_targets.len() > 1 {
                        groups.insert("clients".to_string(), step_targets[1..].to_vec());
                    }
                }
                groups
            };

            let datasite_key = target.clone().unwrap_or_else(|| "local".to_string());

            for input in &module_spec.inputs {
                // Parse as conditional input to check only/without constraints
                let with_value = step.with.get(&input.name);

                let conditional = with_value.and_then(ConditionalInput::from_yaml);

                // Check if this input applies to the current target
                if let Some(ref cond) = conditional {
                    if !cond.applies_to(&datasite_key, &datasite_groups) {
                        // Input doesn't apply to this target, skip it
                        continue;
                    }
                }

                // Get the binding value (either from conditional or simple value)
                let binding = conditional
                    .map(|c| c.value)
                    .or_else(|| with_value.and_then(value_to_string))
                    .ok_or_else(|| {
                        anyhow!(
                            "Step '{}' is missing a binding for input '{}'",
                            step.id,
                            input.name
                        )
                    })?;

                let override_key = (step.id.clone(), input.name.clone());
                let resolved_binding = step_overrides
                    .get(&override_key)
                    .cloned()
                    .unwrap_or(binding);

                // Render template variables in the binding (e.g., {datasite.index}, {datasites})
                let rendered_binding = render_flow_template(
                    &resolved_binding,
                    &datasite_key,
                    &step_targets,
                    &run_id,
                    &spec.name,
                    Some(step_number),
                    Some(&step.id),
                    &spec.vars,
                );

                let resolved_value = resolve_binding(
                    &rendered_binding,
                    &input.raw_type,
                    &spec.inputs,
                    &resolved_inputs,
                    &step_outputs,
                    &step.id,
                    current_datasite,
                    &base_results_dir_abs,
                    &step_output_order,
                    target_data_dir.as_ref(),
                    sandbox_root.as_ref(),
                )?;

                step_args.push(format!("--{}", input.name));
                step_args.push(resolved_value);
            }

            for param in &module_spec.parameters {
                let override_key = (step.id.clone(), param.name.clone());
                if let Some(value) = step_overrides.get(&override_key) {
                    step_args.push("--set".to_string());
                    step_args.push(format!("params.{}={}", param.name, value));
                }
            }

            if !nextflow_passthrough.is_empty() {
                step_args.extend(nextflow_passthrough.clone());
            }

            let step_results_dir = if let Some(site) = &target {
                base_results_dir_abs
                    .join(&step.id)
                    .join(sanitize_identifier(site))
            } else {
                base_results_dir_abs.join(&step.id)
            };
            let step_results_dir_str = step_results_dir.to_string_lossy().to_string();

            let previous_override = env::var("BIOVAULT_DATASITE_OVERRIDE").ok();
            let previous_datasites_override = env::var("BIOVAULT_DATASITES_OVERRIDE").ok();
            let previous_syftbox_data_dir = env::var("SYFTBOX_DATA_DIR").ok();
            let previous_run_id = env::var("BV_RUN_ID").ok();
            let previous_flow_name = env::var("BV_FLOW_NAME").ok();
            match &target {
                Some(site) => env::set_var("BIOVAULT_DATASITE_OVERRIDE", site),
                None => env::remove_var("BIOVAULT_DATASITE_OVERRIDE"),
            }
            if step_targets.is_empty() {
                env::remove_var("BIOVAULT_DATASITES_OVERRIDE");
            } else {
                env::set_var("BIOVAULT_DATASITES_OVERRIDE", step_targets.join(","));
            }
            if let Some(ref data_dir) = target_data_dir {
                env::set_var("SYFTBOX_DATA_DIR", data_dir.to_string_lossy().as_ref());
            }
            // Set flow context for MPC coordination
            env::set_var("BV_RUN_ID", &run_id);
            env::set_var("BV_FLOW_NAME", &spec.name);

            let step_label = match &target {
                Some(site) => format!("{}@{}", step.id, site),
                None => step.id.clone(),
            };
            println!(
                "\n▶️  Running step {} using module {}",
                step_label.bold(),
                module_root
            );

            // Mark step as started in progress log
            if let Some(ref mut progress) = flow_progress {
                progress.mark_step_started(&step.id);
                if let Some(ref shared_base) = progress_write_base {
                    if let Err(e) = write_state(shared_base, &spec.name, &run_id, progress) {
                        if verbose_progress {
                            println!("⚠️  Could not write step start progress: {}", e);
                        }
                    }
                    let _ = append_progress(
                        shared_base,
                        &spec.name,
                        &run_id,
                        "step_started",
                        Some(&step.id),
                        &format!("Starting step {} using module {}", step_label, module_root),
                    );
                }
            }

            append_desktop_log(&format!(
                "[Flow] Running step {} with args: {:?}",
                step_label, step_args
            ));

            let step_result = run_dynamic::execute_dynamic(
                &module_root,
                step_args.clone(),
                dry_run,
                resume,
                Some(step_results_dir_str.clone()),
                run_dynamic::RunSettings::default(),
            )
            .await;

            // Mark step as completed or failed in progress log
            if let Some(ref mut progress) = flow_progress {
                let (event, msg) = match &step_result {
                    Ok(_) => {
                        progress.mark_step_completed(&step.id);
                        ("step_completed", format!("Step {} completed successfully", step_label))
                    }
                    Err(e) => {
                        progress.mark_step_failed(&step.id, &e.to_string());
                        ("step_failed", format!("Step {} failed: {}", step_label, e))
                    }
                };
                if let Some(ref shared_base) = progress_write_base {
                    if let Err(e) = write_state(shared_base, &spec.name, &run_id, progress) {
                        if verbose_progress {
                            println!("⚠️  Could not write step completion progress: {}", e);
                        }
                    }
                    let _ = append_progress(
                        shared_base,
                        &spec.name,
                        &run_id,
                        event,
                        Some(&step.id),
                        &msg,
                    );
                }
            }

            step_result?;

            append_desktop_log(&format!("[Flow] Completed step {}", step_label));

            match previous_override {
                Some(prev) => env::set_var("BIOVAULT_DATASITE_OVERRIDE", prev),
                None => env::remove_var("BIOVAULT_DATASITE_OVERRIDE"),
            }
            match previous_datasites_override {
                Some(prev) => env::set_var("BIOVAULT_DATASITES_OVERRIDE", prev),
                None => env::remove_var("BIOVAULT_DATASITES_OVERRIDE"),
            }
            match previous_syftbox_data_dir {
                Some(prev) => env::set_var("SYFTBOX_DATA_DIR", prev),
                None => env::remove_var("SYFTBOX_DATA_DIR"),
            }
            match previous_run_id {
                Some(prev) => env::set_var("BV_RUN_ID", prev),
                None => env::remove_var("BV_RUN_ID"),
            }
            match previous_flow_name {
                Some(prev) => env::set_var("BV_FLOW_NAME", prev),
                None => env::remove_var("BV_FLOW_NAME"),
            }

            let mut outputs = HashMap::new();
            for output in &module_spec.outputs {
                let path = output
                    .path
                    .as_ref()
                    .map(|p| step_results_dir.join(p))
                    .unwrap_or_else(|| step_results_dir.join(&output.name));
                outputs.insert(output.name.clone(), path.to_string_lossy().to_string());
            }

            // Ensure publish aliases reference concrete paths
            for alias in step.publish.keys() {
                if let Some(path) = outputs.get(alias) {
                    outputs.insert(alias.clone(), path.clone());
                }
            }

            let datasite_key = target.unwrap_or_default();
            if let Some(shared_outputs) = share_locations_by_target.get(&datasite_key) {
                for (name, location) in shared_outputs {
                    outputs.insert(name.clone(), location.syft_url.clone());
                }
            }

            if !step.share.is_empty() && !datasite_key.is_empty() && !dry_run {
                for (share_name, share_spec) in &step.share {
                    let share_location = share_locations_by_target
                        .get(&datasite_key)
                        .and_then(|paths| paths.get(share_name))
                        .ok_or_else(|| {
                            anyhow!(
                                "Shared output '{}' not resolved for step '{}' datasite '{}'",
                                share_name,
                                step.id,
                                datasite_key
                            )
                        })?;
                    let source_name = parse_share_source(&share_spec.source)?;
                    let source_path = outputs.get(&source_name).ok_or_else(|| {
                        anyhow!(
                            "Shared output '{}' references unknown source '{}' in step '{}'",
                            share_name,
                            share_spec.source,
                            step.id
                        )
                    })?;
                    let all_datasites = resolve_all_datasites_from_spec(&spec, &resolved_inputs);
                    create_shared_file(
                        Path::new(source_path),
                        &share_location.local_path,
                        share_spec,
                        &all_datasites,
                        Some(datasite_key.as_str()),
                        sandbox_root.as_ref(),
                    )?;
                }
            }

            step_outputs
                .entry(step.id.clone())
                .or_default()
                .insert(datasite_key, outputs.clone());

            if !step.store.is_empty() {
                let db_conn = db.as_mut().ok_or_else(|| {
                    anyhow!("BioVault database not available for store operations")
                })?;
                for (store_name, store_spec) in &step.store {
                    match store_spec {
                        FlowStoreSpec::Sql(sql) => {
                            if let Some(url) = parse_sql_destination(sql.target.as_deref())? {
                                return Err(anyhow!(
                                    "SQL store '{}' destination '{}' not supported yet (only built-in biovault is available)",
                                    store_name,
                                    url
                                )
                                .into());
                            }
                            let source_path = outputs.get(&sql.source).ok_or_else(|| {
                                anyhow!(
                                    "Store '{}' in step '{}' references unknown output '{}'",
                                    store_name,
                                    step.id,
                                    sql.source
                                )
                            })?;
                            store_sql_output(
                                db_conn,
                                store_name,
                                sql,
                                Path::new(source_path),
                                &run_id,
                            )?;
                        }
                    }
                }
            }
        }

        if !share_locations_by_target.is_empty() {
            let step_entry = step_outputs.entry(step.id.clone()).or_default();
            for (target, share_outputs) in &share_locations_by_target {
                let outputs = step_entry.entry(target.clone()).or_default();
                for (name, location) in share_outputs {
                    outputs
                        .entry(name.clone())
                        .or_insert_with(|| location.syft_url.clone());
                }
            }
        }
    }

    // Mark flow as completed in progress log
    if let Some(ref mut progress) = flow_progress {
        progress.mark_flow_completed();
        if let Some(ref shared_base) = progress_write_base {
            if let Err(e) = write_state(shared_base, &spec.name, &run_id, progress) {
                if verbose_progress {
                    println!("⚠️  Could not write flow completion progress: {}", e);
                }
            }
        }
    }

    println!("\n✅ Flow run completed successfully");
    println!("Results stored under {}", base_results_dir.display());
    Ok(())
}

#[instrument(skip_all, fields(component = "flow", flow = %file), err)]
pub fn validate(file: &str, diagram: bool) -> Result<()> {
    let flow_path = Path::new(file);
    if !flow_path.exists() {
        return Err(anyhow!("Flow file not found: {}", flow_path.display()).into());
    }

    let spec = FlowSpec::load(flow_path)?;
    let db = BioVaultDb::new().ok();
    let validation = validate_internal(flow_path, &spec, db.as_ref())?;

    print_validation(&validation, diagram);

    if validation.has_errors() {
        Err(anyhow!("Flow validation reported errors").into())
    } else {
        println!("\n✅ Flow is valid");
        Ok(())
    }
}

pub fn inspect(file: &str) -> Result<()> {
    let flow_path = Path::new(file);
    if !flow_path.exists() {
        return Err(anyhow!("Flow file not found: {}", flow_path.display()).into());
    }
    let spec = FlowSpec::load(flow_path)?;
    let db = BioVaultDb::new().ok();
    let validation = validate_internal(flow_path, &spec, db.as_ref())?;
    print_steps(&validation);
    Ok(())
}

struct FlowWizard<'a> {
    flow_path: &'a Path,
    db: Option<&'a BioVaultDb>,
    theme: ColorfulTheme,
}

impl<'a> FlowWizard<'a> {
    fn new(flow_path: &'a Path, db: Option<&'a BioVaultDb>) -> Self {
        Self {
            flow_path,
            db,
            theme: ColorfulTheme::default(),
        }
    }

    fn run(&self) -> Result<FlowSpec> {
        println!("🧪 Flow wizard — let’s assemble your flow");

        let spec_name = Input::with_theme(&self.theme)
            .with_prompt("Flow name")
            .validate_with(|input: &String| {
                if input.trim().is_empty() {
                    Err("Flow name cannot be empty")
                } else {
                    Ok(())
                }
            })
            .interact_text()
            .cli_result()?;

        let mut spec = FlowSpec {
            name: spec_name,
            ..FlowSpec::default()
        };

        loop {
            spec = self.add_step_to(spec)?;
            if !Confirm::with_theme(&self.theme)
                .with_prompt("Add another step?")
                .default(true)
                .interact()
                .cli_result()?
            {
                break;
            }
        }

        Ok(spec)
    }

    fn add_step_to(&self, mut spec: FlowSpec) -> Result<FlowSpec> {
        let resolved_steps = resolve_existing_steps(&spec, self.flow_dir(), self.db)?;
        let available_outputs = collect_available_outputs(&resolved_steps);

        let module_choice = self.prompt_module_choice()?;

        let loaded_module = match &module_choice {
            ModuleChoice::Defer => None,
            _ => Some(load_module_spec(&module_choice)?),
        };

        if let Some(module) = &loaded_module {
            println!("\nConfiguring step for module {}", module.summary().bold());
        } else {
            println!("\nConfiguring step without a module (assign later)");
        }

        let default_id = if let Some(module) = &loaded_module {
            generate_default_step_id(module, &spec)
        } else {
            format!("step{}", spec.steps.len() + 1)
        };
        let step_id = Input::with_theme(&self.theme)
            .with_prompt("Step id")
            .default(default_id)
            .validate_with(|value: &String| {
                if value.trim().is_empty() {
                    return Err("Step id cannot be empty");
                }
                Ok(())
            })
            .interact_text()
            .cli_result()?;

        if spec.steps.iter().any(|s| s.id == step_id) {
            return Err(anyhow!("Step id '{}' already exists", step_id).into());
        }

        let mut with_map = BTreeMap::new();
        let mut publish_map = BTreeMap::new();

        if let Some(module) = &loaded_module {
            for input in &module.spec.inputs {
                if let Some(binding) =
                    self.prompt_input_binding(&mut spec, &step_id, input, &available_outputs)?
                {
                    with_map.insert(input.name.clone(), binding);
                }
            }

            for input in &module.spec.inputs {
                with_map.entry(input.name.clone()).or_insert_with(|| {
                    let key = ensure_flow_input(&mut spec, &step_id, input);
                    YamlValue::String(format!("inputs.{}", key))
                });
            }

            publish_map = prompt_publish_map(&self.theme, &module.spec.outputs)?;
        }

        let step = FlowStepSpec {
            id: step_id,
            uses: module_choice.uses_value(),
            where_exec: None,
            foreach: None,
            runs_on: None,
            order: None,
            coordination: None,
            with: with_map,
            publish: publish_map,
            share: BTreeMap::new(),
            store: BTreeMap::new(),
            permissions: Vec::new(),
            barrier: None,
        };

        spec.steps.push(step);
        Ok(spec)
    }

    fn prompt_module_choice(&self) -> Result<ModuleChoice> {
        let mut options = list_registered_modules(self.db);
        options.push(ModuleChoice::EnterPath);
        options.push(ModuleChoice::Defer);

        let items: Vec<String> = options.iter().map(|choice| choice.display()).collect();

        let index = Select::with_theme(&self.theme)
            .with_prompt("Select a module to use")
            .items(&items)
            .default(0)
            .interact()
            .cli_result()?;

        match &options[index] {
            ModuleChoice::EnterPath => {
                let input: String = Input::with_theme(&self.theme)
                    .with_prompt("Module directory path or module.yaml")
                    .interact_text()
                    .cli_result()?;
                let (reference, root) = normalize_module_reference(&input, self.flow_dir())?;
                Ok(ModuleChoice::Path { reference, root })
            }
            ModuleChoice::Defer => Ok(ModuleChoice::Defer),
            other => Ok(other.clone()),
        }
    }

    fn prompt_input_binding(
        &self,
        spec: &mut FlowSpec,
        step_id: &str,
        input: &InputSpec,
        available_outputs: &[CandidateOutput],
    ) -> Result<Option<YamlValue>> {
        println!("\nInput: {} ({})", input.name.bold(), input.raw_type);

        let mut options: Vec<String> = available_outputs
            .iter()
            .map(|candidate| candidate.display_for(input))
            .collect();

        let literal_index = options.len();
        options.push("Enter literal / path".to_string());
        let flow_index = options.len();
        options.push("Use flow input".to_string());

        let selection = Select::with_theme(&self.theme)
            .with_prompt("Binding")
            .items(&options)
            .default(0)
            .interact()
            .cli_result()?;

        if selection < available_outputs.len() {
            let chosen = &available_outputs[selection];
            println!("  → Using {}", chosen.binding.clone().dimmed());
            return Ok(Some(YamlValue::String(chosen.binding.clone())));
        }

        if selection == literal_index {
            let value = Input::with_theme(&self.theme)
                .with_prompt("Enter literal value")
                .interact_text()
                .cli_result()?;
            return Ok(Some(YamlValue::String(value)));
        }

        if selection == flow_index {
            let matching_inputs: Vec<String> = spec
                .inputs
                .iter()
                .filter_map(|(name, existing)| {
                    if types_compatible(existing.raw_type(), &input.raw_type) {
                        Some(name.clone())
                    } else {
                        None
                    }
                })
                .collect();

            let chosen_key = if matching_inputs.is_empty() {
                ensure_flow_input(spec, step_id, input)
            } else {
                let mut flow_options = matching_inputs
                    .iter()
                    .map(|name| format!("Existing: {}", name))
                    .collect::<Vec<_>>();
                flow_options.push(format!("Create new input for '{}'", input.name));

                let selected = Select::with_theme(&self.theme)
                    .with_prompt("Flow input")
                    .items(&flow_options)
                    .default(0)
                    .interact()
                    .cli_result()?;

                if selected < matching_inputs.len() {
                    matching_inputs[selected].clone()
                } else {
                    ensure_flow_input(spec, step_id, input)
                }
            };

            return Ok(Some(YamlValue::String(format!("inputs.{}", chosen_key))));
        }

        Ok(None)
    }

    fn flow_dir(&self) -> &Path {
        self.flow_path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."))
    }
}

fn prompt_publish_map(
    theme: &ColorfulTheme,
    outputs: &[OutputSpec],
) -> Result<BTreeMap<String, String>> {
    if outputs.is_empty() {
        return Ok(BTreeMap::new());
    }

    if !Confirm::with_theme(theme)
        .with_prompt("Customize published outputs (rename/subset)?")
        .default(false)
        .interact()
        .cli_result()?
    {
        return Ok(default_publish_map(outputs));
    }

    let mut map = BTreeMap::new();
    for output in outputs {
        let include = Confirm::with_theme(theme)
            .with_prompt(format!("Publish output '{}' ?", output.name))
            .default(true)
            .interact()
            .cli_result()?;

        if !include {
            continue;
        }

        let alias = Input::with_theme(theme)
            .with_prompt(format!("Alias for '{}'", output.name))
            .default(output.name.clone())
            .interact_text()
            .cli_result()?;

        map.insert(alias, publish_literal(output));
    }

    Ok(map)
}

#[derive(Clone)]
enum ModuleChoice {
    Registered { name: String, path: PathBuf },
    Path { reference: String, root: PathBuf },
    EnterPath,
    Defer,
}

impl ModuleChoice {
    fn display(&self) -> String {
        match self {
            ModuleChoice::Registered { name, path } => {
                format!("{} (registered at {})", name, path.display())
            }
            ModuleChoice::Path { reference, .. } => reference.clone(),
            ModuleChoice::EnterPath => "Enter module path".to_string(),
            ModuleChoice::Defer => "Skip (assign at run time)".to_string(),
        }
    }

    fn uses_value(&self) -> Option<String> {
        match self {
            ModuleChoice::Registered { name, .. } => Some(name.clone()),
            ModuleChoice::Path { reference, .. } => Some(reference.clone()),
            ModuleChoice::EnterPath => None,
            ModuleChoice::Defer => None,
        }
    }
}

struct LoadedModule {
    spec: ModuleSpec,
    root: PathBuf,
}

impl LoadedModule {
    fn summary(&self) -> String {
        format!("{} ({})", self.spec.name, self.root.display())
    }
}

fn load_module_spec(choice: &ModuleChoice) -> Result<LoadedModule> {
    match choice {
        ModuleChoice::Registered { name: _, path } => {
            let spec_path = Path::new(path).join("module.yaml");
            let spec = ModuleSpec::load(&spec_path)?;
            Ok(LoadedModule {
                spec,
                root: path.clone(),
            })
        }
        ModuleChoice::Path { reference: _, root } => {
            let spec_path = Path::new(root).join("module.yaml");
            if !spec_path.exists() {
                return Err(anyhow!("No module.yaml found at {}", spec_path.display()).into());
            }
            let spec = ModuleSpec::load(&spec_path)?;
            Ok(LoadedModule {
                spec,
                root: root.clone(),
            })
        }
        ModuleChoice::EnterPath => Err(anyhow!("Unexpected module choice").into()),
        ModuleChoice::Defer => Err(anyhow!("Deferred module has no spec").into()),
    }
}

fn list_registered_modules(db: Option<&BioVaultDb>) -> Vec<ModuleChoice> {
    if let Some(db) = db {
        match db.list_modules() {
            Ok(modules) => modules
                .into_iter()
                .map(|p| ModuleChoice::Registered {
                    name: p.name,
                    path: PathBuf::from(p.module_path),
                })
                .collect(),
            Err(_) => vec![],
        }
    } else {
        vec![]
    }
}

fn generate_default_step_id(module: &LoadedModule, spec: &FlowSpec) -> String {
    let mut base: String = module
        .spec
        .name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '-'
            }
        })
        .collect();

    while base.contains("--") {
        base = base.replace("--", "-");
    }

    let base = base.trim_matches('-').to_string();

    if base.is_empty() {
        return format!("step{}", spec.steps.len() + 1);
    }

    let mut candidate = base.clone();
    let mut suffix = 2;
    while spec.steps.iter().any(|s| s.id == candidate) {
        candidate = format!("{}-{}", base, suffix);
        suffix += 1;
    }
    candidate
}

struct CandidateOutput {
    binding: String,
    output: OutputSpec,
}

impl CandidateOutput {
    fn display_for(&self, input: &InputSpec) -> String {
        let status = if types_compatible(&input.raw_type, &self.output.raw_type) {
            "match"
        } else {
            "mismatch"
        };
        format!(
            "{} ← {} ({} → {}) [{}]",
            input.name, self.binding, self.output.raw_type, input.raw_type, status
        )
    }
}

fn collect_available_outputs(resolved_steps: &[ResolvedStep]) -> Vec<CandidateOutput> {
    let mut outputs = Vec::new();
    for step in resolved_steps {
        for output in &step.module.spec.outputs {
            outputs.push(CandidateOutput {
                binding: format!("step.{}.outputs.{}", step.step_id, output.name),
                output: output.clone(),
            });
        }
    }
    outputs
}

fn ensure_flow_input(spec: &mut FlowSpec, step_id: &str, input: &InputSpec) -> String {
    if let Some(existing) = spec.inputs.get(&input.name) {
        if existing.raw_type() == input.raw_type {
            return input.name.clone();
        }
    }

    let mut key = input.name.clone();
    let mut suffix = 2;
    while let Some(existing) = spec.inputs.get(&key) {
        if existing.raw_type() == input.raw_type {
            return key;
        }
        key = format!("{}-{}", input.name, suffix);
        suffix += 1;
    }

    // If we had to rename due to type conflict, inform user via console.
    if key != input.name {
        println!(
            "⚠️  Input '{}' on step '{}' conflicts with existing flow input; recorded as '{}'.",
            input.name, step_id, key
        );
    }

    spec.inputs
        .insert(key.clone(), FlowInputSpec::from_type(&input.raw_type));
    key
}

fn is_pending_binding_literal(value: &str) -> bool {
    parse_typed_literal(value)
        .map(|(_, inner)| inner.is_empty())
        .unwrap_or(false)
}

fn is_type_placeholder(binding: &str, expected: &str) -> bool {
    let trimmed = binding.trim();
    if trimmed.contains('(') || trimmed.contains(')') {
        return false;
    }
    crate::module_spec::types_compatible(trimmed, expected)
}

fn normalize_manifest_binding(value: &str) -> (String, bool) {
    // .url_list is the preferred suffix - indicates you're getting a list of URLs
    // .manifest kept for backwards compatibility
    for suffix in [".url_list", ".manifest", ".all", ".paths"] {
        if let Some(base) = value.strip_suffix(suffix) {
            return (base.to_string(), true);
        }
    }
    (value.to_string(), false)
}

fn parse_typed_literal(value: &str) -> Option<(String, String)> {
    let trimmed = value.trim();
    let open = trimmed.find('(')?;
    let close = trimmed.rfind(')')?;
    if close < open {
        return None;
    }
    let type_part = trimmed[..open].trim();
    if type_part.is_empty() {
        return None;
    }
    if !type_part
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '[' || c == ']' || c == '?')
    {
        return None;
    }
    let value_part = trimmed[open + 1..close].trim().to_string();
    Some((type_part.to_string(), value_part))
}

fn publish_literal(output: &OutputSpec) -> String {
    match output.path.as_deref() {
        Some(path) if !path.trim().is_empty() => {
            format!("{}({})", output.raw_type, path.trim())
        }
        _ => output.raw_type.clone(),
    }
}

fn default_publish_map(outputs: &[OutputSpec]) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    for output in outputs {
        map.insert(output.name.clone(), publish_literal(output));
    }
    map
}

fn parse_overrides(
    extra_args: &[String],
    initial_results_dir: Option<String>,
) -> Result<ParseOverridesResult> {
    let mut step_overrides = HashMap::new();
    let mut flow_overrides = HashMap::new();
    let mut results_dir = initial_results_dir;
    let mut nextflow_args = Vec::new();
    let mut i = 0;
    while i < extra_args.len() {
        let arg = &extra_args[i];
        if arg == "--set" {
            if i + 1 >= extra_args.len() {
                return Err(anyhow!("--set requires an argument like step.input=value").into());
            }
            parse_override_pair(&extra_args[i + 1], &mut step_overrides, &mut flow_overrides)?;
            i += 2;
        } else if let Some(kv) = arg.strip_prefix("--set=") {
            parse_override_pair(kv, &mut step_overrides, &mut flow_overrides)?;
            i += 1;
        } else if arg == "--results-dir" {
            if i + 1 >= extra_args.len() {
                return Err(anyhow!("--results-dir requires a value").into());
            }
            results_dir = Some(extra_args[i + 1].clone());
            i += 2;
        } else if let Some(dir) = arg.strip_prefix("--results-dir=") {
            results_dir = Some(dir.to_string());
            i += 1;
        } else if arg == "--" {
            nextflow_args.extend(extra_args[i + 1..].iter().cloned());
            break;
        } else {
            nextflow_args.push(arg.clone());
            i += 1;
        }
    }
    Ok((step_overrides, flow_overrides, results_dir, nextflow_args))
}

fn parse_override_pair(
    pair: &str,
    step_overrides: &mut HashMap<(String, String), String>,
    flow_overrides: &mut HashMap<String, String>,
) -> Result<()> {
    let (key, value) = pair
        .split_once('=')
        .ok_or_else(|| anyhow!("Override must be in the form step.input=value"))?;
    let (step, input) = key
        .split_once('.')
        .ok_or_else(|| anyhow!("Override key must be step.input"))?;
    if step.trim().is_empty() || input.trim().is_empty() {
        return Err(anyhow!("Override key must include both step and input names").into());
    }
    let value = value.trim().to_string();
    if step.trim() == "inputs" {
        flow_overrides.insert(input.trim().to_string(), value);
    } else {
        step_overrides.insert((step.trim().to_string(), input.trim().to_string()), value);
    }
    Ok(())
}

fn resolve_syft_path(
    value: &str,
    syftbox_data_dir: Option<&PathBuf>,
    sandbox_root: Option<&PathBuf>,
) -> Result<String> {
    if !value.starts_with("syft://") {
        return Ok(value.to_string());
    }
    let parsed = SyftURL::parse(value)?;

    if let Some(root) = sandbox_root {
        let path = root
            .join(&parsed.email)
            .join("datasites")
            .join(&parsed.email)
            .join(&parsed.path);
        return Ok(path.to_string_lossy().to_string());
    }

    let data_dir = syftbox_data_dir
        .ok_or_else(|| anyhow!("SYFTBOX_DATA_DIR is required to resolve syft:// outputs"))?;
    let path = data_dir
        .join("datasites")
        .join(&parsed.email)
        .join(&parsed.path);
    Ok(path.to_string_lossy().to_string())
}

fn resolve_binding(
    binding: &str,
    expected_type: &str,
    flow_inputs: &BTreeMap<String, FlowInputSpec>,
    resolved_inputs: &HashMap<String, String>,
    step_outputs: &HashMap<String, HashMap<String, HashMap<String, String>>>,
    current_step_id: &str,
    current_datasite: Option<&str>,
    results_dir: &Path,
    step_output_order: &HashMap<String, Vec<String>>,
    syftbox_data_dir: Option<&PathBuf>,
    sandbox_root: Option<&PathBuf>,
) -> Result<String> {
    if let Some(input_name) = binding.strip_prefix("inputs.") {
        if let Some(value) = resolved_inputs.get(input_name) {
            return Ok(value.clone());
        }
        if let Some(_spec) = flow_inputs.get(input_name) {
            return Err(anyhow!(
                "Flow input '{}' is not set. Provide a value with --set inputs.{}=<value>.",
                input_name,
                input_name
            )
            .into());
        } else {
            return Err(anyhow!(
                "Flow input '{}' referenced in step '{}' is not declared",
                input_name,
                current_step_id
            )
            .into());
        }
    }

    if binding.starts_with("step.") {
        let parts: Vec<&str> = binding.split('.').collect();
        let wants_manifest = parts.len() == 5 && matches!(parts[4], "url_list" | "manifest" | "all" | "paths");
        let valid_namespace = parts.len() >= 3 && matches!(parts[2], "outputs" | "share");
        if (parts.len() != 4 && !wants_manifest) || !valid_namespace {
            return Err(anyhow!(
                "Invalid step reference '{}' in step '{}'. Expected format step.<id>.outputs.<name> or step.<id>.share.<name>.url_list",
                binding,
                current_step_id
            )
            .into());
        }
        let source_step = parts[1];
        let output_name = parts[3];
        let outputs_by_site = step_outputs.get(source_step).ok_or_else(|| {
            anyhow!(
                "Step '{}' references outputs from step '{}' which has not run yet",
                current_step_id,
                source_step
            )
        })?;

        if !wants_manifest {
            if let Some(site) = current_datasite {
                if let Some(outputs) = outputs_by_site.get(site) {
                    if let Some(value) = outputs.get(output_name) {
                        return resolve_syft_path(value, syftbox_data_dir, sandbox_root);
                    }
                }
            }

            if outputs_by_site.len() == 1 {
                if let Some(outputs) = outputs_by_site.values().next() {
                    if let Some(value) = outputs.get(output_name) {
                        return resolve_syft_path(value, syftbox_data_dir, sandbox_root);
                    }
                }
            }
        }

        if outputs_by_site.is_empty() {
            return Err(anyhow!(
                "Step '{}' output '{}' not found when referenced from step '{}'",
                source_step,
                output_name,
                current_step_id
            )
            .into());
        }

        let manifest_dir = results_dir.join("manifests").join(source_step);
        let manifest_path = manifest_dir.join(format!("{}_paths.txt", output_name));
        std::fs::create_dir_all(&manifest_dir).with_context(|| {
            format!(
                "Failed to create manifest directory {}",
                manifest_dir.display()
            )
        })?;

        let mut ordered_sites = step_output_order
            .get(source_step)
            .cloned()
            .unwrap_or_default();
        if ordered_sites.is_empty() {
            ordered_sites = outputs_by_site.keys().cloned().collect();
            ordered_sites.sort();
        }

        let mut lines = Vec::new();
        for site in ordered_sites {
            if let Some(outputs) = outputs_by_site.get(&site) {
                if let Some(value) = outputs.get(output_name) {
                    let resolved = resolve_syft_path(value, syftbox_data_dir, sandbox_root)?;
                    lines.push(format!("{}\t{}", site, resolved));
                }
            }
        }

        if lines.is_empty() {
            return Err(anyhow!(
                "Step '{}' output '{}' not found when referenced from step '{}'",
                source_step,
                output_name,
                current_step_id
            )
            .into());
        }

        std::fs::write(&manifest_path, lines.join("\n")).with_context(|| {
            format!(
                "Failed to write manifest for step '{}' output '{}' to {}",
                source_step,
                output_name,
                manifest_path.display()
            )
        })?;
        return Ok(manifest_path.to_string_lossy().to_string());
    }

    if is_type_placeholder(binding, expected_type) {
        return Err(anyhow!(
            "Input '{}' is still a placeholder. Provide a value like {}(<value>).",
            expected_type,
            expected_type
        )
        .into());
    }

    literal_to_value(binding, expected_type)
}

fn literal_to_value(literal: &str, expected_type: &str) -> Result<String> {
    if let Some((literal_type, literal_value)) = parse_typed_literal(literal) {
        if literal_value.is_empty() {
            return Err(anyhow!(
                "Value for {} cannot be empty. Use {}(<value>).",
                expected_type,
                literal_type
            )
            .into());
        }

        if !types_compatible(expected_type, &literal_type) {
            return Err(anyhow!(
                "Type mismatch. Expected {}, but literal binding declares {}",
                expected_type,
                literal_type
            )
            .into());
        }
        Ok(literal_value)
    } else {
        Ok(literal.to_string())
    }
}

fn sanitize_identifier(name: &str) -> String {
    let mut result = String::new();
    for c in name.chars() {
        if c.is_ascii_alphanumeric() {
            result.push(c.to_ascii_lowercase());
        } else {
            result.push('_');
        }
    }
    while result.contains("__") {
        result = result.replace("__", "_");
    }
    let mut trimmed = result.trim_matches('_').to_string();
    if trimmed.is_empty() {
        trimmed = DEFAULT_TABLE_FALLBACK.to_string();
    }
    if trimmed.chars().next().unwrap().is_ascii_digit() {
        trimmed = format!("{}_{}", DEFAULT_TABLE_FALLBACK, trimmed);
    }
    trimmed
}

fn resolve_step_targets(step: &FlowStepSpec) -> Vec<String> {
    if let Some(runs_on) = &step.runs_on {
        return runs_on.as_vec();
    }
    step.foreach.clone().unwrap_or_default()
}

fn resolve_run_targets(
    step_targets: &[String],
    current_datasite: Option<&str>,
    run_all: bool,
) -> Vec<Option<String>> {
    if step_targets.is_empty() {
        return vec![None];
    }
    if run_all || current_datasite.is_none() {
        return step_targets.iter().cloned().map(Some).collect();
    }
    let current = current_datasite.unwrap();
    if step_targets.iter().any(|site| site == current) {
        vec![Some(current.to_string())]
    } else {
        Vec::new()
    }
}

fn resolve_all_datasites_from_spec(
    spec: &FlowSpec,
    resolved_inputs: &HashMap<String, String>,
) -> Vec<String> {
    if let Some(value) = resolved_inputs.get("datasites") {
        return value.split(',').map(|s| s.trim().to_string()).collect();
    }
    if !spec.datasites.is_empty() {
        return spec.datasites.clone();
    }
    Vec::new()
}

fn expand_permission_patterns(
    patterns: &[String],
    all_datasites: &[String],
    current_datasite: Option<&str>,
) -> Vec<String> {
    let mut result = Vec::new();
    for pattern in patterns {
        let mut trimmed = pattern.trim().to_string();

        // Expand {datasite.current} first
        if let Some(current) = current_datasite {
            trimmed = trimmed.replace("{datasite.current}", current);
        }

        if trimmed.contains("{datasites[") {
            if trimmed.contains("[*]") {
                result.extend(all_datasites.iter().cloned());
            } else if let Some(inner) = trimmed
                .strip_prefix("{datasites[")
                .and_then(|s| s.strip_suffix("]}"))
            {
                if let Ok(index) = inner.trim().parse::<usize>() {
                    if let Some(site) = all_datasites.get(index) {
                        result.push(site.clone());
                    }
                }
            }
        } else {
            result.push(trimmed.to_string());
        }
    }
    result
}

fn resolve_current_datasite() -> Option<String> {
    if let Ok(env_override) = env::var("BIOVAULT_DATASITE_OVERRIDE") {
        if !env_override.trim().is_empty() {
            return Some(env_override.trim().to_string());
        }
    }
    if let Ok(cfg) = crate::config::Config::load() {
        if !cfg.email.trim().is_empty() {
            return Some(cfg.email);
        }
    }
    if let Ok(env_email) = env::var("SYFTBOX_EMAIL") {
        if !env_email.trim().is_empty() {
            return Some(env_email.trim().to_string());
        }
    }
    if let Ok(env_email) = env::var("BIOVAULT_DATASITE") {
        if !env_email.trim().is_empty() {
            return Some(env_email.trim().to_string());
        }
    }
    None
}

fn render_flow_template(
    template: &str,
    current: &str,
    datasites: &[String],
    run_id: &str,
    flow_name: &str,
    step_number: Option<usize>,
    step_id: Option<&str>,
    user_vars: &std::collections::BTreeMap<String, String>,
) -> String {
    let mut rendered = template.to_string();

    // First pass: expand user-defined variables (namespaced as {vars.name})
    // Do multiple passes to handle variables that reference other variables
    for _ in 0..5 {
        let before = rendered.clone();
        for (name, value) in user_vars {
            rendered = rendered.replace(&format!("{{vars.{}}}", name), value);
        }
        if rendered == before {
            break;
        }
    }

    // Second pass: expand built-in variables
    rendered = rendered.replace("{current_datasite}", current);
    rendered = rendered.replace("{datasite.current}", current);
    rendered = rendered.replace("{flow_name}", flow_name);
    rendered = rendered.replace("{run_id}", run_id);
    if let Some(index) = datasites.iter().position(|site| site == current) {
        rendered = rendered.replace("{datasites.index}", &index.to_string());
        rendered = rendered.replace("{datasite.index}", &index.to_string());
    }
    rendered = rendered.replace("{datasites}", &datasites.join(","));
    if let Some(num) = step_number {
        rendered = rendered.replace("{step.number}", &format!("{:02}", num));
    }
    if let Some(id) = step_id {
        rendered = rendered.replace("{step.id}", id);
    }
    rendered
}

struct ShareLocation {
    local_path: PathBuf,
    syft_url: String,
}

fn resolve_share_location(
    url_template: &str,
    current: &str,
    datasites: &[String],
    run_id: &str,
    flow_name: &str,
    data_dir: &Path,
    step_number: usize,
    step_id: &str,
    user_vars: &std::collections::BTreeMap<String, String>,
) -> Result<ShareLocation> {
    // Expand template variables in the URL
    let rendered = render_flow_template(
        url_template,
        current,
        datasites,
        run_id,
        flow_name,
        Some(step_number),
        Some(step_id),
        user_vars,
    );

    // Require syft:// URL scheme
    if !rendered.starts_with("syft://") {
        return Err(anyhow!(
            "Share URL must use syft:// scheme, got: {}",
            rendered
        ).into());
    }

    // Parse and validate using SDK
    let parsed = SyftURL::parse(&rendered)
        .with_context(|| format!("Invalid share URL: {}", rendered))?;

    // Validate no path traversal
    validate_no_path_traversal(&parsed.path)?;

    let local_path = data_dir
        .join("datasites")
        .join(&parsed.email)
        .join(&parsed.path);

    Ok(ShareLocation {
        local_path,
        syft_url: rendered,
    })
}

fn dedupe_emails(values: &[String]) -> Vec<String> {
    let mut set = HashSet::new();
    let mut out = Vec::new();
    for value in values {
        let trimmed = value.trim();
        if trimmed.is_empty() || !set.insert(trimmed.to_string()) {
            continue;
        }
        out.push(trimmed.to_string());
    }
    out
}

/// Validate a local path has no path traversal
fn validate_no_path_traversal(path: &str) -> Result<()> {
    if path.contains("..") {
        return Err(anyhow!("Path traversal (..) not allowed: {}", path).into());
    }
    Ok(())
}

/// Parse share source reference (e.g., "self.outputs.rsids" -> "rsids")
fn parse_share_source(source: &str) -> Result<String> {
    let source = source.trim();
    if let Some(name) = source.strip_prefix("self.outputs.") {
        if name.is_empty() {
            return Err(anyhow!("Invalid share source: {}", source).into());
        }
        return Ok(name.to_string());
    }
    if let Some(name) = source.strip_prefix("outputs.") {
        if name.is_empty() {
            return Err(anyhow!("Invalid share source: {}", source).into());
        }
        return Ok(name.to_string());
    }
    // For backwards compatibility, allow bare names
    Ok(source.to_string())
}

fn execute_permissions_step(
    permissions: &[FlowPermissionSpec],
    target_data_dir: &Path,
    all_datasites: &[String],
    current_datasite: Option<&str>,
    flow_name: &str,
    run_id: &str,
) -> Result<()> {
    use crate::types::{AccessControl, PermissionRule, SyftPermissions};

    for perm_spec in permissions {
        // Expand templates in the URL
        let expanded_url = perm_spec
            .url
            .replace("{flow_name}", flow_name)
            .replace("{run_id}", run_id);
        let expanded_url = if let Some(current) = current_datasite {
            expanded_url.replace("{datasite.current}", current)
        } else {
            expanded_url
        };

        // Parse and validate syft:// URL using the SDK parser
        let parsed = SyftURL::parse(&expanded_url)
            .with_context(|| format!("Invalid permissions URL: {}", expanded_url))?;

        // Validate no path traversal
        validate_no_path_traversal(&parsed.path)?;

        // Build local path from syft URL components
        let full_path = target_data_dir.join("datasites").join(&parsed.email).join(&parsed.path);

        std::fs::create_dir_all(&full_path)
            .with_context(|| format!("Failed to create permissions directory: {}", full_path.display()))?;

        let rules: Vec<PermissionRule> = perm_spec
            .rules
            .iter()
            .map(|rule| {
                let read = expand_permission_patterns(&rule.access.read, all_datasites, current_datasite);
                let write = expand_permission_patterns(&rule.access.write, all_datasites, current_datasite);
                let admin = expand_permission_patterns(&rule.access.admin, all_datasites, current_datasite);
                PermissionRule {
                    pattern: rule.pattern.clone(),
                    access: AccessControl {
                        read: dedupe_emails(&read),
                        write: dedupe_emails(&write),
                        admin: dedupe_emails(&admin),
                    },
                }
            })
            .collect();

        let permissions = SyftPermissions { rules };
        let perms_path = full_path.join("syft.pub.yaml");
        let yaml = serde_yaml::to_string(&permissions)
            .context("Failed to serialize permissions")?;
        std::fs::write(&perms_path, yaml)
            .with_context(|| format!("Failed to write {}", perms_path.display()))?;

        println!("  📝 Created permissions at: {}", perms_path.display());
    }

    Ok(())
}

fn create_shared_file(
    source_path: &Path,
    shared_path: &Path,
    share_spec: &FlowShareSpec,
    all_datasites: &[String],
    current_datasite: Option<&str>,
    _sandbox_root: Option<&PathBuf>,
) -> Result<()> {
    let read = if share_spec.read.is_empty() {
        all_datasites.to_vec()
    } else {
        expand_permission_patterns(&share_spec.read, all_datasites, current_datasite)
    };
    let read = dedupe_emails(&read);
    let write = if share_spec.write.is_empty() {
        read.clone()
    } else {
        expand_permission_patterns(&share_spec.write, all_datasites, current_datasite)
    };
    let write = dedupe_emails(&write);
    let admin = expand_permission_patterns(&share_spec.admin, all_datasites, current_datasite);
    let admin = dedupe_emails(&admin);

    if let Some(parent) = shared_path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create shared directory {}", parent.display()))?;
        let file_pattern = shared_path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("*");
        let permissions = SyftPermissions {
            rules: vec![PermissionRule {
                pattern: file_pattern.to_string(),
                access: AccessControl {
                    read: read.clone(),
                    write,
                    admin,
                },
            }],
        };
        let perms_path = parent.join("syft.pub.yaml");
        let yaml =
            serde_yaml::to_string(&permissions).context("Failed to serialize share permissions")?;
        std::fs::write(&perms_path, yaml)
            .with_context(|| format!("Failed to write {}", perms_path.display()))?;
    }

    std::fs::copy(source_path, shared_path).with_context(|| {
        format!(
            "Failed to copy {} -> {}",
            source_path.display(),
            shared_path.display()
        )
    })?;

    Ok(())
}

fn detect_table_name(store_name: &str, spec: &FlowSqlStoreSpec, run_id: &str) -> String {
    let template = spec
        .table
        .clone()
        .unwrap_or_else(|| format!("{}_{}", store_name, run_id));
    let table_name = template.replace("{run_id}", run_id);
    format!("{}{}", RESULTS_TABLE_PREFIX, table_name)
}

fn detect_format(spec: &FlowSqlStoreSpec, output_path: &Path) -> String {
    if let Some(fmt) = spec.format.as_deref() {
        return fmt.to_ascii_lowercase();
    }
    output_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_else(|| "csv".to_string())
}

fn get_delimiter_for_format(format: &str) -> u8 {
    match format {
        "tsv" => b'\t',
        _ => b',',
    }
}

fn parse_sql_destination(target: Option<&str>) -> Result<Option<String>> {
    let Some(raw) = target.map(|t| t.trim()) else {
        return Ok(None);
    };
    let raw = raw.trim_matches(|c| c == '\'' || c == '"');
    if raw.is_empty() {
        return Ok(None);
    }
    let upper = raw.to_ascii_uppercase();
    if upper == "BIOVAULT" || upper == "SQL" || upper == "SQL()" {
        return Ok(None);
    }
    if upper.starts_with("SQL") {
        let inner = raw
            .strip_prefix("SQL")
            .and_then(|s| s.strip_prefix('('))
            .map(|s| s.trim_end_matches(')'))
            .unwrap_or("")
            .trim();
        if inner.is_empty() {
            return Ok(None);
        }
        if let Some(url_spec) = inner
            .strip_prefix("url:")
            .or_else(|| inner.strip_prefix("url="))
        {
            let url = url_spec.trim().trim_matches(|c| c == '"' || c == '\'');
            if url.is_empty() {
                return Err(anyhow!("SQL destination url cannot be empty").into());
            }
            return Ok(Some(url.to_string()));
        }
        return Err(anyhow!(
            "Unsupported SQL destination syntax '{}'. Use SQL() or SQL(url:<connection>)",
            raw
        )
        .into());
    }
    Err(anyhow!(
        "Unsupported SQL destination '{}'. Use SQL() or SQL(url:<connection>)",
        raw
    )
    .into())
}

fn store_sql_output(
    db: &mut BioVaultDb,
    store_name: &str,
    spec: &FlowSqlStoreSpec,
    output_path: &Path,
    run_id: &str,
) -> Result<()> {
    if !output_path.exists() {
        return Err(anyhow!(
            "Store '{}' refers to missing output file {}",
            store_name,
            output_path.display()
        )
        .into());
    }

    if let Some(url) = parse_sql_destination(spec.target.as_deref())? {
        return Err(anyhow!(
            "SQL store '{}' destination '{}' not supported yet (only built-in biovault is available)",
            store_name,
            url
        )
        .into());
    }

    let table_template = detect_table_name(store_name, spec, run_id);
    let table_identifier = sanitize_identifier(&table_template);
    if table_identifier.is_empty() {
        return Err(anyhow!(
            "Unable to derive a valid table name from '{}'",
            table_template
        )
        .into());
    }

    let format = detect_format(spec, output_path);
    let delimiter = get_delimiter_for_format(&format);

    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(delimiter)
        .from_path(output_path)
        .with_context(|| format!("Failed to open {}", output_path.display()))?;

    let headers = reader
        .headers()
        .with_context(|| format!("Failed to read headers from {}", output_path.display()))?
        .clone();
    if headers.is_empty() {
        return Err(anyhow!(
            "CSV {} has no headers; cannot create SQL table",
            output_path.display()
        )
        .into());
    }

    let participant_column_raw = spec.participant_column.as_deref();
    let participant_column = participant_column_raw.map(|s| s.to_ascii_lowercase());
    let mut participant_column_found = participant_column.is_none();

    let mut used_columns = HashSet::new();
    let mut columns = Vec::new();
    for (idx, header) in headers.iter().enumerate() {
        let mut base = sanitize_identifier(header);
        if let Some(ref target) = participant_column {
            if header.eq_ignore_ascii_case(target) {
                base = "participant_id".to_string();
                participant_column_found = true;
            }
        }
        if base.is_empty() {
            base = format!("{}{}", DEFAULT_COLUMN_PREFIX, idx + 1);
        }
        let mut candidate = base.clone();
        let mut counter = 2;
        while !used_columns.insert(candidate.clone()) {
            candidate = format!("{}_{}", base, counter);
            counter += 1;
        }
        columns.push((candidate, header.to_string(), idx));
    }

    if let Some(raw) = participant_column_raw {
        if !participant_column_found {
            return Err(anyhow!(
                "participant_column '{}' not found in output '{}'",
                raw,
                output_path.display()
            )
            .into());
        }
    }

    let mut create_defs = Vec::new();
    for (sanitized, _, _) in &columns {
        let col_def = format!("\"{}\" TEXT", sanitized);
        create_defs.push(col_def);
    }

    let total_columns = columns.len();
    let tx = db
        .conn
        .transaction()
        .with_context(|| format!("Failed to start SQL transaction for store '{}'", store_name))?;
    if spec.overwrite.unwrap_or(true) {
        tx.execute(
            &format!("DROP TABLE IF EXISTS \"{}\"", table_identifier),
            [],
        )
        .with_context(|| format!("Failed to drop existing table {}", table_identifier))?;
    }
    tx.execute(
        &format!(
            "CREATE TABLE IF NOT EXISTS \"{}\" ({})",
            table_identifier,
            create_defs.join(", ")
        ),
        [],
    )
    .with_context(|| format!("Failed to create table {}", table_identifier))?;

    let column_list = columns
        .iter()
        .map(|(sanitized, _, _)| format!("\"{}\"", sanitized))
        .collect::<Vec<_>>()
        .join(", ");
    let placeholders = (0..total_columns)
        .map(|_| "?")
        .collect::<Vec<_>>()
        .join(", ");
    let insert_sql = format!(
        "INSERT OR REPLACE INTO \"{}\" ({}) VALUES ({})",
        table_identifier, column_list, placeholders
    );

    let mut stmt = tx
        .prepare(&insert_sql)
        .with_context(|| format!("Failed to prepare insert for table {}", table_identifier))?;
    let mut row_count = 0usize;
    for record in reader.records() {
        let record =
            record.with_context(|| format!("Failed to read row from {}", output_path.display()))?;
        let mut values = Vec::with_capacity(total_columns);
        for (_, _, idx) in &columns {
            values.push(record.get(*idx).unwrap_or("").to_string());
        }
        stmt.execute(params_from_iter(values.iter()))
            .with_context(|| format!("Failed to insert row into table {}", table_identifier))?;
        row_count += 1;
    }
    drop(stmt);
    tx.commit()
        .with_context(|| format!("Failed to commit SQL store for table {}", table_identifier))?;

    let db_path = crate::config::get_biovault_home()?.join("biovault.db");
    println!(
        "💾 Stored '{}' output '{}' into table {} (rows: {}).",
        store_name, spec.source, table_identifier, row_count
    );
    println!("    source: {}", output_path.display());
    println!("    database: {}", db_path.display());

    Ok(())
}

struct ResolvedStep {
    step_id: String,
    module: LoadedModule,
    store: BTreeMap<String, FlowStoreSpec>,
}

fn resolve_existing_steps(
    spec: &FlowSpec,
    flow_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<Vec<ResolvedStep>> {
    let mut resolved = Vec::new();
    for step in &spec.steps {
        if let Some(module) = resolve_module(step, flow_dir, db)? {
            resolved.push(ResolvedStep {
                step_id: step.id.clone(),
                module,
                store: step.store.clone(),
            });
        }
    }
    Ok(resolved)
}

fn normalize_module_reference(raw: &str, base_dir: &Path) -> Result<(String, PathBuf)> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(anyhow!("Module path cannot be empty").into());
    }

    let mut reference_path = PathBuf::from(trimmed);
    if reference_path
        .file_name()
        .map(|name| name == "module.yaml")
        .unwrap_or(false)
    {
        reference_path.pop();
    }

    let reference_string = if reference_path.as_os_str().is_empty() {
        ".".to_string()
    } else {
        reference_path.to_string_lossy().to_string()
    };

    let mut absolute = if reference_path.as_os_str().is_empty() {
        base_dir.to_path_buf()
    } else if reference_path.is_absolute() {
        reference_path.clone()
    } else if base_dir.exists() {
        base_dir.join(&reference_path)
    } else {
        env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(&reference_path)
    };

    if !absolute.join("module.yaml").exists() {
        return Err(anyhow!(
            "No module.yaml found at {}",
            absolute.join("module.yaml").display()
        )
        .into());
    }

    // Normalize absolute path by removing trailing components like '.'
    if let Ok(canon) = absolute.canonicalize() {
        absolute = canon;
    }

    Ok((reference_string, absolute))
}

fn resolve_module(
    step: &FlowStepSpec,
    flow_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<Option<LoadedModule>> {
    match step.uses.as_deref().map(str::trim) {
        None | Some("") => Ok(None),
        Some(uses) => resolve_module_by_uses(uses, flow_dir, db).map(Some),
    }
}

fn resolve_module_by_uses(
    uses: &str,
    flow_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<LoadedModule> {
    // Try 1: Resolve as path (relative or absolute)
    if let Ok((reference, root)) = normalize_module_reference(uses, flow_dir) {
        return load_module_spec(&ModuleChoice::Path { reference, root });
    }

    // Try 2: Database lookup with exact name
    if let Some(db) = db {
        if let Ok(Some(module)) = db.get_module(uses) {
            return load_module_spec(&ModuleChoice::Registered {
                name: module.name,
                path: PathBuf::from(module.module_path),
            });
        }
    }

    Err(anyhow!("Unable to resolve module reference '{}'", uses).into())
}

struct FlowValidationResult {
    resolved_steps: Vec<ResolvedStep>,
    unresolved_steps: Vec<UnresolvedStep>,
    bindings: Vec<BindingStatus>,
    errors: Vec<String>,
    warnings: Vec<String>,
    flow_inputs: BTreeMap<String, FlowInputSpec>,
}

impl FlowValidationResult {
    fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

struct BindingStatus {
    step_id: String,
    input: String,
    binding: Option<String>,
    expected: String,
    actual: Option<String>,
    state: BindingState,
    message: Option<String>,
}

enum BindingState {
    Ok,
    Warning,
    Error,
    Pending,
}

struct UnresolvedStep {
    step_id: String,
    reason: String,
}

// Note: validate_module_interface and validate_module_assets functions removed
// They require FlowFileSpec.modules which is not available in FlowSpec
// TODO: Re-add when modules field is added to FlowSpec

// Note: Module interface validation requires access to FlowFileSpec.modules
// which is not available in FlowSpec. This validation is skipped for now.
// TODO: Add modules field to FlowSpec or refactor validation to use FlowFile

fn validate_internal(
    flow_path: &Path,
    spec: &FlowSpec,
    db: Option<&BioVaultDb>,
) -> Result<FlowValidationResult> {
    spec.ensure_unique_step_ids()?;
    let flow_dir = flow_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    let mut resolved_steps = Vec::new();
    let mut unresolved_steps = Vec::new();
    let mut available: HashMap<String, (String, OutputSpec)> = HashMap::new();
    let mut bindings = Vec::new();
    let mut errors = Vec::new();
    let mut warnings = Vec::new();

    for step_spec in &spec.steps {
        match resolve_module(step_spec, flow_dir, db) {
            Ok(Some(module)) => {
                let resolved_step = ResolvedStep {
                    step_id: step_spec.id.clone(),
                    module,
                    store: step_spec.store.clone(),
                };

                let module_spec = &resolved_step.module.spec;

                // Note: Module interface/asset validation skipped - requires FlowFileSpec.modules

                for input in &module_spec.inputs {
                    let raw_binding = step_spec.with.get(&input.name).cloned();
                    let binding_str = raw_binding.as_ref().and_then(value_to_string);

                    if let Some(binding_value) = &binding_str {
                        let (binding_key, binding_is_manifest) =
                            normalize_manifest_binding(binding_value);
                        if let Some((producer_step, output)) = available.get(&binding_key) {
                            let actual_type = if binding_is_manifest {
                                "File"
                            } else {
                                output.raw_type.as_str()
                            };
                            if !types_compatible(&input.raw_type, actual_type) {
                                let msg = format!(
                                    "Type mismatch: expected {} but found {} from step {}",
                                    input.raw_type, actual_type, producer_step
                                );
                                errors.push(msg.clone());
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: Some(actual_type.to_string()),
                                    state: BindingState::Error,
                                    message: Some(msg),
                                });
                                continue;
                            }

                            bindings.push(BindingStatus {
                                step_id: resolved_step.step_id.clone(),
                                input: input.name.clone(),
                                binding: Some(binding_value.clone()),
                                expected: input.raw_type.clone(),
                                actual: Some(actual_type.to_string()),
                                state: BindingState::Ok,
                                message: None,
                            });
                            continue;
                        }

                        if let Some(input_name) = binding_value.strip_prefix("inputs.") {
                            if let Some(spec_input) = spec.inputs.get(input_name) {
                                if !types_compatible(&input.raw_type, spec_input.raw_type()) {
                                    let msg = format!(
                                        "Type mismatch: expected {} but flow input '{}' declares {}",
                                        input.raw_type,
                                        input_name,
                                        spec_input.raw_type()
                                    );
                                    errors.push(msg.clone());
                                    bindings.push(BindingStatus {
                                        step_id: resolved_step.step_id.clone(),
                                        input: input.name.clone(),
                                        binding: Some(binding_value.clone()),
                                        expected: input.raw_type.clone(),
                                        actual: Some(spec_input.raw_type().to_string()),
                                        state: BindingState::Error,
                                        message: Some(msg),
                                    });
                                } else {
                                    bindings.push(BindingStatus {
                                        step_id: resolved_step.step_id.clone(),
                                        input: input.name.clone(),
                                        binding: Some(binding_value.clone()),
                                        expected: input.raw_type.clone(),
                                        actual: None,
                                        state: BindingState::Pending,
                                        message: Some(format!(
                                            "Depends on flow input '{}'. Set it via --set inputs.{}=<value>.",
                                            input_name,
                                            input_name
                                        )),
                                    });
                                }
                            } else {
                                errors.push(format!(
                                    "Binding '{}' for step '{}' input '{}' references unknown flow input",
                                    binding_value,
                                    resolved_step.step_id,
                                    input.name
                                ));
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: None,
                                    state: BindingState::Error,
                                    message: Some("Unknown flow input".to_string()),
                                });
                            }
                        } else if binding_value.starts_with("step.") {
                            errors.push(format!(
                                "Binding '{}' for step '{}' input '{}' references an unknown output",
                                binding_value, resolved_step.step_id, input.name
                            ));
                            bindings.push(BindingStatus {
                                step_id: resolved_step.step_id.clone(),
                                input: input.name.clone(),
                                binding: Some(binding_value.clone()),
                                expected: input.raw_type.clone(),
                                actual: None,
                                state: BindingState::Error,
                                message: Some("Unknown output reference".to_string()),
                            });
                        } else if is_type_placeholder(binding_value, &input.raw_type) {
                            bindings.push(BindingStatus {
                                step_id: resolved_step.step_id.clone(),
                                input: input.name.clone(),
                                binding: Some(binding_value.clone()),
                                expected: input.raw_type.clone(),
                                actual: None,
                                state: BindingState::Pending,
                                message: Some(
                                    "Placeholder – replace with Type(value) when ready".to_string(),
                                ),
                            });
                        } else if is_pending_binding_literal(binding_value) {
                            bindings.push(BindingStatus {
                                step_id: resolved_step.step_id.clone(),
                                input: input.name.clone(),
                                binding: Some(binding_value.clone()),
                                expected: input.raw_type.clone(),
                                actual: None,
                                state: BindingState::Pending,
                                message: Some(
                                    "Pending input – assign via CLI override or flow edit"
                                        .to_string(),
                                ),
                            });
                        } else if let Some((literal_type, literal_value)) =
                            parse_typed_literal(binding_value)
                        {
                            if literal_value.is_empty() {
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: Some(literal_type),
                                    state: BindingState::Pending,
                                    message: Some(
                                        "Provide a value inside parentheses to supply this input"
                                            .to_string(),
                                    ),
                                });
                            } else if !types_compatible(&input.raw_type, &literal_type) {
                                let msg = format!(
                                    "Type mismatch: expected {} but literal binding declared {}",
                                    input.raw_type, literal_type
                                );
                                errors.push(msg.clone());
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: Some(literal_type),
                                    state: BindingState::Error,
                                    message: Some(msg),
                                });
                            } else {
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: Some(literal_type),
                                    state: BindingState::Ok,
                                    message: Some("Literal binding".to_string()),
                                });
                            }
                        } else {
                            warnings.push(format!(
                                "Step '{}' input '{}' uses a literal binding",
                                resolved_step.step_id, input.name
                            ));
                            bindings.push(BindingStatus {
                                step_id: resolved_step.step_id.clone(),
                                input: input.name.clone(),
                                binding: Some(binding_value.clone()),
                                expected: input.raw_type.clone(),
                                actual: None,
                                state: BindingState::Warning,
                                message: Some("Literal value provided".to_string()),
                            });
                        }
                    } else {
                        bindings.push(BindingStatus {
                            step_id: resolved_step.step_id.clone(),
                            input: input.name.clone(),
                            binding: None,
                            expected: input.raw_type.clone(),
                            actual: None,
                            state: BindingState::Pending,
                            message: Some(
                                "Pending input – supply via CLI overrides or a downstream step"
                                    .to_string(),
                            ),
                        });
                    }
                }

                for output in &module_spec.outputs {
                    let key = format!("step.{}.outputs.{}", resolved_step.step_id, output.name);
                    available.insert(key, (resolved_step.step_id.clone(), output.clone()));
                }
                for (share_name, share_spec) in &step_spec.share {
                    let source_name = match parse_share_source(&share_spec.source) {
                        Ok(name) => name,
                        Err(e) => {
                            errors.push(format!(
                                "Share '{}' in step '{}': {}",
                                share_name, step_spec.id, e
                            ));
                            continue;
                        }
                    };
                    let Some(source_output) = module_spec
                        .outputs
                        .iter()
                        .find(|output| output.name == source_name)
                    else {
                        errors.push(format!(
                            "Share '{}' in step '{}' references unknown output '{}'",
                            share_name, step_spec.id, share_spec.source
                        ));
                        continue;
                    };
                    let mut shared_output = source_output.clone();
                    shared_output.name = share_name.clone();
                    // Use 'share' namespace to distinguish from module outputs
                    let key = format!("step.{}.share.{}", resolved_step.step_id, share_name);
                    available.insert(key, (resolved_step.step_id.clone(), shared_output));
                }

                let mut known_outputs: HashSet<String> =
                    module_spec.outputs.iter().map(|o| o.name.clone()).collect();
                for alias in step_spec.publish.keys() {
                    known_outputs.insert(alias.clone());
                }
                for share_name in step_spec.share.keys() {
                    known_outputs.insert(share_name.clone());
                }
                for (store_name, store_spec) in &step_spec.store {
                    match store_spec {
                        FlowStoreSpec::Sql(sql) => {
                            if !known_outputs.contains(&sql.source) {
                                errors.push(format!(
                                    "Store '{}' in step '{}' references unknown output '{}'",
                                    store_name, step_spec.id, sql.source
                                ));
                            }
                        }
                    }
                }

                resolved_steps.push(resolved_step);
            }
            Ok(None) => {
                // Skip adding to unresolved if this is a permissions-only or barrier step
                if step_spec.permissions.is_empty() && step_spec.barrier.is_none() {
                    unresolved_steps.push(UnresolvedStep {
                        step_id: step_spec.id.clone(),
                        reason: "Module not specified (provide via --module or edit flow)".to_string(),
                    });
                }
            }
            Err(err) => {
                errors.push(format!("Step '{}': {}", step_spec.id, err));
            }
        }
    }

    Ok(FlowValidationResult {
        resolved_steps,
        unresolved_steps,
        bindings,
        errors,
        warnings,
        flow_inputs: spec.inputs.clone(),
    })
}

fn print_validation(result: &FlowValidationResult, diagram: bool) {
    if !result.errors.is_empty() {
        println!("\n❌ Errors:");
        for err in &result.errors {
            println!("  - {}", err.red());
        }
    }

    if !result.warnings.is_empty() {
        println!("\n⚠️  Warnings:");
        for warn in &result.warnings {
            println!("  - {}", warn.yellow());
        }
    }

    if !result.flow_inputs.is_empty() {
        println!("\n🔑 Flow inputs:");
        for (name, spec_input) in &result.flow_inputs {
            let mut descriptor = spec_input.raw_type().to_string();
            if let Some(default) = spec_input.default_literal() {
                descriptor.push_str(&format!(" (default: {})", default));
            }
            println!("  - {} : {}", name.bold(), descriptor);
        }
    }

    if !result.unresolved_steps.is_empty() {
        println!("\nℹ️  Pending steps:");
        for step in &result.unresolved_steps {
            println!(
                "  - {} ({})",
                step.step_id.bold(),
                step.reason.as_str().dimmed()
            );
        }
    }

    let pending_bindings: Vec<_> = result
        .bindings
        .iter()
        .filter(|b| matches!(b.state, BindingState::Pending))
        .collect();

    if !pending_bindings.is_empty() {
        println!("\nℹ️  Pending inputs:");
        for binding in pending_bindings {
            println!(
                "  - Step {} input {} awaiting value",
                binding.step_id.bold(),
                binding.input.as_str().bold()
            );
        }
    }

    if diagram {
        print_steps(result);
    }
}

fn print_steps(result: &FlowValidationResult) {
    println!("\n📦 Steps:");
    for step in &result.resolved_steps {
        println!("  • {} → {}", step.step_id.bold(), step.module.summary());
        for binding in result.bindings.iter().filter(|b| b.step_id == step.step_id) {
            let status = match binding.state {
                BindingState::Ok => "ok".green(),
                BindingState::Warning => "warn".yellow(),
                BindingState::Error => "error".red(),
                BindingState::Pending => "pending".cyan(),
            };
            let binding_display = binding
                .binding
                .as_deref()
                .map(|b| {
                    if let Some(name) = b.strip_prefix("inputs.") {
                        format!("inputs.{}", name)
                    } else if is_type_placeholder(b, &binding.expected) {
                        "(pending)".to_string()
                    } else if let Some((_, value)) = parse_typed_literal(b) {
                        if value.is_empty() {
                            b.to_string()
                        } else {
                            value
                        }
                    } else {
                        b.to_string()
                    }
                })
                .unwrap_or_else(|| "(unbound)".to_string());
            let type_info = match &binding.actual {
                Some(actual) => format!("{} → {}", actual, binding.expected.clone()),
                None => binding.expected.clone(),
            };
            println!(
                "      {} ← {} [{}] {}",
                binding.input,
                binding_display.dimmed(),
                status,
                type_info.as_str().dimmed()
            );
            if !matches!(binding.state, BindingState::Ok) {
                if let Some(msg) = &binding.message {
                    println!("        {}", msg.as_str().dimmed());
                }
            }
        }
        if !step.store.is_empty() {
            println!("      store:");
            for (name, spec) in &step.store {
                match spec {
                    FlowStoreSpec::Sql(sql) => {
                        let table_display = sql.table.as_deref().unwrap_or("(auto)");
                        println!(
                            "        {} → sql(table: {}, source: {})",
                            name, table_display, sql.source
                        );
                    }
                }
            }
        }
        println!();
    }

    for pending in &result.unresolved_steps {
        println!(
            "  • {} (pending module) {}",
            pending.step_id.bold(),
            pending.reason.as_str().dimmed()
        );
    }
}

fn resolve_flow_path(file: Option<String>) -> PathBuf {
    file.map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("flow.yaml"))
}

fn types_compatible(expected: &str, actual: &str) -> bool {
    crate::module_spec::types_compatible(expected, actual)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_detect_table_name_default() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let result = detect_table_name("my_store", &spec, "20251023120000");
        assert_eq!(result, "z_results_my_store_20251023120000");
    }

    #[test]
    fn test_detect_table_name_custom_template() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: Some("custom_table_{run_id}".to_string()),
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let result = detect_table_name("my_store", &spec, "20251023120000");
        assert_eq!(result, "z_results_custom_table_20251023120000");
    }

    #[test]
    fn test_detect_table_name_no_run_id_substitution() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: Some("static_table".to_string()),
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let result = detect_table_name("my_store", &spec, "20251023120000");
        assert_eq!(result, "z_results_static_table");
    }

    #[test]
    fn test_detect_table_name_multiple_run_id() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: Some("tbl_{run_id}_v2_{run_id}".to_string()),
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let result = detect_table_name("my_store", &spec, "123");
        assert_eq!(result, "z_results_tbl_123_v2_123");
    }

    #[test]
    fn test_parse_sql_destination_empty() {
        let result = parse_sql_destination(None).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sql_destination_biovault() {
        let result = parse_sql_destination(Some("biovault")).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sql_destination_sql_empty() {
        let result = parse_sql_destination(Some("SQL()")).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sql_destination_sql_case_insensitive() {
        let result = parse_sql_destination(Some("sql()")).unwrap();
        assert_eq!(result, None);
        let result2 = parse_sql_destination(Some("Sql()")).unwrap();
        assert_eq!(result2, None);
    }

    #[test]
    fn test_parse_sql_destination_with_quotes() {
        let result = parse_sql_destination(Some("'SQL()'")).unwrap();
        assert_eq!(result, None);
        let result2 = parse_sql_destination(Some("\"biovault\"")).unwrap();
        assert_eq!(result2, None);
    }

    #[test]
    fn test_parse_sql_destination_url() {
        let result = parse_sql_destination(Some("SQL(url:postgres://localhost/db)")).unwrap();
        assert_eq!(result, Some("postgres://localhost/db".to_string()));
    }

    #[test]
    fn test_parse_sql_destination_url_with_equals() {
        let result = parse_sql_destination(Some("SQL(url=sqlite:///tmp/test.db)")).unwrap();
        assert_eq!(result, Some("sqlite:///tmp/test.db".to_string()));
    }

    #[test]
    fn test_parse_sql_destination_url_with_quotes() {
        let result = parse_sql_destination(Some("SQL(url:'postgres://localhost/db')")).unwrap();
        assert_eq!(result, Some("postgres://localhost/db".to_string()));
    }

    #[test]
    fn test_parse_sql_destination_invalid_syntax() {
        let result = parse_sql_destination(Some("SQL(invalid)"));
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Unsupported SQL destination syntax"));
    }

    #[test]
    fn test_parse_sql_destination_empty_url() {
        let result = parse_sql_destination(Some("SQL(url:)"));
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("url cannot be empty"));
    }

    #[test]
    fn test_parse_sql_destination_unsupported_prefix() {
        let result = parse_sql_destination(Some("POSTGRES()"));
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Unsupported SQL destination"));
    }

    #[test]
    fn test_sanitize_identifier_simple() {
        assert_eq!(sanitize_identifier("my_table"), "my_table");
    }

    #[test]
    fn test_sanitize_identifier_uppercase() {
        assert_eq!(sanitize_identifier("MyTable"), "mytable");
    }

    #[test]
    fn test_sanitize_identifier_spaces() {
        assert_eq!(sanitize_identifier("my table name"), "my_table_name");
    }

    #[test]
    fn test_sanitize_identifier_special_chars() {
        assert_eq!(sanitize_identifier("my-table@name!"), "my_table_name");
    }

    #[test]
    fn test_sanitize_identifier_multiple_underscores() {
        assert_eq!(sanitize_identifier("my___table"), "my_table");
    }

    #[test]
    fn test_sanitize_identifier_leading_trailing_underscores() {
        assert_eq!(sanitize_identifier("_table_"), "table");
    }

    #[test]
    fn test_sanitize_identifier_starts_with_number() {
        assert_eq!(sanitize_identifier("123table"), "t_123table");
    }

    #[test]
    fn test_sanitize_identifier_empty() {
        assert_eq!(sanitize_identifier(""), "t");
    }

    #[test]
    fn test_sanitize_identifier_only_special_chars() {
        assert_eq!(sanitize_identifier("@#$%"), "t");
    }

    #[test]
    fn test_sanitize_identifier_unicode() {
        assert_eq!(sanitize_identifier("my_tåble"), "my_t_ble");
    }

    #[test]
    fn test_detect_format_csv_extension() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let path = PathBuf::from("output.csv");
        assert_eq!(detect_format(&spec, &path), "csv");
    }

    #[test]
    fn test_detect_format_tsv_extension() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let path = PathBuf::from("output.tsv");
        assert_eq!(detect_format(&spec, &path), "tsv");
    }

    #[test]
    fn test_detect_format_uppercase_extension() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let path = PathBuf::from("output.CSV");
        assert_eq!(detect_format(&spec, &path), "csv");
    }

    #[test]
    fn test_detect_format_explicit_format() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: Some("TSV".to_string()),
        };
        let path = PathBuf::from("output.csv");
        assert_eq!(detect_format(&spec, &path), "tsv");
    }

    #[test]
    fn test_detect_format_no_extension() {
        let spec = FlowSqlStoreSpec {
            target: None,
            source: "output".to_string(),
            table: None,
            participant_column: None,
            overwrite: None,
            format: None,
        };
        let path = PathBuf::from("output");
        assert_eq!(detect_format(&spec, &path), "csv");
    }

    #[test]
    fn test_types_compatible_same() {
        assert!(types_compatible("File", "File"));
    }

    #[test]
    fn test_types_compatible_case_insensitive() {
        assert!(types_compatible("File", "file"));
    }

    #[test]
    fn test_types_compatible_optional() {
        assert!(types_compatible("File?", "File"));
        assert!(types_compatible("File", "File?"));
    }

    #[test]
    fn test_types_compatible_different() {
        assert!(!types_compatible("File", "Directory"));
    }

    #[test]
    fn test_types_compatible_map_and_record() {
        assert!(types_compatible(
            "Map[String, Record{bed: File, bim: File, fam: File}]",
            "map[string, record{fam: file, bed: file, bim: file}]"
        ));
        assert!(!types_compatible(
            "Map[String, Record{bed: File, bim: File, fam: File}]",
            "Map[String, File]"
        ));
    }

    #[test]
    fn test_resolve_flow_path_default() {
        assert_eq!(resolve_flow_path(None), PathBuf::from("flow.yaml"));
    }

    #[test]
    fn test_resolve_flow_path_custom() {
        assert_eq!(
            resolve_flow_path(Some("custom.yaml".to_string())),
            PathBuf::from("custom.yaml")
        );
    }

    #[test]
    fn test_get_delimiter_for_tsv() {
        assert_eq!(get_delimiter_for_format("tsv"), b'\t');
    }

    #[test]
    fn test_get_delimiter_for_csv() {
        assert_eq!(get_delimiter_for_format("csv"), b',');
    }

    #[test]
    fn test_get_delimiter_for_unknown() {
        assert_eq!(get_delimiter_for_format("unknown"), b',');
    }
}
