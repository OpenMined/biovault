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

use crate::data::BioVaultDb;
use crate::error::Result;
use crate::pipeline_spec::{
    value_to_string, PipelineInputSpec, PipelineSpec, PipelineSqlStoreSpec, PipelineStepSpec,
    PipelineStoreSpec,
};
use crate::project_spec::{InputSpec, OutputSpec, ProjectSpec};
use anyhow::{anyhow, Context};
use chrono::Utc;
use colored::Colorize;
use csv::ReaderBuilder;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use rusqlite::params_from_iter;
use serde_yaml::Value as YamlValue;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::path::{Path, PathBuf};
use tokio::fs;

use super::run_dynamic;

type StepOverrides = HashMap<(String, String), String>;
type PipelineOverrides = HashMap<String, String>;
type ParseOverridesResult = (
    StepOverrides,
    PipelineOverrides,
    Option<String>,
    Vec<String>,
);

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

pub async fn create(
    file: Option<String>,
    name: Option<String>,
    uses: Option<String>,
    step_id: Option<String>,
) -> Result<()> {
    let pipeline_path = resolve_pipeline_path(file);
    let db = BioVaultDb::new().ok();

    if pipeline_path.exists()
        && !Confirm::with_theme(&ColorfulTheme::default())
            .with_prompt(format!(
                "Pipeline file {} already exists. Overwrite?",
                pipeline_path.display()
            ))
            .default(false)
            .interact()
            .cli_result()?
    {
        println!("✋ Aborting pipeline creation.");
        return Ok(());
    }

    if let Some(uses_value) = uses {
        let pipeline_dir = pipeline_path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));

        let (reference, root) = normalize_project_reference(&uses_value, pipeline_dir)?;
        let choice = ProjectChoice::Path {
            reference: reference.clone(),
            root: root.clone(),
        };
        let project = load_project_spec(&choice)?;

        let mut spec = PipelineSpec::default();
        spec.name = name.unwrap_or_else(|| project.spec.name.clone());

        let generated_id = generate_default_step_id(&project, &spec);
        let step_id_value = step_id.unwrap_or(generated_id);

        let mut with_map = BTreeMap::new();
        for input in &project.spec.inputs {
            let key = ensure_pipeline_input(&mut spec, &step_id_value, input);
            with_map.insert(
                input.name.clone(),
                YamlValue::String(format!("inputs.{}", key)),
            );
        }

        let publish_map = default_publish_map(&project.spec.outputs);

        spec.steps.push(PipelineStepSpec {
            id: step_id_value,
            uses: Some(reference),
            where_exec: None,
            with: with_map,
            publish: publish_map,
            store: BTreeMap::new(),
        });

        spec.save(&pipeline_path)?;

        println!(
            "\n✅ Saved pipeline to {}",
            pipeline_path.display().to_string().bold()
        );

        let validation = validate_internal(&pipeline_path, &spec, db.as_ref())?;
        print_validation(&validation, true);

        return Ok(());
    }

    let wizard = PipelineWizard::new(&pipeline_path, db.as_ref());
    let spec = wizard.run()?;
    spec.save(&pipeline_path)?;

    println!(
        "\n✅ Saved pipeline to {}",
        pipeline_path.display().to_string().bold()
    );

    let validation = validate_internal(&pipeline_path, &spec, db.as_ref())?;
    print_validation(&validation, true);

    Ok(())
}

pub async fn add_step(file: Option<String>) -> Result<()> {
    let pipeline_path = resolve_pipeline_path(file);
    if !pipeline_path.exists() {
        return Err(anyhow!(
            "Pipeline file not found: {}. Run 'bv pipeline create' first.",
            pipeline_path.display()
        )
        .into());
    }

    let mut spec = PipelineSpec::load(&pipeline_path)?;
    let db = BioVaultDb::new().ok();

    let wizard = PipelineWizard::new(&pipeline_path, db.as_ref());
    spec = wizard.add_step_to(spec)?;
    spec.save(&pipeline_path)?;

    println!(
        "\n✅ Updated pipeline at {}",
        pipeline_path.display().to_string().bold()
    );

    let validation = validate_internal(&pipeline_path, &spec, db.as_ref())?;
    print_validation(&validation, true);

    Ok(())
}

pub async fn run_pipeline(
    pipeline_path: &str,
    extra_args: Vec<String>,
    dry_run: bool,
    resume: bool,
    results_dir_override: Option<String>,
) -> Result<()> {
    let path = Path::new(pipeline_path);
    if !path.exists() {
        return Err(anyhow!("Pipeline file not found: {}", pipeline_path).into());
    }

    let (step_overrides, pipeline_overrides, explicit_results_dir, nextflow_passthrough) =
        parse_overrides(&extra_args, results_dir_override.clone())?;

    let spec = PipelineSpec::load(path)?;
    let mut db = BioVaultDb::new().ok();
    let validation = validate_internal(path, &spec, db.as_ref())?;

    if validation.has_errors() {
        print_validation(&validation, true);
        return Err(anyhow!("Pipeline has validation errors. Fix them before running.").into());
    }

    if !validation.unresolved_steps.is_empty() {
        return Err(anyhow!(
            "Pipeline contains steps without a project. Set 'uses' or register the project before running."
        )
        .into());
    }

    let pipeline_dir = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    let run_id = Utc::now().format("%Y%m%d%H%M%S").to_string();
    let run_msg = format!("🆔 Pipeline run {}", run_id);
    println!("{}", run_msg);
    append_desktop_log(&run_msg);

    for key in pipeline_overrides.keys() {
        if !spec.inputs.contains_key(key) {
            return Err(anyhow!("Unknown pipeline input '{}' in --set", key).into());
        }
    }

    let mut resolved_inputs: HashMap<String, String> = HashMap::new();
    for (name, input_spec) in &spec.inputs {
        if let Some(value) = pipeline_overrides.get(name) {
            let resolved = literal_to_value(value, input_spec.raw_type())?;
            resolved_inputs.insert(name.clone(), resolved);
        } else if let Some(default_literal) = input_spec.default_literal() {
            let resolved = literal_to_value(default_literal, input_spec.raw_type())?;
            resolved_inputs.insert(name.clone(), resolved);
        }
    }

    let mut step_outputs: HashMap<String, HashMap<String, String>> = HashMap::new();

    let requested_results_dir = explicit_results_dir.or(results_dir_override);

    let base_results_dir = match requested_results_dir {
        Some(dir) => PathBuf::from(dir),
        None => {
            let mut base = PathBuf::from("results/pipelines");
            base.push(&spec.name);
            base
        }
    };

    if !dry_run {
        fs::create_dir_all(&base_results_dir).await?;
    }

    for step in &spec.steps {
        let project = resolve_project(step, pipeline_dir, db.as_ref())?
            .ok_or_else(|| anyhow!("Step '{}' has no project to run", step.id))?;

        let project_root = project.root.to_string_lossy().to_string();
        let project_spec = &project.spec;

        let mut step_args = Vec::new();

        for input in &project_spec.inputs {
            let binding = step
                .with
                .get(&input.name)
                .and_then(value_to_string)
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

            let resolved_value = resolve_binding(
                &resolved_binding,
                &input.raw_type,
                &spec.inputs,
                &resolved_inputs,
                &step_outputs,
                &step.id,
            )?;

            step_args.push(format!("--{}", input.name));
            step_args.push(resolved_value);
        }

        if !nextflow_passthrough.is_empty() {
            step_args.extend(nextflow_passthrough.clone());
        }

        let step_results_dir = base_results_dir.join(&step.id);
        let step_results_dir_str = step_results_dir.to_string_lossy().to_string();

        println!(
            "\n▶️  Running step {} using project {}",
            step.id.bold(),
            project_root
        );

        append_desktop_log(&format!(
            "[Pipeline] Running step {} with args: {:?}",
            step.id, step_args
        ));
        run_dynamic::execute_dynamic(
            &project_root,
            step_args.clone(),
            dry_run,
            resume,
            Some(step_results_dir_str.clone()),
        )
        .await?;

        append_desktop_log(&format!("[Pipeline] Completed step {}", step.id));

        let mut outputs = HashMap::new();
        for output in &project_spec.outputs {
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

        step_outputs.insert(step.id.clone(), outputs);

        if !step.store.is_empty() {
            let db_conn = db
                .as_mut()
                .ok_or_else(|| anyhow!("BioVault database not available for store operations"))?;
            let step_outputs_ref = step_outputs.get(&step.id).unwrap();
            for (store_name, store_spec) in &step.store {
                match store_spec {
                    PipelineStoreSpec::Sql(sql) => {
                        if let Some(url) = parse_sql_destination(sql.target.as_deref())? {
                            return Err(anyhow!(
                                "SQL store '{}' destination '{}' not supported yet (only built-in biovault is available)",
                                store_name,
                                url
                            )
                            .into());
                        }
                        let source_path = step_outputs_ref.get(&sql.source).ok_or_else(|| {
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

    println!("\n✅ Pipeline run completed successfully");
    println!("Results stored under {}", base_results_dir.display());
    Ok(())
}

pub fn validate(file: &str, diagram: bool) -> Result<()> {
    let pipeline_path = Path::new(file);
    if !pipeline_path.exists() {
        return Err(anyhow!("Pipeline file not found: {}", pipeline_path.display()).into());
    }

    let spec = PipelineSpec::load(pipeline_path)?;
    let db = BioVaultDb::new().ok();
    let validation = validate_internal(pipeline_path, &spec, db.as_ref())?;

    print_validation(&validation, diagram);

    if validation.has_errors() {
        Err(anyhow!("Pipeline validation reported errors").into())
    } else {
        println!("\n✅ Pipeline is valid");
        Ok(())
    }
}

pub fn inspect(file: &str) -> Result<()> {
    let pipeline_path = Path::new(file);
    if !pipeline_path.exists() {
        return Err(anyhow!("Pipeline file not found: {}", pipeline_path.display()).into());
    }
    let spec = PipelineSpec::load(pipeline_path)?;
    let db = BioVaultDb::new().ok();
    let validation = validate_internal(pipeline_path, &spec, db.as_ref())?;
    print_steps(&validation);
    Ok(())
}

struct PipelineWizard<'a> {
    pipeline_path: &'a Path,
    db: Option<&'a BioVaultDb>,
    theme: ColorfulTheme,
}

impl<'a> PipelineWizard<'a> {
    fn new(pipeline_path: &'a Path, db: Option<&'a BioVaultDb>) -> Self {
        Self {
            pipeline_path,
            db,
            theme: ColorfulTheme::default(),
        }
    }

    fn run(&self) -> Result<PipelineSpec> {
        println!("🧪 Pipeline wizard — let’s assemble your pipeline");

        let spec_name = Input::with_theme(&self.theme)
            .with_prompt("Pipeline name")
            .validate_with(|input: &String| {
                if input.trim().is_empty() {
                    Err("Pipeline name cannot be empty")
                } else {
                    Ok(())
                }
            })
            .interact_text()
            .cli_result()?;

        let mut spec = PipelineSpec {
            name: spec_name,
            ..PipelineSpec::default()
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

    fn add_step_to(&self, mut spec: PipelineSpec) -> Result<PipelineSpec> {
        let resolved_steps = resolve_existing_steps(&spec, self.pipeline_dir(), self.db)?;
        let available_outputs = collect_available_outputs(&resolved_steps);

        let project_choice = self.prompt_project_choice()?;

        let loaded_project = match &project_choice {
            ProjectChoice::Defer => None,
            _ => Some(load_project_spec(&project_choice)?),
        };

        if let Some(project) = &loaded_project {
            println!(
                "\nConfiguring step for project {}",
                project.summary().bold()
            );
        } else {
            println!("\nConfiguring step without a project (assign later)");
        }

        let default_id = if let Some(project) = &loaded_project {
            generate_default_step_id(project, &spec)
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

        if let Some(project) = &loaded_project {
            for input in &project.spec.inputs {
                if let Some(binding) =
                    self.prompt_input_binding(&mut spec, &step_id, input, &available_outputs)?
                {
                    with_map.insert(input.name.clone(), binding);
                }
            }

            for input in &project.spec.inputs {
                with_map.entry(input.name.clone()).or_insert_with(|| {
                    let key = ensure_pipeline_input(&mut spec, &step_id, input);
                    YamlValue::String(format!("inputs.{}", key))
                });
            }

            publish_map = prompt_publish_map(&self.theme, &project.spec.outputs)?;
        }

        let step = PipelineStepSpec {
            id: step_id,
            uses: project_choice.uses_value(),
            where_exec: None,
            with: with_map,
            publish: publish_map,
            store: BTreeMap::new(),
        };

        spec.steps.push(step);
        Ok(spec)
    }

    fn prompt_project_choice(&self) -> Result<ProjectChoice> {
        let mut options = list_registered_projects(self.db);
        options.push(ProjectChoice::EnterPath);
        options.push(ProjectChoice::Defer);

        let items: Vec<String> = options.iter().map(|choice| choice.display()).collect();

        let index = Select::with_theme(&self.theme)
            .with_prompt("Select a project to use")
            .items(&items)
            .default(0)
            .interact()
            .cli_result()?;

        match &options[index] {
            ProjectChoice::EnterPath => {
                let input: String = Input::with_theme(&self.theme)
                    .with_prompt("Project directory path or project.yaml")
                    .interact_text()
                    .cli_result()?;
                let (reference, root) = normalize_project_reference(&input, self.pipeline_dir())?;
                Ok(ProjectChoice::Path { reference, root })
            }
            ProjectChoice::Defer => Ok(ProjectChoice::Defer),
            other => Ok(other.clone()),
        }
    }

    fn prompt_input_binding(
        &self,
        spec: &mut PipelineSpec,
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
        let pipeline_index = options.len();
        options.push("Use pipeline input".to_string());

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

        if selection == pipeline_index {
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
                ensure_pipeline_input(spec, step_id, input)
            } else {
                let mut pipeline_options = matching_inputs
                    .iter()
                    .map(|name| format!("Existing: {}", name))
                    .collect::<Vec<_>>();
                pipeline_options.push(format!("Create new input for '{}'", input.name));

                let selected = Select::with_theme(&self.theme)
                    .with_prompt("Pipeline input")
                    .items(&pipeline_options)
                    .default(0)
                    .interact()
                    .cli_result()?;

                if selected < matching_inputs.len() {
                    matching_inputs[selected].clone()
                } else {
                    ensure_pipeline_input(spec, step_id, input)
                }
            };

            return Ok(Some(YamlValue::String(format!("inputs.{}", chosen_key))));
        }

        Ok(None)
    }

    fn pipeline_dir(&self) -> &Path {
        self.pipeline_path
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
enum ProjectChoice {
    Registered { name: String, path: PathBuf },
    Path { reference: String, root: PathBuf },
    EnterPath,
    Defer,
}

impl ProjectChoice {
    fn display(&self) -> String {
        match self {
            ProjectChoice::Registered { name, path } => {
                format!("{} (registered at {})", name, path.display())
            }
            ProjectChoice::Path { reference, .. } => reference.clone(),
            ProjectChoice::EnterPath => "Enter project path".to_string(),
            ProjectChoice::Defer => "Skip (assign at run time)".to_string(),
        }
    }

    fn uses_value(&self) -> Option<String> {
        match self {
            ProjectChoice::Registered { name, .. } => Some(name.clone()),
            ProjectChoice::Path { reference, .. } => Some(reference.clone()),
            ProjectChoice::EnterPath => None,
            ProjectChoice::Defer => None,
        }
    }
}

struct LoadedProject {
    spec: ProjectSpec,
    root: PathBuf,
}

impl LoadedProject {
    fn summary(&self) -> String {
        format!("{} ({})", self.spec.name, self.root.display())
    }
}

fn load_project_spec(choice: &ProjectChoice) -> Result<LoadedProject> {
    match choice {
        ProjectChoice::Registered { name: _, path } => {
            let spec_path = Path::new(path).join("project.yaml");
            let spec = ProjectSpec::load(&spec_path)?;
            Ok(LoadedProject {
                spec,
                root: path.clone(),
            })
        }
        ProjectChoice::Path { reference: _, root } => {
            let spec_path = Path::new(root).join("project.yaml");
            if !spec_path.exists() {
                return Err(anyhow!("No project.yaml found at {}", spec_path.display()).into());
            }
            let spec = ProjectSpec::load(&spec_path)?;
            Ok(LoadedProject {
                spec,
                root: root.clone(),
            })
        }
        ProjectChoice::EnterPath => Err(anyhow!("Unexpected project choice").into()),
        ProjectChoice::Defer => Err(anyhow!("Deferred project has no spec").into()),
    }
}

fn list_registered_projects(db: Option<&BioVaultDb>) -> Vec<ProjectChoice> {
    if let Some(db) = db {
        match db.list_projects() {
            Ok(projects) => projects
                .into_iter()
                .map(|p| ProjectChoice::Registered {
                    name: p.name,
                    path: PathBuf::from(p.project_path),
                })
                .collect(),
            Err(_) => vec![],
        }
    } else {
        vec![]
    }
}

fn generate_default_step_id(project: &LoadedProject, spec: &PipelineSpec) -> String {
    let mut base: String = project
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
        for output in &step.project.spec.outputs {
            outputs.push(CandidateOutput {
                binding: format!("step.{}.outputs.{}", step.step_id, output.name),
                output: output.clone(),
            });
        }
    }
    outputs
}

fn ensure_pipeline_input(spec: &mut PipelineSpec, step_id: &str, input: &InputSpec) -> String {
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
            "⚠️  Input '{}' on step '{}' conflicts with existing pipeline input; recorded as '{}'.",
            input.name, step_id, key
        );
    }

    spec.inputs
        .insert(key.clone(), PipelineInputSpec::from_type(&input.raw_type));
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
    normalize_type(trimmed) == normalize_type(expected)
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
    let mut pipeline_overrides = HashMap::new();
    let mut results_dir = initial_results_dir;
    let mut nextflow_args = Vec::new();
    let mut i = 0;
    while i < extra_args.len() {
        let arg = &extra_args[i];
        if arg == "--set" {
            if i + 1 >= extra_args.len() {
                return Err(anyhow!("--set requires an argument like step.input=value").into());
            }
            parse_override_pair(
                &extra_args[i + 1],
                &mut step_overrides,
                &mut pipeline_overrides,
            )?;
            i += 2;
        } else if let Some(kv) = arg.strip_prefix("--set=") {
            parse_override_pair(kv, &mut step_overrides, &mut pipeline_overrides)?;
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
    Ok((
        step_overrides,
        pipeline_overrides,
        results_dir,
        nextflow_args,
    ))
}

fn parse_override_pair(
    pair: &str,
    step_overrides: &mut HashMap<(String, String), String>,
    pipeline_overrides: &mut HashMap<String, String>,
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
        pipeline_overrides.insert(input.trim().to_string(), value);
    } else {
        step_overrides.insert((step.trim().to_string(), input.trim().to_string()), value);
    }
    Ok(())
}

fn resolve_binding(
    binding: &str,
    expected_type: &str,
    pipeline_inputs: &BTreeMap<String, PipelineInputSpec>,
    resolved_inputs: &HashMap<String, String>,
    step_outputs: &HashMap<String, HashMap<String, String>>,
    current_step_id: &str,
) -> Result<String> {
    if let Some(input_name) = binding.strip_prefix("inputs.") {
        if let Some(value) = resolved_inputs.get(input_name) {
            return Ok(value.clone());
        }
        if let Some(_spec) = pipeline_inputs.get(input_name) {
            return Err(anyhow!(
                "Pipeline input '{}' is not set. Provide a value with --set inputs.{}=<value>.",
                input_name,
                input_name
            )
            .into());
        } else {
            return Err(anyhow!(
                "Pipeline input '{}' referenced in step '{}' is not declared",
                input_name,
                current_step_id
            )
            .into());
        }
    }

    if binding.starts_with("step.") {
        let parts: Vec<&str> = binding.split('.').collect();
        if parts.len() != 4 || parts[2] != "outputs" {
            return Err(anyhow!(
                "Invalid step output reference '{}' in step '{}'. Expected format step.<id>.outputs.<name>",
                binding, current_step_id
            )
            .into());
        }
        let source_step = parts[1];
        let output_name = parts[3];
        let outputs = step_outputs.get(source_step).ok_or_else(|| {
            anyhow!(
                "Step '{}' references outputs from step '{}' which has not run yet",
                current_step_id,
                source_step
            )
        })?;
        let value = outputs.get(output_name).ok_or_else(|| {
            anyhow!(
                "Step '{}' output '{}' not found when referenced from step '{}'",
                source_step,
                output_name,
                current_step_id
            )
        })?;
        return Ok(value.clone());
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

fn detect_table_name(store_name: &str, spec: &PipelineSqlStoreSpec, run_id: &str) -> String {
    let template = spec
        .table
        .clone()
        .unwrap_or_else(|| format!("{}_{}", store_name, run_id));
    let table_name = template.replace("{run_id}", run_id);
    format!("{}{}", RESULTS_TABLE_PREFIX, table_name)
}

fn detect_format(spec: &PipelineSqlStoreSpec, output_path: &Path) -> String {
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
    spec: &PipelineSqlStoreSpec,
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
    project: LoadedProject,
    store: BTreeMap<String, PipelineStoreSpec>,
}

fn resolve_existing_steps(
    spec: &PipelineSpec,
    pipeline_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<Vec<ResolvedStep>> {
    let mut resolved = Vec::new();
    for step in &spec.steps {
        if let Some(project) = resolve_project(step, pipeline_dir, db)? {
            resolved.push(ResolvedStep {
                step_id: step.id.clone(),
                project,
                store: step.store.clone(),
            });
        }
    }
    Ok(resolved)
}

fn normalize_project_reference(raw: &str, base_dir: &Path) -> Result<(String, PathBuf)> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(anyhow!("Project path cannot be empty").into());
    }

    let mut reference_path = PathBuf::from(trimmed);
    if reference_path
        .file_name()
        .map(|name| name == "project.yaml")
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

    if !absolute.join("project.yaml").exists() {
        return Err(anyhow!(
            "No project.yaml found at {}",
            absolute.join("project.yaml").display()
        )
        .into());
    }

    // Normalize absolute path by removing trailing components like '.'
    if let Ok(canon) = absolute.canonicalize() {
        absolute = canon;
    }

    Ok((reference_string, absolute))
}

fn resolve_project(
    step: &PipelineStepSpec,
    pipeline_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<Option<LoadedProject>> {
    match step.uses.as_deref().map(str::trim) {
        None | Some("") => Ok(None),
        Some(uses) => resolve_project_by_uses(uses, pipeline_dir, db).map(Some),
    }
}

fn resolve_project_by_uses(
    uses: &str,
    pipeline_dir: &Path,
    db: Option<&BioVaultDb>,
) -> Result<LoadedProject> {
    // Try 1: Resolve as path (relative or absolute)
    if let Ok((reference, root)) = normalize_project_reference(uses, pipeline_dir) {
        return load_project_spec(&ProjectChoice::Path { reference, root });
    }

    // Try 2: Database lookup with exact name
    if let Some(db) = db {
        if let Ok(Some(project)) = db.get_project(uses) {
            return load_project_spec(&ProjectChoice::Registered {
                name: project.name,
                path: PathBuf::from(project.project_path),
            });
        }
    }

    Err(anyhow!("Unable to resolve project reference '{}'", uses).into())
}

struct PipelineValidationResult {
    resolved_steps: Vec<ResolvedStep>,
    unresolved_steps: Vec<UnresolvedStep>,
    bindings: Vec<BindingStatus>,
    errors: Vec<String>,
    warnings: Vec<String>,
    pipeline_inputs: BTreeMap<String, PipelineInputSpec>,
}

impl PipelineValidationResult {
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

fn validate_internal(
    pipeline_path: &Path,
    spec: &PipelineSpec,
    db: Option<&BioVaultDb>,
) -> Result<PipelineValidationResult> {
    spec.ensure_unique_step_ids()?;
    let pipeline_dir = pipeline_path
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
        match resolve_project(step_spec, pipeline_dir, db) {
            Ok(Some(project)) => {
                let resolved_step = ResolvedStep {
                    step_id: step_spec.id.clone(),
                    project,
                    store: step_spec.store.clone(),
                };

                let project_spec = &resolved_step.project.spec;

                for input in &project_spec.inputs {
                    let raw_binding = step_spec.with.get(&input.name).cloned();
                    let binding_str = raw_binding.as_ref().and_then(value_to_string);

                    if let Some(binding_value) = &binding_str {
                        if let Some((producer_step, output)) = available.get(binding_value) {
                            if !types_compatible(&input.raw_type, &output.raw_type) {
                                let msg = format!(
                                    "Type mismatch: expected {} but found {} from step {}",
                                    input.raw_type, output.raw_type, producer_step
                                );
                                errors.push(msg.clone());
                                bindings.push(BindingStatus {
                                    step_id: resolved_step.step_id.clone(),
                                    input: input.name.clone(),
                                    binding: Some(binding_value.clone()),
                                    expected: input.raw_type.clone(),
                                    actual: Some(output.raw_type.clone()),
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
                                actual: Some(output.raw_type.clone()),
                                state: BindingState::Ok,
                                message: None,
                            });
                            continue;
                        }

                        if let Some(input_name) = binding_value.strip_prefix("inputs.") {
                            if let Some(spec_input) = spec.inputs.get(input_name) {
                                if !types_compatible(&input.raw_type, spec_input.raw_type()) {
                                    let msg = format!(
                                        "Type mismatch: expected {} but pipeline input '{}' declares {}",
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
                                            "Depends on pipeline input '{}'. Set it via --set inputs.{}=<value>.",
                                            input_name,
                                            input_name
                                        )),
                                    });
                                }
                            } else {
                                errors.push(format!(
                                    "Binding '{}' for step '{}' input '{}' references unknown pipeline input",
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
                                    message: Some("Unknown pipeline input".to_string()),
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
                                    "Pending input – assign via CLI override or pipeline edit"
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

                for output in &project_spec.outputs {
                    let key = format!("step.{}.outputs.{}", resolved_step.step_id, output.name);
                    available.insert(key, (resolved_step.step_id.clone(), output.clone()));
                }

                let mut known_outputs: HashSet<String> = project_spec
                    .outputs
                    .iter()
                    .map(|o| o.name.clone())
                    .collect();
                for alias in step_spec.publish.keys() {
                    known_outputs.insert(alias.clone());
                }
                for (store_name, store_spec) in &step_spec.store {
                    match store_spec {
                        PipelineStoreSpec::Sql(sql) => {
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
                unresolved_steps.push(UnresolvedStep {
                    step_id: step_spec.id.clone(),
                    reason: "Project not specified (provide via --project or edit pipeline)"
                        .to_string(),
                });
            }
            Err(err) => {
                errors.push(format!("Step '{}': {}", step_spec.id, err));
            }
        }
    }

    Ok(PipelineValidationResult {
        resolved_steps,
        unresolved_steps,
        bindings,
        errors,
        warnings,
        pipeline_inputs: spec.inputs.clone(),
    })
}

fn print_validation(result: &PipelineValidationResult, diagram: bool) {
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

    if !result.pipeline_inputs.is_empty() {
        println!("\n🔑 Pipeline inputs:");
        for (name, spec_input) in &result.pipeline_inputs {
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

fn print_steps(result: &PipelineValidationResult) {
    println!("\n📦 Steps:");
    for step in &result.resolved_steps {
        println!("  • {} → {}", step.step_id.bold(), step.project.summary());
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
                    PipelineStoreSpec::Sql(sql) => {
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
            "  • {} (pending project) {}",
            pending.step_id.bold(),
            pending.reason.as_str().dimmed()
        );
    }
}

fn resolve_pipeline_path(file: Option<String>) -> PathBuf {
    file.map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("pipeline.yaml"))
}

fn types_compatible(expected: &str, actual: &str) -> bool {
    normalize_type(expected) == normalize_type(actual)
}

fn normalize_type(value: &str) -> String {
    value.trim().trim_end_matches('?').to_ascii_lowercase()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_detect_table_name_default() {
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
        let spec = PipelineSqlStoreSpec {
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
    fn test_normalize_type_simple() {
        assert_eq!(normalize_type("File"), "file");
    }

    #[test]
    fn test_normalize_type_optional() {
        assert_eq!(normalize_type("File?"), "file");
    }

    #[test]
    fn test_normalize_type_with_whitespace() {
        assert_eq!(normalize_type("  File?  "), "file");
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
    fn test_resolve_pipeline_path_default() {
        assert_eq!(resolve_pipeline_path(None), PathBuf::from("pipeline.yaml"));
    }

    #[test]
    fn test_resolve_pipeline_path_custom() {
        assert_eq!(
            resolve_pipeline_path(Some("custom.yaml".to_string())),
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
