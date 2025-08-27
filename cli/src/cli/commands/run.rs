use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

#[derive(Debug, Serialize, Deserialize)]
struct ProjectConfig {
    name: String,
    author: String,
    workflow: String,
    #[serde(default, deserialize_with = "deserialize_string_or_vec")]
    assets: Vec<String>,
    #[serde(default)]
    patients: Vec<String>,
}

fn deserialize_string_or_vec<'de, D>(deserializer: D) -> std::result::Result<Vec<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::{self, Visitor};
    use std::fmt;

    struct StringOrVec;

    impl<'de> Visitor<'de> for StringOrVec {
        type Value = Vec<String>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("string or list of strings")
        }

        fn visit_str<E>(self, value: &str) -> std::result::Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(vec![value.to_string()])
        }

        fn visit_seq<A>(self, mut seq: A) -> std::result::Result<Self::Value, A::Error>
        where
            A: de::SeqAccess<'de>,
        {
            let mut vec = Vec::new();
            while let Some(value) = seq.next_element()? {
                vec.push(value);
            }
            Ok(vec)
        }
    }

    deserializer.deserialize_any(StringOrVec)
}

#[derive(Debug, Serialize, Deserialize)]
struct PatientData {
    ref_version: String,
    #[serde(rename = "ref")]
    ref_path: String,
    ref_index: String,
    aligned: String,
    aligned_index: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct PatientFile {
    patient: std::collections::HashMap<String, PatientData>,
}

pub async fn execute(
    project_folder: &str,
    patient_file: &str,
    patients: Option<Vec<String>>,
    patient: Option<String>,
    all: bool,
    test: bool,
    dry_run: bool,
    with_docker: bool,
    work_dir: Option<String>,
    resume: bool,
) -> Result<()> {
    // Validate project directory
    let project_path = PathBuf::from(project_folder);
    if !project_path.exists() {
        return Err(anyhow::anyhow!("Project folder does not exist: {}", project_folder).into());
    }

    let project_yaml = project_path.join("project.yaml");
    if !project_yaml.exists() {
        return Err(anyhow::anyhow!("project.yaml not found in {}", project_folder).into());
    }

    let workflow_file = project_path
        .join("workflow.nf")
        .canonicalize()
        .with_context(|| format!("Failed to resolve workflow.nf path in {}", project_folder))?;
    if !workflow_file.exists() {
        return Err(anyhow::anyhow!("workflow.nf not found in {}", project_folder).into());
    }

    // Parse project configuration
    let project_content = fs::read_to_string(&project_yaml).with_context(|| {
        format!(
            "Failed to read project.yaml from {}",
            project_yaml.display()
        )
    })?;
    let project_config: ProjectConfig =
        serde_yaml::from_str(&project_content).context("Failed to parse project.yaml")?;

    // Parse patient file
    let patient_file_path = PathBuf::from(patient_file);
    if !patient_file_path.exists() {
        return Err(anyhow::anyhow!("Patient file does not exist: {}", patient_file).into());
    }

    // Get the directory containing the patient file for resolving relative paths
    let patient_dir = patient_file_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Could not determine patient file directory"))?;

    let patient_content = fs::read_to_string(&patient_file_path).with_context(|| {
        format!(
            "Failed to read patient file from {}",
            patient_file_path.display()
        )
    })?;
    let patient_data: PatientFile =
        serde_yaml::from_str(&patient_content).context("Failed to parse patient file")?;

    // Determine which patients to process
    let patients_to_run =
        determine_patients(&patient_data, &project_config, patients, patient, all, test)?;

    if patients_to_run.is_empty() {
        println!("No patients to process");
        return Ok(());
    }

    println!("Processing {} patient(s):", patients_to_run.len());
    for p in &patients_to_run {
        println!("  - {}", p);
    }

    // Get BioVault environment directory
    let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
        PathBuf::from(test_home)
    } else {
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
    };
    let env_dir = home_dir.join(".biovault").join("env").join("default");

    // Check if templates exist
    let template_nf = env_dir.join("template.nf");
    let nextflow_config = env_dir.join("nextflow.config");

    if !template_nf.exists() || !nextflow_config.exists() {
        return Err(
            anyhow::anyhow!("Nextflow templates not found. Please run 'bv init' first").into(),
        );
    }

    // Create temporary directory for execution
    let temp_dir = TempDir::new().context("Failed to create temp directory")?;

    // Copy template.nf and nextflow.config to temp directory
    let temp_template = temp_dir.path().join("template.nf");
    let temp_config = temp_dir.path().join("nextflow.config");
    fs::copy(&template_nf, &temp_template).context("Failed to copy template.nf")?;
    fs::copy(&nextflow_config, &temp_config).context("Failed to copy nextflow.config")?;

    // Determine assets directory
    let assets_dir = if !project_config.assets.is_empty() {
        project_path.join(&project_config.assets[0])
    } else {
        project_path.join("assets")
    };

    // Create assets directory if it doesn't exist
    if !assets_dir.exists() {
        fs::create_dir_all(&assets_dir).with_context(|| {
            format!(
                "Failed to create assets directory: {}",
                assets_dir.display()
            )
        })?;
    }

    // Get absolute path for assets directory
    let assets_dir = assets_dir.canonicalize().with_context(|| {
        format!(
            "Failed to resolve assets directory path: {}",
            assets_dir.display()
        )
    })?;

    // Process each patient
    let mut success_count = 0;
    let mut fail_count = 0;

    for patient_id in patients_to_run {
        println!("\n{}", "=".repeat(60));
        println!("Processing patient: {}", patient_id);
        println!("{}", "=".repeat(60));

        // Get patient data
        let patient_info = patient_data
            .patient
            .get(&patient_id)
            .ok_or_else(|| anyhow::anyhow!("Patient {} not found in patient file", patient_id))?;

        // Create results directory for this patient
        let results_dir = project_path.join("results").join(&patient_id);
        if !results_dir.exists() {
            fs::create_dir_all(&results_dir).with_context(|| {
                format!(
                    "Failed to create results directory: {}",
                    results_dir.display()
                )
            })?;
        }

        // Get absolute path for results directory
        let results_dir = results_dir.canonicalize().with_context(|| {
            format!(
                "Failed to resolve results directory path: {}",
                results_dir.display()
            )
        })?;

        // Resolve patient data paths relative to patient file location
        let ref_path = patient_dir.join(&patient_info.ref_path);
        let ref_index_path = patient_dir.join(&patient_info.ref_index);
        let aligned_path = patient_dir.join(&patient_info.aligned);
        let aligned_index_path = patient_dir.join(&patient_info.aligned_index);

        // Verify that resolved paths exist
        if !ref_path.exists() {
            return Err(anyhow::anyhow!(
                "Reference file not found: {} (resolved to {})",
                patient_info.ref_path,
                ref_path.display()
            )
            .into());
        }
        if !aligned_path.exists() {
            return Err(anyhow::anyhow!(
                "Aligned file not found: {} (resolved to {})",
                patient_info.aligned,
                aligned_path.display()
            )
            .into());
        }

        // Build nextflow command
        let mut cmd = Command::new("nextflow");

        // Set working directory to project directory instead of temp directory
        cmd.current_dir(&project_path);

        cmd.arg("run")
            .arg(&temp_template)
            .arg("--patient_id")
            .arg(&patient_id)
            .arg("--ref_version")
            .arg(&patient_info.ref_version)
            .arg("--ref")
            .arg(ref_path.to_string_lossy().as_ref())
            .arg("--ref_index")
            .arg(ref_index_path.to_string_lossy().as_ref())
            .arg("--aligned")
            .arg(aligned_path.to_string_lossy().as_ref())
            .arg("--aligned_index")
            .arg(aligned_index_path.to_string_lossy().as_ref())
            .arg("--work_flow_file")
            .arg(workflow_file.to_string_lossy().as_ref())
            .arg("--assets_dir")
            .arg(assets_dir.to_string_lossy().as_ref())
            .arg("--results_dir")
            .arg(results_dir.to_string_lossy().as_ref());

        // Add work directory if specified
        if let Some(ref work_dir_path) = work_dir {
            cmd.arg("-work-dir").arg(work_dir_path);
        }

        // Add resume flag if specified
        if resume {
            cmd.arg("-resume");
        }

        // Add Docker flag
        if with_docker {
            cmd.arg("-with-docker");
        }

        // Add config
        cmd.arg("-c").arg(&temp_config);

        // Show command if dry-run or normal execution
        let cmd_str = format!("{:?}", cmd);
        println!("\nCommand: {}", cmd_str);

        if dry_run {
            println!("[DRY RUN] Would execute the above command");
            success_count += 1;
        } else {
            println!("\nExecuting nextflow pipeline...\n");

            // Execute the command
            let output = cmd.output().context("Failed to execute nextflow command")?;

            if output.status.success() {
                println!("✓ Successfully processed patient: {}", patient_id);
                success_count += 1;

                // Print stdout if available
                if !output.stdout.is_empty() {
                    println!("\nNextflow output:");
                    println!("{}", String::from_utf8_lossy(&output.stdout));
                }
            } else {
                println!("✗ Failed to process patient: {}", patient_id);
                fail_count += 1;

                // Print stderr if available
                if !output.stderr.is_empty() {
                    eprintln!("\nNextflow error:");
                    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
                }
            }
        }
    }

    // Print summary
    println!("\n{}", "=".repeat(60));
    println!("Execution Summary");
    println!("{}", "=".repeat(60));
    println!("Successful: {}", success_count);
    println!("Failed: {}", fail_count);
    println!("Total: {}", success_count + fail_count);

    if fail_count > 0 {
        return Err(anyhow::anyhow!("{} patient(s) failed processing", fail_count).into());
    }

    Ok(())
}

fn determine_patients(
    patient_data: &PatientFile,
    project_config: &ProjectConfig,
    patients_arg: Option<Vec<String>>,
    patient_arg: Option<String>,
    all: bool,
    test: bool,
) -> Result<Vec<String>> {
    // Priority 1: Test mode
    if test {
        if !patient_data.patient.contains_key("TEST") {
            return Err(anyhow::anyhow!("TEST patient not found in patient file").into());
        }
        return Ok(vec!["TEST".to_string()]);
    }

    // Priority 2: Command-line patient override (single)
    if let Some(patient_id) = patient_arg {
        if !patient_data.patient.contains_key(&patient_id) {
            return Err(anyhow::anyhow!("Patient {} not found in patient file", patient_id).into());
        }
        return Ok(vec![patient_id]);
    }

    // Priority 3: Command-line patients override (multiple)
    if let Some(patient_list) = patients_arg {
        for patient_id in &patient_list {
            if !patient_data.patient.contains_key(patient_id) {
                return Err(
                    anyhow::anyhow!("Patient {} not found in patient file", patient_id).into(),
                );
            }
        }
        return Ok(patient_list);
    }

    // Priority 4: All patients
    if all {
        return Ok(patient_data.patient.keys().cloned().collect());
    }

    // Priority 5: Project-defined patients
    if !project_config.patients.is_empty() {
        for patient_id in &project_config.patients {
            if !patient_data.patient.contains_key(patient_id) {
                return Err(anyhow::anyhow!(
                    "Patient {} (from project.yaml) not found in patient file",
                    patient_id
                )
                .into());
            }
        }
        return Ok(project_config.patients.clone());
    }

    // No patients specified
    Err(anyhow::anyhow!("No patients specified. Use --patient, --patients, --all, --test flags or define patients in project.yaml").into())
}
