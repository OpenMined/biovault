use crate::error::Error;
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
    participants: Vec<String>,
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
struct ParticipantData {
    ref_version: String,
    #[serde(rename = "ref")]
    ref_path: String,
    ref_index: String,
    aligned: String,
    aligned_index: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct ParticipantFile {
    participant: std::collections::HashMap<String, ParticipantData>,
}

pub struct RunParams {
    pub project_folder: String,
    pub participant_file: String,
    pub participants: Option<Vec<String>>,
    pub participant: Option<String>,
    pub all: bool,
    pub test: bool,
    pub dry_run: bool,
    pub with_docker: bool,
    pub work_dir: Option<String>,
    pub resume: bool,
}

pub async fn execute(params: RunParams) -> anyhow::Result<()> {
    // Validate project directory
    let project_path = PathBuf::from(&params.project_folder);
    if !project_path.exists() {
        return Err(Error::ProjectFolderMissing(params.project_folder.clone()).into());
    }

    let project_yaml = project_path.join("project.yaml");
    if !project_yaml.exists() {
        return Err(Error::ProjectConfigMissing(params.project_folder.clone()).into());
    }

    let workflow_file = project_path
        .join("workflow.nf")
        .canonicalize()
        .with_context(|| {
            format!(
                "Failed to resolve workflow.nf path in {}",
                params.project_folder
            )
        })?;
    if !workflow_file.exists() {
        return Err(Error::WorkflowMissing(params.project_folder.clone()).into());
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

    // Parse participant file
    let participant_file_path = PathBuf::from(&params.participant_file);
    if !participant_file_path.exists() {
        return Err(Error::ParticipantFileMissing(params.participant_file.clone()).into());
    }

    // Get the directory containing the participant file for resolving relative paths
    let participant_dir = participant_file_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Could not determine participant file directory"))?;

    let participant_content = fs::read_to_string(&participant_file_path).with_context(|| {
        format!(
            "Failed to read participant file from {}",
            participant_file_path.display()
        )
    })?;
    let participant_data: ParticipantFile =
        serde_yaml::from_str(&participant_content).context("Failed to parse participant file")?;

    // Determine which participants to process
    let participants_to_run = determine_participants(
        &participant_data,
        &project_config,
        params.participants,
        params.participant,
        params.all,
        params.test,
    )?;

    if participants_to_run.is_empty() {
        println!("No participants to process");
        return Ok(());
    }

    println!("Processing {} participant(s):", participants_to_run.len());
    for p in &participants_to_run {
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
        return Err(Error::TemplatesNotFound.into());
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

    // Process each participant
    let mut success_count = 0;
    let mut fail_count = 0;

    for participant_id in participants_to_run {
        println!("\n{}", "=".repeat(60));
        println!("Processing participant: {}", participant_id);
        println!("{}", "=".repeat(60));

        // Get participant data
        let participant_info = participant_data
            .participant
            .get(&participant_id)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Participant {} not found in participant file",
                    participant_id
                )
            })?;

        // Create results directory for this participant
        let results_dir = project_path.join("results").join(&participant_id);
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

        // Resolve participant data paths relative to participant file location
        let ref_path = participant_dir.join(&participant_info.ref_path);
        let ref_index_path = participant_dir.join(&participant_info.ref_index);
        let aligned_path = participant_dir.join(&participant_info.aligned);
        let aligned_index_path = participant_dir.join(&participant_info.aligned_index);

        // Verify that resolved paths exist and canonicalize them
        if !ref_path.exists() {
            return Err(Error::FileNotFound {
                file: participant_info.ref_path.clone(),
                details: format!("resolved to {}", ref_path.display()),
            }
            .into());
        }
        let ref_path = ref_path.canonicalize().with_context(|| {
            format!(
                "Failed to resolve reference file path: {}",
                ref_path.display()
            )
        })?;

        if !ref_index_path.exists() {
            return Err(Error::FileNotFound {
                file: participant_info.ref_index.clone(),
                details: format!("resolved to {}", ref_index_path.display()),
            }
            .into());
        }
        let ref_index_path = ref_index_path.canonicalize().with_context(|| {
            format!(
                "Failed to resolve reference index path: {}",
                ref_index_path.display()
            )
        })?;

        if !aligned_path.exists() {
            return Err(Error::FileNotFound {
                file: participant_info.aligned.clone(),
                details: format!("resolved to {}", aligned_path.display()),
            }
            .into());
        }
        let aligned_path = aligned_path.canonicalize().with_context(|| {
            format!(
                "Failed to resolve aligned file path: {}",
                aligned_path.display()
            )
        })?;

        if !aligned_index_path.exists() {
            return Err(Error::FileNotFound {
                file: participant_info.aligned_index.clone(),
                details: format!("resolved to {}", aligned_index_path.display()),
            }
            .into());
        }
        let aligned_index_path = aligned_index_path.canonicalize().with_context(|| {
            format!(
                "Failed to resolve aligned index path: {}",
                aligned_index_path.display()
            )
        })?;

        // Build nextflow command
        let mut cmd = Command::new("nextflow");

        // Set working directory to project directory instead of temp directory
        cmd.current_dir(&project_path);

        cmd.arg("run")
            .arg(&temp_template)
            .arg("--participant_id")
            .arg(&participant_id)
            .arg("--ref_version")
            .arg(&participant_info.ref_version)
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
        if let Some(ref work_dir_path) = params.work_dir {
            cmd.arg("-work-dir").arg(work_dir_path);
        }

        // Add resume flag if specified
        if params.resume {
            cmd.arg("-resume");
        }

        // Add Docker flag
        if params.with_docker {
            cmd.arg("-with-docker");
        }

        // Add config
        cmd.arg("-c").arg(&temp_config);

        // Show command if dry-run or normal execution
        let cmd_str = format!("{:?}", cmd);
        println!("\nCommand: {}", cmd_str);

        if params.dry_run {
            println!("[DRY RUN] Would execute the above command");
            success_count += 1;
        } else {
            println!("\nExecuting nextflow pipeline...\n");

            // Execute the command
            let output = cmd.output().context("Failed to execute nextflow command")?;

            if output.status.success() {
                println!("✓ Successfully processed participant: {}", participant_id);
                success_count += 1;

                // Print stdout if available
                if !output.stdout.is_empty() {
                    println!("\nNextflow output:");
                    println!("{}", String::from_utf8_lossy(&output.stdout));
                }
            } else {
                println!("✗ Failed to process participant: {}", participant_id);
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
        return Err(Error::ProcessingFailed(fail_count).into());
    }

    Ok(())
}

fn determine_participants(
    participant_data: &ParticipantFile,
    project_config: &ProjectConfig,
    participants_arg: Option<Vec<String>>,
    participant_arg: Option<String>,
    all: bool,
    test: bool,
) -> anyhow::Result<Vec<String>> {
    // Priority 1: Test mode
    if test {
        if !participant_data.participant.contains_key("TEST") {
            return Err(Error::ParticipantNotFoundInFile("TEST".to_string()).into());
        }
        return Ok(vec!["TEST".to_string()]);
    }

    // Priority 2: Command-line participant override (single)
    if let Some(participant_id) = participant_arg {
        if !participant_data.participant.contains_key(&participant_id) {
            return Err(Error::ParticipantNotFoundInFile(participant_id.clone()).into());
        }
        return Ok(vec![participant_id]);
    }

    // Priority 3: Command-line participants override (multiple)
    if let Some(participant_list) = participants_arg {
        for participant_id in &participant_list {
            if !participant_data.participant.contains_key(participant_id) {
                return Err(Error::ParticipantNotFoundInFile(participant_id.clone()).into());
            }
        }
        return Ok(participant_list);
    }

    // Priority 4: All participants
    if all {
        return Ok(participant_data.participant.keys().cloned().collect());
    }

    // Priority 5: Project-defined participants
    if !project_config.participants.is_empty() {
        for participant_id in &project_config.participants {
            if !participant_data.participant.contains_key(participant_id) {
                return Err(Error::ParticipantNotFoundInFile(format!(
                    "{} (from project.yaml)",
                    participant_id
                ))
                .into());
            }
        }
        return Ok(project_config.participants.clone());
    }

    // No participants specified
    Err(Error::NoParticipantsSpecified.into())
}
