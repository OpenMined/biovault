use crate::cli::download_cache::{
    ChecksumPolicy, ChecksumPolicyType, DownloadCache, DownloadOptions,
};
use crate::cli::syft_url::SyftURL;
use crate::error::Error;
use anyhow::{anyhow, Context};
use colored::Colorize;
use dialoguer::Confirm;
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
// use tempfile::TempDir; // no longer needed since we don't copy templates
use tracing::{debug, info};

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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParticipantData {
    #[serde(default)]
    pub id: String,
    pub ref_version: String,
    #[serde(rename = "ref")]
    pub ref_path: String,
    pub ref_index: String,
    pub aligned: String,
    pub aligned_index: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_b3sum: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_index_b3sum: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_b3sum: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_index_b3sum: Option<String>,
}

pub struct RunParams {
    pub project_folder: String,
    pub participant_source: String,
    pub test: bool,
    pub download: bool,
    pub dry_run: bool,
    pub with_docker: bool,
    pub work_dir: Option<String>,
    pub resume: bool,
}

enum ParticipantSource {
    LocalFile(PathBuf, Option<String>), // path, fragment
    SyftUrl(SyftURL),
    HttpUrl(String),
    SampleDataId(String),
}

impl ParticipantSource {
    fn parse(source: &str) -> anyhow::Result<Self> {
        if source.starts_with("syft://") {
            Ok(ParticipantSource::SyftUrl(SyftURL::parse(source)?))
        } else if source.starts_with("http://") || source.starts_with("https://") {
            Ok(ParticipantSource::HttpUrl(source.to_string()))
        } else {
            // Check if this matches a known sample data ID
            #[derive(serde::Deserialize)]
            struct SampleDataConfig {
                sample_data_urls: std::collections::HashMap<String, serde_yaml::Value>,
            }
            let sample_yaml = include_str!("../../sample_data.yaml");
            if let Ok(cfg) = serde_yaml::from_str::<SampleDataConfig>(sample_yaml) {
                if cfg.sample_data_urls.contains_key(source) {
                    return Ok(ParticipantSource::SampleDataId(source.to_string()));
                }
            }

            // Otherwise treat as local file path
            let (path, fragment) = if let Some(hash_pos) = source.find('#') {
                (
                    source[..hash_pos].to_string(),
                    Some(source[hash_pos + 1..].to_string()),
                )
            } else {
                (source.to_string(), None)
            };
            Ok(ParticipantSource::LocalFile(PathBuf::from(path), fragment))
        }
    }
}

async fn fetch_participant_file(
    source: &ParticipantSource,
) -> anyhow::Result<(String, Option<String>)> {
    match source {
        ParticipantSource::LocalFile(path, fragment) => {
            if !path.exists() {
                return Err(anyhow!("Local file not found: {:?}", path));
            }
            let content = fs::read_to_string(path)
                .with_context(|| format!("Failed to read file: {:?}", path))?;
            Ok((content, fragment.clone()))
        }
        ParticipantSource::SyftUrl(syft_url) => {
            // Convert to HTTP URL and fetch
            // For now, use the first relay server
            let http_url = syft_url.to_http_relay_url("syftbox.net");
            fetch_http_content(&http_url)
                .await
                .map(|content| (content, syft_url.fragment.clone()))
        }
        ParticipantSource::HttpUrl(url) => {
            let (main_url, fragment) = if let Some(hash_pos) = url.find('#') {
                (
                    url[..hash_pos].to_string(),
                    Some(url[hash_pos + 1..].to_string()),
                )
            } else {
                (url.clone(), None)
            };
            fetch_http_content(&main_url)
                .await
                .map(|content| (content, fragment))
        }
        ParticipantSource::SampleDataId(sample_id) => {
            // Ensure sample data is fetched
            crate::cli::commands::sample_data::fetch(Some(vec![sample_id.clone()]), false, true)
                .await?;

            // Load sample_data.yaml to get URLs and compute filenames
            #[derive(serde::Deserialize)]
            struct SampleEntry {
                ref_version: String,
                #[serde(rename = "ref")]
                ref_url: String,
                ref_index: String,
                aligned: serde_yaml::Value,
                aligned_index: String,
            }
            #[derive(serde::Deserialize)]
            struct SampleDataConfig {
                sample_data_urls: std::collections::HashMap<String, SampleEntry>,
            }

            let sample_yaml = include_str!("../../sample_data.yaml");
            let cfg: SampleDataConfig = serde_yaml::from_str(sample_yaml)
                .context("Failed to parse embedded sample data configuration")?;
            let entry = cfg
                .sample_data_urls
                .get(sample_id)
                .ok_or_else(|| anyhow!("Sample data '{}' not found", sample_id))?;

            // Compute local absolute paths under biovault sample data dir
            let biovault_home = crate::config::get_biovault_home()?;
            let sample_data_dir = biovault_home.join("data").join("sample");
            let reference_dir = sample_data_dir.join("reference");
            let participant_dir = sample_data_dir.join(sample_id);

            // Extract filenames from URLs
            fn filename_from_url(url: &str) -> String {
                url.rsplit('/')
                    .next()
                    .unwrap_or("")
                    .split('#')
                    .next()
                    .unwrap()
                    .split('?')
                    .next()
                    .unwrap()
                    .to_string()
            }

            let ref_filename = filename_from_url(&entry.ref_url);
            let ref_index_filename = filename_from_url(&entry.ref_index);

            // Determine aligned file final name
            let aligned_abs_path = match &entry.aligned {
                serde_yaml::Value::String(url) => participant_dir.join(filename_from_url(url)),
                serde_yaml::Value::Sequence(seq) if !seq.is_empty() => {
                    // multiple parts, derive base name from first
                    if let Some(serde_yaml::Value::String(first_url)) = seq.first() {
                        let first_name = filename_from_url(first_url);
                        let base_name = if first_name.ends_with(".tar.gz.aa") {
                            first_name.trim_end_matches(".aa").to_string()
                        } else {
                            first_name
                        };
                        let cram_name = base_name.trim_end_matches(".tar.gz").to_string();
                        participant_dir.join(cram_name)
                    } else {
                        anyhow::bail!("Invalid aligned URL list in sample data");
                    }
                }
                _ => anyhow::bail!("Invalid 'aligned' field in sample data"),
            };

            let aligned_index_abs = if !entry.aligned_index.is_empty() {
                participant_dir.join(filename_from_url(&entry.aligned_index))
            } else {
                PathBuf::new()
            };

            // Build a minimal participants YAML that our existing parser expects
            let mut yaml = String::new();
            yaml.push_str("participants:\n");
            yaml.push_str(&format!("  {}:\n", sample_id));
            yaml.push_str(&format!("    ref_version: {}\n", entry.ref_version));
            yaml.push_str(&format!(
                "    ref: {}\n",
                reference_dir.join(ref_filename).to_string_lossy()
            ));
            yaml.push_str(&format!(
                "    ref_index: {}\n",
                reference_dir.join(ref_index_filename).to_string_lossy()
            ));
            yaml.push_str(&format!(
                "    aligned: {}\n",
                aligned_abs_path.to_string_lossy()
            ));
            if !aligned_index_abs.as_os_str().is_empty() {
                yaml.push_str(&format!(
                    "    aligned_index: {}\n",
                    aligned_index_abs.to_string_lossy()
                ));
            } else {
                yaml.push_str("    aligned_index: \n");
            }

            Ok((yaml, Some(format!("participants.{}", sample_id))))
        }
    }
}

async fn fetch_http_content(url: &str) -> anyhow::Result<String> {
    println!("Fetching participant file from: {}", url.cyan());

    let response = reqwest::get(url)
        .await
        .with_context(|| format!("Failed to fetch URL: {}", url))?;

    if !response.status().is_success() {
        return Err(anyhow!(
            "HTTP request failed with status: {}",
            response.status()
        ));
    }

    response
        .text()
        .await
        .with_context(|| format!("Failed to read response from: {}", url))
}

fn extract_participant_data(
    yaml_content: &str,
    fragment: Option<String>,
    use_mock: bool,
) -> anyhow::Result<(ParticipantData, Option<String>)> {
    let yaml: YamlValue =
        serde_yaml::from_str(yaml_content).with_context(|| "Failed to parse participant YAML")?;

    // Parse fragment to get participant ID
    let participant_id = if let Some(ref frag) = fragment {
        // Expected format: "participants.MADHAVA"
        if frag.starts_with("participants.") {
            frag.strip_prefix("participants.").unwrap().to_string()
        } else {
            return Err(anyhow!(
                "Invalid fragment format. Expected: participants.ID"
            ));
        }
    } else {
        return Err(anyhow!("No participant specified in fragment"));
    };

    // Navigate to the participant
    let participant_yaml = yaml
        .get("participants")
        .and_then(|p| p.get(&participant_id))
        .ok_or_else(|| anyhow!("Participant '{}' not found", participant_id))?;

    if use_mock {
        // Check for mock field
        if let Some(mock_yaml) = participant_yaml.get("mock") {
            // Try to find the mock data key by looking at the ref_version
            // This is a heuristic - we'll use ref_version to determine the mock key
            let mut mock_data: ParticipantData = serde_yaml::from_value(mock_yaml.clone())
                .with_context(|| "Failed to parse mock data")?;
            mock_data.id = participant_id;

            // Determine mock data key based on ref_version
            let mock_key = format!("mock_data_{}", mock_data.ref_version.to_lowercase());

            return Ok((mock_data, Some(mock_key)));
        } else {
            println!(
                "{}",
                "Warning: --test flag set but no mock data available for this participant".yellow()
            );
        }
    }

    // Parse regular participant data
    let mut participant: ParticipantData = serde_yaml::from_value(participant_yaml.clone())
        .with_context(|| format!("Failed to parse participant data for '{}'", participant_id))?;
    participant.id = participant_id;

    Ok((participant, None))
}

async fn ensure_files_exist(
    participant: &ParticipantData,
    auto_download: bool,
    source: &ParticipantSource,
    mock_key: Option<&str>,
) -> anyhow::Result<ParticipantData> {
    let mut local_participant = participant.clone();
    let mut cache = DownloadCache::new(None)?;

    // Get cache directory for checking
    let cache_base = crate::config::get_cache_dir()?;
    let biovault_home = crate::config::get_biovault_home()?;
    let downloads_base = biovault_home.join("data").join("downloads");

    // Create downloads directory based on source
    let participant_downloads_dir = match source {
        ParticipantSource::SyftUrl(syft_url) => {
            // For Syft URLs: downloads/email/participant_id or downloads/email/mock_key
            if let Some(mock_key) = mock_key {
                // For mock data, use the mock key as the directory name
                downloads_base.join(&syft_url.email).join(mock_key)
            } else {
                downloads_base.join(&syft_url.email).join(&participant.id)
            }
        }
        ParticipantSource::HttpUrl(url) => {
            // Try to extract datasite from URL if possible
            if let Ok(parsed_url) = reqwest::Url::parse(url) {
                if let Some(host) = parsed_url.host_str() {
                    if host.contains("syftbox") && url.contains("/datasites/") {
                        // Extract email from URL like https://syftbox.net/datasites/madhava@openmined.org/...
                        if let Some(email_start) = url.find("/datasites/") {
                            let after_datasites = &url[email_start + 11..];
                            if let Some(slash_pos) = after_datasites.find('/') {
                                let email = &after_datasites[..slash_pos];
                                if let Some(mock_key) = mock_key {
                                    downloads_base.join(email).join(mock_key)
                                } else {
                                    downloads_base.join(email).join(&participant.id)
                                }
                            } else {
                                downloads_base.join(&participant.id)
                            }
                        } else {
                            downloads_base.join(&participant.id)
                        }
                    } else {
                        downloads_base.join(&participant.id)
                    }
                } else {
                    downloads_base.join(&participant.id)
                }
            } else {
                downloads_base.join(&participant.id)
            }
        }
        ParticipantSource::LocalFile(path, _) => {
            // For local files, use filename (without extension) and participant_id
            let filename = path.file_stem().and_then(|s| s.to_str()).unwrap_or("local");
            if let Some(mock_key) = mock_key {
                downloads_base.join(filename).join(mock_key)
            } else {
                downloads_base.join(filename).join(&participant.id)
            }
        }
        ParticipantSource::SampleDataId(_) => {
            // Use a dedicated directory for sample data
            downloads_base.join("sample").join(&participant.id)
        }
    };

    fs::create_dir_all(&participant_downloads_dir)?;

    // Helper function to extract filename from URL
    fn extract_filename(url: &str) -> String {
        url.split('/').next_back().unwrap_or("unknown").to_string()
    }

    // List of files to check/download
    let files_to_check = vec![
        (
            "reference",
            participant.ref_path.clone(),
            participant.ref_b3sum.clone(),
        ),
        (
            "reference index",
            participant.ref_index.clone(),
            participant.ref_index_b3sum.clone(),
        ),
        (
            "aligned",
            participant.aligned.clone(),
            participant.aligned_b3sum.clone(),
        ),
        (
            "aligned index",
            participant.aligned_index.clone(),
            participant.aligned_index_b3sum.clone(),
        ),
    ];

    let mut downloads_needed = Vec::new();
    let mut cached_paths: HashMap<String, (PathBuf, String)> = HashMap::new();

    // First check what needs downloading
    for (name, url, b3sum) in &files_to_check {
        // Check if it's a URL
        if url.starts_with("http://") || url.starts_with("https://") {
            // Check cache first if we have a checksum
            if let Some(checksum) = b3sum {
                let cache_path = cache_base.join("by-hash").join(checksum);
                if cache_path.exists() {
                    debug!("Found {} in cache: {:?}", name, cache_path);
                    let filename = extract_filename(url);
                    cached_paths.insert(name.to_string(), (cache_path, filename));
                    continue;
                }
            }
            downloads_needed.push((name.to_string(), url.clone(), b3sum.clone()));
        } else {
            // Local file - check existence and keep path as-is
            if !Path::new(url).exists() {
                return Err(anyhow!("Local file not found: {} at {}", name, url));
            }
            // Keep local paths unchanged
            match *name {
                "reference" => local_participant.ref_path = url.clone(),
                "reference index" => local_participant.ref_index = url.clone(),
                "aligned" => local_participant.aligned = url.clone(),
                "aligned index" => local_participant.aligned_index = url.clone(),
                _ => {}
            }
        }
    }

    // Create symlinks for cached files with proper filenames
    if !cached_paths.is_empty() {
        println!("Using cached files with proper filenames:");
    }

    for (name, (cache_path, filename)) in cached_paths {
        let symlink_path = participant_downloads_dir.join(&filename);

        // Remove existing symlink if it exists
        if symlink_path.exists() || symlink_path.is_symlink() {
            fs::remove_file(&symlink_path).ok();
        }

        // Create symlink to cache
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(&cache_path, &symlink_path)
                .with_context(|| format!("Failed to create symlink for {}", name))?;
        }
        #[cfg(windows)]
        {
            std::os::windows::fs::symlink_file(&cache_path, &symlink_path)
                .with_context(|| format!("Failed to create symlink for {}", name))?;
        }

        println!("  • {} → {}", name.green(), symlink_path.display());

        // Update participant paths with symlink paths
        match name.as_str() {
            "reference" => local_participant.ref_path = symlink_path.to_string_lossy().to_string(),
            "reference index" => {
                local_participant.ref_index = symlink_path.to_string_lossy().to_string()
            }
            "aligned" => local_participant.aligned = symlink_path.to_string_lossy().to_string(),
            "aligned index" => {
                local_participant.aligned_index = symlink_path.to_string_lossy().to_string()
            }
            _ => {}
        }
    }

    if !downloads_needed.is_empty() {
        println!("\nThe following files need to be downloaded:");
        for (name, url, _) in &downloads_needed {
            println!("  - {} from {}", name.cyan(), url);
        }

        let should_download = if auto_download {
            true
        } else {
            Confirm::new()
                .with_prompt("Do you want to download these files?")
                .default(true)
                .interact()?
        };

        if !should_download {
            return Err(anyhow!("File downloads cancelled by user"));
        }

        // Download files and update paths
        for (name, url, b3sum) in &downloads_needed {
            println!("Downloading {}...", name.green());

            // Extract filename from URL
            let filename = extract_filename(url);
            let symlink_path = participant_downloads_dir.join(&filename);

            // Create a temporary path for download
            let temp_dir = tempfile::tempdir()?;
            let temp_path = temp_dir.path().join(&filename);

            // Set up download options
            let checksum_policy = if let Some(hash) = b3sum {
                ChecksumPolicy {
                    policy_type: ChecksumPolicyType::Required,
                    expected_hash: Some(hash.clone()),
                }
            } else {
                ChecksumPolicy {
                    policy_type: ChecksumPolicyType::Optional,
                    expected_hash: None,
                }
            };

            let options = DownloadOptions {
                checksum_policy,
                show_progress: true,
                cache_strategy: Default::default(),
            };

            // Download (will be stored in cache)
            let _downloaded_path = cache.download_with_cache(url, &temp_path, options).await?;

            // After download, the file is in cache. Create symlink with proper filename
            if let Some(hash) = b3sum {
                let cache_path = crate::config::get_cache_dir()?.join("by-hash").join(hash);

                // Remove existing symlink if it exists
                if symlink_path.exists() || symlink_path.is_symlink() {
                    fs::remove_file(&symlink_path).ok();
                }

                // Create symlink to cache
                #[cfg(unix)]
                {
                    std::os::unix::fs::symlink(&cache_path, &symlink_path)
                        .with_context(|| format!("Failed to create symlink for {}", name))?;
                }
                #[cfg(windows)]
                {
                    std::os::windows::fs::symlink_file(&cache_path, &symlink_path)
                        .with_context(|| format!("Failed to create symlink for {}", name))?;
                }

                println!("  • {} → {}", name.green(), symlink_path.display());
            } else {
                // If no checksum, copy the temp file to downloads directory
                fs::copy(&temp_path, &symlink_path)?;
                println!("  • {} → {}", name.green(), symlink_path.display());
            }

            // Update the participant data with symlink paths
            match name.as_str() {
                "reference" => {
                    local_participant.ref_path = symlink_path.to_string_lossy().to_string()
                }
                "reference index" => {
                    local_participant.ref_index = symlink_path.to_string_lossy().to_string()
                }
                "aligned" => local_participant.aligned = symlink_path.to_string_lossy().to_string(),
                "aligned index" => {
                    local_participant.aligned_index = symlink_path.to_string_lossy().to_string()
                }
                _ => {}
            }
        }
    }

    Ok(local_participant)
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
        .context("Failed to resolve workflow.nf path")?;

    if !workflow_file.exists() {
        return Err(Error::WorkflowMissing(params.project_folder.clone()).into());
    }

    // Read project config
    let config_content =
        fs::read_to_string(&project_yaml).context("Failed to read project.yaml")?;
    let config: ProjectConfig =
        serde_yaml::from_str(&config_content).context("Failed to parse project.yaml")?;

    // Parse participant source
    let source = ParticipantSource::parse(&params.participant_source)?;

    // Fetch participant file
    let (yaml_content, fragment) = fetch_participant_file(&source).await?;

    // Extract participant data
    let (mut participant, mock_key) =
        extract_participant_data(&yaml_content, fragment, params.test)?;

    // Ensure all required files exist (download if needed)
    participant =
        ensure_files_exist(&participant, params.download, &source, mock_key.as_deref()).await?;

    // Get BioVault environment directory
    let biovault_home = crate::config::get_biovault_home()?;
    let env_dir = biovault_home.join("env").join("default");

    // Check if templates exist
    let template_nf = env_dir.join("template.nf");
    let nextflow_config = env_dir.join("nextflow.config");

    if !template_nf.exists() || !nextflow_config.exists() {
        return Err(Error::TemplatesNotFound.into());
    }

    // Use templates directly from env dir instead of copying to temp
    let temp_template = template_nf;
    let temp_config = nextflow_config;

    // Determine assets directory
    let assets_dir = if !config.assets.is_empty() {
        project_path.join(&config.assets[0])
    } else {
        project_path.join("assets")
    };

    // Create assets directory if it doesn't exist
    if !assets_dir.exists() {
        fs::create_dir_all(&assets_dir)?;
    }

    // Get absolute path for assets directory
    let assets_dir = assets_dir
        .canonicalize()
        .context("Failed to resolve assets directory path")?;

    // Create results directory for this participant
    let results_dir = project_path.join("results").join(&participant.id);
    if !results_dir.exists() {
        fs::create_dir_all(&results_dir)?;
    }

    // Get absolute path for results directory
    let results_dir = results_dir
        .canonicalize()
        .context("Failed to resolve results directory path")?;

    info!(
        "Running workflow '{}' from project '{}'",
        config.workflow, config.name
    );

    println!("Processing participant: {}", participant.id.cyan());

    // Build Nextflow command
    let mut cmd = Command::new("nextflow");

    // Set working directory to project directory
    cmd.current_dir(&project_path);

    cmd.arg("run")
        .arg(&temp_template)
        .arg("--participant_id")
        .arg(&participant.id)
        .arg("--ref_version")
        .arg(&participant.ref_version)
        .arg("--ref")
        .arg(&participant.ref_path)
        .arg("--ref_index")
        .arg(&participant.ref_index)
        .arg("--aligned")
        .arg(&participant.aligned)
        .arg("--aligned_index")
        .arg(&participant.aligned_index)
        .arg("--work_flow_file")
        .arg(workflow_file.to_string_lossy().as_ref())
        .arg("--assets_dir")
        .arg(assets_dir.to_string_lossy().as_ref())
        .arg("--results_dir")
        .arg(results_dir.to_string_lossy().as_ref());

    if params.resume {
        cmd.arg("-resume");
    }

    if let Some(work_dir) = params.work_dir {
        cmd.arg("-work-dir");
        cmd.arg(work_dir);
    }

    // Docker/Singularity configuration
    if params.with_docker {
        cmd.arg("-with-docker");
    }

    // Add config file
    cmd.arg("-c").arg(&temp_config);

    // Print the command that will be executed
    println!("\n{}", "Nextflow command:".green().bold());

    // Build command string for display
    let mut cmd_str = String::from("nextflow");
    for arg in cmd.get_args() {
        cmd_str.push(' ');
        let arg_str = arg.to_string_lossy();
        // Quote arguments with spaces
        if arg_str.contains(' ') {
            cmd_str.push_str(&format!("'{}'", arg_str));
        } else {
            cmd_str.push_str(&arg_str);
        }
    }
    println!("{}\n", cmd_str.cyan());

    if params.dry_run {
        println!("{}", "[DRY RUN] Would execute the above command".yellow());
        return Ok(());
    }

    // Execute Nextflow
    println!("Executing Nextflow workflow...");
    let status = cmd.status().context("Failed to execute Nextflow")?;

    if !status.success() {
        return Err(anyhow!("Nextflow execution failed"));
    }

    println!("{}", "Workflow completed successfully!".green().bold());
    Ok(())
}
