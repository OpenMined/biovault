use crate::cli::download_cache::{
    ChecksumPolicy, ChecksumPolicyType, DownloadCache, DownloadOptions,
};
use crate::cli::syft_url::SyftURL;
use crate::error::Error;
use anyhow::{anyhow, Context};
use chrono::Local;
use colored::Colorize;
use dialoguer::Confirm;
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::HashMap;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, BufReader, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Duration;
// use tempfile::TempDir; // no longer needed since we don't copy templates
use tracing::{debug, info};

#[derive(Debug, Serialize, Deserialize)]
struct ProjectConfig {
    name: String,
    author: String,
    workflow: String,
    #[serde(default)]
    template: Option<String>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ref_version: Option<String>,
    #[serde(rename = "ref", default, skip_serializing_if = "Option::is_none")]
    pub ref_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ref_index: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned_index: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ref_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ref_index_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned_index_b3sum: Option<String>,
    // SNP data fields
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub snp: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub snp_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uncompress: Option<bool>,
}

fn append_desktop_log_to_path(path: &str, level: &str, message: &str) -> std::io::Result<()> {
    let timestamp = Local::now().format("%Y-%m-%dT%H:%M:%S%:z");
    let line = format!("[{}][{}] {}\n", timestamp, level, message);

    if let Some(parent) = std::path::Path::new(path).parent() {
        fs::create_dir_all(parent)?;
    }

    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    file.write_all(line.as_bytes())
}

fn append_desktop_log(level: &str, message: &str) {
    if let Ok(path) = std::env::var("BIOVAULT_DESKTOP_LOG_FILE") {
        if path.is_empty() {
            return;
        }
        if let Err(err) = append_desktop_log_to_path(&path, level, message) {
            debug!("Failed to append to desktop log {}: {}", path, err);
        }
    }
}

fn sanitize_stream_line(line: &str) -> String {
    let segment = line.rsplit('\r').next().unwrap_or(line);
    strip_ansi_sequences(segment)
}

fn strip_ansi_sequences(input: &str) -> String {
    let mut output = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();

    while let Some(ch) = chars.next() {
        if ch == '\u{001b}' {
            match chars.peek().copied() {
                Some('[') => {
                    chars.next();
                    for c in chars.by_ref() {
                        if ('@'..='~').contains(&c) {
                            break;
                        }
                    }
                }
                Some(']') => {
                    chars.next();
                    while let Some(c) = chars.next() {
                        if c == '\u{0007}' {
                            break;
                        }
                        if c == '\u{001b}' {
                            if matches!(chars.peek(), Some('\\')) {
                                chars.next();
                            }
                            break;
                        }
                    }
                }
                Some(_) => {}
                None => {}
            }
            continue;
        }
        output.push(ch);
    }

    output
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

struct LogTailHandle {
    stop_flag: Arc<AtomicBool>,
    handle: thread::JoinHandle<()>,
}

impl LogTailHandle {
    fn stop(self) {
        self.stop_flag.store(true, Ordering::SeqCst);
        let _ = self.handle.join();
    }
}

fn spawn_nextflow_log_forwarder(log_path: PathBuf) -> Option<LogTailHandle> {
    match std::env::var("BIOVAULT_DESKTOP_LOG_FILE") {
        Ok(path) if !path.is_empty() => path,
        _ => return None,
    };

    let stop_flag = Arc::new(AtomicBool::new(false));
    let thread_flag = stop_flag.clone();

    let handle = thread::spawn(move || {
        let mut offset: u64 = 0;
        while !thread_flag.load(Ordering::SeqCst) {
            forward_nextflow_lines(&log_path, &mut offset);
            thread::sleep(Duration::from_millis(500));
        }
        forward_nextflow_lines(&log_path, &mut offset);
    });

    Some(LogTailHandle { stop_flag, handle })
}

fn forward_nextflow_lines(log_path: &Path, offset: &mut u64) {
    if let Ok(metadata) = fs::metadata(log_path) {
        if metadata.len() < *offset {
            *offset = 0;
        }
    }

    let mut file = match std::fs::File::open(log_path) {
        Ok(file) => file,
        Err(_) => return, // File will appear once Nextflow starts writing.
    };

    if file.seek(SeekFrom::Start(*offset)).is_err() {
        return;
    }

    let mut reader = BufReader::new(file);
    let mut line = String::new();

    loop {
        line.clear();
        match reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(n) => {
                *offset += n as u64;
                let trimmed = line.trim_end_matches(['\r', '\n']);
                if !should_forward_nextflow_line(trimmed) {
                    continue;
                }
                append_desktop_log("INFO", &format!("[Nextflow log] {}", trimmed));
            }
            Err(_) => break,
        }
    }
}

fn should_forward_nextflow_line(line: &str) -> bool {
    if line.is_empty() {
        return false;
    }
    let info = "] INFO ";
    let warn = "] WARN ";
    let error = "] ERROR ";
    (line.contains(info) || line.contains(warn) || line.contains(error))
        && !line.contains("] DEBUG ")
}

pub(crate) fn execute_with_logging(
    mut cmd: Command,
    nextflow_log_path: Option<PathBuf>,
) -> anyhow::Result<std::process::ExitStatus> {
    cmd.stdout(Stdio::piped());
    cmd.stderr(Stdio::piped());

    let mut child = cmd.spawn().context("Failed to execute Nextflow")?;

    let log_path = std::env::var("BIOVAULT_DESKTOP_LOG_FILE").ok();

    let stdout_handle = child.stdout.take().map(|stdout| {
        let log_path = log_path.clone();
        thread::spawn(move || {
            let reader = BufReader::new(stdout);
            let mut last_line: Option<String> = None;
            for line in reader.lines().map_while(Result::ok) {
                println!("{}", line);
                let sanitized = sanitize_stream_line(&line);
                if sanitized.is_empty() {
                    continue;
                }
                if last_line.as_deref() == Some(&sanitized) {
                    continue;
                }
                last_line = Some(sanitized.clone());
                if let Some(ref path) = log_path {
                    let _ = append_desktop_log_to_path(path, "INFO", &sanitized);
                } else {
                    append_desktop_log("INFO", &sanitized);
                }
            }
        })
    });

    let stderr_handle = child.stderr.take().map(|stderr| {
        let log_path = log_path.clone();
        thread::spawn(move || {
            let reader = BufReader::new(stderr);
            let mut last_line: Option<String> = None;
            for line in reader.lines().map_while(Result::ok) {
                eprintln!("{}", line);
                let sanitized = sanitize_stream_line(&line);
                if sanitized.is_empty() {
                    continue;
                }
                if last_line.as_deref() == Some(&sanitized) {
                    continue;
                }
                last_line = Some(sanitized.clone());
                if let Some(ref path) = log_path {
                    let _ = append_desktop_log_to_path(path, "ERROR", &sanitized);
                } else {
                    append_desktop_log("ERROR", &sanitized);
                }
            }
        })
    });

    let log_forwarder = nextflow_log_path.and_then(spawn_nextflow_log_forwarder);

    let status = child
        .wait()
        .context("Failed to wait for Nextflow process")?;

    if let Some(handle) = stdout_handle {
        let _ = handle.join();
    }
    if let Some(handle) = stderr_handle {
        let _ = handle.join();
    }
    if let Some(forwarder) = log_forwarder {
        forwarder.stop();
    }

    Ok(status)
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
    pub template: Option<String>,
    pub results_dir: Option<String>,
    pub nextflow_args: Vec<String>,
}

enum ParticipantSource {
    LocalFile(PathBuf, Option<String>), // path, fragment
    SyftUrl(SyftURL),
    HttpUrl(String),
    SampleDataId(String),
    RegisteredParticipant(String), // participant ID from BioVault home participants.yaml
}

impl ParticipantSource {
    fn parse(source: &str) -> anyhow::Result<Self> {
        if source.starts_with("syft://") {
            Ok(ParticipantSource::SyftUrl(SyftURL::parse(source)?))
        } else if source.starts_with("http://") || source.starts_with("https://") {
            Ok(ParticipantSource::HttpUrl(source.to_string()))
        } else {
            // First check if this is a registered participant (no path separators, no fragment)
            if !source.contains('/') && !source.contains('#') {
                // Try to load registered participants
                if let Ok(participants_path) = crate::config::get_biovault_home() {
                    let participants_file = participants_path.join("participants.yaml");
                    if participants_file.exists() {
                        if let Ok(contents) = fs::read_to_string(&participants_file) {
                            #[derive(serde::Deserialize)]
                            struct ParticipantsFile {
                                participants: std::collections::HashMap<String, serde_yaml::Value>,
                            }
                            if let Ok(parsed) = serde_yaml::from_str::<ParticipantsFile>(&contents)
                            {
                                if parsed.participants.contains_key(source) {
                                    return Ok(ParticipantSource::RegisteredParticipant(
                                        source.to_string(),
                                    ));
                                }
                            }
                        }
                    }
                }
            }

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
    auto_download: bool,
) -> anyhow::Result<(String, Option<String>)> {
    match source {
        ParticipantSource::RegisteredParticipant(participant_id) => {
            // Load the participants file
            let participants_path = crate::config::get_biovault_home()?.join("participants.yaml");
            if !participants_path.exists() {
                return Err(anyhow!("No registered participants found"));
            }
            let content = fs::read_to_string(&participants_path).with_context(|| {
                format!("Failed to read participants file: {:?}", participants_path)
            })?;
            // Return the full file content with the participant ID as fragment
            Ok((content, Some(format!("participants.{}", participant_id))))
        }
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
            // Check if sample data exists locally first
            let biovault_home = crate::config::get_biovault_home()?;
            let participants_file = biovault_home
                .join("data")
                .join("sample")
                .join("participants.yaml");

            // If participants file doesn't exist or auto_download is false, we need to prompt
            if !participants_file.exists() && !auto_download {
                println!("Sample data for '{}' needs to be downloaded.", sample_id);
                println!("This may take some time depending on the file size.");

                let proceed = dialoguer::Confirm::new()
                    .with_prompt("Do you want to download the sample data now?")
                    .default(true)
                    .interact()?;

                if !proceed {
                    return Err(anyhow!("Sample data download cancelled by user"));
                }
            }

            // Fetch sample data (with quiet=false so user sees progress)
            crate::cli::commands::sample_data::fetch(Some(vec![sample_id.clone()]), false, false)
                .await?;

            // Load sample_data.yaml to get URLs and compute filenames
            #[derive(serde::Deserialize)]
            struct PostProcess {
                #[serde(default)]
                #[allow(dead_code)]
                uncompress: Option<bool>,
                #[serde(default)]
                file: Option<String>,
            }
            #[derive(serde::Deserialize)]
            struct SampleEntry {
                #[serde(default)]
                ref_version: Option<String>,
                #[serde(rename = "ref", default)]
                ref_url: Option<String>,
                #[serde(default)]
                ref_index: Option<String>,
                #[serde(default)]
                aligned: Option<serde_yaml::Value>,
                #[serde(default)]
                aligned_index: Option<String>,
                // SNP fields
                #[serde(default)]
                snp: Option<String>,
                #[serde(default)]
                #[allow(dead_code)]
                snp_b3sum: Option<String>,
                #[serde(default)]
                snp_post_process: Option<PostProcess>,
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

            let ref_filename = entry
                .ref_url
                .as_ref()
                .map(|url| filename_from_url(url))
                .unwrap_or_default();
            let ref_index_filename = entry
                .ref_index
                .as_ref()
                .map(|url| filename_from_url(url))
                .unwrap_or_default();

            // Determine aligned file final name
            let aligned_abs_path = match entry.aligned.as_ref() {
                Some(serde_yaml::Value::String(url)) => {
                    participant_dir.join(filename_from_url(url))
                }
                Some(serde_yaml::Value::Sequence(seq)) if !seq.is_empty() => {
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
                None => PathBuf::new(), // No aligned field for SNP data
                _ => anyhow::bail!("Invalid 'aligned' field in sample data"),
            };

            let aligned_index_abs = if let Some(aligned_index) = entry.aligned_index.as_ref() {
                if !aligned_index.is_empty() {
                    participant_dir.join(filename_from_url(aligned_index))
                } else {
                    PathBuf::new()
                }
            } else {
                PathBuf::new()
            };

            // Build a minimal participants YAML that our existing parser expects
            let mut yaml = String::new();
            yaml.push_str("participants:\n");
            yaml.push_str(&format!("  {}:\n", sample_id));
            if let Some(ref_version) = &entry.ref_version {
                yaml.push_str(&format!("    ref_version: {}\n", ref_version));
            }

            // Check if this is SNP data or CRAM data
            if let Some(snp) = &entry.snp {
                // SNP data - point to the specific file if specified in post_process
                let snp_path = if let Some(ref post_process) = entry.snp_post_process {
                    if let Some(ref file) = post_process.file {
                        participant_dir.join(file)
                    } else {
                        // Strip .zip if present
                        participant_dir.join(filename_from_url(snp).replace(".zip", ""))
                    }
                } else {
                    // Strip .zip if present
                    participant_dir.join(filename_from_url(snp).replace(".zip", ""))
                };
                yaml.push_str(&format!("    snp: {}\n", snp_path.to_string_lossy()));
            } else {
                // CRAM data - include reference and alignment files
                if !ref_filename.is_empty() {
                    yaml.push_str(&format!(
                        "    ref: {}\n",
                        reference_dir.join(ref_filename).to_string_lossy()
                    ));
                }
                if !ref_index_filename.is_empty() {
                    yaml.push_str(&format!(
                        "    ref_index: {}\n",
                        reference_dir.join(ref_index_filename).to_string_lossy()
                    ));
                }
                if !aligned_abs_path.as_os_str().is_empty() {
                    yaml.push_str(&format!(
                        "    aligned: {}\n",
                        aligned_abs_path.to_string_lossy()
                    ));
                }
                if !aligned_index_abs.as_os_str().is_empty() {
                    yaml.push_str(&format!(
                        "    aligned_index: {}\n",
                        aligned_index_abs.to_string_lossy()
                    ));
                }
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
            let mock_key = format!(
                "mock_data_{}",
                mock_data
                    .ref_version
                    .as_ref()
                    .unwrap_or(&"unknown".to_string())
                    .to_lowercase()
            );

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
        ParticipantSource::RegisteredParticipant(_) => {
            // For registered participants, use "registered" subdirectory
            downloads_base.join("registered").join(&participant.id)
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
        if let Some(url_str) = url {
            if url_str.starts_with("http://") || url_str.starts_with("https://") {
                // Check cache first if we have a checksum
                if let Some(checksum) = b3sum {
                    let cache_path = cache_base.join("by-hash").join(checksum);
                    if cache_path.exists() {
                        debug!("Found {} in cache: {:?}", name, cache_path);
                        let filename = extract_filename(url_str);
                        cached_paths.insert(name.to_string(), (cache_path, filename));
                        continue;
                    }
                }
                downloads_needed.push((name.to_string(), url.clone(), b3sum.clone()));
            } else {
                // Local file - check existence and keep path as-is
                if !Path::new(url_str).exists() {
                    return Err(anyhow!("Local file not found: {} at {}", name, url_str));
                }
                // Keep local paths unchanged
                match *name {
                    "reference" => local_participant.ref_path = Some(url_str.to_string()),
                    "reference index" => local_participant.ref_index = Some(url_str.to_string()),
                    "aligned" => local_participant.aligned = Some(url_str.to_string()),
                    "aligned index" => local_participant.aligned_index = Some(url_str.to_string()),
                    "snp" => local_participant.snp = Some(url_str.to_string()),
                    _ => {}
                }
            }
        } // Close the if let Some(url_str) = url block
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
            "reference" => {
                local_participant.ref_path = Some(symlink_path.to_string_lossy().to_string())
            }
            "reference index" => {
                local_participant.ref_index = Some(symlink_path.to_string_lossy().to_string())
            }
            "aligned" => {
                local_participant.aligned = Some(symlink_path.to_string_lossy().to_string())
            }
            "aligned index" => {
                local_participant.aligned_index = Some(symlink_path.to_string_lossy().to_string())
            }
            _ => {}
        }
    }

    if !downloads_needed.is_empty() {
        println!("\nThe following files need to be downloaded:");
        for (name, url, _) in &downloads_needed {
            println!(
                "  - {} from {}",
                name.cyan(),
                url.as_ref().unwrap_or(&"unknown".to_string())
            );
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

            let url_str = url
                .as_ref()
                .ok_or_else(|| anyhow!("Missing URL for {}", name))?;

            // Extract filename from URL
            let filename = extract_filename(url_str);
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
            let _downloaded_path = cache
                .download_with_cache(url_str, &temp_path, options)
                .await?;

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
                    local_participant.ref_path = Some(symlink_path.to_string_lossy().to_string())
                }
                "reference index" => {
                    local_participant.ref_index = Some(symlink_path.to_string_lossy().to_string())
                }
                "aligned" => {
                    local_participant.aligned = Some(symlink_path.to_string_lossy().to_string())
                }
                "aligned index" => {
                    local_participant.aligned_index =
                        Some(symlink_path.to_string_lossy().to_string())
                }
                _ => {}
            }
        }
    }

    Ok(local_participant)
}

async fn execute_sheet_workflow(params: &RunParams, config: &ProjectConfig) -> anyhow::Result<()> {
    let project_path = PathBuf::from(&params.project_folder);

    println!("Running sheet-based workflow: {}", config.name.cyan());

    // For sheet template, participant_source is the required path to the samplesheet
    if params.participant_source.is_empty() {
        return Err(anyhow!(
            "Sheet template requires a CSV/TSV file path. Usage: bv run {} <path/to/samplesheet.csv>",
            params.project_folder
        ));
    }

    let samplesheet_path = PathBuf::from(&params.participant_source);
    if !samplesheet_path.exists() {
        return Err(anyhow!(
            "Samplesheet file not found: {}",
            samplesheet_path.display()
        ));
    }

    // Determine assets directory and look for schema.yaml there
    let mut assets_dir = project_path.join("assets");
    if !config.assets.is_empty() {
        let candidate = project_path.join(&config.assets[0]);
        if candidate.is_dir() {
            assets_dir = candidate;
        }
    }

    // Look for schema.yaml in assets directory (preferred) or project root
    let schema_path = if config.assets.contains(&"schema.yaml".to_string()) {
        // schema.yaml is listed in assets, find it
        let schema_in_assets = assets_dir.join("schema.yaml");
        if schema_in_assets.exists() {
            schema_in_assets
        } else {
            project_path.join("schema.yaml")
        }
    } else {
        project_path.join("schema.yaml")
    };

    if !schema_path.exists() {
        println!(
            "{}",
            "Warning: schema.yaml not found. Using defaults.".yellow()
        );
    }

    // Get BioVault environment directory
    let biovault_home = crate::config::get_biovault_home()?;
    let template_name = config
        .template
        .clone()
        .unwrap_or_else(|| "sheet".to_string());
    let env_dir = biovault_home.join("env").join(&template_name);

    // Check if templates exist
    let template_nf = env_dir.join("template.nf");
    let nextflow_config = env_dir.join("nextflow.config");

    if !template_nf.exists() || !nextflow_config.exists() {
        return Err(Error::TemplatesNotFound.into());
    }

    println!("Using sheet template from: {}", env_dir.display());

    // Convert to absolute paths for Nextflow
    let temp_template = template_nf
        .canonicalize()
        .unwrap_or_else(|_| template_nf.clone());
    let temp_config = nextflow_config
        .canonicalize()
        .unwrap_or_else(|_| nextflow_config.clone());

    // Get workflow file
    let workflow_file = project_path
        .join("workflow.nf")
        .canonicalize()
        .context("Failed to resolve workflow.nf path")?;

    // Create assets directory if it doesn't exist (we already determined it above)
    if !assets_dir.exists() {
        fs::create_dir_all(&assets_dir)?;
    }

    let assets_dir = assets_dir
        .canonicalize()
        .context("Failed to resolve assets directory path")?;

    // Create results directory
    let results_base = if let Some(ref custom_dir) = params.results_dir {
        custom_dir.as_str()
    } else if params.test {
        "results-test"
    } else {
        "results-real"
    };
    let results_dir = project_path.join(results_base);
    if !results_dir.exists() {
        fs::create_dir_all(&results_dir)?;
    }

    let results_dir = results_dir
        .canonicalize()
        .context("Failed to resolve results directory path")?;

    info!(
        "Running sheet workflow '{}' from project '{}'",
        config.workflow, config.name
    );

    let nextflow_log_path = project_path.join(".nextflow.log");
    fs::remove_file(&nextflow_log_path).ok();

    // Get configured Nextflow path or use default
    let bv_config = crate::config::get_config().ok();
    let nextflow_cmd = resolve_binary_path(bv_config.as_ref(), "nextflow")
        .unwrap_or_else(|| "nextflow".to_string());

    // Build Nextflow command
    let mut cmd = Command::new(&nextflow_cmd);

    // Set working directory to project directory
    cmd.current_dir(&project_path);

    cmd.arg("-log").arg(&nextflow_log_path);

    cmd.arg("run")
        .arg(&temp_template)
        .arg("--samplesheet")
        .arg(samplesheet_path.canonicalize().unwrap_or(samplesheet_path));

    if schema_path.exists() {
        cmd.arg("--schema_yaml")
            .arg(schema_path.canonicalize().unwrap_or(schema_path));
    }

    cmd.arg("--work_flow_file")
        .arg(workflow_file.to_string_lossy().as_ref())
        .arg("--assets_dir")
        .arg(assets_dir.to_string_lossy().as_ref())
        .arg("--results_dir")
        .arg(results_dir.to_string_lossy().as_ref());

    if params.resume {
        cmd.arg("-resume");
    }

    if let Some(work_dir) = &params.work_dir {
        cmd.arg("-work-dir");
        cmd.arg(work_dir);
    }

    // Docker/Singularity configuration
    // Only add -with-docker if user hasn't provided it in nextflow_args
    let has_docker_arg = params.nextflow_args.iter().any(|arg| {
        arg.starts_with("-with-docker")
            || arg.starts_with("-with-singularity")
            || arg.starts_with("-with-podman")
    });

    if params.with_docker && !has_docker_arg {
        cmd.arg("-with-docker");
    }

    // Add additional Nextflow arguments
    for arg in &params.nextflow_args {
        cmd.arg(arg);
    }

    // Add config file
    cmd.arg("-c").arg(&temp_config);
    cmd.arg("-log").arg(&nextflow_log_path);

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
    append_desktop_log("INFO", &format!("Nextflow command: {}", cmd_str));

    if params.dry_run {
        println!("{}", "[DRY RUN] Would execute the above command".yellow());
        append_desktop_log("INFO", "[DRY RUN] Would execute the above command");
        return Ok(());
    }

    // Execute Nextflow
    println!("Executing Nextflow sheet workflow...");
    append_desktop_log("INFO", "Executing Nextflow sheet workflow...");
    let status =
        execute_with_logging(cmd, Some(nextflow_log_path)).context("Failed to execute Nextflow")?;

    if !status.success() {
        append_desktop_log("ERROR", "Nextflow execution failed");
        return Err(anyhow!("Nextflow execution failed"));
    }

    println!(
        "{}",
        "Sheet workflow completed successfully!".green().bold()
    );
    append_desktop_log("INFO", "Sheet workflow completed successfully!");
    Ok(())
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

    // Check if this is a sheet template project
    let is_sheet_template = config
        .template
        .as_ref()
        .map(|t| t == "sheet")
        .unwrap_or(false);

    if is_sheet_template {
        return execute_sheet_workflow(&params, &config).await;
    }

    // Parse participant source for non-sheet workflows
    let source = ParticipantSource::parse(&params.participant_source)?;

    // Fetch participant file
    let (yaml_content, fragment) = fetch_participant_file(&source, params.download).await?;

    // Extract participant data
    let (mut participant, mock_key) =
        extract_participant_data(&yaml_content, fragment, params.test)?;

    // Ensure all required files exist (download if needed)
    participant =
        ensure_files_exist(&participant, params.download, &source, mock_key.as_deref()).await?;

    // Determine which template to use
    // Priority: CLI flag > project.yaml > default
    let template_name = params
        .template
        .or(config.template.clone())
        .unwrap_or_else(|| "default".to_string());

    // Get BioVault environment directory
    let biovault_home = crate::config::get_biovault_home()?;
    let env_dir = biovault_home.join("env").join(&template_name);

    // Check if templates exist
    let template_nf = env_dir.join("template.nf");
    let nextflow_config = env_dir.join("nextflow.config");

    if !template_nf.exists() || !nextflow_config.exists() {
        return Err(Error::TemplatesNotFound.into());
    }

    println!("Using template: {}", template_name);

    // Use templates directly from env dir instead of copying to temp
    // Convert to absolute paths for Nextflow
    let temp_template = template_nf
        .canonicalize()
        .unwrap_or_else(|_| template_nf.clone());
    let temp_config = nextflow_config
        .canonicalize()
        .unwrap_or_else(|_| nextflow_config.clone());

    // Determine assets directory
    // Prefer the conventional 'assets' folder. If the first assets entry is an existing directory,
    // use it; otherwise do NOT create a directory from a file name like 'eye_color.py'.
    let mut assets_dir = project_path.join("assets");
    if !config.assets.is_empty() {
        let candidate = project_path.join(&config.assets[0]);
        if candidate.is_dir() {
            assets_dir = candidate;
        }
    }

    // Create assets directory if it doesn't exist (only for directories)
    if !assets_dir.exists() {
        fs::create_dir_all(&assets_dir)?;
    }

    // Get absolute path for assets directory
    let assets_dir = assets_dir
        .canonicalize()
        .context("Failed to resolve assets directory path")?;

    // Create results directory for this participant
    // Use results-test for sample data or if test flag is set
    let is_sample_data = matches!(source, ParticipantSource::SampleDataId(_));
    let results_base = if let Some(ref custom_dir) = params.results_dir {
        custom_dir.as_str()
    } else if params.test || is_sample_data {
        "results-test"
    } else {
        "results-real"
    };
    let results_dir = project_path.join(results_base).join(&participant.id);
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

    if is_sample_data {
        println!(
            "Processing sample data participant: {}",
            participant.id.cyan()
        );
    } else {
        println!("Processing participant: {}", participant.id.cyan());
    }

    let nextflow_log_path = project_path.join(".nextflow.log");
    fs::remove_file(&nextflow_log_path).ok();

    // Get configured Nextflow path or use default
    let bv_config = crate::config::get_config().ok();
    let nextflow_cmd = resolve_binary_path(bv_config.as_ref(), "nextflow")
        .unwrap_or_else(|| "nextflow".to_string());

    // Build Nextflow command
    let mut cmd = Command::new(&nextflow_cmd);

    // Set working directory to project directory
    cmd.current_dir(&project_path);

    cmd.arg("-log").arg(&nextflow_log_path);

    cmd.arg("run")
        .arg(&temp_template)
        .arg("--participant_id")
        .arg(&participant.id);

    // Add CRAM-specific parameters if present
    if let Some(ref_version) = &participant.ref_version {
        cmd.arg("--ref_version").arg(ref_version);
    }
    if let Some(ref_path) = &participant.ref_path {
        cmd.arg("--ref").arg(ref_path);
    }
    if let Some(ref_index) = &participant.ref_index {
        cmd.arg("--ref_index").arg(ref_index);
    }
    if let Some(aligned) = &participant.aligned {
        cmd.arg("--aligned").arg(aligned);
    }
    if let Some(aligned_index) = &participant.aligned_index {
        cmd.arg("--aligned_index").arg(aligned_index);
    }

    // Add SNP-specific parameters if present
    if let Some(snp) = &participant.snp {
        // Convert to absolute path if it's a file path
        let snp_path = PathBuf::from(snp);
        let snp_abs = if snp_path.exists() {
            snp_path.canonicalize().unwrap_or(snp_path)
        } else {
            snp_path
        };
        cmd.arg("--snp").arg(snp_abs);
    }

    cmd.arg("--work_flow_file")
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
    // Only add -with-docker if user hasn't provided it in nextflow_args
    let has_docker_arg = params.nextflow_args.iter().any(|arg| {
        arg.starts_with("-with-docker")
            || arg.starts_with("-with-singularity")
            || arg.starts_with("-with-podman")
    });

    if params.with_docker && !has_docker_arg {
        cmd.arg("-with-docker");
    }

    // Add additional Nextflow arguments
    for arg in &params.nextflow_args {
        cmd.arg(arg);
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
    append_desktop_log("INFO", &format!("Nextflow command: {}", cmd_str));

    if params.dry_run {
        println!("{}", "[DRY RUN] Would execute the above command".yellow());
        append_desktop_log("INFO", "[DRY RUN] Would execute the above command");
        return Ok(());
    }

    // Execute Nextflow
    println!("Executing Nextflow workflow...");
    append_desktop_log("INFO", "Executing Nextflow workflow...");
    let status =
        execute_with_logging(cmd, Some(nextflow_log_path)).context("Failed to execute Nextflow")?;

    if !status.success() {
        append_desktop_log("ERROR", "Nextflow execution failed");
        return Err(anyhow!("Nextflow execution failed"));
    }

    println!("{}", "Workflow completed successfully!".green().bold());
    append_desktop_log("INFO", "Workflow completed successfully!");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::cli::download_cache::manifest::Manifest;
    use crate::config::{clear_test_biovault_home, set_test_biovault_home};
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_deserialize_string_as_vec() {
        #[derive(Deserialize)]
        struct TestStruct {
            #[serde(deserialize_with = "deserialize_string_or_vec")]
            assets: Vec<String>,
        }

        let yaml_str = "assets: single_asset";
        let result: TestStruct = serde_yaml::from_str(yaml_str).unwrap();
        assert_eq!(result.assets, vec!["single_asset"]);

        let yaml_list = "assets:\n  - asset1\n  - asset2";
        let result: TestStruct = serde_yaml::from_str(yaml_list).unwrap();
        assert_eq!(result.assets, vec!["asset1", "asset2"]);
    }

    #[test]
    fn test_project_config_deserialize() {
        let yaml = r#"
name: test_project
author: test@example.com
workflow: test.nf
template: test_template
assets: test_asset
participants:
  - p1
  - p2
"#;
        let config: ProjectConfig = serde_yaml::from_str(yaml).unwrap();
        assert_eq!(config.name, "test_project");
        assert_eq!(config.author, "test@example.com");
        assert_eq!(config.workflow, "test.nf");
        assert_eq!(config.template, Some("test_template".to_string()));
        assert_eq!(config.assets, vec!["test_asset"]);
        assert_eq!(config.participants, vec!["p1", "p2"]);
    }

    #[test]
    fn test_participant_data_serialize() {
        let participant = ParticipantData {
            id: "test_id".to_string(),
            ref_version: Some("GRCh38".to_string()),
            ref_path: Some("/path/to/ref.fa".to_string()),
            ref_index: Some("/path/to/ref.fa.fai".to_string()),
            aligned: Some("/path/to/aligned.cram".to_string()),
            aligned_index: Some("/path/to/aligned.cram.crai".to_string()),
            ref_b3sum: None,
            ref_index_b3sum: None,
            aligned_b3sum: None,
            aligned_index_b3sum: None,
            snp: None,
            snp_b3sum: None,
            uncompress: None,
        };

        let yaml = serde_yaml::to_string(&participant).unwrap();
        assert!(yaml.contains("id: test_id"));
        assert!(yaml.contains("ref_version: GRCh38"));
        assert!(yaml.contains("ref: /path/to/ref.fa"));
    }

    #[test]
    fn test_participant_data_snp_variant() {
        let participant = ParticipantData {
            id: "snp_test".to_string(),
            ref_version: None,
            ref_path: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            ref_b3sum: None,
            ref_index_b3sum: None,
            aligned_b3sum: None,
            aligned_index_b3sum: None,
            snp: Some("/path/to/snp.vcf".to_string()),
            snp_b3sum: Some("abc123".to_string()),
            uncompress: Some(true),
        };

        assert_eq!(participant.id, "snp_test");
        assert_eq!(participant.snp, Some("/path/to/snp.vcf".to_string()));
        assert_eq!(participant.snp_b3sum, Some("abc123".to_string()));
        assert_eq!(participant.uncompress, Some(true));
    }

    #[test]
    fn test_run_params_default() {
        let params = RunParams {
            project_folder: "/test".to_string(),
            participant_source: "participants.yaml#TEST".to_string(),
            test: false,
            download: false,
            dry_run: false,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: None,
            results_dir: None,
            nextflow_args: vec![],
        };

        assert_eq!(params.project_folder, "/test");
        assert_eq!(params.participant_source, "participants.yaml#TEST");
        assert!(!params.test);
        assert!(!params.download);
        assert!(!params.with_docker);
        assert!(!params.dry_run);
        assert!(params.work_dir.is_none());
        assert!(!params.resume);
        assert!(params.template.is_none());
        assert!(params.nextflow_args.is_empty());
    }

    #[test]
    fn test_participant_data_clone() {
        let original = ParticipantData {
            id: "clone_test".to_string(),
            ref_version: Some("GRCh37".to_string()),
            ref_path: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            ref_b3sum: None,
            ref_index_b3sum: None,
            aligned_b3sum: None,
            aligned_index_b3sum: None,
            snp: None,
            snp_b3sum: None,
            uncompress: None,
        };

        let cloned = original.clone();
        assert_eq!(cloned.id, original.id);
        assert_eq!(cloned.ref_version, original.ref_version);
    }

    #[test]
    fn test_project_config_minimal() {
        let yaml = r#"
name: minimal
author: user@example.com
workflow: main.nf
"#;
        let config: ProjectConfig = serde_yaml::from_str(yaml).unwrap();
        assert_eq!(config.name, "minimal");
        assert_eq!(config.author, "user@example.com");
        assert_eq!(config.workflow, "main.nf");
        assert!(config.template.is_none());
        assert!(config.assets.is_empty());
        assert!(config.participants.is_empty());
    }

    #[test]
    fn participant_source_parse_variants() {
        // Local file with fragment
        let ps = ParticipantSource::parse("/tmp/p.yml#participants.A").unwrap();
        match ps {
            ParticipantSource::LocalFile(path, frag) => {
                assert!(path.ends_with("p.yml"));
                assert_eq!(frag.as_deref(), Some("participants.A"));
            }
            _ => panic!("expected LocalFile"),
        }

        // HTTP url
        let ps = ParticipantSource::parse("https://example.com/p.yml#participants.A").unwrap();
        match ps {
            ParticipantSource::HttpUrl(u) => assert!(u.starts_with("https://example.com")),
            _ => panic!("expected HttpUrl"),
        }

        // Syft URL
        let ps = ParticipantSource::parse("syft://user@example.com/path#participants.A").unwrap();
        match ps {
            ParticipantSource::SyftUrl(u) => {
                assert_eq!(u.email, "user@example.com");
            }
            _ => panic!("expected SyftUrl"),
        }

        // Sample data id present in embedded config (uses include_str)
        let ps = ParticipantSource::parse("NA06985").unwrap();
        match ps {
            ParticipantSource::SampleDataId(id) => assert_eq!(id, "NA06985"),
            _ => panic!("expected SampleDataId"),
        }
    }

    #[test]
    fn extract_participant_data_happy_and_error_paths() {
        // Happy path
        let yaml = r#"
participants:
  TEST:
    ref_version: GRCh38
"#;
        let (p, mock) = extract_participant_data(yaml, Some("participants.TEST".into()), false)
            .expect("parse ok");
        assert_eq!(p.id, "TEST");
        assert_eq!(p.ref_version.as_deref(), Some("GRCh38"));
        assert!(mock.is_none());

        // Missing fragment
        let err = extract_participant_data(yaml, None, false).unwrap_err();
        assert!(format!("{}", err).contains("No participant specified"));

        // Wrong fragment prefix
        let err = extract_participant_data(yaml, Some("foo.TEST".into()), false).unwrap_err();
        assert!(format!("{}", err).contains("Invalid fragment"));

        // Participant not found
        let err = extract_participant_data(yaml, Some("participants.X".into()), false).unwrap_err();
        assert!(format!("{}", err).contains("not found"));
    }

    #[tokio::test]
    async fn ensure_files_exist_with_local_paths() {
        // Prepare local files
        let td = TempDir::new().unwrap();
        let refp = td.path().join("ref.fa");
        let refi = td.path().join("ref.fa.fai");
        let cram = td.path().join("aln.cram");
        let crai = td.path().join("aln.cram.crai");
        fs::write(&refp, b"ref").unwrap();
        fs::write(&refi, b"idx").unwrap();
        fs::write(&cram, b"cram").unwrap();
        fs::write(&crai, b"crai").unwrap();

        // Participant with local paths
        let p = ParticipantData {
            id: "P1".into(),
            ref_version: Some("GRCh38".into()),
            ref_path: Some(refp.to_string_lossy().to_string()),
            ref_index: Some(refi.to_string_lossy().to_string()),
            aligned: Some(cram.to_string_lossy().to_string()),
            aligned_index: Some(crai.to_string_lossy().to_string()),
            ref_b3sum: None,
            ref_index_b3sum: None,
            aligned_b3sum: None,
            aligned_index_b3sum: None,
            snp: None,
            snp_b3sum: None,
            uncompress: None,
        };

        // Ensure BIOVAULT home is isolated
        let home = TempDir::new().unwrap();
        set_test_biovault_home(home.path());

        let src = ParticipantSource::LocalFile(PathBuf::from("participants.yaml"), None);
        // Point cache dir to a writable temp to avoid platform HOME surprises
        let cache_td = TempDir::new().unwrap();
        std::env::set_var(
            "BIOVAULT_CACHE_DIR",
            cache_td.path().to_string_lossy().to_string(),
        );
        let out = ensure_files_exist(&p, false, &src, None).await.unwrap();
        std::env::remove_var("BIOVAULT_CACHE_DIR");
        assert_eq!(out.ref_path.as_deref(), p.ref_path.as_deref());
        assert_eq!(out.ref_index.as_deref(), p.ref_index.as_deref());
        assert_eq!(out.aligned.as_deref(), p.aligned.as_deref());
        assert_eq!(out.aligned_index.as_deref(), p.aligned_index.as_deref());

        clear_test_biovault_home();
    }

    #[tokio::test]
    async fn execute_dry_run_minimal_project() {
        // Isolate BIOVAULT home and create template files
        let bv_home = TempDir::new().unwrap();
        let env_dir = bv_home.path().join("env").join("test_tpl");
        fs::create_dir_all(&env_dir).unwrap();
        fs::write(env_dir.join("template.nf"), "// template").unwrap();
        fs::write(env_dir.join("nextflow.config"), "// config").unwrap();
        set_test_biovault_home(bv_home.path());

        // Create minimal project
        let proj = TempDir::new().unwrap();
        fs::write(
            proj.path().join("project.yaml"),
            "name: p\nauthor: a\nworkflow: main.nf\ntemplate: test_tpl\n",
        )
        .unwrap();
        fs::write(proj.path().join("workflow.nf"), "// wf").unwrap();
        fs::write(
            proj.path().join("participants.yaml"),
            "participants:\n  X:\n    ref_version: GRCh38\n",
        )
        .unwrap();

        let params = RunParams {
            project_folder: proj.path().to_string_lossy().to_string(),
            participant_source: proj
                .path()
                .join("participants.yaml#participants.X")
                .to_string_lossy()
                .to_string(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: Some("test_tpl".into()),
            results_dir: None,
            nextflow_args: vec![],
        };

        // Use a writable cache dir during test
        // Create the cache directory structure to match what DownloadCache expects
        let cache_td = TempDir::new().unwrap();
        let cache_dir = cache_td.path().join("data").join("cache");
        fs::create_dir_all(&cache_dir).unwrap();
        std::env::set_var(
            "BIOVAULT_CACHE_DIR",
            cache_dir.to_string_lossy().to_string(),
        );
        // Should return Ok before trying to execute nextflow
        execute(params).await.expect("dry-run ok");
        std::env::remove_var("BIOVAULT_CACHE_DIR");

        clear_test_biovault_home();
    }

    #[tokio::test]
    async fn execute_dry_run_with_all_params_and_fields() {
        // Isolate BIOVAULT home and create template files
        let bv_home = TempDir::new().unwrap();
        let env_dir = bv_home.path().join("env").join("full_tpl");
        fs::create_dir_all(&env_dir).unwrap();
        fs::write(env_dir.join("template.nf"), "// template").unwrap();
        fs::write(env_dir.join("nextflow.config"), "// config").unwrap();
        set_test_biovault_home(bv_home.path());

        // Create minimal project with assets dir
        let proj = TempDir::new().unwrap();
        fs::create_dir_all(proj.path().join("assets")).unwrap();
        fs::write(
            proj.path().join("project.yaml"),
            "name: p\nauthor: a\nworkflow: main.nf\ntemplate: full_tpl\nassets: assets\n",
        )
        .unwrap();
        fs::write(proj.path().join("workflow.nf"), "// wf").unwrap();

        // Participant with many optional fields set and corresponding local files
        fs::write(proj.path().join("ref.fa"), b"ref").unwrap();
        fs::write(proj.path().join("ref.fa.fai"), b"idx").unwrap();
        fs::write(proj.path().join("aln.cram"), b"cram").unwrap();
        fs::write(proj.path().join("aln.cram.crai"), b"crai").unwrap();
        fs::write(proj.path().join("snp.vcf"), b"##vcf\n").unwrap();
        // Participant with many optional fields set
        let participants_yaml = format!(
            "participants:\n  Y:\n    ref_version: GRCh38\n    ref: {}\n    ref_index: {}\n    aligned: {}\n    aligned_index: {}\n    snp: {}\n",
            proj.path().join("ref.fa").display(),
            proj.path().join("ref.fa.fai").display(),
            proj.path().join("aln.cram").display(),
            proj.path().join("aln.cram.crai").display(),
            proj.path().join("snp.vcf").display(),
        );
        fs::write(proj.path().join("participants.yaml"), participants_yaml).unwrap();

        let params = RunParams {
            project_folder: proj.path().to_string_lossy().to_string(),
            participant_source: proj
                .path()
                .join("participants.yaml#participants.Y")
                .to_string_lossy()
                .to_string(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: true,
            work_dir: Some("workdir".into()),
            resume: true,
            template: Some("full_tpl".into()),
            results_dir: None,
            nextflow_args: vec![],
        };

        let cache_td = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_CACHE_DIR", cache_td.path());
        execute(params).await.expect("dry-run ok");
        std::env::remove_var("BIOVAULT_CACHE_DIR");

        clear_test_biovault_home();
    }

    // Removed HTTP cache test to avoid network in restricted environments

    #[tokio::test]
    async fn execute_errors_when_paths_missing() {
        // Missing project directory
        let params = RunParams {
            project_folder: "/definitely/not/here".into(),
            participant_source: "participants.yaml#participants.X".into(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: None,
            results_dir: None,
            nextflow_args: vec![],
        };
        assert!(execute(params).await.is_err());

        // Project exists but missing project.yaml
        let proj = TempDir::new().unwrap();
        let params = RunParams {
            project_folder: proj.path().to_string_lossy().to_string(),
            participant_source: "participants.yaml#participants.X".into(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: None,
            results_dir: None,
            nextflow_args: vec![],
        };
        assert!(execute(params).await.is_err());

        // project.yaml present, workflow.nf missing
        fs::write(
            proj.path().join("project.yaml"),
            "name: p\nauthor: a\nworkflow: main.nf\n",
        )
        .unwrap();
        let params = RunParams {
            project_folder: proj.path().to_string_lossy().to_string(),
            participant_source: "participants.yaml#participants.X".into(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: None,
            results_dir: None,
            nextflow_args: vec![],
        };
        assert!(execute(params).await.is_err());

        // workflow present but template missing in env dir -> TemplatesNotFound
        fs::write(proj.path().join("workflow.nf"), "// wf").unwrap();
        // Participants file with minimal entry
        fs::write(
            proj.path().join("participants.yaml"),
            "participants:\n  X:\n    ref_version: GRCh38\n",
        )
        .unwrap();
        // point test home to an empty env dir
        let bv_home = TempDir::new().unwrap();
        set_test_biovault_home(bv_home.path());
        let params = RunParams {
            project_folder: proj.path().to_string_lossy().to_string(),
            participant_source: proj
                .path()
                .join("participants.yaml#participants.X")
                .to_string_lossy()
                .to_string(),
            test: false,
            download: false,
            dry_run: true,
            with_docker: false,
            work_dir: None,
            resume: false,
            template: Some("missing_tpl".into()),
            results_dir: None,
            nextflow_args: vec![],
        };
        assert!(execute(params).await.is_err());

        clear_test_biovault_home();
    }

    #[tokio::test]
    async fn fetch_participant_file_local_missing_errors() {
        let res = fetch_participant_file(
            &ParticipantSource::LocalFile(
                PathBuf::from("/nope/participants.yaml"),
                Some("participants.X".into()),
            ),
            false,
        )
        .await;
        assert!(res.is_err());
    }

    #[test]
    fn extract_participant_data_mock_branch() {
        let yaml = r#"
participants:
  P:
    mock:
      ref_version: GRCh38
      aligned: /tmp/test.cram
      aligned_index: /tmp/test.cram.crai
"#;
        let (p, mock) =
            extract_participant_data(yaml, Some("participants.P".into()), true).unwrap();
        assert_eq!(p.id, "P");
        assert_eq!(p.ref_version.as_deref(), Some("GRCh38"));
        assert_eq!(mock.as_deref(), Some("mock_data_grch38"));
    }
}
