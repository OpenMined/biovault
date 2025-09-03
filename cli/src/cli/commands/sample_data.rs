use crate::cli::download_cache::{
    ChecksumPolicy, ChecksumPolicyType, DownloadCache, DownloadOptions,
};
use crate::error::Error;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use uuid::Uuid;

const SAMPLE_DATA_YAML: &str = include_str!("../../sample_data.yaml");

#[derive(Debug, Serialize, Deserialize)]
struct SampleDataConfig {
    sample_data_urls: HashMap<String, ParticipantData>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ParticipantData {
    ref_version: String,
    #[serde(rename = "ref")]
    ref_url: String,
    ref_index: String,
    aligned: String,
    aligned_index: String,
    #[serde(default)]
    ref_b3sum: String,
    #[serde(default)]
    ref_index_b3sum: String,
    #[serde(default)]
    aligned_b3sum: String,
    #[serde(default)]
    aligned_index_b3sum: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct ParticipantsFile {
    participant: HashMap<String, ParticipantRecord>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ParticipantRecord {
    ref_version: String,
    #[serde(rename = "ref")]
    ref_path: String,
    ref_index: String,
    aligned: String,
    aligned_index: String,
}

pub async fn fetch(participant_ids: Option<Vec<String>>, all: bool) -> anyhow::Result<()> {
    let config: SampleDataConfig = serde_yaml::from_str(SAMPLE_DATA_YAML)
        .context("Failed to parse embedded sample data configuration")?;

    let participants_to_fetch = determine_participants_to_fetch(&config, participant_ids, all)?;

    if participants_to_fetch.is_empty() {
        println!("No participants specified to fetch");
        return Ok(());
    }

    let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
        PathBuf::from(test_home)
    } else {
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
    };

    let sample_data_dir = home_dir.join(".biovault").join("data").join("sample");
    let reference_dir = sample_data_dir.join("reference");

    fs::create_dir_all(&reference_dir).context("Failed to create reference directory")?;

    println!(
        "Fetching sample data for {} participant(s):",
        participants_to_fetch.len()
    );
    for id in &participants_to_fetch {
        println!("  - {}", id);
    }
    println!();

    let participants_file_path = sample_data_dir.join("participants.yaml");
    let mut participants_file = load_or_create_participants_file(&participants_file_path)?;

    // Initialize download cache
    let mut download_cache = DownloadCache::new(None)?;

    for participant_id in participants_to_fetch {
        println!("\n{}", "=".repeat(60));
        println!("Fetching data for participant: {}", participant_id);
        println!("{}", "=".repeat(60));

        let participant_data = config
            .sample_data_urls
            .get(&participant_id)
            .ok_or_else(|| {
                anyhow::anyhow!("Participant {} not found in configuration", participant_id)
            })?;

        let participant_dir = sample_data_dir.join(&participant_id);
        fs::create_dir_all(&participant_dir).context("Failed to create participant directory")?;

        // Extract filenames from URLs
        let ref_filename = extract_filename_from_url(&participant_data.ref_url)?;
        let ref_index_filename = extract_filename_from_url(&participant_data.ref_index)?;
        let aligned_filename = extract_filename_from_url(&participant_data.aligned)?;
        let aligned_index_filename = extract_filename_from_url(&participant_data.aligned_index)?;

        let downloads = vec![
            (
                &participant_data.ref_url,
                reference_dir.join(&ref_filename),
                "Reference genome",
                &participant_data.ref_b3sum,
            ),
            (
                &participant_data.ref_index,
                reference_dir.join(&ref_index_filename),
                "Reference index",
                &participant_data.ref_index_b3sum,
            ),
            (
                &participant_data.aligned,
                participant_dir.join(&aligned_filename),
                "Aligned CRAM",
                &participant_data.aligned_b3sum,
            ),
            (
                &participant_data.aligned_index,
                participant_dir.join(&aligned_index_filename),
                "CRAM index",
                &participant_data.aligned_index_b3sum,
            ),
        ];

        for (url, target_path, description, expected_b3sum) in downloads {
            println!(
                "\n  Processing {}: {}",
                description,
                target_path.file_name().unwrap().to_string_lossy()
            );

            // Set up download options based on whether we have a checksum
            let options = if !expected_b3sum.is_empty() {
                DownloadOptions {
                    checksum_policy: ChecksumPolicy {
                        policy_type: ChecksumPolicyType::Required,
                        expected_hash: Some(expected_b3sum.to_string()),
                    },
                    ..Default::default()
                }
            } else {
                DownloadOptions::default()
            };

            // Download to a temporary location (will be cached)
            let temp_filename = target_path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("download");
            let temp_path =
                std::env::temp_dir().join(format!("bv_{}_{}", temp_filename, Uuid::new_v4()));

            // Use the download cache
            download_cache
                .download_with_cache(url, &temp_path, options)
                .await
                .with_context(|| format!("Failed to download {}", description))?;

            // The file is now in cache, create a symlink to it
            if !expected_b3sum.is_empty() {
                let cache_base = home_dir.join(".biovault").join("data").join("cache");
                let cache_path = cache_base.join("by-hash").join(expected_b3sum);

                // Remove any existing file or symlink at target
                if target_path.exists() || target_path.is_symlink() {
                    fs::remove_file(&target_path).ok();
                }

                // Create symlink to cache
                #[cfg(unix)]
                {
                    std::os::unix::fs::symlink(&cache_path, &target_path)
                        .with_context(|| format!("Failed to create symlink for {}", description))?;
                }
                #[cfg(windows)]
                {
                    std::os::windows::fs::symlink_file(&cache_path, &target_path)
                        .with_context(|| format!("Failed to create symlink for {}", description))?;
                }

                println!("    ✓ Linked to cache (saving disk space)");
            } else {
                // If no checksum, we need to copy the file from temp to target
                // This shouldn't happen in practice since all sample data has checksums
                fs::rename(&temp_path, &target_path)
                    .or_else(|_| fs::copy(&temp_path, &target_path).map(|_| ()))
                    .with_context(|| format!("Failed to move {} to target", description))?;
            }

            // Clean up temp file if it exists
            fs::remove_file(&temp_path).ok();
        }

        let participant_record = ParticipantRecord {
            ref_version: participant_data.ref_version.clone(),
            ref_path: format!("./reference/{}", ref_filename),
            ref_index: format!("./reference/{}", ref_index_filename),
            aligned: format!("./{}/{}", participant_id, aligned_filename),
            aligned_index: format!("./{}/{}", participant_id, aligned_index_filename),
        };

        participants_file
            .participant
            .insert(participant_id.clone(), participant_record);

        save_participants_file(&participants_file_path, &participants_file)?;
        println!("  ✓ Updated participants.yaml");
    }

    println!("\n{}", "=".repeat(60));
    println!("✓ Sample data fetch complete!");
    println!("  Data location: {}", sample_data_dir.display());
    println!("  Participants file: {}", participants_file_path.display());
    println!("{}", "=".repeat(60));

    Ok(())
}

fn extract_filename_from_url(url: &str) -> anyhow::Result<String> {
    Ok(url
        .rsplit('/')
        .next()
        .ok_or_else(|| anyhow::anyhow!("Invalid URL: {}", url))?
        .split('#')
        .next()
        .unwrap()
        .split('?')
        .next()
        .unwrap()
        .to_string())
}

fn determine_participants_to_fetch(
    config: &SampleDataConfig,
    participant_ids: Option<Vec<String>>,
    all: bool,
) -> anyhow::Result<Vec<String>> {
    if all {
        Ok(config.sample_data_urls.keys().cloned().collect())
    } else if let Some(ids) = participant_ids {
        for id in &ids {
            if !config.sample_data_urls.contains_key(id) {
                return Err(Error::ParticipantNotFound(id.clone()).into());
            }
        }
        Ok(ids)
    } else {
        Err(Error::NoParticipantsSpecified.into())
    }
}

fn load_or_create_participants_file(path: &Path) -> anyhow::Result<ParticipantsFile> {
    if path.exists() {
        let content =
            fs::read_to_string(path).context("Failed to read existing participants.yaml")?;
        Ok(serde_yaml::from_str(&content).context("Failed to parse existing participants.yaml")?)
    } else {
        Ok(ParticipantsFile {
            participant: HashMap::new(),
        })
    }
}

fn save_participants_file(path: &Path, participants: &ParticipantsFile) -> anyhow::Result<()> {
    let yaml =
        serde_yaml::to_string(participants).context("Failed to serialize participants data")?;
    fs::write(path, yaml).context("Failed to write participants.yaml")?;
    Ok(())
}

pub async fn list() -> anyhow::Result<()> {
    let config: SampleDataConfig = serde_yaml::from_str(SAMPLE_DATA_YAML)
        .context("Failed to parse embedded sample data configuration")?;

    println!("Available sample data:");
    println!("{}", "=".repeat(60));

    for (participant_id, data) in &config.sample_data_urls {
        println!("\nParticipant ID: {}", participant_id);
        println!("  ref_version: {}", data.ref_version);
        println!("  ref: {}", data.ref_url);
        println!("  ref_index: {}", data.ref_index);
        println!("  aligned: {}", data.aligned);
        println!("  aligned_index: {}", data.aligned_index);
    }

    println!("\n{}", "=".repeat(60));
    println!("Use 'bv sample-data fetch <PARTICIPANT_ID>' to download");
    println!("Use 'bv sample-data fetch --all' to download all samples");

    Ok(())
}
