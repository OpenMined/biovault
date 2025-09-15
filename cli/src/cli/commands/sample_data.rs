use crate::cli::download_cache::{
    ChecksumPolicy, ChecksumPolicyType, DownloadCache, DownloadOptions,
};
use crate::error::Error;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use uuid::Uuid;

const SAMPLE_DATA_YAML: &str = include_str!("../../sample_data.yaml");

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(untagged)]
enum AlignedUrl {
    Single(String),
    Multiple(Vec<String>),
}

impl AlignedUrl {
    fn to_vec(&self) -> Vec<String> {
        match self {
            AlignedUrl::Single(url) => vec![url.clone()],
            AlignedUrl::Multiple(urls) => urls.clone(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(untagged)]
enum AlignedChecksum {
    Single(String),
    Multiple(Vec<String>),
}

impl AlignedChecksum {
    fn to_vec(&self) -> Vec<String> {
        match self {
            AlignedChecksum::Single(checksum) => vec![checksum.clone()],
            AlignedChecksum::Multiple(checksums) => checksums.clone(),
        }
    }
}

impl Default for AlignedChecksum {
    fn default() -> Self {
        AlignedChecksum::Single(String::new())
    }
}

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
    aligned: AlignedUrl,
    aligned_index: String,
    #[serde(default)]
    ref_b3sum: String,
    #[serde(default)]
    ref_index_b3sum: String,
    #[serde(default)]
    aligned_b3sum: AlignedChecksum,
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

    let biovault_home = crate::config::get_biovault_home()?;
    let sample_data_dir = biovault_home.join("data").join("sample");
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

        // Handle aligned URLs (can be single or multiple)
        let aligned_urls = participant_data.aligned.to_vec();
        let aligned_checksums = participant_data.aligned_b3sum.to_vec();
        let aligned_index_filename = if !participant_data.aligned_index.is_empty() {
            extract_filename_from_url(&participant_data.aligned_index)?
        } else {
            String::new()
        };

        // Build downloads list - first add reference files
        let mut downloads = vec![
            (
                participant_data.ref_url.clone(),
                reference_dir.join(&ref_filename),
                "Reference genome".to_string(),
                participant_data.ref_b3sum.clone(),
            ),
            (
                participant_data.ref_index.clone(),
                reference_dir.join(&ref_index_filename),
                "Reference index".to_string(),
                participant_data.ref_index_b3sum.clone(),
            ),
        ];

        // Handle aligned files - might be split into multiple parts
        if !aligned_urls.is_empty() && !aligned_urls[0].is_empty() {
            if aligned_urls.len() > 1 {
                // Multiple parts - need to download all and combine
                let mut part_files = Vec::new();
                for (i, url) in aligned_urls.iter().enumerate() {
                    let part_filename = extract_filename_from_url(url)?;
                    let checksum = if i < aligned_checksums.len() {
                        aligned_checksums[i].clone()
                    } else {
                        String::new()
                    };
                    part_files.push(part_filename.clone());
                    downloads.push((
                        url.clone(),
                        participant_dir.join(&part_filename),
                        format!("Aligned CRAM part {}/{}", i + 1, aligned_urls.len()),
                        checksum,
                    ));
                }
            } else {
                // Single aligned file
                let aligned_filename = extract_filename_from_url(&aligned_urls[0])?;
                downloads.push((
                    aligned_urls[0].clone(),
                    participant_dir.join(&aligned_filename),
                    "Aligned CRAM".to_string(),
                    if !aligned_checksums.is_empty() {
                        aligned_checksums[0].clone()
                    } else {
                        String::new()
                    },
                ));
            }
        }

        // Add aligned index if present
        if !participant_data.aligned_index.is_empty() {
            downloads.push((
                participant_data.aligned_index.clone(),
                participant_dir.join(&aligned_index_filename),
                "CRAM index".to_string(),
                participant_data.aligned_index_b3sum.clone(),
            ));
        };

        for (url, target_path, description, expected_b3sum) in &downloads {
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
                let cache_base = crate::config::get_cache_dir()?;
                let cache_path = cache_base.join("by-hash").join(expected_b3sum);

                // Remove any existing file or symlink at target
                if target_path.exists() || target_path.is_symlink() {
                    fs::remove_file(target_path).ok();
                }

                // Create symlink to cache
                #[cfg(unix)]
                {
                    std::os::unix::fs::symlink(&cache_path, target_path)
                        .with_context(|| format!("Failed to create symlink for {}", description))?;
                }
                #[cfg(windows)]
                {
                    std::os::windows::fs::symlink_file(&cache_path, target_path)
                        .with_context(|| format!("Failed to create symlink for {}", description))?;
                }

                println!("    ✓ Linked to cache (saving disk space)");
            } else {
                // If no checksum, we need to copy the file from temp to target
                // This shouldn't happen in practice since all sample data has checksums
                fs::rename(&temp_path, target_path)
                    .or_else(|_| fs::copy(&temp_path, target_path).map(|_| ()))
                    .with_context(|| format!("Failed to move {} to target", description))?;
            }

            // Clean up temp file if it exists
            fs::remove_file(&temp_path).ok();
        }

        // Process split tar files if needed
        let final_aligned_filename = if aligned_urls.len() > 1 {
            // We have multiple parts, need to combine them
            println!("\n  Combining split archive parts...");

            // Determine the base filename (remove .aa, .ab suffixes)
            let first_part = extract_filename_from_url(&aligned_urls[0])?;
            let base_name = if first_part.ends_with(".tar.gz.aa") {
                first_part.trim_end_matches(".aa")
            } else {
                &first_part
            };

            // Create a temp directory for extraction
            let temp_extract_dir =
                std::env::temp_dir().join(format!("bv_extract_{}", Uuid::new_v4()));
            fs::create_dir_all(&temp_extract_dir)
                .context("Failed to create temp extraction dir")?;

            let combined_tar_path = temp_extract_dir.join(base_name);
            let mut combined_file = fs::File::create(&combined_tar_path)
                .context("Failed to create combined tar file")?;

            // Combine all parts
            for url in &aligned_urls {
                let part_filename = extract_filename_from_url(url)?;
                let part_path = participant_dir.join(&part_filename);

                if part_path.exists() {
                    let mut part_file =
                        fs::File::open(&part_path).context("Failed to open part file")?;
                    std::io::copy(&mut part_file, &mut combined_file)
                        .context("Failed to copy part to combined file")?;
                }
            }

            println!("    ✓ Combined archive parts into {}", base_name);

            // Extract the tar.gz to get the final CRAM file
            println!("  Extracting archive...");
            let output = std::process::Command::new("tar")
                .args(["xzf", base_name])
                .current_dir(&temp_extract_dir)
                .output()
                .context("Failed to extract tar archive")?;

            if !output.status.success() {
                anyhow::bail!(
                    "Failed to extract archive: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
            }

            // The extracted file should be the CRAM file
            let cram_name = base_name.trim_end_matches(".tar.gz");
            let extracted_cram_path = temp_extract_dir.join(cram_name);

            // Check if we have the expected checksum for the final CRAM
            // Note: For split files, we might want to add a field for the final CRAM checksum
            // For now, we'll compute it and cache it
            println!("  Computing checksum for extracted CRAM...");
            let cram_hash = crate::cli::download_cache::calculate_blake3(&extracted_cram_path)
                .context("Failed to compute CRAM checksum")?;
            println!("    CRAM blake3: {}", cram_hash);

            // Cache the extracted CRAM file
            let cache_base = crate::config::get_cache_dir()?;
            let cache_path = cache_base.join("by-hash").join(&cram_hash);

            if !cache_path.exists() {
                fs::create_dir_all(cache_path.parent().unwrap())?;
                fs::rename(&extracted_cram_path, &cache_path)
                    .or_else(|_| -> std::io::Result<()> {
                        fs::copy(&extracted_cram_path, &cache_path)?;
                        fs::remove_file(&extracted_cram_path)?;
                        Ok(())
                    })
                    .context("Failed to move CRAM to cache")?;
                println!("    ✓ Cached extracted CRAM file");
            } else {
                println!("    ✓ CRAM already in cache");
                fs::remove_file(&extracted_cram_path).ok();
            }

            // Create symlink in participant directory
            let target_cram_path = participant_dir.join(cram_name);
            if target_cram_path.exists() || target_cram_path.is_symlink() {
                fs::remove_file(&target_cram_path).ok();
            }

            #[cfg(unix)]
            {
                std::os::unix::fs::symlink(&cache_path, target_cram_path)
                    .context("Failed to create symlink for extracted CRAM")?;
            }
            #[cfg(windows)]
            {
                std::os::windows::fs::symlink_file(&cache_path, target_cram_path)
                    .context("Failed to create symlink for extracted CRAM")?;
            }

            // Clean up tar file, parts, and temp directory
            fs::remove_file(&combined_tar_path).ok();
            for url in &aligned_urls {
                let part_filename = extract_filename_from_url(url)?;
                let part_path = participant_dir.join(&part_filename);
                fs::remove_file(&part_path).ok();
            }
            fs::remove_dir_all(&temp_extract_dir).ok();

            println!("    ✓ Extracted and cached {}", cram_name);
            cram_name.to_string()
        } else if !aligned_urls.is_empty() && !aligned_urls[0].is_empty() {
            extract_filename_from_url(&aligned_urls[0])?
        } else {
            String::new()
        };

        let participant_record = ParticipantRecord {
            ref_version: participant_data.ref_version.clone(),
            ref_path: format!("./reference/{}", ref_filename),
            ref_index: format!("./reference/{}", ref_index_filename),
            aligned: if !final_aligned_filename.is_empty() {
                format!("./{}/{}", participant_id, final_aligned_filename)
            } else {
                String::new()
            },
            aligned_index: if !aligned_index_filename.is_empty() {
                format!("./{}/{}", participant_id, aligned_index_filename)
            } else {
                String::new()
            },
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

        // Handle aligned URLs which can be single or multiple
        match &data.aligned {
            AlignedUrl::Single(url) => {
                println!("  aligned: {}", url);
            }
            AlignedUrl::Multiple(urls) => {
                if urls.len() == 1 {
                    println!("  aligned: {}", urls[0]);
                } else {
                    println!("  aligned: {} parts", urls.len());
                    for (i, url) in urls.iter().enumerate() {
                        println!("    part {}: {}", i + 1, url);
                    }
                }
            }
        }

        println!("  aligned_index: {}", data.aligned_index);
    }

    println!("\n{}", "=".repeat(60));
    println!("Use 'bv sample-data fetch <PARTICIPANT_ID>' to download");
    println!("Use 'bv sample-data fetch --all' to download all samples");

    Ok(())
}
