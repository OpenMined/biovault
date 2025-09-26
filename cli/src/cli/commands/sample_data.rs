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
struct PostProcess {
    #[serde(default)]
    uncompress: Option<bool>,
    #[serde(default)]
    file: Option<String>,
    #[serde(default)]
    extract: Option<String>,
    #[serde(default)]
    rename: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ParticipantData {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    ref_version: Option<String>,
    #[serde(rename = "ref", default, skip_serializing_if = "Option::is_none")]
    ref_url: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    ref_index: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    aligned: Option<AlignedUrl>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    aligned_index: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    ref_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    ref_index_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    aligned_b3sum: Option<AlignedChecksum>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    aligned_index_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    snp: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    snp_b3sum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    snp_post_process: Option<PostProcess>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ParticipantsFile {
    participant: HashMap<String, ParticipantRecord>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ParticipantRecord {
    #[serde(skip_serializing_if = "Option::is_none")]
    ref_version: Option<String>,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    ref_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ref_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    aligned: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    aligned_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    snp: Option<String>,
}

pub async fn fetch(
    participant_ids: Option<Vec<String>>,
    all: bool,
    quiet: bool,
) -> anyhow::Result<()> {
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

    if !quiet {
        println!(
            "Fetching sample data for {} participant(s):",
            participants_to_fetch.len()
        );
        for id in &participants_to_fetch {
            println!("  - {}", id);
        }
        println!();
    }

    let participants_file_path = sample_data_dir.join("participants.yaml");
    let mut participants_file = load_or_create_participants_file(&participants_file_path)?;

    // Initialize download cache
    let mut download_cache = DownloadCache::new(None)?;

    for participant_id in participants_to_fetch {
        if !quiet {
            println!("\n{}", "=".repeat(60));
            println!("Fetching data for participant: {}", participant_id);
            println!("{}", "=".repeat(60));
        }

        let participant_data = config
            .sample_data_urls
            .get(&participant_id)
            .ok_or_else(|| {
                anyhow::anyhow!("Participant {} not found in configuration", participant_id)
            })?;

        let participant_dir = sample_data_dir.join(&participant_id);
        fs::create_dir_all(&participant_dir).context("Failed to create participant directory")?;

        // Check if this is SNP data or CRAM data
        if let Some(snp_url) = &participant_data.snp {
            // Handle SNP data
            let snp_filename = extract_filename_from_url(snp_url)?;
            let snp_checksum = participant_data
                .snp_b3sum
                .as_ref()
                .unwrap_or(&String::new())
                .clone();

            let downloads = vec![(
                snp_url.clone(),
                participant_dir.join(&snp_filename),
                "SNP data".to_string(),
                snp_checksum,
            )];

            for (url, target_path, description, expected_b3sum) in &downloads {
                if !quiet {
                    println!(
                        "\n  Processing {}: {}",
                        description,
                        target_path.file_name().unwrap().to_string_lossy()
                    );
                }

                // Set up download options based on whether we have a checksum
                let mut options = if !expected_b3sum.is_empty() {
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
                options.show_progress = !quiet;

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

                // Handle archive extraction if post_process specifies it
                let should_uncompress = participant_data
                    .snp_post_process
                    .as_ref()
                    .and_then(|pp| pp.uncompress)
                    .unwrap_or(false)
                    || snp_filename.ends_with(".zip");

                if should_uncompress {
                    if !quiet {
                        println!("  Extracting SNP archive...");
                    }

                    // Create temp extraction directory
                    let temp_extract_dir =
                        std::env::temp_dir().join(format!("bv_snp_extract_{}", Uuid::new_v4()));
                    fs::create_dir_all(&temp_extract_dir)
                        .context("Failed to create temp extraction dir")?;

                    // Extract the zip file
                    let output = std::process::Command::new("unzip")
                        .args(["-q", "-o", temp_path.to_str().unwrap()])
                        .current_dir(&temp_extract_dir)
                        .output()
                        .context("Failed to extract SNP zip archive")?;

                    if !output.status.success() {
                        anyhow::bail!(
                            "Failed to extract SNP archive: {}",
                            String::from_utf8_lossy(&output.stderr)
                        );
                    }

                    // Move extracted files to participant directory
                    let entries = fs::read_dir(&temp_extract_dir)
                        .context("Failed to read extracted files")?;
                    for entry in entries {
                        let entry = entry?;
                        let file_name = entry.file_name();
                        let source_path = entry.path();
                        let target_file_path = participant_dir.join(&file_name);

                        // Move or copy file
                        fs::rename(&source_path, &target_file_path)
                            .or_else(|_| -> std::io::Result<()> {
                                fs::copy(&source_path, &target_file_path)?;
                                fs::remove_file(&source_path)?;
                                Ok(())
                            })
                            .with_context(|| {
                                format!("Failed to move extracted file: {:?}", file_name)
                            })?;

                        if !quiet {
                            println!("    ✓ Extracted: {}", file_name.to_string_lossy());
                        }
                    }

                    // Clean up
                    fs::remove_dir_all(&temp_extract_dir).ok();
                    fs::remove_file(&temp_path).ok();
                } else if !expected_b3sum.is_empty() {
                    // Non-zip file with checksum - create symlink to cache
                    let cache_base = crate::config::get_cache_dir()?;
                    let cache_path = cache_base.join("by-hash").join(expected_b3sum);

                    // Remove any existing file or symlink at target
                    if target_path.exists() || target_path.is_symlink() {
                        fs::remove_file(target_path).ok();
                    }

                    // Create symlink to cache
                    #[cfg(unix)]
                    {
                        std::os::unix::fs::symlink(&cache_path, target_path).with_context(
                            || format!("Failed to create symlink for {}", description),
                        )?;
                    }
                    #[cfg(windows)]
                    {
                        std::os::windows::fs::symlink_file(&cache_path, target_path).with_context(
                            || format!("Failed to create symlink for {}", description),
                        )?;
                    }

                    if !quiet {
                        println!("    ✓ Linked to cache (saving disk space)");
                    }

                    // Clean up temp file
                    fs::remove_file(&temp_path).ok();
                } else {
                    // No checksum, copy file to target
                    fs::rename(&temp_path, target_path)
                        .or_else(|_| fs::copy(&temp_path, target_path).map(|_| ()))
                        .with_context(|| format!("Failed to move {} to target", description))?;
                    fs::remove_file(&temp_path).ok();
                }
            }

            // Update participants file for SNP data
            // Include the specific file if specified in post_process
            let snp_path = if let Some(ref post_process) = participant_data.snp_post_process {
                if let Some(ref file) = post_process.file {
                    format!("./{}/{}", participant_id, file)
                } else {
                    format!("./{}", participant_id)
                }
            } else {
                format!("./{}", participant_id)
            };

            let participant_record = ParticipantRecord {
                ref_version: None,
                ref_path: None,
                ref_index: None,
                aligned: None,
                aligned_index: None,
                snp: Some(snp_path),
            };

            participants_file
                .participant
                .insert(participant_id.clone(), participant_record);

            save_participants_file(&participants_file_path, &participants_file)?;
            if !quiet {
                println!("  ✓ Updated participants.yaml");
            }

            continue; // Skip the rest of the CRAM processing
        }

        // Extract filenames from URLs for CRAM data
        let ref_filename = extract_filename_from_url(
            participant_data
                .ref_url
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("Missing ref URL for CRAM data"))?,
        )?;
        let ref_index_filename = extract_filename_from_url(
            participant_data
                .ref_index
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("Missing ref_index URL for CRAM data"))?,
        )?;

        // Handle aligned URLs (can be single or multiple)
        let aligned_urls = participant_data
            .aligned
            .as_ref()
            .map(|a| a.to_vec())
            .unwrap_or_default();
        let aligned_checksums = participant_data
            .aligned_b3sum
            .as_ref()
            .map(|a| a.to_vec())
            .unwrap_or_default();
        let aligned_index_filename = if let Some(aligned_index) = &participant_data.aligned_index {
            if !aligned_index.is_empty() {
                extract_filename_from_url(aligned_index)?
            } else {
                String::new()
            }
        } else {
            String::new()
        };

        // Build downloads list - first add reference files
        let mut downloads = vec![
            (
                participant_data.ref_url.clone().unwrap_or_default(),
                reference_dir.join(&ref_filename),
                "Reference genome".to_string(),
                participant_data.ref_b3sum.clone().unwrap_or_default(),
            ),
            (
                participant_data.ref_index.clone().unwrap_or_default(),
                reference_dir.join(&ref_index_filename),
                "Reference index".to_string(),
                participant_data.ref_index_b3sum.clone().unwrap_or_default(),
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
        if let Some(aligned_index) = &participant_data.aligned_index {
            if !aligned_index.is_empty() {
                downloads.push((
                    aligned_index.clone(),
                    participant_dir.join(&aligned_index_filename),
                    "CRAM index".to_string(),
                    participant_data
                        .aligned_index_b3sum
                        .clone()
                        .unwrap_or_default(),
                ));
            }
        };

        for (url, target_path, description, expected_b3sum) in &downloads {
            if !quiet {
                println!(
                    "\n  Processing {}: {}",
                    description,
                    target_path.file_name().unwrap().to_string_lossy()
                );
            }

            // Set up download options based on whether we have a checksum
            let mut options = if !expected_b3sum.is_empty() {
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
            options.show_progress = !quiet;

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

                if !quiet {
                    println!("    ✓ Linked to cache (saving disk space)");
                }
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
            if !quiet {
                println!("\n  Combining split archive parts...");
            }

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

            if !quiet {
                println!("    ✓ Combined archive parts into {}", base_name);
            }

            // Extract the tar.gz to get the final CRAM file
            if !quiet {
                println!("  Extracting archive...");
            }
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
            if !quiet {
                println!("  Computing checksum for extracted CRAM...");
            }
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
                if !quiet {
                    println!("    ✓ Cached extracted CRAM file");
                }
            } else {
                if !quiet {
                    println!("    ✓ CRAM already in cache");
                }
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
            ref_path: Some(format!("./reference/{}", ref_filename)),
            ref_index: Some(format!("./reference/{}", ref_index_filename)),
            aligned: if !final_aligned_filename.is_empty() {
                Some(format!("./{}/{}", participant_id, final_aligned_filename))
            } else {
                None
            },
            aligned_index: if !aligned_index_filename.is_empty() {
                Some(format!("./{}/{}", participant_id, aligned_index_filename))
            } else {
                None
            },
            snp: None,
        };

        participants_file
            .participant
            .insert(participant_id.clone(), participant_record);

        save_participants_file(&participants_file_path, &participants_file)?;
        if !quiet {
            println!("  ✓ Updated participants.yaml");
        }
    }

    if !quiet {
        println!("\n{}", "=".repeat(60));
        println!("✓ Sample data fetch complete!");
        println!("  Data location: {}", sample_data_dir.display());
        println!("  Participants file: {}", participants_file_path.display());
        println!("{}", "=".repeat(60));
    }

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

        // Check if this is SNP data or CRAM data
        if let Some(snp_url) = &data.snp {
            println!("  Type: SNP");
            println!("  snp: {}", snp_url);
            if let Some(snp_b3sum) = &data.snp_b3sum {
                println!("  snp_b3sum: {}", snp_b3sum);
            }
        } else {
            // CRAM data
            if let Some(ref_version) = &data.ref_version {
                println!("  ref_version: {}", ref_version);
            }
            if let Some(ref_url) = &data.ref_url {
                println!("  ref: {}", ref_url);
            }
            if let Some(ref_index) = &data.ref_index {
                println!("  ref_index: {}", ref_index);
            }

            // Handle aligned URLs which can be single or multiple
            if let Some(aligned) = &data.aligned {
                match aligned {
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
            }

            if let Some(aligned_index) = &data.aligned_index {
                println!("  aligned_index: {}", aligned_index);
            }
        }
    }

    println!("\n{}", "=".repeat(60));
    println!("Use 'bv sample-data fetch <PARTICIPANT_ID>' to download");
    println!("Use 'bv sample-data fetch --all' to download all samples");

    Ok(())
}

#[cfg(test)]
mod list_tests {
    use super::*;

    #[tokio::test]
    async fn sample_list_runs() {
        // Should print available sample info without error
        list().await.unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn filename_extraction_handles_fragments_and_queries() {
        assert_eq!(
            extract_filename_from_url("https://x/y/file.txt").unwrap(),
            "file.txt"
        );
        assert_eq!(
            extract_filename_from_url("https://x/y/file.txt?download=1").unwrap(),
            "file.txt"
        );
        assert_eq!(
            extract_filename_from_url("https://x/y/file.txt#frag").unwrap(),
            "file.txt"
        );
    }

    #[test]
    fn determine_participants_selection_rules() {
        let mut map = std::collections::HashMap::new();
        map.insert(
            "A".to_string(),
            ParticipantData {
                ref_version: None,
                ref_url: None,
                ref_index: None,
                aligned: None,
                aligned_index: None,
                ref_b3sum: None,
                ref_index_b3sum: None,
                aligned_b3sum: None,
                aligned_index_b3sum: None,
                snp: None,
                snp_b3sum: None,
                snp_post_process: None,
            },
        );
        let cfg = SampleDataConfig {
            sample_data_urls: map,
        };
        // all=true picks all keys
        let mut ids = determine_participants_to_fetch(&cfg, None, true).unwrap();
        ids.sort();
        assert_eq!(ids, vec!["A".to_string()]);

        // explicit list must exist
        let ids2 = determine_participants_to_fetch(&cfg, Some(vec!["A".into()]), false).unwrap();
        assert_eq!(ids2, vec!["A".to_string()]);
        // non-existent yields error
        assert!(determine_participants_to_fetch(&cfg, Some(vec!["B".into()]), false).is_err());
        // none provided and all=false -> error
        assert!(determine_participants_to_fetch(&cfg, None, false).is_err());
    }

    #[test]
    fn load_or_create_and_save_participants_file() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("participants.yaml");
        // create new empty
        let pf = load_or_create_participants_file(&p).unwrap();
        assert!(pf.participant.is_empty());
        // save and load back
        save_participants_file(&p, &pf).unwrap();
        let pf2 = load_or_create_participants_file(&p).unwrap();
        assert!(pf2.participant.is_empty());
    }

    #[test]
    fn aligned_url_and_checksum_to_vec() {
        let u1 = AlignedUrl::Single("http://x/file.cram".into());
        assert_eq!(u1.to_vec(), vec!["http://x/file.cram".to_string()]);
        let u2 = AlignedUrl::Multiple(vec!["u1".into(), "u2".into()]);
        assert_eq!(u2.to_vec(), vec!["u1".to_string(), "u2".to_string()]);

        let c1 = AlignedChecksum::Single("abc".into());
        assert_eq!(c1.to_vec(), vec!["abc".to_string()]);
        let c2 = AlignedChecksum::Multiple(vec!["a".into(), "b".into()]);
        assert_eq!(c2.to_vec(), vec!["a".to_string(), "b".to_string()]);

        // Default is a single empty string
        let def = AlignedChecksum::default();
        assert_eq!(def.to_vec(), vec![String::new()]);
    }
}
