use crate::Result;
use anyhow::Context;
use blake3;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use tokio::fs::File;
use tokio::io::AsyncWriteExt;

const SAMPLE_DATA_YAML: &str = include_str!("../../sample_data.yaml");

#[derive(Debug, Serialize, Deserialize)]
struct SampleDataConfig {
    sample_data_urls: HashMap<String, PatientData>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct PatientData {
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
struct PatientsFile {
    patient: HashMap<String, PatientRecord>,
}

#[derive(Debug, Serialize, Deserialize)]
struct PatientRecord {
    ref_version: String,
    #[serde(rename = "ref")]
    ref_path: String,
    ref_index: String,
    aligned: String,
    aligned_index: String,
}

pub async fn fetch(patient_ids: Option<Vec<String>>, all: bool) -> Result<()> {
    let config: SampleDataConfig = serde_yaml::from_str(SAMPLE_DATA_YAML)
        .context("Failed to parse embedded sample data configuration")?;

    let patients_to_fetch = determine_patients_to_fetch(&config, patient_ids, all)?;

    if patients_to_fetch.is_empty() {
        println!("No patients specified to fetch");
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
        "Fetching sample data for {} patient(s):",
        patients_to_fetch.len()
    );
    for id in &patients_to_fetch {
        println!("  - {}", id);
    }
    println!();

    let patients_file_path = sample_data_dir.join("patients.yaml");
    let mut patients_file = load_or_create_patients_file(&patients_file_path)?;

    for patient_id in patients_to_fetch {
        println!("\n{}", "=".repeat(60));
        println!("Fetching data for patient: {}", patient_id);
        println!("{}", "=".repeat(60));

        let patient_data = config
            .sample_data_urls
            .get(&patient_id)
            .ok_or_else(|| anyhow::anyhow!("Patient {} not found in configuration", patient_id))?;

        let patient_dir = sample_data_dir.join(&patient_id);
        fs::create_dir_all(&patient_dir).context("Failed to create patient directory")?;

        // Extract filenames from URLs
        let ref_filename = extract_filename_from_url(&patient_data.ref_url)?;
        let ref_index_filename = extract_filename_from_url(&patient_data.ref_index)?;
        let aligned_filename = extract_filename_from_url(&patient_data.aligned)?;
        let aligned_index_filename = extract_filename_from_url(&patient_data.aligned_index)?;

        let downloads = vec![
            (
                &patient_data.ref_url,
                reference_dir.join(&ref_filename),
                "Reference genome",
                &patient_data.ref_b3sum,
            ),
            (
                &patient_data.ref_index,
                reference_dir.join(&ref_index_filename),
                "Reference index",
                &patient_data.ref_index_b3sum,
            ),
            (
                &patient_data.aligned,
                patient_dir.join(&aligned_filename),
                "Aligned CRAM",
                &patient_data.aligned_b3sum,
            ),
            (
                &patient_data.aligned_index,
                patient_dir.join(&aligned_index_filename),
                "CRAM index",
                &patient_data.aligned_index_b3sum,
            ),
        ];

        for (url, target_path, description, expected_b3sum) in downloads {
            if target_path.exists() {
                println!("  ✓ {} already exists at:", description);
                println!("    {}", target_path.display());
                if !expected_b3sum.is_empty() {
                    print!("    Verifying BLAKE3 checksum... ");
                    std::io::stdout().flush().unwrap();
                    match verify_file_checksum(&target_path, expected_b3sum) {
                        Ok(true) => println!("✓ Valid"),
                        Ok(false) => {
                            println!("✗ Invalid checksum!");
                            println!("    Expected: {}", expected_b3sum);
                            if let Ok(actual) = calculate_blake3(&target_path) {
                                println!("    Actual:   {}", actual);
                            }
                            return Err(anyhow::anyhow!(
                                "Checksum verification failed for {}",
                                description
                            )
                            .into());
                        }
                        Err(e) => {
                            println!("⚠ Could not verify: {}", e);
                        }
                    }
                }
            } else {
                println!("  ↓ Downloading {} to:", description);
                println!("    {}", target_path.display());
                download_file(url, &target_path)
                    .await
                    .with_context(|| format!("Failed to download {}", description))?;
                println!("    ✓ Download complete");

                if !expected_b3sum.is_empty() {
                    print!("    Verifying BLAKE3 checksum... ");
                    std::io::stdout().flush().unwrap();
                    match verify_file_checksum(&target_path, expected_b3sum) {
                        Ok(true) => println!("✓ Valid"),
                        Ok(false) => {
                            println!("✗ Invalid checksum!");
                            println!("    Expected: {}", expected_b3sum);
                            if let Ok(actual) = calculate_blake3(&target_path) {
                                println!("    Actual:   {}", actual);
                            }
                            // Delete the corrupted file
                            let _ = fs::remove_file(&target_path);
                            return Err(anyhow::anyhow!(
                                "Downloaded file has invalid checksum for {}",
                                description
                            )
                            .into());
                        }
                        Err(e) => {
                            println!("⚠ Could not verify: {}", e);
                        }
                    }
                }
            }
        }

        let patient_record = PatientRecord {
            ref_version: patient_data.ref_version.clone(),
            ref_path: format!("./reference/{}", ref_filename),
            ref_index: format!("./reference/{}", ref_index_filename),
            aligned: format!("./{}/{}", patient_id, aligned_filename),
            aligned_index: format!("./{}/{}", patient_id, aligned_index_filename),
        };

        patients_file
            .patient
            .insert(patient_id.clone(), patient_record);

        save_patients_file(&patients_file_path, &patients_file)?;
        println!("  ✓ Updated patients.yaml");
    }

    println!("\n{}", "=".repeat(60));
    println!("✓ Sample data fetch complete!");
    println!("  Data location: {}", sample_data_dir.display());
    println!("  Patients file: {}", patients_file_path.display());
    println!("{}", "=".repeat(60));

    Ok(())
}

fn extract_filename_from_url(url: &str) -> Result<String> {
    url.split('/')
        .last()
        .ok_or_else(|| anyhow::anyhow!("Invalid URL format: {}", url))
        .map(|s| s.to_string())
        .map_err(Into::into)
}

fn calculate_blake3(path: &Path) -> Result<String> {
    // For large files, blake3::Hasher::update_rayon provides parallel hashing
    let mut file = fs::File::open(path)
        .with_context(|| format!("Failed to open file for checksum: {}", path.display()))?;

    // Get file size to decide strategy
    let metadata = file.metadata()?;
    let file_size = metadata.len();

    if file_size > 100 * 1024 * 1024 {
        // For files > 100MB, use parallel hashing
        let mut hasher = blake3::Hasher::new();
        let mut buffer = vec![0; 64 * 1024 * 1024]; // 64MB buffer for parallel processing

        loop {
            let bytes_read = file
                .read(&mut buffer)
                .with_context(|| format!("Failed to read file for checksum: {}", path.display()))?;

            if bytes_read == 0 {
                break;
            }

            // Use rayon-enabled parallel update for large chunks
            hasher.update_rayon(&buffer[..bytes_read]);
        }

        Ok(hasher.finalize().to_hex().to_string())
    } else {
        // For smaller files, use regular update
        let mut hasher = blake3::Hasher::new();
        let mut buffer = vec![0; 8 * 1024 * 1024]; // 8MB buffer

        loop {
            let bytes_read = file
                .read(&mut buffer)
                .with_context(|| format!("Failed to read file for checksum: {}", path.display()))?;

            if bytes_read == 0 {
                break;
            }

            hasher.update(&buffer[..bytes_read]);
        }

        Ok(hasher.finalize().to_hex().to_string())
    }
}

fn verify_file_checksum(path: &Path, expected_b3sum: &str) -> Result<bool> {
    let actual = calculate_blake3(path)?;
    Ok(actual == expected_b3sum)
}

fn determine_patients_to_fetch(
    config: &SampleDataConfig,
    patient_ids: Option<Vec<String>>,
    all: bool,
) -> Result<Vec<String>> {
    if all {
        Ok(config.sample_data_urls.keys().cloned().collect())
    } else if let Some(ids) = patient_ids {
        for id in &ids {
            if !config.sample_data_urls.contains_key(id) {
                return Err(anyhow::anyhow!(
                    "Patient {} not found in sample data configuration",
                    id
                )
                .into());
            }
        }
        Ok(ids)
    } else {
        Err(anyhow::anyhow!("No patients specified. Use patient IDs or --all flag").into())
    }
}

fn load_or_create_patients_file(path: &Path) -> Result<PatientsFile> {
    if path.exists() {
        let content = fs::read_to_string(path).context("Failed to read existing patients.yaml")?;
        Ok(serde_yaml::from_str(&content).context("Failed to parse existing patients.yaml")?)
    } else {
        Ok(PatientsFile {
            patient: HashMap::new(),
        })
    }
}

fn save_patients_file(path: &Path, patients: &PatientsFile) -> Result<()> {
    let yaml = serde_yaml::to_string(patients).context("Failed to serialize patients data")?;
    fs::write(path, yaml).context("Failed to write patients.yaml")?;
    Ok(())
}

async fn download_file(url: &str, target_path: &Path) -> Result<()> {
    let client = reqwest::Client::builder()
        .timeout(std::time::Duration::from_secs(3600))
        .build()
        .context("Failed to create HTTP client")?;

    let response = client
        .get(url)
        .send()
        .await
        .context("Failed to send request")?;

    if !response.status().is_success() {
        return Err(
            anyhow::anyhow!("HTTP request failed with status: {}", response.status()).into(),
        );
    }

    let total_size = response.content_length().unwrap_or(0);

    let pb =
        if total_size > 0 {
            let pb = ProgressBar::new(total_size);
            pb.set_style(ProgressStyle::default_bar()
            .template("    [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
            .expect("Failed to set progress bar template")
            .progress_chars("#>-"));
            Some(pb)
        } else {
            println!("    Downloading (size unknown)...");
            None
        };

    let mut file = File::create(target_path)
        .await
        .context("Failed to create target file")?;

    let mut downloaded = 0u64;
    let mut stream = response.bytes_stream();

    while let Some(chunk) = futures_util::StreamExt::next(&mut stream).await {
        let chunk = chunk.context("Failed to read chunk")?;
        file.write_all(&chunk)
            .await
            .context("Failed to write chunk to file")?;

        downloaded += chunk.len() as u64;
        if let Some(ref pb) = pb {
            pb.set_position(downloaded);
        }
    }

    if let Some(pb) = pb {
        pb.finish_and_clear();
    }

    Ok(())
}

pub async fn list() -> Result<()> {
    let config: SampleDataConfig = serde_yaml::from_str(SAMPLE_DATA_YAML)
        .context("Failed to parse embedded sample data configuration")?;

    println!("Available sample data:");
    println!("{}", "=".repeat(60));

    for (patient_id, data) in &config.sample_data_urls {
        println!("\nPatient ID: {}", patient_id);
        println!("  ref_version: {}", data.ref_version);
        println!("  ref: {}", data.ref_url);
        println!("  ref_index: {}", data.ref_index);
        println!("  aligned: {}", data.aligned);
        println!("  aligned_index: {}", data.aligned_index);
    }

    println!("\n{}", "=".repeat(60));
    println!("Use 'bv sample-data fetch <PATIENT_ID>' to download");
    println!("Use 'bv sample-data fetch --all' to download all samples");

    Ok(())
}
