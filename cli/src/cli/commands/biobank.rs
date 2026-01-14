use anyhow::{anyhow, Context, Result};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
use serde_yaml::Value as YamlValue;
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use tracing::debug;

use crate::cli::commands::participant::{Participant, ParticipantsFile};
use crate::cli::syft_url::SyftURL;
use crate::config::{get_config, Config};
use crate::error::Error;
use crate::syftbox::storage::SyftBoxStorage;

const PUBLIC_TEMPLATE_PATH: &str = include_str!("../../participants.public.template.yaml");

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PublicParticipant {
    pub id: String,
    pub url: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_version: Option<String>,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    pub reference: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snp: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mock: Option<YamlValue>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PublicParticipantsFile {
    pub private_url: String,
    #[serde(default)]
    pub datasite: String,
    #[serde(default = "default_http_relay_servers")]
    pub http_relay_servers: Vec<String>,
    #[serde(default)]
    pub public_url: String,
    pub participants: HashMap<String, PublicParticipant>,
    #[serde(flatten)]
    pub mock_data: HashMap<String, YamlValue>,
}

fn default_http_relay_servers() -> Vec<String> {
    vec!["syftbox.net".to_string()]
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftBoxConfig {
    pub data_dir: String,
    #[serde(flatten)]
    pub other: HashMap<String, JsonValue>,
}

fn get_syftbox_config_path() -> Result<PathBuf> {
    let config = get_config()?;
    let syftbox_config_path = config
        .syftbox_config
        .ok_or_else(|| anyhow!("SyftBox config path not found in BioVault config"))?;

    let path = PathBuf::from(&syftbox_config_path);
    if !path.exists() {
        return Err(anyhow!(
            "SyftBox config file not found at: {}",
            syftbox_config_path
        ));
    }

    Ok(path)
}

fn get_syftbox_data_dir() -> Result<PathBuf> {
    let syftbox_config_path = get_syftbox_config_path()?;
    let contents = fs::read_to_string(&syftbox_config_path)
        .with_context(|| format!("Failed to read SyftBox config at {:?}", syftbox_config_path))?;

    let config: SyftBoxConfig =
        serde_json::from_str(&contents).with_context(|| "Failed to parse SyftBox config JSON")?;

    let data_dir = PathBuf::from(&config.data_dir);
    if !data_dir.exists() {
        return Err(anyhow!(
            "SyftBox data directory not found at: {}",
            config.data_dir
        ));
    }

    Ok(data_dir)
}

fn syft_storage() -> Result<SyftBoxStorage> {
    let data_dir = get_syftbox_data_dir()?;
    Ok(SyftBoxStorage::new(&data_dir))
}

fn read_text_with_optional_storage(path: &Path) -> Result<String> {
    if let Ok(data_dir) = get_syftbox_data_dir() {
        let storage = SyftBoxStorage::new(&data_dir);
        if path.starts_with(&data_dir) {
            let bytes = storage
                .read_plaintext_file(path)
                .with_context(|| format!("Failed to read {:?}", path))?;
            return Ok(String::from_utf8(bytes)?);
        }
    }
    fs::read_to_string(path).with_context(|| format!("Failed to read {:?}", path))
}

fn get_datasite_path(email: &str) -> Result<PathBuf> {
    let data_dir = get_syftbox_data_dir()?;
    let datasite_path = data_dir.join("datasites").join(email);

    if !datasite_path.exists() {
        return Err(anyhow!(
            "Datasite path does not exist at: {:?}\nPlease ensure SyftBox is properly configured for email: {}",
            datasite_path,
            email
        ));
    }

    Ok(datasite_path)
}

fn get_public_participants_path(email: &str) -> Result<PathBuf> {
    let datasite_path = get_datasite_path(email)?;
    Ok(datasite_path
        .join("public")
        .join("biovault")
        .join("participants.yaml"))
}

fn load_private_participants() -> Result<ParticipantsFile> {
    ParticipantsFile::load()
}

fn load_public_participants(email: &str) -> Result<PublicParticipantsFile> {
    let path = get_public_participants_path(email)?;
    let storage = syft_storage()?;

    if !path.exists() {
        // Create with defaults if doesn't exist
        let template: YamlValue = serde_yaml::from_str(PUBLIC_TEMPLATE_PATH)?;
        let mock_data = if let YamlValue::Mapping(map) = template {
            map.into_iter()
                .map(|(k, v)| {
                    let key = k.as_str().unwrap_or_default().to_string();
                    (key, v)
                })
                .collect()
        } else {
            HashMap::new()
        };

        Ok(PublicParticipantsFile {
            private_url: format!("syft://{}/private/biovault/participants.yaml", email),
            datasite: email.to_string(),
            http_relay_servers: vec!["syftbox.net".to_string()],
            public_url: format!("syft://{}/public/biovault/participants.yaml", email),
            participants: HashMap::new(),
            mock_data,
        })
    } else {
        let bytes = storage
            .read_plaintext_file(&path)
            .with_context(|| format!("Failed to read public participants file at {:?}", path))?;
        let contents = String::from_utf8(bytes)?;
        let mut parsed: PublicParticipantsFile = serde_yaml::from_str(&contents)
            .with_context(|| "Failed to parse public participants YAML")?;

        // Fill in missing fields for backward compatibility
        if parsed.datasite.is_empty() {
            parsed.datasite = email.to_string();
        }
        if parsed.public_url.is_empty() {
            parsed.public_url = format!("syft://{}/public/biovault/participants.yaml", email);
        }
        // http_relay_servers already has a default via serde

        Ok(parsed)
    }
}

fn save_public_participants(email: &str, file: &PublicParticipantsFile) -> Result<()> {
    let path = get_public_participants_path(email)?;
    let parent = path.parent().ok_or_else(|| anyhow!("Invalid path"))?;
    let storage = syft_storage()?;
    storage
        .ensure_dir(parent)
        .with_context(|| format!("Failed to create directory {:?}", parent))?;

    // Build YAML content manually to support proper anchor references
    let mut yaml_content = String::new();

    // Add datasite
    yaml_content.push_str(&format!("datasite: {}\n", email));

    // Add http_relay_servers
    yaml_content.push_str("http_relay_servers:\n");
    for server in &file.http_relay_servers {
        yaml_content.push_str(&format!("  - {}\n", server));
    }

    // Add public_url
    yaml_content.push_str(&format!(
        "public_url: \"syft://{}/public/biovault/participants.yaml\"\n",
        email
    ));

    // Add private_url
    yaml_content.push_str(&format!(
        "private_url: \"syft://{}/private/biovault/participants.yaml\"\n\n",
        email
    ));

    // Add mock data templates with anchors FIRST (before references)
    yaml_content.push_str(PUBLIC_TEMPLATE_PATH);
    yaml_content.push('\n');

    // Add participants section AFTER anchors are defined
    yaml_content.push_str("participants:\n");
    for (id, participant) in &file.participants {
        yaml_content.push_str(&format!("  {}:\n", id));
        yaml_content.push_str(&format!("    id: {}\n", participant.id));
        yaml_content.push_str(&format!("    url: \"{}\"\n", participant.url));

        // Check if it's an SNP participant
        if let Some(snp) = &participant.snp {
            yaml_content.push_str(&format!("    snp: \"{}\"\n", snp));
            // Add mock field for SNP
            yaml_content.push_str("    mock: *mock_data_snp\n");
        } else {
            // Default CRAM/BAM participant
            if let Some(ref_version) = &participant.ref_version {
                yaml_content.push_str(&format!("    ref_version: {}\n", ref_version));
            }
            if let Some(reference) = &participant.reference {
                yaml_content.push_str(&format!("    ref: \"{}\"\n", reference));
            }
            if let Some(ref_index) = &participant.ref_index {
                yaml_content.push_str(&format!("    ref_index: \"{}\"\n", ref_index));
            }
            if let Some(aligned) = &participant.aligned {
                yaml_content.push_str(&format!("    aligned: \"{}\"\n", aligned));
            }
            if let Some(aligned_index) = &participant.aligned_index {
                yaml_content.push_str(&format!("    aligned_index: \"{}\"\n", aligned_index));
            }
            // Add mock field as a proper YAML anchor reference (no quotes)
            // Only add mock field if we have corresponding mock data in the template
            if let Some(ref_version) = &participant.ref_version {
                if ref_version.to_lowercase() == "grch38" {
                    yaml_content.push_str("    mock: *mock_data_grch38\n");
                }
            }
            // Note: GRCh37 mock data not included yet
        }
    }

    storage
        .write_plaintext_file(&path, yaml_content.as_bytes(), true)
        .with_context(|| format!("Failed to write public participants file at {:?}", path))?;

    Ok(())
}

pub async fn publish(
    participant_id: Option<String>,
    all: bool,
    http_relay_servers: Option<Vec<String>>,
) -> Result<()> {
    let config = get_config()?;
    let email = &config.email;

    println!("{}", "Publishing participants to SyftBox...".green().bold());

    // Validate SyftBox configuration
    let datasite_path = get_datasite_path(email)?;
    debug!("Using datasite path: {:?}", datasite_path);

    let private_participants = load_private_participants()?;
    let mut public_participants = load_public_participants(email)?;

    // Update http_relay_servers if provided
    if let Some(servers) = http_relay_servers {
        public_participants.http_relay_servers = servers;
    }

    let participants_to_publish: Vec<(String, Participant)> = if all {
        private_participants.participants.into_iter().collect()
    } else if let Some(id) = participant_id {
        let participant = private_participants
            .participants
            .get(&id)
            .ok_or_else(|| anyhow!("Participant '{}' not found", id))?;
        vec![(id, participant.clone())]
    } else {
        return Err(anyhow!("Must specify either --participant_id or --all"));
    };

    if participants_to_publish.is_empty() {
        println!("{}", "No participants to publish.".yellow());
        return Ok(());
    }

    for (id, participant) in participants_to_publish {
        let public_participant = if participant.snp.is_some() {
            // SNP participant
            PublicParticipant {
                id: id.clone(),
                url: format!("{{root.private_url}}#participants.{}", id),
                ref_version: None,
                reference: None,
                ref_index: None,
                aligned: None,
                aligned_index: None,
                snp: Some("{url}.snp".to_string()),
                mock: None, // Handled manually in save_public_participants
            }
        } else {
            // Default CRAM/BAM participant
            PublicParticipant {
                id: id.clone(),
                url: format!("{{root.private_url}}#participants.{}", id),
                ref_version: participant.ref_version.clone(),
                reference: Some("{url}.ref".to_string()),
                ref_index: Some("{url}.ref_index".to_string()),
                aligned: Some("{url}.aligned".to_string()),
                aligned_index: Some("{url}.aligned_index".to_string()),
                snp: None,
                mock: None, // Handled manually in save_public_participants
            }
        };

        public_participants
            .participants
            .insert(id.clone(), public_participant);
        println!("✓ Published participant: {}", id.cyan());
    }

    save_public_participants(email, &public_participants)?;

    let path = get_public_participants_path(email)?;
    println!(
        "{}",
        format!("✓ Successfully published to: {:?}", path)
            .green()
            .bold()
    );

    Ok(())
}

pub async fn unpublish(participant_id: Option<String>, all: bool) -> Result<()> {
    let config = get_config()?;
    let email = &config.email;

    println!(
        "{}",
        "Unpublishing participants from SyftBox...".green().bold()
    );

    // Validate SyftBox configuration
    let datasite_path = get_datasite_path(email)?;
    debug!("Using datasite path: {:?}", datasite_path);

    let mut public_participants = load_public_participants(email)?;

    if all {
        let count = public_participants.participants.len();
        public_participants.participants.clear();
        println!("✓ Unpublished {} participants", count);
    } else if let Some(id) = participant_id {
        if public_participants.participants.remove(&id).is_some() {
            println!("✓ Unpublished participant: {}", id.cyan());
        } else {
            return Err(anyhow!("Participant '{}' not found in public file", id));
        }
    } else {
        return Err(anyhow!("Must specify either --participant_id or --all"));
    }

    save_public_participants(email, &public_participants)?;

    let path = get_public_participants_path(email)?;
    println!(
        "{}",
        format!("✓ Successfully updated: {:?}", path).green().bold()
    );

    Ok(())
}

#[derive(Debug, Clone)]
#[allow(dead_code)] // TODO: Will be used when biobank read commands are implemented
pub enum ReadMode {
    Normal,
    Mock,
}

#[allow(dead_code)] // TODO: Will be used when biobank read commands are implemented
pub fn read_participant_with_mode(
    participant_file: &Path,
    participant_id: &str,
    mode: ReadMode,
) -> Result<HashMap<String, String>> {
    let contents = read_text_with_optional_storage(participant_file)
        .with_context(|| format!("Failed to read participant file at {:?}", participant_file))?;

    let yaml: YamlValue =
        serde_yaml::from_str(&contents).with_context(|| "Failed to parse participant YAML")?;

    let mut result = HashMap::new();

    // Extract participant data
    let participant = yaml
        .get("participants")
        .and_then(|p| p.get(participant_id))
        .ok_or_else(|| anyhow!("Participant '{}' not found", participant_id))?;

    match mode {
        ReadMode::Normal => {
            // Return fields as-is with URL references
            if let YamlValue::Mapping(map) = participant {
                for (key, value) in map {
                    if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                        result.insert(k.to_string(), v.to_string());
                    }
                }
            }
        }
        ReadMode::Mock => {
            // Check for mock field and apply mock data if present
            if let Some(mock_value) = participant.get("mock") {
                // Check if it's already dereferenced (anchor resolved to mapping)
                if let YamlValue::Mapping(mock_map) = mock_value {
                    // Use the mock data values directly
                    for (key, value) in mock_map {
                        if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                            result.insert(k.to_string(), v.to_string());
                        }
                    }
                    // Add ID field back
                    if let Some(id) = participant.get("id").and_then(|i| i.as_str()) {
                        result.insert("id".to_string(), id.to_string());
                    }
                } else if let Some(mock_ref) = mock_value.as_str() {
                    // It's a string reference like "*mock_data_grch38", look it up
                    let mock_key = mock_ref.trim_start_matches('*');

                    if let Some(YamlValue::Mapping(mock_map)) = yaml.get(mock_key) {
                        for (key, value) in mock_map {
                            if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                                result.insert(k.to_string(), v.to_string());
                            }
                        }
                        // Add ID field back
                        if let Some(id) = participant.get("id").and_then(|i| i.as_str()) {
                            result.insert("id".to_string(), id.to_string());
                        }
                    } else {
                        // Couldn't find mock data, return normal fields
                        if let YamlValue::Mapping(map) = participant {
                            for (key, value) in map {
                                if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                                    result.insert(k.to_string(), v.to_string());
                                }
                            }
                        }
                    }
                } else {
                    // Unexpected mock value type, return normal fields
                    if let YamlValue::Mapping(map) = participant {
                        for (key, value) in map {
                            if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                                result.insert(k.to_string(), v.to_string());
                            }
                        }
                    }
                }
            } else {
                // No mock data, return normal
                if let YamlValue::Mapping(map) = participant {
                    for (key, value) in map {
                        if let (Some(k), Some(v)) = (key.as_str(), value.as_str()) {
                            result.insert(k.to_string(), v.to_string());
                        }
                    }
                }
            }
        }
    }

    Ok(result)
}

// Biobank list command structures and implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiobankFile {
    #[serde(default)]
    pub datasite: String,
    #[serde(default = "default_http_relay_servers")]
    pub http_relay_servers: Vec<String>,
    #[serde(default)]
    pub public_url: String,
    pub participants: BTreeMap<String, BiobankParticipant>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiobankParticipant {
    pub id: String,
    pub url: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_version: Option<String>,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    pub reference: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snp: Option<String>,
}

#[derive(Debug)]
pub struct BioBankInfo {
    pub email: String,
    pub public_url: String,
    pub http_relay_servers: Vec<String>,
    pub participants: BTreeMap<String, BiobankParticipant>,
}

pub async fn list(override_path: Option<PathBuf>) -> crate::error::Result<()> {
    let config = Config::load()?;

    let data_dir = if let Some(path) = override_path {
        path
    } else {
        config.get_syftbox_data_dir()?
    };
    let storage = SyftBoxStorage::new(&data_dir);

    let datasites_dir = data_dir.join("datasites");

    if !datasites_dir.exists() {
        return Err(Error::DatasitesDirMissing(
            datasites_dir.display().to_string(),
        ));
    }

    let biobanks = find_biobanks(&storage, &datasites_dir)?;

    if biobanks.is_empty() {
        println!("No biobanks found in {}", datasites_dir.display());
        return Ok(());
    }

    println!("BioBanks (sorted alphabetically)\n");

    for biobank in biobanks {
        println!("Datasite: {}", biobank.email);
        println!("Syft URL: {}", biobank.public_url);
        println!("HTTP Relay URLs:");

        for server in &biobank.http_relay_servers {
            if let Ok(syft_url) = SyftURL::parse(&biobank.public_url) {
                let http_url = syft_url.to_http_relay_url(server);
                println!("  - {}", http_url);
            } else {
                println!("  - Error: Invalid Syft URL format");
            }
        }

        println!("\n\nNumber of Participants: {}", biobank.participants.len());
        println!("Participants:\n");

        for (participant_id, participant) in &biobank.participants {
            let participant_type = if participant.snp.is_some() {
                "SNP".to_string()
            } else if let Some(ref ref_version) = participant.ref_version {
                ref_version.clone()
            } else {
                "Unknown".to_string()
            };
            println!("{} ({})", participant_id, participant_type);

            // Try to parse the participant URL
            if let Ok(syft_url) = SyftURL::parse(&participant.url) {
                println!("Syft URL: {}", syft_url);
            } else {
                // If parsing fails, construct from public_url
                if let Ok(base_url) = SyftURL::parse(&biobank.public_url) {
                    let participant_url =
                        base_url.with_fragment(format!("participants.{}", participant_id));
                    println!("Syft URL: {}", participant_url);
                } else {
                    println!("Syft URL: Error - Invalid URL format");
                }
            }
            println!();
        }

        println!("--------------------------------");
    }

    Ok(())
}

fn find_biobanks(storage: &SyftBoxStorage, datasites_dir: &Path) -> Result<Vec<BioBankInfo>> {
    let mut biobanks = Vec::new();
    let entries = storage.list_dir(datasites_dir)?;
    for entry in entries {
        if !entry.is_dir() {
            continue;
        }
        let Some(email) = entry.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        let participants_path = entry
            .join("public")
            .join("biovault")
            .join("participants.yaml");
        if !participants_path.exists() {
            continue;
        }

        match load_biobank_file(storage, &participants_path, email) {
            Ok(info) => biobanks.push(info),
            Err(e) => {
                eprintln!(
                    "Invalid file at path: {} ({})",
                    participants_path.display(),
                    e
                );
            }
        }
    }

    biobanks.sort_by(|a, b| a.email.cmp(&b.email));

    Ok(biobanks)
}

fn load_biobank_file(storage: &SyftBoxStorage, path: &Path, email: &str) -> Result<BioBankInfo> {
    let content = if storage.contains(path) {
        let bytes = storage
            .read_plaintext_file(path)
            .with_context(|| format!("Failed to read participants file: {}", path.display()))?;
        String::from_utf8(bytes).with_context(|| format!("Invalid UTF-8 in {}", path.display()))?
    } else {
        fs::read_to_string(path)
            .with_context(|| format!("Failed to read participants file: {}", path.display()))?
    };

    let mut biobank_file: BiobankFile = serde_yaml::from_str(&content)
        .with_context(|| format!("Failed to parse participants file: {}", path.display()))?;

    // Fill in missing fields for backward compatibility
    if biobank_file.datasite.is_empty() {
        biobank_file.datasite = email.to_string();
    }
    if biobank_file.public_url.is_empty() {
        biobank_file.public_url = format!("syft://{}/public/biovault/participants.yaml", email);
    }

    // Ensure participants have proper URL format with dot notation
    for (id, participant) in biobank_file.participants.iter_mut() {
        if participant.url.is_empty() || participant.url.contains("#participants/") {
            // Convert old slash format to dot notation or set default
            participant.url = format!(
                "syft://{}/public/biovault/participants.yaml#participants.{}",
                email, id
            );
        }
    }

    Ok(BioBankInfo {
        email: biobank_file.datasite.clone(),
        public_url: biobank_file.public_url,
        http_relay_servers: biobank_file.http_relay_servers,
        participants: biobank_file.participants,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::syftbox::SyftBoxApp;
    use tempfile::TempDir;

    // Helper function to create a test config
    fn create_test_config(syftbox_data_dir: &Path, email: &str, test_home: &Path) -> Result<()> {
        let config_dir = test_home.join(".biovault");
        fs::create_dir_all(&config_dir)?;

        let syftbox_config_path = syftbox_data_dir.join("syftbox_config.json");
        let syftbox_config = SyftBoxConfig {
            data_dir: syftbox_data_dir.to_str().unwrap().to_string(),
            other: HashMap::new(),
        };
        fs::write(
            &syftbox_config_path,
            serde_json::to_string(&syftbox_config)?,
        )?;

        let config = Config {
            email: email.to_string(),
            syftbox_config: Some(syftbox_config_path.to_str().unwrap().to_string()),
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let config_path = config_dir.join("config.yaml");
        config.save(&config_path)?;

        Ok(())
    }

    // Helper to create test participants file
    fn create_test_participants_file(test_home: &Path) -> Result<()> {
        let participants_dir = test_home.join(".biovault");
        fs::create_dir_all(&participants_dir)?;

        let participants = ParticipantsFile {
            participants: vec![
                (
                    "TEST1".to_string(),
                    Participant {
                        id: "TEST1".to_string(),
                        ref_version: Some("GRCh38".to_string()),
                        r#ref: Some("/data/ref/grch38.fa".to_string()),
                        ref_index: Some("/data/ref/grch38.fa.fai".to_string()),
                        aligned: Some("/data/aligned/test1.cram".to_string()),
                        aligned_index: Some("/data/aligned/test1.cram.crai".to_string()),
                        snp: None,
                    },
                ),
                (
                    "TEST2".to_string(),
                    Participant {
                        id: "TEST2".to_string(),
                        ref_version: Some("GRCh38".to_string()),
                        r#ref: Some("/data/ref/grch38_2.fa".to_string()),
                        ref_index: Some("/data/ref/grch38_2.fa.fai".to_string()),
                        aligned: Some("/data/aligned/test2.bam".to_string()),
                        aligned_index: Some("/data/aligned/test2.bam.bai".to_string()),
                        snp: None,
                    },
                ),
            ]
            .into_iter()
            .collect(),
        };

        let participants_path = participants_dir.join("participants.yaml");
        fs::write(&participants_path, serde_yaml::to_string(&participants)?)?;

        Ok(())
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_publish_and_unpublish() {
        let temp_dir = TempDir::new().unwrap();
        let email = "test@example.com";

        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));
        let syftbox_config_path = temp_dir.path().join("syftbox_config.json");
        crate::config::set_test_config(crate::config::Config {
            email: email.to_string(),
            syftbox_config: Some(syftbox_config_path.to_string_lossy().to_string()),
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        });

        let app = SyftBoxApp::new(temp_dir.path(), email, "biovault").unwrap();
        let data_dir = app.data_dir.clone();

        create_test_config(temp_dir.path(), email, temp_dir.path()).unwrap();
        create_test_participants_file(temp_dir.path()).unwrap();

        let read_public = |path: &Path| -> String {
            String::from_utf8(
                app.storage
                    .read_plaintext_file(path)
                    .expect("read participants"),
            )
            .expect("utf8")
        };

        let public_path = data_dir
            .join("datasites")
            .join(email)
            .join("public")
            .join("biovault")
            .join("participants.yaml");

        let result = publish(Some("TEST1".to_string()), false, None).await;
        assert!(result.is_ok());

        let mut public_content = read_public(&public_path);
        assert!(public_content.contains("TEST1"));
        assert!(public_content.contains("mock: *mock_data_grch38"));

        let result = publish(None, true, None).await;
        assert!(result.is_ok());

        public_content = read_public(&public_path);
        assert!(public_content.contains("TEST1"));
        assert!(public_content.contains("TEST2"));

        let custom_servers = vec![
            "relay1.example.com".to_string(),
            "relay2.example.com".to_string(),
        ];
        let result = publish(Some("TEST1".to_string()), false, Some(custom_servers)).await;
        assert!(result.is_ok());

        public_content = read_public(&public_path);
        assert!(public_content
            .contains("http_relay_servers:\n  - relay1.example.com\n  - relay2.example.com"));

        let result = unpublish(Some("TEST1".to_string()), false).await;
        assert!(result.is_ok());

        public_content = read_public(&public_path);
        assert!(!public_content.contains("TEST1:"));
        assert!(public_content.contains("TEST2"));

        let result = unpublish(None, true).await;
        assert!(result.is_ok());

        public_content = read_public(&public_path);
        assert!(!public_content.contains("TEST1:"));
        assert!(!public_content.contains("TEST2:"));

        crate::config::clear_test_config();
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_read_participant_with_published_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("public_participants.yaml");

        // Create a published participants file with mock data
        // Note: YAML anchors must be defined before they are referenced
        let yaml_content = r#"private_url: "syft://test@example.com/private/biovault/participants.yaml"

# Mock data definitions (anchors must be defined before references)
mock_data_grch38: &mock_data_grch38
  ref_version: GRCh38
  ref: https://example.com/grch38.fa
  ref_index: https://example.com/grch38.fa.fai
  aligned: https://example.com/test.cram
  aligned_index: https://example.com/test.cram.crai

participants:
  TEST1:
    id: TEST1
    url: "{root.private_url}#participants.TEST1"
    ref_version: GRCh38
    ref: "{url}.ref"
    ref_index: "{url}.ref_index"
    aligned: "{url}.aligned"
    aligned_index: "{url}.aligned_index"
    mock: *mock_data_grch38
"#;
        fs::write(&file_path, yaml_content).unwrap();

        // Test normal mode - should return template references
        let result = read_participant_with_mode(&file_path, "TEST1", ReadMode::Normal).unwrap();
        assert_eq!(result.get("id").unwrap(), "TEST1");
        assert_eq!(result.get("ref_version").unwrap(), "GRCh38");
        assert_eq!(result.get("ref").unwrap(), "{url}.ref");
        assert_eq!(
            result.get("url").unwrap(),
            "{root.private_url}#participants.TEST1"
        );

        // Test mock mode - should return mock data values
        let result = read_participant_with_mode(&file_path, "TEST1", ReadMode::Mock).unwrap();
        assert_eq!(result.get("id").unwrap(), "TEST1");
        assert_eq!(result.get("ref_version").unwrap(), "GRCh38");
        // When there's a mock field, it should substitute with the mock data values
        assert_eq!(result.get("ref").unwrap(), "https://example.com/grch38.fa");
        assert_eq!(
            result.get("aligned").unwrap(),
            "https://example.com/test.cram"
        );
    }

    #[test]
    fn test_read_participant_normal_mode() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test.yaml");

        let yaml_content = r#"
participants:
  TEST:
    id: TEST
    url: "{root.private_url}#participants.TEST"
    ref_version: GRCh38
    ref: "{url}.ref"
    ref_index: "{url}.ref_index"
    aligned: "{url}.aligned"
    aligned_index: "{url}.aligned_index"
"#;
        fs::write(&file_path, yaml_content).unwrap();

        let result = read_participant_with_mode(&file_path, "TEST", ReadMode::Normal).unwrap();

        assert_eq!(result.get("id").unwrap(), "TEST");
        assert_eq!(result.get("ref_version").unwrap(), "GRCh38");
        assert_eq!(result.get("ref").unwrap(), "{url}.ref");
        assert_eq!(
            result.get("url").unwrap(),
            "{root.private_url}#participants.TEST"
        );
    }

    #[test]
    fn test_read_participant_mock_mode() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test.yaml");

        let yaml_content = r#"
participants:
  TEST:
    id: TEST
    url: "{root.private_url}#participants.TEST"
    ref_version: GRCh38
    ref: "{url}.ref"
    mock: "*mock_data_grch38"

mock_data_grch38:
  ref_version: GRCh38
  ref: https://example.com/reference.fa
  ref_index: https://example.com/reference.fa.fai
  aligned: https://example.com/aligned.cram
  aligned_index: https://example.com/aligned.cram.crai
"#;
        fs::write(&file_path, yaml_content).unwrap();

        let result = read_participant_with_mode(&file_path, "TEST", ReadMode::Mock).unwrap();

        assert_eq!(result.get("id").unwrap(), "TEST");
        assert_eq!(result.get("ref_version").unwrap(), "GRCh38");
        assert_eq!(
            result.get("ref").unwrap(),
            "https://example.com/reference.fa"
        );
        assert_eq!(
            result.get("aligned").unwrap(),
            "https://example.com/aligned.cram"
        );
    }

    #[test]
    fn test_read_participant_mock_mode_no_mock_field() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test.yaml");

        let yaml_content = r#"
participants:
  TEST:
    id: TEST
    ref_version: GRCh38
    ref: "/path/to/reference.fa"
    aligned: "/path/to/aligned.cram"
"#;
        fs::write(&file_path, yaml_content).unwrap();

        let result = read_participant_with_mode(&file_path, "TEST", ReadMode::Mock).unwrap();

        assert_eq!(result.get("id").unwrap(), "TEST");
        assert_eq!(result.get("ref_version").unwrap(), "GRCh38");
        assert_eq!(result.get("ref").unwrap(), "/path/to/reference.fa");
        assert_eq!(result.get("aligned").unwrap(), "/path/to/aligned.cram");
    }

    #[tokio::test]
    #[serial_test::serial]
    async fn test_list_biobanks() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let storage = SyftBoxStorage::new(temp_dir.path());
        let datasites_dir = temp_dir.path().join("datasites");
        storage.ensure_dir(&datasites_dir)?;

        create_test_biobank(
            &storage,
            &datasites_dir,
            "alice@example.com",
            vec![("PARTICIPANT1", "GRCh38"), ("PARTICIPANT2", "GRCh37")],
        )?;

        create_test_biobank(
            &storage,
            &datasites_dir,
            "bob@example.org",
            vec![("TEST", "GRCh38")],
        )?;

        create_test_biobank(
            &storage,
            &datasites_dir,
            "charlie@example.net",
            vec![
                ("SAMPLE1", "GRCh38"),
                ("SAMPLE2", "GRCh38"),
                ("SAMPLE3", "GRCh37"),
            ],
        )?;

        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault"));

        let config_dir = temp_dir.path().join(".biovault");
        fs::create_dir_all(&config_dir)?;

        let config = Config {
            email: "test@example.com".to_string(),
            syftbox_config: Some(
                temp_dir
                    .path()
                    .join(".syftbox/config.json")
                    .to_string_lossy()
                    .to_string(),
            ),
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        config.save(config_dir.join("config.yaml"))?;

        let syftbox_dir = temp_dir.path().join(".syftbox");
        fs::create_dir_all(&syftbox_dir)?;

        let syftbox_config = serde_json::json!({
            "data_dir": temp_dir.path().to_string_lossy()
        });
        fs::write(
            syftbox_dir.join("config.json"),
            serde_json::to_string_pretty(&syftbox_config)?,
        )?;

        let result = list(Some(temp_dir.path().to_path_buf())).await;
        assert!(result.is_ok());

        crate::config::clear_test_biovault_home();

        Ok(())
    }

    fn create_test_biobank(
        storage: &SyftBoxStorage,
        datasites_dir: &Path,
        email: &str,
        participants: Vec<(&str, &str)>,
    ) -> Result<()> {
        let biobank_dir = datasites_dir
            .join("example.com")
            .join(email)
            .join("public")
            .join("biovault");
        storage.ensure_dir(&biobank_dir)?;

        let mut participants_map = BTreeMap::new();
        for (id, ref_version) in participants {
            participants_map.insert(
                id.to_string(),
                BiobankParticipant {
                    id: id.to_string(),
                    url: format!(
                        "syft://{}/public/biovault/participants.yaml#participants.{}",
                        email, id
                    ),
                    ref_version: Some(ref_version.to_string()),
                    reference: None,
                    ref_index: None,
                    aligned: None,
                    aligned_index: None,
                    snp: None,
                },
            );
        }

        let biobank_file = BiobankFile {
            datasite: email.to_string(),
            http_relay_servers: vec!["syftbox.net".to_string()],
            public_url: format!("syft://{}/public/biovault/participants.yaml", email),
            participants: participants_map,
        };

        let yaml_content = serde_yaml::to_string(&biobank_file)?;
        storage.write_plaintext_file(
            &biobank_dir.join("participants.yaml"),
            yaml_content.as_bytes(),
            true,
        )?;

        Ok(())
    }
}
