use anyhow::{anyhow, Context, Result};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
use serde_yaml::Value as YamlValue;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use tracing::debug;

use crate::cli::commands::participant::{Participant, ParticipantsFile};
use crate::config::get_config;

const PUBLIC_TEMPLATE_PATH: &str = include_str!("../../participants.public.template.yaml");

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PublicParticipant {
    pub id: String,
    pub url: String,
    pub ref_version: String,
    #[serde(rename = "ref")]
    pub reference: String,
    pub ref_index: String,
    pub aligned: String,
    pub aligned_index: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mock: Option<YamlValue>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PublicParticipantsFile {
    pub private_url: String,
    pub participants: HashMap<String, PublicParticipant>,
    #[serde(flatten)]
    pub mock_data: HashMap<String, YamlValue>,
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
            participants: HashMap::new(),
            mock_data,
        })
    } else {
        let contents = fs::read_to_string(&path)
            .with_context(|| format!("Failed to read public participants file at {:?}", path))?;
        let parsed: PublicParticipantsFile = serde_yaml::from_str(&contents)
            .with_context(|| "Failed to parse public participants YAML")?;
        Ok(parsed)
    }
}

fn save_public_participants(email: &str, file: &PublicParticipantsFile) -> Result<()> {
    let path = get_public_participants_path(email)?;
    let parent = path.parent().ok_or_else(|| anyhow!("Invalid path"))?;
    fs::create_dir_all(parent)
        .with_context(|| format!("Failed to create directory {:?}", parent))?;

    // Build YAML content manually to support proper anchor references
    let mut yaml_content = String::new();

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
        yaml_content.push_str(&format!("    ref_version: {}\n", participant.ref_version));
        yaml_content.push_str(&format!("    ref: \"{}\"\n", participant.reference));
        yaml_content.push_str(&format!("    ref_index: \"{}\"\n", participant.ref_index));
        yaml_content.push_str(&format!("    aligned: \"{}\"\n", participant.aligned));
        yaml_content.push_str(&format!(
            "    aligned_index: \"{}\"\n",
            participant.aligned_index
        ));

        // Add mock field as a proper YAML anchor reference (no quotes)
        // Only add mock field if we have corresponding mock data in the template
        if participant.ref_version.to_lowercase() == "grch38" {
            yaml_content.push_str("    mock: *mock_data_grch38\n");
        }
        // Note: GRCh37 mock data not included yet
    }

    fs::write(&path, yaml_content)
        .with_context(|| format!("Failed to write public participants file at {:?}", path))?;

    Ok(())
}

pub async fn publish(participant_id: Option<String>, all: bool) -> Result<()> {
    let config = get_config()?;
    let email = &config.email;

    println!("{}", "Publishing participants to SyftBox...".green().bold());

    // Validate SyftBox configuration
    let datasite_path = get_datasite_path(email)?;
    debug!("Using datasite path: {:?}", datasite_path);

    let private_participants = load_private_participants()?;
    let mut public_participants = load_public_participants(email)?;

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
        let public_participant = PublicParticipant {
            id: id.clone(),
            url: format!("{{root.private_url}}#participants.{}", id),
            ref_version: participant.ref_version.clone(),
            reference: "{url}.ref".to_string(),
            ref_index: "{url}.ref_index".to_string(),
            aligned: "{url}.aligned".to_string(),
            aligned_index: "{url}.aligned_index".to_string(),
            mock: None, // Handled manually in save_public_participants
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
    let contents = fs::read_to_string(participant_file)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use tempfile::TempDir;

    // Helper function to create a test config
    fn create_test_config(syftbox_data_dir: &Path, email: &str) -> Result<()> {
        let config_dir = dirs::home_dir()
            .ok_or_else(|| anyhow!("Could not determine home directory"))?
            .join(".biovault");
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
        };

        let config_path = config_dir.join("config.yaml");
        fs::write(&config_path, serde_yaml::to_string(&config)?)?;

        Ok(())
    }

    // Helper to create test participants file
    fn create_test_participants_file() -> Result<()> {
        let home = dirs::home_dir().ok_or_else(|| anyhow!("Could not determine home directory"))?;
        let participants_dir = home.join(".biovault");
        fs::create_dir_all(&participants_dir)?;

        let participants = ParticipantsFile {
            participants: vec![
                (
                    "TEST1".to_string(),
                    Participant {
                        id: "TEST1".to_string(),
                        ref_version: "GRCh38".to_string(),
                        r#ref: "/data/ref/grch38.fa".to_string(),
                        ref_index: "/data/ref/grch38.fa.fai".to_string(),
                        aligned: "/data/aligned/test1.cram".to_string(),
                        aligned_index: "/data/aligned/test1.cram.crai".to_string(),
                    },
                ),
                (
                    "TEST2".to_string(),
                    Participant {
                        id: "TEST2".to_string(),
                        ref_version: "GRCh38".to_string(),
                        r#ref: "/data/ref/grch38_2.fa".to_string(),
                        ref_index: "/data/ref/grch38_2.fa.fai".to_string(),
                        aligned: "/data/aligned/test2.bam".to_string(),
                        aligned_index: "/data/aligned/test2.bam.bai".to_string(),
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
    async fn test_publish_and_unpublish() {
        let temp_dir = TempDir::new().unwrap();
        let email = "test@example.com";

        // Setup test environment
        let datasites_dir = temp_dir.path().join("datasites").join(email);
        fs::create_dir_all(&datasites_dir).unwrap();

        create_test_config(temp_dir.path(), email).unwrap();
        create_test_participants_file().unwrap();

        // Test publishing a single participant
        let result = publish(Some("TEST1".to_string()), false).await;
        assert!(result.is_ok());

        // Check that the public file was created
        let public_path = datasites_dir
            .join("public")
            .join("biovault")
            .join("participants.yaml");
        assert!(public_path.exists());

        // Read and verify the public file content
        let public_content = fs::read_to_string(&public_path).unwrap();
        assert!(public_content.contains("TEST1"));
        assert!(public_content.contains("mock: *mock_data_grch38"));
        assert!(public_content.contains("mock_data_grch38: &mock_data_grch38"));
        assert!(public_content.contains("{root.private_url}#participants.TEST1"));
        assert!(public_content.contains("{url}.ref"));

        // Test publishing all participants
        let result = publish(None, true).await;
        assert!(result.is_ok());

        let public_content = fs::read_to_string(&public_path).unwrap();
        assert!(public_content.contains("TEST1"));
        assert!(public_content.contains("TEST2"));
        // Both TEST1 and TEST2 use GRCh38, so both should have mock references
        assert!(public_content.contains("mock: *mock_data_grch38"));
        assert!(public_content.contains("mock_data_grch38: &mock_data_grch38"));

        // Test unpublishing a single participant
        let result = unpublish(Some("TEST1".to_string()), false).await;
        assert!(result.is_ok());

        let public_content = fs::read_to_string(&public_path).unwrap();
        assert!(!public_content.contains("TEST1:"));
        assert!(public_content.contains("TEST2"));

        // Test unpublishing all
        let result = unpublish(None, true).await;
        assert!(result.is_ok());

        let public_content = fs::read_to_string(&public_path).unwrap();
        assert!(!public_content.contains("TEST1:"));
        assert!(!public_content.contains("TEST2:"));

        // Clean up
        let config_dir = dirs::home_dir().unwrap().join(".biovault");
        let _ = fs::remove_file(config_dir.join("config.yaml"));
        let _ = fs::remove_file(config_dir.join("participants.yaml"));
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
}
