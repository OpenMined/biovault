use crate::config::Config;
use crate::Result;
use anyhow::Context;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use tracing::debug;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Participants {
    pub participants: BTreeMap<String, Participant>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Participant {
    pub ref_version: String,
    #[serde(rename = "ref")]
    pub reference: Option<String>,
    pub ref_index: Option<String>,
    pub aligned: Option<String>,
    pub aligned_index: Option<String>,
}

#[derive(Debug)]
pub struct BioBankInfo {
    pub email: String,
    pub participants: Participants,
}

pub async fn list(override_path: Option<PathBuf>) -> Result<()> {
    let config = Config::load()?;

    let data_dir = if let Some(path) = override_path {
        path
    } else {
        config.get_syftbox_data_dir()?
    };

    let datasites_dir = data_dir.join("datasites");

    if !datasites_dir.exists() {
        return Err(anyhow::anyhow!(
            "Datasites directory not found at: {}",
            datasites_dir.display()
        )
        .into());
    }

    let biobanks = find_biobanks(&datasites_dir)?;

    if biobanks.is_empty() {
        println!("No biobanks found in {}", datasites_dir.display());
        return Ok(());
    }

    println!("BioBanks (sorted alphabetically)\n");

    for biobank in biobanks {
        println!("--------------------------------");
        println!("{}", biobank.email);
        println!(
            "Number of Participants: {}",
            biobank.participants.participants.len()
        );
        println!("\nParticipants:");

        for (participant_id, participant) in &biobank.participants.participants {
            println!("{} ({})", participant_id, participant.ref_version);
            println!(
                "url: syft://{}/private/biovault/participants.yaml#participants/{}",
                biobank.email, participant_id
            );
            println!();
        }
    }

    println!("--------------------------------");

    Ok(())
}

fn find_biobanks(datasites_dir: &Path) -> Result<Vec<BioBankInfo>> {
    let mut biobanks = Vec::new();

    let pattern = format!(
        "{}/**/**/public/biobank/participants.yaml",
        datasites_dir.display()
    );
    debug!("Searching for biobank files with pattern: {}", pattern);

    for entry in glob::glob(&pattern)? {
        match entry {
            Ok(path) => {
                debug!("Found participants file: {}", path.display());

                if let Some(email) = extract_email_from_path(&path, datasites_dir) {
                    match load_participants(&path) {
                        Ok(participants) => {
                            biobanks.push(BioBankInfo {
                                email,
                                participants,
                            });
                        }
                        Err(e) => {
                            eprintln!(
                                "Warning: Failed to load participants from {}: {}",
                                path.display(),
                                e
                            );
                        }
                    }
                } else {
                    eprintln!(
                        "Warning: Could not extract email from path: {}",
                        path.display()
                    );
                }
            }
            Err(e) => {
                eprintln!("Warning: Error accessing path: {}", e);
            }
        }
    }

    biobanks.sort_by(|a, b| a.email.cmp(&b.email));

    Ok(biobanks)
}

fn extract_email_from_path(path: &Path, datasites_dir: &Path) -> Option<String> {
    let relative_path = path.strip_prefix(datasites_dir).ok()?;
    let components: Vec<&str> = relative_path
        .components()
        .filter_map(|c| c.as_os_str().to_str())
        .collect();

    if components.len() >= 2 {
        Some(components[1].to_string())
    } else {
        None
    }
}

fn load_participants(path: &Path) -> Result<Participants> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Failed to read participants file: {}", path.display()))?;

    let participants: Participants = serde_yaml::from_str(&content)
        .with_context(|| format!("Failed to parse participants file: {}", path.display()))?;

    Ok(participants)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[tokio::test]
    async fn test_list_biobanks() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let datasites_dir = temp_dir.path().join("datasites");

        create_test_biobank(
            &datasites_dir,
            "alice@example.com",
            vec![("PATIENT1", "GRCh38"), ("PATIENT2", "GRCh37")],
        )?;

        create_test_biobank(&datasites_dir, "bob@example.org", vec![("TEST", "GRCh38")])?;

        create_test_biobank(
            &datasites_dir,
            "charlie@example.net",
            vec![
                ("SAMPLE1", "GRCh38"),
                ("SAMPLE2", "GRCh38"),
                ("SAMPLE3", "GRCh37"),
            ],
        )?;

        std::env::set_var("BIOVAULT_TEST_HOME", temp_dir.path());

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

        std::env::remove_var("BIOVAULT_TEST_HOME");

        Ok(())
    }

    fn create_test_biobank(
        datasites_dir: &Path,
        email: &str,
        participants: Vec<(&str, &str)>,
    ) -> Result<()> {
        let biobank_dir = datasites_dir
            .join("example.com")
            .join(email)
            .join("public")
            .join("biobank");
        fs::create_dir_all(&biobank_dir)?;

        let mut participants_map = BTreeMap::new();
        for (id, ref_version) in participants {
            participants_map.insert(
                id.to_string(),
                Participant {
                    ref_version: ref_version.to_string(),
                    reference: None,
                    ref_index: None,
                    aligned: None,
                    aligned_index: None,
                },
            );
        }

        let participants_data = Participants {
            participants: participants_map,
        };

        let yaml_content = serde_yaml::to_string(&participants_data)?;
        fs::write(biobank_dir.join("participants.yaml"), yaml_content)?;

        Ok(())
    }
}
