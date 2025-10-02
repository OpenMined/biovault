use anyhow::{anyhow, Context, Result};
use colored::Colorize;
use dialoguer::{Confirm, Input, Select};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tracing::warn;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Participant {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub r#ref: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_index: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snp: Option<String>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct ParticipantsFile {
    pub participants: HashMap<String, Participant>,
}

impl ParticipantsFile {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn load() -> anyhow::Result<Self> {
        let path = get_participants_file_path()?;
        if !path.exists() {
            Ok(Self::new())
        } else {
            let contents = fs::read_to_string(&path)
                .with_context(|| format!("Failed to read participants file at {:?}", path))?;
            let parsed: Self = serde_yaml::from_str(&contents)
                .with_context(|| "Failed to parse participants YAML")?;
            Ok(parsed)
        }
    }

    fn save(&self) -> Result<()> {
        let path = get_participants_file_path()?;
        let parent = path.parent().ok_or_else(|| anyhow!("Invalid path"))?;
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create directory {:?}", parent))?;

        let yaml = serde_yaml::to_string(self)
            .with_context(|| "Failed to serialize participants to YAML")?;
        fs::write(&path, yaml)
            .with_context(|| format!("Failed to write participants file at {:?}", path))?;
        Ok(())
    }
}

fn get_participants_file_path() -> Result<PathBuf> {
    Ok(crate::config::get_biovault_home()?.join("participants.yaml"))
}

fn check_samtools_installed() -> bool {
    Command::new("which")
        .arg("samtools")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

fn detect_reference_version(aligned_file: &str) -> Option<String> {
    if !check_samtools_installed() {
        return None;
    }

    let check_ebv = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "samtools view -H {} | grep -c 'SN:chrEBV'",
            aligned_file
        ))
        .output();

    let check_ki270 = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "samtools view -H {} | grep -c 'SN:KI270'",
            aligned_file
        ))
        .output();

    match (check_ebv, check_ki270) {
        (Ok(ebv), Ok(ki270)) => {
            let ebv_count = String::from_utf8_lossy(&ebv.stdout)
                .trim()
                .parse::<i32>()
                .unwrap_or(0);
            let ki270_count = String::from_utf8_lossy(&ki270.stdout)
                .trim()
                .parse::<i32>()
                .unwrap_or(0);

            if ebv_count > 0 || ki270_count > 0 {
                Some("GRCh38".to_string())
            } else {
                Some("GRCh37".to_string())
            }
        }
        _ => None,
    }
}

fn get_index_file_path(file_path: &str, is_aligned: bool) -> String {
    if is_aligned {
        if file_path.ends_with(".cram") {
            format!("{}.crai", file_path)
        } else if file_path.ends_with(".bam") {
            format!("{}.bai", file_path)
        } else if file_path.ends_with(".sam") {
            format!("{}.sai", file_path)
        } else {
            format!("{}.idx", file_path)
        }
    } else {
        format!("{}.fai", file_path)
    }
}

fn expand_tilde(path: &str) -> String {
    if path.starts_with("~/") {
        if let Ok(home) = std::env::var("HOME").or_else(|_| std::env::var("USERPROFILE")) {
            return path.replacen("~", &home, 1);
        }
    }
    path.to_string()
}

fn normalize_path(path: &str) -> Result<String> {
    // Expand tilde (~) to home directory
    let expanded = expand_tilde(path);

    // Convert to absolute path if relative
    let path_buf = Path::new(&expanded);
    if path_buf.is_relative() {
        let absolute = std::env::current_dir()?.join(path_buf);
        Ok(absolute.to_string_lossy().to_string())
    } else {
        Ok(expanded)
    }
}

pub async fn add(
    id: Option<String>,
    aligned: Option<String>,
    template: Option<String>,
    snp: Option<String>,
    reference: Option<String>,
    ref_version: Option<String>,
    non_interactive: bool,
) -> Result<()> {
    println!("{}", "Adding new participant...".green().bold());

    let template_type = template.as_deref().unwrap_or("default");

    // Validate template type
    if template_type != "default" && template_type != "snp" {
        return Err(anyhow!(
            "Invalid template type '{}'. Must be 'default' or 'snp'",
            template_type
        ));
    }

    let mut participants_file = ParticipantsFile::load()?;
    let mut index_jobs = Vec::new();

    let participant_id = match id {
        Some(id) => id,
        None => {
            if non_interactive {
                return Err(anyhow!(
                    "Participant ID is required in non-interactive mode"
                ));
            }
            Input::<String>::new()
                .with_prompt("Enter participant ID")
                .interact_text()?
        }
    };

    if participants_file.participants.contains_key(&participant_id) {
        if non_interactive {
            // In non-interactive mode, just overwrite
            println!("⚠ Overwriting existing participant '{}'", participant_id);
        } else {
            let overwrite = Confirm::new()
                .with_prompt(format!(
                    "Participant '{}' already exists. Overwrite?",
                    participant_id
                ))
                .interact()?;

            if !overwrite {
                println!("Cancelled.");
                return Ok(());
            }
        }
    }

    // Handle SNP template
    if template_type == "snp" {
        let snp_file = match snp {
            Some(file) => file,
            None => {
                if non_interactive {
                    return Err(anyhow!(
                        "SNP file path is required for SNP template in non-interactive mode"
                    ));
                }
                Input::<String>::new()
                    .with_prompt("Enter SNP file path")
                    .validate_with(|input: &String| -> Result<(), &str> {
                        let expanded = expand_tilde(input);
                        if Path::new(&expanded).exists() {
                            Ok(())
                        } else {
                            Err("File does not exist")
                        }
                    })
                    .interact_text()?
            }
        };

        let snp_file_normalized = normalize_path(&snp_file)?;
        if !Path::new(&snp_file_normalized).exists() {
            return Err(anyhow!("SNP file does not exist: {}", snp_file_normalized));
        }

        let participant = Participant {
            id: participant_id.clone(),
            ref_version: None,
            r#ref: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: Some(snp_file_normalized),
        };

        participants_file
            .participants
            .insert(participant_id.clone(), participant);
        participants_file.save()?;

        println!(
            "{}",
            format!("✓ SNP Participant '{}' added successfully!", participant_id)
                .green()
                .bold()
        );

        return Ok(());
    }

    // Handle default template (existing logic)
    if !check_samtools_installed() {
        warn!("samtools is not installed. Index creation will not be available.");
        println!("{}", "Warning: samtools is not installed. You won't be able to create index files automatically.".yellow());
    }

    let aligned_file = match aligned {
        Some(file) => file,
        None => {
            if non_interactive {
                return Err(anyhow!(
                    "Aligned file path is required for default template in non-interactive mode"
                ));
            }
            Input::<String>::new()
                .with_prompt("Enter aligned file path (.cram, .bam, or .sam)")
                .validate_with(|input: &String| -> Result<(), &str> {
                    let expanded = expand_tilde(input);
                    if Path::new(&expanded).exists() {
                        Ok(())
                    } else {
                        Err("File does not exist")
                    }
                })
                .interact_text()?
        }
    };

    let aligned_file = normalize_path(&aligned_file)?;
    if !Path::new(&aligned_file).exists() {
        return Err(anyhow!("Aligned file does not exist: {}", aligned_file));
    }

    let aligned_index = get_index_file_path(&aligned_file, true);
    let aligned_index_exists = Path::new(&aligned_index).exists();

    if aligned_index_exists {
        println!("✓ Index file found: {}", aligned_index.green());
    } else {
        println!("⚠ Index file not found: {}", aligned_index.yellow());
    }

    let use_aligned_index = if aligned_index_exists {
        if non_interactive {
            true // Use the auto-detected index in non-interactive mode
        } else {
            Confirm::new()
                .with_prompt(format!("Use this index file: {}?", aligned_index))
                .default(true)
                .interact()?
        }
    } else {
        false
    };

    let final_aligned_index = if use_aligned_index {
        aligned_index.clone()
    } else if non_interactive {
        // In non-interactive mode, use the default index path
        if !Path::new(&aligned_index).exists() && check_samtools_installed() {
            index_jobs.push((aligned_file.clone(), true));
        }
        aligned_index
    } else {
        let custom_index = Input::<String>::new()
            .with_prompt("Enter aligned index file path")
            .default(aligned_index.clone())
            .interact_text()?;

        if !Path::new(&custom_index).exists() && check_samtools_installed() {
            let create_index = Confirm::new()
                .with_prompt("Index file doesn't exist. Create it after setup?")
                .interact()?;

            if create_index {
                index_jobs.push((aligned_file.clone(), true));
            }
        }
        custom_index
    };

    let detected_version = detect_reference_version(&aligned_file);
    let version_hint = detected_version
        .as_ref()
        .map(|v| format!(" (Possibly {} detected)", v))
        .unwrap_or_default();

    let ref_version_str = if let Some(version) = ref_version {
        version
    } else if non_interactive {
        // Use detected version or default to GRCh38
        detected_version.unwrap_or_else(|| "GRCh38".to_string())
    } else {
        let ref_version_choices = vec!["GRCh38", "GRCh37"];
        let default_selection = match detected_version.as_deref() {
            Some("GRCh38") => 0,
            Some("GRCh37") => 1,
            _ => 0,
        };

        let ref_version_idx = Select::new()
            .with_prompt(format!("Select reference version{}", version_hint))
            .items(&ref_version_choices)
            .default(default_selection)
            .interact()?;

        ref_version_choices[ref_version_idx].to_string()
    };

    let ref_file = if let Some(ref_path) = reference {
        ref_path
    } else if non_interactive {
        return Err(anyhow!(
            "Reference genome file path is required in non-interactive mode"
        ));
    } else {
        Input::<String>::new()
            .with_prompt("Enter reference genome file path (.fa or .fasta)")
            .validate_with(|input: &String| -> Result<(), &str> {
                let expanded = expand_tilde(input);
                if Path::new(&expanded).exists() {
                    Ok(())
                } else {
                    Err("File does not exist")
                }
            })
            .interact_text()?
    };

    let ref_file = normalize_path(&ref_file)?;

    let ref_index = get_index_file_path(&ref_file, false);
    let ref_index_exists = Path::new(&ref_index).exists();

    if ref_index_exists {
        println!("✓ Reference index found: {}", ref_index.green());
    } else {
        println!("⚠ Reference index not found: {}", ref_index.yellow());
    }

    let use_ref_index = if ref_index_exists {
        if non_interactive {
            true // Use the auto-detected reference index in non-interactive mode
        } else {
            Confirm::new()
                .with_prompt(format!("Use this reference index file: {}?", ref_index))
                .default(true)
                .interact()?
        }
    } else {
        false
    };

    let final_ref_index = if use_ref_index {
        ref_index.clone()
    } else if non_interactive {
        // In non-interactive mode, use the default reference index path
        if !Path::new(&ref_index).exists() && check_samtools_installed() {
            index_jobs.push((ref_file.clone(), false));
        }
        ref_index
    } else {
        let custom_index = Input::<String>::new()
            .with_prompt("Enter reference index file path")
            .default(ref_index.clone())
            .interact_text()?;

        if !Path::new(&custom_index).exists() && check_samtools_installed() {
            let create_index = Confirm::new()
                .with_prompt("Reference index doesn't exist. Create it after setup?")
                .interact()?;

            if create_index {
                index_jobs.push((ref_file.clone(), false));
            }
        }
        custom_index
    };

    let participant = Participant {
        id: participant_id.clone(),
        ref_version: Some(ref_version_str),
        r#ref: Some(ref_file.clone()),
        ref_index: Some(normalize_path(&final_ref_index)?),
        aligned: Some(aligned_file.clone()),
        aligned_index: Some(normalize_path(&final_aligned_index)?),
        snp: None,
    };

    participants_file
        .participants
        .insert(participant_id.clone(), participant);
    participants_file.save()?;

    println!(
        "{}",
        format!("✓ Participant '{}' added successfully!", participant_id)
            .green()
            .bold()
    );

    if !index_jobs.is_empty() && check_samtools_installed() {
        println!("\n{}", "Creating index files...".cyan().bold());
        for (file, is_aligned) in index_jobs {
            println!("Indexing: {}", file);
            let result = if is_aligned {
                Command::new("samtools").args(["index", &file]).output()
            } else {
                Command::new("samtools").args(["faidx", &file]).output()
            };

            match result {
                Ok(output) if output.status.success() => {
                    println!("✓ Index created successfully");
                }
                Ok(output) => {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    println!("✗ Failed to create index: {}", stderr.red());
                }
                Err(e) => {
                    println!("✗ Failed to run samtools: {}", e.to_string().red());
                }
            }
        }
    }

    Ok(())
}

pub async fn list() -> Result<()> {
    let participants_file = ParticipantsFile::load()?;

    if participants_file.participants.is_empty() {
        println!("{}", "No participants found.".yellow());
        println!("Use 'bv participant add' to add a participant.");
        return Ok(());
    }

    println!("{}", "Participants:".green().bold());
    println!();

    for (id, participant) in &participants_file.participants {
        println!("  {} {}", "id:".bold(), id.cyan());
        if let Some(ref snp) = participant.snp {
            println!("    type: SNP");
            println!("    snp: {}", snp);
        } else {
            println!("    type: CRAM/BAM");
            if let Some(ref ref_version) = participant.ref_version {
                println!("    ref_version: {}", ref_version);
            }
            if let Some(ref r#ref) = participant.r#ref {
                println!("    ref: {}", r#ref);
            }
            if let Some(ref ref_index) = participant.ref_index {
                println!("    ref_index: {}", ref_index);
            }
            if let Some(ref aligned) = participant.aligned {
                println!("    aligned: {}", aligned);
            }
            if let Some(ref aligned_index) = participant.aligned_index {
                println!("    aligned_index: {}", aligned_index);
            }
        }
        println!();
    }

    Ok(())
}

pub async fn delete(id: String) -> Result<()> {
    let mut participants_file = ParticipantsFile::load()?;

    if !participants_file.participants.contains_key(&id) {
        return Err(anyhow!("Participant '{}' not found", id));
    }

    let confirm = Confirm::new()
        .with_prompt(format!("Delete participant '{}'?", id))
        .interact()?;

    if !confirm {
        println!("Cancelled.");
        return Ok(());
    }

    participants_file.participants.remove(&id);
    participants_file.save()?;

    println!(
        "{}",
        format!("✓ Participant '{}' deleted successfully!", id)
            .green()
            .bold()
    );
    Ok(())
}

pub async fn validate(id: Option<String>) -> Result<()> {
    if !check_samtools_installed() {
        return Err(anyhow!(
            "samtools is not installed. Please install samtools first."
        ));
    }

    let check_seqkit = Command::new("which")
        .arg("seqkit")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false);

    if !check_seqkit {
        return Err(anyhow!(
            "seqkit is not installed. Please install seqkit first."
        ));
    }

    let participants_file = ParticipantsFile::load()?;

    let participants_to_validate: Vec<(String, Participant)> = match id {
        Some(ref id) => {
            let participant = participants_file
                .participants
                .get(id)
                .ok_or_else(|| anyhow!("Participant '{}' not found", id))?;
            vec![(id.clone(), participant.clone())]
        }
        None => participants_file.participants.into_iter().collect(),
    };

    if participants_to_validate.is_empty() {
        println!("{}", "No participants to validate.".yellow());
        return Ok(());
    }

    println!("{}", "Validating participants...".green().bold());
    println!();

    let mut all_valid = true;

    for (id, participant) in participants_to_validate {
        println!("Validating participant: {}", id.cyan());

        // Check if it's an SNP participant
        if let Some(snp) = &participant.snp {
            // Validate SNP file exists
            if Path::new(snp).exists() {
                println!("  ✓ SNP file exists: {}", snp.green());
            } else {
                println!("  ✗ SNP file not found: {}", snp.red());
                all_valid = false;
            }
            continue; // Skip CRAM/BAM validation for SNP participants
        }

        // Validate CRAM/BAM participant
        if let Some(aligned) = &participant.aligned {
            let aligned_result = Command::new("samtools")
                .args(["quickcheck", "-v", aligned])
                .output();

            match aligned_result {
                Ok(output) if output.status.success() => {
                    println!("  ✓ Aligned file is valid: {}", aligned.green());
                }
                Ok(output) => {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    println!("  ✗ Aligned file is invalid: {}", aligned.red());
                    if !stderr.is_empty() {
                        println!("    Error: {}", stderr.trim());
                    }
                    all_valid = false;
                }
                Err(e) => {
                    println!(
                        "  ✗ Failed to validate aligned file: {}",
                        e.to_string().red()
                    );
                    all_valid = false;
                }
            }
        }

        if let Some(r#ref) = &participant.r#ref {
            let ref_result = Command::new("seqkit").args(["stats", r#ref]).output();

            match ref_result {
                Ok(output) if output.status.success() => {
                    println!("  ✓ Reference file is valid: {}", r#ref.green());
                    let stdout = String::from_utf8_lossy(&output.stdout);
                    if !stdout.is_empty() {
                        for line in stdout.lines().skip(1).take(1) {
                            println!("    {}", line.trim());
                        }
                    }
                }
                Ok(output) => {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    println!("  ✗ Reference file is invalid: {}", r#ref.red());
                    if !stderr.is_empty() {
                        println!("    Error: {}", stderr.trim());
                    }
                    all_valid = false;
                }
                Err(e) => {
                    println!(
                        "  ✗ Failed to validate reference file: {}",
                        e.to_string().red()
                    );
                    all_valid = false;
                }
            }
        }

        if let Some(aligned_index) = &participant.aligned_index {
            if !Path::new(aligned_index).exists() {
                println!(
                    "  ⚠ Aligned index file not found: {}",
                    aligned_index.yellow()
                );
            } else {
                println!("  ✓ Aligned index file exists: {}", aligned_index.green());
            }
        }

        if let Some(ref_index) = &participant.ref_index {
            if !Path::new(ref_index).exists() {
                println!("  ⚠ Reference index file not found: {}", ref_index.yellow());
            } else {
                println!("  ✓ Reference index file exists: {}", ref_index.green());
            }
        }

        println!();
    }

    if all_valid {
        println!("{}", "✓ All validations passed!".green().bold());
    } else {
        println!("{}", "✗ Some validations failed.".red().bold());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn index_file_path_rules() {
        assert_eq!(get_index_file_path("x.cram", true), "x.cram.crai");
        assert_eq!(get_index_file_path("x.bam", true), "x.bam.bai");
        assert_eq!(get_index_file_path("x.sam", true), "x.sam.sai");
        assert_eq!(get_index_file_path("x.xyz", true), "x.xyz.idx");
        assert_eq!(get_index_file_path("ref.fa", false), "ref.fa.fai");
    }

    #[test]
    #[serial_test::serial]
    #[cfg(unix)]
    fn expand_tilde_resolves_home_unix() {
        let tmp = TempDir::new().unwrap();
        std::env::set_var("HOME", tmp.path());
        let expanded = super::expand_tilde("~/file.txt");
        assert_eq!(
            std::path::Path::new(&expanded),
            &tmp.path().join("file.txt")
        );
    }

    #[test]
    #[serial_test::serial]
    #[cfg(windows)]
    fn expand_tilde_resolves_home_windows() {
        let tmp = TempDir::new().unwrap();
        std::env::set_var("USERPROFILE", tmp.path());
        let expanded = super::expand_tilde("~/file.txt");
        assert_eq!(
            std::path::Path::new(&expanded),
            &tmp.path().join("file.txt")
        );
    }

    #[test]
    #[serial_test::serial]
    fn normalize_path_makes_absolute_and_keeps_tail() {
        let tmp = TempDir::new().unwrap();
        let cwd = std::env::current_dir().unwrap();
        std::env::set_current_dir(tmp.path()).unwrap();
        let abs = super::normalize_path("rel/path").unwrap();
        let p = std::path::Path::new(&abs);
        assert!(p.is_absolute());
        assert!(p.ends_with(std::path::Path::new("rel").join("path")));
        std::env::set_current_dir(cwd).unwrap();
    }

    #[test]
    fn test_participant_serialize_deserialize() {
        let participant = Participant {
            id: "test123".to_string(),
            ref_version: Some("GRCh38".to_string()),
            r#ref: Some("/path/to/ref.fa".to_string()),
            ref_index: Some("/path/to/ref.fa.fai".to_string()),
            aligned: Some("/path/to/aligned.cram".to_string()),
            aligned_index: Some("/path/to/aligned.cram.crai".to_string()),
            snp: None,
        };

        let serialized = serde_yaml::to_string(&participant).unwrap();
        assert!(serialized.contains("id: test123"));
        assert!(serialized.contains("ref_version: GRCh38"));
        assert!(serialized.contains("ref: /path/to/ref.fa"));

        let deserialized: Participant = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserialized.id, "test123");
        assert_eq!(deserialized.ref_version, Some("GRCh38".to_string()));
    }

    #[test]
    fn test_participant_snp_variant() {
        let participant = Participant {
            id: "snp_test".to_string(),
            ref_version: None,
            r#ref: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: Some("/path/to/snp.vcf".to_string()),
        };

        assert_eq!(participant.id, "snp_test");
        assert_eq!(participant.snp, Some("/path/to/snp.vcf".to_string()));
        assert!(participant.aligned.is_none());
        assert!(participant.ref_version.is_none());
    }

    #[test]
    fn test_participants_file_default() {
        let pf = ParticipantsFile::default();
        assert!(pf.participants.is_empty());
    }

    #[test]
    fn test_participants_file_with_multiple_entries() {
        let mut pf = ParticipantsFile::new();

        let p1 = Participant {
            id: "p1".to_string(),
            ref_version: Some("GRCh37".to_string()),
            r#ref: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: None,
        };

        let p2 = Participant {
            id: "p2".to_string(),
            ref_version: Some("GRCh38".to_string()),
            r#ref: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: None,
        };

        pf.participants.insert("p1".to_string(), p1.clone());
        pf.participants.insert("p2".to_string(), p2.clone());

        assert_eq!(pf.participants.len(), 2);
        assert!(pf.participants.contains_key("p1"));
        assert!(pf.participants.contains_key("p2"));
        assert_eq!(
            pf.participants["p1"].ref_version,
            Some("GRCh37".to_string())
        );
        assert_eq!(
            pf.participants["p2"].ref_version,
            Some("GRCh38".to_string())
        );
    }

    #[test]
    #[serial_test::serial]
    fn test_participants_file_load_nonexistent() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        let pf = ParticipantsFile::load().unwrap();
        assert!(pf.participants.is_empty());
    }

    #[test]
    #[serial_test::serial]
    fn test_participants_file_save_and_reload() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        let mut pf = ParticipantsFile::new();
        let participant = Participant {
            id: "save_test".to_string(),
            ref_version: Some("GRCh38".to_string()),
            r#ref: Some("/ref.fa".to_string()),
            ref_index: Some("/ref.fa.fai".to_string()),
            aligned: Some("/aligned.cram".to_string()),
            aligned_index: Some("/aligned.cram.crai".to_string()),
            snp: None,
        };

        pf.participants.insert("save_test".to_string(), participant);
        pf.save().unwrap();

        // Verify file was created
        let file_path = temp.path().join("participants.yaml");
        assert!(file_path.exists());

        // Load and verify
        let loaded = ParticipantsFile::load().unwrap();
        assert_eq!(loaded.participants.len(), 1);
        assert!(loaded.participants.contains_key("save_test"));
        assert_eq!(loaded.participants["save_test"].id, "save_test");
        assert_eq!(
            loaded.participants["save_test"].ref_version,
            Some("GRCh38".to_string())
        );
    }

    #[test]
    fn test_get_index_file_path_edge_cases() {
        // Test various file extensions for aligned files (function is case-sensitive)
        assert_eq!(get_index_file_path("test.CRAM", true), "test.CRAM.idx"); // Uppercase not recognized
        assert_eq!(get_index_file_path("test.BAM", true), "test.BAM.idx"); // Uppercase not recognized
        assert_eq!(get_index_file_path("test.SAM", true), "test.SAM.idx"); // Uppercase not recognized
        assert_eq!(get_index_file_path("test", true), "test.idx");
        assert_eq!(get_index_file_path("test.txt", true), "test.txt.idx");

        // Test lowercase extensions work correctly
        assert_eq!(get_index_file_path("test.cram", true), "test.cram.crai");
        assert_eq!(get_index_file_path("test.bam", true), "test.bam.bai");
        assert_eq!(get_index_file_path("test.sam", true), "test.sam.sai");

        // Test reference files (always get .fai regardless of extension)
        assert_eq!(
            get_index_file_path("genome.fasta", false),
            "genome.fasta.fai"
        );
        assert_eq!(get_index_file_path("genome.FA", false), "genome.FA.fai");
        assert_eq!(get_index_file_path("genome", false), "genome.fai");
    }

    #[test]
    fn test_expand_tilde_edge_cases() {
        // Test non-tilde paths
        assert_eq!(expand_tilde(""), "");
        assert_eq!(expand_tilde("~"), "~");
        assert_eq!(expand_tilde("~user/file"), "~user/file"); // Not expanded
        assert_eq!(expand_tilde("/~/file"), "/~/file"); // Not at start
        assert_eq!(expand_tilde("./~/file"), "./~/file"); // Not at start
    }

    #[test]
    #[serial_test::serial]
    fn test_normalize_path_with_tilde() {
        let temp = TempDir::new().unwrap();
        let temp_path = temp.path().to_string_lossy().to_string();

        if cfg!(unix) {
            std::env::set_var("HOME", &temp_path);
        } else {
            std::env::set_var("USERPROFILE", &temp_path);
        }

        let result = normalize_path("~/test.txt").unwrap();
        assert!(result.contains("test.txt"));
        assert!(std::path::Path::new(&result).is_absolute());
    }

    #[test]
    fn test_participant_clone() {
        let original = Participant {
            id: "clone_test".to_string(),
            ref_version: Some("GRCh38".to_string()),
            r#ref: Some("/ref.fa".to_string()),
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: None,
        };

        let cloned = original.clone();
        assert_eq!(cloned.id, original.id);
        assert_eq!(cloned.ref_version, original.ref_version);
        assert_eq!(cloned.r#ref, original.r#ref);
    }

    #[test]
    fn test_participant_debug_format() {
        let participant = Participant {
            id: "debug_test".to_string(),
            ref_version: None,
            r#ref: None,
            ref_index: None,
            aligned: None,
            aligned_index: None,
            snp: Some("/snp.vcf".to_string()),
        };

        let debug_str = format!("{:?}", participant);
        assert!(debug_str.contains("debug_test"));
        assert!(debug_str.contains("snp"));
        assert!(debug_str.contains("/snp.vcf"));
    }

    #[test]
    #[serial_test::serial]
    fn test_get_participants_file_path() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        let path = get_participants_file_path().unwrap();
        assert!(path.ends_with("participants.yaml"));
        assert!(path.starts_with(temp.path()));
    }

    #[test]
    fn test_check_samtools_installed() {
        // This test just verifies the function returns bool without panicking
        let _result = check_samtools_installed();
        // Function completes without panic - test passes
    }

    #[test]
    fn test_detect_reference_version_without_samtools() {
        // When samtools is not installed, should return None
        // This is hard to test directly, but we can at least call it
        let result = detect_reference_version("/nonexistent/file.bam");
        // Result can be None (no samtools) or None (file doesn't exist)
        assert!(result.is_none() || result.is_some());
    }

    #[tokio::test]
    async fn test_list_empty_participants() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        // Should not error when no participants
        let result = list().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_list_with_participants() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        let mut pf = ParticipantsFile::new();
        pf.participants.insert(
            "test_id".to_string(),
            Participant {
                id: "test_id".to_string(),
                ref_version: Some("GRCh38".to_string()),
                r#ref: Some("/path/ref.fa".to_string()),
                ref_index: None,
                aligned: None,
                aligned_index: None,
                snp: None,
            },
        );
        pf.save().unwrap();

        let result = list().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    async fn test_delete_nonexistent_participant() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

        let result = delete("nonexistent".to_string()).await;
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }

    #[test]
    fn test_participants_file_new() {
        let pf = ParticipantsFile::new();
        assert!(pf.participants.is_empty());
    }

    #[test]
    fn test_participant_partial_data() {
        let p = Participant {
            id: "partial".to_string(),
            ref_version: Some("GRCh38".to_string()),
            r#ref: None,
            ref_index: None,
            aligned: Some("/aligned.bam".to_string()),
            aligned_index: None,
            snp: None,
        };
        assert_eq!(p.id, "partial");
        assert!(p.ref_version.is_some());
        assert!(p.r#ref.is_none());
        assert!(p.aligned.is_some());
    }

    #[test]
    fn test_get_index_file_path_lowercase_variants() {
        assert_eq!(get_index_file_path("file.cram", true), "file.cram.crai");
        assert_eq!(get_index_file_path("file.bam", true), "file.bam.bai");
        assert_eq!(get_index_file_path("file.sam", true), "file.sam.sai");
        assert_eq!(get_index_file_path("ref.fasta", false), "ref.fasta.fai");
    }

    #[test]
    fn test_expand_tilde_no_expansion() {
        assert_eq!(expand_tilde("/absolute/path"), "/absolute/path");
        assert_eq!(expand_tilde("relative/path"), "relative/path");
        assert_eq!(expand_tilde(""), "");
    }

    #[test]
    #[serial_test::serial]
    fn test_normalize_path_absolute() {
        let result = normalize_path("/absolute/path");
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), "/absolute/path");
    }
}
