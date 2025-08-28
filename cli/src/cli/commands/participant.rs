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
    pub ref_version: String,
    pub r#ref: String,
    pub ref_index: String,
    pub aligned: String,
    pub aligned_index: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ParticipantsFile {
    pub participants: HashMap<String, Participant>,
}

impl ParticipantsFile {
    fn new() -> Self {
        Self {
            participants: HashMap::new(),
        }
    }

    fn load() -> Result<Self> {
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
    let home = dirs::home_dir().ok_or_else(|| anyhow!("Could not determine home directory"))?;
    Ok(home.join(".biovault").join("participants.yaml"))
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
        .arg(format!("samtools view -H {} | grep -c 'SN:chrEBV'", aligned_file))
        .output();

    let check_ki270 = Command::new("sh")
        .arg("-c")
        .arg(format!("samtools view -H {} | grep -c 'SN:KI270'", aligned_file))
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

fn normalize_path(path: &str) -> Result<String> {
    let path_buf = PathBuf::from(path);
    if path_buf.is_absolute() {
        Ok(path.to_string())
    } else {
        Ok(path.to_string())
    }
}

pub async fn add(id: Option<String>, aligned: Option<String>) -> Result<()> {
    println!("{}", "Adding new participant...".green().bold());
    
    if !check_samtools_installed() {
        warn!("samtools is not installed. Index creation will not be available.");
        println!("{}", "Warning: samtools is not installed. You won't be able to create index files automatically.".yellow());
    }

    let mut participants_file = ParticipantsFile::load()?;
    let mut index_jobs = Vec::new();

    let participant_id = match id {
        Some(id) => id,
        None => Input::<String>::new()
            .with_prompt("Enter participant ID")
            .interact_text()?,
    };

    if participants_file.participants.contains_key(&participant_id) {
        let overwrite = Confirm::new()
            .with_prompt(format!("Participant '{}' already exists. Overwrite?", participant_id))
            .interact()?;
        
        if !overwrite {
            println!("Cancelled.");
            return Ok(());
        }
    }

    let aligned_file = match aligned {
        Some(file) => file,
        None => Input::<String>::new()
            .with_prompt("Enter aligned file path (.cram, .bam, or .sam)")
            .validate_with(|input: &String| -> Result<(), &str> {
                if Path::new(input).exists() {
                    Ok(())
                } else {
                    Err("File does not exist")
                }
            })
            .interact_text()?,
    };

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
        Confirm::new()
            .with_prompt(format!("Use this index file: {}?", aligned_index))
            .default(true)
            .interact()?
    } else {
        false
    };

    let final_aligned_index = if use_aligned_index {
        aligned_index.clone()
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
    
    let ref_version = ref_version_choices[ref_version_idx].to_string();

    let ref_file = Input::<String>::new()
        .with_prompt("Enter reference genome file path (.fa or .fasta)")
        .validate_with(|input: &String| -> Result<(), &str> {
            if Path::new(input).exists() {
                Ok(())
            } else {
                Err("File does not exist")
            }
        })
        .interact_text()?;

    let ref_index = get_index_file_path(&ref_file, false);
    let ref_index_exists = Path::new(&ref_index).exists();
    
    if ref_index_exists {
        println!("✓ Reference index found: {}", ref_index.green());
    } else {
        println!("⚠ Reference index not found: {}", ref_index.yellow());
    }

    let use_ref_index = if ref_index_exists {
        Confirm::new()
            .with_prompt(format!("Use this reference index file: {}?", ref_index))
            .default(true)
            .interact()?
    } else {
        false
    };

    let final_ref_index = if use_ref_index {
        ref_index.clone()
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
        ref_version,
        r#ref: normalize_path(&ref_file)?,
        ref_index: normalize_path(&final_ref_index)?,
        aligned: normalize_path(&aligned_file)?,
        aligned_index: normalize_path(&final_aligned_index)?,
    };

    participants_file.participants.insert(participant_id.clone(), participant);
    participants_file.save()?;

    println!("{}", format!("✓ Participant '{}' added successfully!", participant_id).green().bold());

    if !index_jobs.is_empty() && check_samtools_installed() {
        println!("\n{}", "Creating index files...".cyan().bold());
        for (file, is_aligned) in index_jobs {
            println!("Indexing: {}", file);
            let result = if is_aligned {
                Command::new("samtools")
                    .args(&["index", &file])
                    .output()
            } else {
                Command::new("samtools")
                    .args(&["faidx", &file])
                    .output()
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
        println!("    ref_version: {}", participant.ref_version);
        println!("    ref: {}", participant.r#ref);
        println!("    ref_index: {}", participant.ref_index);
        println!("    aligned: {}", participant.aligned);
        println!("    aligned_index: {}", participant.aligned_index);
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
    
    println!("{}", format!("✓ Participant '{}' deleted successfully!", id).green().bold());
    Ok(())
}

pub async fn validate(id: Option<String>) -> Result<()> {
    if !check_samtools_installed() {
        return Err(anyhow!("samtools is not installed. Please install samtools first."));
    }

    let check_seqkit = Command::new("which")
        .arg("seqkit")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false);

    if !check_seqkit {
        return Err(anyhow!("seqkit is not installed. Please install seqkit first."));
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
        
        let aligned_result = Command::new("samtools")
            .args(&["quickcheck", "-v", &participant.aligned])
            .output();
        
        match aligned_result {
            Ok(output) if output.status.success() => {
                println!("  ✓ Aligned file is valid: {}", participant.aligned.green());
            }
            Ok(output) => {
                let stderr = String::from_utf8_lossy(&output.stderr);
                println!("  ✗ Aligned file is invalid: {}", participant.aligned.red());
                if !stderr.is_empty() {
                    println!("    Error: {}", stderr.trim());
                }
                all_valid = false;
            }
            Err(e) => {
                println!("  ✗ Failed to validate aligned file: {}", e.to_string().red());
                all_valid = false;
            }
        }

        let ref_result = Command::new("seqkit")
            .args(&["stats", &participant.r#ref])
            .output();
        
        match ref_result {
            Ok(output) if output.status.success() => {
                println!("  ✓ Reference file is valid: {}", participant.r#ref.green());
                let stdout = String::from_utf8_lossy(&output.stdout);
                if !stdout.is_empty() {
                    for line in stdout.lines().skip(1).take(1) {
                        println!("    {}", line.trim());
                    }
                }
            }
            Ok(output) => {
                let stderr = String::from_utf8_lossy(&output.stderr);
                println!("  ✗ Reference file is invalid: {}", participant.r#ref.red());
                if !stderr.is_empty() {
                    println!("    Error: {}", stderr.trim());
                }
                all_valid = false;
            }
            Err(e) => {
                println!("  ✗ Failed to validate reference file: {}", e.to_string().red());
                all_valid = false;
            }
        }

        if !Path::new(&participant.aligned_index).exists() {
            println!("  ⚠ Aligned index file not found: {}", participant.aligned_index.yellow());
        } else {
            println!("  ✓ Aligned index file exists: {}", participant.aligned_index.green());
        }

        if !Path::new(&participant.ref_index).exists() {
            println!("  ⚠ Reference index file not found: {}", participant.ref_index.yellow());
        } else {
            println!("  ✓ Reference index file exists: {}", participant.ref_index.green());
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