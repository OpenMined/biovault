use crate::cli::syft_url::SyftURL;
use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{ProjectYaml, SyftPermissions};
use anyhow::Context;
use chrono::Local;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

pub async fn submit(project_path: String, destination: String) -> Result<()> {
    let config = Config::load()?;

    // Parse destination - could be email or full syft URL
    let (datasite_email, participant_url) = if destination.starts_with("syft://") {
        // Full syft URL provided
        let dest_url = SyftURL::parse(&destination)?;
        let email = dest_url.email.clone();
        let participant = if dest_url.fragment.is_some() {
            Some(destination.clone())
        } else {
            None
        };
        (email, participant)
    } else if destination.contains('@') {
        // Just an email/datasite provided
        (destination.clone(), None)
    } else {
        return Err(Error::from(anyhow::anyhow!(
            "Invalid destination: must be either an email address or a syft:// URL"
        )));
    };

    // Determine project directory
    let project_dir = PathBuf::from(&project_path);
    let project_dir = if project_dir.is_relative() && project_path == "." {
        std::env::current_dir()?
    } else {
        project_dir
    };

    if !project_dir.exists() {
        return Err(Error::from(anyhow::anyhow!(
            "Project directory not found: {}",
            project_dir.display()
        )));
    }

    let project_yaml_path = project_dir.join("project.yaml");
    if !project_yaml_path.exists() {
        return Err(Error::from(anyhow::anyhow!(
            "project.yaml not found in: {}",
            project_dir.display()
        )));
    }

    // Load and validate project
    let mut project = ProjectYaml::from_file(&project_yaml_path)?;

    // Override security-sensitive fields
    project.author = config.email.clone();
    project.datasites = Some(vec![datasite_email.clone()]);

    // Handle participants from destination URL
    project.participants = participant_url.clone().map(|url| vec![url]);

    // Hash workflow.nf
    let workflow_path = project_dir.join("workflow.nf");
    if !workflow_path.exists() {
        return Err(Error::from(anyhow::anyhow!(
            "workflow.nf not found in project directory"
        )));
    }
    let workflow_hash = hash_file(&workflow_path)?;

    // Collect and hash asset files
    let assets_dir = project_dir.join("assets");
    let mut asset_files = Vec::new();
    let mut b3_hashes = HashMap::new();

    // Add workflow hash
    b3_hashes.insert("workflow.nf".to_string(), workflow_hash.clone());

    if assets_dir.exists() && assets_dir.is_dir() {
        for entry in WalkDir::new(&assets_dir)
            .min_depth(1)
            .follow_links(false)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
        {
            let path = entry.path();
            let relative_path = path
                .strip_prefix(&project_dir)
                .unwrap_or(path)
                .to_string_lossy()
                .to_string();

            asset_files.push(relative_path.clone());
            let file_hash = hash_file(path)?;
            b3_hashes.insert(relative_path, file_hash);
        }
    }

    project.assets = if asset_files.is_empty() {
        None
    } else {
        Some(asset_files)
    };
    project.b3_hashes = Some(b3_hashes);

    // Calculate project hash from workflow.nf and assets only (deterministic)
    let mut hash_content = String::new();
    hash_content.push_str(&project.name);
    hash_content.push_str(&workflow_hash);
    if let Some(ref hashes) = project.b3_hashes {
        let mut sorted_hashes: Vec<_> = hashes.iter().collect();
        sorted_hashes.sort_by_key(|k| k.0);
        for (file, hash) in sorted_hashes {
            hash_content.push_str(file);
            hash_content.push_str(hash);
        }
    }
    let project_hash = blake3::hash(hash_content.as_bytes()).to_hex().to_string();

    // Create submission folder name
    let date_str = Local::now().format("%Y-%m-%d").to_string();
    let short_hash = &project_hash[0..8];
    let submission_folder_name = format!("{}-{}-{}", project.name, date_str, short_hash);

    // Create submission path
    let submissions_path = config.get_shared_submissions_path()?;
    fs::create_dir_all(&submissions_path)?;

    let submission_path = submissions_path.join(&submission_folder_name);

    // Check if already submitted
    if submission_path.exists() {
        // Verify the hash matches using the same deterministic method
        let existing_project_yaml = submission_path.join("project.yaml");
        if existing_project_yaml.exists() {
            let existing_project = ProjectYaml::from_file(&existing_project_yaml)?;

            // Calculate existing project hash using same method
            let mut existing_hash_content = String::new();
            existing_hash_content.push_str(&existing_project.name);

            // Get workflow hash from existing submission
            if let Some(ref existing_hashes) = existing_project.b3_hashes {
                if let Some(existing_workflow_hash) = existing_hashes.get("workflow.nf") {
                    existing_hash_content.push_str(existing_workflow_hash);

                    let mut sorted_hashes: Vec<_> = existing_hashes.iter().collect();
                    sorted_hashes.sort_by_key(|k| k.0);
                    for (file, hash) in sorted_hashes {
                        existing_hash_content.push_str(file);
                        existing_hash_content.push_str(hash);
                    }
                }
            }

            let existing_hash = blake3::hash(existing_hash_content.as_bytes())
                .to_hex()
                .to_string();

            if existing_hash == project_hash {
                println!("⚠️  This exact project version has already been submitted.");
                println!("   Location: {}", submission_path.display());
                println!("   Hash: {}", short_hash);
                return Ok(());
            }
        }

        return Err(Error::from(anyhow::anyhow!(
            "A submission with this name and date already exists but with different content: {}",
            submission_path.display()
        )));
    }

    // Create submission directory
    fs::create_dir_all(&submission_path)?;

    // Copy project files
    copy_project_files(&project_dir, &submission_path)?;

    // Save updated project.yaml
    project.save(&submission_path.join("project.yaml"))?;

    // Create permissions file
    let permissions = SyftPermissions::new_for_datasite(&datasite_email);
    let permissions_path = submission_path.join("syft.pub.yaml");
    permissions.save(&permissions_path)?;

    println!("✓ Project submitted successfully!");
    println!("  Name: {}", project.name);
    println!("  To: {}", datasite_email);
    if let Some(participants) = &project.participants {
        println!("  Participants: {}", participants.join(", "));
    }
    println!("  Location: {}", submission_path.display());
    println!("  Hash: {}", short_hash);

    Ok(())
}

fn hash_file(path: &Path) -> Result<String> {
    let content =
        fs::read(path).with_context(|| format!("Failed to read file: {}", path.display()))?;
    Ok(blake3::hash(&content).to_hex().to_string())
}

fn copy_project_files(src: &Path, dest: &Path) -> Result<()> {
    // Copy workflow.nf
    let src_workflow = src.join("workflow.nf");
    let dest_workflow = dest.join("workflow.nf");
    fs::copy(&src_workflow, &dest_workflow)
        .with_context(|| "Failed to copy workflow.nf".to_string())?;

    // Copy assets directory if it exists
    let src_assets = src.join("assets");
    if src_assets.exists() && src_assets.is_dir() {
        let dest_assets = dest.join("assets");
        fs::create_dir_all(&dest_assets)?;

        for entry in WalkDir::new(&src_assets)
            .min_depth(1)
            .follow_links(false)
            .into_iter()
            .filter_map(|e| e.ok())
        {
            let src_path = entry.path();
            let relative_path = src_path.strip_prefix(&src_assets).unwrap();
            let dest_path = dest_assets.join(relative_path);

            if entry.file_type().is_dir() {
                fs::create_dir_all(&dest_path)?;
            } else {
                if let Some(parent) = dest_path.parent() {
                    fs::create_dir_all(parent)?;
                }
                fs::copy(src_path, &dest_path)
                    .with_context(|| format!("Failed to copy file: {}", src_path.display()))?;
            }
        }
    }

    Ok(())
}
