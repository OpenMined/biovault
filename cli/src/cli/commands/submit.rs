use crate::cli::commands::messages::init_message_system;
use crate::cli::syft_url::SyftURL;
use crate::config::Config;
use crate::error::{Error, Result};
use crate::messages::{Message, MessageType};
use crate::syftbox::storage::{SyftBoxStorage, WritePolicy};
#[cfg(test)]
use crate::syftbox::SyftBoxApp;
use crate::types::{ProjectYaml, SyftPermissions};
use anyhow::Context;
use chrono::Local;
use dialoguer::{Confirm, Editor};
use serde::de::DeserializeOwned;
use serde::Serialize;
use serde_json::json;
use serde_yaml;
use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

pub async fn submit(
    project_path: String,
    destination: String,
    non_interactive: bool,
    force: bool,
) -> Result<()> {
    let config = Config::load()?;

    let storage = syftbox_storage(&config)?;

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

    // Check if this is a self-submission (for testing)
    let is_self_submission = datasite_email == config.email;
    if is_self_submission {
        println!("🧪 Self-submission detected (testing mode)");
    }

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
    let mut datasites = project.datasites.clone().unwrap_or_default();
    if datasites.is_empty() {
        datasites.push(datasite_email.clone());
    } else if !datasites.contains(&datasite_email) {
        datasites.push(datasite_email.clone());
    }
    project.datasites = Some(datasites);

    // Handle participants from destination URL
    if let Some(ref url) = participant_url {
        project.participants = Some(vec![url.clone()]);
    }

    // Hash workflow file
    let workflow_path = project_dir.join(&project.workflow);
    if !workflow_path.exists() {
        return Err(Error::from(anyhow::anyhow!(
            "Workflow file '{}' not found in project directory",
            project.workflow
        )));
    }
    let workflow_hash = hash_file(&workflow_path)?;

    // Collect and hash asset files
    let assets_dir = project_dir.join("assets");
    let mut asset_files = Vec::new();
    let mut b3_hashes = HashMap::new();

    // Add workflow hash
    b3_hashes.insert(project.workflow.clone(), workflow_hash.clone());

    if assets_dir.exists() && assets_dir.is_dir() {
        for entry in WalkDir::new(&assets_dir)
            .min_depth(1)
            .follow_links(false)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
        {
            let path = entry.path();
            // Relative path from project root (for hashes)
            let relative_path = path
                .strip_prefix(&project_dir)
                .unwrap_or(path)
                .to_string_lossy()
                .to_string();
            // Relative path from assets dir (for YAML assets list)
            let assets_rel = path
                .strip_prefix(&assets_dir)
                .unwrap_or(path)
                .to_string_lossy()
                .to_string();

            asset_files.push(assets_rel);
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

    // Calculate project hash from workflow and assets only (deterministic)
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
    storage.ensure_dir(&submissions_path)?;

    let submission_path = submissions_path.join(&submission_folder_name);

    // Check if already submitted
    if submission_path.exists() {
        // Verify the hash matches using the same deterministic method
        let existing_project_yaml = submission_path.join("project.yaml");
        if existing_project_yaml.exists() {
            let existing_project =
                read_yaml_from_storage::<ProjectYaml>(&storage, &existing_project_yaml)?;

            // Calculate existing project hash using same method
            let mut existing_hash_content = String::new();
            existing_hash_content.push_str(&existing_project.name);

            // Get workflow hash from existing submission
            if let Some(ref existing_hashes) = existing_project.b3_hashes {
                let workflow_key = if existing_project.workflow.trim().is_empty() {
                    "workflow.nf"
                } else {
                    existing_project.workflow.as_str()
                };
                if let Some(existing_workflow_hash) = existing_hashes
                    .get(workflow_key)
                    .or_else(|| existing_hashes.get("workflow.nf"))
                {
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
                // Check if we're submitting to a different recipient
                let existing_permissions_path = submission_path.join("syft.pub.yaml");
                if existing_permissions_path.exists() && datasite_email != config.email {
                    println!(
                        "📝 Updating submission with new recipient: {}",
                        datasite_email
                    );

                    let mut existing_project =
                        read_yaml_from_storage::<ProjectYaml>(&storage, &existing_project_yaml)?;
                    let mut share_datasites =
                        existing_project.datasites.clone().unwrap_or_default();
                    if share_datasites.is_empty() {
                        share_datasites = project
                            .datasites
                            .clone()
                            .unwrap_or_else(|| vec![datasite_email.clone()]);
                    }
                    if !share_datasites.contains(&datasite_email) {
                        share_datasites.push(datasite_email.clone());
                    }
                    share_datasites = dedupe_recipients(share_datasites);

                    let mut permissions = SyftPermissions::new_for_datasites(&share_datasites);
                    permissions.add_rule(
                        "results/**/*",
                        share_datasites.clone(),
                        share_datasites.clone(),
                    );
                    write_yaml_to_storage(
                        &storage,
                        &existing_permissions_path,
                        &permissions,
                        None,
                        true,
                    )?;

                    let recipients = recipients_from_permissions(
                        &storage,
                        &existing_permissions_path,
                        &config.email,
                    )?;

                    copy_project_files(
                        &project_dir,
                        &submission_path,
                        &storage,
                        &recipients,
                        &project.workflow,
                    )?;

                    existing_project.datasites = Some(share_datasites);

                    // Handle participants from destination URL
                    if let Some(ref new_participant) = participant_url {
                        if let Some(ref mut participants) = existing_project.participants {
                            if !participants.contains(new_participant) {
                                participants.push(new_participant.clone());
                            }
                        } else {
                            existing_project.participants = Some(vec![new_participant.clone()]);
                        }
                    }

                    write_yaml_to_storage(
                        &storage,
                        &existing_project_yaml,
                        &existing_project,
                        Some(recipients),
                        true,
                    )?;

                    println!("✓ Permissions updated for existing submission");
                    println!("  Location: {}", submission_path.display());

                    // Send the message to the new recipient
                    return send_project_message(
                        &config,
                        &project,
                        &datasite_email,
                        &submission_path,
                        &submission_folder_name,
                        &project_hash,
                        &date_str,
                        non_interactive,
                        is_self_submission,
                    );
                } else if !force {
                    println!("⚠️  This exact project version has already been submitted.");
                    println!("   Location: {}", submission_path.display());
                    println!("   Hash: {}", short_hash);
                    println!("   Use --force to resend the message.");
                    return Ok(());
                } else {
                    println!(
                        "ℹ️  Project already submitted, but --force flag used. Sending message..."
                    );
                    println!("   Location: {}", submission_path.display());
                    println!("   Hash: {}", short_hash);

                    // Check if destination has changed
                    let existing_datasites = existing_project.datasites.clone().unwrap_or_default();
                    let destination_changed = !existing_datasites.contains(&datasite_email);

                    if destination_changed {
                        // Update project.yaml and syft.pub.yaml with new destination
                        println!("✓ Updating destination to: {}", datasite_email);

                        // Update project.yaml with new datasite
                        let mut updated_project = existing_project.clone();
                        let mut updated_datasites =
                            updated_project.datasites.clone().unwrap_or_default();
                        if updated_datasites.is_empty() {
                            updated_datasites = project.datasites.clone().unwrap_or_default();
                        }
                        if updated_datasites.is_empty() {
                            updated_datasites.push(datasite_email.clone());
                        } else if !updated_datasites.contains(&datasite_email) {
                            updated_datasites.push(datasite_email.clone());
                        }
                        updated_datasites = dedupe_recipients(updated_datasites);
                        updated_project.datasites = Some(updated_datasites.clone());
                        if let Some(ref url) = participant_url {
                            updated_project.participants = Some(vec![url.clone()]);
                        }
                        let mut recipients = updated_datasites.clone();
                        recipients.push(config.email.clone());
                        let recipients = dedupe_recipients(recipients);
                        write_yaml_to_storage(
                            &storage,
                            &existing_project_yaml,
                            &updated_project,
                            Some(recipients.clone()),
                            true,
                        )?;

                        // Update syft.pub.yaml
                        let existing_permissions_path = submission_path.join("syft.pub.yaml");
                        let mut updated_permissions =
                            SyftPermissions::new_for_datasites(&updated_datasites);
                        updated_permissions.add_rule(
                            "results/**/*",
                            updated_datasites.clone(),
                            updated_datasites.clone(),
                        );
                        write_yaml_to_storage(
                            &storage,
                            &existing_permissions_path,
                            &updated_permissions,
                            None,
                            true,
                        )?;

                        copy_project_files(
                            &project_dir,
                            &submission_path,
                            &storage,
                            &recipients,
                            &project.workflow,
                        )?;
                    }

                    // Send new request message for same project
                    return send_project_message(
                        &config,
                        &project,
                        &datasite_email,
                        &submission_path,
                        &submission_folder_name,
                        &project_hash,
                        &date_str,
                        non_interactive,
                        is_self_submission,
                    );
                }
            } else {
                // Project has changed
                if !force {
                    println!("⚠️  A submission exists with the same name but different content.");
                    println!("   Existing: {}", submission_path.display());
                    println!("   Use --force to create a new version.");
                    return Ok(());
                }

                // Create new submission with updated timestamp
                let new_date_str = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
                let new_submission_folder =
                    format!("{}-{}-{}", project.name, new_date_str, short_hash);
                let new_submission_path = submissions_path.join(&new_submission_folder);

                return create_and_submit_project(
                    &config,
                    &storage,
                    &project,
                    &datasite_email,
                    &project_dir,
                    &new_submission_path,
                    &new_submission_folder,
                    &project_hash,
                    &new_date_str,
                    participant_url,
                    non_interactive,
                );
            }
        }
    }

    // New submission - create and submit
    create_and_submit_project(
        &config,
        &storage,
        &project,
        &datasite_email,
        &project_dir,
        &submission_path,
        &submission_folder_name,
        &project_hash,
        &date_str,
        participant_url,
        non_interactive,
    )
}

fn hash_file(path: &Path) -> Result<String> {
    let content =
        fs::read(path).with_context(|| format!("Failed to read file: {}", path.display()))?;
    Ok(blake3::hash(&content).to_hex().to_string())
}

fn syftbox_storage(config: &Config) -> Result<SyftBoxStorage> {
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(SyftBoxStorage::new(&data_dir))
}

fn read_yaml_from_storage<T: DeserializeOwned>(storage: &SyftBoxStorage, path: &Path) -> Result<T> {
    let bytes = storage
        .read_plaintext_file(path)
        .with_context(|| format!("Failed to read {:?}", path))?;
    Ok(serde_yaml::from_slice(&bytes).with_context(|| format!("Failed to parse {:?}", path))?)
}

fn write_yaml_to_storage<T: Serialize>(
    storage: &SyftBoxStorage,
    path: &Path,
    value: &T,
    recipients: Option<Vec<String>>,
    overwrite: bool,
) -> Result<()> {
    let yaml = serde_yaml::to_string(value)?;
    let policy = if let Some(recipients) = recipients {
        WritePolicy::Envelope {
            recipients,
            hint: path
                .file_name()
                .map(|name| name.to_string_lossy().into_owned()),
        }
    } else {
        WritePolicy::Plaintext
    };
    storage.write_with_shadow(path, yaml.as_bytes(), policy, overwrite)?;
    Ok(())
}

fn dedupe_recipients(mut recipients: Vec<String>) -> Vec<String> {
    let mut set = BTreeSet::new();
    for recipient in recipients.drain(..) {
        let trimmed = recipient.trim();
        if !trimmed.is_empty() {
            set.insert(trimmed.to_string());
        }
    }
    set.into_iter().collect()
}

fn recipients_from_permissions(
    storage: &SyftBoxStorage,
    permissions_path: &Path,
    sender: &str,
) -> Result<Vec<String>> {
    let permissions: SyftPermissions = read_yaml_from_storage(storage, permissions_path)?;
    let mut recipients = vec![sender.to_string()];
    for rule in permissions.rules {
        recipients.extend(rule.access.read);
        recipients.extend(rule.access.write);
    }
    Ok(dedupe_recipients(recipients))
}

fn copy_project_files(
    src: &Path,
    dest: &Path,
    storage: &SyftBoxStorage,
    recipients: &[String],
    workflow_name: &str,
) -> Result<()> {
    storage.ensure_dir(dest)?;
    // Copy workflow file
    let src_workflow = src.join(workflow_name);
    let dest_workflow = dest.join(workflow_name);
    let workflow_bytes = fs::read(&src_workflow)
        .with_context(|| format!("Failed to read {} from project", workflow_name))?;
    let workflow_policy = WritePolicy::Envelope {
        recipients: recipients.to_vec(),
        hint: Some(workflow_name.to_string()),
    };
    storage.write_with_shadow(&dest_workflow, &workflow_bytes, workflow_policy, true)?;

    // Copy assets directory if it exists
    let src_assets = src.join("assets");
    if src_assets.exists() && src_assets.is_dir() {
        let dest_assets = dest.join("assets");
        storage.ensure_dir(&dest_assets)?;

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
                storage.ensure_dir(&dest_path)?;
            } else {
                let bytes = fs::read(src_path)
                    .with_context(|| format!("Failed to read file: {}", src_path.display()))?;
                let hint = dest_path
                    .strip_prefix(dest)
                    .ok()
                    .and_then(|p| p.to_str())
                    .map(|s| format!("assets/{s}"));
                let asset_policy = WritePolicy::Envelope {
                    recipients: recipients.to_vec(),
                    hint,
                };
                storage.write_with_shadow(&dest_path, &bytes, asset_policy, true)?;
            }
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn create_and_submit_project(
    config: &Config,
    storage: &SyftBoxStorage,
    project: &ProjectYaml,
    datasite_email: &str,
    project_dir: &Path,
    submission_path: &Path,
    submission_folder_name: &str,
    project_hash: &str,
    date_str: &str,
    participant_url: Option<String>,
    non_interactive: bool,
) -> Result<()> {
    // Create submission directory
    storage.ensure_dir(submission_path)?;

    let share_datasites = dedupe_recipients(
        project
            .datasites
            .clone()
            .unwrap_or_else(|| vec![datasite_email.to_string()]),
    );
    let mut recipients = share_datasites.clone();
    recipients.push(config.email.clone());
    let recipients = dedupe_recipients(recipients);

    // Copy project files
    copy_project_files(
        project_dir,
        submission_path,
        storage,
        &recipients,
        &project.workflow,
    )?;

    // Save updated project.yaml
    let mut final_project = project.clone();
    final_project.participants = participant_url.map(|url| vec![url]);
    let project_yaml_path = submission_path.join("project.yaml");
    write_yaml_to_storage(
        storage,
        &project_yaml_path,
        &final_project,
        Some(recipients),
        true,
    )?;

    // Create permissions file
    let mut permissions = SyftPermissions::new_for_datasites(&share_datasites);
    permissions.add_rule(
        "results/**/*",
        share_datasites.clone(),
        share_datasites.clone(),
    );
    let permissions_path = submission_path.join("syft.pub.yaml");
    write_yaml_to_storage(storage, &permissions_path, &permissions, None, true)?;

    println!("✓ Project submitted successfully!");
    println!("  Name: {}", project.name);
    println!("  To: {}", datasite_email);
    if let Some(participants) = &final_project.participants {
        println!("  Participants: {}", participants.join(", "));
    }
    println!("  Location: {}", submission_path.display());
    let short_hash = &project_hash[0..8];
    println!("  Hash: {}", short_hash);

    // Send the project message
    // Determine if this is a self-submission
    let is_self_submission = datasite_email == config.email;

    send_project_message(
        config,
        &final_project,
        datasite_email,
        submission_path,
        submission_folder_name,
        project_hash,
        date_str,
        non_interactive,
        is_self_submission,
    )
}

#[allow(clippy::too_many_arguments)]
fn send_project_message(
    config: &Config,
    project: &ProjectYaml,
    datasite_email: &str,
    submission_path: &Path,
    submission_folder_name: &str,
    project_hash: &str,
    date_str: &str,
    non_interactive: bool,
    is_self_submission: bool,
) -> Result<()> {
    // Build a syft:// URL for the saved submission so the receiver can locate it
    let datasite_root = config.get_datasite_path()?;
    let rel_from_datasite = submission_path
        .strip_prefix(&datasite_root)
        .unwrap_or(submission_path)
        .to_string_lossy()
        .to_string();
    let submission_syft_url = if is_self_submission {
        // For self-submission, the URL should point to own datasite
        format!("syft://{}/{}", config.email, rel_from_datasite)
    } else {
        format!("syft://{}/{}", config.email, rel_from_datasite)
    };

    // Prepare default message body and allow user to override
    let default_body = "I would like to run the following project.".to_string();
    let mut body = if non_interactive {
        // In non-interactive mode, use default message
        default_body.clone()
    } else {
        let use_custom = Confirm::new()
            .with_prompt("Write a custom message body?")
            .default(false)
            .interact()
            .unwrap_or(false);

        if use_custom {
            match Editor::new().edit(&default_body) {
                Ok(Some(content)) if !content.trim().is_empty() => content,
                _ => default_body.clone(),
            }
        } else {
            default_body.clone()
        }
    };

    // Add handy paths for the recipient to copy/paste
    let sender_local_path = submission_path.to_string_lossy().to_string();
    let receiver_local_path_template = format!(
        "$SYFTBOX_DATA_DIR/datasites/{}/shared/biovault/submissions/{}",
        config.email, submission_folder_name
    );
    body.push_str(&format!(
        "\n\nSubmission location references:\n- syft URL: {}\n- Sender local path: {}\n- Receiver local path (template): {}\n",
        submission_syft_url, sender_local_path, receiver_local_path_template
    ));

    // Construct metadata for the message
    let metadata = json!({
        "project": project,
        "project_location": submission_syft_url,
        // Receiver guidance about which of their participants to use
        "participants": "With your participants: ALL",
        // Date component used in the submission folder name
        "date": date_str,
        // Explicit list of asset files (if any)
        "assets": project.assets.clone().unwrap_or_default(),
        // Helpful paths for receiver tooling
        "sender_local_path": sender_local_path,
        "receiver_local_path_template": receiver_local_path_template,
    });

    // Initialize messaging system and send a project message
    let (db, sync) = init_message_system(config)?;

    let mut msg = Message::new(config.email.clone(), datasite_email.to_string(), body);
    msg.subject = Some(format!("Project Request - {}", project.name));
    msg.message_type = MessageType::Project {
        project_name: project.name.clone(),
        submission_id: submission_folder_name.to_string(),
        files_hash: Some(project_hash.to_string()),
    };
    msg.metadata = Some(metadata);

    db.insert_message(&msg)?;

    if is_self_submission {
        // For self-submission, write directly to our own inbox
        println!("🧪 Writing message directly to own inbox for testing");
        // The sync system will handle this correctly as a self-message
    }

    // Try to send immediately; if offline, it will remain queued locally
    sync.send_message(&msg.id)
        .context("failed to send project message")?;

    println!("✉️  Project message prepared for {}", datasite_email);
    if let Some(subj) = &msg.subject {
        println!("  Subject: {}", subj);
    }
    println!("  Submission URL: {}", submission_syft_url);

    Ok(())
}

#[allow(dead_code)]
fn update_permissions_for_new_recipient(
    storage: &SyftBoxStorage,
    permissions_path: &Path,
    new_recipient: &str,
) -> Result<()> {
    let mut permissions: SyftPermissions = read_yaml_from_storage(storage, permissions_path)?;

    // Add new recipient to all rules
    for rule in &mut permissions.rules {
        if !rule.access.read.contains(&new_recipient.to_string()) {
            rule.access.read.push(new_recipient.to_string());
        }

        // For results pattern, also add write permission
        if rule.pattern == "results/**/*" && !rule.access.write.contains(&new_recipient.to_string())
        {
            rule.access.write.push(new_recipient.to_string());
        }
    }

    // Save updated permissions
    write_yaml_to_storage(storage, permissions_path, &permissions, None, true)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::commands::messages::get_message_db_path;
    use crate::config::{
        clear_test_biovault_home, clear_test_syftbox_data_dir, set_test_biovault_home,
        set_test_syftbox_data_dir,
    };
    use crate::messages::{MessageDb, MessageType};
    use serial_test::serial;
    use std::collections::HashMap;
    use std::time::Duration;
    use tempfile::TempDir;

    fn write_project_yaml(path: &Path, datasites: Option<Vec<String>>) {
        let project = ProjectYaml {
            name: "demo-project".into(),
            author: "author@example.com".into(),
            datasites,
            participants: None,
            workflow: "workflow.nf".into(),
            template: None,
            assets: Some(vec!["assets/data.csv".into()]),
            inputs: None,
            outputs: None,
            b3_hashes: None,
            extra: HashMap::new(),
        };
        project.save(&path.to_path_buf()).unwrap();
    }

    fn setup_submit_env(temp: &TempDir, config_email: &str) -> (PathBuf, SyftBoxApp, PathBuf) {
        let bv_home = temp.path().join(".biovault");
        fs::create_dir_all(&bv_home).unwrap();
        set_test_biovault_home(&bv_home);

        let config = Config {
            email: config_email.into(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };
        config.save(bv_home.join("config.yaml")).unwrap();

        let syft_dir = temp.path().join("syft_data");
        fs::create_dir_all(&syft_dir).unwrap();
        set_test_syftbox_data_dir(&syft_dir);
        let app = SyftBoxApp::new(&syft_dir, config_email, "biovault").unwrap();

        let project_dir = temp.path().join("project");
        fs::create_dir_all(project_dir.join("assets")).unwrap();
        fs::write(project_dir.join("workflow.nf"), b"process MAIN {}").unwrap();
        fs::write(project_dir.join("assets/data.csv"), b"data").unwrap();
        write_project_yaml(&project_dir.join("project.yaml"), None);

        (project_dir, app, bv_home)
    }

    fn list_submissions(app: &SyftBoxApp) -> Vec<PathBuf> {
        let submissions_dir = app
            .data_dir
            .join("datasites")
            .join(&app.email)
            .join("shared")
            .join("biovault")
            .join("submissions");
        app.storage.list_dir(&submissions_dir).unwrap_or_default()
    }

    #[test]
    fn hash_file_computes_blake3() {
        let tmp = TempDir::new().unwrap();
        let f = tmp.path().join("a.txt");
        fs::write(&f, b"content").unwrap();
        let h = hash_file(&f).unwrap();
        assert_eq!(h.len(), 64);
    }

    #[test]
    fn copy_project_files_copies_workflow_and_assets() {
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("src");
        let dest = tmp.path().join("datasites/dest");
        fs::create_dir_all(&dest).unwrap();
        fs::create_dir_all(src.join("assets/nested")).unwrap();
        fs::write(src.join("workflow.nf"), b"wf").unwrap();
        fs::write(src.join("assets/nested/file.bin"), b"x").unwrap();
        let storage = SyftBoxStorage::new(tmp.path());
        let recipients = vec!["peer@example.com".to_string()];

        copy_project_files(&src, &dest, &storage, &recipients, "workflow.nf").unwrap();

        assert!(dest.join("workflow.nf").exists());
        assert!(dest.join("assets/nested/file.bin").exists());

        // No assets dir case should still succeed
        let src2 = tmp.path().join("src2");
        fs::create_dir_all(&src2).unwrap();
        fs::write(src2.join("workflow.nf"), b"wf").unwrap();
        let dest2 = tmp.path().join("datasites/dest2");
        fs::create_dir_all(&dest2).unwrap();
        copy_project_files(&src2, &dest2, &storage, &recipients, "workflow.nf").unwrap();
    }

    #[test]
    fn test_hash_file_with_different_content() {
        let tmp = TempDir::new().unwrap();

        let file1 = tmp.path().join("file1.txt");
        fs::write(&file1, b"content1").unwrap();
        let hash1 = hash_file(&file1).unwrap();

        let file2 = tmp.path().join("file2.txt");
        fs::write(&file2, b"content2").unwrap();
        let hash2 = hash_file(&file2).unwrap();

        // Different content should produce different hashes
        assert_ne!(hash1, hash2);

        // Same content should produce same hash
        let file3 = tmp.path().join("file3.txt");
        fs::write(&file3, b"content1").unwrap();
        let hash3 = hash_file(&file3).unwrap();
        assert_eq!(hash1, hash3);
    }

    #[test]
    fn test_hash_file_empty_file() {
        let tmp = TempDir::new().unwrap();
        let empty_file = tmp.path().join("empty.txt");
        fs::write(&empty_file, b"").unwrap();
        let hash = hash_file(&empty_file).unwrap();
        // Blake3 hash of empty string is a known value
        assert_eq!(hash.len(), 64);
    }

    #[test]
    fn test_copy_project_files_missing_workflow() {
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("src");
        let dest = tmp.path().join("dest");
        fs::create_dir_all(&src).unwrap();
        fs::create_dir_all(&dest).unwrap();
        // No workflow.nf file

        let storage = SyftBoxStorage::new(tmp.path());
        let recipients = vec!["peer@example.com".to_string()];
        let result = copy_project_files(&src, &dest, &storage, &recipients, "workflow.nf");
        assert!(result.is_err());
    }

    #[test]
    fn test_copy_project_files_with_symlinks() {
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("src");
        let dest = tmp.path().join("datasites/dest");

        fs::create_dir_all(&src).unwrap();
        fs::create_dir_all(&dest).unwrap();
        fs::create_dir_all(src.join("assets")).unwrap();

        fs::write(src.join("workflow.nf"), b"workflow").unwrap();
        fs::write(src.join("assets/real_file.txt"), b"real").unwrap();

        // Create a symlink (this test will only work on Unix-like systems)
        #[cfg(unix)]
        {
            use std::os::unix::fs::symlink;
            let _ = symlink(
                src.join("assets/real_file.txt"),
                src.join("assets/link_file.txt"),
            );
        }

        let storage = SyftBoxStorage::new(tmp.path());
        let recipients = vec!["peer@example.com".to_string()];
        let result = copy_project_files(&src, &dest, &storage, &recipients, "workflow.nf");
        assert!(result.is_ok());
        assert!(dest.join("workflow.nf").exists());
    }

    #[test]
    fn test_copy_project_files_deeply_nested() {
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("src");
        let dest = tmp.path().join("datasites/dest");

        fs::create_dir_all(&dest).unwrap();
        fs::create_dir_all(src.join("assets/a/b/c/d")).unwrap();
        fs::write(src.join("workflow.nf"), b"wf").unwrap();
        fs::write(src.join("assets/a/b/c/d/deep.txt"), b"deep").unwrap();

        let storage = SyftBoxStorage::new(tmp.path());
        let recipients = vec!["peer@example.com".to_string()];
        copy_project_files(&src, &dest, &storage, &recipients, "workflow.nf").unwrap();

        assert!(dest.join("workflow.nf").exists());
        assert!(dest.join("assets/a/b/c/d/deep.txt").exists());

        let content = fs::read(dest.join("assets/a/b/c/d/deep.txt")).unwrap();
        assert_eq!(content, b"deep");
    }

    #[test]
    fn test_copy_project_files_preserves_content() {
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("src");
        let dest = tmp.path().join("datasites/dest");

        fs::create_dir_all(&dest).unwrap();
        fs::create_dir_all(src.join("assets")).unwrap();

        let workflow_content = b"nextflow.enable.dsl=2\nworkflow { }";
        let asset_content = b"important data";

        fs::write(src.join("workflow.nf"), workflow_content).unwrap();
        fs::write(src.join("assets/data.csv"), asset_content).unwrap();

        let storage = SyftBoxStorage::new(tmp.path());
        let recipients = vec!["peer@example.com".to_string()];
        copy_project_files(&src, &dest, &storage, &recipients, "workflow.nf").unwrap();

        let copied_workflow = fs::read(dest.join("workflow.nf")).unwrap();
        let copied_asset = fs::read(dest.join("assets/data.csv")).unwrap();

        assert_eq!(copied_workflow, workflow_content);
        assert_eq!(copied_asset, asset_content);
    }

    #[test]
    fn test_hash_file_nonexistent() {
        let result = hash_file(Path::new("/nonexistent/file.txt"));
        assert!(result.is_err());
    }

    #[test]
    fn update_permissions_adds_new_recipient_without_duplicates() {
        let tmp = TempDir::new().unwrap();
        let datasites_dir = tmp.path().join("datasites");
        fs::create_dir_all(&datasites_dir).unwrap();
        let path = datasites_dir.join("syft.pub.yaml");
        let mut original = SyftPermissions::new_for_datasite("original@example.com");
        original.add_rule(
            "results/**/*",
            vec!["original@example.com".to_string()],
            vec!["original@example.com".to_string()],
        );
        original.save(&path.to_path_buf()).unwrap();

        let storage = SyftBoxStorage::new(tmp.path());
        update_permissions_for_new_recipient(&storage, &path, "colleague@example.com").unwrap();

        let updated: SyftPermissions =
            serde_yaml::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

        // Every rule should grant read access to both recipients
        for rule in &updated.rules {
            assert!(rule
                .access
                .read
                .contains(&"original@example.com".to_string()));
            assert!(rule
                .access
                .read
                .contains(&"colleague@example.com".to_string()));
        }

        // Results rule should also grant write access
        let results_rule = updated
            .rules
            .iter()
            .find(|r| r.pattern == "results/**/*")
            .expect("results rule present");
        assert!(results_rule
            .access
            .write
            .contains(&"colleague@example.com".to_string()));

        // Running again should not duplicate entries
        update_permissions_for_new_recipient(&storage, &path, "colleague@example.com").unwrap();
        let deduped: SyftPermissions =
            serde_yaml::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
        let reads: Vec<_> = deduped
            .rules
            .iter()
            .flat_map(|r| r.access.read.iter())
            .collect();
        assert_eq!(
            reads
                .iter()
                .filter(|v| v.as_str() == "colleague@example.com")
                .count(),
            deduped.rules.len()
        );
    }

    #[test]
    #[serial]
    fn send_project_message_writes_db_and_rpc_request() {
        let tmp = TempDir::new().unwrap();
        let bv_home = tmp.path().join(".biovault");
        fs::create_dir_all(&bv_home).unwrap();
        set_test_biovault_home(&bv_home);

        let syft_dir = tmp.path().join("syft_data");
        fs::create_dir_all(&syft_dir).unwrap();
        set_test_syftbox_data_dir(&syft_dir);

        let config = Config {
            email: "me@example.com".into(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let app = SyftBoxApp::new(&syft_dir, &config.email, "biovault").unwrap();
        let datasite_path = app.data_dir.join("datasites").join(&config.email);

        let submission_folder = "demo-20240101-deadbeef";
        let submission_path = datasite_path
            .join("shared")
            .join("biovault")
            .join("submissions")
            .join(submission_folder);
        app.storage.ensure_dir(&submission_path).unwrap();

        let mut hashes = HashMap::new();
        hashes.insert("workflow.nf".to_string(), "abc123".to_string());
        let project = ProjectYaml {
            name: "Demo Project".into(),
            author: config.email.clone(),
            datasites: Some(vec![config.email.clone()]),
            participants: Some(vec!["participant1".into()]),
            workflow: "workflow.nf".into(),
            template: None,
            assets: Some(vec!["assets/data.csv".into()]),
            inputs: None,
            outputs: None,
            b3_hashes: Some(hashes),
            extra: HashMap::new(),
        };

        write_yaml_to_storage(
            &app.storage,
            &submission_path.join("project.yaml"),
            &project,
            Some(vec![config.email.clone()]),
            true,
        )
        .unwrap();
        let perms = SyftPermissions::new_for_datasite(&config.email);
        write_yaml_to_storage(
            &app.storage,
            &submission_path.join("syft.pub.yaml"),
            &perms,
            None,
            true,
        )
        .unwrap();

        send_project_message(
            &config,
            &project,
            &config.email,
            &submission_path,
            submission_folder,
            "deadbeefcafebabe",
            "2024-01-01",
            true,
            true,
        )
        .unwrap();

        let db_path = get_message_db_path(&config).unwrap();
        let db = MessageDb::new(&db_path).unwrap();
        let messages = db.list_messages(None).unwrap();
        assert_eq!(messages.len(), 1);
        let msg = &messages[0];
        assert_eq!(msg.to, config.email);
        assert!(matches!(msg.message_type, MessageType::Project { .. }));

        let metadata = msg.metadata.as_ref().expect("metadata present");
        let project_location = metadata
            .get("project_location")
            .and_then(|v| v.as_str())
            .unwrap();
        assert!(project_location.ends_with(submission_folder));
        assert!(project_location.starts_with(&format!("syft://{}", config.email)));

        // Ensure RPC request file written for the recipient
        let rpc_dir = app.register_endpoint("/message").unwrap();
        let rpc_entries = app.storage.list_dir(&rpc_dir).unwrap_or_default();
        assert!(!rpc_entries.is_empty());

        clear_test_syftbox_data_dir();
        clear_test_biovault_home();
    }

    #[tokio::test]
    #[serial]
    async fn submit_updates_existing_submission_for_new_recipient() {
        let temp = TempDir::new().unwrap();
        let config_email = "me@example.com";
        let (project_dir, app, _bv_home) = setup_submit_env(&temp, config_email);

        submit(
            project_dir.to_string_lossy().to_string(),
            "friend@example.com".into(),
            true,
            false,
        )
        .await
        .unwrap();

        let before_dirs = list_submissions(&app);
        assert_eq!(before_dirs.len(), 1);
        let submission_path = before_dirs[0].clone();

        submit(
            project_dir.to_string_lossy().to_string(),
            "colleague@example.com".into(),
            true,
            false,
        )
        .await
        .unwrap();

        let after_dirs = list_submissions(&app);
        assert_eq!(after_dirs.len(), 1);
        assert_eq!(after_dirs[0], submission_path);

        let permissions: SyftPermissions = serde_yaml::from_str(
            &app.storage
                .read_plaintext_string(&submission_path.join("syft.pub.yaml"))
                .unwrap(),
        )
        .unwrap();
        for rule in &permissions.rules {
            assert!(rule
                .access
                .read
                .contains(&"colleague@example.com".to_string()));
        }

        let project = ProjectYaml::from_file(&submission_path.join("project.yaml")).unwrap();
        let datasites = project.datasites.unwrap();
        assert_eq!(datasites.len(), 2);
        assert!(datasites.contains(&"friend@example.com".to_string()));
        assert!(datasites.contains(&"colleague@example.com".to_string()));

        let db_path = get_message_db_path(&Config {
            email: config_email.into(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        })
        .unwrap();
        let db = MessageDb::new(&db_path).unwrap();
        assert_eq!(db.list_messages(None).unwrap().len(), 2);

        clear_test_syftbox_data_dir();
        clear_test_biovault_home();
    }

    #[tokio::test]
    #[serial]
    async fn submit_creates_new_version_when_forced_after_changes() {
        let temp = TempDir::new().unwrap();
        let config_email = "force@example.com";
        let (project_dir, app, _bv_home) = setup_submit_env(&temp, config_email);

        submit(
            project_dir.to_string_lossy().to_string(),
            "partner@example.com".into(),
            true,
            false,
        )
        .await
        .unwrap();

        let initial_dirs = list_submissions(&app);
        assert_eq!(initial_dirs.len(), 1);

        // Change workflow to alter project hash
        tokio::time::sleep(Duration::from_millis(5)).await; // ensure timestamp difference if needed
        fs::write(
            project_dir.join("workflow.nf"),
            b"process MAIN { echo 'changed' }",
        )
        .unwrap();
        write_project_yaml(&project_dir.join("project.yaml"), None);

        submit(
            project_dir.to_string_lossy().to_string(),
            "partner@example.com".into(),
            true,
            false,
        )
        .await
        .unwrap();

        // Changed content produces a new submission folder automatically
        let after_no_force = list_submissions(&app);
        assert!(after_no_force.len() >= 2);

        let latest_submission = after_no_force
            .iter()
            .max_by_key(|p| p.metadata().ok().and_then(|m| m.modified().ok()))
            .cloned()
            .unwrap_or_else(|| after_no_force.last().unwrap().clone());
        let latest_name = latest_submission.file_name().unwrap().to_string_lossy();
        assert!(latest_name.contains('-'));

        submit(
            project_dir.to_string_lossy().to_string(),
            "partner@example.com".into(),
            true,
            true,
        )
        .await
        .unwrap();

        let final_dirs = list_submissions(&app);
        assert_eq!(final_dirs.len(), after_no_force.len());

        clear_test_syftbox_data_dir();
        clear_test_biovault_home();
    }
}
