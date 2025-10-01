use crate::cli::commands::run::{execute as run_execute, RunParams};
use crate::cli::syft_url::SyftURL;
use crate::config::Config;
use crate::messages::{Message, MessageDb, MessageSync};
use crate::types::ProjectYaml;
use crate::types::SyftPermissions;
use anyhow::Result;
use colored::Colorize;
use dialoguer::{Confirm, Input, Select};
use serde_json::json;
use std::fs;
use std::path::Path;
use std::path::PathBuf;

const MESSAGE_ENDPOINT: &str = "/message";

/// Clean up stale database locks
pub fn cleanup_locks(config: &Config, all: bool) -> Result<()> {
    use crate::messages::MessageDb;

    if all {
        println!("🧹 Scanning all BioVault virtualenvs for stale locks...");
        // This would be more complex - scan all possible BioVault installations
        println!("⚠️  --all mode not yet implemented. Cleaning current environment only.");
    }

    let db_path = get_message_db_path(config)?;
    let cleaned = MessageDb::clean_stale_lock(&db_path)?;

    if !cleaned {
        println!("✅ No stale locks found");
    }

    Ok(())
}

/// Expand environment variables in text (specifically $SYFTBOX_DATA_DIR)
fn expand_env_vars_in_text(text: &str) -> Result<String> {
    let mut result = text.to_string();

    // Expand $SYFTBOX_DATA_DIR
    if let Ok(data_dir) = std::env::var("SYFTBOX_DATA_DIR") {
        result = result.replace("$SYFTBOX_DATA_DIR", &data_dir);
    }

    Ok(result)
}

/// Get the path to the message database
pub fn get_message_db_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = config.get_biovault_dir()?;
    let db_path = biovault_dir.join("data").join("messages.db");

    // Ensure the data directory exists
    if let Some(parent) = db_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    Ok(db_path)
}

#[cfg(test)]
mod tests_fast_helpers {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn get_message_db_path_creates_parent_dir() {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path().join(".bvtest"));
        let cfg = Config {
            email: "e@example".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };
        let path = get_message_db_path(&cfg).unwrap();
        // Parent dir should exist now
        assert!(path.parent().unwrap().exists());
        assert!(path.ends_with("messages.db"));
        crate::config::clear_test_biovault_home();
    }
}

/// Initialize the message system
pub fn init_message_system(config: &Config) -> Result<(MessageDb, MessageSync)> {
    let db_path = get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    let data_dir = config.get_syftbox_data_dir()?;
    let app = crate::syftbox::SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
    app.register_endpoint(MESSAGE_ENDPOINT)?;

    let sync = MessageSync::new(&db_path, app)?;

    println!("BioVault messaging initialized for {}", config.email);

    Ok((db, sync))
}

/// Send a message
pub fn send_message(
    config: &Config,
    recipient: &str,
    body: &str,
    subject: Option<&str>,
) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync first to check for any pending ACKs
    let _ = sync.sync_quiet();

    // Create the message
    let mut msg = Message::new(
        config.email.clone(),
        recipient.to_string(),
        body.to_string(),
    );

    if let Some(subj) = subject {
        msg.subject = Some(subj.to_string());
    }

    // Save to local database
    db.insert_message(&msg)?;

    // Send via RPC
    sync.send_message(&msg.id)?;

    println!("✉️  Message sent to {}", recipient);
    if let Some(subj) = &msg.subject {
        println!("   Subject: {}", subj);
    }

    Ok(())
}

/// Reply to a message
pub fn reply_message(config: &Config, message_id: &str, body: &str) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync first to ensure we have the latest messages
    let _ = sync.sync_quiet();

    // Get the original message
    let original = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    // Create reply
    let reply = Message::reply_to(&original, config.email.clone(), body.to_string());

    // Save to local database
    db.insert_message(&reply)?;

    // Send via RPC
    sync.send_message(&reply.id)?;

    println!("↩️  Reply sent to {}", reply.to);

    Ok(())
}

/// Delete a message
pub fn delete_message(config: &Config, message_id: &str) -> Result<()> {
    let (db, _) = init_message_system(config)?;

    // First get the message to ensure it exists and get the full ID
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    // Now delete with the full ID
    db.delete_message(&msg.id)?;

    println!(
        "🗑️  Message deleted: {} ({})",
        &msg.id[..8],
        msg.display_subject()
    );

    Ok(())
}

/// List messages
pub fn list_messages(
    config: &Config,
    unread_only: bool,
    sent_only: bool,
    projects_only: bool,
) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync to get latest messages and show notification if new
    let (_new_msg_ids, count) = sync.sync_quiet()?;
    if count > 0 {
        println!("🆕 {} new message(s) received", count);
    }

    let mut messages = if unread_only {
        db.list_unread_messages()?
    } else if sent_only {
        db.list_sent_messages(Some(50))?
    } else {
        db.list_messages(Some(50))?
    };

    // Filter for projects if requested
    if projects_only {
        messages.retain(|msg| {
            matches!(
                msg.message_type,
                crate::messages::MessageType::Project { .. }
            )
        });
    }

    if messages.is_empty() {
        if unread_only {
            println!("No unread messages");
        } else {
            println!("No messages");
        }
        return Ok(());
    }

    println!("\n📬 Messages:");
    println!("─────────────");

    for msg in messages {
        let status_icon = match msg.status {
            crate::messages::MessageStatus::Draft => "📝",
            crate::messages::MessageStatus::Sent => "📤",
            crate::messages::MessageStatus::Received => "📥",
            crate::messages::MessageStatus::Read => "👁️",
            crate::messages::MessageStatus::Deleted => "🗑️",
            crate::messages::MessageStatus::Archived => "📁",
        };

        println!("\n{} [{}]", status_icon, &msg.id[..8]);
        println!("  From: {}", msg.from);
        println!("  To: {}", msg.to);
        println!("  Subject: {}", msg.display_subject());
        // Convert to local time
        let local_time = msg.created_at.with_timezone(&chrono::Local);
        println!("  Date: {}", local_time.format("%Y-%m-%d %H:%M:%S %Z"));

        // Show first 100 chars of body
        let preview = if msg.body.len() > 100 {
            format!("{}...", &msg.body[..100])
        } else {
            msg.body.clone()
        };
        println!("  Body: {}", preview);

        // If this is a project message with metadata, try quick verification
        if let Some(meta) = &msg.metadata {
            if msg.message_type.to_string() == "project" {
                match verify_project_from_metadata(config, meta) {
                    Ok((true, note)) => println!("  Project Verify: OK{}", note),
                    Ok((false, note)) => println!("  Project Verify: FAIL{}", note),
                    Err(e) => println!("  Project Verify: UNVERIFIED ({})", e),
                }
            }
        }

        if msg.parent_id.is_some() {
            println!("  ↩️  Reply to: {}", msg.parent_id.as_ref().unwrap());
        }
    }

    Ok(())
}

/// Read a specific message
pub async fn read_message(config: &Config, message_id: &str) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync first in case there are new messages
    let _ = sync.sync_quiet();

    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    // Mark as read if it was received
    if msg.status == crate::messages::MessageStatus::Received {
        db.mark_as_read(message_id)?;
    }

    println!("\n📧 Message Details");
    println!("═══════════════════");
    println!("ID: {}", msg.id);
    println!("From: {}", msg.from);
    println!("To: {}", msg.to);
    println!("Subject: {}", msg.display_subject());
    let local_time = msg.created_at.with_timezone(&chrono::Local);
    println!("Date: {}", local_time.format("%Y-%m-%d %H:%M:%S %Z"));

    if let Some(parent_id) = &msg.parent_id {
        println!("Reply to: {}", parent_id);
    }

    if let Some(thread_id) = &msg.thread_id {
        println!("Thread: {}", thread_id);
    }

    println!("\nBody:");
    println!("─────");

    // Expand environment variables in message body
    let expanded_body = expand_env_vars_in_text(&msg.body)?;
    println!("{}", expanded_body);

    // If this is a project message with metadata, attempt verification and show details
    if let Some(meta) = &msg.metadata {
        if msg.message_type.to_string() == "project" {
            println!("\nProject Verification:");
            println!("──────────────────────");
            match verify_project_from_metadata(config, meta) {
                Ok((true, note)) => println!("Status: OK{}", note),
                Ok((false, note)) => println!("Status: FAIL{}", note),
                Err(e) => println!("Status: UNVERIFIED ({})", e),
            }

            println!("\nDetails:");
            println!("────────");
            if let Some(loc) = meta.get("project_location").and_then(|v| v.as_str()) {
                println!("Project location: {}", loc.cyan());
                if let Ok(p) = resolve_syft_url_to_path(config, loc) {
                    println!("Local path: {}", p.display());
                }
            }
            if let Some(date) = meta.get("date").and_then(|v| v.as_str()) {
                println!("Date: {}", date);
            }
            if let Some(project) = meta.get("project") {
                if let Some(participants) = project.get("participants").and_then(|v| v.as_array()) {
                    if !participants.is_empty() {
                        let parts: Vec<String> = participants
                            .iter()
                            .filter_map(|p| p.as_str().map(|s| s.to_string()))
                            .collect();
                        println!("Desired participants: {}", parts.join(", "));
                    }
                }
                if let Some(assets) = project.get("assets").and_then(|v| v.as_array()) {
                    if !assets.is_empty() {
                        println!("Assets:");
                        for a in assets {
                            if let Some(s) = a.as_str() {
                                println!("  - {}", s);
                            }
                        }
                    }
                }
            }

            if let Some(status) = meta.get("remote_status").and_then(|v| v.as_str()) {
                println!("Remote status: {}", status);
                if let Some(reason) = meta.get("remote_reason").and_then(|v| v.as_str()) {
                    if !reason.is_empty() {
                        println!("Reason: {}", reason);
                    }
                }
                if status == "approved" {
                    if let Some(path_str) = meta.get("results_path").and_then(|v| v.as_str()) {
                        let results = std::path::PathBuf::from(path_str);
                        println!("Results location: {}", results.display());
                        if results.exists() {
                            println!("Results tree:");
                            print_dir_tree(&results, 3)?;
                        }
                    } else if let Some(loc) = meta.get("project_location").and_then(|v| v.as_str())
                    {
                        if let Ok(root) = resolve_syft_url_to_path(config, loc) {
                            let results = root.join("results");
                            println!("Results location: {}", results.display());
                            if results.exists() {
                                println!("Results tree:");
                                print_dir_tree(&results, 3)?;
                            }
                        }
                    }
                }
            }

            // Offer actions to the recipient (show regardless of read status)
            if msg.to == config.email {
                println!("\nActions:");
                println!("────────");
                let actions = vec![
                    "Reject",
                    "Review",
                    "Approve (run if needed, release results)",
                    "Run on test data",
                    "Run on real data",
                    "Back",
                ];
                let choice = Select::new()
                    .with_prompt("Choose an action")
                    .items(&actions)
                    .default(5)
                    .interact_opt()?;

                if let Some(idx) = choice {
                    match idx {
                        0 => reject_project(config, &msg)?,
                        1 => review_project(config, &msg)?,
                        2 => approve_project(config, &msg).await?,
                        3 => run_project_test(config, &msg).await?,
                        4 => run_project_real(config, &msg).await?,
                        _ => {}
                    }
                }
            }

            // Sender-side archive action after approval to revoke write and mark done
            if msg.from == config.email {
                println!("\nSender Actions:");
                println!("───────────────");
                let actions = vec!["Archive (finalize and revoke write)", "Back"];
                let choice = Select::new()
                    .with_prompt("Choose an action")
                    .items(&actions)
                    .default(1)
                    .interact_opt()?;
                if let Some(0) = choice {
                    archive_project(config, &msg)?;
                }
            }
        }
        // If this is a status update reply, show status and results info directly
        if let crate::messages::MessageType::Request { request_type, .. } = &msg.message_type {
            if request_type == "status" {
                println!("\nStatus Update:");
                println!("──────────────");
                if let Some(update) = meta.get("status_update") {
                    if let Some(status) = update.get("status").and_then(|v| v.as_str()) {
                        println!("Status: {}", status);
                    }
                    if let Some(reason) = update.get("reason").and_then(|v| v.as_str()) {
                        if !reason.is_empty() {
                            println!("Note: {}", reason);
                        }
                    }
                }
                if let Some(path_str) = meta.get("results_path").and_then(|v| v.as_str()) {
                    let results = std::path::PathBuf::from(path_str);
                    println!("Results location: {}", results.display());
                    if results.exists() {
                        println!("Results tree:");
                        print_dir_tree(&results, 3)?;
                    }
                }
            }
        }
    }

    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub enum ProjectAction {
    Reject,
    Review,
    Approve,
    RunTest,
    RunReal,
}

/// Public entrypoint so other commands (like inbox) can trigger project triage actions
pub async fn perform_project_action(
    config: &Config,
    message_id: &str,
    action: ProjectAction,
) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    match action {
        ProjectAction::Reject => reject_project(config, &msg)?,
        ProjectAction::Review => review_project(config, &msg)?,
        ProjectAction::Approve => approve_project(config, &msg).await?,
        ProjectAction::RunTest => run_project_test(config, &msg).await?,
        ProjectAction::RunReal => run_project_real(config, &msg).await?,
    }
    Ok(())
}

/// Process a project message non-interactively (for automated testing)
pub async fn process_project_message(
    config: &Config,
    message_id: &str,
    test: bool,
    _real: bool,
    participant: Option<String>,
    approve: bool,
    _non_interactive: bool, // Currently always non-interactive
) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    // Verify it's a project message
    if !matches!(
        msg.message_type,
        crate::messages::MessageType::Project { .. }
    ) {
        return Err(anyhow::anyhow!(
            "Message {} is not a project message",
            message_id
        ));
    }

    // Build the project copy in private directory
    let dest = build_run_project_copy(config, &msg)?;

    // Determine participant source
    let participant_source = if let Some(ref p) = participant {
        // Normalize the participant source if it's just an ID
        if !test {
            // For real data, normalize the participant ID to full path
            normalize_participant_source_for_real(p)?
        } else {
            // For test data, just use the ID as-is (sample data)
            p.clone()
        }
    } else {
        return Err(anyhow::anyhow!(
            "No participant specified. Please provide a participant source with --participant"
        ));
    };

    // Run the project
    use crate::cli::commands::run::{execute as run_execute, RunParams};

    let result = run_execute(RunParams {
        project_folder: dest.to_string_lossy().to_string(),
        participant_source: participant_source.clone(),
        test,
        download: true,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: Some("snp".to_string()),
        results_dir: None,
    })
    .await;

    match result {
        Ok(_) => {
            println!(
                "✓ Project processed successfully with participant: {}",
                participant_source
            );

            // Show results location
            let results_base = if test { "results-test" } else { "results-real" };
            let results_dir = dest.join(results_base).join(&participant_source);
            if results_dir.exists() {
                println!("Results saved to: {}", results_dir.display());
            }

            // Approve if requested
            if approve {
                approve_project_non_interactive(config, &msg).await?;
                println!("✓ Project approved and results sent to sender");
            }
        }
        Err(e) => {
            eprintln!("✗ Failed to process project: {}", e);
            return Err(e);
        }
    }

    Ok(())
}

/// Archive a message (for non-interactive use)
pub fn archive_message(config: &Config, message_id: &str) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    archive_project(config, &msg)?;
    Ok(())
}

fn archive_project(config: &Config, msg: &Message) -> anyhow::Result<()> {
    // Update syft.pub.yaml by removing the results write rule
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let project_location = meta
        .get("project_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing project_location"))?;
    let root = resolve_syft_url_to_path(config, project_location)?;
    let perm_path = root.join("syft.pub.yaml");
    if perm_path.exists() {
        let content = fs::read_to_string(&perm_path)?;
        let mut perms: SyftPermissions = serde_yaml::from_str(&content)?;
        perms.rules.retain(|r| r.pattern != "results/**/*");
        perms.save(&perm_path)?;
        println!("Revoked write permissions to results for the recipient.");
    } else {
        println!(
            "Warning: permission file not found at {}",
            perm_path.display()
        );
    }

    // Mark the message as archived
    let (db, _) = init_message_system(config)?;
    if let Some(mut full) = db.get_message(&msg.id)? {
        full.status = crate::messages::MessageStatus::Archived;
        db.update_message(&full)?;
        println!("Message archived.");
    }
    Ok(())
}

fn sender_project_root(
    config: &Config,
    meta: &serde_json::Value,
) -> anyhow::Result<(PathBuf, String)> {
    let project_location = meta
        .get("project_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing project_location"))?;
    let path = resolve_syft_url_to_path(config, project_location)?;
    let folder = Path::new(project_location)
        .components()
        .next_back()
        .and_then(|c| match c {
            std::path::Component::Normal(os) => os.to_str(),
            _ => None,
        })
        .unwrap_or("submission")
        .to_string();
    Ok((path, folder))
}

fn receiver_private_submissions_path(config: &Config) -> anyhow::Result<PathBuf> {
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(data_dir
        .join("datasites")
        .join(&config.email)
        .join("private")
        .join("app_data")
        .join("biovault")
        .join("submissions"))
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> anyhow::Result<()> {
    fs::create_dir_all(dst)?;
    for entry in walkdir::WalkDir::new(src)
        .into_iter()
        .filter_map(|e| e.ok())
    {
        let rel = entry.path().strip_prefix(src).unwrap();
        let out = dst.join(rel);
        if entry.file_type().is_dir() {
            fs::create_dir_all(&out)?;
        } else {
            // Skip any syft.pub.yaml files anywhere in the tree
            if entry
                .path()
                .file_name()
                .and_then(|n| n.to_str())
                .map(|n| n == "syft.pub.yaml")
                .unwrap_or(false)
            {
                continue;
            }
            if let Some(parent) = out.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::copy(entry.path(), &out)?;
        }
    }
    Ok(())
}

fn build_run_project_copy(config: &Config, msg: &Message) -> anyhow::Result<PathBuf> {
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, folder_name) = sender_project_root(config, meta)?;
    let dest_root = receiver_private_submissions_path(config)?;
    let dest_path = dest_root.join(&folder_name);
    if !dest_path.exists() {
        copy_dir_recursive(&sender_root, &dest_path)?;
    }
    Ok(dest_path)
}

fn prompt_participant_source(default_source: &str) -> anyhow::Result<String> {
    let input: String = Input::new()
        .with_prompt("Participant source (syft://, file.yaml#fragment, or sample ID)")
        .default(default_source.to_string())
        .interact_text()?;
    Ok(input)
}

fn normalize_participant_source_for_real(input: &str) -> anyhow::Result<String> {
    // If it already looks like a URL/path/fragment, keep as-is
    if input.contains("://")
        || input.contains('/')
        || input.contains(".yaml")
        || input.contains('#')
    {
        return Ok(input.to_string());
    }
    // Otherwise treat as participant ID under local participants.yaml
    let biovault_home = crate::config::get_biovault_home()?;
    let participants_file = biovault_home.join("participants.yaml");
    let file_str = participants_file.to_string_lossy();
    Ok(format!("{}#participants.{}", file_str, input))
}

fn send_status_ack_with_meta(
    config: &Config,
    original: &Message,
    status: &str,
    body: Option<String>,
    extra_metadata: Option<serde_json::Value>,
) -> anyhow::Result<()> {
    let (db, sync) = init_message_system(config)?;

    let mut reply = Message::reply_to(
        original,
        config.email.clone(),
        body.unwrap_or_else(|| status.to_string()),
    );
    reply.subject = Some(format!("Project {}", status));
    reply.message_type = crate::messages::MessageType::Request {
        request_type: "status".to_string(),
        params: None,
    };
    let mut md = json!({
        "status_update": {
            "message_id": original.id,
            "status": status,
            "reason": reply.body,
        }
    });
    if let Some(extra) = extra_metadata {
        if let Some(obj) = md.as_object_mut() {
            if let Some(extra_obj) = extra.as_object() {
                for (k, v) in extra_obj.iter() {
                    obj.insert(k.clone(), v.clone());
                }
            }
        }
    }
    reply.metadata = Some(md);

    db.insert_message(&reply)?;
    let _ = sync.send_message(&reply.id);
    Ok(())
}

fn reject_project(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let custom = Confirm::new()
        .with_prompt("Add a rejection reason?")
        .default(true)
        .interact()
        .unwrap_or(true);
    let body = if custom {
        Some(
            Input::new()
                .with_prompt("Reason")
                .default("Sorry your request has been rejected.".to_string())
                .interact_text()?,
        )
    } else {
        Some("Sorry your request has been rejected.".to_string())
    };
    send_status_ack_with_meta(config, msg, "rejected", body, None)?;
    println!("{}", "Sent rejection to sender.".yellow());
    Ok(())
}

fn review_project(config: &Config, msg: &Message) -> anyhow::Result<()> {
    send_status_ack_with_meta(
        config,
        msg,
        "reviewing",
        Some("Request is under review.".to_string()),
        None,
    )?;
    println!("{}", "Marked as reviewing and notified sender.".yellow());
    Ok(())
}

async fn run_project_test(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_project_copy(config, msg)?;
    let source = prompt_participant_source("NA06985")?;
    let source_for_run = source.clone();
    run_execute(RunParams {
        project_folder: dest.to_string_lossy().to_string(),
        participant_source: source_for_run,
        test: true,
        download: true,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: None,
        results_dir: None,
    })
    .await?;
    println!("{}", "Test run completed.".green());
    // Show results directory and a tree of contents
    print_results_location_and_tree(&dest, &source, true)?;
    Ok(())
}

async fn run_project_real(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_project_copy(config, msg)?;
    let biovault_home = crate::config::get_biovault_home()?;
    let participants_file = biovault_home.join("participants.yaml");
    let default_source = if participants_file.exists() {
        format!("{}#participants.ALL", participants_file.to_string_lossy())
    } else {
        "participants.yaml#participants.ALL".to_string()
    };
    let raw = prompt_participant_source(&default_source)?;
    let source = normalize_participant_source_for_real(&raw)?;
    let source_for_run = source.clone();
    run_execute(RunParams {
        project_folder: dest.to_string_lossy().to_string(),
        participant_source: source_for_run,
        test: false,
        download: false,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: None,
        results_dir: None,
    })
    .await?;
    println!("{}", "Real data run completed.".green());
    // Show results directory and a tree of contents
    print_results_location_and_tree(&dest, &source, false)?;
    Ok(())
}

/// Non-interactive version of approve_project for automated testing
async fn approve_project_non_interactive(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_project_copy(config, msg)?;
    let results_dir = dest.join("results-real");
    let needs_run = !results_dir.exists()
        || fs::read_dir(&results_dir)
            .map(|mut i| i.next().is_none())
            .unwrap_or(true);
    if needs_run {
        println!("No results found. Running on real data before approval...");
        run_project_real(config, msg).await?;
    }

    // Release results to sender shared location
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, _folder) = sender_project_root(config, meta)?;
    let sender_results = sender_root.join("results");
    copy_dir_recursive(&results_dir, &sender_results)?;

    // Get project details for the message
    let project_location = meta
        .get("project_location")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown");

    // Check if results exist and get basic info
    let results_info = if sender_results.exists() {
        let mut info = Vec::new();
        if let Ok(entries) = fs::read_dir(&sender_results) {
            for entry in entries.flatten() {
                if entry.file_type().map(|ft| ft.is_dir()).unwrap_or(false) {
                    info.push(format!("{}/", entry.file_name().to_string_lossy()));
                } else {
                    info.push(entry.file_name().to_string_lossy().to_string());
                }
            }
        }
        if info.is_empty() {
            "No results found".to_string()
        } else {
            format!("Results:\n  {}", info.join("\n  "))
        }
    } else {
        "Results directory not found".to_string()
    };

    // Create relative path for results location
    let data_dir = config.get_syftbox_data_dir()?;
    let relative_results = if sender_results.starts_with(&data_dir) {
        format!(
            "$SYFTBOX_DATA_DIR/{}",
            sender_results.strip_prefix(&data_dir)?.to_string_lossy()
        )
    } else {
        sender_results.to_string_lossy().to_string()
    };

    // Create the default approval message with results location
    let body = format!(
        "Your project has been approved.\n\nProject location: {}\nResults location: {}\n\n{}",
        project_location, relative_results, results_info
    );

    // Include results_path in the status update so sender can see it directly
    let extra = json!({ "results_path": relative_results });
    send_status_ack_with_meta(config, msg, "approved", Some(body), Some(extra))?;
    println!(
        "{}",
        "Approved, results released, and sender notified.".green()
    );
    Ok(())
}

async fn approve_project(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_project_copy(config, msg)?;
    let results_dir = dest.join("results-real");
    let needs_run = !results_dir.exists()
        || fs::read_dir(&results_dir)
            .map(|mut i| i.next().is_none())
            .unwrap_or(true);
    if needs_run {
        println!("No results found. Running on real data before approval...");
        run_project_real(config, msg).await?;
    }

    // Release results to sender shared location
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, _folder) = sender_project_root(config, meta)?;
    let sender_results = sender_root.join("results");
    copy_dir_recursive(&results_dir, &sender_results)?;

    // Get project details for the message
    let project_location = meta
        .get("project_location")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown");

    // Check if results exist and get basic info
    let results_info = if sender_results.exists() {
        let mut info = Vec::new();
        if let Ok(entries) = fs::read_dir(&sender_results) {
            for entry in entries.flatten() {
                if entry.file_type().map(|ft| ft.is_dir()).unwrap_or(false) {
                    info.push(format!("{}/", entry.file_name().to_string_lossy()));
                } else {
                    info.push(entry.file_name().to_string_lossy().to_string());
                }
            }
        }
        if info.is_empty() {
            "No results found".to_string()
        } else {
            format!("Results:\n  {}", info.join("\n  "))
        }
    } else {
        "Results directory not found".to_string()
    };

    // Create relative path for results location
    let data_dir = config.get_syftbox_data_dir()?;
    let relative_results = if sender_results.starts_with(&data_dir) {
        format!(
            "$SYFTBOX_DATA_DIR/{}",
            sender_results.strip_prefix(&data_dir)?.to_string_lossy()
        )
    } else {
        sender_results.to_string_lossy().to_string()
    };

    // Create the default approval message with results location
    let default_message = format!(
        "Your project has been approved.\n\nProject location: {}\nResults location: {}\n\n{}",
        project_location, relative_results, results_info
    );

    // Optional approval message
    let add_note = Confirm::new()
        .with_prompt("Add a message to approval?")
        .default(false)
        .interact()
        .unwrap_or(false);
    let body = if add_note {
        Some(
            Input::new()
                .with_prompt("Message")
                .default(default_message.clone())
                .interact_text()?,
        )
    } else {
        Some(default_message)
    };

    // Include results_path in the status update so sender can see it directly
    let extra = json!({ "results_path": relative_results });
    send_status_ack_with_meta(config, msg, "approved", body, Some(extra))?;
    println!(
        "{}",
        "Approved, results released, and sender notified.".green()
    );
    Ok(())
}

fn extract_participant_id(source: &str) -> Option<String> {
    if let Some(pos) = source.find('#') {
        let frag = &source[pos + 1..];
        if let Some(rest) = frag.strip_prefix("participants.") {
            return Some(rest.to_string());
        }
    }
    // If it doesn't look like a path/URL and has no fragment, treat the whole string as the ID
    if !(source.contains("://")
        || source.contains('/')
        || source.contains(".yaml")
        || source.contains('#'))
    {
        return Some(source.to_string());
    }
    None
}

fn print_results_location_and_tree(
    project_root: &Path,
    participant_source: &str,
    is_test: bool,
) -> anyhow::Result<()> {
    let id = extract_participant_id(participant_source).unwrap_or_else(|| "ALL".to_string());
    let base = if is_test {
        "results-test"
    } else {
        "results-real"
    };
    let results_dir = project_root.join(base).join(&id);
    println!("Results location: {}", results_dir.display());
    if results_dir.exists() {
        println!("Results tree:");
        print_dir_tree(&results_dir, 3)?;
    } else {
        println!("Results folder not found yet at {}", results_dir.display());
    }
    Ok(())
}

fn print_dir_tree(root: &Path, max_depth: usize) -> anyhow::Result<()> {
    fn walk(dir: &Path, depth: usize, max_depth: usize) -> anyhow::Result<()> {
        if depth > max_depth {
            return Ok(());
        }
        let mut entries: Vec<_> = fs::read_dir(dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .collect();
        entries.sort();
        for path in entries {
            let indent = "  ".repeat(depth.saturating_sub(0));
            let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if path.is_dir() {
                println!("{}{}/", indent, name);
                walk(&path, depth + 1, max_depth)?;
            } else {
                println!("{}{}", indent, name);
            }
        }
        Ok(())
    }
    println!("{}/", root.display());
    walk(root, 1, max_depth)
}

/// View a message thread
pub fn view_thread(config: &Config, thread_id: &str) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync to get latest messages in thread
    let _ = sync.sync_quiet();

    let messages = db.get_thread_messages(thread_id)?;

    if messages.is_empty() {
        println!("No messages found in thread: {}", thread_id);
        return Ok(());
    }

    println!("\n💬 Thread: {}", thread_id);
    println!("═══════════════════════════");

    for msg in messages {
        let local_time = msg.created_at.with_timezone(&chrono::Local);
        println!("\n[{}] {}", local_time.format("%Y-%m-%d %H:%M"), msg.from);

        if let Some(subj) = &msg.subject {
            println!("Subject: {}", subj);
        }

        // Expand environment variables in message body
        let expanded_body = expand_env_vars_in_text(&msg.body)?;
        println!("{}", expanded_body);
        println!("─────────────────────");
    }

    Ok(())
}

/// Sync messages (check for new incoming and update ACKs)
pub fn sync_messages(config: &Config) -> Result<()> {
    let (_, sync) = init_message_system(config)?;

    println!("🔄 Syncing messages...");
    sync.sync()?;

    Ok(())
}

/// Resolve a syft:// URL to a local filesystem path within the SyftBox data dir
fn resolve_syft_url_to_path(config: &Config, url: &str) -> anyhow::Result<PathBuf> {
    let parsed = SyftURL::parse(url)?;
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(data_dir
        .join("datasites")
        .join(parsed.email)
        .join(parsed.path))
}

/// Verify a project message using embedded metadata
/// Returns (is_ok, extra_note)
fn verify_project_from_metadata(
    config: &Config,
    metadata: &serde_json::Value,
) -> anyhow::Result<(bool, String)> {
    // Extract the embedded project.yaml
    let project_val = metadata
        .get("project")
        .ok_or_else(|| anyhow::anyhow!("missing project metadata"))?;

    let project: ProjectYaml = serde_json::from_value(project_val.clone())
        .map_err(|e| anyhow::anyhow!("invalid project metadata: {}", e))?;

    // Extract project location syft URL
    let project_location = metadata
        .get("project_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing project_location syft URL"))?;

    let root = resolve_syft_url_to_path(config, project_location)?;
    if !root.exists() {
        return Ok((false, format!(" (missing path: {})", root.display())));
    }

    // Need b3_hashes to verify
    let Some(expected_hashes) = project.b3_hashes.clone() else {
        return Ok((false, " (no hashes provided)".to_string()));
    };

    // Verify each expected file hash
    let mut mismatches: Vec<String> = Vec::new();
    for (rel, expected) in expected_hashes.iter() {
        let file_path = root.join(rel);
        if !file_path.exists() {
            mismatches.push(format!("missing: {}", rel));
            continue;
        }
        let bytes = std::fs::read(&file_path)
            .map_err(|e| anyhow::anyhow!("failed to read {}: {}", file_path.display(), e))?;
        let got = blake3::hash(&bytes).to_hex().to_string();
        if got != *expected {
            mismatches.push(format!("mismatch: {}", rel));
        }
    }

    if mismatches.is_empty() {
        Ok((true, String::new()))
    } else {
        let preview = if mismatches.len() > 3 {
            format!(" ({} issues; first: {})", mismatches.len(), mismatches[0])
        } else {
            format!(" ({})", mismatches.join(", "))
        };
        Ok((false, preview))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use tempfile::TempDir;

    fn create_test_config() -> Config {
        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        }
    }

    #[test]
    fn test_init_message_system() -> Result<()> {
        let temp_dir = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(temp_dir.path());
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault_test"));
        let config = create_test_config();

        // Initialize the message system
        let db_path = get_message_db_path(&config)?;
        let db = MessageDb::new(&db_path)?;

        // Test that we can list messages (should be empty in fresh test DB)
        let messages = db.list_messages(None)?;
        assert_eq!(messages.len(), 0);

        Ok(())
    }

    #[test]
    fn test_message_crud() -> Result<()> {
        let temp_dir = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(temp_dir.path());
        crate::config::set_test_biovault_home(temp_dir.path().join(".biovault_test"));
        let config = create_test_config();

        // Initialize just the database, not the full system (to avoid sync)
        let db_path = get_message_db_path(&config)?;
        let db = MessageDb::new(&db_path)?;

        // Create a message
        let msg = Message::new(
            "test@example.com".to_string(),
            "recipient@example.com".to_string(),
            "Test message body".to_string(),
        );

        // Insert
        db.insert_message(&msg)?;

        // Read
        let retrieved = db.get_message(&msg.id)?;
        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().body, "Test message body");

        // List
        let messages = db.list_messages(None)?;
        assert_eq!(messages.len(), 1);

        // Delete
        db.delete_message(&msg.id)?;

        // Verify it's marked as deleted (soft delete)
        let deleted_msg = db.get_message(&msg.id)?;
        assert!(deleted_msg.is_some());
        assert_eq!(
            deleted_msg.unwrap().status,
            crate::messages::MessageStatus::Deleted
        );

        // Verify it doesn't show in normal list
        let messages_after_delete = db.list_messages(None)?;
        assert_eq!(messages_after_delete.len(), 0);

        Ok(())
    }

    #[test]
    fn normalize_participant_source_for_real_behaviour() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));

        // Existing URL/path strings are preserved
        assert_eq!(
            normalize_participant_source_for_real("syft://x/y#z")?,
            "syft://x/y#z"
        );
        assert_eq!(
            normalize_participant_source_for_real("/abs/path/participants.yaml#p.ID")?,
            "/abs/path/participants.yaml#p.ID"
        );

        // Bare ID resolves to participants.yaml under BIOVAULT home
        let got = normalize_participant_source_for_real("TESTID")?;
        assert!(got.ends_with("participants.yaml#participants.TESTID"));
        Ok(())
    }

    #[test]
    fn resolve_and_paths_and_copy_dir() -> Result<()> {
        let tmp = TempDir::new()?;
        // SyftBox data dir and config email
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = Config {
            email: "u@example.com".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };

        // Build a fake project tree under datasites/u@example.com/app_data/biovault/submissions/proj1
        let root = tmp
            .path()
            .join("datasites/u@example.com/app_data/biovault/submissions/proj1");
        std::fs::create_dir_all(root.join("a/b")).unwrap();
        std::fs::write(root.join("a/b/file.txt"), b"hi").unwrap();
        // This file should be skipped when copying
        std::fs::write(root.join("a/syft.pub.yaml"), b"rules:").unwrap();

        // sender_project_root from metadata
        let loc = format!(
            "syft://u@example.com/app_data/biovault/submissions/{}",
            "proj1"
        );
        let meta = json!({
            "project_location": loc
        });
        let (sender_root, folder) = sender_project_root(&cfg, &meta)?;
        assert_eq!(folder, "proj1");
        assert!(sender_root.ends_with("proj1"));

        // receiver path
        let recv = receiver_private_submissions_path(&cfg)?;
        assert!(recv.ends_with("datasites/u@example.com/private/app_data/biovault/submissions"));

        // Copy behaviour: skip syft.pub.yaml and recreate tree
        let dest = recv.join(&folder);
        copy_dir_recursive(&sender_root, &dest)?;
        assert!(dest.join("a/b/file.txt").exists());
        assert!(!dest.join("a/syft.pub.yaml").exists());

        Ok(())
    }

    #[test]
    fn build_run_project_copy_creates_dest_once() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = Config {
            email: "u@example.com".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };

        // Make sender tree and one file
        let proj = tmp
            .path()
            .join("datasites/u@example.com/app_data/biovault/submissions/projX");
        std::fs::create_dir_all(proj.join("dir")).unwrap();
        std::fs::write(proj.join("dir/f.txt"), b"x").unwrap();
        let meta =
            json!({"project_location": "syft://u@example.com/app_data/biovault/submissions/projX"});
        let mut msg = Message::new("u@example.com".into(), "v@example.com".into(), "b".into());
        msg.metadata = Some(meta);

        let dest1 = build_run_project_copy(&cfg, &msg)?;
        assert!(dest1.join("dir/f.txt").exists());
        // Call again; should not error and keep same path
        let dest2 = build_run_project_copy(&cfg, &msg)?;
        assert_eq!(dest1, dest2);
        Ok(())
    }

    #[test]
    fn verify_project_from_metadata_ok_and_fail() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = Config {
            email: "u@example.com".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };

        // Build project root with two files and compute blake3
        let root = tmp
            .path()
            .join("datasites/u@example.com/app_data/biovault/submissions/projZ");
        std::fs::create_dir_all(&root).unwrap();
        let f1 = root.join("a.txt");
        let f2 = root.join("b.txt");
        std::fs::write(&f1, b"A").unwrap();
        std::fs::write(&f2, b"B").unwrap();
        let h1 = blake3::hash(&std::fs::read(&f1).unwrap())
            .to_hex()
            .to_string();
        let h2 = blake3::hash(&std::fs::read(&f2).unwrap())
            .to_hex()
            .to_string();

        // Construct embedded project.yaml equivalent in metadata
        let project = json!({
            "name": "n", "author":"a", "workflow":"w",
            "b3_hashes": {"a.txt": h1, "b.txt": h2}
        });
        let meta_ok = json!({
            "project": project,
            "project_location": "syft://u@example.com/app_data/biovault/submissions/projZ"
        });
        let (ok, note) = verify_project_from_metadata(&cfg, &meta_ok)?;
        assert!(ok, "expected OK, got note: {}", note);

        // Now break one file and expect failure with note
        std::fs::write(&f2, b"BROKEN").unwrap();
        let (ok2, note2) = verify_project_from_metadata(&cfg, &meta_ok)?;
        assert!(!ok2);
        assert!(!note2.is_empty());
        Ok(())
    }

    #[test]
    fn list_messages_displays_without_actions() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = Config {
            email: "me@example.com".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };
        // Init DB only
        let db_path = get_message_db_path(&cfg)?;
        let db = MessageDb::new(&db_path)?;
        // Insert a simple draft message from me -> other (so no interactive recipient actions)
        let m = Message::new(
            "me@example.com".into(),
            "you@example.com".into(),
            "body".into(),
        );
        db.insert_message(&m)?;
        // Should render list without interacting
        list_messages(&cfg, false, false, false)?;
        Ok(())
    }

    #[tokio::test]
    async fn read_message_marks_received_as_read() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = Config {
            email: "me@example.com".into(),
            syftbox_config: None,
            version: None,

            binary_paths: None,
        };
        let db_path = get_message_db_path(&cfg)?;
        let db = MessageDb::new(&db_path)?;

        let mut m = Message::new(
            "you@example.com".into(),
            "me@example.com".into(),
            "hi".into(),
        );
        m.status = crate::messages::MessageStatus::Received;
        db.insert_message(&m)?;

        read_message(&cfg, &m.id).await?;
        let updated = db.get_message(&m.id)?.unwrap();
        assert_eq!(updated.status, crate::messages::MessageStatus::Read);
        Ok(())
    }

    #[test]
    fn send_and_delete_message_flow() -> Result<()> {
        let tmp = TempDir::new()?;
        // Point syftbox data dir and biovault home to temp
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = create_test_config();

        // Send a message to another recipient
        super::send_message(&cfg, "you@example.com", "hello", Some("subj"))?;

        // Validate DB has a sent message
        let db_path = get_message_db_path(&cfg)?;
        let db = MessageDb::new(&db_path)?;
        let msgs = db.list_sent_messages(None)?;
        assert!(!msgs.is_empty());
        let first = &msgs[0];
        assert_eq!(first.status, crate::messages::MessageStatus::Sent);
        assert_eq!(first.to, "you@example.com");

        // Delete the message (soft delete)
        super::delete_message(&cfg, &first.id)?;
        let after = db.get_message(&first.id)?.unwrap();
        assert_eq!(after.status, crate::messages::MessageStatus::Deleted);
        Ok(())
    }

    #[test]
    fn reply_message_creates_response() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = create_test_config();
        // Initialize DB and insert an incoming message to reply to
        let db_path = get_message_db_path(&cfg)?;
        let db = MessageDb::new(&db_path)?;
        let mut original = Message::new("alice@example.com".into(), cfg.email.clone(), "hi".into());
        original.status = crate::messages::MessageStatus::Received;
        db.insert_message(&original)?;

        super::reply_message(&cfg, &original.id, "re: hi")?;

        // Should have at least one sent message now
        let sent = db.list_sent_messages(None)?;
        assert!(!sent.is_empty());
        Ok(())
    }

    #[test]
    fn list_messages_smoke() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = create_test_config();

        // With empty DB
        super::list_messages(&cfg, false, false, false)?;
        super::list_messages(&cfg, true, false, false)?;
        Ok(())
    }
}
