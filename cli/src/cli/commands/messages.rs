use crate::cli::commands::run::{execute as run_execute, RunParams};
use crate::cli::commands::run_dynamic;
use crate::cli::syft_url::SyftURL;
use crate::config::Config;
use crate::data::{self, BioVaultDb};
use crate::messages::{Message, MessageDb, MessageSync};
use crate::module_spec::ModuleSpec;
use crate::subscriptions;
use crate::syftbox::storage::SyftBoxStorage;
use crate::types::ModuleYaml;
use crate::types::SyftPermissions;
use anyhow::{Context as _, Result};
use colored::Colorize;
use csv::Writer;
use dialoguer::{Confirm, Input, Select};
use serde_json::json;
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};
use std::thread;
use std::time::{Duration, Instant};

const MESSAGE_ENDPOINT: &str = "/message";

fn init_printed_emails() -> &'static Mutex<HashSet<String>> {
    static SET: OnceLock<Mutex<HashSet<String>>> = OnceLock::new();
    SET.get_or_init(|| Mutex::new(HashSet::new()))
}

fn should_log_init_line() -> bool {
    matches!(
        std::env::var("BIOVAULT_MESSAGE_INIT_LOG")
            .unwrap_or_default()
            .to_ascii_lowercase()
            .as_str(),
        "1" | "true" | "yes" | "on"
    )
}

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

fn syftbox_storage(config: &Config) -> Result<SyftBoxStorage> {
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(SyftBoxStorage::new(&data_dir))
}

#[cfg(test)]
mod tests_fast_helpers {
    use super::*;
    use std::sync::{Mutex, OnceLock};
    use tempfile::TempDir;

    fn env_lock() -> &'static Mutex<()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    #[test]
    fn get_message_db_path_creates_parent_dir() {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path().join(".bvtest"));
        let cfg = Config {
            email: "e@example".into(),
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
        let path = get_message_db_path(&cfg).unwrap();
        // Parent dir should exist now
        assert!(path.parent().unwrap().exists());
        assert!(path.ends_with("messages.db"));
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_expand_env_vars_in_text_with_syftbox() {
        let _guard = env_lock().lock().unwrap();
        std::env::set_var("SYFTBOX_DATA_DIR", "/test/data");
        let result = expand_env_vars_in_text("Path: $SYFTBOX_DATA_DIR/file").unwrap();
        assert_eq!(result, "Path: /test/data/file");
        std::env::remove_var("SYFTBOX_DATA_DIR");
    }

    #[test]
    fn test_expand_env_vars_in_text_no_env() {
        let _guard = env_lock().lock().unwrap();
        std::env::remove_var("SYFTBOX_DATA_DIR");
        let result = expand_env_vars_in_text("Plain text").unwrap();
        assert_eq!(result, "Plain text");
    }

    #[test]
    fn test_cleanup_locks_no_stale_locks() {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = Config {
            email: "test@example.com".into(),
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
        // Create db first
        let _ = get_message_db_path(&cfg).unwrap();
        // Should not panic
        let result = cleanup_locks(&cfg, false);
        assert!(result.is_ok());
        crate::config::clear_test_biovault_home();
    }
}

/// Initialize the message system
pub fn init_message_system(config: &Config) -> Result<(MessageDb, MessageSync)> {
    let (db, sync) = build_message_system(config)?;

    if should_log_init_line() {
        let mut seen = init_printed_emails()
            .lock()
            .expect("message init tracker poisoned");
        if seen.insert(config.email.clone()) {
            println!("BioVault messaging initialized for {}", config.email);
        }
    }

    Ok((db, sync))
}

/// Initialize the message system without printing status output
pub fn init_message_system_quiet(config: &Config) -> Result<(MessageDb, MessageSync)> {
    build_message_system(config)
}

fn build_message_system(config: &Config) -> Result<(MessageDb, MessageSync)> {
    let db_path = get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    let data_dir = config.get_syftbox_data_dir()?;
    let app = crate::syftbox::SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
    app.register_endpoint(MESSAGE_ENDPOINT)?;

    let sync = MessageSync::new(&db_path, app)?;

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
    modules_only: bool,
    json_output: bool,
) -> Result<()> {
    let (db, sync) = if json_output {
        init_message_system_quiet(config)?
    } else {
        init_message_system(config)?
    };

    // Quietly sync to get latest messages and show notification if new
    let (_new_msg_ids, count) = sync.sync_quiet()?;
    if count > 0 && !json_output {
        println!("🆕 {} new message(s) received", count);
    }

    let mut messages = if unread_only {
        db.list_unread_messages()?
    } else if sent_only {
        db.list_sent_messages(Some(50))?
    } else {
        db.list_messages(Some(50))?
    };

    // Filter for modules if requested
    if modules_only {
        messages.retain(|msg| {
            matches!(
                msg.message_type,
                crate::messages::MessageType::Module { .. }
            )
        });
    }

    if json_output {
        let json_messages: Vec<_> = messages
            .into_iter()
            .map(|msg| {
                serde_json::json!({
                    "id": msg.id,
                    "from": msg.from,
                    "to": msg.to,
                    "subject": msg.display_subject(),
                    "status": msg.status.to_string(),
                    "message_type": msg.message_type.to_string(),
                    "created_at": msg.created_at.to_rfc3339(),
                })
            })
            .collect();

        println!("{}", serde_json::to_string_pretty(&json_messages)?);
        return Ok(());
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

        // If this is a module message with metadata, try quick verification
        if let Some(meta) = &msg.metadata {
            if msg.message_type.to_string() == "module" {
                match verify_module_from_metadata(config, meta) {
                    Ok((true, note)) => println!("  Module Verify: OK{}", note),
                    Ok((false, note)) => println!("  Module Verify: FAIL{}", note),
                    Err(e) => println!("  Module Verify: UNVERIFIED ({})", e),
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
pub async fn read_message(config: &Config, message_id: &str, non_interactive: bool) -> Result<()> {
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

    // If this is a module message with metadata, attempt verification and show details
    if let Some(meta) = &msg.metadata {
        if msg.message_type.to_string() == "module" {
            println!("\nModule Verification:");
            println!("──────────────────────");
            match verify_module_from_metadata(config, meta) {
                Ok((true, note)) => println!("Status: OK{}", note),
                Ok((false, note)) => println!("Status: FAIL{}", note),
                Err(e) => println!("Status: UNVERIFIED ({})", e),
            }

            println!("\nDetails:");
            println!("────────");
            if let Some(loc) = meta.get("module_location").and_then(|v| v.as_str()) {
                println!("Module location: {}", loc.cyan());
                if let Ok(p) = resolve_syft_url_to_path(config, loc) {
                    println!("Local path: {}", p.display());
                }
            }
            if let Some(date) = meta.get("date").and_then(|v| v.as_str()) {
                println!("Date: {}", date);
            }
            if let Some(module) = meta.get("module") {
                if let Some(participants) = module.get("participants").and_then(|v| v.as_array()) {
                    if !participants.is_empty() {
                        let parts: Vec<String> = participants
                            .iter()
                            .filter_map(|p| p.as_str().map(|s| s.to_string()))
                            .collect();
                        println!("Desired participants: {}", parts.join(", "));
                    }
                }
                if let Some(assets) = module.get("assets").and_then(|v| v.as_array()) {
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
                    let storage_opt = syftbox_storage(config).ok();
                    if let Some(path_str) = meta.get("results_path").and_then(|v| v.as_str()) {
                        let results = std::path::PathBuf::from(path_str);
                        println!("Results location: {}", results.display());
                        if results.exists() {
                            println!("Results tree:");
                            print_dir_tree(storage_opt.as_ref(), &results, 3)?;
                        }
                    } else if let Some(loc) = meta.get("module_location").and_then(|v| v.as_str()) {
                        if let Ok(root) = resolve_syft_url_to_path(config, loc) {
                            let results = root.join("results");
                            println!("Results location: {}", results.display());
                            if results.exists() {
                                println!("Results tree:");
                                print_dir_tree(storage_opt.as_ref(), &results, 3)?;
                            }
                        }
                    }
                }
            }

            // Offer actions to the recipient (show regardless of read status)
            if msg.to == config.email && !non_interactive {
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
                        0 => reject_module(config, &msg)?,
                        1 => review_module(config, &msg)?,
                        2 => approve_module(config, &msg).await?,
                        3 => run_module_test(config, &msg).await?,
                        4 => run_module_real(config, &msg).await?,
                        _ => {}
                    }
                }
            }

            // Sender-side archive action after approval to revoke write and mark done
            if msg.from == config.email && !non_interactive {
                println!("\nSender Actions:");
                println!("───────────────");
                let actions = vec!["Archive (finalize and revoke write)", "Back"];
                let choice = Select::new()
                    .with_prompt("Choose an action")
                    .items(&actions)
                    .default(1)
                    .interact_opt()?;
                if let Some(0) = choice {
                    archive_module(config, &msg)?;
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
                        let storage_opt = syftbox_storage(config).ok();
                        print_dir_tree(storage_opt.as_ref(), &results, 3)?;
                    }
                }
            }
        }
    }

    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub enum ModuleAction {
    Reject,
    Review,
    Approve,
    RunTest,
    RunReal,
}

/// Public entrypoint so other commands (like inbox) can trigger module triage actions
pub async fn perform_module_action(
    config: &Config,
    message_id: &str,
    action: ModuleAction,
) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    if let Some(meta) = msg.metadata.as_ref() {
        if let Some(loc) = meta.get("module_location").and_then(|v| v.as_str()) {
            if let Err(err) = add_subscription_for_syft_url(config, loc) {
                eprintln!("⚠️  Warning: failed to add subscription: {}", err);
            }
        }
    }

    match action {
        ModuleAction::Reject => reject_module(config, &msg)?,
        ModuleAction::Review => review_module(config, &msg)?,
        ModuleAction::Approve => approve_module(config, &msg).await?,
        ModuleAction::RunTest => run_module_test(config, &msg).await?,
        ModuleAction::RunReal => run_module_real(config, &msg).await?,
    }
    Ok(())
}

/// Process a module message non-interactively (for automated testing)
pub async fn process_module_message(
    config: &Config,
    message_id: &str,
    test: bool,
    _real: bool,
    participant: Option<String>,
    approve: bool,
    _non_interactive: bool, // Currently always non-interactive
    sync: bool,
    sync_timeout: Option<u64>,
) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    if let Some(meta) = msg.metadata.as_ref() {
        if let Some(loc) = meta.get("module_location").and_then(|v| v.as_str()) {
            if let Err(err) = add_subscription_for_syft_url(config, loc) {
                eprintln!("⚠️  Warning: failed to add subscription: {}", err);
            }
            if sync {
                // NOTE: This waits on files declared in module metadata. Flow request
                // submissions can include multiple modules/assets outside this manifest,
                // so they may still race until flow-level manifests are added.
                wait_for_submission_sync(config, loc, meta, sync_timeout.unwrap_or(60))?;
            }
        }
    }

    // Verify it's a module message
    if !matches!(
        msg.message_type,
        crate::messages::MessageType::Module { .. }
    ) {
        return Err(anyhow::anyhow!(
            "Message {} is not a module message",
            message_id
        ));
    }

    // Build the module copy in private directory
    let dest = build_run_module_copy(config, &msg)?;

    // Load spec to determine runtime handling (dynamic vs legacy)
    let module_spec_path = dest.join("module.yaml");
    let spec = ModuleSpec::load(&module_spec_path)?;

    match spec.runtime.as_deref() {
        Some("nextflow") | Some("dynamic-nextflow") => {
            return process_dynamic_module_message(
                config,
                &msg,
                &dest,
                participant.clone(),
                test,
                approve,
            )
            .await;
        }
        Some("shell") => {
            return process_shell_module_message(config, &msg, &dest, test, approve).await;
        }
        _ => {}
    }

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

    // Run the module
    use crate::cli::commands::run::{execute as run_execute, RunParams};

    let result = run_execute(RunParams {
        module_folder: dest.to_string_lossy().to_string(),
        participant_source: participant_source.clone(),
        test,
        download: true,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: Some("snp".to_string()),
        results_dir: None,
        nextflow_args: vec![],
    })
    .await;

    match result {
        Ok(_) => {
            println!(
                "✓ Module processed successfully with participant: {}",
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
                approve_module_non_interactive(config, &msg, test).await?;
                println!("✓ Module approved and results sent to sender");
            }
        }
        Err(e) => {
            eprintln!("✗ Failed to process module: {}", e);
            return Err(e);
        }
    }

    Ok(())
}

fn wait_for_submission_sync(
    config: &Config,
    url: &str,
    metadata: &serde_json::Value,
    timeout_secs: u64,
) -> anyhow::Result<()> {
    let root = resolve_syft_url_to_path(config, url)?;
    let module_yaml = root.join("module.yaml");
    let expected_paths = expected_submission_files(metadata);
    let start = Instant::now();

    while start.elapsed() < Duration::from_secs(timeout_secs) {
        if module_yaml.exists() && expected_paths.iter().all(|rel| root.join(rel).exists()) {
            return Ok(());
        }
        thread::sleep(Duration::from_secs(1));
    }

    anyhow::bail!(
        "Timed out waiting for submission sync (missing files under {} after {}s)",
        root.display(),
        timeout_secs
    )
}

fn expected_submission_files(metadata: &serde_json::Value) -> Vec<String> {
    let mut files: Vec<String> = Vec::new();

    if let Some(module) = metadata.get("module") {
        if let Some(hashes) = module.get("b3_hashes").and_then(|v| v.as_object()) {
            for key in hashes.keys() {
                if !key.trim().is_empty() {
                    files.push(key.to_string());
                }
            }
        }
    }

    if !files.iter().any(|p| p == "module.yaml") {
        files.push("module.yaml".to_string());
    }

    files.sort();
    files.dedup();
    files
}

fn add_subscription_for_syft_url(config: &Config, url: &str) -> anyhow::Result<()> {
    let parsed = SyftURL::parse(url)?;
    let data_dir = config.get_syftbox_data_dir()?;
    let syftsub_path = data_dir.join(".data").join("syft.sub.yaml");

    let mut cfg =
        subscriptions::load(&syftsub_path).unwrap_or_else(|_| subscriptions::default_config());

    let mut path = parsed.path.trim_start_matches('/').to_string();
    let path_obj = Path::new(&path);
    let last_segment_is_file = path_obj
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| !ext.is_empty())
        .unwrap_or(false);
    if last_segment_is_file {
        if let Some(parent) = path_obj.parent() {
            path = parent.to_string_lossy().to_string();
        }
    }
    if path.is_empty() {
        return Ok(());
    }
    let rule_path = format!("{}/**", path.trim_end_matches('/'));

    let exists = cfg.rules.iter().any(|rule| {
        rule.action == subscriptions::Action::Allow
            && rule
                .datasite
                .as_deref()
                .map(|ds| ds.eq_ignore_ascii_case(&parsed.email))
                .unwrap_or(false)
            && rule.path == rule_path
    });

    if !exists {
        cfg.rules.push(subscriptions::Rule {
            action: subscriptions::Action::Allow,
            datasite: Some(parsed.email.clone()),
            path: rule_path.clone(),
        });
        subscriptions::save(&syftsub_path, &cfg)?;
        println!(
            "🔧 Added subscription: datasite={} path={} (file: {})",
            parsed.email,
            rule_path,
            syftsub_path.display()
        );
    } else {
        println!(
            "ℹ️  Subscription already present: datasite={} path={} (file: {})",
            parsed.email,
            rule_path,
            syftsub_path.display()
        );
    }

    Ok(())
}

/// Archive a message (for non-interactive use)
pub fn archive_message(config: &Config, message_id: &str) -> anyhow::Result<()> {
    let (db, _sync) = init_message_system(config)?;
    let msg = db
        .get_message(message_id)?
        .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

    archive_module(config, &msg)?;
    Ok(())
}

fn archive_module(config: &Config, msg: &Message) -> anyhow::Result<()> {
    // Update syft.pub.yaml by removing the results write rule
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let module_location = meta
        .get("module_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing module_location"))?;
    let root = resolve_syft_url_to_path(config, module_location)?;
    let perm_path = root.join("syft.pub.yaml");
    if perm_path.exists() {
        let storage = syftbox_storage(config)?;
        let bytes = storage.read_plaintext_file(&perm_path)?;
        let mut perms: SyftPermissions = serde_yaml::from_slice(&bytes)?;
        perms.rules.retain(|r| r.pattern != "results/**/*");
        let yaml = serde_yaml::to_string(&perms)?;
        storage.write_plaintext_file(&perm_path, yaml.as_bytes(), true)?;
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

fn sender_module_root(
    config: &Config,
    meta: &serde_json::Value,
) -> anyhow::Result<(PathBuf, String)> {
    let module_location = meta
        .get("module_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing module_location"))?;
    let path = resolve_syft_url_to_path(config, module_location)?;
    let folder = Path::new(module_location)
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

    // Check if this is a normal SyftBox root (has datasites/ folder inside)
    let real_data_dir = if data_dir.join("datasites").exists() {
        // Normal case: data_dir/datasites/email exists, so data_dir is the root
        data_dir
    } else if data_dir.components().any(|c| c.as_os_str() == "datasites")
        && data_dir.to_string_lossy().contains(&config.email)
    {
        // Edge case: SYFTBOX_DATA_DIR is pointing to datasite itself (no datasites/ folder inside)
        // Walk up to find the parent that doesn't contain "datasites"
        let mut parent = data_dir.clone();
        while parent.components().any(|c| c.as_os_str() == "datasites") {
            if let Some(p) = parent.parent() {
                parent = p.to_path_buf();
            } else {
                break;
            }
        }
        parent
    } else {
        // Fallback: treat as normal case
        data_dir
    };

    // Private is at data_dir root level
    Ok(real_data_dir
        .join("private")
        .join("app_data")
        .join("biovault")
        .join("submissions"))
}

fn copy_dir_recursive(config: &Config, src: &Path, dst: &Path) -> anyhow::Result<()> {
    let storage = syftbox_storage(config)?;
    let skip_syft_pub = |path: &Path| {
        path.file_name()
            .and_then(|n| n.to_str())
            .map(|n| n == "syft.pub.yaml")
            .unwrap_or(false)
    };

    // Determine recipient from destination path (e.g., "client1@sandbox.local")
    let data_dir = config.get_syftbox_data_dir()?;
    let datasites_dir = data_dir.join("datasites");
    let recipient = if dst.starts_with(&datasites_dir) {
        dst.strip_prefix(&datasites_dir)
            .ok()
            .and_then(|p| p.components().next())
            .and_then(|c| c.as_os_str().to_str())
            .map(|s| s.to_string())
    } else {
        None
    };

    for entry in walkdir::WalkDir::new(src)
        .into_iter()
        .filter_map(|e| e.ok())
    {
        let rel = entry
            .path()
            .strip_prefix(src)
            .with_context(|| format!("Failed to relativize {:?}", entry.path()))?;
        let out = dst.join(rel);

        if entry.file_type().is_dir() {
            ensure_dir_with_storage(&storage, &out)?;
            continue;
        }

        if skip_syft_pub(entry.path()) {
            continue;
        }

        if let Some(parent) = out.parent() {
            ensure_dir_with_storage(&storage, parent)?;
        }

        let data = if storage.contains(entry.path()) {
            // Use shadow-aware read to populate shadow folder when decrypting
            storage
                .read_with_shadow(entry.path())
                .with_context(|| format!("Failed to read {:?}", entry.path()))?
        } else {
            fs::read(entry.path()).with_context(|| format!("Failed to read {:?}", entry.path()))?
        };

        if storage.contains(&out) {
            // If we have a recipient, encrypt the file using shadow folder pattern
            if let Some(ref recip) = recipient {
                storage
                    .write_encrypted_with_shadow(&out, &data, vec![recip.clone()], None, true)
                    .with_context(|| format!("Failed to write encrypted {:?}", out))?;
            } else {
                // Fallback to plaintext if recipient unknown
                storage
                    .write_plaintext_file(&out, &data, true)
                    .with_context(|| format!("Failed to write {:?}", out))?;
            }
        } else {
            fs::write(&out, &data).with_context(|| format!("Failed to write {:?}", out))?;
        }
    }

    Ok(())
}

fn ensure_dir_with_storage(storage: &SyftBoxStorage, dir: &Path) -> anyhow::Result<()> {
    if storage.contains(dir) {
        storage.ensure_dir(dir)?;
    } else {
        fs::create_dir_all(dir)?;
    }
    Ok(())
}

fn list_dir_any(storage: Option<&SyftBoxStorage>, dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    if let Some(storage) = storage {
        if storage.contains(dir) {
            return storage.list_dir(dir);
        }
    }
    if !dir.exists() {
        return Ok(Vec::new());
    }
    Ok(fs::read_dir(dir)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect())
}

fn build_run_module_copy(config: &Config, msg: &Message) -> anyhow::Result<PathBuf> {
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, folder_name) = sender_module_root(config, meta)?;
    let dest_root = receiver_private_submissions_path(config)?;
    let dest_path = dest_root.join(&folder_name);

    // Always copy if key files are missing (handles incomplete/failed copies)
    let module_yaml = dest_path.join("module.yaml");
    let assets_dir = dest_path.join("assets");
    let assets_missing = !assets_dir.exists()
        || fs::read_dir(&assets_dir)
            .map(|mut entries| entries.next().is_none())
            .unwrap_or(true);
    if !module_yaml.exists() || assets_missing {
        copy_dir_recursive(config, &sender_root, &dest_path)?;
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
    reply.subject = Some(format!("Module {}", status));
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

fn reject_module(config: &Config, msg: &Message) -> anyhow::Result<()> {
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

fn review_module(config: &Config, msg: &Message) -> anyhow::Result<()> {
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

async fn run_module_test(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_module_copy(config, msg)?;
    let run_dir = prepare_run_directory(config, &dest, &msg.id)?;

    if let Some(invocation) = try_prepare_dynamic_run(config, &run_dir, true)? {
        run_dynamic::execute_dynamic(
            run_dir.to_string_lossy().as_ref(),
            invocation.args.clone(),
            false,
            false,
            Some(invocation.results_dir.clone()),
            run_dynamic::RunSettings::default(),
        )
        .await?;
        print_dynamic_results(&run_dir, &invocation.results_dir)?;
        return Ok(());
    }

    let source = prompt_participant_source("NA06985")?;
    let source_for_run = source.clone();
    run_execute(RunParams {
        module_folder: run_dir.to_string_lossy().to_string(),
        participant_source: source_for_run,
        test: true,
        download: true,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: None,
        results_dir: None,
        nextflow_args: vec![],
    })
    .await?;
    println!("{}", "Test run completed.".green());
    print_results_location_and_tree(config, &run_dir, &source, true)?;
    Ok(())
}

async fn run_module_real(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_module_copy(config, msg)?;
    let run_dir = prepare_run_directory(config, &dest, &msg.id)?;

    if let Some(invocation) = try_prepare_dynamic_run(config, &run_dir, false)? {
        run_dynamic::execute_dynamic(
            run_dir.to_string_lossy().as_ref(),
            invocation.args.clone(),
            false,
            false,
            Some(invocation.results_dir.clone()),
            run_dynamic::RunSettings::default(),
        )
        .await?;
        print_dynamic_results(&run_dir, &invocation.results_dir)?;

        let source_results = run_dir.join(&invocation.results_dir);
        let dest_results = dest.join(&invocation.results_dir);
        if source_results.exists() {
            copy_dir_recursive(config, &source_results, &dest_results)?;
        }
        return Ok(());
    }

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
        module_folder: run_dir.to_string_lossy().to_string(),
        participant_source: source_for_run,
        test: false,
        download: false,
        dry_run: false,
        with_docker: false,
        work_dir: None,
        resume: false,
        template: None,
        results_dir: None,
        nextflow_args: vec![],
    })
    .await?;
    println!("{}", "Real data run completed.".green());
    print_results_location_and_tree(config, &run_dir, &source, false)?;

    let results_dir = run_dir.join("results-real");
    let dest_results = dest.join("results-real");
    if results_dir.exists() {
        copy_dir_recursive(config, &results_dir, &dest_results)?;
    }
    Ok(())
}

async fn process_shell_module_message(
    config: &Config,
    msg: &Message,
    dest: &Path,
    test: bool,
    approve: bool,
) -> anyhow::Result<()> {
    let run_dir = prepare_run_directory(config, dest, &msg.id)?;
    let results_dir = if test { "results-test" } else { "results-real" };

    run_dynamic::execute_dynamic(
        run_dir.to_string_lossy().as_ref(),
        Vec::new(),
        false,
        false,
        Some(results_dir.to_string()),
        run_dynamic::RunSettings::default(),
    )
    .await?;

    let source_results = run_dir.join(results_dir);
    let dest_results = dest.join(results_dir);
    if source_results.exists() {
        copy_dir_recursive(config, &source_results, &dest_results)?;
    }

    println!("✓ Module processed successfully (shell template)");
    if dest_results.exists() {
        println!("Results saved to: {}", dest_results.display());
    }

    if approve {
        approve_module_non_interactive(config, msg, test).await?;
        println!("✓ Module approved and results sent to sender");
    }

    Ok(())
}

async fn process_dynamic_module_message(
    config: &Config,
    msg: &Message,
    dest: &Path,
    participant: Option<String>,
    test: bool,
    approve: bool,
) -> anyhow::Result<()> {
    let participant_hint = participant
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or("ALL");

    let run_dir = prepare_run_directory(config, dest, &msg.id)?;
    let invocation = prepare_dynamic_run_non_interactive(config, &run_dir, participant_hint, test)?;

    run_dynamic::execute_dynamic(
        run_dir.to_string_lossy().as_ref(),
        invocation.args.clone(),
        false,
        false,
        Some(invocation.results_dir.clone()),
        run_dynamic::RunSettings::default(),
    )
    .await?;

    let source_results = run_dir.join(&invocation.results_dir);
    let dest_results = dest.join(&invocation.results_dir);
    if source_results.exists() {
        copy_dir_recursive(config, &source_results, &dest_results)?;
    }

    let participant_label = Path::new(&invocation.results_dir)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(participant_hint);

    println!(
        "✓ Module processed successfully with participant: {}",
        participant_label
    );
    let results_path = dest.join(&invocation.results_dir);
    if results_path.exists() {
        println!("Results saved to: {}", results_path.display());
    }

    if approve {
        approve_module_non_interactive(config, msg, test).await?;
        println!("✓ Module approved and results sent to sender");
    }

    Ok(())
}

fn prepare_dynamic_run_non_interactive(
    config: &Config,
    module_dir: &Path,
    participant_hint: &str,
    test: bool,
) -> anyhow::Result<DynamicRunInvocation> {
    let spec_path = module_dir.join("module.yaml");
    let spec = ModuleSpec::load(&spec_path)?;

    if participant_hint.eq_ignore_ascii_case("ALL") || participant_hint.eq_ignore_ascii_case("AUTO")
    {
        ensure_catalog_ready_for_all(config)?;
    }

    let mut args = Vec::new();
    let mut samplesheet_value: Option<String> = None;
    let mut needs_data_dir = false;

    for input in &spec.inputs {
        match input.name.as_str() {
            "samplesheet" => {
                let path = resolve_samplesheet_path(config, module_dir, participant_hint)?;
                args.push("--set".to_string());
                args.push(format!("inputs.{}={}", input.name, path));
                samplesheet_value = Some(path);
            }
            "data_dir" => {
                needs_data_dir = true;
            }
            other => {
                let is_optional = input.raw_type.trim_end().ends_with('?');
                if !is_optional {
                    return Err(anyhow::anyhow!(
                        "No automatic binding available for required input '{}'",
                        other
                    ));
                }
            }
        }
    }

    let samplesheet_path = samplesheet_value.ok_or_else(|| {
        anyhow::anyhow!("Dynamic module requires a 'samplesheet' input but none could be inferred")
    })?;

    if needs_data_dir {
        let data_dir = infer_data_dir_path(config, &samplesheet_path)?;
        args.push("--set".to_string());
        args.push(format!("inputs.data_dir={}", data_dir.to_string_lossy()));
    }

    let participant_label = derive_participant_label(participant_hint, &samplesheet_path);
    let results_dir = if test {
        format!("results-test/{}", participant_label)
    } else {
        format!("results-real/{}", participant_label)
    };

    Ok(DynamicRunInvocation { args, results_dir })
}

fn resolve_samplesheet_path(
    config: &Config,
    module_dir: &Path,
    participant_hint: &str,
) -> anyhow::Result<String> {
    if participant_hint.eq_ignore_ascii_case("ALL") || participant_hint.eq_ignore_ascii_case("AUTO")
    {
        let generated = auto_generate_samplesheet(config, module_dir)?;
        return Ok(canonicalize_string(PathBuf::from(generated)));
    }

    if participant_hint.starts_with("syft://") {
        let path = resolve_syft_url_to_path(config, participant_hint)?;
        if !path.exists() {
            return Err(anyhow::anyhow!(
                "Resolved samplesheet not found at {}",
                path.display()
            ));
        }
        return Ok(canonicalize_string(path));
    }

    let direct = PathBuf::from(participant_hint);
    let candidate = if direct.is_absolute() && direct.exists() {
        direct
    } else {
        let module_candidate = module_dir.join(&direct);
        if module_candidate.exists() {
            module_candidate
        } else {
            let data_root = config.get_syftbox_data_dir()?;
            let data_candidate = data_root.join(&direct);
            if data_candidate.exists() {
                data_candidate
            } else {
                direct
            }
        }
    };

    if !candidate.exists() {
        return Err(anyhow::anyhow!(
            "Samplesheet not found: {}",
            participant_hint
        ));
    }

    Ok(canonicalize_string(candidate))
}

fn infer_data_dir_path(config: &Config, samplesheet_path: &str) -> anyhow::Result<PathBuf> {
    let sheet_path = PathBuf::from(samplesheet_path);
    if let Some(inferred) = infer_data_dir_from_samplesheet(&sheet_path) {
        return Ok(canonicalize_pathbuf(inferred));
    }

    let mut candidates: Vec<PathBuf> = Vec::new();
    if let Some(parent) = sheet_path.parent() {
        candidates.push(parent.to_path_buf());
    }
    let data_root = config.get_syftbox_data_dir()?;
    let fixtures = data_root.join("genotype_fixtures");
    candidates.push(fixtures);
    candidates.push(data_root);

    for candidate in candidates {
        if candidate.exists() {
            return Ok(canonicalize_pathbuf(candidate));
        }
    }

    Ok(sheet_path
        .parent()
        .map(|p| canonicalize_pathbuf(p.to_path_buf()))
        .unwrap_or_else(|| PathBuf::from(samplesheet_path)))
}

fn infer_data_dir_from_samplesheet(path: &Path) -> Option<PathBuf> {
    let mut reader = csv::Reader::from_path(path).ok()?;
    let headers = reader.headers().ok()?.clone();
    let mut column_index: Option<usize> = None;
    for candidate in ["genotype_file_path", "genotype_file", "file_path", "path"] {
        if let Some(idx) = headers.iter().position(|h| h == candidate) {
            column_index = Some(idx);
            break;
        }
    }
    let col_idx = column_index?;

    for result in reader.records() {
        let record = result.ok()?;
        let value = record.get(col_idx)?.trim();
        if value.is_empty() {
            continue;
        }
        let candidate = PathBuf::from(value);
        if candidate.is_absolute() {
            return candidate.parent().map(|p| p.to_path_buf());
        }
    }
    None
}

fn derive_participant_label(participant_hint: &str, samplesheet_path: &str) -> String {
    if participant_hint.eq_ignore_ascii_case("ALL") || participant_hint.eq_ignore_ascii_case("AUTO")
    {
        return "ALL".to_string();
    }
    if let Some(id) = extract_participant_id(participant_hint) {
        return id;
    }
    Path::new(samplesheet_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("run")
        .to_string()
}

fn canonicalize_string(path: PathBuf) -> String {
    canonicalize_pathbuf(path).to_string_lossy().to_string()
}

fn canonicalize_pathbuf(path: PathBuf) -> PathBuf {
    path.canonicalize().unwrap_or(path)
}

fn ensure_catalog_ready_for_all(_config: &Config) -> anyhow::Result<()> {
    let db = BioVaultDb::new()?;
    if !data::list_files(&db, None, None, false, Some(1))?.is_empty() {
        return Ok(());
    }
    Err(anyhow::anyhow!(
        "No cataloged files found. Import genotype data first with 'bv files import'."
    ))
}

/// Non-interactive version of approve_module for automated testing
async fn approve_module_non_interactive(
    config: &Config,
    msg: &Message,
    test: bool,
) -> anyhow::Result<()> {
    let dest = build_run_module_copy(config, msg)?;
    let storage = syftbox_storage(config)?;
    let results_dir = if test {
        // In test mode, use test results for approval
        dest.join("results-test")
    } else {
        dest.join("results-real")
    };

    // Check if results exist using filesystem directly (more robust than storage.list_dir)
    let needs_run = !results_dir.exists() || {
        std::fs::read_dir(&results_dir)
            .map(|mut entries| entries.next().is_none())
            .unwrap_or(true)
    };

    if needs_run && !test {
        // Only run real data if not in test mode
        println!("No results found. Running on real data before approval...");
        run_module_real(config, msg).await?;
    } else if needs_run && test {
        return Err(anyhow::anyhow!(
            "No test results found. Run with --test first before approving."
        ));
    }

    // Release results to sender shared location
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, _folder) = sender_module_root(config, meta)?;
    let sender_results = sender_root.join("results");
    copy_dir_recursive(config, &results_dir, &sender_results)?;

    // Get module details for the message
    let module_location = meta
        .get("module_location")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown");

    // Check if results exist and get basic info
    let results_info = if sender_results.exists() {
        let mut info = Vec::new();
        if let Ok(entries) = storage.list_dir(&sender_results) {
            for entry in entries {
                if entry.is_dir() {
                    if let Some(name) = entry.file_name() {
                        info.push(format!("{}/", name.to_string_lossy()));
                    }
                } else if let Some(name) = entry.file_name() {
                    info.push(name.to_string_lossy().to_string());
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
        "Your module has been approved.\n\nModule location: {}\nResults location: {}\n\n{}",
        module_location, relative_results, results_info
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

async fn approve_module(config: &Config, msg: &Message) -> anyhow::Result<()> {
    let dest = build_run_module_copy(config, msg)?;
    let results_dir = dest.join("results-real");
    let storage = syftbox_storage(config)?;
    let needs_run = !results_dir.exists()
        || storage
            .list_dir(&results_dir)
            .map(|entries| entries.is_empty())
            .unwrap_or(true);
    if needs_run {
        println!("No results found. Running on real data before approval...");
        run_module_real(config, msg).await?;
    }

    // Release results to sender shared location
    let meta = msg
        .metadata
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing metadata"))?;
    let (sender_root, _folder) = sender_module_root(config, meta)?;
    let sender_results = sender_root.join("results");
    copy_dir_recursive(config, &results_dir, &sender_results)?;

    // Get module details for the message
    let module_location = meta
        .get("module_location")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown");

    // Check if results exist and get basic info
    let results_info = if sender_results.exists() {
        let mut info = Vec::new();
        if let Ok(entries) = storage.list_dir(&sender_results) {
            for entry in entries {
                if entry.is_dir() {
                    if let Some(name) = entry.file_name() {
                        info.push(format!("{}/", name.to_string_lossy()));
                    }
                } else if let Some(name) = entry.file_name() {
                    info.push(name.to_string_lossy().to_string());
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
        "Your module has been approved.\n\nModule location: {}\nResults location: {}\n\n{}",
        module_location, relative_results, results_info
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

struct DynamicRunInvocation {
    args: Vec<String>,
    results_dir: String,
}

fn try_prepare_dynamic_run(
    config: &Config,
    module_dir: &Path,
    is_test: bool,
) -> anyhow::Result<Option<DynamicRunInvocation>> {
    let spec_path = module_dir.join("module.yaml");
    if !spec_path.exists() {
        return Ok(None);
    }

    let spec = match ModuleSpec::load(&spec_path) {
        Ok(spec) => spec,
        Err(err) => {
            println!(
                "⚠ Failed to load module spec at {}: {}",
                spec_path.display(),
                err
            );
            return Ok(None);
        }
    };

    if spec.inputs.is_empty() {
        return Ok(None);
    }

    println!(
        "Module '{}' requires {} input(s).",
        spec.name,
        spec.inputs.len()
    );

    let mut args = Vec::new();
    for input in &spec.inputs {
        let mut prompt = format!("Enter value for input '{}'", input.name);
        if let Some(desc) = &input.description {
            if !desc.trim().is_empty() {
                prompt.push_str(&format!(" ({})", desc.trim()));
            }
        }
        prompt.push(':');

        let input_value: String = Input::new().with_prompt(prompt).interact_text()?;
        let trimmed = input_value.trim();

        let resolved = if trimmed.is_empty() {
            if input.name == "samplesheet" {
                auto_generate_samplesheet(config, module_dir)?
            } else {
                trimmed.to_string()
            }
        } else if input.name == "samplesheet"
            && matches!(trimmed.to_ascii_uppercase().as_str(), "ALL" | "AUTO")
        {
            auto_generate_samplesheet(config, module_dir)?
        } else {
            trimmed.to_string()
        };

        if resolved.is_empty() {
            println!(
                "⚠ Input '{}' left empty; skipping dynamic run fallback.",
                input.name
            );
            return Ok(None);
        }

        println!("  • {} = {}", input.name, resolved);
        args.push("--set".to_string());
        args.push(format!("inputs.{}={}", input.name, resolved));
    }

    let results_dir = if is_test {
        "results-test".to_string()
    } else {
        "results-real".to_string()
    };

    let preview = format_bv_run_command(module_dir, &results_dir, &args);
    println!("Nextflow command:\n  {}", preview);

    Ok(Some(DynamicRunInvocation { args, results_dir }))
}

fn format_bv_run_command(module_dir: &Path, results_dir: &str, args: &[String]) -> String {
    let module_path = module_dir.to_string_lossy().into_owned();
    let mut tokens = vec![
        "bv".to_string(),
        "run".to_string(),
        shell_escape(&module_path),
        "--results-dir".to_string(),
        shell_escape(results_dir),
    ];
    for arg in args {
        tokens.push(shell_escape(arg));
    }
    tokens.join(" ")
}

fn shell_escape(arg: &str) -> String {
    if arg
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || "-._/@=+:".contains(c))
    {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', "'\\''"))
    }
}

fn print_dynamic_results(module_root: &Path, results_dir: &str) -> anyhow::Result<()> {
    let results_path = module_root.join(results_dir);
    println!("Results location: {}", results_path.display());
    if results_path.exists() {
        println!("Results tree:");
        print_dir_tree(None, &results_path, 3)?;
    } else {
        println!("Results folder not found yet at {}", results_path.display());
    }
    Ok(())
}

fn prepare_run_directory(
    config: &Config,
    source: &Path,
    message_id: &str,
) -> anyhow::Result<PathBuf> {
    let base = config.get_biovault_dir()?;
    let run_root = base.join("runs");
    fs::create_dir_all(&run_root)?;
    let run_dir = run_root.join(message_id);
    if run_dir.exists() {
        fs::remove_dir_all(&run_dir)?;
    }
    copy_dir_recursive(config, source, &run_dir)?;
    Ok(run_dir)
}

fn auto_generate_samplesheet(_config: &Config, module_dir: &Path) -> anyhow::Result<String> {
    let db = BioVaultDb::new()?;
    let mut stmt = db.conn.prepare(
        "SELECT f.file_path, p.participant_id
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
         ORDER BY f.file_path",
    )?;
    let mut rows = stmt.query([])?;
    let mut entries: Vec<(String, String)> = Vec::new();
    while let Some(row) = rows.next()? {
        let path: String = row.get(0)?;
        let participant: Option<String> = row.get(1)?;
        entries.push((participant.unwrap_or_default(), path));
    }
    if entries.is_empty() {
        anyhow::bail!(
            "No cataloged files found. Import genotype data first with 'bv files import'."
        );
    }

    let inputs_dir = module_dir.join("inputs");
    fs::create_dir_all(&inputs_dir)?;
    let sheet_path = inputs_dir.join("auto_samplesheet.csv");
    let mut writer = Writer::from_path(&sheet_path)?;
    writer.write_record(["participant_id", "genotype_file_path"])?;
    for (participant, file_path) in entries {
        let pid = if participant.trim().is_empty() {
            Path::new(&file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string()
        } else {
            participant
        };
        writer.write_record([pid, file_path])?;
    }
    writer.flush()?;
    println!("  ✓ Generated samplesheet at {}", sheet_path.display());
    Ok(sheet_path.to_string_lossy().to_string())
}

fn print_results_location_and_tree(
    config: &Config,
    module_root: &Path,
    participant_source: &str,
    is_test: bool,
) -> anyhow::Result<()> {
    let id = extract_participant_id(participant_source).unwrap_or_else(|| "ALL".to_string());
    let base = if is_test {
        "results-test"
    } else {
        "results-real"
    };
    let storage = syftbox_storage(config).ok();
    let results_dir = module_root.join(base).join(&id);
    let base_dir = module_root.join(base);
    println!("Results location: {}", results_dir.display());
    if results_dir.exists() {
        println!("Results tree:");
        print_dir_tree(storage.as_ref(), &results_dir, 3)?;
    } else if base_dir.exists() {
        println!(
            "Results folder for participant not found; showing {} instead",
            base_dir.display()
        );
        print_dir_tree(storage.as_ref(), &base_dir, 3)?;
    } else {
        println!("Results folder not found yet at {}", results_dir.display());
    }
    Ok(())
}

fn print_dir_tree(
    storage: Option<&SyftBoxStorage>,
    root: &Path,
    max_depth: usize,
) -> anyhow::Result<()> {
    fn walk(
        storage: Option<&SyftBoxStorage>,
        dir: &Path,
        depth: usize,
        max_depth: usize,
    ) -> anyhow::Result<()> {
        if depth > max_depth {
            return Ok(());
        }
        let mut entries = list_dir_any(storage, dir)?;
        entries.sort();
        for path in entries {
            let indent = "  ".repeat(depth.saturating_sub(0));
            let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if path.is_dir() {
                println!("{}{}/", indent, name);
                walk(storage, &path, depth + 1, max_depth)?;
            } else {
                println!("{}{}", indent, name);
            }
        }
        Ok(())
    }
    println!("{}/", root.display());
    walk(storage, root, 1, max_depth)
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
pub fn sync_messages(config: &Config, no_cleanup: bool) -> Result<()> {
    let (_, sync) = init_message_system(config)?;

    println!("🔄 Syncing messages...");
    if no_cleanup {
        println!("🐛 Debug mode: keeping consumed messages (--no-cleanup)");
    }
    sync.sync(no_cleanup)?;

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

/// Verify a module message using embedded metadata
/// Returns (is_ok, extra_note)
fn verify_module_from_metadata(
    config: &Config,
    metadata: &serde_json::Value,
) -> anyhow::Result<(bool, String)> {
    // Extract the embedded module.yaml
    let module_val = metadata
        .get("module")
        .ok_or_else(|| anyhow::anyhow!("missing module metadata"))?;

    let module: ModuleYaml = serde_json::from_value(module_val.clone())
        .map_err(|e| anyhow::anyhow!("invalid module metadata: {}", e))?;

    // Extract module location syft URL
    let module_location = metadata
        .get("module_location")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("missing module_location syft URL"))?;

    let root = resolve_syft_url_to_path(config, module_location)?;
    if !root.exists() {
        return Ok((false, format!(" (missing path: {})", root.display())));
    }

    let storage = syftbox_storage(config)?;

    // Need b3_hashes to verify
    let Some(expected_hashes) = module.b3_hashes.clone() else {
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
        let bytes = storage
            .read_plaintext_file(&file_path)
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
    use crate::syftbox::SyftBoxApp;
    use serde_json::json;
    use tempfile::TempDir;

    fn create_test_config() -> Config {
        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
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
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
        };

        let app = SyftBoxApp::new(tmp.path(), &cfg.email, "biovault").unwrap();
        let root = app.app_data_dir.join("submissions").join("proj1");
        app.storage.ensure_dir(&root.join("a/b")).unwrap();
        app.storage
            .write_plaintext_file(&root.join("a/b/file.txt"), b"hi", true)
            .unwrap();
        app.storage
            .write_plaintext_file(&root.join("a/syft.pub.yaml"), b"rules:", true)
            .unwrap();

        // sender_module_root from metadata
        let loc = format!(
            "syft://u@example.com/app_data/biovault/submissions/{}",
            "proj1"
        );
        let meta = json!({
            "module_location": loc
        });
        let (sender_root, folder) = sender_module_root(&cfg, &meta)?;
        assert_eq!(folder, "proj1");
        assert!(sender_root.ends_with("proj1"));

        // receiver path (should be at root level, not inside datasite)
        let recv = receiver_private_submissions_path(&cfg)?;
        assert!(recv.ends_with("private/app_data/biovault/submissions"));
        assert!(!recv
            .to_string_lossy()
            .contains("datasites/u@example.com/private"));

        // Copy behaviour: skip syft.pub.yaml and recreate tree
        let dest = recv.join(&folder);
        copy_dir_recursive(&cfg, &sender_root, &dest)?;
        assert!(dest.join("a/b/file.txt").exists());
        assert!(!dest.join("a/syft.pub.yaml").exists());

        Ok(())
    }

    #[test]
    fn build_run_module_copy_creates_dest_once() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = Config {
            email: "u@example.com".into(),
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

        // Make sender tree and one file
        let app = SyftBoxApp::new(tmp.path(), &cfg.email, "biovault").unwrap();
        let proj = app.app_data_dir.join("submissions").join("projX");
        app.storage.ensure_dir(&proj.join("dir")).unwrap();
        app.storage
            .write_plaintext_file(&proj.join("dir/f.txt"), b"x", true)
            .unwrap();
        let meta =
            json!({"module_location": "syft://u@example.com/app_data/biovault/submissions/projX"});
        let mut msg = Message::new("u@example.com".into(), "v@example.com".into(), "b".into());
        msg.metadata = Some(meta);

        let dest1 = build_run_module_copy(&cfg, &msg)?;
        assert!(dest1.join("dir/f.txt").exists());
        // Call again; should not error and keep same path
        let dest2 = build_run_module_copy(&cfg, &msg)?;
        assert_eq!(dest1, dest2);
        Ok(())
    }

    #[test]
    fn verify_module_from_metadata_ok_and_fail() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = Config {
            email: "u@example.com".into(),
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

        // Build module root with two files and compute blake3
        let app = SyftBoxApp::new(tmp.path(), &cfg.email, "biovault").unwrap();
        let root = app.app_data_dir.join("submissions").join("projZ");
        app.storage.ensure_dir(&root).unwrap();
        let f1 = root.join("a.txt");
        let f2 = root.join("b.txt");
        app.storage.write_plaintext_file(&f1, b"A", true).unwrap();
        app.storage.write_plaintext_file(&f2, b"B", true).unwrap();
        let h1 = blake3::hash(&app.storage.read_plaintext_file(&f1).unwrap())
            .to_hex()
            .to_string();
        let h2 = blake3::hash(&app.storage.read_plaintext_file(&f2).unwrap())
            .to_hex()
            .to_string();

        // Construct embedded module.yaml equivalent in metadata
        let module = serde_json::to_value(ModuleYaml {
            name: "n".into(),
            version: "0.1.0".into(),
            author: "a".into(),
            datasites: None,
            participants: None,
            workflow: "w".into(),
            template: None,
            assets: None,
            parameters: None,
            inputs: None,
            outputs: None,
            b3_hashes: Some(std::collections::HashMap::from([
                ("a.txt".to_string(), h1.clone()),
                ("b.txt".to_string(), h2.clone()),
            ])),
            extra: std::collections::HashMap::new(),
        })?;
        let meta_ok = json!({
            "module": module,
            "module_location": "syft://u@example.com/app_data/biovault/submissions/projZ"
        });
        let (ok, note) = verify_module_from_metadata(&cfg, &meta_ok)?;
        assert!(ok, "expected OK, got note: {}", note);

        // Now break one file and expect failure with note
        app.storage
            .write_plaintext_file(&f2, b"BROKEN", true)
            .unwrap();
        let (ok2, note2) = verify_module_from_metadata(&cfg, &meta_ok)?;
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
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
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
        list_messages(&cfg, false, false, false, false)?;
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
            syftbox_credentials: None,
            agent_bridge_enabled: None,
            agent_bridge_port: None,
            agent_bridge_http_port: None,
            agent_bridge_token: None,
            agent_bridge_blocklist: None,
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

        read_message(&cfg, &m.id, false).await?;
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
        super::list_messages(&cfg, false, false, false, false)?;
        super::list_messages(&cfg, true, false, false, false)?;
        Ok(())
    }

    #[test]
    fn test_view_thread_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = create_test_config();

        // View non-existent thread should not error
        super::view_thread(&cfg, "nonexistent-thread")?;
        Ok(())
    }

    #[test]
    fn test_sync_messages() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let cfg = create_test_config();

        // Should not error even without messages
        super::sync_messages(&cfg, false)?;
        Ok(())
    }

    #[test]
    fn test_resolve_syft_url_to_path() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = create_test_config();

        let url = "syft://user@example.com/path/to/file";
        let path = super::resolve_syft_url_to_path(&cfg, url)?;

        assert!(path.to_string_lossy().contains("datasites"));
        assert!(path.to_string_lossy().contains("user@example.com"));
        assert!(path.to_string_lossy().contains("path/to/file"));
        Ok(())
    }

    #[test]
    fn test_get_message_db_path() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_biovault_home(tmp.path());
        let cfg = create_test_config();

        let path = get_message_db_path(&cfg)?;
        assert!(path.to_string_lossy().contains("messages.db"));
        assert!(path.parent().unwrap().exists());
        Ok(())
    }

    #[test]
    fn test_cleanup_locks() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_biovault_home(tmp.path());
        let cfg = create_test_config();

        let result = cleanup_locks(&cfg, false);
        assert!(result.is_ok());
        Ok(())
    }

    #[test]
    fn test_expand_env_vars_in_text_multiple() -> Result<()> {
        std::env::set_var("SYFTBOX_DATA_DIR", "/data");
        let text = "$SYFTBOX_DATA_DIR/a and $SYFTBOX_DATA_DIR/b";
        let result = expand_env_vars_in_text(text)?;
        assert_eq!(result, "/data/a and /data/b");
        std::env::remove_var("SYFTBOX_DATA_DIR");
        Ok(())
    }

    #[test]
    fn test_receiver_private_submissions_path() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = create_test_config();

        let path = receiver_private_submissions_path(&cfg)?;
        assert!(path.to_string_lossy().contains("private"));
        assert!(path.to_string_lossy().contains("submissions"));
        Ok(())
    }

    #[test]
    fn test_copy_dir_recursive_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        let src = tmp.path().join("src");
        let dst = tmp.path().join("dst");
        std::fs::create_dir(&src)?;

        crate::config::set_test_syftbox_data_dir(tmp.path());
        let cfg = create_test_config();
        copy_dir_recursive(&cfg, &src, &dst)?;
        crate::config::clear_test_syftbox_data_dir();
        assert!(dst.exists());
        Ok(())
    }

    #[test]
    fn test_normalize_participant_source_for_real_url() -> Result<()> {
        let tmp = TempDir::new()?;
        crate::config::set_test_biovault_home(tmp.path());

        let url = "syft://user@example.com/file#id";
        let result = normalize_participant_source_for_real(url)?;
        assert_eq!(result, url);
        Ok(())
    }
}
