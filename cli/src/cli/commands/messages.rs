use crate::config::Config;
use crate::messages::{Message, MessageDb, MessageSync};
use anyhow::Result;
use std::path::PathBuf;

const MESSAGE_ENDPOINT: &str = "/message";

/// Get the path to the message database
fn get_message_db_path(config: &Config) -> Result<PathBuf> {
    let biovault_dir = config.get_biovault_dir()?;
    let db_path = biovault_dir.join("data").join("messages.db");

    // Ensure the data directory exists
    if let Some(parent) = db_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    Ok(db_path)
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
pub fn list_messages(config: &Config, unread_only: bool) -> Result<()> {
    let (db, sync) = init_message_system(config)?;

    // Quietly sync to get latest messages and show notification if new
    let (_new_msg_ids, count) = sync.sync_quiet()?;
    if count > 0 {
        println!("🆕 {} new message(s) received", count);
    }

    let messages = if unread_only {
        db.list_unread_messages()?
    } else {
        db.list_messages(Some(50))?
    };

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

        if msg.parent_id.is_some() {
            println!("  ↩️  Reply to: {}", msg.parent_id.as_ref().unwrap());
        }
    }

    Ok(())
}

/// Read a specific message
pub fn read_message(config: &Config, message_id: &str) -> Result<()> {
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
    println!("{}", msg.body);

    Ok(())
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

        println!("{}", msg.body);
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn create_test_config() -> Config {
        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
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
}
