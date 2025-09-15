use crate::config::Config;
use crate::messages::{MessageDb, MessageType};
use anyhow::Result;
use dialoguer::{theme::ColorfulTheme, Select};
use std::io::{self, Write};

/// Filter options for listing messages
pub struct ListFilters {
    pub sent: bool,
    pub all: bool,
    pub unread: bool,
    pub projects: bool,
    pub message_type: Option<String>,
    pub from: Option<String>,
    pub search: Option<String>,
}

/// Display messages in the inbox with optional filters
pub fn list(config: &Config, filters: ListFilters) -> Result<()> {
    let db_path = super::messages::get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    // Apply filters
    let messages = if filters.all {
        db.list_messages(Some(100))?
    } else if filters.sent {
        db.list_sent_messages(Some(100))?
    } else if filters.unread {
        db.list_unread_messages()?
    } else if filters.projects {
        db.list_messages_by_type("project", Some(100))?
    } else if let Some(ref search_term) = filters.search {
        db.search_messages(search_term, Some(100))?
    } else if let Some(ref msg_type) = filters.message_type {
        db.list_messages_by_type(msg_type, Some(100))?
    } else {
        // Default to inbox
        db.list_inbox_messages(Some(100))?
    };

    // Filter by sender if specified
    let messages = if let Some(ref sender) = filters.from {
        messages
            .into_iter()
            .filter(|m| m.from.contains(sender))
            .collect()
    } else {
        messages
    };

    if messages.is_empty() {
        println!("No messages found");
        return Ok(());
    }

    // Display messages
    println!("\n📬 Messages ({} total):", messages.len());
    println!("═══════════════════════════════════════");

    for (i, msg) in messages.iter().enumerate() {
        let status_icon = match msg.status {
            crate::messages::MessageStatus::Draft => "📝",
            crate::messages::MessageStatus::Sent => "📤",
            crate::messages::MessageStatus::Received => "📥",
            crate::messages::MessageStatus::Read => "👁️",
            crate::messages::MessageStatus::Deleted => "🗑️",
            crate::messages::MessageStatus::Archived => "📁",
        };

        let type_icon = match &msg.message_type {
            MessageType::Text => "✉️",
            MessageType::Project { .. } => "📦",
            MessageType::Request { .. } => "🔔",
        };

        println!(
            "\n{}. {} {} [{}]",
            i + 1,
            status_icon,
            type_icon,
            &msg.id[..8]
        );
        println!("   From: {}", msg.from);
        println!("   To: {}", msg.to);

        if let Some(ref subject) = msg.subject {
            if !subject.is_empty() {
                println!("   Subject: {}", subject);
            }
        }

        let local_time = msg.created_at.with_timezone(&chrono::Local);
        println!("   Date: {}", local_time.format("%Y-%m-%d %H:%M"));

        // Show preview of body
        let preview_len = 80;
        let preview = if msg.body.len() > preview_len {
            format!("{}...", &msg.body[..preview_len])
        } else {
            msg.body.clone()
        };
        println!("   {}", preview);
    }

    println!("\n─────────────────────────────────────");
    println!("Use 'bv inbox -i' for interactive mode");

    Ok(())
}

/// Interactive mode for inbox
pub fn interactive(config: &Config, initial_view: Option<String>) -> Result<()> {
    let db_path = super::messages::get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    // Sync messages first
    let sync = super::messages::init_message_system(config)?.1;
    let _ = sync.sync_quiet();

    let mut current_view = initial_view.unwrap_or_else(|| "inbox".to_string());

    loop {
        // Clear and reset cursor position for each loop iteration
        print!("\x1B[2J\x1B[1;1H");
        io::stdout().flush()?;

        // Get messages based on current view
        let messages = match current_view.as_str() {
            "inbox" => db.list_inbox_messages(Some(50))?,
            "sent" => db.list_sent_messages(Some(50))?,
            "all" => db.list_messages(Some(50))?,
            "unread" => db.list_unread_messages()?,
            "projects" => db.list_messages_by_type("project", Some(50))?,
            _ => db.list_inbox_messages(Some(50))?,
        };

        // Print header (ASCII only to avoid emoji width issues)
        println!("======================================================");
        println!(
            "BioVault Inbox - {} ({} messages)",
            current_view.to_uppercase(),
            messages.len()
        );
        println!("======================================================");
        println!();

        // Build display options
        let mut display_options = Vec::new();

        // Add messages if any
        if messages.is_empty() {
            display_options.push("(No messages in this view)".to_string());
        } else {
            for msg in &messages {
                // Build concise ASCII-only line to prevent width mis-calculation and wrapping
                let status = match msg.status {
                    crate::messages::MessageStatus::Draft => "DRAFT",
                    crate::messages::MessageStatus::Sent => "SENT",
                    crate::messages::MessageStatus::Received => "RECV",
                    crate::messages::MessageStatus::Read => "READ",
                    crate::messages::MessageStatus::Deleted => "DEL",
                    crate::messages::MessageStatus::Archived => "ARCH",
                };

                let from_or_to = if current_view == "sent" {
                    format!("To:{}", msg.to)
                } else {
                    format!("From:{}", msg.from)
                };

                let subject = msg.subject.as_deref().unwrap_or("(No Subject)");
                let mut option = format!("[{status}] {from_or_to} - {subject}");
                // Truncate aggressively to avoid line wrapping that causes scroll glitches
                const MAX_ITEM_LEN: usize = 80;
                if option.len() > MAX_ITEM_LEN {
                    option.truncate(MAX_ITEM_LEN);
                }
                display_options.push(option);
            }
        }

        // Always add menu actions at the bottom
        display_options.push("------------------------------------------------------".to_string());
        display_options.push("Change View".to_string());
        display_options.push("Sync Messages".to_string());
        display_options.push("Quit".to_string());

        // Limit display height to prevent scrolling issues
        // Keep it small to avoid the scrolling problem
        let selection = Select::with_theme(&ColorfulTheme::default())
            // Keep prompt ASCII and short to avoid wrapping
            .with_prompt("Select (Arrows/Enter, Esc to quit)")
            .default(0)
            .items(&display_options)
            .max_length(6) // Very conservative to prevent scrolling
            .interact_opt()?;

        match selection {
            None => {
                // User pressed Esc/Ctrl+C
                break;
            }
            Some(idx) => {
                // Calculate where we are in the menu
                let separator_idx = if messages.is_empty() {
                    1
                } else {
                    messages.len()
                };
                let change_view_idx = separator_idx + 1;
                let sync_idx = separator_idx + 2;
                let quit_idx = separator_idx + 3;

                if messages.is_empty() && idx == 0 {
                    // "(No messages in this view)" selected - do nothing
                    continue;
                } else if !messages.is_empty() && idx < messages.len() {
                    // Message selected - show actions
                    let msg = &messages[idx];
                    message_actions(config, &db, msg)?;
                } else if idx == separator_idx {
                    // Separator line - do nothing
                    continue;
                } else if idx == change_view_idx {
                    // Change view
                    let views = ["inbox", "sent", "all", "unread", "projects"];
                    // ASCII-only names to avoid width calc issues
                    let view_names = ["Inbox", "Sent", "All Messages", "Unread", "Projects"];

                    println!("\nSelect view:");
                    let view_selection = Select::with_theme(&ColorfulTheme::default())
                        .with_prompt("Choose a view")
                        .default(0)
                        .items(&view_names)
                        .interact_opt()?;

                    if let Some(selection) = view_selection {
                        current_view = views[selection].to_string();
                    }
                } else if idx == sync_idx {
                    // Sync
                    println!("\nSyncing messages...");
                    sync.sync()?;
                    println!("Sync complete! Press Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                } else if idx == quit_idx {
                    // Quit
                    break;
                }
            }
        }
    }

    Ok(())
}

/// Handle actions for a selected message
/// Returns false if user selected "Back", true otherwise
fn message_actions(
    config: &Config,
    db: &MessageDb,
    msg: &crate::messages::Message,
) -> Result<bool> {
    let actions = vec![
        "Read",
        "Reply",
        "Delete",
        "Mark as Read/Unread",
        "Back to Messages",
    ];

    let selection = Select::with_theme(&ColorfulTheme::default())
        .with_prompt("Select action")
        .default(0)
        .items(&actions)
        .interact_opt()?;

    match selection {
        Some(0) => {
            // Read
            super::messages::read_message(config, &msg.id)?;
            println!("\nPress Enter to continue...");
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            Ok(false) // Continue to message list
        }
        Some(1) => {
            // Reply
            println!("\nEnter your reply (press Enter when done):");
            print!("> ");
            io::stdout().flush()?;
            let mut reply_body = String::new();
            io::stdin().read_line(&mut reply_body)?;

            if !reply_body.trim().is_empty() {
                super::messages::reply_message(config, &msg.id, reply_body.trim())?;
                println!("✅ Reply sent!");
            } else {
                println!("Reply cancelled (empty message)");
            }

            println!("Press Enter to continue...");
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            Ok(false)
        }
        Some(2) => {
            // Delete
            super::messages::delete_message(config, &msg.id)?;
            println!("🗑️ Message deleted. Press Enter to continue...");
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            Ok(false)
        }
        Some(3) => {
            // Mark as read/unread
            if msg.status == crate::messages::MessageStatus::Received {
                db.mark_as_read(&msg.id)?;
                println!("✅ Marked as read");
            } else if msg.status == crate::messages::MessageStatus::Read {
                println!("ℹ️ Message is already read (mark as unread not implemented yet)");
            } else {
                println!("ℹ️ Cannot change read status for this message type");
            }
            println!("Press Enter to continue...");
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            Ok(false)
        }
        Some(4) | None => {
            // Back to messages or ESC pressed
            Ok(false)
        }
        _ => Ok(false),
    }
}
