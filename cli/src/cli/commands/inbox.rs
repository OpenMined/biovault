use crate::config::Config;
use crate::messages::{MessageDb, MessageType};
use anyhow::Result;
use dialoguer::{theme::ColorfulTheme, Select};
use std::io::Read;
use std::time::{Duration, Instant};

/// Expand environment variables in text (specifically $SYFTBOX_DATA_DIR)
fn expand_env_vars_in_text(text: &str) -> Result<String> {
    let mut result = text.to_string();

    // Expand $SYFTBOX_DATA_DIR
    if let Ok(data_dir) = std::env::var("SYFTBOX_DATA_DIR") {
        result = result.replace("$SYFTBOX_DATA_DIR", &data_dir);
    }

    Ok(result)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Key {
    Up,
    Down,
    Enter,
    Esc,
    Char(char),
}

fn enable_raw_mode_cmd() -> Result<()> {
    // Use `stty` to disable canonical mode and echo
    std::process::Command::new("stty")
        .arg("-icanon")
        .arg("-echo")
        .arg("min")
        .arg("1")
        .arg("time")
        .arg("0")
        .status()?;
    Ok(())
}

fn disable_raw_mode_cmd() -> Result<()> {
    std::process::Command::new("stty").arg("sane").status()?;
    Ok(())
}

fn read_key() -> Result<Key> {
    let mut stdin = io::stdin();
    let mut buf = [0u8; 1];
    stdin.read_exact(&mut buf)?;
    match buf[0] {
        b'\n' | b'\r' => Ok(Key::Enter),
        0x1B => {
            // Escape sequence for arrows: ESC [ A/B
            let mut seq = [0u8; 2];
            if stdin.read_exact(&mut seq).is_ok() && seq[0] == b'[' {
                return match seq[1] {
                    b'A' => Ok(Key::Up),
                    b'B' => Ok(Key::Down),
                    _ => Ok(Key::Esc),
                };
            }
            Ok(Key::Esc)
        }
        b => Ok(Key::Char(b as char)),
    }
}

fn read_key_with_timeout(timeout_ms: u64) -> Result<Option<Key>> {
    use std::sync::mpsc;
    use std::thread;

    let (tx, rx) = mpsc::channel();

    thread::spawn(move || {
        if let Ok(key) = read_key() {
            let _ = tx.send(key);
        }
    });

    match rx.recv_timeout(Duration::from_millis(timeout_ms)) {
        Ok(key) => Ok(Some(key)),
        Err(mpsc::RecvTimeoutError::Timeout) => Ok(None),
        Err(mpsc::RecvTimeoutError::Disconnected) => Ok(None),
    }
}
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
        let expanded_body = expand_env_vars_in_text(&msg.body)?;
        let preview = if expanded_body.len() > preview_len {
            format!("{}...", &expanded_body[..preview_len])
        } else {
            expanded_body
        };
        println!("   {}", preview);
    }

    println!("\n-------------------------------------");
    println!("Tip: use 'bv inbox --plain' for non-interactive output");

    Ok(())
}

/// Interactive mode for inbox
pub async fn interactive(config: &Config, initial_view: Option<String>) -> Result<()> {
    let db_path = super::messages::get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    // Check if daemon is running - if so, don't sync manually
    let daemon_running = super::daemon::is_daemon_running(config).unwrap_or(false);
    if !daemon_running {
        // Sync messages first only if daemon is not running
        let sync = super::messages::init_message_system(config)?.1;
        let _ = sync.sync_quiet();
    }

    let mut current_view = initial_view.unwrap_or_else(|| "inbox".to_string());

    // Use raw mode and a simple key-driven UI (via stty)
    enable_raw_mode_cmd()?;
    let mut selected: usize = 0; // index into current messages; extra index for Quit
    let mut last_refresh = Instant::now();
    let mut last_message_count = 0;
    let refresh_interval = Duration::from_secs(5); // Refresh every 5 seconds when daemon is active

    loop {
        // Load messages for current view
        let messages = match current_view.as_str() {
            "inbox" => db.list_inbox_messages(Some(200))?,
            "sent" => db.list_sent_messages(Some(200))?,
            "all" => db.list_messages(Some(200))?,
            "unread" => db.list_unread_messages()?,
            "projects" => db.list_messages_by_type("project", Some(200))?,
            _ => db.list_inbox_messages(Some(200))?,
        };

        // Check if new messages arrived (for notification)
        let new_messages_arrived = messages.len() > last_message_count && last_message_count > 0;
        if new_messages_arrived && daemon_running {
            // Reset selection to top to highlight new messages
            selected = 0;
        }
        last_message_count = messages.len();

        // Bound selection to range [0 .. messages.len()] where last is Quit
        if selected > messages.len() {
            selected = messages.len();
        }

        // Render screen
        print!("\x1B[2J\x1B[1;1H");
        io::stdout().flush()?;
        println!("======================================================");
        println!(
            "BioVault Inbox - {} ({} messages){}",
            current_view.to_uppercase(),
            messages.len(),
            if new_messages_arrived && daemon_running { " 🆕" } else { "" }
        );

        // Show daemon status with auto-refresh indicator
        let daemon_status = if daemon_running {
            let time_since_refresh = Instant::now().duration_since(last_refresh);
            let refresh_countdown = refresh_interval.as_secs() - time_since_refresh.as_secs().min(refresh_interval.as_secs());
            format!("🤖 Daemon: ACTIVE | Auto-refresh in {}s", refresh_countdown)
        } else {
            "⚠️  Daemon: STOPPED".to_string()
        };
        println!("Status: {}", daemon_status);
        println!("======================================================");
        println!("(Press '?' for shortcuts)");
        println!("------------------------------------------------------");

        if messages.is_empty() {
            println!("(No messages in this view)");
        } else {
            for (i, msg) in messages.iter().enumerate() {
                let status = match msg.status {
                    crate::messages::MessageStatus::Draft => "DRAFT",
                    crate::messages::MessageStatus::Sent => "SENT",
                    crate::messages::MessageStatus::Received => "RECV",
                    crate::messages::MessageStatus::Read => "READ",
                    crate::messages::MessageStatus::Deleted => "DEL",
                    crate::messages::MessageStatus::Archived => "ARCH",
                };
                let who = if current_view == "sent" {
                    &msg.to
                } else {
                    &msg.from
                };
                let subject = msg.subject.as_deref().unwrap_or("(No Subject)");
                let mut line = format!(
                    "{} {status} {} - {}",
                    if i == selected { ">" } else { " " },
                    who,
                    subject
                );
                if line.len() > 80 {
                    line.truncate(80);
                }
                println!("{}", line);
            }
        }

        // Quit item at the bottom
        println!(
            "{} Quit",
            if selected == messages.len() { ">" } else { " " }
        );

        // Wait for a key event (with auto-refresh if daemon is running)
        let key_opt = if daemon_running {
            // Check if it's time to refresh
            if Instant::now().duration_since(last_refresh) >= refresh_interval {
                last_refresh = Instant::now();
                // Force a refresh by continuing the loop
                read_key_with_timeout(100)?
            } else {
                // Poll for key with short timeout
                read_key_with_timeout(500)?
            }
        } else {
            // Blocking read when daemon is not running
            Some(read_key()?)
        };

        // If no key was pressed and daemon is running, continue loop for refresh
        let key = match key_opt {
            Some(k) => k,
            None => {
                if daemon_running {
                    continue; // This will reload messages and refresh the display
                } else {
                    continue;
                }
            }
        };

        match key {
            Key::Char('q') | Key::Esc => {
                disable_raw_mode_cmd()?;
                break;
            }
            Key::Char('?') | Key::Char('h') | Key::Char('H') => {
                print!("\x1B[2J\x1B[1;1H");
                println!("Shortcuts:\n  ? / h : Help\n  n     : New Message\n  s     : Sync Messages\n  v     : Change View (menu)\n  q / Esc: Quit\n  1..5  : Tabs (Inbox, Sent, All, Unread, Projects)\n\nArrows to move, Enter to open.");
                println!("\nPress any key to return...");
                let _ = read_key();
            }
            Key::Char('n') | Key::Char('N') => {
                disable_raw_mode_cmd()?;
                compose_new_message(config)?;
                enable_raw_mode_cmd()?;
            }
            Key::Char('s') | Key::Char('S') => {
                disable_raw_mode_cmd()?;
                if !daemon_running {
                    println!("\nSyncing messages...");
                    let sync = super::messages::init_message_system(config)?.1;
                    let _ = sync.sync();
                    println!("Sync complete. Press Enter...");
                } else {
                    println!("\n⚠️  Daemon is running - sync happens automatically. Press Enter...");
                }
                let mut t = String::new();
                io::stdin().read_line(&mut t).ok();
                enable_raw_mode_cmd()?;
            }
            Key::Char('v') | Key::Char('V') => {
                disable_raw_mode_cmd()?;
                let views = ["inbox", "sent", "all", "unread", "projects"];
                let view_names = ["Inbox", "Sent", "All Messages", "Unread", "Projects"];
                println!("\nSelect view:");
                let view_selection = Select::with_theme(&ColorfulTheme::default())
                    .with_prompt("Choose a view")
                    .default(0)
                    .items(&view_names)
                    .interact_opt()?;
                if let Some(sel) = view_selection {
                    current_view = views[sel].to_string();
                    selected = 0;
                }
                enable_raw_mode_cmd()?;
            }
            Key::Char('1') => {
                current_view = "inbox".to_string();
                selected = 0;
            }
            Key::Char('2') => {
                current_view = "sent".to_string();
                selected = 0;
            }
            Key::Char('3') => {
                current_view = "all".to_string();
                selected = 0;
            }
            Key::Char('4') => {
                current_view = "unread".to_string();
                selected = 0;
            }
            Key::Char('5') => {
                current_view = "projects".to_string();
                selected = 0;
            }
            Key::Up => {
                selected = selected.saturating_sub(1);
            }
            Key::Down => {
                if selected < messages.len() {
                    selected += 1;
                }
            }
            Key::Enter => {
                if selected == messages.len() {
                    disable_raw_mode_cmd()?;
                    break;
                }
                if !messages.is_empty() {
                    let msg = &messages[selected];
                    disable_raw_mode_cmd()?;
                    let _ = message_actions(config, &db, msg).await?;
                    enable_raw_mode_cmd()?;
                }
            }
            Key::Char(_) => {}
        }
    }

    Ok(())
}

/// Handle actions for a selected message
/// Returns false if user selected "Back", true otherwise
async fn message_actions(
    config: &Config,
    db: &MessageDb,
    msg: &crate::messages::Message,
) -> Result<bool> {
    let mut actions = vec!["Read", "Reply", "Delete", "Mark as Read/Unread"];
    // If project message addressed to this user, add triage actions inline
    if let crate::messages::MessageType::Project { .. } = msg.message_type {
        if msg.to == config.email {
            actions.push("Run on test data");
            actions.push("Run on real data");
            actions.push("Reject");
            actions.push("Review");
            actions.push("Approve");
        }
        if msg.from == config.email {
            actions.push("Archive (finalize and revoke write)");
        }
    }
    actions.push("Back to Messages");
    // Capture dynamic action indexes for later dispatch
    let idx_run_test = actions.iter().position(|s| *s == "Run on test data");
    let idx_run_real = actions.iter().position(|s| *s == "Run on real data");
    let idx_reject = actions.iter().position(|s| *s == "Reject");
    let idx_review = actions.iter().position(|s| *s == "Review");
    let idx_approve = actions.iter().position(|s| *s == "Approve");
    let idx_archive = actions
        .iter()
        .position(|s| *s == "Archive (finalize and revoke write)");

    let selection = Select::with_theme(&ColorfulTheme::default())
        .with_prompt("Select action")
        .default(0)
        .items(&actions)
        .interact_opt()?;

    match selection {
        Some(0) => {
            // Read
            super::messages::read_message(config, &msg.id).await?;
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
        // Dispatch dynamic actions if present
        Some(idx) => {
            use super::messages::{perform_project_action, ProjectAction};
            if let Some(i) = idx_run_test {
                if idx == i {
                    perform_project_action(config, &msg.id, ProjectAction::RunTest).await?;
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            if let Some(i) = idx_run_real {
                if idx == i {
                    perform_project_action(config, &msg.id, ProjectAction::RunReal).await?;
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            if let Some(i) = idx_reject {
                if idx == i {
                    perform_project_action(config, &msg.id, ProjectAction::Reject).await?;
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            if let Some(i) = idx_review {
                if idx == i {
                    perform_project_action(config, &msg.id, ProjectAction::Review).await?;
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            if let Some(i) = idx_approve {
                if idx == i {
                    perform_project_action(config, &msg.id, ProjectAction::Approve).await?;
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            if let Some(i) = idx_archive {
                if idx == i {
                    super::messages::read_message(config, &msg.id).await?; // Archive via read view
                    println!("\nPress Enter to continue...");
                    let mut input = String::new();
                    io::stdin().read_line(&mut input)?;
                    return Ok(false);
                }
            }
            Ok(false)
        }
        None => Ok(false),
    }
}

/// Compose and send a new message interactively
fn compose_new_message(config: &Config) -> Result<()> {
    println!("\nCompose New Message");
    println!("--------------------");

    // Recipient
    print!("Recipient email: ");
    io::stdout().flush()?;
    let mut recipient = String::new();
    io::stdin().read_line(&mut recipient)?;
    let recipient = recipient.trim().to_string();
    if recipient.is_empty() {
        println!("Cancelled (no recipient)");
        println!("Press Enter to continue...");
        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        return Ok(());
    }

    // Subject (optional)
    print!("Subject (optional): ");
    io::stdout().flush()?;
    let mut subject = String::new();
    io::stdin().read_line(&mut subject)?;
    let subject = subject.trim().to_string();
    let subject_opt = if subject.is_empty() {
        None
    } else {
        Some(subject.as_str())
    };

    // Body (single line for simplicity)
    println!("Body (single line, press Enter to finish):");
    print!("> ");
    io::stdout().flush()?;
    let mut body = String::new();
    io::stdin().read_line(&mut body)?;
    let body = body.trim();
    if body.is_empty() {
        println!("Cancelled (empty body)");
        println!("Press Enter to continue...");
        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        return Ok(());
    }

    super::messages::send_message(config, &recipient, body, subject_opt)?;
    println!("\nMessage sent. Press Enter to continue...");
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_list_filters_default() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: false,
            message_type: None,
            from: None,
            search: None,
        };
        assert!(!filters.sent);
        assert!(!filters.all);
        assert!(!filters.unread);
        assert!(!filters.projects);
        assert!(filters.message_type.is_none());
        assert!(filters.from.is_none());
        assert!(filters.search.is_none());
    }

    #[test]
    fn test_key_enum_equality() {
        assert_eq!(Key::Up, Key::Up);
        assert_eq!(Key::Down, Key::Down);
        assert_eq!(Key::Enter, Key::Enter);
        assert_eq!(Key::Esc, Key::Esc);
        assert_eq!(Key::Char('a'), Key::Char('a'));
        assert_ne!(Key::Char('a'), Key::Char('b'));
        assert_ne!(Key::Up, Key::Down);
    }

    #[test]
    fn test_key_enum_debug() {
        assert_eq!(format!("{:?}", Key::Up), "Up");
        assert_eq!(format!("{:?}", Key::Down), "Down");
        assert_eq!(format!("{:?}", Key::Enter), "Enter");
        assert_eq!(format!("{:?}", Key::Esc), "Esc");
        assert_eq!(format!("{:?}", Key::Char('x')), "Char('x')");
    }

    #[test]
    fn test_list_filters_with_values() {
        let filters = ListFilters {
            sent: true,
            all: false,
            unread: true,
            projects: false,
            message_type: Some("project".to_string()),
            from: Some("test@example.com".to_string()),
            search: Some("search term".to_string()),
        };
        assert!(filters.sent);
        assert!(!filters.all);
        assert!(filters.unread);
        assert!(!filters.projects);
        assert_eq!(filters.message_type, Some("project".to_string()));
        assert_eq!(filters.from, Some("test@example.com".to_string()));
        assert_eq!(filters.search, Some("search term".to_string()));
    }

    #[test]
    fn test_key_char_creation() {
        let key = Key::Char('a');
        match key {
            Key::Char(c) => assert_eq!(c, 'a'),
            _ => panic!("Expected Char variant"),
        }
    }

    #[test]
    fn test_key_copy_trait() {
        let key1 = Key::Up;
        let key2 = key1; // This works because Key implements Copy
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_list_filters_sent_filter() {
        let filters = ListFilters {
            sent: true,
            all: false,
            unread: false,
            projects: false,
            message_type: None,
            from: None,
            search: None,
        };
        assert!(filters.sent);
        assert!(!filters.all);
    }

    #[test]
    fn test_list_filters_search_filter() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: false,
            message_type: None,
            from: None,
            search: Some("test query".to_string()),
        };
        assert_eq!(filters.search, Some("test query".to_string()));
    }

    #[test]
    fn test_list_filters_message_type_filter() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: false,
            message_type: Some("request".to_string()),
            from: None,
            search: None,
        };
        assert_eq!(filters.message_type, Some("request".to_string()));
    }

    #[test]
    fn test_list_filters_from_filter() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: false,
            message_type: None,
            from: Some("user@example.com".to_string()),
            search: None,
        };
        assert_eq!(filters.from, Some("user@example.com".to_string()));
    }

    #[test]
    fn test_list_filters_projects_filter() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: true,
            message_type: None,
            from: None,
            search: None,
        };
        assert!(filters.projects);
        assert!(!filters.sent);
    }

    #[test]
    fn test_list_filters_unread_filter() {
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: true,
            projects: false,
            message_type: None,
            from: None,
            search: None,
        };
        assert!(filters.unread);
        assert!(!filters.all);
    }

    #[test]
    fn test_list_filters_all_filter() {
        let filters = ListFilters {
            sent: false,
            all: true,
            unread: false,
            projects: false,
            message_type: None,
            from: None,
            search: None,
        };
        assert!(filters.all);
        assert!(!filters.sent);
    }

    #[test]
    fn test_list_filters_combined() {
        let filters = ListFilters {
            sent: true,
            all: false,
            unread: true,
            projects: false,
            message_type: Some("system".to_string()),
            from: Some("admin@test.com".to_string()),
            search: Some("important".to_string()),
        };
        assert!(filters.sent);
        assert!(filters.unread);
        assert_eq!(filters.message_type, Some("system".to_string()));
        assert_eq!(filters.from, Some("admin@test.com".to_string()));
        assert_eq!(filters.search, Some("important".to_string()));
    }

    fn test_config(tmp: &std::path::Path) -> Config {
        crate::config::set_test_biovault_home(tmp.join(".bv"));
        crate::config::set_test_syftbox_data_dir(tmp);
        Config {
            email: "me@example.com".into(),
            syftbox_config: None,
            version: None,
        }
    }

    #[test]
    fn list_handles_empty_db() {
        let tmp = tempfile::TempDir::new().unwrap();
        let cfg = test_config(tmp.path());
        let filters = ListFilters {
            sent: false,
            all: false,
            unread: false,
            projects: false,
            message_type: None,
            from: None,
            search: None,
        };
        list(&cfg, filters).unwrap();
    }

    #[test]
    fn list_filters_sent_unread_projects_search_from() {
        use crate::cli::commands::messages::get_message_db_path;
        use crate::messages::{Message, MessageStatus, MessageType};

        let tmp = tempfile::TempDir::new().unwrap();
        let cfg = test_config(tmp.path());
        let db_path = get_message_db_path(&cfg).unwrap();
        let db = crate::messages::MessageDb::new(&db_path).unwrap();

        // Sent/draft message from me
        let mut m1 = Message::new(
            "me@example.com".into(),
            "you@example.com".into(),
            "hello world".into(),
        );
        m1.status = MessageStatus::Sent;
        db.insert_message(&m1).unwrap();

        // Unread (received) message to me, from alice
        let mut m2 = Message::new(
            "alice@example.com".into(),
            "me@example.com".into(),
            "body 2".into(),
        );
        m2.status = MessageStatus::Received;
        m2.subject = Some("greetings".into());
        db.insert_message(&m2).unwrap();

        // Project-type message (counted in projects filter)
        let mut m3 = Message::new(
            "bob@example.com".into(),
            "me@example.com".into(),
            "proj body".into(),
        );
        m3.message_type = MessageType::Project {
            project_name: "P".into(),
            submission_id: "S".into(),
            files_hash: None,
        };
        m3.status = MessageStatus::Sent;
        db.insert_message(&m3).unwrap();

        // sent filter
        list(
            &cfg,
            ListFilters {
                sent: true,
                all: false,
                unread: false,
                projects: false,
                message_type: None,
                from: None,
                search: None,
            },
        )
        .unwrap();
        // unread filter
        list(
            &cfg,
            ListFilters {
                sent: false,
                all: false,
                unread: true,
                projects: false,
                message_type: None,
                from: None,
                search: None,
            },
        )
        .unwrap();
        // projects filter
        list(
            &cfg,
            ListFilters {
                sent: false,
                all: false,
                unread: false,
                projects: true,
                message_type: None,
                from: None,
                search: None,
            },
        )
        .unwrap();
        // message_type filter explicit
        list(
            &cfg,
            ListFilters {
                sent: false,
                all: false,
                unread: false,
                projects: false,
                message_type: Some("project".into()),
                from: None,
                search: None,
            },
        )
        .unwrap();
        // search filter by subject/body
        list(
            &cfg,
            ListFilters {
                sent: false,
                all: false,
                unread: false,
                projects: false,
                message_type: None,
                from: None,
                search: Some("greet".into()),
            },
        )
        .unwrap();
        // from filter narrows by sender
        list(
            &cfg,
            ListFilters {
                sent: false,
                all: true,
                unread: false,
                projects: false,
                message_type: None,
                from: Some("alice".into()),
                search: None,
            },
        )
        .unwrap();
    }
}
