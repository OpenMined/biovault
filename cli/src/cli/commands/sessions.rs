//! CLI commands for managing collaborative Jupyter sessions
//!
//! Sessions allow data scientists to create collaborative Jupyter environments
//! with peer invitations for secure data sharing.

use crate::config::Config;
use crate::data::BioVaultDb;
use crate::messages::{Message, MessageDb, MessageStatus};
use anyhow::{Context, Result};
use chrono::Utc;
use clap::Subcommand;
use colored::Colorize;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum SessionsCommands {
    #[command(about = "Create a new collaborative session")]
    Create {
        #[arg(help = "Session name")]
        name: String,

        #[arg(long, help = "Peer email to invite")]
        peer: Option<String>,

        #[arg(long, help = "Session description")]
        description: Option<String>,
    },

    #[command(about = "List all sessions")]
    List {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Show session details")]
    Show {
        #[arg(help = "Session ID")]
        session_id: String,

        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Invite a peer to an existing session")]
    Invite {
        #[arg(help = "Session ID")]
        session_id: String,

        #[arg(help = "Peer email to invite")]
        peer: String,
    },

    #[command(about = "List pending session invitations")]
    Invitations {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Accept a session invitation")]
    Accept {
        #[arg(help = "Session ID from invitation")]
        session_id: String,
    },

    #[command(about = "Reject a session invitation")]
    Reject {
        #[arg(help = "Session ID from invitation")]
        session_id: String,

        #[arg(long, help = "Reason for rejection")]
        reason: Option<String>,
    },

    #[command(about = "Send a chat message in a session")]
    Chat {
        #[arg(help = "Session ID")]
        session_id: String,

        #[arg(help = "Message text")]
        message: String,
    },

    #[command(about = "Delete a session")]
    Delete {
        #[arg(help = "Session ID")]
        session_id: String,

        #[arg(long, help = "Skip confirmation prompt")]
        yes: bool,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Session {
    pub id: i64,
    pub session_id: String,
    pub name: String,
    pub description: Option<String>,
    pub session_path: String,
    pub owner: String,
    pub peer: Option<String>,
    pub role: String,
    pub status: String,
    pub created_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionInvitation {
    pub session_id: String,
    pub requester: String,
    pub target: String,
    pub session_name: Option<String>,
    pub description: Option<String>,
    pub created_at: String,
    pub status: String,
}

/// Handle sessions subcommands
pub async fn handle(command: SessionsCommands, config: &Config) -> Result<()> {
    match command {
        SessionsCommands::Create {
            name,
            peer,
            description,
        } => create_session(config, &name, peer, description),
        SessionsCommands::List { json } => list_sessions(config, json),
        SessionsCommands::Show { session_id, json } => show_session(config, &session_id, json),
        SessionsCommands::Invite { session_id, peer } => invite_peer(config, &session_id, &peer),
        SessionsCommands::Invitations { json } => list_invitations(config, json),
        SessionsCommands::Accept { session_id } => accept_invitation(config, &session_id),
        SessionsCommands::Reject { session_id, reason } => {
            reject_invitation(config, &session_id, reason)
        }
        SessionsCommands::Chat {
            session_id,
            message,
        } => send_chat_message(config, &session_id, &message),
        SessionsCommands::Delete { session_id, yes } => delete_session(config, &session_id, yes),
    }
}

fn generate_session_id() -> String {
    let mut rng = rand::rng();
    let bytes: [u8; 6] = rng.random();
    hex::encode(bytes)
}

fn get_sessions_dir(config: &Config) -> Result<PathBuf> {
    let biovault_home = config.get_biovault_dir()?;
    Ok(biovault_home.join("sessions"))
}

fn get_rpc_session_dir(config: &Config, target_email: &str) -> Result<PathBuf> {
    let data_dir = config.get_syftbox_data_dir()?;
    Ok(data_dir
        .join("datasites")
        .join(target_email)
        .join("app_data")
        .join("biovault")
        .join("rpc")
        .join("session"))
}

fn get_my_rpc_session_dir(config: &Config) -> Result<PathBuf> {
    get_rpc_session_dir(config, &config.email)
}

/// Create a new session
fn create_session(
    config: &Config,
    name: &str,
    peer: Option<String>,
    description: Option<String>,
) -> Result<()> {
    let db = BioVaultDb::new()?;
    let session_id = generate_session_id();
    let sessions_dir = get_sessions_dir(config)?;
    let session_path = sessions_dir.join(&session_id);
    let owner = config.email.clone();

    // Create session directory
    fs::create_dir_all(&session_path).context("Failed to create session directory")?;
    fs::create_dir_all(session_path.join("data")).context("Failed to create data directory")?;

    let session_path_str = session_path.to_string_lossy().to_string();

    // Insert into database
    db.connection().execute(
        "INSERT INTO sessions (session_id, name, description, session_path, owner, peer, role, status)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, 'owner', 'active')",
        rusqlite::params![&session_id, name, &description, &session_path_str, &owner, &peer],
    )?;

    // Write session config file
    let session_config = serde_json::json!({
        "session_id": &session_id,
        "name": name,
        "owner": &owner,
        "peer": &peer,
        "created_at": Utc::now().to_rfc3339(),
    });
    let config_path = session_path.join("session.json");
    fs::write(&config_path, serde_json::to_string_pretty(&session_config)?)?;

    println!("\n✅ Session created: {}", name.green());
    println!("   📁 ID: {}", session_id.cyan());
    println!("   📂 Path: {}", session_path_str);

    // Send invitation if peer specified
    if let Some(peer_email) = &peer {
        send_session_invitation(config, &session_id, name, &owner, peer_email, &description)?;
        println!("   📨 Invitation sent to: {}", peer_email.cyan());
    }

    println!(
        "\nStart Jupyter with: {}",
        format!("bv jupyter start {}", session_path_str).cyan()
    );

    Ok(())
}

/// Send a session invitation via RPC and messaging
fn send_session_invitation(
    config: &Config,
    session_id: &str,
    session_name: &str,
    owner: &str,
    peer_email: &str,
    description: &Option<String>,
) -> Result<()> {
    // Create RPC folder structure for the peer
    let rpc_path = get_rpc_session_dir(config, peer_email)?;
    fs::create_dir_all(&rpc_path)?;

    let invitation = serde_json::json!({
        "session_id": session_id,
        "requester": owner,
        "target": peer_email,
        "session_name": session_name,
        "description": description,
        "created_at": Utc::now().to_rfc3339(),
        "message": format!("{} invites you to a BioVault session", owner),
        "status": "pending"
    });

    let request_file = rpc_path.join(format!("{}.request", session_id));
    fs::write(&request_file, serde_json::to_string_pretty(&invitation)?)?;

    // Also send via messaging system
    let db_path = crate::cli::commands::messages::get_message_db_path(config)?;
    let db = MessageDb::new(&db_path)?;

    let metadata = serde_json::json!({
        "session_invite": {
            "session_id": session_id,
            "session_name": session_name,
            "from": owner,
            "description": description,
            "created_at": Utc::now().to_rfc3339(),
        }
    });

    let mut msg = Message::new(
        owner.to_string(),
        peer_email.to_string(),
        format!(
            "{} invited you to the session \"{}\" (ID: {}).",
            owner, session_name, session_id
        ),
    );
    msg.subject = Some(format!("Session Invite: {}", session_name));
    msg.metadata = Some(metadata);
    msg.status = MessageStatus::Sent;

    db.insert_message(&msg)?;

    // Try to send via RPC sync
    let data_dir = config.get_syftbox_data_dir()?;
    let app = crate::syftbox::SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
    let sync = crate::messages::MessageSync::new(&db_path, app)?;
    let _ = sync.send_message(&msg.id);

    Ok(())
}

/// List all sessions
fn list_sessions(_config: &Config, json_output: bool) -> Result<()> {
    let db = BioVaultDb::new()?;

    let sessions: Vec<Session> = db
        .connection()
        .prepare(
            "SELECT id, session_id, name, description, session_path, owner, peer, role, status, created_at
             FROM sessions ORDER BY created_at DESC",
        )?
        .query_map([], |row| {
            Ok(Session {
                id: row.get(0)?,
                session_id: row.get(1)?,
                name: row.get(2)?,
                description: row.get(3)?,
                session_path: row.get(4)?,
                owner: row.get(5)?,
                peer: row.get(6)?,
                role: row.get(7)?,
                status: row.get(8)?,
                created_at: row.get(9)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;

    if json_output {
        println!("{}", serde_json::to_string_pretty(&sessions)?);
        return Ok(());
    }

    if sessions.is_empty() {
        println!("No sessions yet.");
        println!(
            "\nCreate one with: {}",
            "bv session create \"My Session\"".cyan()
        );
        return Ok(());
    }

    println!("\n🔬 {} Session(s)", "Your".bold());
    println!("═══════════════════════════════════════════════════════════════\n");

    for session in &sessions {
        let status_icon = match session.status.as_str() {
            "active" => "🟢",
            "pending" => "🟡",
            "closed" => "⚫",
            _ => "⚪",
        };

        let role_badge = if session.role == "owner" {
            "(owner)".dimmed()
        } else {
            "(peer)".blue()
        };

        println!(
            "  {} {} {} [{}]",
            status_icon,
            session.name.green(),
            role_badge,
            session.session_id.cyan()
        );

        if let Some(peer) = &session.peer {
            println!("     👤 with: {}", peer);
        }
        if let Some(desc) = &session.description {
            if !desc.is_empty() {
                println!("     📝 {}", desc.dimmed());
            }
        }
        println!();
    }

    println!("Total: {} session(s)", sessions.len());

    Ok(())
}

/// Show session details
fn show_session(_config: &Config, session_id: &str, json_output: bool) -> Result<()> {
    let db = BioVaultDb::new()?;

    let session: Session = db
        .connection()
        .query_row(
            "SELECT id, session_id, name, description, session_path, owner, peer, role, status, created_at
             FROM sessions WHERE session_id = ?1",
            [session_id],
            |row| {
                Ok(Session {
                    id: row.get(0)?,
                    session_id: row.get(1)?,
                    name: row.get(2)?,
                    description: row.get(3)?,
                    session_path: row.get(4)?,
                    owner: row.get(5)?,
                    peer: row.get(6)?,
                    role: row.get(7)?,
                    status: row.get(8)?,
                    created_at: row.get(9)?,
                })
            },
        )
        .context(format!("Session not found: {}", session_id))?;

    if json_output {
        println!("{}", serde_json::to_string_pretty(&session)?);
        return Ok(());
    }

    println!("\n🔬 Session: {}", session.name.green());
    println!("═══════════════════════════════════════════════════════════════");
    println!("   ID:          {}", session.session_id.cyan());
    println!("   Status:      {}", session.status);
    println!("   Role:        {}", session.role);
    println!("   Owner:       {}", session.owner);
    println!(
        "   Peer:        {}",
        session.peer.as_deref().unwrap_or("(none)")
    );
    println!("   Path:        {}", session.session_path);
    println!("   Created:     {}", session.created_at);

    if let Some(desc) = &session.description {
        if !desc.is_empty() {
            println!("   Description: {}", desc);
        }
    }

    println!(
        "\nStart Jupyter: {}",
        format!("bv jupyter start {}", session.session_path).cyan()
    );

    Ok(())
}

/// Invite a peer to an existing session
fn invite_peer(config: &Config, session_id: &str, peer: &str) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get session details
    let (name, owner, description): (String, String, Option<String>) = db
        .connection()
        .query_row(
            "SELECT name, owner, description FROM sessions WHERE session_id = ?1",
            [session_id],
            |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
        )
        .context(format!("Session not found: {}", session_id))?;

    // Update session with peer
    db.connection().execute(
        "UPDATE sessions SET peer = ?1, updated_at = CURRENT_TIMESTAMP WHERE session_id = ?2",
        rusqlite::params![peer, session_id],
    )?;

    // Send invitation
    send_session_invitation(config, session_id, &name, &owner, peer, &description)?;

    println!("\n📨 Invitation sent to: {}", peer.green());
    println!("   Session: {} [{}]", name, session_id.cyan());

    Ok(())
}

/// List pending invitations
fn list_invitations(config: &Config, json_output: bool) -> Result<()> {
    let my_rpc_dir = get_my_rpc_session_dir(config)?;
    let mut invitations = Vec::new();

    // Check RPC folder for .request files
    if my_rpc_dir.exists() {
        for entry in fs::read_dir(&my_rpc_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.extension().map(|e| e == "request").unwrap_or(false) {
                if let Ok(content) = fs::read_to_string(&path) {
                    if let Ok(invite) = serde_json::from_str::<SessionInvitation>(&content) {
                        if invite.status == "pending"
                            && !session_exists(config, &invite.session_id)?
                        {
                            invitations.push(invite);
                        }
                    }
                }
            }
        }
    }

    // Also check messaging system for session invites
    let db_path = crate::cli::commands::messages::get_message_db_path(config)?;
    if let Ok(db) = MessageDb::new(&db_path) {
        if let Ok(messages) = db.list_messages(None) {
            for msg in messages {
                if msg.to != config.email {
                    continue;
                }
                if !matches!(msg.status, MessageStatus::Received | MessageStatus::Read) {
                    continue;
                }
                if let Some(meta) = msg.metadata {
                    if let Some(invite) = meta.get("session_invite") {
                        if let Some(sid) = invite.get("session_id").and_then(|v| v.as_str()) {
                            if session_exists(config, sid)? {
                                continue;
                            }
                            let session_name = invite
                                .get("session_name")
                                .and_then(|v| v.as_str())
                                .map(String::from);
                            let description = invite
                                .get("description")
                                .and_then(|v| v.as_str())
                                .map(String::from);
                            let created_at = invite
                                .get("created_at")
                                .and_then(|v| v.as_str())
                                .map(String::from)
                                .unwrap_or_else(|| msg.created_at.to_rfc3339());
                            let requester = invite
                                .get("from")
                                .and_then(|v| v.as_str())
                                .unwrap_or(&msg.from)
                                .to_string();

                            // Avoid duplicates
                            if !invitations.iter().any(|i| i.session_id == sid) {
                                invitations.push(SessionInvitation {
                                    session_id: sid.to_string(),
                                    requester,
                                    target: config.email.clone(),
                                    session_name,
                                    description,
                                    created_at,
                                    status: "pending".to_string(),
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    invitations.sort_by(|a, b| b.created_at.cmp(&a.created_at));

    if json_output {
        println!("{}", serde_json::to_string_pretty(&invitations)?);
        return Ok(());
    }

    if invitations.is_empty() {
        println!("No pending session invitations.");
        return Ok(());
    }

    println!(
        "\n📬 {} Pending Invitation(s)",
        invitations.len().to_string().yellow()
    );
    println!("═══════════════════════════════════════════════════════════════\n");

    for invite in &invitations {
        let name = invite.session_name.as_deref().unwrap_or("Unnamed Session");
        println!("  📋 {} [{}]", name.green(), invite.session_id.cyan());
        println!("     From: {}", invite.requester);
        if let Some(desc) = &invite.description {
            if !desc.is_empty() {
                println!("     📝 {}", desc.dimmed());
            }
        }
        println!(
            "     Accept: {}",
            format!("bv session accept {}", invite.session_id).cyan()
        );
        println!();
    }

    Ok(())
}

fn session_exists(_config: &Config, session_id: &str) -> Result<bool> {
    let db = BioVaultDb::new()?;
    let exists: Option<i64> = db
        .connection()
        .query_row(
            "SELECT id FROM sessions WHERE session_id = ?1 LIMIT 1",
            [session_id],
            |row| row.get(0),
        )
        .ok();
    Ok(exists.is_some())
}

/// Accept a session invitation
fn accept_invitation(config: &Config, session_id: &str) -> Result<()> {
    let my_rpc_dir = get_my_rpc_session_dir(config)?;
    let request_file = my_rpc_dir.join(format!("{}.request", session_id));

    if !request_file.exists() {
        anyhow::bail!(
            "Invitation not found: {}. Use 'bv session invitations' to see pending invites.",
            session_id
        );
    }

    let content = fs::read_to_string(&request_file)?;
    let mut invitation: SessionInvitation = serde_json::from_str(&content)?;

    if invitation.session_name.is_none() {
        invitation.session_name = Some(format!("Session with {}", invitation.requester));
    }

    // Update invitation status
    let mut updated = invitation.clone();
    updated.status = "accepted".to_string();
    fs::write(&request_file, serde_json::to_string_pretty(&updated)?)?;

    // Create local session
    let sessions_dir = get_sessions_dir(config)?;
    let session_path = sessions_dir.join(session_id);
    fs::create_dir_all(&session_path)?;
    fs::create_dir_all(session_path.join("data"))?;

    let session_name = invitation
        .session_name
        .clone()
        .unwrap_or_else(|| format!("Session with {}", invitation.requester));

    let db = BioVaultDb::new()?;
    db.connection().execute(
        "INSERT INTO sessions (session_id, name, description, session_path, owner, peer, role, status)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, 'peer', 'active')",
        rusqlite::params![
            session_id,
            &session_name,
            &invitation.description,
            session_path.to_string_lossy().to_string(),
            &config.email,
            &invitation.requester,
        ],
    )?;

    // Write session config
    let session_config = serde_json::json!({
        "session_id": session_id,
        "name": &session_name,
        "owner": &config.email,
        "peer": &invitation.requester,
        "created_at": Utc::now().to_rfc3339(),
    });
    fs::write(
        session_path.join("session.json"),
        serde_json::to_string_pretty(&session_config)?,
    )?;

    // Send acceptance response via RPC
    let requester_rpc = get_rpc_session_dir(config, &invitation.requester)?;
    fs::create_dir_all(&requester_rpc)?;

    let response = serde_json::json!({
        "session_id": session_id,
        "status": "accepted",
        "accepted_at": Utc::now().to_rfc3339(),
        "responder": &config.email,
        "session_name": &session_name,
    });
    fs::write(
        requester_rpc.join(format!("{}.response", session_id)),
        serde_json::to_string_pretty(&response)?,
    )?;

    println!(
        "\n✅ Accepted invitation from {}",
        invitation.requester.green()
    );
    println!("   Session: {} [{}]", session_name, session_id.cyan());
    println!(
        "\nStart Jupyter: {}",
        format!("bv jupyter start {}", session_path.display()).cyan()
    );

    Ok(())
}

/// Reject a session invitation
fn reject_invitation(config: &Config, session_id: &str, reason: Option<String>) -> Result<()> {
    let my_rpc_dir = get_my_rpc_session_dir(config)?;
    let request_file = my_rpc_dir.join(format!("{}.request", session_id));

    if !request_file.exists() {
        anyhow::bail!("Invitation not found: {}", session_id);
    }

    let content = fs::read_to_string(&request_file)?;
    let invitation: SessionInvitation = serde_json::from_str(&content)?;

    // Update invitation status
    let mut updated = invitation.clone();
    updated.status = "rejected".to_string();
    fs::write(&request_file, serde_json::to_string_pretty(&updated)?)?;

    // Send rejection response
    let requester_rpc = get_rpc_session_dir(config, &invitation.requester)?;
    fs::create_dir_all(&requester_rpc)?;

    let response = serde_json::json!({
        "session_id": session_id,
        "status": "rejected",
        "rejected_at": Utc::now().to_rfc3339(),
        "reason": reason,
        "responder": &config.email,
    });
    fs::write(
        requester_rpc.join(format!("{}.response", session_id)),
        serde_json::to_string_pretty(&response)?,
    )?;

    println!(
        "\n❌ Rejected invitation from {}",
        invitation.requester.yellow()
    );

    Ok(())
}

/// Send a chat message in a session
fn send_chat_message(config: &Config, session_id: &str, message: &str) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get session details
    let (name, peer, role): (String, Option<String>, String) = db
        .connection()
        .query_row(
            "SELECT name, peer, role FROM sessions WHERE session_id = ?1",
            [session_id],
            |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
        )
        .context(format!("Session not found: {}", session_id))?;

    let recipient = if role == "owner" {
        peer.ok_or_else(|| anyhow::anyhow!("No peer set for this session"))?
    } else {
        // Get owner from session
        db.connection().query_row(
            "SELECT owner FROM sessions WHERE session_id = ?1",
            [session_id],
            |row| row.get(0),
        )?
    };

    // Send via messaging system
    let db_path = crate::cli::commands::messages::get_message_db_path(config)?;
    let msg_db = MessageDb::new(&db_path)?;

    let metadata = serde_json::json!({
        "session_chat": {
            "session_id": session_id,
            "session_name": name,
            "from": &config.email,
            "created_at": Utc::now().to_rfc3339(),
        }
    });

    let mut msg = Message::new(config.email.clone(), recipient.clone(), message.to_string());
    msg.subject = Some(format!("Session: {}", name));
    msg.metadata = Some(metadata);
    msg.status = MessageStatus::Sent;

    msg_db.insert_message(&msg)?;

    // Try to send via RPC
    let data_dir = config.get_syftbox_data_dir()?;
    let app = crate::syftbox::SyftBoxApp::new(&data_dir, &config.email, "biovault")?;
    let sync = crate::messages::MessageSync::new(&db_path, app)?;
    let _ = sync.send_message(&msg.id);

    println!(
        "💬 Message sent to {} in session {}",
        recipient.green(),
        name.cyan()
    );

    Ok(())
}

/// Delete a session
fn delete_session(_config: &Config, session_id: &str, skip_confirm: bool) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get session details
    let session_path: String = db
        .connection()
        .query_row(
            "SELECT session_path FROM sessions WHERE session_id = ?1",
            [session_id],
            |row| row.get(0),
        )
        .context(format!("Session not found: {}", session_id))?;

    if !skip_confirm {
        use dialoguer::Confirm;
        let confirmed = Confirm::new()
            .with_prompt(format!(
                "Delete session {}?\nThe session folder will be preserved at: {}",
                session_id, session_path
            ))
            .default(false)
            .interact()
            .unwrap_or(false);

        if !confirmed {
            println!("Cancelled.");
            return Ok(());
        }
    }

    // Delete from database
    db.connection()
        .execute("DELETE FROM sessions WHERE session_id = ?1", [session_id])?;

    println!("\n🗑️  Deleted session: {}", session_id.yellow());
    println!("   Files preserved at: {}", session_path);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn create_test_config(temp: &TempDir) -> Config {
        crate::config::set_test_syftbox_data_dir(temp.path());
        crate::config::set_test_biovault_home(temp.path().join(".biovault"));

        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
        }
    }

    #[test]
    fn test_generate_session_id() {
        let id1 = generate_session_id();
        let id2 = generate_session_id();

        assert_eq!(id1.len(), 12); // 6 bytes = 12 hex chars
        assert_ne!(id1, id2);
    }

    #[test]
    fn test_list_sessions_empty() {
        let temp = TempDir::new().unwrap();
        let config = create_test_config(&temp);

        // Should not error on empty database
        let result = list_sessions(&config, false);
        assert!(result.is_ok());

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }
}
