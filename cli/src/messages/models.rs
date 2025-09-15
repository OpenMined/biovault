use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Message {
    pub id: String,
    pub thread_id: Option<String>,
    pub parent_id: Option<String>,

    pub from: String, // Email for now, but kept as string for future
    pub to: String,

    pub subject: Option<String>, // Optional, defaults to "(No Subject)"
    pub body: String,

    pub status: MessageStatus,
    pub sync_status: SyncStatus,

    pub created_at: DateTime<Utc>,
    pub sent_at: Option<DateTime<Utc>>,
    pub received_at: Option<DateTime<Utc>>,
    pub read_at: Option<DateTime<Utc>>,

    pub rpc_request_id: Option<String>,
    pub rpc_ack_status: Option<i32>,
    pub rpc_ack_at: Option<DateTime<Utc>>,
}

impl Message {
    pub fn new(from: String, to: String, body: String) -> Self {
        let id = Uuid::new_v4().to_string();
        Self {
            id: id.clone(),
            thread_id: Some(id.clone()), // New messages start their own thread
            parent_id: None,
            from,
            to,
            subject: None,
            body,
            status: MessageStatus::Draft,
            sync_status: SyncStatus::Local,
            created_at: Utc::now(),
            sent_at: None,
            received_at: None,
            read_at: None,
            rpc_request_id: None,
            rpc_ack_status: None,
            rpc_ack_at: None,
        }
    }

    pub fn reply_to(original: &Message, from: String, body: String) -> Self {
        Self {
            id: Uuid::new_v4().to_string(),
            thread_id: original.thread_id.clone(),
            parent_id: Some(original.id.clone()),
            from,
            to: original.from.clone(), // Reply goes back to sender
            subject: None,             // Will be handled in display
            body,
            status: MessageStatus::Draft,
            sync_status: SyncStatus::Local,
            created_at: Utc::now(),
            sent_at: None,
            received_at: None,
            read_at: None,
            rpc_request_id: None,
            rpc_ack_status: None,
            rpc_ack_at: None,
        }
    }

    pub fn display_subject(&self) -> String {
        self.subject
            .clone()
            .unwrap_or_else(|| "(No Subject)".to_string())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum MessageStatus {
    Draft,
    Sent,
    Received,
    Read,
    Deleted,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SyncStatus {
    Local,
    Syncing,
    Synced,
    Failed,
}

impl std::fmt::Display for MessageStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MessageStatus::Draft => write!(f, "draft"),
            MessageStatus::Sent => write!(f, "sent"),
            MessageStatus::Received => write!(f, "received"),
            MessageStatus::Read => write!(f, "read"),
            MessageStatus::Deleted => write!(f, "deleted"),
        }
    }
}

impl std::fmt::Display for SyncStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SyncStatus::Local => write!(f, "local"),
            SyncStatus::Syncing => write!(f, "syncing"),
            SyncStatus::Synced => write!(f, "synced"),
            SyncStatus::Failed => write!(f, "failed"),
        }
    }
}
