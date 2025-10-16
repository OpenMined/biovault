use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
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

    pub message_type: MessageType,
    pub metadata: Option<JsonValue>,

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
            message_type: MessageType::Text,
            metadata: None,
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
            message_type: MessageType::Text,
            metadata: None,
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub enum MessageType {
    #[default]
    Text,
    Project {
        project_name: String,
        submission_id: String,
        files_hash: Option<String>,
    },
    Request {
        request_type: String,
        params: Option<JsonValue>,
    },
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ThreadFilter {
    All,
    Inbox,
    Sent,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MessageThreadSummary {
    pub thread_id: String,
    pub subject: String,
    pub participants: Vec<String>,
    pub last_message_preview: String,
    pub last_message_at: DateTime<Utc>,
    pub last_message_id: String,
    pub last_message_status: MessageStatus,
    pub unread_count: usize,
    pub total_messages: usize,
    pub has_project: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub project_name: Option<String>,
}

impl MessageThreadSummary {
    pub fn last_message_preview(&self) -> &str {
        &self.last_message_preview
    }
}

impl std::fmt::Display for MessageType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MessageType::Text => write!(f, "text"),
            MessageType::Project { .. } => write!(f, "project"),
            MessageType::Request { .. } => write!(f, "request"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum MessageStatus {
    Draft,
    Sent,
    Received,
    Read,
    Deleted,
    Archived,
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
            MessageStatus::Archived => write!(f, "archived"),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn messages_new_defaults_and_reply_to() {
        let m = Message::new("a@x".into(), "b@y".into(), "hi".into());
        assert_eq!(m.from, "a@x");
        assert_eq!(m.to, "b@y");
        assert_eq!(m.status, MessageStatus::Draft);
        assert_eq!(m.sync_status, SyncStatus::Local);
        assert_eq!(m.display_subject(), "(No Subject)");

        let r = Message::reply_to(&m, "b@y".into(), "re".into());
        assert_eq!(r.to, "a@x");
        assert_eq!(r.thread_id, m.thread_id);
        assert_eq!(r.parent_id.as_deref(), Some(m.id.as_str()));
    }

    #[test]
    fn messages_display_impls() {
        let t1 = MessageType::Text;
        assert_eq!(t1.to_string(), "text");
        let t2 = MessageType::Project {
            project_name: "p".into(),
            submission_id: "s".into(),
            files_hash: None,
        };
        assert_eq!(t2.to_string(), "project");
        let t3 = MessageType::Request {
            request_type: "x".into(),
            params: None,
        };
        assert_eq!(t3.to_string(), "request");

        assert_eq!(MessageStatus::Draft.to_string(), "draft");
        assert_eq!(MessageStatus::Sent.to_string(), "sent");
        assert_eq!(SyncStatus::Synced.to_string(), "synced");
    }
}
