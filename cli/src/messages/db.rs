use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use rusqlite::{params, Connection, Row};
use std::path::Path;

use super::models::{Message, MessageStatus, SyncStatus};

pub struct MessageDb {
    conn: Connection,
}

impl MessageDb {
    pub fn new(db_path: &Path) -> Result<Self> {
        let conn = Connection::open(db_path)
            .with_context(|| format!("Failed to open database at {:?}", db_path))?;

        // Create tables if they don't exist
        conn.execute(
            "CREATE TABLE IF NOT EXISTS messages (
                id TEXT PRIMARY KEY,
                thread_id TEXT,
                parent_id TEXT,

                from_address TEXT NOT NULL,
                to_address TEXT NOT NULL,

                subject TEXT,
                body TEXT NOT NULL,

                status TEXT NOT NULL,
                sync_status TEXT NOT NULL,

                created_at TEXT NOT NULL,
                sent_at TEXT,
                received_at TEXT,
                read_at TEXT,

                rpc_request_id TEXT,
                rpc_ack_status INTEGER,
                rpc_ack_at TEXT
            )",
            [],
        )?;

        // Create indexes
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_thread_id ON messages(thread_id)",
            [],
        )?;
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_from_address ON messages(from_address)",
            [],
        )?;
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_to_address ON messages(to_address)",
            [],
        )?;
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_status ON messages(status)",
            [],
        )?;
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_created_at ON messages(created_at)",
            [],
        )?;

        Ok(Self { conn })
    }

    pub fn insert_message(&self, msg: &Message) -> Result<()> {
        self.conn.execute(
            "INSERT INTO messages (
                id, thread_id, parent_id,
                from_address, to_address,
                subject, body,
                status, sync_status,
                created_at, sent_at, received_at, read_at,
                rpc_request_id, rpc_ack_status, rpc_ack_at
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16)",
            params![
                msg.id,
                msg.thread_id,
                msg.parent_id,
                msg.from,
                msg.to,
                msg.subject,
                msg.body,
                msg.status.to_string(),
                msg.sync_status.to_string(),
                msg.created_at.to_rfc3339(),
                msg.sent_at.map(|dt| dt.to_rfc3339()),
                msg.received_at.map(|dt| dt.to_rfc3339()),
                msg.read_at.map(|dt| dt.to_rfc3339()),
                msg.rpc_request_id,
                msg.rpc_ack_status,
                msg.rpc_ack_at.map(|dt| dt.to_rfc3339()),
            ],
        )?;
        Ok(())
    }

    pub fn update_message(&self, msg: &Message) -> Result<()> {
        self.conn.execute(
            "UPDATE messages SET
                thread_id = ?2,
                parent_id = ?3,
                from_address = ?4,
                to_address = ?5,
                subject = ?6,
                body = ?7,
                status = ?8,
                sync_status = ?9,
                created_at = ?10,
                sent_at = ?11,
                received_at = ?12,
                read_at = ?13,
                rpc_request_id = ?14,
                rpc_ack_status = ?15,
                rpc_ack_at = ?16
            WHERE id = ?1",
            params![
                msg.id,
                msg.thread_id,
                msg.parent_id,
                msg.from,
                msg.to,
                msg.subject,
                msg.body,
                msg.status.to_string(),
                msg.sync_status.to_string(),
                msg.created_at.to_rfc3339(),
                msg.sent_at.map(|dt| dt.to_rfc3339()),
                msg.received_at.map(|dt| dt.to_rfc3339()),
                msg.read_at.map(|dt| dt.to_rfc3339()),
                msg.rpc_request_id,
                msg.rpc_ack_status,
                msg.rpc_ack_at.map(|dt| dt.to_rfc3339()),
            ],
        )?;
        Ok(())
    }

    pub fn get_message(&self, id: &str) -> Result<Option<Message>> {
        // First try exact match
        let mut stmt = self.conn.prepare("SELECT * FROM messages WHERE id = ?1")?;

        let mut rows = stmt.query(params![id])?;
        if let Some(row) = rows.next()? {
            return Ok(Some(Self::row_to_message(row)?));
        }

        // If not found and ID is short (partial), try prefix match
        if id.len() < 36 {
            // UUID is 36 chars with dashes
            let mut stmt = self
                .conn
                .prepare("SELECT * FROM messages WHERE id LIKE ?1 || '%'")?;

            let mut rows = stmt.query(params![id])?;
            let mut matches = Vec::new();
            while let Some(row) = rows.next()? {
                matches.push(Self::row_to_message(row)?);
            }

            if matches.len() == 1 {
                return Ok(Some(matches.into_iter().next().unwrap()));
            } else if matches.len() > 1 {
                return Err(anyhow::anyhow!(
                    "Multiple messages found with ID prefix '{}'. Please use a longer prefix.",
                    id
                ));
            }
        }

        Ok(None)
    }

    pub fn list_messages(&self, limit: Option<usize>) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!("SELECT * FROM messages WHERE status != 'deleted' ORDER BY created_at DESC LIMIT {}", limit)
        } else {
            "SELECT * FROM messages WHERE status != 'deleted' ORDER BY created_at DESC".to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn list_unread_messages(&self) -> Result<Vec<Message>> {
        let mut stmt = self
            .conn
            .prepare("SELECT * FROM messages WHERE status = 'received' ORDER BY created_at DESC")?;

        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn get_thread_messages(&self, thread_id: &str) -> Result<Vec<Message>> {
        let mut stmt = self.conn.prepare(
            "SELECT * FROM messages WHERE thread_id = ?1 AND status != 'deleted' ORDER BY created_at ASC"
        )?;

        let rows = stmt.query_map(params![thread_id], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn mark_as_read(&self, id: &str) -> Result<()> {
        self.conn.execute(
            "UPDATE messages SET status = 'read', read_at = ?2 WHERE id = ?1",
            params![id, Utc::now().to_rfc3339()],
        )?;
        Ok(())
    }

    pub fn mark_as_sent(&self, id: &str) -> Result<()> {
        self.conn.execute(
            "UPDATE messages SET status = 'sent', sent_at = ?2 WHERE id = ?1",
            params![id, Utc::now().to_rfc3339()],
        )?;
        Ok(())
    }

    pub fn delete_message(&self, id: &str) -> Result<()> {
        self.conn.execute(
            "UPDATE messages SET status = 'deleted' WHERE id = ?1",
            params![id],
        )?;
        Ok(())
    }

    pub fn list_deleted_messages(&self, limit: Option<usize>) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!(
                "SELECT * FROM messages WHERE status = 'deleted' ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT * FROM messages WHERE status = 'deleted' ORDER BY created_at DESC".to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn permanently_delete(&self, id: &str) -> Result<()> {
        self.conn
            .execute("DELETE FROM messages WHERE id = ?1", params![id])?;
        Ok(())
    }

    fn row_to_message(row: &Row) -> Result<Message> {
        Ok(Message {
            id: row.get(0)?,
            thread_id: row.get(1)?,
            parent_id: row.get(2)?,
            from: row.get(3)?,
            to: row.get(4)?,
            subject: row.get(5)?,
            body: row.get(6)?,
            status: Self::parse_status(&row.get::<_, String>(7)?)?,
            sync_status: Self::parse_sync_status(&row.get::<_, String>(8)?)?,
            created_at: DateTime::parse_from_rfc3339(&row.get::<_, String>(9)?)?
                .with_timezone(&Utc),
            sent_at: row
                .get::<_, Option<String>>(10)?
                .and_then(|s| DateTime::parse_from_rfc3339(&s).ok())
                .map(|dt| dt.with_timezone(&Utc)),
            received_at: row
                .get::<_, Option<String>>(11)?
                .and_then(|s| DateTime::parse_from_rfc3339(&s).ok())
                .map(|dt| dt.with_timezone(&Utc)),
            read_at: row
                .get::<_, Option<String>>(12)?
                .and_then(|s| DateTime::parse_from_rfc3339(&s).ok())
                .map(|dt| dt.with_timezone(&Utc)),
            rpc_request_id: row.get(13)?,
            rpc_ack_status: row.get(14)?,
            rpc_ack_at: row
                .get::<_, Option<String>>(15)?
                .and_then(|s| DateTime::parse_from_rfc3339(&s).ok())
                .map(|dt| dt.with_timezone(&Utc)),
        })
    }

    fn parse_status(s: &str) -> Result<MessageStatus> {
        match s {
            "draft" => Ok(MessageStatus::Draft),
            "sent" => Ok(MessageStatus::Sent),
            "received" => Ok(MessageStatus::Received),
            "read" => Ok(MessageStatus::Read),
            "deleted" => Ok(MessageStatus::Deleted),
            _ => anyhow::bail!("Unknown message status: {}", s),
        }
    }

    fn parse_sync_status(s: &str) -> Result<SyncStatus> {
        match s {
            "local" => Ok(SyncStatus::Local),
            "syncing" => Ok(SyncStatus::Syncing),
            "synced" => Ok(SyncStatus::Synced),
            "failed" => Ok(SyncStatus::Failed),
            _ => anyhow::bail!("Unknown sync status: {}", s),
        }
    }
}
