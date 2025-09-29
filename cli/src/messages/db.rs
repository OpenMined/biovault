use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use rusqlite::{params, Connection, Row};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

use super::models::{Message, MessageStatus, SyncStatus};

#[derive(Serialize, Deserialize)]
struct LockInfo {
    pid: u32,
    timestamp: DateTime<Utc>,
}

pub struct DatabaseLock {
    lock_path: PathBuf,
}

static LOCK_REGISTRY: OnceLock<Mutex<HashMap<PathBuf, usize>>> = OnceLock::new();

impl DatabaseLock {
    fn registry() -> &'static Mutex<HashMap<PathBuf, usize>> {
        LOCK_REGISTRY.get_or_init(|| Mutex::new(HashMap::new()))
    }

    fn increment_if_held(lock_path: &Path) -> bool {
        let registry = Self::registry();
        let mut map = registry.lock().unwrap();
        if let Some(count) = map.get_mut(lock_path) {
            *count += 1;
            return true;
        }
        false
    }

    fn increment_with_owner(lock_path: &Path) {
        let registry = Self::registry();
        let mut map = registry.lock().unwrap();
        let entry = map.entry(lock_path.to_path_buf()).or_insert(1);
        *entry += 1;
    }

    fn register_owner(lock_path: &Path) {
        let registry = Self::registry();
        let mut map = registry.lock().unwrap();
        let entry = map.entry(lock_path.to_path_buf()).or_insert(0);
        if *entry == 0 {
            *entry = 1;
        }
    }

    fn release_handle(lock_path: &Path) -> bool {
        let registry = Self::registry();
        let mut map = registry.lock().unwrap();
        if let Some(count) = map.get_mut(lock_path) {
            if *count > 1 {
                *count -= 1;
                return false;
            }
            map.remove(lock_path);
            return true;
        }
        true
    }

    fn acquire(db_path: &Path, timeout: Duration) -> Result<Self> {
        let lock_path = db_path.with_extension("lock");
        let start = Instant::now();
        let current_pid = std::process::id();

        loop {
            // If this process already holds the lock, allow re-entrant access
            if Self::increment_if_held(&lock_path) {
                return Ok(DatabaseLock {
                    lock_path: lock_path.clone(),
                });
            }

            // Check if lock file exists
            if lock_path.exists() {
                // Try to read existing lock info
                if let Ok(content) = fs::read_to_string(&lock_path) {
                    if let Ok(lock_info) = serde_json::from_str::<LockInfo>(&content) {
                        if lock_info.pid == current_pid {
                            Self::increment_with_owner(&lock_path);
                            return Ok(DatabaseLock {
                                lock_path: lock_path.clone(),
                            });
                        }
                        let age = Utc::now() - lock_info.timestamp;

                        // Check if process is still alive first (more important than age)
                        let mut should_remove = false;
                        let mut process_exists = true;

                        #[cfg(unix)]
                        {
                            unsafe {
                                let result = libc::kill(lock_info.pid as i32, 0);
                                if result != 0 {
                                    // Process doesn't exist, remove stale lock
                                    should_remove = true;
                                    process_exists = false;
                                }
                            }
                        }

                        // Also check if lock is stale (older than 30 seconds for faster recovery)
                        if age.num_seconds() > 30 {
                            should_remove = true;
                        }

                        if should_remove {
                            let _ = fs::remove_file(&lock_path);
                            // Continue to retry immediately after removing stale lock
                            continue;
                        } else if !process_exists {
                            // Process doesn't exist but we couldn't remove - try harder
                            let _ = fs::remove_file(&lock_path);
                            continue;
                        } else {
                            // Lock is valid - wait or timeout
                            if start.elapsed() > timeout {
                                return Err(anyhow::anyhow!(
                                    "Database lock timeout - held by active PID {} since {}",
                                    lock_info.pid,
                                    lock_info.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
                                ));
                            }
                            std::thread::sleep(Duration::from_millis(50));
                            continue;
                        }
                    } else {
                        // Invalid lock file format, remove it
                        let _ = fs::remove_file(&lock_path);
                    }
                }
            }

            // Try to create lock file atomically
            let lock_info = LockInfo {
                pid: current_pid,
                timestamp: Utc::now(),
            };

            // Try to create the lock file atomically using temp file + rename
            let temp_lock_path = lock_path.with_extension("tmp");
            match fs::write(&temp_lock_path, serde_json::to_string(&lock_info)?) {
                Ok(_) => {
                    // Try atomic rename - this is the true lock acquisition
                    match fs::rename(&temp_lock_path, &lock_path) {
                        Ok(_) => {
                            Self::register_owner(&lock_path);
                            return Ok(DatabaseLock {
                                lock_path: lock_path.clone(),
                            });
                        }
                        Err(_) => {
                            // Another process got the lock first, clean up our temp file
                            let _ = fs::remove_file(&temp_lock_path);
                            // Loop will check stale lock on next iteration
                        }
                    }
                }
                Err(_) => {
                    // Could not create temp file, wait and retry
                }
            }

            // If we get here, lock acquisition failed - wait and retry
            if start.elapsed() > timeout {
                // Try to provide better error info
                if let Ok(content) = fs::read_to_string(&lock_path) {
                    if let Ok(info) = serde_json::from_str::<LockInfo>(&content) {
                        return Err(anyhow::anyhow!(
                            "Database lock timeout - held by PID {} since {}",
                            info.pid,
                            info.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
                        ));
                    }
                }
                return Err(anyhow::anyhow!("Database lock timeout"));
            }
            std::thread::sleep(Duration::from_millis(50));
        }
    }
}

impl Drop for DatabaseLock {
    fn drop(&mut self) {
        if DatabaseLock::release_handle(&self.lock_path) {
            // Remove the lock file when the last holder in this process releases it
            let _ = fs::remove_file(&self.lock_path);
        }
    }
}

pub struct MessageDb {
    conn: Connection,
    _lock: Option<DatabaseLock>,
}

impl MessageDb {
    pub fn new(db_path: &Path) -> Result<Self> {
        Self::new_with_timeout(db_path, Duration::from_secs(30))
    }

    pub fn new_with_timeout(db_path: &Path, timeout: Duration) -> Result<Self> {
        // In tests, use a shorter timeout but still fail if can't acquire lock
        let lock = if cfg!(test) {
            Some(DatabaseLock::acquire(db_path, Duration::from_secs(10))?)
        } else {
            Some(DatabaseLock::acquire(db_path, timeout)?)
        };

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

                message_type TEXT DEFAULT 'text',
                metadata TEXT,

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

        // Migration: Add message_type column if it doesn't exist
        let has_message_type: bool = conn
            .prepare("SELECT COUNT(*) FROM pragma_table_info('messages') WHERE name='message_type'")
            .and_then(|mut stmt| stmt.query_row([], |row| row.get(0)))
            .unwrap_or(false);

        if !has_message_type {
            conn.execute(
                "ALTER TABLE messages ADD COLUMN message_type TEXT DEFAULT 'text'",
                [],
            )?;
        }

        // Migration: Add metadata column if it doesn't exist
        let has_metadata: bool = conn
            .prepare("SELECT COUNT(*) FROM pragma_table_info('messages') WHERE name='metadata'")
            .and_then(|mut stmt| stmt.query_row([], |row| row.get(0)))
            .unwrap_or(false);

        if !has_metadata {
            conn.execute("ALTER TABLE messages ADD COLUMN metadata TEXT", [])?;
        }

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
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_message_type ON messages(message_type)",
            [],
        )?;

        Ok(Self { conn, _lock: lock })
    }

    pub fn insert_message(&self, msg: &Message) -> Result<()> {
        let metadata_json = msg
            .metadata
            .as_ref()
            .map(serde_json::to_string)
            .transpose()?;

        self.conn.execute(
            "INSERT INTO messages (
                id, thread_id, parent_id,
                from_address, to_address,
                subject, body,
                message_type, metadata,
                status, sync_status,
                created_at, sent_at, received_at, read_at,
                rpc_request_id, rpc_ack_status, rpc_ack_at
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18)",
            params![
                msg.id,
                msg.thread_id,
                msg.parent_id,
                msg.from,
                msg.to,
                msg.subject,
                msg.body,
                msg.message_type.to_string(),
                metadata_json,
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
        let metadata_json = msg
            .metadata
            .as_ref()
            .map(serde_json::to_string)
            .transpose()?;

        self.conn.execute(
            "UPDATE messages SET
                thread_id = ?2,
                parent_id = ?3,
                from_address = ?4,
                to_address = ?5,
                subject = ?6,
                body = ?7,
                message_type = ?8,
                metadata = ?9,
                status = ?10,
                sync_status = ?11,
                created_at = ?12,
                sent_at = ?13,
                received_at = ?14,
                read_at = ?15,
                rpc_request_id = ?16,
                rpc_ack_status = ?17,
                rpc_ack_at = ?18
            WHERE id = ?1",
            params![
                msg.id,
                msg.thread_id,
                msg.parent_id,
                msg.from,
                msg.to,
                msg.subject,
                msg.body,
                msg.message_type.to_string(),
                metadata_json,
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
        let mut stmt = self.conn.prepare(
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE id = ?1",
        )?;

        let mut rows = stmt.query(params![id])?;
        if let Some(row) = rows.next()? {
            return Ok(Some(Self::row_to_message(row)?));
        }

        // If not found and ID is short (partial), try prefix match
        if id.len() < 36 {
            // UUID is 36 chars with dashes
            let mut stmt = self.conn.prepare(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE id LIKE ?1 || '%'",
            )?;

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
            format!(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                            status, sync_status, created_at, sent_at, received_at, read_at,
                            rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                     FROM messages WHERE status != 'deleted' ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE status != 'deleted' ORDER BY created_at DESC"
                .to_string()
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
        let mut stmt = self.conn.prepare(
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE status = 'received' ORDER BY created_at DESC",
        )?;

        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn get_thread_messages(&self, thread_id: &str) -> Result<Vec<Message>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE thread_id = ?1 AND status != 'deleted' ORDER BY created_at ASC",
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
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE status = 'deleted' ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE status = 'deleted' ORDER BY created_at DESC"
                .to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn list_inbox_messages(&self, limit: Option<usize>) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE status IN ('received', 'read') ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE status IN ('received', 'read') ORDER BY created_at DESC"
                .to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn list_sent_messages(&self, limit: Option<usize>) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE status IN ('sent', 'draft') ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE status IN ('sent', 'draft') ORDER BY created_at DESC"
                .to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map([], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn list_messages_by_type(
        &self,
        message_type: &str,
        limit: Option<usize>,
    ) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE message_type = ?1 AND status != 'deleted' ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE message_type = ?1 AND status != 'deleted' ORDER BY created_at DESC"
                .to_string()
        };

        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map(params![message_type], |row| Ok(Self::row_to_message(row)))?;

        let mut messages = Vec::new();
        for row in rows {
            messages.push(row??);
        }
        Ok(messages)
    }

    pub fn search_messages(&self, search: &str, limit: Option<usize>) -> Result<Vec<Message>> {
        let query = if let Some(limit) = limit {
            format!(
                "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                        status, sync_status, created_at, sent_at, received_at, read_at,
                        rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
                 FROM messages WHERE (body LIKE ?1 OR subject LIKE ?1) AND status != 'deleted' ORDER BY created_at DESC LIMIT {}",
                limit
            )
        } else {
            "SELECT id, thread_id, parent_id, from_address, to_address, subject, body,
                    status, sync_status, created_at, sent_at, received_at, read_at,
                    rpc_request_id, rpc_ack_status, rpc_ack_at, message_type, metadata
             FROM messages WHERE (body LIKE ?1 OR subject LIKE ?1) AND status != 'deleted' ORDER BY created_at DESC".to_string()
        };

        let search_pattern = format!("%{}%", search);
        let mut stmt = self.conn.prepare(&query)?;
        let rows = stmt.query_map(params![search_pattern], |row| Ok(Self::row_to_message(row)))?;

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
        // With our standardized SELECT columns, the order is:
        // 0: id, 1: thread_id, 2: parent_id, 3: from_address, 4: to_address,
        // 5: subject, 6: body, 7: status, 8: sync_status, 9: created_at,
        // 10: sent_at, 11: received_at, 12: read_at, 13: rpc_request_id,
        // 14: rpc_ack_status, 15: rpc_ack_at, 16: message_type, 17: metadata

        // Parse message_type
        let type_str: Option<String> = row.get(16)?;
        let message_type = type_str
            .map(|s| Self::parse_message_type(&s))
            .transpose()?
            .unwrap_or(crate::messages::MessageType::Text);

        // Parse metadata JSON
        let metadata_str: Option<String> = row.get(17)?;
        let metadata = metadata_str
            .and_then(|s| if s.is_empty() { None } else { Some(s) })
            .map(|s| serde_json::from_str(&s))
            .transpose()?;

        Ok(Message {
            id: row.get(0)?,
            thread_id: row.get(1)?,
            parent_id: row.get(2)?,
            from: row.get(3)?,
            to: row.get(4)?,
            subject: row.get(5)?,
            body: row.get(6)?,
            message_type,
            metadata,
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

    fn parse_message_type(s: &str) -> Result<crate::messages::MessageType> {
        use crate::messages::MessageType;
        match s {
            "text" => Ok(MessageType::Text),
            "project" => Ok(MessageType::Project {
                project_name: String::new(),
                submission_id: String::new(),
                files_hash: None,
            }),
            "request" => Ok(MessageType::Request {
                request_type: String::new(),
                params: None,
            }),
            _ => Ok(MessageType::Text), // Default to text for backwards compatibility
        }
    }

    fn parse_status(s: &str) -> Result<MessageStatus> {
        match s {
            "draft" => Ok(MessageStatus::Draft),
            "sent" => Ok(MessageStatus::Sent),
            "received" => Ok(MessageStatus::Received),
            "read" => Ok(MessageStatus::Read),
            "deleted" => Ok(MessageStatus::Deleted),
            "archived" => Ok(MessageStatus::Archived),
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

    /// Clean up stale database locks
    pub fn clean_stale_lock(db_path: &Path) -> Result<bool> {
        let lock_path = db_path.with_extension("lock");

        if !lock_path.exists() {
            return Ok(false);
        }

        // Try to read lock info
        match fs::read_to_string(&lock_path) {
            Ok(content) => {
                if let Ok(lock_info) = serde_json::from_str::<LockInfo>(&content) {
                    let age = Utc::now() - lock_info.timestamp;

                    // Check if lock is stale (older than 60 seconds)
                    if age.num_seconds() > 60 {
                        fs::remove_file(&lock_path)?;
                        println!(
                            "🧹 Removed stale lock from PID {} (age: {}s)",
                            lock_info.pid,
                            age.num_seconds()
                        );
                        return Ok(true);
                    }

                    // Check if process is still alive
                    #[cfg(unix)]
                    {
                        unsafe {
                            let result = libc::kill(lock_info.pid as i32, 0);
                            if result != 0 {
                                fs::remove_file(&lock_path)?;
                                println!("🧹 Removed stale lock from dead PID {}", lock_info.pid);
                                return Ok(true);
                            }
                        }
                    }

                    println!(
                        "🔒 Lock is held by active PID {} since {}",
                        lock_info.pid,
                        lock_info.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
                    );
                    Ok(false)
                } else {
                    // Invalid lock file format, remove it
                    fs::remove_file(&lock_path)?;
                    println!("🧹 Removed corrupted lock file");
                    Ok(true)
                }
            }
            Err(_) => {
                // Can't read lock file, remove it
                fs::remove_file(&lock_path)?;
                println!("🧹 Removed unreadable lock file");
                Ok(true)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::messages::models::{Message, MessageStatus};
    use tempfile::TempDir;

    fn new_db(tmp: &TempDir) -> MessageDb {
        let db_path = tmp.path().join("msgs.sqlite");
        MessageDb::new(&db_path).expect("create db")
    }

    #[test]
    fn reentrant_lock_allows_multiple_connections() {
        let tmp = TempDir::new().unwrap();
        let db_path = tmp.path().join("msgs.sqlite");
        let lock_path = db_path.with_extension("lock");

        let db1 = MessageDb::new(&db_path).expect("first db");
        assert!(lock_path.exists(), "lock file should be created");

        let db2 = MessageDb::new(&db_path).expect("second db");
        assert!(
            lock_path.exists(),
            "lock file should persist while both handles are active"
        );

        drop(db1);
        assert!(
            lock_path.exists(),
            "lock should remain until last handle drops"
        );

        drop(db2);
        assert!(
            !lock_path.exists(),
            "lock file should be removed after final handle"
        );
    }

    #[test]
    fn insert_get_and_list() {
        let tmp = TempDir::new().unwrap();
        let db = new_db(&tmp);

        let m = Message::new("a@x".into(), "b@y".into(), "hello world".into());
        db.insert_message(&m).unwrap();

        // Exact get
        let got = db.get_message(&m.id).unwrap().expect("found");
        assert_eq!(got.body, "hello world");

        // Prefix get
        let prefix = &m.id[..8];
        let got2 = db.get_message(prefix).unwrap().expect("prefix found");
        assert_eq!(got2.id, m.id);

        // list_messages (non-deleted)
        let all = db.list_messages(None).unwrap();
        assert_eq!(all.len(), 1);

        // list_sent_messages includes drafts and sent
        let sent = db.list_sent_messages(None).unwrap();
        assert_eq!(sent.len(), 1);
    }

    #[test]
    fn thread_and_status_updates() {
        let tmp = TempDir::new().unwrap();
        let db = new_db(&tmp);

        let m1 = Message::new("a@x".into(), "b@y".into(), "m1".into());
        db.insert_message(&m1).unwrap();

        // Reply creates same thread_id
        let m2 = Message::reply_to(&m1, "b@y".into(), "m2".into());
        db.insert_message(&m2).unwrap();

        let thread = db
            .get_thread_messages(m1.thread_id.as_ref().unwrap())
            .unwrap();
        assert_eq!(thread.len(), 2);
        assert_eq!(thread[0].id, m1.id);
        assert_eq!(thread[1].id, m2.id);

        // Mark as sent and read, and query inbox
        db.mark_as_sent(&m1.id).unwrap();

        // To make it show in inbox, set to received and update, then mark as read
        let mut m2_upd = m2.clone();
        m2_upd.status = MessageStatus::Received;
        db.update_message(&m2_upd).unwrap();
        let unread = db.list_unread_messages().unwrap();
        assert_eq!(unread.len(), 1);
        assert_eq!(unread[0].id, m2.id);

        db.mark_as_read(&m2.id).unwrap();
        let inbox = db.list_inbox_messages(None).unwrap();
        assert_eq!(inbox.len(), 1);
        assert_eq!(inbox[0].status, MessageStatus::Read);
    }

    #[test]
    fn delete_and_search_and_permanent_delete() {
        let tmp = TempDir::new().unwrap();
        let db = new_db(&tmp);

        let mut m = Message::new("a@x".into(), "b@y".into(), "find me body".into());
        m.subject = Some("about rust".into());
        db.insert_message(&m).unwrap();

        // Search by body
        assert_eq!(db.search_messages("find me", None).unwrap().len(), 1);
        // Search by subject
        assert_eq!(db.search_messages("rust", None).unwrap().len(), 1);
        assert_eq!(db.search_messages("nomatch", None).unwrap().len(), 0);

        // Soft delete
        db.delete_message(&m.id).unwrap();
        assert_eq!(db.list_deleted_messages(None).unwrap().len(), 1);
        assert_eq!(db.list_messages(None).unwrap().len(), 0);

        // Permanent delete
        db.permanently_delete(&m.id).unwrap();
        assert!(db.get_message(&m.id).unwrap().is_none());
    }

    #[test]
    fn list_messages_by_type_filters_correctly() {
        let tmp = TempDir::new().unwrap();
        let db = new_db(&tmp);

        let mut p = Message::new("a@x".into(), "b@y".into(), "p".into());
        p.message_type = crate::messages::MessageType::Project {
            project_name: "n".into(),
            submission_id: "s".into(),
            files_hash: None,
        };
        db.insert_message(&p).unwrap();

        let mut r = Message::new("a@x".into(), "b@y".into(), "r".into());
        r.message_type = crate::messages::MessageType::Request {
            request_type: "t".into(),
            params: None,
        };
        db.insert_message(&r).unwrap();

        let projects = db.list_messages_by_type("project", None).unwrap();
        assert_eq!(projects.len(), 1);
        assert_eq!(projects[0].id, p.id);

        let requests = db.list_messages_by_type("request", None).unwrap();
        assert_eq!(requests.len(), 1);
        assert_eq!(requests[0].id, r.id);
    }
}
