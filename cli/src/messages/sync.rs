use anyhow::{Context, Result};
use chrono::Utc;
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::endpoint::Endpoint;
use crate::syftbox::storage::WritePolicy;
use crate::syftbox::types::{RpcRequest, RpcResponse};

use super::db::MessageDb;
use super::models::{Message, MessageStatus, MessageType, SyncStatus};

#[derive(Debug, Serialize, Deserialize)]
struct MessagePayload {
    message_id: String,
    thread_id: Option<String>,
    parent_id: Option<String>,
    from: String,
    to: String,
    subject: Option<String>,
    body: String,
    message_type: String,
    metadata: Option<serde_json::Value>,
    created_at: String,
}

pub struct MessageSync {
    db: MessageDb,
    app: SyftBoxApp,
}

impl MessageSync {
    pub fn new(db_path: &Path, app: SyftBoxApp) -> Result<Self> {
        let db = MessageDb::new(db_path)?;
        Ok(Self { db, app })
    }

    /// Send a message via RPC
    pub fn send_message(&self, message_id: &str) -> Result<()> {
        let mut msg = self
            .db
            .get_message(message_id)?
            .ok_or_else(|| anyhow::anyhow!("Message not found: {}", message_id))?;

        // Create RPC request with message payload
        let payload = MessagePayload {
            message_id: msg.id.clone(),
            thread_id: msg.thread_id.clone(),
            parent_id: msg.parent_id.clone(),
            from: msg.from.clone(),
            to: msg.to.clone(),
            subject: msg.subject.clone(),
            body: msg.body.clone(),
            message_type: msg.message_type.to_string(),
            metadata: msg.metadata.clone(),
            created_at: msg.created_at.to_rfc3339(),
        };

        let rpc_request = RpcRequest::new(
            msg.from.clone(),
            format!(
                "syft://{}@openmined.org/app_data/biovault/rpc/message",
                msg.to
            ),
            "POST".to_string(),
            serde_json::to_vec(&payload)?,
        );

        // Update message with RPC tracking
        msg.rpc_request_id = Some(rpc_request.id.clone());
        msg.sync_status = SyncStatus::Syncing;
        msg.status = MessageStatus::Sent;
        msg.sent_at = Some(Utc::now());
        self.db.update_message(&msg)?;

        // Write request directly to recipient's datasite RPC folder
        let recipient_rpc_dir = self
            .app
            .data_dir
            .join("datasites")
            .join(&msg.to)
            .join("app_data")
            .join("biovault")
            .join("rpc")
            .join("message");

        self.app
            .storage
            .ensure_dir(&recipient_rpc_dir)
            .with_context(|| format!("Failed to prepare RPC dir {:?}", recipient_rpc_dir))?;

        let request_filename = format!("{}.request", rpc_request.id);
        let request_path = recipient_rpc_dir.join(request_filename);

        let write_policy = WritePolicy::Envelope {
            recipients: vec![msg.to.clone()],
            hint: Some(format!("message-{}", rpc_request.id)),
        };

        self.app
            .storage
            .write_json_with_shadow(&request_path, &rpc_request, write_policy, true)?;

        println!("📤 Message sent to {}", msg.to);
        println!("Request written to: {:?}", request_path);
        Ok(())
    }

    /// Check for incoming messages
    pub fn check_incoming(&self, no_cleanup: bool) -> Result<Vec<String>> {
        let endpoint = Endpoint::new(&self.app, "/message")?;
        let requests = endpoint.check_requests()?;

        let mut new_message_ids = Vec::new();

        for (request_path, rpc_request) in requests {
            // Parse the message payload
            let body_bytes = rpc_request
                .decode_body()
                .context("Failed to decode request body")?;
            let payload: MessagePayload =
                serde_json::from_slice(&body_bytes).context("Failed to parse message payload")?;

            // Create local message from received payload
            let mut msg = Message::new(
                payload.from.clone(),
                self.app.email.clone(), // We are the recipient
                payload.body,
            );

            msg.id = payload.message_id;
            msg.thread_id = payload.thread_id;
            msg.parent_id = payload.parent_id;
            msg.subject = payload.subject;
            // Parse message type from string
            msg.message_type = match payload.message_type.as_str() {
                "project" => MessageType::Project {
                    project_name: String::new(),
                    submission_id: String::new(),
                    files_hash: None,
                },
                "request" => MessageType::Request {
                    request_type: String::new(),
                    params: None,
                },
                _ => MessageType::Text,
            };
            msg.metadata = payload.metadata;
            msg.status = MessageStatus::Received;
            msg.sync_status = SyncStatus::Synced;
            msg.received_at = Some(Utc::now());
            msg.created_at =
                chrono::DateTime::parse_from_rfc3339(&payload.created_at)?.with_timezone(&Utc);

            // Check if message already exists before inserting
            if self.db.get_message(&msg.id)?.is_none() {
                // Save to database
                self.db.insert_message(&msg)?;
                new_message_ids.push(msg.id.clone());
            }

            // Process status-update requests to update original message metadata
            // Status update: detect via metadata.status_update and update original message metadata
            if let Some(meta) = &msg.metadata {
                if let Some(update) = meta.get("status_update") {
                    if let (Some(orig_id), Some(status)) = (
                        update.get("message_id").and_then(|v| v.as_str()),
                        update.get("status").and_then(|v| v.as_str()),
                    ) {
                        if let Some(mut orig) = self.db.get_message(orig_id)? {
                            // Work on a cloned metadata value to avoid moving out of 'orig'
                            let mut md = orig.metadata.clone().unwrap_or(serde_json::json!({}));
                            // Ensure object shape
                            if !md.is_object() {
                                md = serde_json::json!({});
                            }
                            if let Some(obj) = md.as_object_mut() {
                                obj.insert(
                                    "remote_status".to_string(),
                                    serde_json::Value::String(status.to_string()),
                                );
                                if let Some(reason) = update.get("reason").and_then(|v| v.as_str())
                                {
                                    obj.insert(
                                        "remote_reason".to_string(),
                                        serde_json::Value::String(reason.to_string()),
                                    );
                                }
                                if let Some(results_path) =
                                    meta.get("results_path").and_then(|v| v.as_str())
                                {
                                    obj.insert(
                                        "results_path".to_string(),
                                        serde_json::Value::String(results_path.to_string()),
                                    );
                                }
                            }
                            orig.metadata = Some(md);
                            self.db.update_message(&orig)?;
                        }
                    }
                }
            }

            // Send ACK response
            let ack_response = RpcResponse::new(
                &rpc_request,
                self.app.email.clone(),
                200,
                b"Message received".to_vec(),
            );
            endpoint.send_response(&request_path, &rpc_request, &ack_response, no_cleanup)?;

            // Silent - will be shown when listing messages
        }

        Ok(new_message_ids)
    }

    /// Check for ACK responses to sent messages
    pub fn check_acks(&self, no_cleanup: bool) -> Result<()> {
        let endpoint = Endpoint::new(&self.app, "/message")?;
        let responses = endpoint.check_responses()?;

        for (response_path, rpc_response) in responses {
            // The response ID matches the request ID
            let request_id = &rpc_response.id;

            // Find message by RPC request ID
            let messages = self.db.list_messages(None)?;
            for mut msg in messages {
                if msg.rpc_request_id.as_ref() == Some(request_id) {
                    // Update message with ACK info
                    msg.rpc_ack_status = Some(rpc_response.status_code as i32);
                    msg.rpc_ack_at = Some(Utc::now());
                    msg.sync_status = if rpc_response.status_code == 200 {
                        SyncStatus::Synced
                    } else {
                        SyncStatus::Failed
                    };
                    self.db.update_message(&msg)?;

                    // Silent - confirmation tracked in database
                    break;
                }
            }

            // Clean up response file (unless --no-cleanup mode)
            if !no_cleanup {
                endpoint.cleanup_response(&response_path)?;
            }
        }

        Ok(())
    }

    /// Full sync cycle: check incoming, send pending, check acks (verbose)
    pub fn sync(&self, no_cleanup: bool) -> Result<()> {
        // Check for new incoming messages
        let new_messages = self.check_incoming(no_cleanup)?;
        if !new_messages.is_empty() {
            println!("📬 {} new message(s) received", new_messages.len());
        }

        // Check for ACK responses
        self.check_acks(no_cleanup)?;

        Ok(())
    }

    /// Silent sync - returns results without printing
    pub fn sync_quiet(&self) -> Result<(Vec<String>, usize)> {
        // Check for new incoming messages (always cleanup in quiet mode - used internally)
        let new_messages = self.check_incoming(false)?;

        // Check for ACK responses
        self.check_acks(false)?;

        Ok((new_messages.clone(), new_messages.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::messages::models::Message;
    use crate::syftbox::storage::SyftBoxStorage;
    use tempfile::TempDir;

    // Helper to list response files under an endpoint path
    fn list_response_files(
        storage: &SyftBoxStorage,
        dir: &std::path::Path,
    ) -> Vec<std::path::PathBuf> {
        storage
            .list_dir(dir)
            .unwrap_or_default()
            .into_iter()
            .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("response"))
            .collect()
    }

    #[test]
    fn send_receive_and_ack_flow_via_fs() {
        let tmp = TempDir::new().unwrap();
        let data_dir = tmp.path();

        // Both parties share the same syftbox data_dir (as in a sync-mounted folder), with different datasites
        let app_sender = SyftBoxApp::new(data_dir, "alice@example.com", "biovault").unwrap();
        let app_recipient = SyftBoxApp::new(data_dir, "bob@example.com", "biovault").unwrap();

        let db_sender = tmp.path().join("sender.sqlite");
        let db_recipient = tmp.path().join("recipient.sqlite");
        let ms_sender = MessageSync::new(&db_sender, app_sender.clone()).unwrap();
        let ms_recipient = MessageSync::new(&db_recipient, app_recipient.clone()).unwrap();

        // Insert a draft message in sender DB targeting bob
        let m = Message::new(
            "alice@example.com".into(),
            "bob@example.com".into(),
            "hello".into(),
        );
        ms_sender.db.insert_message(&m).unwrap();

        // Send the message: should create request file under bob's endpoint and update sender DB fields
        ms_sender.send_message(&m.id).unwrap();

        // Recipient checks incoming; should ingest and write an ACK response in bob's endpoint dir
        let new_msgs = ms_recipient.check_incoming(false).unwrap();
        assert_eq!(new_msgs.len(), 1);

        // Simulate syft sync by copying the response from bob's endpoint to alice's endpoint
        let recipient_ep = app_recipient.endpoint_path("/message");
        let sender_ep = app_sender.endpoint_path("/message");
        app_sender.storage.ensure_dir(&sender_ep).unwrap();
        for resp in list_response_files(&app_recipient.storage, &recipient_ep) {
            let file_name = resp.file_name().unwrap();
            let dest = sender_ep.join(file_name);
            app_sender.storage.copy_raw_file(&resp, &dest).unwrap();
        }

        // Sender checks ACKs and should update message sync status
        ms_sender.check_acks(false).unwrap();
    }

    #[test]
    fn ack_failure_sets_failed_status() {
        let tmp = TempDir::new().unwrap();
        let data_dir = tmp.path();
        let app_sender = SyftBoxApp::new(data_dir, "alice@example.com", "biovault").unwrap();
        let db_sender = tmp.path().join("sender.sqlite");
        let ms_sender = MessageSync::new(&db_sender, app_sender.clone()).unwrap();

        // Insert a message and mark it as sent with an RPC request id
        let m = Message::new(
            "alice@example.com".into(),
            "bob@example.com".into(),
            "hello".into(),
        );
        ms_sender.db.insert_message(&m).unwrap();
        ms_sender.send_message(&m.id).unwrap();

        // Create a failure response file for the same request id
        let endpoint = Endpoint::new(&app_sender, "/message").unwrap();
        let req_id = ms_sender
            .db
            .get_message(&m.id)
            .unwrap()
            .unwrap()
            .rpc_request_id
            .unwrap();
        // Build a dummy RpcResponse with non-200 code
        let dummy_req = RpcRequest::new(
            "x".into(),
            app_sender.build_syft_url("/message"),
            "POST".into(),
            b"{}".to_vec(),
        );
        let mut resp = RpcResponse::new(&dummy_req, "bob@example.com".into(), 500, b"err".to_vec());
        // Force response id to match the request id we need to ack
        resp.id = req_id.clone();
        let resp_path = endpoint.path.join(format!("{}.response", req_id));
        app_sender
            .storage
            .write_plaintext_file(
                &resp_path,
                serde_json::to_string_pretty(&resp).unwrap().as_bytes(),
                true,
            )
            .unwrap();

        // Process ACKs and verify status is Failed
        ms_sender.check_acks(false).unwrap();
        let updated = ms_sender.db.get_message(&m.id).unwrap().unwrap();
        assert_eq!(updated.sync_status, SyncStatus::Failed);
        assert!(updated.rpc_ack_at.is_some());
    }
}
