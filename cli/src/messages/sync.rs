use anyhow::{Context, Result};
use chrono::Utc;
use reqwest::blocking::Client as BlockingClient;
use reqwest::Url;
use serde::{Deserialize, Serialize};
use std::path::Path;
use tracing::instrument;

use crate::config::Config;
use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::endpoint::Endpoint;
use crate::syftbox::storage::WritePolicy;
use crate::syftbox::types::{RpcRequest, RpcResponse};
use syftbox_sdk::trace_context;
use syftbox_sdk::{has_syc_magic, parse_envelope};

use super::db::MessageDb;
use super::models::{
    DecryptionFailureReason, FailedMessage, Message, MessageStatus, MessageType, SyncStatus,
};

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

#[derive(Debug, Deserialize)]
struct SendHandlerResponse {
    request_id: String,
}

impl MessageSync {
    pub fn new(db_path: &Path, app: SyftBoxApp) -> Result<Self> {
        let db = MessageDb::new(db_path)?;
        Ok(Self { db, app })
    }

    fn syftbox_auth_enabled() -> bool {
        match std::env::var("SYFTBOX_AUTH_ENABLED") {
            Ok(v) => {
                let v = v.trim().to_lowercase();
                !(v.is_empty() || v == "0" || v == "false" || v == "no")
            }
            Err(_) => true,
        }
    }

    fn send_via_send_handler(&self, from: &str, to: &str, body: &[u8]) -> Result<String> {
        let server_url = std::env::var("SYFTBOX_SERVER_URL")
            .context("SYFTBOX_SERVER_URL is required to send messages via the server")?;
        let server_url = server_url.trim().trim_end_matches('/');
        if server_url.is_empty() {
            anyhow::bail!("SYFTBOX_SERVER_URL is empty");
        }

        // In dev mode (SYFTBOX_AUTH_ENABLED=0/false/no), SyftBox can accept unauthenticated send handler requests.
        // In auth-enabled mode, require a token so the send handler has auth.
        let auth_enabled = Self::syftbox_auth_enabled();
        let access_token = if auth_enabled {
            Some(
                Config::load()
                    .ok()
                    .and_then(|c| c.syftbox_credentials)
                    .and_then(|c| c.access_token)
                    .filter(|t| !t.trim().is_empty())
                    .ok_or_else(|| {
                        anyhow::anyhow!(
                            "SyftBox access token is missing. Re-authenticate in Settings → SyftBox."
                        )
                    })?,
            )
        } else {
            None
        };

        let target_syft_url = format!("syft://{}/app_data/biovault/rpc/message", to);
        let mut url = Url::parse(&format!("{}/api/v1/send/msg", server_url))
            .context("Failed to build send handler URL")?;
        url.query_pairs_mut()
            .append_pair("x-syft-url", &target_syft_url)
            .append_pair("x-syft-from", from);

        let client = BlockingClient::builder()
            .build()
            .context("Failed to build HTTP client")?;

        let mut req = client
            .post(url)
            .header(reqwest::header::CONTENT_TYPE, "application/json")
            .body(body.to_vec());
        if let Some(token) = access_token {
            req = req.bearer_auth(token);
        }
        let resp = req
            .send()
            .context("Failed to send message via send handler")?;

        let status = resp.status();
        let resp_text = resp
            .text()
            .context("Failed to read send handler response body")?;

        if status.is_client_error() || status.is_server_error() {
            anyhow::bail!("Send handler error {}: {}", status, resp_text);
        }

        let parsed: SendHandlerResponse = serde_json::from_str(&resp_text)
            .with_context(|| format!("Failed to parse send handler response: {}", resp_text))?;
        Ok(parsed.request_id)
    }

    /// Send a message via RPC
    #[instrument(skip(self), fields(component = "messaging", message_id = %message_id), err)]
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

        let payload_bytes = serde_json::to_vec(&payload)?;

        // Mark as sent (we'll set rpc_request_id after successful send)
        msg.sync_status = SyncStatus::Syncing;
        msg.status = MessageStatus::Sent;
        msg.sent_at = Some(Utc::now());
        self.db.update_message(&msg)?;

        // Default transport: write request directly into the recipient's RPC folder.
        // SyftBox sync will upload this write if ACL permits it, delivering the request to the peer.
        //
        // Optional override: set BIOVAULT_USE_SEND_HANDLER=1 to use the SyftBox HTTP send handler instead.
        let prefer_send_handler = matches!(
            std::env::var("BIOVAULT_USE_SEND_HANDLER")
                .ok()
                .as_deref()
                .map(str::trim)
                .map(str::to_lowercase)
                .as_deref(),
            Some("1" | "true" | "yes")
        );

        let mut mode = "filesystem";
        let send_result: Result<String> = if prefer_send_handler {
            mode = "send-handler";
            self.send_via_send_handler(&msg.from, &msg.to, &payload_bytes)
        } else {
            let mut rpc_request = RpcRequest::new(
                msg.from.clone(),
                format!("syft://{}/app_data/biovault/rpc/message", msg.to),
                "POST".to_string(),
                payload_bytes.clone(),
            );

            // Inject trace context for distributed tracing
            trace_context::inject_trace_context(&mut rpc_request);

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

            let request_path = recipient_rpc_dir.join(format!("{}.request", rpc_request.id));

            let write_policy = WritePolicy::Envelope {
                recipients: vec![msg.to.clone()],
                hint: Some(format!("message-{}", rpc_request.id)),
            };

            self.app.storage.write_json_with_shadow(
                &request_path,
                &rpc_request,
                write_policy,
                true,
            )?;

            Ok(rpc_request.id)
        };

        match send_result {
            Ok(request_id) => {
                msg.rpc_request_id = Some(request_id.clone());
                self.db.update_message(&msg)?;
                println!("📤 Message sent to {} ({})", msg.to, mode);
                Ok(())
            }
            Err(e) => {
                msg.sync_status = SyncStatus::Failed;
                let _ = self.db.update_message(&msg);
                Err(e)
            }
        }
    }

    /// Check for incoming messages
    #[instrument(skip(self), fields(component = "messaging"), err)]
    pub fn check_incoming(&self, no_cleanup: bool) -> Result<Vec<String>> {
        let (new_message_ids, _new_failed) = self.check_incoming_with_failures(no_cleanup)?;
        Ok(new_message_ids)
    }

    /// Check for incoming messages and also capture decryption failures
    /// Returns (new_message_ids, new_failed_count)
    pub fn check_incoming_with_failures(&self, no_cleanup: bool) -> Result<(Vec<String>, usize)> {
        let endpoint = Endpoint::new(&self.app, "/message")?;
        let (requests, failures) = endpoint.check_requests_with_failures()?;

        let mut new_message_ids = Vec::new();
        let mut new_failed_count = 0;

        // Process successful decryptions
        for (request_path, rpc_request) in requests {
            // Log trace context from incoming request for distributed tracing
            if let Some(traceparent) = rpc_request.traceparent() {
                if let Some((trace_id, parent_span_id, _flags)) =
                    trace_context::parse_traceparent(traceparent)
                {
                    tracing::info!(
                        trace_id = %trace_id,
                        parent_span_id = %parent_span_id,
                        "Processing request with trace context"
                    );
                }
            }

            let body_bytes = match rpc_request.decode_body() {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("Failed to decode request body: {}", e);
                    continue;
                }
            };
            let payload: MessagePayload = match serde_json::from_slice(&body_bytes) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("Failed to parse message payload: {}", e);
                    continue;
                }
            };

            let mut msg = Message::new(payload.from.clone(), self.app.email.clone(), payload.body);

            msg.id = payload.message_id;
            msg.thread_id = payload.thread_id;
            msg.parent_id = payload.parent_id;
            msg.subject = payload.subject;
            msg.message_type = match payload.message_type.as_str() {
                "module" => MessageType::Module {
                    module_name: String::new(),
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
            if let Ok(dt) = chrono::DateTime::parse_from_rfc3339(&payload.created_at) {
                msg.created_at = dt.with_timezone(&Utc);
            }

            if self.db.get_message(&msg.id)?.is_none() {
                self.db.insert_message(&msg)?;
                new_message_ids.push(msg.id.clone());
            }

            // If we successfully decrypted, remove any previous failed record for this RPC ID
            let _ = self.db.delete_failed_message_by_rpc_id(&rpc_request.id);

            // Send ACK and clean up
            let ack_response = RpcResponse::new(
                &rpc_request,
                self.app.email.clone(),
                200,
                b"Message received".to_vec(),
            );
            let _ = endpoint.send_response(&request_path, &rpc_request, &ack_response, no_cleanup);
        }

        // Process failures - extract envelope metadata and store
        for (request_path, error_msg, raw_bytes) in failures {
            // Extract RPC request ID from filename (UUID.request -> UUID)
            let rpc_request_id = request_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();

            // Skip if we already have a failed record for this RPC ID
            if self
                .db
                .get_failed_message_by_rpc_id(&rpc_request_id)?
                .is_some()
            {
                continue;
            }

            eprintln!(
                "⚠️  Failed to decrypt/parse incoming message (rpc_id={}): {}",
                rpc_request_id, error_msg
            );

            // Try to extract envelope metadata
            let (
                sender_identity,
                sender_fingerprint,
                recipient_fingerprint,
                filename_hint,
                failure_reason,
            ) = self.extract_envelope_metadata(&raw_bytes, &error_msg);

            let mut failed = FailedMessage::new(
                request_path.to_string_lossy().to_string(),
                rpc_request_id,
                sender_identity,
                sender_fingerprint,
                failure_reason,
                error_msg,
            );
            failed.recipient_identity = Some(self.app.email.clone());
            failed.recipient_fingerprint = recipient_fingerprint;
            failed.filename_hint = filename_hint;

            self.db.insert_failed_message(&failed)?;
            new_failed_count += 1;
        }

        Ok((new_message_ids, new_failed_count))
    }

    /// Extract envelope metadata from raw bytes for failure reporting
    fn extract_envelope_metadata(
        &self,
        raw_bytes: &[u8],
        error_msg: &str,
    ) -> (
        String,
        String,
        Option<String>,
        Option<String>,
        DecryptionFailureReason,
    ) {
        // Default values if we can't parse the envelope
        let mut sender_identity = "unknown".to_string();
        let mut sender_fingerprint = "unknown".to_string();
        let mut recipient_fingerprint = None;
        let mut filename_hint = None;

        // Determine failure reason from error message
        let failure_reason = if error_msg.contains("sender bundle not cached")
            || error_msg.contains("no cached bundle")
        {
            DecryptionFailureReason::SenderBundleNotCached
        } else if error_msg.contains("recipient")
            || error_msg.contains("not addressed")
            || error_msg.contains("wrong recipient")
        {
            DecryptionFailureReason::WrongRecipient
        } else if error_msg.contains("decrypt") || error_msg.contains("cipher") {
            DecryptionFailureReason::DecryptionFailed
        } else if error_msg.contains("fingerprint") || error_msg.contains("key mismatch") {
            DecryptionFailureReason::RecipientKeyMismatch
        } else {
            DecryptionFailureReason::Other(error_msg.to_string())
        };

        // Try to parse the envelope
        if has_syc_magic(raw_bytes) {
            if let Ok(parsed) = parse_envelope(raw_bytes) {
                sender_identity = parsed.prelude.sender.identity.clone();
                sender_fingerprint = parsed.prelude.sender.ik_fingerprint.clone();

                // Get recipient info from first wrapping
                if let Some(wrapping) = parsed.prelude.wrappings.first() {
                    if let Some(ref ri) = wrapping.recipient_identity {
                        // We could use this if needed
                        let _ = ri;
                    }
                }

                // Get recipient fingerprint from first recipient
                if let Some(recipient) = parsed.prelude.recipients.first() {
                    if let Some(ref spk_fp) = recipient.spk_fingerprint {
                        recipient_fingerprint = Some(spk_fp.clone());
                    }
                }

                // Get filename hint if available
                if let Some(ref meta) = parsed.prelude.public_meta {
                    filename_hint = meta.filename_hint.clone();
                }
            }
        }

        (
            sender_identity,
            sender_fingerprint,
            recipient_fingerprint,
            filename_hint,
            failure_reason,
        )
    }

    /// Check for ACK responses to sent messages
    #[instrument(skip(self), fields(component = "messaging"), err)]
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
    #[instrument(skip(self), fields(component = "messaging"), err)]
    pub fn sync(&self, no_cleanup: bool) -> Result<()> {
        let (new_message_ids, new_failed_count) = self.check_incoming_with_failures(no_cleanup)?;
        if !new_message_ids.is_empty() {
            println!("📬 {} new message(s) received", new_message_ids.len());
        }
        if new_failed_count > 0 {
            eprintln!(
                "⚠️  Recorded {} message decryption failure(s) (see failed_messages table).",
                new_failed_count
            );
        }
        self.check_acks(no_cleanup)?;
        Ok(())
    }

    /// Silent sync - returns results without printing
    /// Returns (new_message_ids, new_message_count)
    pub fn sync_quiet(&self) -> Result<(Vec<String>, usize)> {
        let (new_messages, _count, _failed_count) = self.sync_quiet_with_failures()?;
        Ok((new_messages.clone(), new_messages.len()))
    }

    /// Silent sync with failure tracking - returns (new_message_ids, new_message_count, failed_count)
    pub fn sync_quiet_with_failures(&self) -> Result<(Vec<String>, usize, usize)> {
        let (new_messages, failed_count) = self.check_incoming_with_failures(false)?;
        self.check_acks(false)?;
        Ok((new_messages.clone(), new_messages.len(), failed_count))
    }

    /// List failed messages from the database
    pub fn list_failed_messages(&self, include_dismissed: bool) -> Result<Vec<FailedMessage>> {
        self.db.list_failed_messages(include_dismissed)
    }

    /// Count non-dismissed failed messages
    pub fn count_failed_messages(&self) -> Result<usize> {
        self.db.count_failed_messages()
    }

    /// Dismiss a failed message
    pub fn dismiss_failed_message(&self, id: &str) -> Result<bool> {
        self.db.dismiss_failed_message(id)
    }

    /// Delete a failed message
    pub fn delete_failed_message(&self, id: &str) -> Result<bool> {
        self.db.delete_failed_message(id)
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

        // Get the request ID that was assigned during send
        let req_id = ms_sender
            .db
            .get_message(&m.id)
            .unwrap()
            .unwrap()
            .rpc_request_id
            .unwrap();

        // Create a failure response file in the recipient's endpoint folder (where the request was sent).
        // check_responses scans other people's folders to find responses to our requests.
        let recipient_endpoint_path = app_sender
            .data_dir
            .join("datasites")
            .join(&m.to)
            .join("app_data")
            .join("biovault")
            .join("rpc")
            .join("message");

        // Build a dummy RpcResponse with non-200 code
        let dummy_req = RpcRequest::new(
            "x".into(),
            app_sender.build_syft_url("/message"),
            "POST".into(),
            b"{}".to_vec(),
        );
        let mut resp = RpcResponse::new(&dummy_req, m.to.clone(), 500, b"err".to_vec());
        // Force response id to match the request id we need to ack
        resp.id = req_id.clone();
        let resp_path = recipient_endpoint_path.join(format!("{}.response", req_id));
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
