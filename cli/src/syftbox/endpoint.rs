use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::storage::WritePolicy;
use crate::syftbox::types::{RpcRequest, RpcResponse};
use anyhow::{anyhow, bail, Context, Result};
use std::path::{Path, PathBuf};

/// Represents an RPC endpoint
pub struct Endpoint {
    pub name: String,
    pub path: PathBuf,
    #[allow(dead_code)]
    app: SyftBoxApp,
}

impl Endpoint {
    /// Create a new endpoint for a SyftBox app
    pub fn new(app: &SyftBoxApp, name: &str) -> Result<Self> {
        let path = app.register_endpoint(name)?;

        Ok(Self {
            name: name.to_string(),
            path,
            app: app.clone(),
        })
    }

    /// Check for request files in this endpoint (your own folder only)
    /// Requests are written by others TO you in datasites/your-email/app_data/biovault/rpc/message/
    pub fn check_requests(&self) -> Result<Vec<(PathBuf, RpcRequest)>> {
        // Only check your own endpoint folder - this is where others write requests TO you
        self.check_requests_in_folder(&self.path)
    }

    /// Helper to check for requests in a specific folder
    fn check_requests_in_folder(&self, folder: &Path) -> Result<Vec<(PathBuf, RpcRequest)>> {
        let mut requests = Vec::new();

        // Return empty if folder doesn't exist
        if !folder.exists() {
            return Ok(requests);
        }

        for entry in std::fs::read_dir(folder)? {
            let entry = entry?;
            let path = entry.path();

            // Check if it's a .request file
            if path.extension().and_then(|s| s.to_str()) == Some("request") {
                match self.read_request(&path) {
                    Ok(request) => {
                        requests.push((path, request));
                    }
                    Err(e) => {
                        eprintln!("Failed to read request file {:?}: {}", path, e);
                    }
                }
            }
        }

        Ok(requests)
    }

    /// Read a request file
    fn read_request(&self, path: &Path) -> Result<RpcRequest> {
        self.app
            .storage
            .read_json_with_shadow(path)
            .with_context(|| format!("Failed to parse request JSON from: {:?}", path))
    }

    /// Send a response for a request
    pub fn send_response(
        &self,
        request_path: &Path,
        request: &RpcRequest,
        response: &RpcResponse,
        no_cleanup: bool,
    ) -> Result<()> {
        if request_path.is_dir() {
            bail!(
                "request path {:?} is a directory, expected file",
                request_path
            );
        }

        // Get the request ID from the filename (UUID.request -> UUID)
        let request_filename = request_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("Invalid request filename"))?;

        // Create response file path next to the request (UUID.response in same directory)
        // This ensures the sender can find the response in their own app_data folder
        let response_filename = format!("{}.response", request_filename);
        let request_dir = request_path
            .parent()
            .ok_or_else(|| anyhow!("Request path has no parent"))?;
        let response_path = request_dir.join(response_filename);

        eprintln!("DEBUG send_response:");
        eprintln!("  request_path: {:?}", request_path);
        eprintln!("  response_path: {:?}", response_path);

        let write_policy = WritePolicy::Envelope {
            recipients: vec![request.sender.clone()],
            hint: Some(format!("rpc-response-{}", response.id)),
        };

        self.app
            .storage
            .write_json_with_shadow(&response_path, response, write_policy, true)?;

        // Delete the request file (unless --no-cleanup mode)
        if !no_cleanup {
            self.app.storage.remove_path(request_path)?;
        }

        println!("Sent response to: {:?}", response_path);
        Ok(())
    }

    /// Create a request file (for sending requests to others)
    pub fn create_request(&self, request: &RpcRequest) -> Result<PathBuf> {
        let request_filename = format!("{}.request", request.id);
        let request_path = self.path.join(request_filename);

        let write_policy = match extract_email_from_syft_url(&request.url) {
            Some(email) => WritePolicy::Envelope {
                recipients: vec![email],
                hint: Some(format!("rpc-request-{}", request.id)),
            },
            None => WritePolicy::Plaintext,
        };

        self.app
            .storage
            .write_json_with_shadow(&request_path, request, write_policy, true)?;

        println!("Created request: {:?}", request_path);
        Ok(request_path)
    }

    /// Check for response files (when we sent a request and are waiting for response)
    /// Responses are in OTHER identity folders (where you wrote requests), NOT your own folder
    /// You write requests to datasites/their-email/, they write responses back there
    pub fn check_responses(&self) -> Result<Vec<(PathBuf, RpcResponse)>> {
        let mut responses = Vec::new();

        // Get the datasites root directory
        let datasites_root = self.app.data_dir.join("datasites");
        if !datasites_root.exists() {
            return Ok(responses);
        }

        // Canonicalize paths to handle symlinks (e.g., /var -> /private/var on macOS)
        let canonical_datasites_root = datasites_root
            .canonicalize()
            .unwrap_or_else(|_| datasites_root.clone());
        let canonical_path = self
            .path
            .canonicalize()
            .unwrap_or_else(|_| self.path.clone());

        // Get the endpoint-relative path (e.g., "app_data/biovault/rpc/message")
        let endpoint_suffix = canonical_path
            .strip_prefix(&canonical_datasites_root)
            .ok()
            .and_then(|p| p.strip_prefix(self.app.email.as_str()).ok());

        if endpoint_suffix.is_none() {
            // Path structure unexpected - can't determine suffix
            return Ok(responses);
        }

        let endpoint_suffix = endpoint_suffix.unwrap();

        // Scan all identity folders
        // Responses are where you wrote requests (in other people's folders)
        // Also check your own folder for test/development scenarios
        for entry in std::fs::read_dir(&datasites_root)? {
            let entry = entry?;
            let identity_path = entry.path();

            if !identity_path.is_dir() {
                continue;
            }

            // Check the endpoint folder for this identity
            let endpoint_folder = identity_path.join(endpoint_suffix);
            if endpoint_folder.exists() {
                responses.extend(self.check_responses_in_folder(&endpoint_folder)?);
            }
        }

        Ok(responses)
    }

    /// Helper to check for responses in a specific folder
    fn check_responses_in_folder(&self, folder: &Path) -> Result<Vec<(PathBuf, RpcResponse)>> {
        let mut responses = Vec::new();

        // Return empty if folder doesn't exist
        if !folder.exists() {
            return Ok(responses);
        }

        for entry in std::fs::read_dir(folder)? {
            let entry = entry?;
            let path = entry.path();

            // Check if it's a .response file
            if path.extension().and_then(|s| s.to_str()) == Some("response") {
                match self.read_response(&path) {
                    Ok(response) => {
                        responses.push((path, response));
                    }
                    Err(e) => {
                        eprintln!("Failed to read response file {:?}: {}", path, e);
                    }
                }
            }
        }

        Ok(responses)
    }

    /// Read a response file
    fn read_response(&self, path: &Path) -> Result<RpcResponse> {
        self.app
            .storage
            .read_json_with_shadow(path)
            .with_context(|| format!("Failed to parse response JSON from: {:?}", path))
    }

    /// Clean up a response file after processing
    pub fn cleanup_response(&self, response_path: &Path) -> Result<()> {
        self.app.storage.remove_path(response_path)?;
        Ok(())
    }
}

fn extract_email_from_syft_url(url: &str) -> Option<String> {
    let trimmed = url.strip_prefix("syft://")?;
    trimmed.split('/').next().map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::syftbox::types::RpcRequest;
    use base64::Engine as _;
    use tempfile::TempDir;

    #[test]
    fn test_endpoint_request_response_flow() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;

        let endpoint = Endpoint::new(&app, "/message")?;

        // Create a request
        let request = RpcRequest::new(
            "sender@example.com".to_string(),
            app.build_syft_url("/message"),
            "POST".to_string(),
            b"Hello".to_vec(),
        );

        let request_path = endpoint.create_request(&request)?;
        assert!(request_path.exists());

        // Check for requests
        let requests = endpoint.check_requests()?;
        assert_eq!(requests.len(), 1);

        // Send a response
        let (req_path, req) = &requests[0];
        let response = RpcResponse::ok_json(
            req,
            "test@example.com".to_string(),
            &serde_json::json!({"reply": "Hi!"}),
        )
        .expect("ok_json");

        endpoint.send_response(req_path, req, &response, false)?;

        // Verify request was deleted and response created
        assert!(!request_path.exists());

        let responses = endpoint.check_responses()?;
        assert_eq!(responses.len(), 1);

        Ok(())
    }

    #[test]
    fn malformed_request_and_response_cleanup_and_error_paths() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;
        let endpoint = Endpoint::new(&app, "/message")?;

        // 1) Malformed request file should be skipped without panic
        let bad_req_path = endpoint.path.join("bad.request");
        std::fs::write(&bad_req_path, "{not json}")?;
        let reqs = endpoint.check_requests()?;
        assert!(reqs.is_empty());

        // Also a directory named *.request triggers read error path
        let dir_req_path = endpoint.path.join("dir.request");
        std::fs::create_dir_all(&dir_req_path)?;
        let reqs2 = endpoint.check_requests()?;
        assert!(reqs2.is_empty());

        // 2) Create a valid response file and then cleanup
        let req = RpcRequest::new(
            "sender@example.com".to_string(),
            app.build_syft_url("/message"),
            "POST".to_string(),
            b"Hello".to_vec(),
        );
        // Write a minimal response JSON with required fields
        let resp_path = endpoint.path.join(format!("{}.response", req.id));
        std::fs::write(
            &resp_path,
            serde_json::json!({
                "id": req.id,
                "sender": "r@example.com",
                "url": req.url,
                "headers": {},
                "created": chrono::Utc::now().to_rfc3339(),
                "expires": chrono::Utc::now().to_rfc3339(),
                "status_code": 200,
                "body": base64::engine::general_purpose::STANDARD.encode(b"{}"),
            })
            .to_string(),
        )?;

        let resps = endpoint.check_responses()?;
        assert_eq!(resps.len(), 1);
        let (resp_file, _resp) = &resps[0];
        endpoint.cleanup_response(resp_file)?;
        assert!(!resp_file.exists());

        // 3) Invalid response JSON should be skipped
        let bad_resp_path = endpoint.path.join("oops.response");
        std::fs::write(&bad_resp_path, "not json")?;
        let resps2 = endpoint.check_responses()?;
        // It should skip invalid file and not error
        assert!(resps2.is_empty());

        // 5) Early-return branches when endpoint directory is missing
        // Remove the endpoint directory and ensure both checks return empty
        // Handle path existence check robustly (may not exist due to cleanup or symlink resolution)
        if endpoint.path.exists() {
            // Try to remove the directory, but don't fail if it's already gone
            let _ = std::fs::remove_dir_all(&endpoint.path);
        }
        let none_reqs = endpoint.check_requests()?;
        let none_resps = endpoint.check_responses()?;
        assert!(none_reqs.is_empty());
        assert!(none_resps.is_empty());

        // 6) Unrelated extension files are ignored
        std::fs::create_dir_all(&endpoint.path)?;
        std::fs::write(endpoint.path.join("note.txt"), b"ignore")?;
        let ig_reqs = endpoint.check_requests()?;
        let ig_resps = endpoint.check_responses()?;
        assert!(ig_reqs.is_empty());
        assert!(ig_resps.is_empty());

        // 4) send_response with invalid request path should return error
        let bogus = std::path::Path::new("");
        let ok_resp =
            crate::syftbox::types::RpcResponse::error(&req, "r@example.com".into(), 500, "err");
        // when request metadata is missing, send_response should fail
        assert!(endpoint
            .send_response(bogus, &req, &ok_resp, false)
            .is_err());

        // 7) send_response with directory-as-request should fail on remove_file
        let dir_as_req = endpoint.path.join("dirlike.request");
        std::fs::create_dir_all(&dir_as_req)?;
        // This will succeed writing response but fail deleting the dir
        let result = endpoint.send_response(&dir_as_req, &req, &ok_resp, false);
        assert!(result.is_err());

        // 8) create_request fails when a directory with same name exists
        let clash = RpcRequest::new(
            "sender@example.com".to_string(),
            app.build_syft_url("/message"),
            "POST".to_string(),
            b"Hi".to_vec(),
        );
        let clash_dir = endpoint.path.join(format!("{}.request", clash.id));
        std::fs::create_dir_all(&clash_dir)?;
        let create_err = endpoint.create_request(&clash);
        assert!(create_err.is_err());

        Ok(())
    }
}
