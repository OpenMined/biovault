use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::types::{RpcRequest, RpcResponse};
use anyhow::{Context, Result};
use std::fs;
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

    /// Check for request files in this endpoint
    pub fn check_requests(&self) -> Result<Vec<(PathBuf, RpcRequest)>> {
        let mut requests = Vec::new();

        if !self.path.exists() {
            return Ok(requests);
        }

        for entry in fs::read_dir(&self.path)? {
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
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read request file: {:?}", path))?;

        let request: RpcRequest = serde_json::from_str(&content)
            .with_context(|| format!("Failed to parse request JSON from: {:?}", path))?;

        Ok(request)
    }

    /// Send a response for a request
    pub fn send_response(&self, request_path: &Path, response: &RpcResponse) -> Result<()> {
        // Get the request ID from the filename (UUID.request -> UUID)
        let request_filename = request_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow::anyhow!("Invalid request filename"))?;

        // Create response file path (UUID.response)
        let response_filename = format!("{}.response", request_filename);
        let response_path = self.path.join(response_filename);

        // Write the response file
        let response_json =
            serde_json::to_string_pretty(&response).context("Failed to serialize response")?;

        fs::write(&response_path, response_json)
            .with_context(|| format!("Failed to write response file: {:?}", response_path))?;

        // Delete the request file
        fs::remove_file(request_path)
            .with_context(|| format!("Failed to delete request file: {:?}", request_path))?;

        println!("Sent response to: {:?}", response_path);
        Ok(())
    }

    /// Create a request file (for sending requests to others)
    pub fn create_request(&self, request: &RpcRequest) -> Result<PathBuf> {
        let request_filename = format!("{}.request", request.id);
        let request_path = self.path.join(request_filename);

        let request_json =
            serde_json::to_string_pretty(&request).context("Failed to serialize request")?;

        fs::write(&request_path, request_json)
            .with_context(|| format!("Failed to write request file: {:?}", request_path))?;

        println!("Created request: {:?}", request_path);
        Ok(request_path)
    }

    /// Check for response files (when we sent a request and are waiting for response)
    pub fn check_responses(&self) -> Result<Vec<(PathBuf, RpcResponse)>> {
        let mut responses = Vec::new();

        if !self.path.exists() {
            return Ok(responses);
        }

        for entry in fs::read_dir(&self.path)? {
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
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read response file: {:?}", path))?;

        let response: RpcResponse = serde_json::from_str(&content)
            .with_context(|| format!("Failed to parse response JSON from: {:?}", path))?;

        Ok(response)
    }

    /// Clean up a response file after processing
    pub fn cleanup_response(&self, response_path: &Path) -> Result<()> {
        fs::remove_file(response_path)
            .with_context(|| format!("Failed to delete response file: {:?}", response_path))?;
        Ok(())
    }
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

        endpoint.send_response(req_path, &response)?;

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
        std::fs::remove_dir_all(&endpoint.path)?;
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
        assert!(endpoint.send_response(bogus, &ok_resp).is_err());

        // 7) send_response with directory-as-request should fail on remove_file
        let dir_as_req = endpoint.path.join("dirlike.request");
        std::fs::create_dir_all(&dir_as_req)?;
        // This will succeed writing response but fail deleting the dir
        let result = endpoint.send_response(&dir_as_req, &ok_resp);
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
