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
        )?;

        endpoint.send_response(req_path, &response)?;

        // Verify request was deleted and response created
        assert!(!request_path.exists());

        let responses = endpoint.check_responses()?;
        assert_eq!(responses.len(), 1);

        Ok(())
    }
}
