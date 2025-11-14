use crate::syftbox::app::SyftBoxApp;
use crate::syftbox::endpoint::Endpoint;
use crate::syftbox::types::{RpcRequest, RpcResponse};
use anyhow::Result;
use std::path::{Path, PathBuf};

/// Check for requests across all endpoints of an app
pub fn check_all_requests(app: &SyftBoxApp) -> Result<Vec<(String, PathBuf, RpcRequest)>> {
    let mut all_requests = Vec::new();

    let endpoints = app.list_endpoints()?;
    for endpoint_name in endpoints {
        let endpoint = Endpoint::new(app, &endpoint_name)?;
        let requests = endpoint.check_requests()?;

        for (path, request) in requests {
            all_requests.push((endpoint_name.clone(), path, request));
        }
    }

    Ok(all_requests)
}

/// Check for requests on a specific endpoint
pub fn check_requests(app: &SyftBoxApp, endpoint_name: &str) -> Result<Vec<(PathBuf, RpcRequest)>> {
    let endpoint = Endpoint::new(app, endpoint_name)?;
    endpoint.check_requests()
}

/// Send a response to a request
pub fn send_response(
    app: &SyftBoxApp,
    endpoint_name: &str,
    request_path: &Path,
    request: &RpcRequest,
    response: &RpcResponse,
) -> Result<()> {
    let endpoint = Endpoint::new(app, endpoint_name)?;
    endpoint.send_response(request_path, request, response)
}

/// Create and send a request to another datasite
pub fn send_request(
    app: &SyftBoxApp,
    endpoint_name: &str,
    request: &RpcRequest,
) -> Result<PathBuf> {
    let endpoint = Endpoint::new(app, endpoint_name)?;
    endpoint.create_request(request)
}

/// Check for responses to requests we've sent
pub fn check_responses(
    app: &SyftBoxApp,
    endpoint_name: &str,
) -> Result<Vec<(PathBuf, RpcResponse)>> {
    let endpoint = Endpoint::new(app, endpoint_name)?;
    endpoint.check_responses()
}

/// Helper to process a request and automatically send a response
pub fn process_request<F>(
    app: &SyftBoxApp,
    endpoint_name: &str,
    request_path: &Path,
    request: &RpcRequest,
    handler: F,
) -> Result<()>
where
    F: FnOnce(&RpcRequest) -> Result<RpcResponse>,
{
    match handler(request) {
        Ok(response) => {
            send_response(app, endpoint_name, request_path, request, &response)?;
            Ok(())
        }
        Err(e) => {
            // Send error response
            let error_response = RpcResponse::error(
                request,
                app.email.clone(),
                500,
                &format!("Error processing request: {}", e),
            );
            send_response(app, endpoint_name, request_path, request, &error_response)?;
            Err(e)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_check_all_requests() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;

        // Register multiple endpoints
        app.register_endpoint("/message")?;
        app.register_endpoint("/data")?;

        // Create requests in different endpoints
        let req1 = RpcRequest::new(
            "sender1@example.com".to_string(),
            app.build_syft_url("/message"),
            "POST".to_string(),
            b"Message 1".to_vec(),
        );

        let req2 = RpcRequest::new(
            "sender2@example.com".to_string(),
            app.build_syft_url("/data"),
            "GET".to_string(),
            b"".to_vec(),
        );

        send_request(&app, "/message", &req1)?;
        send_request(&app, "/data", &req2)?;

        // Check all requests
        let all_requests = check_all_requests(&app)?;
        assert_eq!(all_requests.len(), 2);

        Ok(())
    }

    #[test]
    fn test_process_request() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;

        let endpoint_name = "/echo";
        app.register_endpoint(endpoint_name)?;

        // Create a request
        let request = RpcRequest::new(
            "sender@example.com".to_string(),
            app.build_syft_url(endpoint_name),
            "POST".to_string(),
            serde_json::to_vec(&serde_json::json!({"message": "Hello"}))?,
        );

        let request_path = send_request(&app, endpoint_name, &request)?;

        // Get the request
        let requests = check_requests(&app, endpoint_name)?;
        assert_eq!(requests.len(), 1);

        let (req_path, req) = &requests[0];

        // Process the request with a handler
        process_request(&app, endpoint_name, req_path, req, |request| {
            let body: serde_json::Value = request.body_as_json()?;
            RpcResponse::ok_json(
                request,
                app.email.clone(),
                &serde_json::json!({
                    "echo": body["message"],
                    "processed": true
                }),
            )
        })?;

        // Verify request was processed and response created
        assert!(!request_path.exists());

        let responses = check_responses(&app, endpoint_name)?;
        assert_eq!(responses.len(), 1);

        Ok(())
    }
}
