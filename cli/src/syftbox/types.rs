use anyhow::{anyhow, Result};
use base64::{engine::general_purpose::STANDARD as BASE64, Engine as _};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use uuid::Uuid;

/// RPC Request structure matching the SyftBox format
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RpcRequest {
    pub id: String,
    pub sender: String,
    pub url: String,
    #[serde(default)]
    pub headers: RpcHeaders,
    pub created: DateTime<Utc>,
    pub expires: DateTime<Utc>,
    pub method: String,
    pub body: String, // Base64 encoded
}

impl RpcRequest {
    /// Create a new RPC request
    pub fn new(sender: String, url: String, method: String, body: Vec<u8>) -> Self {
        let now = Utc::now();
        let id = Uuid::new_v4().to_string();

        Self {
            id,
            sender,
            url,
            headers: RpcHeaders::default(),
            created: now,
            expires: now + chrono::Duration::days(1),
            method,
            body: BASE64.encode(body),
        }
    }

    /// Decode the base64 body to bytes
    pub fn decode_body(&self) -> Result<Vec<u8>> {
        BASE64
            .decode(&self.body)
            .map_err(|e| anyhow!("Failed to decode body: {}", e))
    }

    /// Decode the body as a UTF-8 string
    pub fn body_as_string(&self) -> Result<String> {
        let bytes = self.decode_body()?;
        String::from_utf8(bytes).map_err(|e| anyhow!("Failed to decode body as UTF-8: {}", e))
    }

    /// Decode the body as JSON
    pub fn body_as_json<T: for<'de> Deserialize<'de>>(&self) -> Result<T> {
        let bytes = self.decode_body()?;
        serde_json::from_slice(&bytes).map_err(|e| anyhow!("Failed to parse body as JSON: {}", e))
    }
}

/// RPC Response structure matching the SyftBox format
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RpcResponse {
    pub id: String,
    pub sender: String,
    pub url: String,
    #[serde(default)]
    pub headers: RpcHeaders,
    pub created: DateTime<Utc>,
    pub expires: DateTime<Utc>,
    pub status_code: u16,
    pub body: String, // Base64 encoded
}

impl RpcResponse {
    /// Create a new RPC response for a request
    pub fn new(request: &RpcRequest, sender: String, status_code: u16, body: Vec<u8>) -> Self {
        let now = Utc::now();

        Self {
            id: request.id.clone(),
            sender,
            url: request.url.clone(),
            headers: RpcHeaders::default(),
            created: now,
            expires: now + chrono::Duration::days(1),
            status_code,
            body: BASE64.encode(body),
        }
    }

    /// Create a simple OK response with JSON body
    pub fn ok_json<T: Serialize>(request: &RpcRequest, sender: String, body: &T) -> Result<Self> {
        let json_bytes = serde_json::to_vec(body)?;
        let json_len = json_bytes.len();
        let mut response = Self::new(request, sender, 200, json_bytes);
        response
            .headers
            .insert("content-type".to_string(), "application/json".to_string());
        response
            .headers
            .insert("content-length".to_string(), json_len.to_string());
        Ok(response)
    }

    /// Create an error response
    pub fn error(request: &RpcRequest, sender: String, status_code: u16, message: &str) -> Self {
        let body = serde_json::json!({
            "error": message
        });
        let json_bytes = serde_json::to_vec(&body).unwrap_or_default();
        let json_len = json_bytes.len();
        let mut response = Self::new(request, sender, status_code, json_bytes);
        response
            .headers
            .insert("content-type".to_string(), "application/json".to_string());
        response
            .headers
            .insert("content-length".to_string(), json_len.to_string());
        response
    }

    /// Decode the base64 body to bytes
    pub fn decode_body(&self) -> Result<Vec<u8>> {
        BASE64
            .decode(&self.body)
            .map_err(|e| anyhow!("Failed to decode body: {}", e))
    }

    /// Decode the body as a UTF-8 string
    pub fn body_as_string(&self) -> Result<String> {
        let bytes = self.decode_body()?;
        String::from_utf8(bytes).map_err(|e| anyhow!("Failed to decode body as UTF-8: {}", e))
    }
}

/// RPC Headers - a simple HashMap wrapper
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct RpcHeaders(HashMap<String, String>);

impl RpcHeaders {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(&mut self, key: String, value: String) -> Option<String> {
        self.0.insert(key.to_lowercase(), value)
    }

    pub fn get(&self, key: &str) -> Option<&String> {
        self.0.get(&key.to_lowercase())
    }

    pub fn contains_key(&self, key: &str) -> bool {
        self.0.contains_key(&key.to_lowercase())
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &String)> {
        self.0.iter()
    }
}

impl From<HashMap<String, String>> for RpcHeaders {
    fn from(map: HashMap<String, String>) -> Self {
        Self(map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_request_body_encoding() {
        let request = RpcRequest::new(
            "test@example.com".to_string(),
            "syft://test@example.com/app_data/test/rpc/message".to_string(),
            "POST".to_string(),
            b"Hello, World!".to_vec(),
        );

        assert_eq!(request.body_as_string().unwrap(), "Hello, World!");
    }

    #[test]
    fn test_response_json() {
        let request = RpcRequest::new(
            "test@example.com".to_string(),
            "syft://test@example.com/app_data/test/rpc/message".to_string(),
            "POST".to_string(),
            b"{}".to_vec(),
        );

        let response = RpcResponse::ok_json(
            &request,
            "responder@example.com".to_string(),
            &serde_json::json!({"message": "Hello"}),
        )
        .unwrap();

        assert_eq!(response.status_code, 200);
        assert_eq!(
            response.headers.get("content-type").unwrap(),
            "application/json"
        );
    }
}
