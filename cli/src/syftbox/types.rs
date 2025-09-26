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

    #[test]
    fn request_body_as_json_and_error_paths() {
        // Success: parse JSON body
        let req_ok = RpcRequest::new(
            "a@b".into(),
            "syft://a@b/app_data/x/rpc/y".into(),
            "POST".into(),
            serde_json::to_vec(&serde_json::json!({"k":"v"})).unwrap(),
        );
        let v: serde_json::Value = req_ok.body_as_json().unwrap();
        assert_eq!(v["k"], "v");

        // Also parse into a concrete type to exercise another monomorphization
        #[derive(Deserialize, Debug)]
        struct K {
            k: String,
        }
        let req_typed = RpcRequest::new(
            "a@b".into(),
            "syft://a@b/app_data/x/rpc/y".into(),
            "POST".into(),
            serde_json::to_vec(&serde_json::json!({"k":"typed"})).unwrap(),
        );
        let kt: K = req_typed.body_as_json().unwrap();
        assert_eq!(kt.k, "typed");

        // Invalid base64 causes decode_body error
        let mut req_bad =
            RpcRequest::new("s".into(), "syft://s/p".into(), "GET".into(), b"x".to_vec());
        req_bad.body = "not_base64!!".into();
        assert!(req_bad.decode_body().is_err());

        // Valid base64 but invalid UTF-8 -> body_as_string error
        let bytes = vec![0xff, 0xff, 0xff];
        let req_bad_utf8 = RpcRequest::new("s".into(), "syft://s/p".into(), "POST".into(), bytes);
        assert!(req_bad_utf8.body_as_string().is_err());

        // Valid UTF-8 but invalid JSON -> body_as_json error
        let req_not_json = RpcRequest::new(
            "s".into(),
            "syft://s/p".into(),
            "POST".into(),
            b"not json".to_vec(),
        );
        let bad: Result<serde_json::Value> = req_not_json.body_as_json();
        assert!(bad.is_err());
    }

    #[test]
    fn response_ok_and_error_and_decoding() {
        let req = RpcRequest::new(
            "me@x".into(),
            "syft://me@x/app/rpc/ep".into(),
            "POST".into(),
            b"{}".to_vec(),
        );
        let resp_ok =
            RpcResponse::ok_json(&req, "me@x".into(), &serde_json::json!({"a":1})).unwrap();
        assert_eq!(resp_ok.status_code, 200);
        // Body decodes to JSON text string
        let s = resp_ok.body_as_string().unwrap();
        assert!(s.contains("\"a\":1"));

        let resp_err = RpcResponse::error(&req, "me@x".into(), 503, "unavailable");
        assert_eq!(resp_err.status_code, 503);
        let body = resp_err.body_as_string().unwrap();
        assert!(body.contains("error"));
        assert!(resp_err.headers.contains_key("content-length"));

        // Corrupt base64 results in decode error
        let mut resp_bad = resp_ok.clone();
        resp_bad.body = "!not base64!".into();
        assert!(resp_bad.decode_body().is_err());

        // Call new() directly to ensure coverage
        let direct = RpcResponse::new(&req, "me@x".into(), 201, b"{}".to_vec());
        assert_eq!(direct.status_code, 201);
    }

    #[test]
    fn headers_insert_get_contains_iter_and_from() {
        let mut h = RpcHeaders::new();
        // insert lowercases key for storage
        h.insert("Content-Type".into(), "application/json".into());
        assert_eq!(h.get("content-type"), Some(&"application/json".to_string()));
        assert!(h.contains_key("CONTENT-TYPE"));
        // Iter exposes items
        assert!(h.iter().count() >= 1);

        // Overwrite returns previous value via Option
        let prev = h.insert("content-type".into(), "text/plain".into());
        assert_eq!(prev, Some("application/json".to_string()));
        assert_eq!(h.get("CONTENT-TYPE").unwrap(), "text/plain");

        // From<HashMap> preserves original keys; get() with different case may not find it
        let mut map = std::collections::HashMap::new();
        map.insert("X-Thing".to_string(), "1".to_string());
        let h2 = RpcHeaders::from(map);
        // contains_key uses lowercase lookup, original key is not normalized
        assert!(!h2.contains_key("x-thing"));
        // But iter still shows the value, covering iteration path
        let items: Vec<_> = h2.iter().collect();
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].1, &"1".to_string());
        // Missing key -> None
        assert!(h2.get("missing").is_none());
    }
}
