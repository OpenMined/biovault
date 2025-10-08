use anyhow::Result;
use serde::Serialize;

/// Wrapper for all CLI JSON responses with API versioning
#[derive(Serialize)]
pub struct CliResponse<T> {
    pub api_version: String,
    pub data: T,
}

impl<T: Serialize> CliResponse<T> {
    pub fn new(data: T) -> Self {
        Self {
            api_version: env!("CARGO_PKG_VERSION").to_string(),
            data,
        }
    }

    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    pub fn to_json_compact(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct TestData {
        message: String,
        count: u32,
    }

    #[test]
    fn test_response_serialization() {
        let data = TestData {
            message: "test".to_string(),
            count: 42,
        };
        let response = CliResponse::new(data);

        let json = response.to_json().unwrap();
        assert!(json.contains("api_version"));
        assert!(json.contains("test"));
        assert!(json.contains("42"));
    }

    #[test]
    fn test_response_has_version() {
        let data = TestData {
            message: "test".to_string(),
            count: 1,
        };
        let response = CliResponse::new(data);

        assert!(!response.api_version.is_empty());
    }
}
