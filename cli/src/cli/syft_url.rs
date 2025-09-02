use anyhow::{anyhow, Result};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct SyftURL {
    pub email: String,
    pub path: String,
    pub fragment: Option<String>,
}

impl SyftURL {
    #[allow(dead_code)]
    pub fn new(email: impl Into<String>, path: impl Into<String>) -> Self {
        Self {
            email: email.into(),
            path: path.into(),
            fragment: None,
        }
    }

    pub fn with_fragment(mut self, fragment: impl Into<String>) -> Self {
        self.fragment = Some(fragment.into());
        self
    }

    pub fn parse(url: &str) -> Result<Self> {
        if !url.starts_with("syft://") {
            return Err(anyhow!("Invalid SyftURL: must start with 'syft://'"));
        }

        let without_scheme = &url[7..]; // Remove "syft://"

        // Split by fragment if present
        let (main_part, fragment) = if let Some(hash_pos) = without_scheme.find('#') {
            let fragment = without_scheme[hash_pos + 1..].to_string();
            (&without_scheme[..hash_pos], Some(fragment))
        } else {
            (without_scheme, None)
        };

        // Find the first slash to separate email from path
        let slash_pos = main_part
            .find('/')
            .ok_or_else(|| anyhow!("Invalid SyftURL: missing path separator"))?;

        let email = main_part[..slash_pos].to_string();
        let path = main_part[slash_pos + 1..].to_string();

        Ok(Self {
            email,
            path,
            fragment,
        })
    }

    pub fn to_http_relay_url(&self, relay_server: &str) -> String {
        let relay_server = relay_server
            .trim_start_matches("https://")
            .trim_start_matches("http://");

        let mut url = format!(
            "https://{}/datasites/{}/{}",
            relay_server, self.email, self.path
        );
        if let Some(ref fragment) = self.fragment {
            url.push('#');
            url.push_str(fragment);
        }
        url
    }

    #[allow(dead_code)]
    pub fn from_http_relay_url(url: &str, relay_server: &str) -> Result<Self> {
        let relay_server = relay_server
            .trim_start_matches("https://")
            .trim_start_matches("http://");

        let prefix = format!("https://{}/datasites/", relay_server);
        let alt_prefix = format!("http://{}/datasites/", relay_server);

        let without_prefix = if url.starts_with(&prefix) {
            &url[prefix.len()..]
        } else if url.starts_with(&alt_prefix) {
            &url[alt_prefix.len()..]
        } else {
            return Err(anyhow!(
                "Invalid HTTP relay URL: doesn't match expected format"
            ));
        };

        // Split by fragment if present
        let (main_part, fragment) = if let Some(hash_pos) = without_prefix.find('#') {
            let fragment = without_prefix[hash_pos + 1..].to_string();
            (&without_prefix[..hash_pos], Some(fragment))
        } else {
            (without_prefix, None)
        };

        // Find the first slash to separate email from path
        let slash_pos = main_part
            .find('/')
            .ok_or_else(|| anyhow!("Invalid HTTP relay URL: missing path separator"))?;

        let email = main_part[..slash_pos].to_string();
        let path = main_part[slash_pos + 1..].to_string();

        Ok(Self {
            email,
            path,
            fragment,
        })
    }
}

impl fmt::Display for SyftURL {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut url = format!("syft://{}/{}", self.email, self.path);
        if let Some(ref fragment) = self.fragment {
            url.push('#');
            url.push_str(fragment);
        }
        write!(f, "{}", url)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_syft_url_construction() {
        let url = SyftURL::new("user@example.com", "public/data/file.yaml");
        assert_eq!(
            url.to_string(),
            "syft://user@example.com/public/data/file.yaml"
        );

        let url_with_fragment = url.with_fragment("participants.TEST1");
        assert_eq!(
            url_with_fragment.to_string(),
            "syft://user@example.com/public/data/file.yaml#participants.TEST1"
        );
    }

    #[test]
    fn test_syft_url_parsing() {
        let url = SyftURL::parse("syft://user@example.com/public/data/file.yaml").unwrap();
        assert_eq!(url.email, "user@example.com");
        assert_eq!(url.path, "public/data/file.yaml");
        assert_eq!(url.fragment, None);

        let url_with_fragment =
            SyftURL::parse("syft://user@example.com/public/data/file.yaml#participants.TEST1")
                .unwrap();
        assert_eq!(url_with_fragment.email, "user@example.com");
        assert_eq!(url_with_fragment.path, "public/data/file.yaml");
        assert_eq!(
            url_with_fragment.fragment,
            Some("participants.TEST1".to_string())
        );
    }

    #[test]
    fn test_http_relay_conversion() {
        let url = SyftURL::new("user@example.com", "public/data/file.yaml")
            .with_fragment("participants.TEST1");

        let http_url = url.to_http_relay_url("syftbox.net");
        assert_eq!(http_url, "https://syftbox.net/datasites/user@example.com/public/data/file.yaml#participants.TEST1");

        let parsed_back = SyftURL::from_http_relay_url(&http_url, "syftbox.net").unwrap();
        assert_eq!(parsed_back, url);
    }
}
