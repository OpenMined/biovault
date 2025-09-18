use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub version: String,
    pub downloads: Downloads,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Downloads {
    pub by_url: HashMap<String, UrlEntry>,
    pub by_hash: HashMap<String, HashEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UrlEntry {
    pub current_hash: String,
    pub last_checked: DateTime<Utc>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_modified: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub etag: Option<String>,
    pub history: Vec<HistoryEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HistoryEntry {
    pub hash: String,
    pub first_seen: DateTime<Utc>,
    pub size: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HashEntry {
    pub urls: Vec<String>,
    pub size: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mime_type: Option<String>,
    pub created: DateTime<Utc>,
    pub last_accessed: DateTime<Utc>,
    pub reference_count: usize,
}

impl Default for Manifest {
    fn default() -> Self {
        Self::new()
    }
}

impl Manifest {
    pub fn new() -> Self {
        Self {
            version: "1.0.0".to_string(),
            downloads: Downloads {
                by_url: HashMap::new(),
                by_hash: HashMap::new(),
            },
        }
    }

    pub fn load(path: &Path) -> Result<Self> {
        if !path.exists() {
            return Ok(Self::new());
        }

        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read manifest from {}", path.display()))?;

        serde_yaml::from_str(&content)
            .with_context(|| format!("Failed to parse manifest from {}", path.display()))
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        // Ensure parent directory exists
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Failed to create directory for manifest at {}",
                    path.display()
                )
            })?;
        }

        let yaml = serde_yaml::to_string(self).context("Failed to serialize manifest")?;

        fs::write(path, yaml)
            .with_context(|| format!("Failed to write manifest to {}", path.display()))?;

        Ok(())
    }

    pub fn get_by_url(&self, url: &str) -> Option<&UrlEntry> {
        self.downloads.by_url.get(url)
    }

    pub fn get_by_hash(&self, hash: &str) -> Option<&HashEntry> {
        self.downloads.by_hash.get(hash)
    }

    pub fn add_download(
        &mut self,
        url: String,
        hash: String,
        size: u64,
        etag: Option<String>,
        last_modified: Option<String>,
    ) {
        let now = Utc::now();

        // Check if URL entry exists and hash changed
        let needs_history_update = self
            .downloads
            .by_url
            .get(&url)
            .map(|e| e.current_hash != hash)
            .unwrap_or(false);

        if needs_history_update {
            // Add current hash to history before updating
            if let Some(url_entry) = self.downloads.by_url.get(&url) {
                let history_entry = HistoryEntry {
                    hash: url_entry.current_hash.clone(),
                    first_seen: url_entry.last_checked,
                    size: self
                        .downloads
                        .by_hash
                        .get(&url_entry.current_hash)
                        .map(|e| e.size)
                        .unwrap_or(0),
                };

                // Update with history
                let mut new_entry = url_entry.clone();
                new_entry.history.push(history_entry);
                new_entry.current_hash = hash.clone();
                new_entry.last_checked = now;
                new_entry.etag = etag.clone();
                new_entry.last_modified = last_modified.clone();

                self.downloads.by_url.insert(url.clone(), new_entry);
            }
        } else {
            // Update or create URL entry
            self.downloads
                .by_url
                .entry(url.clone())
                .and_modify(|e| {
                    e.last_checked = now;
                    e.etag = etag.clone();
                    e.last_modified = last_modified.clone();
                })
                .or_insert_with(|| UrlEntry {
                    current_hash: hash.clone(),
                    last_checked: now,
                    last_modified,
                    etag,
                    history: vec![],
                });
        }

        // Update or create hash entry
        let hash_entry = self
            .downloads
            .by_hash
            .entry(hash)
            .or_insert_with(|| HashEntry {
                urls: vec![],
                size,
                mime_type: None,
                created: now,
                last_accessed: now,
                reference_count: 0,
            });

        if !hash_entry.urls.contains(&url) {
            hash_entry.urls.push(url);
        }
        hash_entry.last_accessed = now;
    }

    pub fn update_last_accessed(&mut self, hash: &str) {
        if let Some(entry) = self.downloads.by_hash.get_mut(hash) {
            entry.last_accessed = Utc::now();
        }
    }

    pub fn increment_reference_count(&mut self, hash: &str) {
        if let Some(entry) = self.downloads.by_hash.get_mut(hash) {
            entry.reference_count += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_manifest_round_trip() {
        let temp_dir = TempDir::new().unwrap();
        let manifest_path = temp_dir.path().join("manifest.yaml");

        let mut manifest = Manifest::new();
        manifest.add_download(
            "https://example.com/file.txt".to_string(),
            "abc123".to_string(),
            1024,
            Some("\"etag123\"".to_string()),
            Some("2024-01-01T00:00:00Z".to_string()),
        );

        manifest.save(&manifest_path).unwrap();
        let loaded = Manifest::load(&manifest_path).unwrap();

        assert_eq!(loaded.downloads.by_url.len(), 1);
        assert_eq!(loaded.downloads.by_hash.len(), 1);
        assert!(loaded.get_by_url("https://example.com/file.txt").is_some());
        assert!(loaded.get_by_hash("abc123").is_some());
    }

    #[test]
    fn manifest_load_default_and_save_parent_dirs() {
        let temp_dir = TempDir::new().unwrap();
        let nonexist = temp_dir.path().join("no/such/manifest.yaml");
        let m = Manifest::load(&nonexist).unwrap();
        assert_eq!(m.downloads.by_url.len(), 0);
        m.save(&nonexist).unwrap();
        assert!(nonexist.exists());
    }

    #[test]
    fn manifest_add_download_history_and_updates() {
        let mut m = Manifest::new();
        m.add_download(
            "https://example.com/f".into(),
            "hash1".into(),
            10,
            Some("etag1".into()),
            Some("2024-01-01".into()),
        );
        m.add_download(
            "https://example.com/f".into(),
            "hash2".into(),
            10,
            Some("etag2".into()),
            Some("2024-01-02".into()),
        );

        let u = m.get_by_url("https://example.com/f").unwrap();
        assert_eq!(u.current_hash, "hash2");
        assert!(!u.history.is_empty());

        let h2 = m.get_by_hash("hash2").unwrap();
        assert!(h2.urls.contains(&"https://example.com/f".into()));

        m.update_last_accessed("hash2");
        m.increment_reference_count("hash2");
        assert!(m.get_by_hash("hash2").unwrap().reference_count >= 1);
    }
}
