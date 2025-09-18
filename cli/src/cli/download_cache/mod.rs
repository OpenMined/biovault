use anyhow::{Context, Result};
use blake3;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Read;
use std::path::Path;

pub mod downloader;
pub mod manifest;

pub use downloader::{DownloadCache, DownloadOptions};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChecksumPolicy {
    #[serde(rename = "type")]
    pub policy_type: ChecksumPolicyType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_hash: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChecksumPolicyType {
    Required,  // Must have expected_hash and it must match
    Preferred, // If expected_hash provided, warn if different but continue
    Optional,  // No hash validation required
}

#[derive(Debug, Clone)]
pub struct CacheStrategy {
    pub check_remote: bool,
}

impl Default for CacheStrategy {
    fn default() -> Self {
        Self { check_remote: true }
    }
}

pub fn calculate_blake3(path: &Path) -> Result<String> {
    let mut file = fs::File::open(path)
        .with_context(|| format!("Failed to open file for checksum: {}", path.display()))?;

    let metadata = file.metadata()?;
    let file_size = metadata.len();

    if file_size > 100 * 1024 * 1024 {
        // For files > 100MB, use parallel hashing
        let mut hasher = blake3::Hasher::new();
        let mut buffer = vec![0; 64 * 1024 * 1024]; // 64MB buffer

        loop {
            let bytes_read = file.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            hasher.update_rayon(&buffer[..bytes_read]);
        }

        Ok(hasher.finalize().to_hex().to_string())
    } else {
        // For smaller files, use regular update
        let mut hasher = blake3::Hasher::new();
        let mut buffer = vec![0; 8 * 1024 * 1024]; // 8MB buffer

        loop {
            let bytes_read = file.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            hasher.update(&buffer[..bytes_read]);
        }

        Ok(hasher.finalize().to_hex().to_string())
    }
}

/// Cross-platform file linking/copying
/// On Unix: creates hard link or symlink
/// On Windows: copies the file
pub fn link_or_copy(source: &Path, target: &Path) -> Result<()> {
    // Ensure target directory exists
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent).with_context(|| {
            format!("Failed to create parent directory for {}", target.display())
        })?;
    }

    // Remove target if it exists
    if target.exists() {
        fs::remove_file(target)
            .with_context(|| format!("Failed to remove existing file at {}", target.display()))?;
    }

    #[cfg(unix)]
    {
        // Try hard link first (more efficient, works across bind mounts)
        if fs::hard_link(source, target).is_ok() {
            return Ok(());
        }

        // Fall back to symlink
        use std::os::unix::fs::symlink;
        symlink(source, target).with_context(|| {
            format!(
                "Failed to create symlink from {} to {}",
                source.display(),
                target.display()
            )
        })?;
    }

    #[cfg(windows)]
    {
        // On Windows, just copy the file
        fs::copy(source, target).with_context(|| {
            format!(
                "Failed to copy file from {} to {}",
                source.display(),
                target.display()
            )
        })?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_calculate_blake3() {
        let temp_dir = TempDir::new().unwrap();
        let test_file = temp_dir.path().join("test.txt");
        fs::write(&test_file, b"Hello, World!").unwrap();

        let hash = calculate_blake3(&test_file).unwrap();
        assert!(!hash.is_empty());
        assert_eq!(hash.len(), 64); // BLAKE3 produces 32 bytes = 64 hex chars
    }

    #[test]
    fn test_link_or_copy() {
        let temp_dir = TempDir::new().unwrap();
        let source = temp_dir.path().join("source.txt");
        let target = temp_dir.path().join("subdir/target.txt");

        fs::write(&source, b"test content").unwrap();

        link_or_copy(&source, &target).unwrap();

        assert!(target.exists());
        assert_eq!(fs::read(&target).unwrap(), b"test content");
    }

    #[test]
    fn test_cache_strategy_default_true() {
        let s = CacheStrategy::default();
        assert!(s.check_remote);
    }

    #[test]
    fn test_calculate_blake3_open_error() {
        let missing = Path::new("/definitely/not/here.bin");
        let err = calculate_blake3(missing).err().unwrap();
        let msg = format!("{}", err);
        assert!(msg.contains("Failed to open file for checksum"));
    }

    #[test]
    fn test_calculate_blake3_large_parallel() {
        // Create a file > 100MB to hit the parallel hashing path
        let temp_dir = TempDir::new().unwrap();
        let big = temp_dir.path().join("big.bin");
        let mut f = std::fs::File::create(&big).unwrap();
        // Write 101 MiB in 8 MiB chunks
        let chunk = vec![0u8; 8 * 1024 * 1024];
        let mut written = 0usize;
        while written < 101 * 1024 * 1024 {
            let to_write = std::cmp::min(chunk.len(), 101 * 1024 * 1024 - written);
            use std::io::Write;
            f.write_all(&chunk[..to_write]).unwrap();
            written += to_write;
        }
        drop(f);
        let h = calculate_blake3(&big).unwrap();
        assert_eq!(h.len(), 64);
    }

    #[test]
    fn test_link_or_copy_parent_dir_create_error() {
        let temp_dir = TempDir::new().unwrap();
        let source = temp_dir.path().join("src.txt");
        std::fs::write(&source, b"x").unwrap();
        // Create a file where the parent dir should be
        let parent_file = temp_dir.path().join("subdir");
        std::fs::write(&parent_file, b"not a dir").unwrap();
        let target = parent_file.join("target.txt");
        let err = link_or_copy(&source, &target).err().unwrap();
        let msg = format!("{}", err);
        assert!(msg.contains("Failed to create parent directory"));
    }

    #[test]
    fn test_link_or_copy_remove_existing_file_error() {
        let temp_dir = TempDir::new().unwrap();
        let source = temp_dir.path().join("src2.txt");
        std::fs::write(&source, b"y").unwrap();
        let target = temp_dir.path().join("exist/target.txt");
        // Create parent dir and then create a directory at the target path
        std::fs::create_dir_all(target.parent().unwrap()).unwrap();
        std::fs::create_dir(&target).unwrap();
        let err = link_or_copy(&source, &target).err().unwrap();
        let msg = format!("{}", err);
        assert!(msg.contains("Failed to remove existing file"));
    }
}
