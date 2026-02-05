use super::manifest::Manifest;
use super::{calculate_blake3, link_or_copy, CacheStrategy, ChecksumPolicy, ChecksumPolicyType};
use anyhow::{anyhow, Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::header::{CONTENT_LENGTH, ETAG, LAST_MODIFIED, RANGE};
use reqwest::StatusCode;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tokio::fs::{File, OpenOptions};
use tokio::io::{AsyncSeekExt, AsyncWriteExt};

pub struct DownloadCache {
    cache_dir: PathBuf,
    manifest_path: PathBuf,
    manifest: Manifest,
}

#[derive(Clone)]
pub struct DownloadOptions {
    pub checksum_policy: ChecksumPolicy,
    pub cache_strategy: CacheStrategy,
    pub show_progress: bool,
    pub progress_callback: Option<Arc<dyn Fn(u64, u64) + Send + Sync>>,
}

impl std::fmt::Debug for DownloadOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DownloadOptions")
            .field("checksum_policy", &self.checksum_policy)
            .field("cache_strategy", &self.cache_strategy)
            .field("show_progress", &self.show_progress)
            .field("progress_callback", &self.progress_callback.is_some())
            .finish()
    }
}

impl Default for DownloadOptions {
    fn default() -> Self {
        Self {
            checksum_policy: ChecksumPolicy {
                policy_type: ChecksumPolicyType::Optional,
                expected_hash: None,
            },
            cache_strategy: CacheStrategy::default(),
            show_progress: true,
            progress_callback: None,
        }
    }
}

impl DownloadCache {
    pub fn new(cache_dir: Option<PathBuf>) -> Result<Self> {
        let cache_dir = if let Some(dir) = cache_dir {
            dir
        } else {
            // Use the shared cache directory
            crate::config::get_cache_dir()?
        };

        fs::create_dir_all(&cache_dir).with_context(|| {
            format!("Failed to create cache directory: {}", cache_dir.display())
        })?;

        let manifest_path = cache_dir
            .parent()
            .ok_or_else(|| anyhow!("Invalid cache directory"))?
            .join("manifest.yaml");

        let manifest = Manifest::load(&manifest_path)?;

        Ok(Self {
            cache_dir,
            manifest_path,
            manifest,
        })
    }

    pub async fn download_with_cache(
        &mut self,
        url: &str,
        target_path: &Path,
        options: DownloadOptions,
    ) -> Result<PathBuf> {
        let verbose = options.show_progress;
        if verbose {
            println!("📦 Processing: {}", url);
        }

        // Check if we have the expected hash in cache
        if let Some(ref expected_hash) = options.checksum_policy.expected_hash {
            if let Some(_hash_entry) = self.manifest.get_by_hash(expected_hash) {
                let cache_path = self.cache_dir.join("by-hash").join(expected_hash);
                if cache_path.exists() {
                    if verbose {
                        println!("  ✓ Found in cache (hash: {}...)", &expected_hash[..8]);
                    }
                    if verbose {
                        println!("  ✓ Cache hit! Using cached version");
                    }

                    // Update last accessed time
                    self.manifest.update_last_accessed(expected_hash);
                    self.manifest.increment_reference_count(expected_hash);
                    self.save_manifest()?;

                    // Link or copy to target
                    link_or_copy(&cache_path, target_path)?;
                    return Ok(target_path.to_path_buf());
                }
            }
        }

        // Check if URL was downloaded before
        let cached_info = self.manifest.get_by_url(url).map(|e| {
            (
                e.current_hash.clone(),
                e.etag.clone(),
                e.last_modified.clone(),
            )
        });

        let should_revalidate =
            if let Some((current_hash, cached_etag, cached_last_modified)) = cached_info {
                let cache_path = self.cache_dir.join("by-hash").join(&current_hash);

                if cache_path.exists() {
                    if options.cache_strategy.check_remote {
                        if verbose {
                            println!("  ↻ Checking if remote file has changed...");
                        }

                        // Perform HEAD request to check if file changed
                        let client = reqwest::Client::new();
                        let head_response = client.head(url).send().await?;

                        let remote_etag = head_response
                            .headers()
                            .get(ETAG)
                            .and_then(|v| v.to_str().ok())
                            .map(|s| s.to_string());

                        let remote_last_modified = head_response
                            .headers()
                            .get(LAST_MODIFIED)
                            .and_then(|v| v.to_str().ok())
                            .map(|s| s.to_string());

                        // Check if file has changed
                        let has_changed = match (&cached_etag, &remote_etag) {
                            (Some(cached), Some(remote)) if cached == remote => false,
                            _ => match (&cached_last_modified, &remote_last_modified) {
                                (Some(cached), Some(remote)) if cached == remote => false,
                                _ => true, // When in doubt, re-download
                            },
                        };

                        if !has_changed {
                            if verbose {
                                println!("  ✓ Remote file unchanged");
                            }
                            if verbose {
                                println!(
                                    "  ✓ Cache hit! Using cached version (hash: {}...)",
                                    &current_hash[..8]
                                );
                            }

                            // Update manifest
                            self.manifest.update_last_accessed(&current_hash);
                            self.manifest.increment_reference_count(&current_hash);
                            self.save_manifest()?;

                            // Link or copy to target
                            link_or_copy(&cache_path, target_path)?;
                            return Ok(target_path.to_path_buf());
                        } else {
                            println!("  ⚠ Remote file has changed, will re-download");
                            true
                        }
                    } else {
                        // Not checking remote, use cached version
                        if verbose {
                            println!(
                                "  ✓ Cache hit! Using cached version (hash: {}...)",
                                &current_hash[..8]
                            );
                        }

                        self.manifest.update_last_accessed(&current_hash);
                        self.manifest.increment_reference_count(&current_hash);
                        self.save_manifest()?;

                        link_or_copy(&cache_path, target_path)?;
                        return Ok(target_path.to_path_buf());
                    }
                } else {
                    true // Cache file missing, need to download
                }
            } else {
                true // Never downloaded before
            };

        if should_revalidate {
            if verbose {
                println!("  ↓ Cache miss - downloading file...");
            }

            // Download to temporary location using URL hash for resumability
            // This allows resuming interrupted downloads since the temp filename is predictable
            let url_hash = blake3::hash(url.as_bytes()).to_hex().to_string();
            let temp_path = self
                .cache_dir
                .join("downloads")
                .join(format!("{}.partial", &url_hash[..16]));
            fs::create_dir_all(temp_path.parent().unwrap())?;

            let (etag, last_modified) = self
                .download_file(
                    url,
                    &temp_path,
                    options.show_progress,
                    options.progress_callback.clone(),
                )
                .await?;

            // Calculate hash
            if verbose {
                print!("  ⟳ Computing BLAKE3 checksum... ");
                std::io::Write::flush(&mut std::io::stdout())?;
            }
            let actual_hash = calculate_blake3(&temp_path)?;
            if verbose {
                println!("done ({}...)", &actual_hash[..8]);
            }

            // Validate checksum if required
            match options.checksum_policy.policy_type {
                ChecksumPolicyType::Required => {
                    let expected = options
                        .checksum_policy
                        .expected_hash
                        .ok_or_else(|| anyhow!("Required checksum not provided"))?;
                    if actual_hash != expected {
                        fs::remove_file(&temp_path)?;
                        return Err(anyhow!(
                            "Checksum mismatch! Expected: {}..., Got: {}...",
                            &expected[..8],
                            &actual_hash[..8]
                        ));
                    }
                    if verbose {
                        println!("  ✓ Checksum verified");
                    }
                }
                ChecksumPolicyType::Preferred => {
                    if let Some(expected) = &options.checksum_policy.expected_hash {
                        if &actual_hash != expected {
                            if verbose {
                                println!("  ⚠ Warning: Checksum mismatch!");
                            }
                            if verbose {
                                println!("    Expected: {}...", &expected[..8]);
                            }
                            if verbose {
                                println!("    Got:      {}...", &actual_hash[..8]);
                            }
                            if verbose {
                                println!("  ⚠ Continuing anyway (policy: preferred)");
                            }
                        } else if verbose {
                            println!("  ✓ Checksum verified");
                        }
                    }
                }
                ChecksumPolicyType::Optional => {
                    // No validation needed
                }
            }

            // Move to cache
            let cache_path = self.cache_dir.join("by-hash").join(&actual_hash);
            fs::create_dir_all(cache_path.parent().unwrap())?;

            if cache_path.exists() {
                // File already in cache (possibly from another URL)
                if verbose {
                    println!("  ✓ File already in cache with same content");
                }
                fs::remove_file(&temp_path)?;
            } else {
                fs::rename(&temp_path, &cache_path)?;
                if verbose {
                    println!("  ✓ Added to cache");
                }
            }

            // Update manifest
            let file_size = fs::metadata(&cache_path)?.len();
            self.manifest.add_download(
                url.to_string(),
                actual_hash.clone(),
                file_size,
                etag,
                last_modified,
            );
            self.manifest.increment_reference_count(&actual_hash);
            self.save_manifest()?;

            // Link or copy to target
            link_or_copy(&cache_path, target_path)?;
            if verbose {
                println!("  ✓ File ready at: {}", target_path.display());
            }

            Ok(target_path.to_path_buf())
        } else {
            Ok(target_path.to_path_buf())
        }
    }

    async fn download_file(
        &self,
        url: &str,
        target_path: &Path,
        show_progress: bool,
        progress_callback: Option<Arc<dyn Fn(u64, u64) + Send + Sync>>,
    ) -> Result<(Option<String>, Option<String>)> {
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(3600))
            .build()?;

        // Check if partial download exists
        let existing_size = if target_path.exists() {
            fs::metadata(target_path).map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        // First, do a HEAD request to get total size and check if server supports Range
        let head_response = client.head(url).send().await?;
        if !head_response.status().is_success() {
            return Err(anyhow!(
                "HTTP HEAD request failed: {}",
                head_response.status()
            ));
        }

        let total_size = head_response
            .headers()
            .get(CONTENT_LENGTH)
            .and_then(|v| v.to_str().ok())
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);

        let accepts_ranges = head_response
            .headers()
            .get("accept-ranges")
            .and_then(|v| v.to_str().ok())
            .map(|s| s != "none")
            .unwrap_or(false);

        // Extract headers for caching
        let etag = head_response
            .headers()
            .get(ETAG)
            .and_then(|v| v.to_str().ok())
            .map(|s| s.to_string());

        let last_modified = head_response
            .headers()
            .get(LAST_MODIFIED)
            .and_then(|v| v.to_str().ok())
            .map(|s| s.to_string());

        // Determine if we can resume
        let resume_from = if existing_size > 0 && accepts_ranges && total_size > 0 {
            if existing_size >= total_size {
                // File already complete
                if show_progress {
                    println!(
                        "    ✓ File already fully downloaded ({} bytes)",
                        existing_size
                    );
                }
                return Ok((etag, last_modified));
            }
            if show_progress {
                println!(
                    "    ↻ Resuming download from {} / {} bytes ({:.1}%)",
                    existing_size,
                    total_size,
                    (existing_size as f64 / total_size as f64) * 100.0
                );
            }
            Some(existing_size)
        } else if existing_size > 0 && !accepts_ranges {
            if show_progress {
                println!("    ⚠ Server doesn't support resume, starting fresh");
            }
            // Remove partial file
            let _ = fs::remove_file(target_path);
            None
        } else {
            None
        };

        // Build the GET request with optional Range header
        let mut request = client.get(url);
        if let Some(offset) = resume_from {
            request = request.header(RANGE, format!("bytes={}-", offset));
        }

        let response = request.send().await?;
        let status = response.status();

        // Handle response status
        if status == StatusCode::RANGE_NOT_SATISFIABLE {
            // Range not satisfiable - file might be complete or changed
            if show_progress {
                println!("    ⚠ Range not satisfiable, starting fresh");
            }
            let _ = fs::remove_file(target_path);
            // Retry without range
            return self
                .download_file_fresh(
                    url,
                    target_path,
                    show_progress,
                    progress_callback,
                    total_size,
                )
                .await;
        }

        if !status.is_success() {
            return Err(anyhow!("HTTP request failed: {}", status));
        }

        // Determine actual start position and total
        let (start_pos, content_length) = if status == StatusCode::PARTIAL_CONTENT {
            // Server returned partial content
            let content_len = response
                .content_length()
                .unwrap_or(total_size - resume_from.unwrap_or(0));
            (resume_from.unwrap_or(0), content_len)
        } else {
            // Server returned full content (200 OK) - start from beginning
            if resume_from.is_some() {
                // We asked for range but got full content - remove partial and start fresh
                let _ = fs::remove_file(target_path);
            }
            (0, response.content_length().unwrap_or(total_size))
        };

        let display_total = if start_pos > 0 {
            start_pos + content_length
        } else {
            content_length
        };

        let pb = if show_progress && display_total > 0 {
            let pb = ProgressBar::new(display_total);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("    [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                    .expect("Failed to set progress bar template")
                    .progress_chars("#>-"),
            );
            if start_pos > 0 {
                pb.set_position(start_pos);
            }
            Some(pb)
        } else if show_progress {
            println!("    Downloading (size unknown)...");
            None
        } else {
            None
        };

        // Open file for writing - append if resuming, create if new
        let mut file = if start_pos > 0 {
            let mut f = OpenOptions::new().write(true).open(target_path).await?;
            f.seek(std::io::SeekFrom::End(0)).await?;
            f
        } else {
            File::create(target_path).await?
        };

        let mut downloaded = start_pos;
        let mut stream = response.bytes_stream();

        while let Some(chunk) = futures_util::StreamExt::next(&mut stream).await {
            let chunk = chunk?;
            file.write_all(&chunk).await?;
            downloaded += chunk.len() as u64;

            if let Some(ref pb) = pb {
                pb.set_position(downloaded);
            }
            if let Some(ref cb) = progress_callback {
                cb(downloaded, display_total);
            }
        }

        // Ensure all data is flushed
        file.flush().await?;

        if let Some(pb) = pb {
            pb.finish_and_clear();
        }

        Ok((etag, last_modified))
    }

    /// Download file from scratch without attempting to resume
    async fn download_file_fresh(
        &self,
        url: &str,
        target_path: &Path,
        show_progress: bool,
        progress_callback: Option<Arc<dyn Fn(u64, u64) + Send + Sync>>,
        total_size: u64,
    ) -> Result<(Option<String>, Option<String>)> {
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(3600))
            .build()?;

        let response = client.get(url).send().await?;

        if !response.status().is_success() {
            return Err(anyhow!("HTTP request failed: {}", response.status()));
        }

        let etag = response
            .headers()
            .get(ETAG)
            .and_then(|v| v.to_str().ok())
            .map(|s| s.to_string());

        let last_modified = response
            .headers()
            .get(LAST_MODIFIED)
            .and_then(|v| v.to_str().ok())
            .map(|s| s.to_string());

        let content_length = response.content_length().unwrap_or(total_size);

        let pb = if show_progress && content_length > 0 {
            let pb = ProgressBar::new(content_length);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("    [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                    .expect("Failed to set progress bar template")
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else {
            None
        };

        let mut file = File::create(target_path).await?;
        let mut downloaded = 0u64;
        let mut stream = response.bytes_stream();

        while let Some(chunk) = futures_util::StreamExt::next(&mut stream).await {
            let chunk = chunk?;
            file.write_all(&chunk).await?;
            downloaded += chunk.len() as u64;

            if let Some(ref pb) = pb {
                pb.set_position(downloaded);
            }
            if let Some(ref cb) = progress_callback {
                cb(downloaded, content_length);
            }
        }

        file.flush().await?;

        if let Some(pb) = pb {
            pb.finish_and_clear();
        }

        Ok((etag, last_modified))
    }

    fn save_manifest(&self) -> Result<()> {
        self.manifest.save(&self.manifest_path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::download_cache::manifest::Manifest;
    use tempfile::TempDir;

    fn write_manifest(parent_dir: &Path, m: &Manifest) {
        let path = parent_dir.join("manifest.yaml");
        m.save(&path).unwrap();
    }

    #[test]
    fn test_download_options_default() {
        let opts = DownloadOptions::default();
        assert_eq!(
            opts.checksum_policy.policy_type,
            ChecksumPolicyType::Optional
        );
        assert!(opts.checksum_policy.expected_hash.is_none());
        assert!(opts.cache_strategy.check_remote);
        assert!(opts.show_progress);
    }

    #[test]
    fn test_download_options_clone() {
        let opts1 = DownloadOptions {
            checksum_policy: ChecksumPolicy {
                policy_type: ChecksumPolicyType::Required,
                expected_hash: Some("hash123".to_string()),
            },
            cache_strategy: CacheStrategy {
                check_remote: false,
            },
            show_progress: false,
            progress_callback: None,
        };
        let opts2 = opts1.clone();
        assert_eq!(
            opts2.checksum_policy.policy_type,
            ChecksumPolicyType::Required
        );
        assert_eq!(
            opts2.checksum_policy.expected_hash,
            Some("hash123".to_string())
        );
        assert!(!opts2.cache_strategy.check_remote);
        assert!(!opts2.show_progress);
    }

    #[test]
    fn test_download_cache_new_with_custom_dir() {
        let tmp = TempDir::new().unwrap();
        let cache_dir = tmp.path().join("custom_cache");
        let dc = DownloadCache::new(Some(cache_dir.clone())).unwrap();
        assert_eq!(dc.cache_dir, cache_dir);
        assert!(cache_dir.exists());
    }

    #[test]
    fn test_download_cache_manifest_path() {
        let tmp = TempDir::new().unwrap();
        let cache_dir = tmp.path().join("cache/subdir");
        let dc = DownloadCache::new(Some(cache_dir.clone())).unwrap();
        // Manifest path should be set correctly
        assert!(dc.manifest_path.to_string_lossy().contains("manifest.yaml"));
        assert_eq!(dc.cache_dir, cache_dir);
    }

    #[tokio::test]
    async fn cache_hit_by_expected_hash_links_file() {
        let tmp = TempDir::new().unwrap();
        let cache_dir = tmp.path().join("cache");
        std::fs::create_dir_all(&cache_dir).unwrap();

        // Prepare cached file under by-hash
        let hash = "abc123".to_string();
        let by_hash = cache_dir.join("by-hash");
        std::fs::create_dir_all(&by_hash).unwrap();
        let cached_file = by_hash.join(&hash);
        std::fs::write(&cached_file, b"cached").unwrap();

        // Manifest with hash entry (url not needed for this branch)
        let mut m = Manifest::new();
        m.add_download("https://example/file".into(), hash.clone(), 6, None, None);
        write_manifest(tmp.path(), &m);

        let mut dc = DownloadCache::new(Some(cache_dir.clone())).unwrap();

        let target = tmp.path().join("out/target.txt");
        let opts = DownloadOptions {
            checksum_policy: ChecksumPolicy {
                policy_type: ChecksumPolicyType::Required,
                expected_hash: Some(hash.clone()),
            },
            cache_strategy: CacheStrategy { check_remote: true },
            show_progress: false,
            progress_callback: None,
        };

        let res = dc
            .download_with_cache("https://irrelevant", &target, opts)
            .await
            .unwrap();
        assert_eq!(res, target);
        assert_eq!(std::fs::read(&target).unwrap(), b"cached");
    }

    #[tokio::test]
    async fn cache_hit_by_url_without_remote_check_links_file() {
        let tmp = TempDir::new().unwrap();
        let cache_dir = tmp.path().join("cache");
        std::fs::create_dir_all(&cache_dir).unwrap();
        let url = "https://example.com/file".to_string();
        let hash = "deadbeef".to_string();

        // Prepare cached file and manifest entry for URL
        let by_hash = cache_dir.join("by-hash");
        std::fs::create_dir_all(&by_hash).unwrap();
        let cached_file = by_hash.join(&hash);
        std::fs::write(&cached_file, b"via-url").unwrap();

        let mut m = Manifest::new();
        m.add_download(
            url.clone(),
            hash.clone(),
            7,
            Some("etag".into()),
            Some("lm".into()),
        );
        write_manifest(tmp.path(), &m);

        let mut dc = DownloadCache::new(Some(cache_dir.clone())).unwrap();
        let target = tmp.path().join("out2/target.txt");
        let opts = DownloadOptions {
            checksum_policy: ChecksumPolicy {
                policy_type: ChecksumPolicyType::Optional,
                expected_hash: None,
            },
            cache_strategy: CacheStrategy {
                check_remote: false,
            },
            show_progress: false,
            progress_callback: None,
        };

        let res = dc.download_with_cache(&url, &target, opts).await.unwrap();
        assert_eq!(res, target);
        assert_eq!(std::fs::read(&target).unwrap(), b"via-url");
    }

    #[tokio::test]
    #[cfg_attr(not(feature = "slow-tests"), ignore = "slow (network error path)")]
    async fn cache_miss_attempts_download_and_errors() {
        let tmp = TempDir::new().unwrap();
        let cache_dir = tmp.path().join("cache");
        std::fs::create_dir_all(&cache_dir).unwrap();

        // Manifest with a URL entry but missing cache file to force revalidation/download
        let mut m = Manifest::new();
        m.add_download(
            "http://127.0.0.1:9/nonexistent".into(),
            "missinghash".into(),
            0,
            Some("etag".into()),
            Some("lm".into()),
        );
        write_manifest(tmp.path(), &m);

        let mut dc = DownloadCache::new(Some(cache_dir.clone())).unwrap();
        let target = tmp.path().join("out/target.txt");
        let opts = DownloadOptions {
            checksum_policy: ChecksumPolicy {
                policy_type: ChecksumPolicyType::Optional,
                expected_hash: None,
            },
            cache_strategy: CacheStrategy { check_remote: true },
            show_progress: false,
            progress_callback: None,
        };

        // With network restricted and URL unreachable, this should error
        let res = dc
            .download_with_cache("http://127.0.0.1:9/nonexistent", &target, opts)
            .await;
        assert!(res.is_err());
    }
}
