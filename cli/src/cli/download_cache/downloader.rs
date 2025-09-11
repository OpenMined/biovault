use super::manifest::Manifest;
use super::{calculate_blake3, link_or_copy, CacheStrategy, ChecksumPolicy, ChecksumPolicyType};
use anyhow::{anyhow, Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::header::{ETAG, LAST_MODIFIED};
use std::fs;
use std::path::{Path, PathBuf};
use tokio::fs::File;
use tokio::io::AsyncWriteExt;

pub struct DownloadCache {
    cache_dir: PathBuf,
    manifest_path: PathBuf,
    manifest: Manifest,
}

#[derive(Debug, Clone)]
pub struct DownloadOptions {
    pub checksum_policy: ChecksumPolicy,
    pub cache_strategy: CacheStrategy,
    pub show_progress: bool,
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
        println!("📦 Processing: {}", url);

        // Check if we have the expected hash in cache
        if let Some(ref expected_hash) = options.checksum_policy.expected_hash {
            if let Some(_hash_entry) = self.manifest.get_by_hash(expected_hash) {
                let cache_path = self.cache_dir.join("by-hash").join(expected_hash);
                if cache_path.exists() {
                    println!("  ✓ Found in cache (hash: {}...)", &expected_hash[..8]);
                    println!("  ✓ Cache hit! Using cached version");

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
                        println!("  ↻ Checking if remote file has changed...");

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
                            println!("  ✓ Remote file unchanged");
                            println!(
                                "  ✓ Cache hit! Using cached version (hash: {}...)",
                                &current_hash[..8]
                            );

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
                        println!(
                            "  ✓ Cache hit! Using cached version (hash: {}...)",
                            &current_hash[..8]
                        );

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
            println!("  ↓ Cache miss - downloading file...");

            // Download to temporary location
            let temp_path = self
                .cache_dir
                .join("downloads")
                .join(format!("{}.tmp", uuid::Uuid::new_v4()));
            fs::create_dir_all(temp_path.parent().unwrap())?;

            let (etag, last_modified) = self
                .download_file(url, &temp_path, options.show_progress)
                .await?;

            // Calculate hash
            print!("  ⟳ Computing BLAKE3 checksum... ");
            std::io::Write::flush(&mut std::io::stdout())?;
            let actual_hash = calculate_blake3(&temp_path)?;
            println!("done ({}...)", &actual_hash[..8]);

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
                    println!("  ✓ Checksum verified");
                }
                ChecksumPolicyType::Preferred => {
                    if let Some(expected) = &options.checksum_policy.expected_hash {
                        if &actual_hash != expected {
                            println!("  ⚠ Warning: Checksum mismatch!");
                            println!("    Expected: {}...", &expected[..8]);
                            println!("    Got:      {}...", &actual_hash[..8]);
                            println!("  ⚠ Continuing anyway (policy: preferred)");
                        } else {
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
                println!("  ✓ File already in cache with same content");
                fs::remove_file(&temp_path)?;
            } else {
                fs::rename(&temp_path, &cache_path)?;
                println!("  ✓ Added to cache");
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
            println!("  ✓ File ready at: {}", target_path.display());

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
    ) -> Result<(Option<String>, Option<String>)> {
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(3600))
            .build()?;

        let response = client.get(url).send().await?;

        if !response.status().is_success() {
            return Err(anyhow!("HTTP request failed: {}", response.status()));
        }

        // Extract headers
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

        let total_size = response.content_length().unwrap_or(0);

        let pb = if show_progress && total_size > 0 {
            let pb = ProgressBar::new(total_size);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("    [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                    .expect("Failed to set progress bar template")
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else if show_progress {
            println!("    Downloading (size unknown)...");
            None
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
        }

        if let Some(pb) = pb {
            pb.finish_and_clear();
        }

        Ok((etag, last_modified))
    }

    fn save_manifest(&self) -> Result<()> {
        self.manifest.save(&self.manifest_path)
    }
}
