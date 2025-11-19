use crate::config;
use crate::syftbox::syc;
use anyhow::{anyhow, bail, Context, Result};
use serde::de::DeserializeOwned;
use serde::Serialize;
use serde_json;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use syft_crypto_protocol::datasite::{
    bytes::{read_bytes, write_bytes, BytesReadOpts, BytesWriteOpts, BytesWriteOutcome},
    context::{ensure_vault_layout, resolve_vault, AppContext},
};
use tracing::warn;
use walkdir::WalkDir;

#[derive(Debug, Clone)]
pub struct SyftBoxStorage {
    root: PathBuf,
    backend: StorageBackend,
}

#[derive(Debug, Clone, Copy)]
pub enum ReadPolicy {
    AllowPlaintext,
    RequireEnvelope,
}

#[derive(Debug, Clone)]
pub enum WritePolicy {
    Plaintext,
    Envelope {
        recipients: Vec<String>,
        hint: Option<String>,
    },
}

#[derive(Debug, Clone)]
enum StorageBackend {
    SyctCrypto(SyctCryptoBackend),
    PlainFs,
}

#[derive(Debug, Clone)]
struct SyctCryptoBackend {
    context: AppContext,
}

impl SyftBoxStorage {
    pub fn new(root: &Path) -> Self {
        let disable_crypto = env::var("BIOVAULT_DISABLE_SYC")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
            || cfg!(test);

        let canonical_root = root.canonicalize().unwrap_or_else(|_| root.to_path_buf());
        let encrypted_root = syc::resolve_encrypted_root(&canonical_root);

        let backend = if disable_crypto {
            StorageBackend::PlainFs
        } else {
            SyctCryptoBackend::new(&encrypted_root)
                .map(StorageBackend::SyctCrypto)
                .unwrap_or_else(|err| {
                    warn!(
                        "syc backend unavailable ({}); falling back to plaintext filesystem access",
                        err
                    );
                    StorageBackend::PlainFs
                })
        };

        let uses_crypto = matches!(backend, StorageBackend::SyctCrypto(_));
        eprintln!(
            "SyftBoxStorage initialized: root={:?}, crypto={}",
            encrypted_root, uses_crypto
        );

        Self {
            root: encrypted_root,
            backend,
        }
    }

    pub fn uses_crypto(&self) -> bool {
        matches!(self.backend, StorageBackend::SyctCrypto(_))
    }

    pub fn write_plaintext_file(
        &self,
        absolute_path: &Path,
        data: &[u8],
        overwrite: bool,
    ) -> Result<()> {
        self.write_with_policy(absolute_path, data, WritePolicy::Plaintext, overwrite)?;
        Ok(())
    }

    pub fn write_encrypted_file(
        &self,
        absolute_path: &Path,
        data: &[u8],
        recipients: Vec<String>,
        hint: Option<String>,
        overwrite: bool,
    ) -> Result<()> {
        self.write_with_policy(
            absolute_path,
            data,
            WritePolicy::Envelope { recipients, hint },
            overwrite,
        )?;
        Ok(())
    }

    /// Write encrypted file using shadow folder pattern:
    /// 1. Write plaintext to shadow folder
    /// 2. Encrypt from shadow → datasites
    pub fn write_encrypted_with_shadow(
        &self,
        datasite_path: &Path,
        data: &[u8],
        recipients: Vec<String>,
        hint: Option<String>,
        overwrite: bool,
    ) -> Result<()> {
        match &self.backend {
            StorageBackend::SyctCrypto(backend) => {
                // Get relative path from datasite root
                let relative = self.relative_from_root(datasite_path)?;

                // Write plaintext to shadow folder
                use syft_crypto_protocol::datasite::context::resolve_shadow_path;
                let shadow_path = resolve_shadow_path(&backend.context, &relative);
                if let Some(parent) = shadow_path.parent() {
                    fs::create_dir_all(parent)
                        .with_context(|| format!("failed to create shadow parent {:?}", parent))?;
                }
                fs::write(&shadow_path, data).with_context(|| {
                    format!("failed to write plaintext to shadow {:?}", shadow_path)
                })?;

                // Encrypt from shadow to datasites using syc file encrypt pattern
                use syft_crypto_protocol::datasite::context::resolve_identity;
                use syft_crypto_protocol::datasite::crypto::{
                    encrypt_envelope_for_recipient, load_private_keys_for_identity,
                    resolve_recipient_bundle,
                };

                if recipients.len() != 1 {
                    bail!("exactly one recipient is supported for encryption");
                }

                let sender_identity = resolve_identity(None, &backend.context.vault_path)?;
                let recipient_identity = recipients[0].clone();
                let sender_keys =
                    load_private_keys_for_identity(&backend.context, &sender_identity)?;
                let recipient_bundle = resolve_recipient_bundle(
                    &backend.context,
                    &sender_keys,
                    &sender_identity,
                    &recipient_identity,
                )?;

                let plaintext = fs::read(&shadow_path)?;
                let envelope = encrypt_envelope_for_recipient(
                    &sender_identity,
                    &sender_keys,
                    &recipient_identity,
                    &recipient_bundle,
                    &plaintext,
                    hint.as_deref(),
                )?;

                // Write encrypted to datasites
                use syft_crypto_protocol::datasite::context::atomic_write;
                atomic_write(datasite_path, &envelope)?;

                Ok(())
            }
            StorageBackend::PlainFs => {
                // Fallback to direct write for PlainFs backend
                self.write_plaintext_file(datasite_path, data, overwrite)
            }
        }
    }

    pub fn read_plaintext_file(&self, absolute_path: &Path) -> Result<Vec<u8>> {
        self.read_with_policy(absolute_path, ReadPolicy::AllowPlaintext)
    }

    /// Read encrypted file using shadow folder pattern:
    /// 1. Decrypt from datasites → shadow
    /// 2. Read plaintext from shadow
    pub fn read_with_shadow(&self, datasite_path: &Path) -> Result<Vec<u8>> {
        match &self.backend {
            StorageBackend::SyctCrypto(backend) => {
                // Get relative path from datasite root
                let relative = self.relative_from_root(datasite_path)?;
                use syft_crypto_protocol::datasite::context::resolve_shadow_path;
                let shadow_path = resolve_shadow_path(&backend.context, &relative);
                eprintln!("DEBUG read_with_shadow:");
                eprintln!("  datasite: {:?}", datasite_path);
                eprintln!("  relative: {:?}", relative);
                eprintln!("  shadow: {:?}", shadow_path);
                eprintln!("  shadow_root: {:?}", backend.context.shadow_root);

                // If shadow copy exists and is newer, use it
                if shadow_path.exists() && datasite_path.exists() {
                    let shadow_meta = fs::metadata(&shadow_path)?;
                    let datasite_meta = fs::metadata(datasite_path)?;
                    if shadow_meta.modified()? >= datasite_meta.modified()? {
                        return fs::read(&shadow_path).with_context(|| {
                            format!("failed to read from shadow {:?}", shadow_path)
                        });
                    }
                }

                // Decrypt from datasites to shadow using syc file decrypt pattern
                use syft_crypto_protocol::datasite::context::resolve_identity;
                use syft_crypto_protocol::datasite::crypto::{
                    decrypt_envelope_for_recipient, load_private_keys_for_identity,
                    parse_optional_envelope, resolve_sender_bundle_for_decrypt,
                };

                let bytes = fs::read(datasite_path)?;

                // Try to decrypt - ONLY cache if we get plaintext
                let envelope = parse_optional_envelope(&bytes)?;

                let plaintext = if let Some(envelope) = envelope {
                    // File is encrypted, MUST decrypt successfully to cache
                    let identity = resolve_identity(None, &backend.context.vault_path)?;
                    let recipient_keys =
                        load_private_keys_for_identity(&backend.context, &identity)?;
                    let sender_bundle =
                        resolve_sender_bundle_for_decrypt(&backend.context, &envelope)?;
                    decrypt_envelope_for_recipient(
                        &identity,
                        &recipient_keys,
                        &sender_bundle,
                        &envelope,
                    )?
                } else {
                    // File is plaintext
                    bytes
                };

                // Cache ONLY successfully decrypted plaintext
                if let Some(parent) = shadow_path.parent() {
                    match fs::create_dir_all(parent) {
                        Ok(_) => eprintln!("  ✓ Created shadow parent: {:?}", parent),
                        Err(e) => {
                            eprintln!("  ✗ Failed to create shadow parent: {:?}: {}", parent, e)
                        }
                    }
                }
                match fs::write(&shadow_path, &plaintext) {
                    Ok(_) => {
                        eprintln!("✓ Cached PLAINTEXT to shadow: {:?}", shadow_path);
                        // Sanity check: ensure we didn't cache encrypted data
                        if plaintext.len() >= 4 && &plaintext[0..4] == b"SYC1" {
                            panic!(
                                "BUG: Cached encrypted data to shadow folder! Path: {:?}",
                                shadow_path
                            );
                        }
                    }
                    Err(e) => eprintln!("✗ Failed to cache shadow: {:?}: {}", shadow_path, e),
                }

                Ok(plaintext)
            }
            StorageBackend::PlainFs => {
                // Fallback to direct read for PlainFs backend
                self.read_plaintext_file(datasite_path)
            }
        }
    }

    pub fn read_plaintext_string(&self, absolute_path: &Path) -> Result<String> {
        let bytes = self.read_plaintext_file(absolute_path)?;
        let content = String::from_utf8(bytes)
            .with_context(|| format!("failed to decode utf-8 from {:?}", absolute_path))?;
        Ok(content)
    }

    pub fn read_json<T: DeserializeOwned>(
        &self,
        absolute_path: &Path,
        policy: ReadPolicy,
    ) -> Result<T> {
        let bytes = self.read_with_policy(absolute_path, policy)?;
        serde_json::from_slice(&bytes)
            .with_context(|| format!("failed to parse JSON from {:?}", absolute_path))
    }

    /// Read JSON using shadow folder pattern for encrypted files
    pub fn read_json_with_shadow<T: DeserializeOwned>(&self, absolute_path: &Path) -> Result<T> {
        let bytes = self.read_with_shadow(absolute_path)?;
        serde_json::from_slice(&bytes)
            .with_context(|| format!("failed to parse JSON from {:?}", absolute_path))
    }

    /// Write with shadow folder pattern - creates both encrypted file and plaintext shadow
    /// This makes unencrypted/ a true mirror of datasites/
    pub fn write_with_shadow(
        &self,
        absolute_path: &Path,
        plaintext: &[u8],
        policy: WritePolicy,
        overwrite: bool,
    ) -> Result<BytesWriteOutcome> {
        // Check if we'll need to create a shadow (before policy is moved)
        let needs_shadow = matches!(&policy, WritePolicy::Envelope { .. });

        // First write the encrypted file (or plaintext if policy is Plaintext)
        let outcome = self.write_with_policy(absolute_path, plaintext, policy, overwrite)?;

        // If crypto is enabled and we wrote an encrypted file, also cache the plaintext to shadow
        if let StorageBackend::SyctCrypto(ref backend) = self.backend {
            if needs_shadow {
                // Calculate shadow path
                if let Ok(relative) = absolute_path.strip_prefix(&self.root) {
                    let shadow_path = backend.context.shadow_root.join(relative);

                    eprintln!("DEBUG write_with_shadow:");
                    eprintln!("  datasite: {:?}", absolute_path);
                    eprintln!("  relative: {:?}", relative);
                    eprintln!("  shadow: {:?}", shadow_path);
                    eprintln!("  shadow_root: {:?}", backend.context.shadow_root);

                    // Create parent directory for shadow
                    if let Some(shadow_parent) = shadow_path.parent() {
                        if let Err(e) = fs::create_dir_all(shadow_parent) {
                            eprintln!(
                                "✗ Failed to create shadow parent: {:?}: {}",
                                shadow_parent, e
                            );
                        } else {
                            eprintln!("  ✓ Created shadow parent: {:?}", shadow_parent);
                        }
                    }

                    // Write plaintext to shadow
                    match fs::write(&shadow_path, plaintext) {
                        Ok(_) => {
                            eprintln!("✓ Cached PLAINTEXT to shadow (write): {:?}", shadow_path);
                            // Sanity check: ensure we didn't cache encrypted data
                            if plaintext.len() >= 4 && &plaintext[0..4] == b"SYC1" {
                                panic!("BUG: Tried to cache encrypted data to shadow folder during write! Path: {:?}", shadow_path);
                            }
                        }
                        Err(e) => {
                            eprintln!("✗ Failed to cache shadow (write): {:?}: {}", shadow_path, e)
                        }
                    }
                }
            }
        }

        Ok(outcome)
    }

    /// Write JSON with shadow folder pattern - creates both encrypted file and plaintext shadow
    pub fn write_json_with_shadow<T: Serialize>(
        &self,
        absolute_path: &Path,
        value: &T,
        policy: WritePolicy,
        overwrite: bool,
    ) -> Result<()> {
        let data = serde_json::to_vec_pretty(value)
            .with_context(|| format!("failed to serialize JSON for {:?}", absolute_path))?;
        self.write_with_shadow(absolute_path, &data, policy, overwrite)?;
        Ok(())
    }

    /// Write JSON using the old method (no shadow) - deprecated, use write_json_with_shadow instead
    pub fn write_json<T: Serialize>(
        &self,
        absolute_path: &Path,
        value: &T,
        policy: WritePolicy,
        overwrite: bool,
    ) -> Result<()> {
        let data = serde_json::to_vec_pretty(value)
            .with_context(|| format!("failed to serialize JSON for {:?}", absolute_path))?;
        self.write_with_policy(absolute_path, &data, policy, overwrite)?;
        Ok(())
    }

    pub fn remove_path(&self, absolute_path: &Path) -> Result<()> {
        self.ensure_within_root(absolute_path)?;

        // Remove the main file/directory if it exists
        if absolute_path.exists() {
            if absolute_path.is_dir() {
                fs::remove_dir_all(absolute_path)
                    .with_context(|| format!("failed to remove directory {:?}", absolute_path))?;
            } else {
                fs::remove_file(absolute_path)
                    .with_context(|| format!("failed to remove file {:?}", absolute_path))?;
            }
        }

        // If crypto is enabled, also remove the corresponding shadow file/directory
        if let StorageBackend::SyctCrypto(ref backend) = self.backend {
            // Calculate shadow path: strip datasites root, add to shadow root
            if let Ok(relative) = absolute_path.strip_prefix(&self.root) {
                let shadow_path = backend.context.shadow_root.join(relative);

                if shadow_path.exists() {
                    if shadow_path.is_dir() {
                        fs::remove_dir_all(&shadow_path).with_context(|| {
                            format!("failed to remove shadow directory {:?}", shadow_path)
                        })?;
                        eprintln!("✓ Removed shadow directory: {:?}", shadow_path);
                    } else {
                        fs::remove_file(&shadow_path).with_context(|| {
                            format!("failed to remove shadow file {:?}", shadow_path)
                        })?;
                        eprintln!("✓ Removed shadow file: {:?}", shadow_path);
                    }
                }
            }
        }

        Ok(())
    }

    pub fn ensure_dir(&self, dir: &Path) -> Result<()> {
        self.ensure_within_root(dir)?;
        fs::create_dir_all(dir).with_context(|| format!("failed to ensure directory {:?}", dir))?;
        Ok(())
    }

    pub fn contains(&self, path: &Path) -> bool {
        path.starts_with(&self.root)
    }

    pub fn list_dir(&self, dir: &Path) -> Result<Vec<PathBuf>> {
        self.ensure_within_root(dir)?;
        if !dir.exists() {
            return Ok(Vec::new());
        }
        let mut entries = Vec::new();
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            entries.push(entry.path());
        }
        Ok(entries)
    }

    pub fn path_exists(&self, absolute_path: &Path) -> Result<bool> {
        self.ensure_within_root(absolute_path)?;
        Ok(absolute_path.exists())
    }

    pub fn copy_raw_file(&self, src: &Path, dst: &Path) -> Result<()> {
        self.ensure_within_root(src)?;
        self.ensure_within_root(dst)?;
        if let Some(parent) = dst.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to ensure parent {:?}", parent))?;
        }
        fs::copy(src, dst).with_context(|| format!("failed to copy {:?} -> {:?}", src, dst))?;
        Ok(())
    }

    pub fn copy_tree<F>(
        &self,
        src: &Path,
        dst: &Path,
        mut skip: F,
        write_policy: WritePolicy,
    ) -> Result<()>
    where
        F: FnMut(&Path) -> bool,
    {
        self.ensure_within_root(src)?;
        self.ensure_within_root(dst)?;
        self.ensure_dir(dst)?;

        for entry in WalkDir::new(src).into_iter().filter_map(|e| e.ok()) {
            let rel = entry
                .path()
                .strip_prefix(src)
                .with_context(|| format!("failed to relativize {:?}", entry.path()))?;
            let target = dst.join(rel);

            if entry.file_type().is_dir() {
                self.ensure_dir(&target)?;
                continue;
            }

            if skip(entry.path()) {
                continue;
            }

            if let Some(parent) = target.parent() {
                self.ensure_dir(parent)?;
            }

            let data = self.read_plaintext_file(entry.path())?;
            self.write_with_policy(&target, &data, write_policy.clone(), true)?;
        }

        Ok(())
    }

    fn write_with_policy(
        &self,
        absolute_path: &Path,
        data: &[u8],
        policy: WritePolicy,
        overwrite: bool,
    ) -> Result<BytesWriteOutcome> {
        let relative = self.relative_from_root(absolute_path)?;
        let (recipients, plaintext, hint) = match policy {
            WritePolicy::Plaintext => (Vec::new(), true, None),
            WritePolicy::Envelope { recipients, hint } => (recipients, false, hint),
        };
        let opts = BytesWriteOpts {
            relative,
            recipients,
            plaintext,
            overwrite,
            hint,
        };
        match &self.backend {
            StorageBackend::SyctCrypto(backend) => write_bytes(&backend.context, &opts, data)
                .map_err(|err| anyhow!("syc write failed: {err}")),
            StorageBackend::PlainFs => {
                if absolute_path.exists() && !overwrite {
                    bail!(
                        "refusing to overwrite existing path {:?} without explicit permission",
                        absolute_path
                    );
                }
                if let Some(parent) = absolute_path.parent() {
                    fs::create_dir_all(parent)
                        .with_context(|| format!("failed to ensure parent {:?}", parent))?;
                }
                fs::write(absolute_path, data)
                    .with_context(|| format!("failed to write {:?}", absolute_path))?;
                Ok(BytesWriteOutcome {
                    destination: absolute_path.to_path_buf(),
                    bytes_written: data.len(),
                    encrypted: false,
                })
            }
        }
    }

    fn read_with_policy(&self, absolute_path: &Path, policy: ReadPolicy) -> Result<Vec<u8>> {
        match &self.backend {
            StorageBackend::SyctCrypto(backend) => {
                let relative = self.relative_from_root(absolute_path)?;
                let opts = BytesReadOpts {
                    relative,
                    identity: None,
                    require_envelope: matches!(policy, ReadPolicy::RequireEnvelope),
                };
                let result = read_bytes(&backend.context, &opts)
                    .map_err(|err| anyhow!("syc read failed: {err}"))?;
                Ok(result.plaintext)
            }
            StorageBackend::PlainFs => fs::read(absolute_path)
                .with_context(|| format!("failed to read {:?}", absolute_path)),
        }
    }

    /// Canonicalize path to handle symlinks (e.g., /var -> /private/var on macOS)
    /// For non-existent paths, find first existing ancestor, canonicalize it, and rebuild path
    fn canonicalize_for_comparison(&self, path: &Path) -> PathBuf {
        if path.exists() {
            return path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
        }

        // Find first existing ancestor
        let mut current = path;
        let mut components = Vec::new();

        while !current.exists() {
            if let Some(name) = current.file_name() {
                components.push(name);
            }
            if let Some(parent) = current.parent() {
                current = parent;
            } else {
                // No existing ancestor found, return as-is
                return path.to_path_buf();
            }
        }

        // Canonicalize the existing ancestor
        let canonical_base = current
            .canonicalize()
            .unwrap_or_else(|_| current.to_path_buf());

        // Rebuild path from canonical base
        components.reverse();
        let mut result = canonical_base;
        for component in components {
            result = result.join(component);
        }
        result
    }

    fn relative_from_root(&self, absolute: &Path) -> Result<PathBuf> {
        let canonical_absolute = self.canonicalize_for_comparison(absolute);

        canonical_absolute
            .strip_prefix(&self.root)
            .map(|p| p.to_path_buf())
            .map_err(|_| {
                anyhow!(
                    "path {:?} is outside of SyftBox root {:?}",
                    canonical_absolute,
                    self.root
                )
            })
    }

    fn ensure_within_root(&self, absolute: &Path) -> Result<()> {
        let canonical_absolute = self.canonicalize_for_comparison(absolute);

        if canonical_absolute.starts_with(&self.root) {
            Ok(())
        } else {
            Err(anyhow!(
                "path {:?} is outside of SyftBox root {:?}",
                canonical_absolute,
                self.root
            ))
        }
    }
}

impl SyctCryptoBackend {
    fn new(root: &Path) -> Result<Self> {
        let vault_override = env::var_os("SYC_VAULT").map(PathBuf::from);
        let default_vault = {
            let home = config::get_biovault_home()
                .context("failed to determine BioVault home for Syft Crypto vault")?;
            syc::vault_path_for_home(&home)
        };
        let vault_path = resolve_vault(Some(
            vault_override.unwrap_or_else(|| default_vault.clone()),
        ));
        ensure_vault_layout(&vault_path)
            .map_err(|err| anyhow!("failed to prepare syc vault: {err}"))?;

        let data_root = root.canonicalize().unwrap_or_else(|_| root.to_path_buf());
        let encrypted_root = syc::resolve_encrypted_root(&data_root);
        let shadow_root = syc::shadow_root_for_data_root(&encrypted_root);
        fs::create_dir_all(&shadow_root)
            .with_context(|| format!("failed to prepare shadow root {:?}", shadow_root))?;

        Ok(Self {
            context: AppContext {
                vault_path,
                data_root: encrypted_root,
                shadow_root,
            },
        })
    }
}
