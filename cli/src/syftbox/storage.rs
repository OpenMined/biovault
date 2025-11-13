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

        let backend = if disable_crypto {
            StorageBackend::PlainFs
        } else {
            SyctCryptoBackend::new(root)
                .map(StorageBackend::SyctCrypto)
                .unwrap_or_else(|err| {
                    warn!(
                        "syc backend unavailable ({}); falling back to plaintext filesystem access",
                        err
                    );
                    StorageBackend::PlainFs
                })
        };

        Self {
            root: root.to_path_buf(),
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

    pub fn read_plaintext_file(&self, absolute_path: &Path) -> Result<Vec<u8>> {
        self.read_with_policy(absolute_path, ReadPolicy::AllowPlaintext)
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
        if absolute_path.exists() {
            if absolute_path.is_dir() {
                fs::remove_dir_all(absolute_path)
                    .with_context(|| format!("failed to remove directory {:?}", absolute_path))?;
            } else {
                fs::remove_file(absolute_path)
                    .with_context(|| format!("failed to remove file {:?}", absolute_path))?;
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
        let relative = self.relative_from_root(absolute_path)?;
        match &self.backend {
            StorageBackend::SyctCrypto(backend) => {
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

    fn relative_from_root(&self, absolute: &Path) -> Result<PathBuf> {
        absolute
            .strip_prefix(&self.root)
            .map(|p| p.to_path_buf())
            .map_err(|_| {
                anyhow!(
                    "path {:?} is outside of SyftBox root {:?}",
                    absolute,
                    self.root
                )
            })
    }

    fn ensure_within_root(&self, absolute: &Path) -> Result<()> {
        if absolute.starts_with(&self.root) {
            Ok(())
        } else {
            Err(anyhow!(
                "path {:?} is outside of SyftBox root {:?}",
                absolute,
                self.root
            ))
        }
    }
}

impl SyctCryptoBackend {
    fn new(root: &Path) -> Result<Self> {
        let vault_path = resolve_vault(None);
        ensure_vault_layout(&vault_path)
            .map_err(|err| anyhow!("failed to prepare syc vault: {err}"))?;

        let data_root = root.canonicalize().unwrap_or_else(|_| root.to_path_buf());
        let shadow_root = root.join(".syc-shadow");
        fs::create_dir_all(&shadow_root)
            .with_context(|| format!("failed to prepare shadow root {:?}", shadow_root))?;

        Ok(Self {
            context: AppContext {
                vault_path,
                data_root,
                shadow_root,
            },
        })
    }
}
