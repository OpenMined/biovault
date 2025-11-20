use crate::config::Config;
use anyhow::{anyhow, Context, Result};
use serde::Serialize;
use serde_json::Value;
use std::fs;
use std::path::{Path, PathBuf};
use syft_crypto_protocol::datasite::context::{ensure_vault_layout, sanitize_identity};
use syft_crypto_protocol::datasite::crypto::{parse_public_bundle, PublicBundleInfo};
use syft_crypto_protocol::identity::generate_identity_material;

const VAULT_DIR_NAME: &str = ".syc";
const CONFIG_DIR: &str = "config";
const DATASITE_JSON: &str = "datasite.json";
const PUBLIC_DID_RELATIVE: &str = "public/crypto/did.json";
const SHADOW_DIR_NAME: &str = "unencrypted";

pub fn vault_path_for_home(home: &Path) -> PathBuf {
    home.parent()
        .map(|parent| parent.join(VAULT_DIR_NAME))
        .unwrap_or_else(|| home.join(VAULT_DIR_NAME))
}

pub fn shadow_root_for_data_root(data_root: &Path) -> PathBuf {
    let base = if data_root
        .file_name()
        .map(|n| n == "datasites")
        .unwrap_or(false)
    {
        data_root.parent().unwrap_or(data_root)
    } else {
        data_root
    };
    base.join(SHADOW_DIR_NAME)
}

pub fn resolve_encrypted_root(data_root: &Path) -> PathBuf {
    if data_root
        .file_name()
        .map(|n| n == "datasites")
        .unwrap_or(false)
    {
        data_root.to_path_buf()
    } else {
        let candidate = data_root.join("datasites");
        if candidate.exists() {
            candidate
        } else {
            data_root.to_path_buf()
        }
    }
}

#[derive(Debug, Clone)]
pub struct IdentityProvisioningOutcome {
    pub identity: String,
    pub generated: bool,
    pub recovery_mnemonic: Option<String>,
    pub vault_path: PathBuf,
    pub bundle_path: PathBuf,
    pub public_bundle_path: PathBuf,
}

#[derive(Serialize)]
struct DatasiteConfigFile<'a> {
    encrypted_root: &'a str,
    shadow_root: &'a str,
}

pub fn provision_local_identity(
    config: &Config,
    biovault_dir: &Path,
) -> Result<IdentityProvisioningOutcome> {
    let data_root = config.get_syftbox_data_dir().context(
        "BioVault init requires a SyftBox data directory. Ensure you're inside a SyftBox environment or set SYFTBOX_DATA_DIR.",
    )?;

    let data_root = data_root
        .canonicalize()
        .unwrap_or_else(|_| data_root.clone());
    let encrypted_root = resolve_encrypted_root(&data_root);
    let shadow_root = shadow_root_for_data_root(&encrypted_root);
    fs::create_dir_all(&shadow_root).with_context(|| {
        format!(
            "failed to create shadow directory: {}",
            shadow_root.display()
        )
    })?;

    let vault_path = vault_path_for_home(biovault_dir);
    ensure_vault_layout(&vault_path).map_err(|err| {
        anyhow!(
            "failed to prepare Syft Crypto vault at {}: {err}",
            vault_path.display()
        )
    })?;

    write_datasite_config(&vault_path, &encrypted_root, &shadow_root)?;

    let identity = config.email.trim();
    let slug = sanitize_identity(identity);
    let key_path = vault_path.join("keys").join(format!("{slug}.key"));
    let bundle_path = vault_path.join("bundles").join(format!("{slug}.json"));

    let mut generated = false;
    let mut recovery_mnemonic = None;
    let mut public_bundle_value = None;

    if !key_path.exists() || !bundle_path.exists() {
        let generated_identity = generate_identity_material(identity)?;
        fs::write(&key_path, &generated_identity.key_file)
            .with_context(|| format!("failed to write private key file: {}", key_path.display()))?;
        fs::write(
            &bundle_path,
            serde_json::to_vec_pretty(&generated_identity.public_bundle)?,
        )
        .with_context(|| format!("failed to write bundle file: {}", bundle_path.display()))?;
        generated = true;
        recovery_mnemonic = Some(generated_identity.recovery_key_mnemonic.clone());
        public_bundle_value = Some(generated_identity.public_bundle);
    }

    let public_bundle = match public_bundle_value {
        Some(bundle) => bundle,
        None => {
            let contents = fs::read_to_string(&bundle_path).with_context(|| {
                format!("failed to read bundle file: {}", bundle_path.display())
            })?;
            serde_json::from_str(&contents).context("failed to parse existing bundle JSON")?
        }
    };

    let public_bundle_path = export_public_bundle(identity, &public_bundle, &encrypted_root)?;

    Ok(IdentityProvisioningOutcome {
        identity: identity.to_string(),
        generated,
        recovery_mnemonic,
        vault_path,
        bundle_path,
        public_bundle_path,
    })
}

pub fn import_public_bundle(
    config: &Config,
    bundle_path: &Path,
    expected_identity: Option<&str>,
) -> Result<PublicBundleInfo> {
    let biovault_home = crate::config::get_biovault_home()?;
    let vault_path = vault_path_for_home(&biovault_home);
    ensure_vault_layout(&vault_path)
        .map_err(|err| anyhow!("failed to prepare Syft Crypto vault: {err}"))?;

    let bundle_bytes = fs::read(bundle_path)
        .with_context(|| format!("failed to read bundle at {}", bundle_path.display()))?;
    let bundle_str = String::from_utf8(bundle_bytes)
        .with_context(|| format!("bundle at {} is not valid UTF-8", bundle_path.display()))?;
    let bundle: Value = serde_json::from_str(&bundle_str)
        .with_context(|| format!("failed to parse bundle JSON at {}", bundle_path.display()))?;
    let parsed = parse_public_bundle_from_str(&bundle_str)?;

    if let Some(expected) = expected_identity {
        if parsed.identity != expected {
            return Err(anyhow!(
                "bundle identity mismatch: expected {}, found {}",
                expected,
                parsed.identity
            ));
        }
    }

    let slug = sanitize_identity(&parsed.identity);
    let target = vault_path.join("bundles").join(format!("{slug}.json"));
    fs::create_dir_all(target.parent().unwrap())?;
    fs::write(&target, serde_json::to_vec_pretty(&bundle)?)
        .with_context(|| format!("failed to write bundle to {}", target.display()))?;

    // If this is our own identity, ensure public bundle copy is refreshed.
    if parsed.identity == config.email {
        if let Ok(data_root) = config.get_syftbox_data_dir() {
            let _ = export_public_bundle(&parsed.identity, &bundle, &data_root);
        }
    }

    Ok(parsed)
}

fn write_datasite_config(
    vault_path: &Path,
    encrypted_root: &Path,
    shadow_root: &Path,
) -> Result<()> {
    let config_dir = vault_path.join(CONFIG_DIR);
    fs::create_dir_all(&config_dir)?;
    let json = DatasiteConfigFile {
        encrypted_root: &encrypted_root.to_string_lossy(),
        shadow_root: &shadow_root.to_string_lossy(),
    };
    let payload = serde_json::to_string_pretty(&json)?;
    fs::write(config_dir.join(DATASITE_JSON), payload)?;
    Ok(())
}

fn export_public_bundle(identity: &str, bundle: &Value, data_root: &Path) -> Result<PathBuf> {
    let public_dir = data_root.join(identity).join(PUBLIC_DID_RELATIVE);
    if let Some(parent) = public_dir.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create directory: {}", parent.display()))?;
    }
    fs::write(&public_dir, serde_json::to_vec_pretty(bundle)?)
        .with_context(|| format!("failed to export bundle to {}", public_dir.display()))?;
    Ok(public_dir)
}

fn parse_public_bundle_from_str(body: &str) -> Result<PublicBundleInfo> {
    parse_public_bundle(body).map_err(|err| anyhow!("invalid bundle: {err}"))
}
