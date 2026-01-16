use crate::config::Config;
use crate::syftbox::syc::{
    import_public_bundle, parse_public_bundle_file, provision_local_identity_with_options,
    restore_identity_from_mnemonic,
};
use crate::Result;
use anyhow::{anyhow, Context};
use clap::Subcommand;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};
use syftbox_sdk::PublicBundleInfo;

#[derive(Subcommand, Debug)]
pub enum KeyCommands {
    /// Generate a key for an identity if missing; prints recovery mnemonic once.
    Generate {
        #[arg(long, help = "Identity email (defaults to config email)")]
        email: Option<String>,
        #[arg(
            long,
            help = "Override SyftBox data dir (defaults to BioVault SyftBox config data_dir)"
        )]
        data_dir: Option<PathBuf>,
        #[arg(
            long,
            help = "Override Syft Crypto vault path (defaults to datasites/.syc)"
        )]
        vault: Option<PathBuf>,
        #[arg(long, help = "Force overwrite existing key material")]
        force: bool,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// Import a public bundle into the vault (TOFU-checked).
    Import {
        #[arg(help = "Path to public bundle (e.g., datasites/<email>/public/crypto/did.json)")]
        bundle: PathBuf,
        #[arg(
            long,
            help = "Expected identity inside bundle; default inferred from bundle"
        )]
        email: Option<String>,
        #[arg(long, help = "Skip TOFU fingerprint check if existing bundle differs")]
        ignore_tofu: bool,
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Override data dir (used for exported copy refresh)")]
        data_dir: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// Export the public bundle for an identity.
    Export {
        #[arg(long, help = "Identity email (defaults to config email)")]
        email: Option<String>,
        #[arg(
            long,
            help = "Output path (defaults to datasites/<id>/public/crypto/did.json)"
        )]
        output: Option<PathBuf>,
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Override data dir")]
        data_dir: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// Restore key material from a BIP-39 mnemonic.
    Restore {
        #[arg(long, help = "Identity email")]
        email: String,
        #[arg(long, help = "BIP-39 recovery phrase")]
        mnemonic: String,
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Override data dir")]
        data_dir: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// Wipe key and bundle for an identity.
    Wipe {
        #[arg(long, help = "Identity email (defaults to config email)")]
        email: Option<String>,
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Override data dir")]
        data_dir: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// Show key/bundle status for an identity.
    Status {
        #[arg(long, help = "Identity email (defaults to config email)")]
        email: Option<String>,
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Override data dir")]
        data_dir: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
    /// List identities present in the vault.
    VaultList {
        #[arg(long, help = "Override vault path")]
        vault: Option<PathBuf>,
        #[arg(long, help = "Output JSON instead of human-readable text")]
        json: bool,
    },
}

pub async fn handle(command: KeyCommands, config: &Config) -> Result<()> {
    match command {
        KeyCommands::Generate {
            email,
            data_dir,
            vault,
            force,
            json,
        } => {
            let email = resolve_email(email.as_deref(), config)?;
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;
            let outcome =
                provision_local_identity_with_options(&email, &data_root, Some(&vault_path), force)
                    .with_context(|| format!("failed to generate identity {email}"))?;
            let bundle = parse_public_bundle_file(&outcome.public_bundle_path)?;
            let result = GenerateResult {
                identity: bundle.identity.clone(),
                fingerprint: bundle.fingerprint.clone(),
                vault: outcome.vault_path.to_string_lossy().to_string(),
                bundle_path: outcome.bundle_path.to_string_lossy().to_string(),
                export_path: outcome.public_bundle_path.to_string_lossy().to_string(),
                mnemonic: outcome.recovery_mnemonic.clone(),
            };
            if json {
                print_json(&result)?;
            } else {
                println!("✓ Identity: {}", bundle.identity);
                println!("  Fingerprint: {}", bundle.fingerprint);
                println!("  Vault: {}", outcome.vault_path.display());
                println!("  Bundle: {}", outcome.bundle_path.display());
                println!("  Exported: {}", outcome.public_bundle_path.display());
                if let Some(mnemonic) = outcome.recovery_mnemonic {
                    println!("  Recovery mnemonic (store securely, shown once): {mnemonic}");
                }
            }
        }
        KeyCommands::Import {
            bundle,
            email,
            ignore_tofu,
            vault,
            data_dir,
            json,
        } => {
            let info = parse_public_bundle_file(&bundle)?;
            let expected_identity = email.as_deref().unwrap_or(&info.identity);
            if info.identity != expected_identity {
                return Err(anyhow!(
                    "bundle identity mismatch: expected {}, found {}",
                    expected_identity,
                    info.identity
                )
                .into());
            }
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;

            if let Some(existing) = load_existing_bundle(&vault_path, &info.identity)? {
                if existing.fingerprint != info.fingerprint && !ignore_tofu {
                    return Err(anyhow!(
                        "TOFU violation: bundle for {} already present with fingerprint {} (incoming {}) – rerun with --ignore-tofu to overwrite",
                        info.identity,
                        existing.fingerprint,
                        info.fingerprint
                    )
                    .into());
                }
            }

            let parsed = import_public_bundle(
                &bundle,
                Some(expected_identity),
                &vault_path,
                Some(&data_root),
                Some(config.email.as_str()),
            )
            .with_context(|| format!("failed to import bundle from {}", bundle.display()))?;
            let result = ImportResult {
                identity: parsed.identity.clone(),
                fingerprint: parsed.fingerprint.clone(),
                vault: vault_path.to_string_lossy().to_string(),
                bundle_path: vault_path
                    .join("bundles")
                    .join(format!(
                        "{}.json",
                        syftbox_sdk::sanitize_identity(&parsed.identity)
                    ))
                    .to_string_lossy()
                    .to_string(),
            };
            if json {
                print_json(&result)?;
            } else {
                println!("✓ Imported bundle for {}", parsed.identity);
                println!("  Fingerprint: {}", parsed.fingerprint);
                println!("  Vault: {}", vault_path.display());
            }
        }
        KeyCommands::Export {
            email,
            output,
            vault,
            data_dir,
            json,
        } => {
            let email = resolve_email(email.as_deref(), config)?;
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;
            let bundle_path = output.unwrap_or_else(|| {
                data_root
                    .join("datasites")
                    .join(&email)
                    .join("public/crypto/did.json")
            });
            let info = load_existing_bundle(&vault_path, &email)?.ok_or_else(|| {
                crate::error::Error::Anyhow(anyhow!("no bundle found in vault for {email}"))
            })?;
            if let Some(parent) = bundle_path.parent() {
                fs::create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }
            fs::write(
                &bundle_path,
                serde_json::to_vec_pretty(&info.value)
                    .context("failed to serialise public bundle")?,
            )
            .with_context(|| format!("failed to write {}", bundle_path.display()))?;
            let result = ExportResult {
                identity: email.clone(),
                fingerprint: info.fingerprint,
                output: bundle_path.to_string_lossy().to_string(),
                vault_bundle: vault_path
                    .join("bundles")
                    .join(format!("{}.json", syftbox_sdk::sanitize_identity(&email)))
                    .to_string_lossy()
                    .to_string(),
            };
            if json {
                print_json(&result)?;
            } else {
                println!(
                    "✓ Exported bundle for {} -> {}",
                    email,
                    bundle_path.display()
                );
            }
        }
        KeyCommands::Restore {
            email,
            mnemonic,
            vault,
            data_dir,
            json,
        } => {
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;
            let outcome =
                restore_identity_from_mnemonic(&email, &mnemonic, &data_root, Some(&vault_path))
                    .with_context(|| format!("failed to restore identity {}", email))?;
            let bundle = parse_public_bundle_file(&outcome.public_bundle_path)?;
            let result = GenerateResult {
                identity: bundle.identity.clone(),
                fingerprint: bundle.fingerprint.clone(),
                vault: outcome.vault_path.to_string_lossy().to_string(),
                bundle_path: outcome.bundle_path.to_string_lossy().to_string(),
                export_path: outcome.public_bundle_path.to_string_lossy().to_string(),
                mnemonic: outcome.recovery_mnemonic.clone(),
            };
            if json {
                print_json(&result)?;
            } else {
                println!("✓ Restored identity {}", bundle.identity);
                println!("  Fingerprint: {}", bundle.fingerprint);
                println!("  Vault: {}", outcome.vault_path.display());
                println!("  Bundle: {}", outcome.bundle_path.display());
            }
        }
        KeyCommands::Wipe {
            email,
            vault,
            data_dir,
            json,
        } => {
            let email = resolve_email(email.as_deref(), config)?;
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;
            let slug = syftbox_sdk::sanitize_identity(&email);
            let key_path = vault_path.join("keys").join(format!("{slug}.key"));
            let bundle_path = vault_path.join("bundles").join(format!("{slug}.json"));
            let export_path = resolve_export_path(&data_root, &email);
            let mut removed_any = false;
            for path in [&key_path, &bundle_path, &export_path] {
                if path.exists() {
                    fs::remove_file(path)
                        .with_context(|| format!("failed to remove {}", path.display()))?;
                    removed_any = true;
                }
            }
            let result = WipeResult {
                identity: email.clone(),
                removed: removed_any,
                key_path: key_path.to_string_lossy().to_string(),
                bundle_path: bundle_path.to_string_lossy().to_string(),
                export_path: export_path.to_string_lossy().to_string(),
            };
            if json {
                print_json(&result)?;
            } else if removed_any {
                println!("✓ Wiped key material for {}", email);
            } else {
                println!("ℹ️ No key material found for {}", email);
            }
        }
        KeyCommands::Status {
            email,
            vault,
            data_dir,
            json,
        } => {
            let email = resolve_email(email.as_deref(), config)?;
            let (data_root, vault_path) =
                resolve_paths(config, data_dir.as_deref(), vault.as_deref())?;
            let bundle_path = vault_path
                .join("bundles")
                .join(format!("{}.json", syftbox_sdk::sanitize_identity(&email)));
            let export_path = resolve_export_path(&data_root, &email);
            let existing = load_existing_bundle(&vault_path, &email)?;

            let mut export_fp = None;
            let mut export_matches = None;
            if let Some(info) = existing.as_ref() {
                if export_path.exists() {
                    if let Ok(einfo) = parse_public_bundle_file(&export_path) {
                        export_matches = Some(einfo.fingerprint == info.fingerprint);
                        export_fp = Some(einfo.fingerprint);
                    }
                }
            }
            let result = StatusResult {
                identity: email.clone(),
                vault: vault_path.to_string_lossy().to_string(),
                bundle: bundle_path.to_string_lossy().to_string(),
                export: export_path.to_string_lossy().to_string(),
                vault_fingerprint: existing.as_ref().map(|i| i.fingerprint.clone()),
                export_fingerprint: export_fp,
                export_matches,
            };
            if json {
                print_json(&result)?;
            } else {
                println!("Identity: {}", email);
                println!("Vault: {}", vault_path.display());
                println!("Bundle: {}", bundle_path.display());
                println!("Export: {}", export_path.display());
                match existing {
                    Some(info) => {
                        println!("  Vault bundle fingerprint: {}", info.fingerprint);
                        if let Some(fp) = result.export_fingerprint {
                            println!(
                                "  Export fingerprint: {} (matches vault: {})",
                                fp,
                                result.export_matches.unwrap_or(false)
                            );
                        } else if export_path.exists() {
                            println!("  Export unreadable");
                        } else {
                            println!("  Export missing");
                        }
                    }
                    None => {
                        println!("  No bundle found in vault");
                    }
                }
            }
        }
        KeyCommands::VaultList { vault, json } => {
            let (vault_path, _) = resolve_vault_only(config, vault.as_deref())?;
            let keys_dir = vault_path.join("keys");
            if !keys_dir.exists() {
                let result = VaultListResult {
                    vault: vault_path.to_string_lossy().to_string(),
                    entries: Vec::new(),
                };
                if json {
                    print_json(&result)?;
                } else {
                    println!("(no vault at {})", vault_path.display());
                }
                return Ok(());
            }
            let mut entries = Vec::new();
            for entry in fs::read_dir(&keys_dir)? {
                let entry = entry?;
                if entry.file_type()?.is_file() {
                    if let Some(id) = entry
                        .file_name()
                        .to_str()
                        .and_then(|s| s.strip_suffix(".key"))
                    {
                        let slug = syftbox_sdk::sanitize_identity(id);
                        let bundle_path = vault_path.join("bundles").join(format!("{slug}.json"));
                        let fingerprint = parse_public_bundle_file(&bundle_path)
                            .ok()
                            .map(|info| info.fingerprint);
                        entries.push(VaultEntry {
                            identity: id.to_string(),
                            bundle: bundle_path.to_string_lossy().to_string(),
                            fingerprint,
                        });
                    }
                }
            }
            let result = VaultListResult {
                vault: vault_path.to_string_lossy().to_string(),
                entries,
            };
            if json {
                print_json(&result)?;
            } else if result.entries.is_empty() {
                println!("(vault empty at {})", vault_path.display());
            } else {
                println!("Vault at {}", vault_path.display());
                for e in &result.entries {
                    if let Some(fp) = &e.fingerprint {
                        println!("  - {} ({})", e.identity, fp);
                    } else {
                        println!("  - {}", e.identity);
                    }
                }
            }
        }
    }
    Ok(())
}

fn resolve_email<'a>(email: Option<&'a str>, config: &'a Config) -> Result<String> {
    if let Some(e) = email {
        return Ok(e.to_string());
    }
    if !config.email.trim().is_empty() {
        return Ok(config.email.clone());
    }
    Err(anyhow!("email is required").into())
}

fn resolve_paths(
    config: &Config,
    data_override: Option<&Path>,
    vault_override: Option<&Path>,
) -> Result<(PathBuf, PathBuf)> {
    let data_root = if let Some(dir) = data_override {
        dir.to_path_buf()
    } else {
        config.get_syftbox_data_dir()?
    };
    let encrypted_root = syftbox_sdk::syftbox::syc::resolve_encrypted_root(&data_root);
    let vault_path = resolve_vault_default(vault_override)?;
    Ok((encrypted_root, vault_path))
}

fn resolve_vault_only(
    config: &Config,
    vault_override: Option<&Path>,
) -> Result<(PathBuf, PathBuf)> {
    let data_root = config.get_syftbox_data_dir()?;
    let encrypted_root = syftbox_sdk::syftbox::syc::resolve_encrypted_root(&data_root);
    let vault_path = resolve_vault_default(vault_override)?;
    Ok((vault_path, encrypted_root))
}

fn load_existing_bundle(vault_path: &Path, identity: &str) -> Result<Option<PublicBundleInfo>> {
    let slug = syftbox_sdk::sanitize_identity(identity);
    let bundle_path = vault_path.join("bundles").join(format!("{slug}.json"));
    if !bundle_path.exists() {
        return Ok(None);
    }
    let info = parse_public_bundle_file(&bundle_path)?;
    Ok(Some(info))
}

fn resolve_vault_default(vault_override: Option<&Path>) -> Result<PathBuf> {
    if let Some(v) = vault_override {
        return Ok(v.to_path_buf());
    }
    if let Some(env_vault) = std::env::var_os("SYC_VAULT") {
        return Ok(PathBuf::from(env_vault));
    }
    Ok(crate::config::resolve_syc_vault_path()?)
}

/// Resolve the vault path using the same logic as key commands.
/// This ensures consistent vault location between `bv key` commands and message sync.
pub fn resolve_vault_for_config(_config: &Config) -> Result<PathBuf> {
    resolve_vault_default(None)
}

fn resolve_export_path(data_root: &Path, identity: &str) -> PathBuf {
    let base = if data_root
        .file_name()
        .map(|n| n == "datasites")
        .unwrap_or(false)
    {
        data_root.to_path_buf()
    } else {
        data_root.join("datasites")
    };
    base.join(identity)
        .join("public")
        .join("crypto")
        .join("did.json")
}

fn print_json<T: Serialize>(value: &T) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(value)?);
    Ok(())
}

#[derive(Serialize)]
struct GenerateResult {
    identity: String,
    fingerprint: String,
    vault: String,
    bundle_path: String,
    export_path: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    mnemonic: Option<String>,
}

#[derive(Serialize)]
struct ImportResult {
    identity: String,
    fingerprint: String,
    vault: String,
    bundle_path: String,
}

#[derive(Serialize)]
struct ExportResult {
    identity: String,
    fingerprint: String,
    output: String,
    vault_bundle: String,
}

#[derive(Serialize)]
struct WipeResult {
    identity: String,
    removed: bool,
    key_path: String,
    bundle_path: String,
    export_path: String,
}

#[derive(Serialize, Clone)]
struct StatusResult {
    identity: String,
    vault: String,
    bundle: String,
    export: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    vault_fingerprint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    export_fingerprint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    export_matches: Option<bool>,
}

#[derive(Serialize)]
struct VaultListResult {
    vault: String,
    entries: Vec<VaultEntry>,
}

#[derive(Serialize)]
struct VaultEntry {
    identity: String,
    bundle: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    fingerprint: Option<String>,
}
