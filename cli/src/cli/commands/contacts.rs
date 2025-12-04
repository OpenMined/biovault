//! CLI commands for managing contacts (trusted public key bundles)
//!
//! Contacts are peers whose public keys you've imported, allowing you to send
//! encrypted messages and collaborate with them.

use crate::config::Config;
use anyhow::{Context, Result};
use clap::Subcommand;
use colored::Colorize;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum ContactsCommands {
    #[command(about = "List your trusted contacts")]
    List {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Scan network for discoverable peers")]
    Scan {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Import a contact's public key from their datasite")]
    Import {
        #[arg(help = "Email/identity of the contact to import")]
        identity: String,
    },

    #[command(about = "Remove a contact from your trusted list")]
    Remove {
        #[arg(help = "Email/identity of the contact to remove")]
        identity: String,

        #[arg(long, help = "Skip confirmation prompt")]
        yes: bool,
    },

    #[command(about = "Trust a contact's changed key (re-import)")]
    Trust {
        #[arg(help = "Email/identity of the contact to trust")]
        identity: String,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactInfo {
    pub identity: String,
    pub fingerprint: String,
    pub bundle_path: String,
    #[serde(default)]
    pub has_changed: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local_fingerprint: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscoveredPeer {
    pub identity: String,
    pub fingerprint: String,
    pub did_path: String,
    pub is_imported: bool,
    #[serde(default)]
    pub has_changed: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScanResult {
    pub contacts: Vec<DiscoveredPeer>,
    pub discovered: Vec<DiscoveredPeer>,
    pub current_identity: String,
}

/// Handle contacts subcommands
pub async fn handle(command: ContactsCommands, config: &Config) -> Result<()> {
    match command {
        ContactsCommands::List { json } => list_contacts(config, json),
        ContactsCommands::Scan { json } => scan_network(config, json),
        ContactsCommands::Import { identity } => import_contact(config, &identity),
        ContactsCommands::Remove { identity, yes } => remove_contact(config, &identity, yes),
        ContactsCommands::Trust { identity } => trust_changed_key(config, &identity),
    }
}

/// Get paths for vault and datasites
fn resolve_paths(config: &Config) -> Result<(PathBuf, PathBuf)> {
    let data_dir = config.get_syftbox_data_dir()?;

    // Vault is in .biovault under data dir
    let vault_path = data_dir.join(".biovault").join("vault");

    // Datasites dir
    let datasites_dir = if data_dir
        .file_name()
        .map(|n| n == "datasites")
        .unwrap_or(false)
    {
        data_dir.clone()
    } else {
        data_dir.join("datasites")
    };

    Ok((datasites_dir, vault_path))
}

/// List all trusted contacts
fn list_contacts(config: &Config, json_output: bool) -> Result<()> {
    let (_, vault_path) = resolve_paths(config)?;
    let bundles_dir = vault_path.join("bundles");

    let mut contacts = Vec::new();

    if bundles_dir.exists() {
        for entry in std::fs::read_dir(&bundles_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.extension().map(|e| e == "json").unwrap_or(false) {
                if let Ok(info) = crate::syftbox::syc::parse_public_bundle_file(&path) {
                    contacts.push(ContactInfo {
                        identity: info.identity,
                        fingerprint: info.fingerprint,
                        bundle_path: path.to_string_lossy().to_string(),
                        has_changed: false,
                        local_fingerprint: None,
                    });
                }
            }
        }
    }

    contacts.sort_by(|a, b| a.identity.to_lowercase().cmp(&b.identity.to_lowercase()));

    if json_output {
        println!("{}", serde_json::to_string_pretty(&contacts)?);
        return Ok(());
    }

    if contacts.is_empty() {
        println!("No trusted contacts yet.");
        println!(
            "\nUse {} to discover peers on the network.",
            "bv contacts scan".cyan()
        );
        println!(
            "Use {} to add a contact.",
            "bv contacts import <email>".cyan()
        );
        return Ok(());
    }

    println!("\n👥 {} Trusted Contacts", "Your".bold());
    println!("═══════════════════════════════════════════════════════════════\n");

    for contact in &contacts {
        let short_fp = &contact.fingerprint[..16.min(contact.fingerprint.len())];
        println!("  📧 {}", contact.identity.green());
        println!("     🔑 {}...", short_fp.dimmed());
        println!();
    }

    println!("Total: {} contact(s)", contacts.len());

    Ok(())
}

/// Scan network for discoverable peers
fn scan_network(config: &Config, json_output: bool) -> Result<()> {
    let current_email = config.email.clone();
    let (datasites_dir, vault_path) = resolve_paths(config)?;
    let bundles_dir = vault_path.join("bundles");

    let mut contacts = Vec::new();
    let mut discovered = Vec::new();

    let current_slug = syftbox_sdk::sanitize_identity(&current_email);

    if datasites_dir.exists() {
        for entry in std::fs::read_dir(&datasites_dir)? {
            let entry = entry?;
            let datasite_path = entry.path();

            if !datasite_path.is_dir() {
                continue;
            }

            let did_path = datasite_path.join("public").join("crypto").join("did.json");

            if !did_path.exists() {
                continue;
            }

            if let Ok(remote_info) = crate::syftbox::syc::parse_public_bundle_file(&did_path) {
                let slug = syftbox_sdk::sanitize_identity(&remote_info.identity);

                // Skip current identity
                if slug == current_slug {
                    continue;
                }

                let local_bundle_path = bundles_dir.join(format!("{slug}.json"));
                let is_imported = local_bundle_path.exists();

                let has_changed = if is_imported {
                    match crate::syftbox::syc::parse_public_bundle_file(&local_bundle_path) {
                        Ok(local_info) => local_info.fingerprint != remote_info.fingerprint,
                        Err(_) => false,
                    }
                } else {
                    false
                };

                let peer = DiscoveredPeer {
                    identity: remote_info.identity,
                    fingerprint: remote_info.fingerprint,
                    did_path: did_path.to_string_lossy().to_string(),
                    is_imported,
                    has_changed,
                };

                if is_imported {
                    contacts.push(peer);
                } else {
                    discovered.push(peer);
                }
            }
        }
    }

    contacts.sort_by(|a, b| a.identity.to_lowercase().cmp(&b.identity.to_lowercase()));
    discovered.sort_by(|a, b| a.identity.to_lowercase().cmp(&b.identity.to_lowercase()));

    let result = ScanResult {
        contacts,
        discovered,
        current_identity: current_email.clone(),
    };

    if json_output {
        println!("{}", serde_json::to_string_pretty(&result)?);
        return Ok(());
    }

    println!(
        "\n🔍 {} {}",
        "Network Scan for".bold(),
        current_email.cyan()
    );
    println!("═══════════════════════════════════════════════════════════════\n");

    // Show contacts with key changes
    let changed: Vec<_> = result.contacts.iter().filter(|c| c.has_changed).collect();
    if !changed.is_empty() {
        println!("⚠️  {} Key(s) Changed:", "Warning:".yellow().bold());
        for peer in changed {
            println!("   {} {}", "⚡".yellow(), peer.identity.yellow());
            println!(
                "      Run: {} to trust the new key",
                format!("bv contacts trust {}", peer.identity).cyan()
            );
        }
        println!();
    }

    // Show trusted contacts
    println!(
        "✅ {} Trusted Contact(s):",
        result.contacts.len().to_string().green()
    );
    if result.contacts.is_empty() {
        println!("   (none)");
    } else {
        for peer in &result.contacts {
            let status = if peer.has_changed {
                "⚠️ key changed".yellow().to_string()
            } else {
                "✓".green().to_string()
            };
            println!("   {} {} {}", status, peer.identity, "".dimmed());
        }
    }
    println!();

    // Show discoverable peers
    println!(
        "🌐 {} Discoverable Peer(s):",
        result.discovered.len().to_string().blue()
    );
    if result.discovered.is_empty() {
        println!("   (none found)");
    } else {
        for peer in &result.discovered {
            let short_fp = &peer.fingerprint[..16.min(peer.fingerprint.len())];
            println!("   📧 {}", peer.identity);
            println!("      🔑 {}...", short_fp.dimmed());
        }
    }
    println!();

    if !result.discovered.is_empty() {
        println!(
            "💡 To add a contact: {}",
            "bv contacts import <email>".cyan()
        );
    }

    Ok(())
}

/// Import a contact's public key bundle
fn import_contact(config: &Config, identity: &str) -> Result<()> {
    let (datasites_dir, vault_path) = resolve_paths(config)?;
    let bundles_dir = vault_path.join("bundles");

    // Ensure bundles directory exists
    if !bundles_dir.exists() {
        std::fs::create_dir_all(&bundles_dir).context("Failed to create bundles directory")?;
    }

    let did_path = datasites_dir
        .join(identity)
        .join("public")
        .join("crypto")
        .join("did.json");

    if !did_path.exists() {
        anyhow::bail!(
            "DID not found for {}. Make sure the peer's datasite is synced.",
            identity
        );
    }

    let remote_info = crate::syftbox::syc::parse_public_bundle_file(&did_path)
        .context("Failed to parse DID file")?;

    let slug = syftbox_sdk::sanitize_identity(&remote_info.identity);
    let local_bundle_path = bundles_dir.join(format!("{slug}.json"));

    let was_existing = local_bundle_path.exists();

    std::fs::copy(&did_path, &local_bundle_path).context("Failed to copy public key bundle")?;

    let action = if was_existing { "Updated" } else { "Imported" };
    let short_fp = &remote_info.fingerprint[..16.min(remote_info.fingerprint.len())];

    println!("\n✅ {} contact: {}", action, remote_info.identity.green());
    println!("   🔑 Fingerprint: {}...", short_fp);
    println!("   📁 Saved to: {}", local_bundle_path.display());
    println!();
    println!(
        "You can now send encrypted messages with: {}",
        format!("bv message send {} \"Hello!\"", identity).cyan()
    );

    Ok(())
}

/// Remove a contact from trusted list
fn remove_contact(config: &Config, identity: &str, skip_confirm: bool) -> Result<()> {
    let (_, vault_path) = resolve_paths(config)?;
    let bundles_dir = vault_path.join("bundles");

    let slug = syftbox_sdk::sanitize_identity(identity);
    let bundle_path = bundles_dir.join(format!("{slug}.json"));

    if !bundle_path.exists() {
        anyhow::bail!("Contact '{}' not found in your trusted list.", identity);
    }

    if !skip_confirm {
        use dialoguer::Confirm;
        let confirmed = Confirm::new()
            .with_prompt(format!(
                "Remove {} from your contacts?\nYou won't be able to send encrypted messages until you re-add them.",
                identity
            ))
            .default(false)
            .interact()
            .unwrap_or(false);

        if !confirmed {
            println!("Cancelled.");
            return Ok(());
        }
    }

    std::fs::remove_file(&bundle_path).context("Failed to remove contact bundle")?;

    println!("\n🗑️  Removed contact: {}", identity.yellow());
    println!(
        "   To re-add: {}",
        format!("bv contacts import {}", identity).cyan()
    );

    Ok(())
}

/// Trust a changed key by re-importing
fn trust_changed_key(config: &Config, identity: &str) -> Result<()> {
    println!("🔄 Re-importing key for {}...", identity);
    import_contact(config, identity)?;
    println!("\n✅ New key trusted for {}", identity.green());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn create_test_config(temp: &TempDir) -> Config {
        crate::config::set_test_syftbox_data_dir(temp.path());
        crate::config::set_test_biovault_home(temp.path().join(".biovault"));

        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
            version: None,
            binary_paths: None,
            syftbox_credentials: None,
        }
    }

    #[test]
    fn test_list_contacts_empty() {
        let temp = TempDir::new().unwrap();
        let config = create_test_config(&temp);

        // Should not error on empty bundles
        let result = list_contacts(&config, false);
        assert!(result.is_ok());

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_scan_network_empty() {
        let temp = TempDir::new().unwrap();
        let config = create_test_config(&temp);

        // Should not error on empty datasites
        let result = scan_network(&config, false);
        assert!(result.is_ok());

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_import_contact_not_found() {
        let temp = TempDir::new().unwrap();
        let config = create_test_config(&temp);

        // Should fail gracefully for non-existent contact
        let result = import_contact(&config, "nonexistent@example.com");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("DID not found"));

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_remove_contact_not_found() {
        let temp = TempDir::new().unwrap();
        let config = create_test_config(&temp);

        // Should fail gracefully for non-existent contact
        let result = remove_contact(&config, "nonexistent@example.com", true);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));

        crate::config::clear_test_syftbox_data_dir();
        crate::config::clear_test_biovault_home();
    }
}
