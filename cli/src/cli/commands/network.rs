//! CLI commands for network discovery and status
//!
//! Network commands provide visibility into the SyftBox network,
//! showing discovered datasites, published datasets, and connectivity status.

use crate::config::Config;
use anyhow::{Context, Result};
use clap::Subcommand;
use colored::Colorize;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum NetworkCommands {
    #[command(about = "Scan network for published datasets")]
    Datasets {
        #[arg(long, help = "Output as JSON")]
        json: bool,

        #[arg(long, help = "Only show datasets from trusted contacts")]
        trusted_only: bool,

        #[arg(long, help = "Only show your own datasets")]
        own_only: bool,
    },

    #[command(about = "Show network status and connectivity")]
    Status {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Check if a specific dataset is publicly visible")]
    CheckDataset {
        #[arg(help = "Dataset name to check")]
        name: String,

        #[arg(long, help = "Output as JSON")]
        json: bool,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscoveredDataset {
    pub name: String,
    pub owner: String,
    pub description: Option<String>,
    pub version: Option<String>,
    pub schema: Option<String>,
    pub public_url: Option<String>,
    pub asset_count: usize,
    pub has_mock_data: bool,
    pub is_trusted: bool,
    pub is_own: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkStatus {
    pub current_identity: String,
    pub datasites_count: usize,
    pub trusted_contacts_count: usize,
    pub own_datasets_count: usize,
    pub network_datasets_count: usize,
    pub syftbox_data_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetVisibilityCheck {
    pub name: String,
    pub is_published: bool,
    pub has_manifest: bool,
    pub in_index: bool,
    pub has_mock_data: bool,
    pub public_url: Option<String>,
    pub issues: Vec<String>,
    pub suggestions: Vec<String>,
}

/// Handle network subcommands
pub async fn handle(command: NetworkCommands, config: &Config) -> Result<()> {
    match command {
        NetworkCommands::Datasets {
            json,
            trusted_only,
            own_only,
        } => scan_datasets(config, json, trusted_only, own_only),
        NetworkCommands::Status { json } => show_status(config, json),
        NetworkCommands::CheckDataset { name, json } => check_dataset(config, &name, json),
    }
}

/// Scan network for published datasets
fn scan_datasets(
    config: &Config,
    json_output: bool,
    trusted_only: bool,
    own_only: bool,
) -> Result<()> {
    let current_email = config.email.clone();
    let data_dir = config.get_syftbox_data_dir()?;
    let datasites_dir = data_dir.join("datasites");

    // Resolve vault path for trusted contacts
    let vault_path = std::env::var_os("SYC_VAULT")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            dirs::home_dir()
                .map(|h| h.join(".syc"))
                .unwrap_or_else(|| PathBuf::from(".syc"))
        });
    let bundles_dir = vault_path.join("bundles");

    let mut datasets = Vec::new();
    let current_slug = syftbox_sdk::sanitize_identity(&current_email);

    if !datasites_dir.exists() {
        if json_output {
            println!("{}", serde_json::to_string_pretty(&datasets)?);
        } else {
            println!("No datasites directory found. SyftBox may not be initialized.");
            println!("\nRun {} to initialize.", "bv init".cyan());
        }
        return Ok(());
    }

    let entries =
        std::fs::read_dir(&datasites_dir).context("Failed to read datasites directory")?;

    for entry in entries.flatten() {
        let datasite_path = entry.path();
        if !datasite_path.is_dir() {
            continue;
        }

        let owner = datasite_path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();

        let owner_slug = syftbox_sdk::sanitize_identity(&owner);
        let is_own = owner_slug == current_slug;
        let is_trusted = is_own || bundles_dir.join(format!("{}.json", owner_slug)).exists();

        // Apply filters
        if own_only && !is_own {
            continue;
        }
        if trusted_only && !is_trusted {
            continue;
        }

        // Look for datasets.yaml index
        let index_path = datasite_path
            .join("public")
            .join("biovault")
            .join("datasets.yaml");

        if !index_path.exists() {
            continue;
        }

        // Parse the datasets index
        let index_bytes = match std::fs::read(&index_path) {
            Ok(b) => b,
            Err(_) => continue,
        };

        let index: super::datasets::DatasetIndex = match serde_yaml::from_slice(&index_bytes) {
            Ok(i) => i,
            Err(_) => continue,
        };

        // Load each dataset's manifest
        for resource in index.resources {
            let dataset_dir = datasite_path
                .join("public")
                .join("biovault")
                .join("datasets")
                .join(&resource.name);
            let manifest_path = dataset_dir.join("dataset.yaml");

            if !manifest_path.exists() {
                continue;
            }

            let manifest_bytes = match std::fs::read(&manifest_path) {
                Ok(b) => b,
                Err(_) => continue,
            };

            let manifest: super::datasets::DatasetManifest =
                match serde_yaml::from_slice(&manifest_bytes) {
                    Ok(m) => m,
                    Err(_) => continue,
                };

            // Check for mock data
            let has_mock_data = manifest.assets.values().any(|asset| asset.mock.is_some());

            datasets.push(DiscoveredDataset {
                name: manifest.name.clone(),
                owner: owner.clone(),
                description: manifest.description.clone(),
                version: manifest.version.clone(),
                schema: manifest.schema.clone(),
                public_url: manifest.public_url.clone(),
                asset_count: manifest.assets.len(),
                has_mock_data,
                is_trusted,
                is_own,
            });
        }
    }

    // Sort: own first, then trusted, then by name
    datasets.sort_by(|a, b| match (a.is_own, b.is_own) {
        (true, false) => std::cmp::Ordering::Less,
        (false, true) => std::cmp::Ordering::Greater,
        _ => match (a.is_trusted, b.is_trusted) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => a.name.to_lowercase().cmp(&b.name.to_lowercase()),
        },
    });

    if json_output {
        println!("{}", serde_json::to_string_pretty(&datasets)?);
        return Ok(());
    }

    if datasets.is_empty() {
        println!("No datasets found on the network.");
        println!("\n{}", "To publish a dataset:".bold());
        println!(
            "  1. Create a dataset in the Data tab or with: {}",
            "bv dataset create".cyan()
        );
        println!(
            "  2. Publish it with: {}",
            "bv dataset publish --name <name>".cyan()
        );
        println!("\n{}", "For others to see your datasets:".bold());
        println!("  • Your datasite must be synced via SyftBox");
        println!("  • They must have your datasite in their datasites folder");
        return Ok(());
    }

    println!("\n📊 {} Datasets on Network", "Discovered".bold());
    println!("═══════════════════════════════════════════════════════════════\n");

    let own_datasets: Vec<_> = datasets.iter().filter(|d| d.is_own).collect();
    let trusted_datasets: Vec<_> = datasets
        .iter()
        .filter(|d| !d.is_own && d.is_trusted)
        .collect();
    let other_datasets: Vec<_> = datasets
        .iter()
        .filter(|d| !d.is_own && !d.is_trusted)
        .collect();

    if !own_datasets.is_empty() {
        println!("  {} {}", "📁".bold(), "Your Datasets".green().bold());
        for dataset in &own_datasets {
            print_dataset(dataset);
        }
        println!();
    }

    if !trusted_datasets.is_empty() {
        println!(
            "  {} {}",
            "🤝".bold(),
            "From Trusted Contacts".cyan().bold()
        );
        for dataset in &trusted_datasets {
            print_dataset(dataset);
        }
        println!();
    }

    if !other_datasets.is_empty() {
        println!("  {} {}", "🌐".bold(), "From Other Peers".dimmed().bold());
        for dataset in &other_datasets {
            print_dataset(dataset);
        }
        println!();
    }

    println!("─────────────────────────────────────────────────────────────────");
    println!(
        "Total: {} dataset(s) | {} yours | {} from trusted | {} from others",
        datasets.len(),
        own_datasets.len(),
        trusted_datasets.len(),
        other_datasets.len()
    );

    if !other_datasets.is_empty() {
        println!(
            "\n💡 {} to encrypt messages to untrusted peers.",
            "Import their contact first".yellow()
        );
        println!("   Use: {}", "bv contacts import <email>".cyan());
    }

    Ok(())
}

fn print_dataset(dataset: &DiscoveredDataset) {
    let mock_indicator = if dataset.has_mock_data {
        "📋".to_string()
    } else {
        "  ".to_string()
    };

    let trust_indicator = if dataset.is_own {
        "✓".green().to_string()
    } else if dataset.is_trusted {
        "🔐".to_string()
    } else {
        "?".dimmed().to_string()
    };

    println!(
        "     {} {} {} ({})",
        trust_indicator,
        dataset.name.bold(),
        mock_indicator,
        dataset.owner.dimmed()
    );

    if let Some(desc) = &dataset.description {
        if !desc.is_empty() {
            let short_desc = if desc.len() > 60 {
                format!("{}...", &desc[..57])
            } else {
                desc.clone()
            };
            println!("        {}", short_desc.dimmed());
        }
    }

    println!(
        "        {} asset(s) | v{}",
        dataset.asset_count,
        dataset.version.as_deref().unwrap_or("1.0.0")
    );
}

/// Show network status
fn show_status(config: &Config, json_output: bool) -> Result<()> {
    let current_email = config.email.clone();
    let data_dir = config.get_syftbox_data_dir()?;
    let datasites_dir = data_dir.join("datasites");

    let vault_path = std::env::var_os("SYC_VAULT")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            dirs::home_dir()
                .map(|h| h.join(".syc"))
                .unwrap_or_else(|| PathBuf::from(".syc"))
        });
    let bundles_dir = vault_path.join("bundles");

    let mut datasites_count = 0;
    let mut trusted_count = 0;
    let mut own_datasets = 0;
    let mut network_datasets = 0;

    // Count trusted contacts
    if bundles_dir.exists() {
        for entry in std::fs::read_dir(&bundles_dir)
            .into_iter()
            .flatten()
            .flatten()
        {
            if entry
                .path()
                .extension()
                .map(|e| e == "json")
                .unwrap_or(false)
            {
                trusted_count += 1;
            }
        }
    }

    // Count datasites and datasets
    if datasites_dir.exists() {
        let current_slug = syftbox_sdk::sanitize_identity(&current_email);

        for entry in std::fs::read_dir(&datasites_dir)
            .into_iter()
            .flatten()
            .flatten()
        {
            let datasite_path = entry.path();
            if !datasite_path.is_dir() {
                continue;
            }
            datasites_count += 1;

            let owner = datasite_path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
                .to_string();
            let owner_slug = syftbox_sdk::sanitize_identity(&owner);
            let is_own = owner_slug == current_slug;

            let index_path = datasite_path
                .join("public")
                .join("biovault")
                .join("datasets.yaml");

            if index_path.exists() {
                if let Ok(bytes) = std::fs::read(&index_path) {
                    if let Ok(index) =
                        serde_yaml::from_slice::<super::datasets::DatasetIndex>(&bytes)
                    {
                        if is_own {
                            own_datasets += index.resources.len();
                        } else {
                            network_datasets += index.resources.len();
                        }
                    }
                }
            }
        }
    }

    let status = NetworkStatus {
        current_identity: current_email.clone(),
        datasites_count,
        trusted_contacts_count: trusted_count,
        own_datasets_count: own_datasets,
        network_datasets_count: network_datasets,
        syftbox_data_dir: data_dir.to_string_lossy().to_string(),
    };

    if json_output {
        println!("{}", serde_json::to_string_pretty(&status)?);
        return Ok(());
    }

    println!("\n🌐 {} Status", "Network".bold());
    println!("═══════════════════════════════════════════════════════════════\n");

    println!("  {} {}", "Identity:".bold(), current_email.green());
    println!("  {} {}", "Data Dir:".bold(), data_dir.display());
    println!();

    println!("  {} {}", "📁 Datasites synced:".bold(), datasites_count);
    println!("  {} {}", "🤝 Trusted contacts:".bold(), trusted_count);
    println!("  {} {}", "📊 Your datasets:".bold(), own_datasets);
    println!("  {} {}", "🌐 Network datasets:".bold(), network_datasets);

    println!("\n─────────────────────────────────────────────────────────────────");

    if trusted_count == 0 {
        println!(
            "\n💡 {} to collaborate securely.",
            "Add trusted contacts".yellow()
        );
        println!("   Scan for peers: {}", "bv contacts scan".cyan());
        println!("   Import contact: {}", "bv contacts import <email>".cyan());
    }

    if own_datasets == 0 {
        println!(
            "\n💡 {} to share with others.",
            "Publish a dataset".yellow()
        );
        println!("   Create dataset: {}", "bv dataset create".cyan());
        println!(
            "   Publish it:     {}",
            "bv dataset publish --name <name>".cyan()
        );
    }

    Ok(())
}

/// Check if a specific dataset is publicly visible
fn check_dataset(config: &Config, name: &str, json_output: bool) -> Result<()> {
    let email = &config.email;
    let data_dir = config.get_syftbox_data_dir()?;

    let public_dir = data_dir
        .join("datasites")
        .join(email)
        .join("public")
        .join("biovault");

    let index_path = public_dir.join("datasets.yaml");
    let dataset_dir = public_dir.join("datasets").join(name);
    let manifest_path = dataset_dir.join("dataset.yaml");

    let mut check = DatasetVisibilityCheck {
        name: name.to_string(),
        is_published: false,
        has_manifest: false,
        in_index: false,
        has_mock_data: false,
        public_url: None,
        issues: Vec::new(),
        suggestions: Vec::new(),
    };

    // Check if manifest exists
    if manifest_path.exists() {
        check.has_manifest = true;

        if let Ok(bytes) = std::fs::read(&manifest_path) {
            if let Ok(manifest) = serde_yaml::from_slice::<super::datasets::DatasetManifest>(&bytes)
            {
                check.public_url = manifest.public_url.clone();
                check.has_mock_data = manifest.assets.values().any(|a| a.mock.is_some());

                if !check.has_mock_data {
                    check
                        .issues
                        .push("No mock data configured for assets".to_string());
                    check.suggestions.push(
                        "Add mock data to let others preview your dataset without accessing real data"
                            .to_string(),
                    );
                }
            }
        }
    } else {
        check
            .issues
            .push(format!("Manifest not found at {}", manifest_path.display()));
        check
            .suggestions
            .push(format!("Run: bv dataset publish --name {}", name));
    }

    // Check if in index
    if index_path.exists() {
        if let Ok(bytes) = std::fs::read(&index_path) {
            if let Ok(index) = serde_yaml::from_slice::<super::datasets::DatasetIndex>(&bytes) {
                check.in_index = index.resources.iter().any(|r| r.name == name);
            }
        }

        if !check.in_index {
            check
                .issues
                .push("Dataset not listed in datasets.yaml index".to_string());
            check.suggestions.push(format!(
                "Run: bv dataset publish --name {} (this updates the index)",
                name
            ));
        }
    } else {
        check
            .issues
            .push("datasets.yaml index file not found".to_string());
        check.suggestions.push(format!(
            "Run: bv dataset publish --name {} (this creates the index)",
            name
        ));
    }

    check.is_published = check.has_manifest && check.in_index;

    if json_output {
        println!("{}", serde_json::to_string_pretty(&check)?);
        return Ok(());
    }

    println!("\n🔍 {} Visibility Check", "Dataset".bold());
    println!("═══════════════════════════════════════════════════════════════\n");

    println!("  {} {}", "Dataset:".bold(), name.cyan());

    if check.is_published {
        println!(
            "  {} {}",
            "Status:".bold(),
            "✅ Published & Visible".green()
        );
    } else {
        println!("  {} {}", "Status:".bold(), "❌ Not Publicly Visible".red());
    }

    println!();
    println!(
        "  {} {}",
        "Has manifest:".bold(),
        bool_icon(check.has_manifest)
    );
    println!("  {} {}", "In index:    ".bold(), bool_icon(check.in_index));
    println!(
        "  {} {}",
        "Has mock data:".bold(),
        bool_icon(check.has_mock_data)
    );

    if let Some(url) = &check.public_url {
        println!("  {} {}", "Public URL:  ".bold(), url.dimmed());
    }

    if !check.issues.is_empty() {
        println!("\n  {} {}", "⚠️ ".bold(), "Issues:".yellow().bold());
        for issue in &check.issues {
            println!("     • {}", issue);
        }
    }

    if !check.suggestions.is_empty() {
        println!("\n  {} {}", "💡".bold(), "Suggestions:".cyan().bold());
        for suggestion in &check.suggestions {
            println!("     • {}", suggestion);
        }
    }

    if check.is_published {
        println!(
            "\n  {} Others can now discover this dataset via the Network tab.",
            "✨".bold()
        );
    }

    Ok(())
}

fn bool_icon(value: bool) -> String {
    if value {
        "✓".green().to_string()
    } else {
        "✗".red().to_string()
    }
}
