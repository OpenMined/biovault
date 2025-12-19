use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use uuid::Uuid;

use crate::config::Config;
use crate::data::{
    datasets::{build_manifest_from_db, get_dataset_with_assets, list_datasets, upsert_dataset},
    BioVaultDb,
};
use crate::syftbox::storage::SyftBoxStorage;

const DEFAULT_SCHEMA: &str = "net.biovault.datasets:1.0.0";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetIndex {
    pub email: String,
    pub resources: Vec<DatasetResource>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DatasetResource {
    pub name: String,
    pub path: String,
    pub schema: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DatasetAsset {
    #[serde(default)]
    pub id: Option<String>,
    #[serde(default, rename = "type")]
    pub kind: Option<String>,
    #[serde(default)]
    pub url: Option<String>,
    /// For twin_list: { url, type }. For single: string reference like "{url}.private"
    #[serde(default)]
    pub private: Option<YamlValue>,
    #[serde(default)]
    pub mock: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mappings: Option<DatasetAssetMapping>,
    #[serde(flatten, default)]
    pub extra: BTreeMap<String, YamlValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DatasetAssetMappingEndpoint {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_file_id: Option<i64>,
    /// For twin_list: array of entry objects with id, file_path, participant_id, etc.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entries: Option<Vec<YamlValue>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DatasetAssetMapping {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub private: Option<DatasetAssetMappingEndpoint>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mock: Option<DatasetAssetMappingEndpoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DatasetManifest {
    pub name: String,
    #[serde(default)]
    pub description: Option<String>,
    #[serde(default)]
    pub author: Option<String>,
    #[serde(default)]
    pub schema: Option<String>,
    #[serde(default)]
    pub version: Option<String>,
    #[serde(default)]
    pub http_relay_servers: Vec<String>,
    #[serde(default)]
    pub public_url: Option<String>,
    #[serde(default)]
    pub private_url: Option<String>,
    #[serde(default)]
    pub assets: BTreeMap<String, DatasetAsset>,
    #[serde(flatten, default)]
    pub extra: BTreeMap<String, YamlValue>,
}

pub async fn publish(
    manifest_path: Option<String>,
    dataset_name: Option<String>,
    copy_mock: bool,
) -> Result<()> {
    let config = Config::load()?;
    let email = &config.email;
    let manifest_dir: PathBuf;

    let mut manifest: DatasetManifest = if let Some(path) = manifest_path {
        let manifest_path = PathBuf::from(path);
        manifest_dir = manifest_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));

        let raw = fs::read_to_string(&manifest_path)
            .with_context(|| format!("Failed to read manifest at {}", manifest_path.display()))?;
        serde_yaml::from_str(&raw)
            .with_context(|| format!("Failed to parse manifest {}", manifest_path.display()))?
    } else if let Some(name) = dataset_name {
        let db = BioVaultDb::new()?;
        let Some((dataset, assets)) = get_dataset_with_assets(&db, &name)? else {
            return Err(anyhow!("Dataset '{}' not found in database", name));
        };
        manifest_dir = PathBuf::from(".");
        build_manifest_from_db(&dataset, &assets)
    } else {
        return Err(anyhow!(
            "Provide either a manifest path or --name to publish from the database"
        ));
    };

    if manifest.name.trim().is_empty() {
        return Err(anyhow!("Dataset manifest is missing a name"));
    }

    // Defaults and normalization
    if manifest.author.is_none() {
        manifest.author = Some(email.clone());
    }
    if manifest.schema.is_none() {
        manifest.schema = Some(DEFAULT_SCHEMA.to_string());
    }
    if manifest.version.is_none() {
        manifest.version = Some("1.0.0".to_string());
    }
    if manifest.http_relay_servers.is_empty() {
        manifest.http_relay_servers = vec!["syftbox.net".to_string()];
    }

    let dataset_name = manifest.name.clone();
    manifest.public_url.get_or_insert_with(|| {
        format!(
            "syft://{}/public/biovault/datasets/{}/dataset.yaml",
            email, dataset_name
        )
    });
    manifest.private_url.get_or_insert_with(|| {
        format!(
            "syft://{}/private/biovault/datasets/{}/dataset.yaml",
            email, dataset_name
        )
    });

    for (asset_key, asset) in manifest.assets.iter_mut() {
        if asset.id.is_none() {
            asset.id = Some(Uuid::new_v4().to_string());
        }
        if asset.kind.is_none() {
            asset.kind = Some("file".to_string());
        }
        if asset.url.is_none() {
            asset.url = Some(format!("{{root.private_url}}#assets.{}", asset_key));
        }
    }

    // Persist into DB and rebuild manifest from DB state
    let mut db = BioVaultDb::new()?;
    upsert_dataset(&mut db, &manifest)?;

    let (dataset_row, asset_rows) = get_dataset_with_assets(&db, &manifest.name)?
        .ok_or_else(|| anyhow!("Dataset '{}' not found after upsert", manifest.name))?;
    manifest = build_manifest_from_db(&dataset_row, &asset_rows);
    let dataset_name = manifest.name.clone();

    let data_dir = config.get_syftbox_data_dir()?;
    let storage = SyftBoxStorage::new(&data_dir);
    let public_dir = data_dir
        .join("datasites")
        .join(email)
        .join("public")
        .join("biovault")
        .join("datasets")
        .join(&dataset_name);

    storage
        .ensure_dir(&public_dir)
        .with_context(|| format!("Failed to create {}", public_dir.display()))?;

    // Persist dataset manifest to public and private locations
    // Strip internal-only hints before publishing YAML
    let mut manifest_for_write = manifest.clone();
    for asset in manifest_for_write.assets.values_mut() {
        asset.extra.remove("mock_source_path");
        asset.mappings = None;

        // For twin_list assets, strip internal fields from mock entries
        if let Some(mock_val) = asset.mock.as_mut() {
            if let Some(mock_map) = mock_val.as_mapping_mut() {
                if let Some(entries_val) =
                    mock_map.get_mut(YamlValue::String("entries".to_string()))
                {
                    if let Some(entries_seq) = entries_val.as_sequence_mut() {
                        for entry in entries_seq.iter_mut() {
                            if let Some(entry_map) = entry.as_mapping_mut() {
                                // Remove internal fields (not for public YAML)
                                entry_map.remove(YamlValue::String("file_path".to_string()));
                                entry_map.remove(YamlValue::String("source_path".to_string()));
                                entry_map.remove(YamlValue::String("db_file_id".to_string()));
                            }
                        }
                    }
                }
            }
        }

        // For twin_list assets, strip entries from private field (keep only type and url)
        // Private entries must never appear in the public YAML
        if let Some(priv_val) = asset.private.as_mut() {
            if let Some(priv_map) = priv_val.as_mapping_mut() {
                // Remove entries array - this contains private file paths
                priv_map.remove(YamlValue::String("entries".to_string()));
                // Also remove any internal fields that shouldn't be public
                priv_map.remove(YamlValue::String("file_path".to_string()));
                priv_map.remove(YamlValue::String("db_file_id".to_string()));
            }
        }
    }
    manifest_for_write.extra.remove("mock_source_path");
    // Remove top-level file lists (used internally for file copying)
    manifest_for_write.extra.remove("mock_files");
    manifest_for_write.extra.remove("private_files");

    let yaml = serde_yaml::to_string(&manifest_for_write)?;
    let public_manifest_path = public_dir.join("dataset.yaml");
    storage.write_plaintext_file(&public_manifest_path, yaml.as_bytes(), true)?;

    let mut mapping_updates: Vec<(String, String)> = Vec::new();

    // Optionally copy mock artifacts alongside public manifest
    if copy_mock {
        let assets_dir = public_dir.join("assets");
        storage
            .ensure_dir(&assets_dir)
            .with_context(|| format!("Failed to create {}", assets_dir.display()))?;

        // Collect expected filenames from manifest to clean up stale files later
        let mut expected_files: std::collections::HashSet<String> =
            std::collections::HashSet::new();

        // First pass: collect all expected filenames
        for (asset_key, asset) in &manifest.assets {
            // CSV file for this asset
            expected_files.insert(format!("{}.csv", asset_key));

            // For twin_list, collect all mock entry filenames
            if asset.kind.as_deref() == Some("twin_list") {
                if let Some(mock_value) = &asset.mock {
                    let entries = mock_value
                        .as_mapping()
                        .and_then(|m| m.get(YamlValue::String("entries".to_string())))
                        .and_then(|v| v.as_sequence());
                    if let Some(entry_list) = entries {
                        for entry in entry_list {
                            if let Some(entry_map) = entry.as_mapping() {
                                // Get filename from URL or source_path
                                let filename = entry_map
                                    .get(YamlValue::String("url".to_string()))
                                    .and_then(|v| v.as_str())
                                    .and_then(|url| url.split('/').next_back())
                                    .or_else(|| {
                                        entry_map
                                            .get(YamlValue::String("source_path".to_string()))
                                            .and_then(|v| v.as_str())
                                            .and_then(|p| std::path::Path::new(p).file_name())
                                            .and_then(|n| n.to_str())
                                    });
                                if let Some(f) = filename {
                                    expected_files.insert(f.to_string());
                                }
                            }
                        }
                    }
                }
            }
        }

        // Clean up stale files in assets directory
        if assets_dir.exists() {
            if let Ok(entries) = fs::read_dir(&assets_dir) {
                for entry in entries.flatten() {
                    let filename = entry.file_name().to_string_lossy().to_string();
                    if !expected_files.contains(&filename) {
                        if let Err(e) = fs::remove_file(entry.path()) {
                            eprintln!("⚠️  Failed to remove stale file {}: {}", filename, e);
                        } else {
                            println!("  🗑️  Removed stale asset: {}", filename);
                        }
                    }
                }
            }
        }

        for (asset_key, asset) in &manifest.assets {
            if let (Some(priv_url), Some(priv_path)) = (
                manifest.private_url.as_ref(),
                resolve_private_source(&db, asset),
            ) {
                let private_fragment = format!("{}#assets.{}", priv_url, asset_key);
                mapping_updates.push((private_fragment, priv_path));
            }

            // Prefer an internal source hint to avoid copying from published URLs
            // Use mappings.mock -> (db_file_id or file_path) if available; fallback to mock field
            let mock_source = if let Some(mapping) = &asset.mappings {
                if let Some(ep) = &mapping.mock {
                    if let Some(id) = ep.db_file_id {
                        resolve_file_path_by_id(&db, id)
                            .ok()
                            .or_else(|| ep.file_path.clone())
                    } else {
                        ep.file_path.clone()
                    }
                } else {
                    None
                }
            } else if let Some(src) = asset
                .extra
                .get("mock_source_path")
                .and_then(|v| v.as_str())
                .map(|s| s.to_string())
            {
                Some(src)
            } else if let Some(mock_value) = &asset.mock {
                extract_mock_path(mock_value)
            } else {
                None
            };

            // Handle twin_list assets - copy each entry file and generate CSV
            if asset.kind.as_deref() == Some("twin_list") {
                if let Some(mock_value) = &asset.mock {
                    // Extract entries from mock object: { type, url, entries: [...] }
                    let entries = mock_value
                        .as_mapping()
                        .and_then(|m| m.get(YamlValue::String("entries".to_string())))
                        .and_then(|v| v.as_sequence());

                    if let Some(entry_list) = entries {
                        let mut csv_rows: Vec<String> = Vec::new();
                        csv_rows.push("participant_id,file".to_string());

                        for entry in entry_list {
                            let entry_map = match entry.as_mapping() {
                                Some(m) => m,
                                None => continue,
                            };

                            // Get participant_id (optional)
                            let participant_id = entry_map
                                .get(YamlValue::String("participant_id".to_string()))
                                .and_then(|v| v.as_str())
                                .unwrap_or("");

                            // Get source_path directly from entry (set by frontend)
                            let source_path_str = entry_map
                                .get(YamlValue::String("source_path".to_string()))
                                .and_then(|v| v.as_str());

                            // Extract filename from URL for CSV
                            let url = entry_map
                                .get(YamlValue::String("url".to_string()))
                                .and_then(|v| v.as_str())
                                .unwrap_or("");

                            let filename = url.split('/').next_back().unwrap_or("");

                            if let Some(src_path) = source_path_str {
                                if src_path.starts_with("http://")
                                    || src_path.starts_with("https://")
                                    || src_path.starts_with("syft://")
                                {
                                    // Already a URL, just use filename in CSV
                                    csv_rows.push(format!("{},{}", participant_id, filename));
                                    continue;
                                }

                                let source_path = PathBuf::from(src_path);
                                if source_path.exists() {
                                    let actual_filename = source_path
                                        .file_name()
                                        .and_then(|n| n.to_str())
                                        .unwrap_or(filename);
                                    let dest = assets_dir.join(actual_filename);
                                    if let Err(e) = fs::copy(&source_path, &dest) {
                                        eprintln!(
                                            "⚠️  Failed to copy {} to {}: {}",
                                            source_path.display(),
                                            dest.display(),
                                            e
                                        );
                                    } else {
                                        println!("  📁 Copied {} to assets/", actual_filename);
                                    }
                                    csv_rows
                                        .push(format!("{},{}", participant_id, actual_filename));
                                } else {
                                    eprintln!("⚠️  Mock file not found: {}", source_path.display());
                                }
                            } else if !filename.is_empty() {
                                // No source_path but we have a URL, entry is already published
                                csv_rows.push(format!("{},{}", participant_id, filename));
                            }
                        }

                        // Write CSV file
                        if csv_rows.len() > 1 {
                            let csv_path = assets_dir.join(format!("{}.csv", asset_key));
                            let csv_content = csv_rows.join("\n");
                            storage.write_plaintext_file(
                                &csv_path,
                                csv_content.as_bytes(),
                                true,
                            )?;
                            println!(
                                "  📄 Generated {} with {} entries",
                                csv_path.display(),
                                csv_rows.len() - 1
                            );
                        }
                    }
                }
                continue; // Skip single-file handling for twin_list
            }

            if let Some(src) = mock_source {
                if src.starts_with("http://")
                    || src.starts_with("https://")
                    || src.starts_with("syft://")
                {
                    continue;
                }
                let source_path = if Path::new(&src).is_relative() {
                    manifest_dir.join(&src)
                } else {
                    PathBuf::from(&src)
                };

                if !source_path.exists() {
                    eprintln!(
                        "⚠️  Mock for asset '{}' not found at {}",
                        asset_key,
                        source_path.display()
                    );
                    continue;
                }

                let dest = assets_dir.join(
                    source_path
                        .file_name()
                        .ok_or_else(|| anyhow!("Invalid mock file name"))?,
                );
                fs::copy(&source_path, &dest).with_context(|| {
                    format!(
                        "Failed to copy mock {} to {}",
                        source_path.display(),
                        dest.display()
                    )
                })?;
            }
        }
    } else {
        // Even when not copying mocks, still capture private mappings for local lookup
        for (asset_key, asset) in &manifest.assets {
            if let (Some(priv_url), Some(priv_path)) = (
                manifest.private_url.as_ref(),
                resolve_private_source(&db, asset),
            ) {
                let private_fragment = format!("{}#assets.{}", priv_url, asset_key);
                mapping_updates.push((private_fragment, priv_path));
            }
        }
    }

    if !mapping_updates.is_empty() {
        update_local_mappings(mapping_updates)?;
    }

    // Update datasets.yaml index
    let index_path = data_dir
        .join("datasites")
        .join(email)
        .join("public")
        .join("biovault")
        .join("datasets.yaml");
    if let Some(parent) = index_path.parent() {
        storage
            .ensure_dir(parent)
            .with_context(|| format!("Failed to create {}", parent.display()))?;
    }
    let mut index = match storage.read_plaintext_file(&index_path) {
        Ok(bytes) => serde_yaml::from_slice::<DatasetIndex>(&bytes).unwrap_or(DatasetIndex {
            email: email.clone(),
            resources: Vec::new(),
        }),
        Err(_) => DatasetIndex {
            email: email.clone(),
            resources: Vec::new(),
        },
    };

    let resource = DatasetResource {
        name: dataset_name.clone(),
        path: format!("./datasets/{}/dataset.yaml", dataset_name),
        schema: manifest
            .schema
            .clone()
            .unwrap_or_else(|| DEFAULT_SCHEMA.to_string()),
        version: manifest.version.clone(),
    };

    index.resources.retain(|r| r.name != dataset_name);
    index.resources.push(resource);
    index.resources.sort_by(|a, b| a.name.cmp(&b.name));

    let index_yaml = serde_yaml::to_string(&index)?;
    storage.write_plaintext_file(&index_path, index_yaml.as_bytes(), true)?;

    println!(
        "{}",
        format!(
            "✓ Published dataset '{}' to SyftBox\n  public: {}",
            manifest.name,
            public_manifest_path.display()
        )
        .green()
        .bold()
    );

    Ok(())
}

pub async fn list() -> Result<()> {
    let db = BioVaultDb::new()?;
    let datasets = list_datasets(&db)?;

    if datasets.is_empty() {
        println!("{}", "No datasets found.".yellow());
        return Ok(());
    }

    println!("{}", "Datasets".bold());
    println!();
    for ds in &datasets {
        println!("• {}", ds.name.cyan());
        println!("  version: {}", ds.version);
        println!("  author: {}", ds.author);
        if let Some(desc) = &ds.description {
            println!("  desc: {}", desc);
        }
        if let Some(pub_url) = &ds.public_url {
            println!("  public_url: {}", pub_url);
        }
        println!();
    }

    Ok(())
}

fn extract_mock_path(value: &YamlValue) -> Option<String> {
    match value {
        YamlValue::String(s) => Some(s.clone()),
        YamlValue::Mapping(map) => map
            .get(YamlValue::String("mock".to_string()))
            .and_then(|v| v.as_str())
            .map(|s| s.to_string()),
        _ => None,
    }
}

#[allow(clippy::too_many_arguments)]
pub async fn init(
    name: String,
    private_path: Option<String>,
    private_id: Option<i64>,
    mock_path: Option<String>,
    mock_id: Option<i64>,
    description: Option<String>,
    author: Option<String>,
    version: Option<String>,
    asset_key: Option<String>,
    asset_type: Option<String>,
    output: Option<String>,
) -> Result<()> {
    let config = Config::load()?;
    let email = &config.email;
    let name = name.trim();
    if name.is_empty() {
        return Err(anyhow!("Dataset name cannot be empty"));
    }

    let db = BioVaultDb::new()?;

    let resolved_private = match (private_id, private_path) {
        (Some(id), _) => resolve_file_path(&db, id)?,
        (None, Some(path)) => path,
        (None, None) => {
            return Err(anyhow!(
                "Provide either --private-id or a private file path argument"
            ))
        }
    };

    let resolved_mock = if let Some(id) = mock_id {
        Some(resolve_file_path(&db, id)?)
    } else {
        mock_path
    };

    let asset_key = asset_key
        .or_else(|| {
            Path::new(&resolved_private)
                .file_stem()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string())
        })
        .unwrap_or_else(|| "asset".to_string());

    let asset_kind = asset_type.unwrap_or_else(|| "twin".to_string());

    let mut manifest = DatasetManifest {
        name: name.to_string(),
        description,
        author: Some(author.unwrap_or_else(|| email.clone())),
        schema: Some(DEFAULT_SCHEMA.to_string()),
        version: Some(version.unwrap_or_else(|| "1.0.0".to_string())),
        http_relay_servers: vec!["syftbox.net".to_string()],
        public_url: Some(format!(
            "syft://{}/public/biovault/datasets/{}/dataset.yaml",
            email, name
        )),
        private_url: Some(format!(
            "syft://{}/private/biovault/datasets/{}/dataset.yaml",
            email, name
        )),
        ..Default::default()
    };

    // Private always points to the private syft URL suffix for this asset
    let private_field = Some(YamlValue::String("{url}.private".to_string()));

    // Mock: point to public asset location under this dataset if provided
    let mock_field = resolved_mock
        .as_ref()
        .and_then(|path| Path::new(path).file_name().and_then(|f| f.to_str()))
        .map(|filename| {
            format!(
                "syft://{}/public/biovault/datasets/{}/assets/{}",
                email, name, filename
            )
        });

    let asset = DatasetAsset {
        id: Some(Uuid::new_v4().to_string()),
        kind: Some(asset_kind),
        url: Some(format!("{{root.private_url}}#assets.{}", asset_key)),
        private: private_field,
        mock: mock_field.clone().map(YamlValue::String),
        mappings: Some(DatasetAssetMapping {
            private: Some(DatasetAssetMappingEndpoint {
                file_path: if private_id.is_none() {
                    Some(resolved_private.clone())
                } else {
                    None
                },
                db_file_id: private_id,
                entries: None,
            }),
            mock: resolved_mock.clone().map(|p| DatasetAssetMappingEndpoint {
                file_path: if mock_id.is_none() { Some(p) } else { None },
                db_file_id: mock_id,
                entries: None,
            }),
        }),
        extra: BTreeMap::new(),
    };

    manifest.assets.insert(asset_key.clone(), asset);

    let out_path = if let Some(out) = output {
        PathBuf::from(out)
    } else {
        PathBuf::from(format!("./{name}.dataset.yaml"))
    };

    let mut db = BioVaultDb::new()?;
    upsert_dataset(&mut db, &manifest)?;

    let yaml = serde_yaml::to_string(&manifest)?;
    fs::write(&out_path, yaml)
        .with_context(|| format!("Failed to write {}", out_path.display()))?;

    println!(
        "{}",
        format!(
            "✓ Generated dataset manifest at {}\n  name: {}\n  asset key: {}\n  private: {}\n  mock: {}",
            out_path.display(),
            manifest.name,
            asset_key,
            "{url}.private",
            mock_field
                .as_deref()
                .unwrap_or("syft://<email>/public/biovault/datasets/<name>/assets/<mock>")
        )
        .green()
        .bold()
    );

    println!(
        "{}",
        "Edit the manifest to refine metadata or add more assets, then run: bv datasets publish <path> --copy-mock"
            .dimmed()
    );

    Ok(())
}

fn resolve_file_path(db: &BioVaultDb, file_id: i64) -> Result<String> {
    let path: String = db
        .connection()
        .query_row(
            "SELECT file_path FROM files WHERE id = ?1",
            [file_id],
            |row| row.get(0),
        )
        .with_context(|| format!("File with id {} not found in catalog", file_id))?;
    Ok(path)
}

fn resolve_file_path_by_id(db: &BioVaultDb, file_id: i64) -> Result<String> {
    resolve_file_path(db, file_id)
}

fn resolve_private_source(db: &BioVaultDb, asset: &DatasetAsset) -> Option<String> {
    let mapping = asset.mappings.as_ref()?;
    if let Some(ep) = &mapping.private {
        if let Some(id) = ep.db_file_id {
            return resolve_file_path_by_id(db, id)
                .ok()
                .or_else(|| ep.file_path.clone());
        }
        if ep.file_path.is_some() {
            return ep.file_path.clone();
        }
    }
    None
}

#[derive(Debug, Default, Serialize, Deserialize)]
#[allow(dead_code)]
struct LocalMappingFile {
    #[serde(default)]
    mappings: BTreeMap<String, String>,
}

fn update_local_mappings(entries: Vec<(String, String)>) -> Result<()> {
    if entries.is_empty() {
        return Ok(());
    }

    crate::data::datasets::update_local_mappings(&entries)
}
