use crate::cli::commands::datasets::{DatasetAsset, DatasetManifest};
use crate::config::get_biovault_home;
use crate::data::db::BioVaultDb;
use anyhow::{Context, Result};
use rusqlite::OptionalExtension;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DatasetRecord {
    pub id: i64,
    pub name: String,
    pub version: String,
    pub author: String,
    pub description: Option<String>,
    pub schema: String,
    pub public_url: Option<String>,
    pub private_url: Option<String>,
    pub http_relay_servers: Vec<String>,
    pub extra: serde_json::Value,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DatasetAssetRecord {
    pub asset_key: String,
    pub asset_uuid: String,
    pub kind: String,
    pub url: String,
    pub private_ref: Option<String>,
    pub mock_ref: Option<String>,
    pub private_file_id: Option<i64>,
    pub mock_file_id: Option<i64>,
    pub private_path: Option<String>,
    pub mock_path: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
struct LocalMappingFile {
    #[serde(default)]
    mappings: BTreeMap<String, String>,
}

/// Remove a dataset (and its assets via cascade) by name. Returns rows deleted (0 or 1).
pub fn delete_dataset(db: &BioVaultDb, name: &str) -> Result<usize> {
    let affected = db
        .conn
        .execute("DELETE FROM datasets WHERE name = ?1", [name])?;
    Ok(affected)
}

pub fn upsert_dataset(db: &mut BioVaultDb, manifest: &DatasetManifest) -> Result<i64> {
    let tx = db.conn.unchecked_transaction()?;

    let http_relay_servers_json = serde_json::to_string(&manifest.http_relay_servers)
        .context("Failed to serialize http_relay_servers")?;
    let extra_json = serde_json::to_value(&manifest.extra).unwrap_or(serde_json::Value::Null);

    let dataset_id: i64 = tx
        .query_row(
            "SELECT id FROM datasets WHERE name = ?1",
            [&manifest.name],
            |row| row.get(0),
        )
        .optional()?
        .unwrap_or_else(|| {
            tx.execute(
                "INSERT INTO datasets (name, version, author, description, schema, public_url, private_url, http_relay_servers, extra)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
                rusqlite::params![
                    manifest.name,
                    manifest.version.clone().unwrap_or_else(|| "1.0.0".into()),
                    manifest.author.clone().unwrap_or_else(|| "".into()),
                    manifest.description,
                    manifest.schema.clone().unwrap_or_else(|| "".into()),
                    manifest.public_url.clone(),
                    manifest.private_url.clone(),
                    http_relay_servers_json,
                    serde_json::to_string(&extra_json).unwrap_or_else(|_| "{}".into()),
                ],
            )
            .expect("insert dataset");
            tx.last_insert_rowid()
        });

    if dataset_id != 0 {
        tx.execute(
            "UPDATE datasets SET version = ?1, author = ?2, description = ?3, schema = ?4, public_url = ?5,
             private_url = ?6, http_relay_servers = ?7, extra = ?8, updated_at = CURRENT_TIMESTAMP WHERE id = ?9",
            rusqlite::params![
                manifest.version.clone().unwrap_or_else(|| "1.0.0".into()),
                manifest.author.clone().unwrap_or_else(|| "".into()),
                manifest.description,
                manifest.schema.clone().unwrap_or_else(|| "".into()),
                manifest.public_url.clone(),
                manifest.private_url.clone(),
                http_relay_servers_json,
                serde_json::to_string(&extra_json).unwrap_or_else(|_| "{}".into()),
                dataset_id
            ],
        )?;
    }

    // Replace assets
    tx.execute(
        "DELETE FROM dataset_assets WHERE dataset_id = ?1",
        [dataset_id],
    )?;

    for (asset_key, asset) in &manifest.assets {
        let extra_json = serde_json::to_value(&asset.extra).unwrap_or(serde_json::Value::Null);
        tx.execute(
            "INSERT INTO dataset_assets
                (dataset_id, asset_key, asset_uuid, kind, url, private_ref, mock_ref, extra, private_file_id, mock_file_id, private_path, mock_path)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12)",
            rusqlite::params![
                dataset_id,
                asset_key,
                asset.id.clone().unwrap_or_else(|| "".into()),
                asset.kind.clone().unwrap_or_else(|| "file".into()),
                asset.url.clone().unwrap_or_else(|| "".into()),
                asset
                    .private
                    .clone()
                    .and_then(|v| {
                        // For single assets: extract string directly
                        // For twin_list: serialize the object to JSON
                        if let Some(s) = v.as_str() {
                            Some(s.to_string())
                        } else {
                            serde_json::to_string(&v).ok()
                        }
                    }),
                asset
                    .mock
                    .clone()
                    .and_then(|v| {
                        // For single assets: extract string directly
                        // For twin_list: serialize the object to JSON
                        if let Some(s) = v.as_str() {
                            Some(s.to_string())
                        } else {
                            serde_json::to_string(&v).ok()
                        }
                    }),
                serde_json::to_string(&extra_json).unwrap_or_else(|_| "{}".into()),
                asset
                    .mappings
                    .as_ref()
                    .and_then(|m| m.private.as_ref())
                    .and_then(|p| p.db_file_id),
                asset
                    .mappings
                    .as_ref()
                    .and_then(|m| m.mock.as_ref())
                    .and_then(|p| p.db_file_id),
                asset
                    .mappings
                    .as_ref()
                    .and_then(|m| m.private.as_ref())
                    .and_then(|p| p.file_path.clone()),
                asset
                    .mappings
                    .as_ref()
                    .and_then(|m| m.mock.as_ref())
                    .and_then(|p| p.file_path.clone()),
            ],
        )?;
    }

    tx.commit()?;
    Ok(dataset_id)
}

pub fn get_dataset_with_assets(
    db: &BioVaultDb,
    name: &str,
) -> Result<Option<(DatasetRecord, Vec<DatasetAssetRecord>)>> {
    let conn = db.connection();
    let dataset: Option<DatasetRecord> = conn
        .query_row(
            "SELECT id, name, version, author, description, schema, public_url, private_url, http_relay_servers, extra
             FROM datasets WHERE name = ?1",
            [name],
            |row| {
                let http_json: String = row.get(8)?;
                let http_relay_servers: Vec<String> =
                    serde_json::from_str(&http_json).unwrap_or_default();
                let extra_json: String = row.get(9).unwrap_or_else(|_| "{}".into());
                Ok(DatasetRecord {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    version: row.get(2)?,
                    author: row.get(3)?,
                    description: row.get(4)?,
                    schema: row.get(5)?,
                    public_url: row.get(6)?,
                    private_url: row.get(7)?,
                    http_relay_servers,
                    extra: serde_json::from_str(&extra_json).unwrap_or_else(|_| serde_json::json!({})),
                })
            },
        )
        .optional()?;

    let Some(dataset) = dataset else {
        return Ok(None);
    };

    let mut stmt = conn.prepare(
        "SELECT asset_key, asset_uuid, kind, url, private_ref, mock_ref, private_file_id, mock_file_id, private_path, mock_path
         FROM dataset_assets WHERE dataset_id = ?1",
    )?;
    let assets_iter = stmt.query_map([dataset.id], |row| {
        Ok(DatasetAssetRecord {
            asset_key: row.get(0)?,
            asset_uuid: row.get(1)?,
            kind: row.get(2)?,
            url: row.get(3)?,
            private_ref: row.get(4)?,
            mock_ref: row.get(5)?,
            private_file_id: row.get::<_, Option<i64>>(6).ok().flatten().or_else(|| {
                row.get::<_, Option<String>>(6)
                    .ok()
                    .flatten()
                    .and_then(|s| s.parse().ok())
            }),
            mock_file_id: row.get::<_, Option<i64>>(7).ok().flatten().or_else(|| {
                row.get::<_, Option<String>>(7)
                    .ok()
                    .flatten()
                    .and_then(|s| s.parse().ok())
            }),
            private_path: row.get(8).ok(),
            mock_path: row.get(9).ok(),
        })
    })?;
    let mut assets = Vec::new();
    for a in assets_iter {
        assets.push(a?);
    }

    Ok(Some((dataset, assets)))
}

pub fn list_datasets(db: &BioVaultDb) -> Result<Vec<DatasetRecord>> {
    let conn = db.connection();
    let mut stmt = conn.prepare(
        "SELECT id, name, version, author, description, schema, public_url, private_url, http_relay_servers, extra
         FROM datasets ORDER BY name",
    )?;
    let rows = stmt.query_map([], |row| {
        let http_json: String = row.get(8)?;
        let http_relay_servers: Vec<String> = serde_json::from_str(&http_json).unwrap_or_default();
        let extra_json: String = row.get(9).unwrap_or_else(|_| "{}".into());
        Ok(DatasetRecord {
            id: row.get(0)?,
            name: row.get(1)?,
            version: row.get(2)?,
            author: row.get(3)?,
            description: row.get(4)?,
            schema: row.get(5)?,
            public_url: row.get(6)?,
            private_url: row.get(7)?,
            http_relay_servers,
            extra: serde_json::from_str(&extra_json).unwrap_or_else(|_| serde_json::json!({})),
        })
    })?;
    let mut out = Vec::new();
    for row in rows {
        out.push(row?);
    }
    Ok(out)
}

/// List all datasets with their assets for UI consumption.
pub fn list_datasets_with_assets(
    db: &BioVaultDb,
) -> Result<Vec<(DatasetRecord, Vec<DatasetAssetRecord>)>> {
    let datasets = list_datasets(db)?;
    let mut out = Vec::with_capacity(datasets.len());
    for ds in datasets {
        let mut stmt = db.conn.prepare(
            "SELECT asset_key, asset_uuid, kind, url, private_ref, mock_ref, private_file_id, mock_file_id, private_path, mock_path
             FROM dataset_assets WHERE dataset_id = ?1",
        )?;
        let assets_iter = stmt.query_map([ds.id], |row| {
            Ok(DatasetAssetRecord {
                asset_key: row.get(0)?,
                asset_uuid: row.get(1)?,
                kind: row.get(2)?,
                url: row.get(3)?,
                private_ref: row.get(4)?,
                mock_ref: row.get(5)?,
                private_file_id: row.get(6).ok(),
                mock_file_id: row.get(7).ok(),
                private_path: row.get(8).ok(),
                mock_path: row.get(9).ok(),
            })
        })?;

        let mut assets = Vec::new();
        for asset in assets_iter {
            assets.push(asset?);
        }

        out.push((ds.clone(), assets));
    }

    Ok(out)
}

/// Persist private URL -> local path mappings to BIOVAULT_HOME/mapping.yaml
pub fn update_local_mappings(entries: &[(String, String)]) -> Result<()> {
    if entries.is_empty() {
        return Ok(());
    }

    let mapping_path = get_biovault_home()?.join("mapping.yaml");
    let mut mapping = if mapping_path.exists() {
        let raw = fs::read_to_string(&mapping_path)
            .with_context(|| format!("Failed to read {}", mapping_path.display()))?;
        serde_yaml::from_str::<LocalMappingFile>(&raw).unwrap_or_default()
    } else {
        LocalMappingFile::default()
    };

    for (public_ref, private_path) in entries {
        mapping
            .mappings
            .insert(public_ref.clone(), private_path.clone());
    }

    let yaml = serde_yaml::to_string(&mapping)?;
    if let Some(parent) = mapping_path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create {}", parent.display()))?;
    }
    fs::write(&mapping_path, yaml)
        .with_context(|| format!("Failed to write {}", mapping_path.display()))?;

    Ok(())
}

/// Load mappings from BIOVAULT_HOME/mapping.yaml
fn load_local_mappings() -> Result<BTreeMap<String, String>> {
    let mapping_path = get_biovault_home()?.join("mapping.yaml");
    if mapping_path.exists() {
        let raw = fs::read_to_string(&mapping_path)
            .with_context(|| format!("Failed to read {}", mapping_path.display()))?;
        let mapping: LocalMappingFile = serde_yaml::from_str(&raw).unwrap_or_default();
        Ok(mapping.mappings)
    } else {
        Ok(BTreeMap::new())
    }
}

/// Resolve a syft:// or file:// URL to a local filesystem path.
///
/// For file:// URLs: strips prefix and returns the path directly
/// For public syft:// URLs (containing "/public/"): resolves directly to the datasite path
/// For private syft:// URLs with #fragment: looks up in mapping.yaml
///
/// Example URLs:
/// - `file:///path/to/file.txt` → /path/to/file.txt
/// - `syft://user@example.com/public/biovault/datasets/test/assets/file.txt` → datasite path
/// - `syft://user@example.com/private/biovault/datasets/test/dataset.yaml#assets.asset_1.private` → mapping.yaml lookup
pub fn resolve_syft_url(data_dir: &std::path::Path, url: &str) -> Result<std::path::PathBuf> {
    use syftbox_sdk::SyftURL;

    // Handle file:// URLs - just strip the prefix
    if url.starts_with("file://") {
        let path = url.strip_prefix("file://").unwrap_or(url);
        return Ok(std::path::PathBuf::from(path));
    }

    let parsed = SyftURL::parse(url)
        .map_err(|e| anyhow::anyhow!("Failed to parse syft URL '{}': {}", url, e))?;

    // Check if URL has a fragment (indicates private reference needing mapping lookup)
    if let Some(ref _fragment) = parsed.fragment {
        // Look up in mapping.yaml using the full URL with fragment
        let mappings = load_local_mappings()?;
        if let Some(local_path) = mappings.get(url) {
            return Ok(std::path::PathBuf::from(local_path));
        }
        // If not in mappings, fall through to direct resolution
    }

    // Direct resolution: syft://email/path → data_dir/datasites/email/path
    Ok(data_dir
        .join("datasites")
        .join(&parsed.email)
        .join(&parsed.path))
}

/// Batch resolve multiple syft:// URLs to local paths
pub fn resolve_syft_urls(
    data_dir: &std::path::Path,
    urls: &[String],
) -> Result<Vec<(String, Option<String>)>> {
    let mut results = Vec::with_capacity(urls.len());
    for url in urls {
        match resolve_syft_url(data_dir, url) {
            Ok(path) => {
                if path.exists() {
                    results.push((url.clone(), Some(path.to_string_lossy().to_string())));
                } else {
                    results.push((url.clone(), None));
                }
            }
            Err(_) => {
                results.push((url.clone(), None));
            }
        }
    }
    Ok(results)
}

/// Convert DB rows into manifest structs so publish can emit YAML
pub fn build_manifest_from_db(
    dataset: &DatasetRecord,
    assets: &[DatasetAssetRecord],
) -> DatasetManifest {
    let mut manifest = DatasetManifest {
        name: dataset.name.clone(),
        version: Some(dataset.version.clone()),
        author: Some(dataset.author.clone()),
        description: dataset.description.clone(),
        schema: Some(dataset.schema.clone()),
        public_url: dataset.public_url.clone(),
        private_url: dataset.private_url.clone(),
        http_relay_servers: dataset.http_relay_servers.clone(),
        extra: serde_json::from_value(dataset.extra.clone()).unwrap_or_default(),
        ..Default::default()
    };

    for a in assets {
        let asset = DatasetAsset {
            id: Some(a.asset_uuid.clone()),
            kind: Some(a.kind.clone()),
            url: Some(a.url.clone()),
            private: a.private_ref.clone().map(|s| {
                // Try to parse as JSON (for twin_list objects), fallback to string
                serde_json::from_str::<serde_yaml::Value>(&s)
                    .unwrap_or_else(|_| serde_yaml::Value::String(s))
            }),
            mock: a.mock_ref.clone().map(|s| {
                // Try to parse as JSON (for twin_list objects), fallback to string
                serde_json::from_str::<serde_yaml::Value>(&s)
                    .unwrap_or_else(|_| serde_yaml::Value::String(s))
            }),
            mappings: Some(crate::cli::commands::datasets::DatasetAssetMapping {
                private: Some(
                    crate::cli::commands::datasets::DatasetAssetMappingEndpoint {
                        file_path: a.private_path.clone(),
                        db_file_id: a.private_file_id,
                        entries: None,
                    },
                ),
                mock: Some(
                    crate::cli::commands::datasets::DatasetAssetMappingEndpoint {
                        file_path: a.mock_path.clone(),
                        db_file_id: a.mock_file_id,
                        entries: None,
                    },
                ),
            }),
            ..Default::default()
        };
        manifest.assets.insert(a.asset_key.clone(), asset);
    }

    manifest
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::commands::datasets::{
        DatasetAssetMapping, DatasetAssetMappingEndpoint, DatasetManifest,
    };
    use crate::config;
    use tempfile::TempDir;

    #[test]
    fn upsert_and_load_persists_asset_links() {
        let temp = TempDir::new().unwrap();
        config::set_test_biovault_home(temp.path());

        let mut db = BioVaultDb::new().unwrap();

        let mut manifest = DatasetManifest {
            name: "single_cell".to_string(),
            author: Some("test@example.com".to_string()),
            schema: Some("net.biovault.datasets:1.0.0".to_string()),
            version: Some("1.0.0".to_string()),
            http_relay_servers: vec!["syftbox.net".to_string()],
            ..Default::default()
        };

        let asset = crate::cli::commands::datasets::DatasetAsset {
            id: Some("asset-uuid".to_string()),
            kind: Some("twin".to_string()),
            url: Some("{root.private_url}#assets.sc_rnaseq".to_string()),
            private: Some(serde_yaml::Value::String("{url}.private".to_string())),
            mock: Some(serde_yaml::Value::String("syft://public/mock".to_string())),
            mappings: Some(DatasetAssetMapping {
                private: Some(DatasetAssetMappingEndpoint {
                    file_path: Some("/path/private.h5ad".to_string()),
                    db_file_id: Some(2),
                    entries: None,
                }),
                mock: Some(DatasetAssetMappingEndpoint {
                    file_path: Some("/path/mock.h5ad".to_string()),
                    db_file_id: Some(1),
                    entries: None,
                }),
            }),
            ..Default::default()
        };

        manifest.assets.insert("sc_rnaseq".to_string(), asset);

        // Upsert and reload
        upsert_dataset(&mut db, &manifest).unwrap();
        let (ds, assets) = get_dataset_with_assets(&db, "single_cell")
            .unwrap()
            .expect("dataset present");

        // Verify asset link columns round-trip
        assert_eq!(assets.len(), 1);
        let a = &assets[0];
        assert_eq!(a.private_file_id, Some(2));
        assert_eq!(a.mock_file_id, Some(1));
        assert_eq!(a.private_path.as_deref(), Some("/path/private.h5ad"));
        assert_eq!(a.mock_path.as_deref(), Some("/path/mock.h5ad"));

        // Build manifest from DB and ensure mappings exist for internal use
        let rebuilt = build_manifest_from_db(&ds, &assets);
        let rebuilt_asset = rebuilt.assets.get("sc_rnaseq").unwrap();
        let mapping = rebuilt_asset.mappings.as_ref().unwrap();
        assert_eq!(mapping.private.as_ref().unwrap().db_file_id, Some(2));
        assert_eq!(mapping.mock.as_ref().unwrap().db_file_id, Some(1));

        config::clear_test_biovault_home();
    }
}
