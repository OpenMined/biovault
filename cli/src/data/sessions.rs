//! Database operations for session datasets
//!
//! Manages the association between sessions and datasets for collaborative analysis.

use crate::data::BioVaultDb;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// A dataset associated with a session
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionDataset {
    pub id: i64,
    pub session_id: String,
    pub dataset_public_url: String,
    pub dataset_owner: String,
    pub dataset_name: String,
    pub role: String, // 'shared' or 'yours'
    pub created_at: String,
}

/// Request to associate a dataset with a session
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AddSessionDatasetRequest {
    pub session_id: String,
    pub dataset_public_url: String,
    pub dataset_owner: String,
    pub dataset_name: String,
    pub role: Option<String>, // defaults to 'shared'
}

/// Add a dataset to a session
pub fn add_session_dataset(db: &BioVaultDb, request: &AddSessionDatasetRequest) -> Result<i64> {
    let role = request.role.as_deref().unwrap_or("shared");

    db.connection().execute(
        "INSERT INTO session_datasets (session_id, dataset_public_url, dataset_owner, dataset_name, role)
         VALUES (?1, ?2, ?3, ?4, ?5)
         ON CONFLICT(session_id, dataset_public_url) DO UPDATE SET
             role = excluded.role",
        rusqlite::params![
            &request.session_id,
            &request.dataset_public_url,
            &request.dataset_owner,
            &request.dataset_name,
            role,
        ],
    )?;

    Ok(db.connection().last_insert_rowid())
}

/// Get all datasets associated with a session
pub fn get_session_datasets(db: &BioVaultDb, session_id: &str) -> Result<Vec<SessionDataset>> {
    let mut stmt = db.connection().prepare(
        "SELECT id, session_id, dataset_public_url, dataset_owner, dataset_name, role, created_at
         FROM session_datasets
         WHERE session_id = ?1
         ORDER BY created_at ASC",
    )?;

    let datasets = stmt
        .query_map([session_id], |row| {
            Ok(SessionDataset {
                id: row.get(0)?,
                session_id: row.get(1)?,
                dataset_public_url: row.get(2)?,
                dataset_owner: row.get(3)?,
                dataset_name: row.get(4)?,
                role: row.get(5)?,
                created_at: row.get(6)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()
        .context("Failed to query session datasets")?;

    Ok(datasets)
}

/// Remove a dataset from a session
pub fn remove_session_dataset(
    db: &BioVaultDb,
    session_id: &str,
    dataset_public_url: &str,
) -> Result<bool> {
    let affected = db.connection().execute(
        "DELETE FROM session_datasets WHERE session_id = ?1 AND dataset_public_url = ?2",
        rusqlite::params![session_id, dataset_public_url],
    )?;

    Ok(affected > 0)
}

/// Remove all datasets from a session
pub fn clear_session_datasets(db: &BioVaultDb, session_id: &str) -> Result<usize> {
    let affected = db.connection().execute(
        "DELETE FROM session_datasets WHERE session_id = ?1",
        [session_id],
    )?;

    Ok(affected)
}

/// Check if a session has any associated datasets
pub fn session_has_datasets(db: &BioVaultDb, session_id: &str) -> Result<bool> {
    let count: i64 = db.connection().query_row(
        "SELECT COUNT(*) FROM session_datasets WHERE session_id = ?1",
        [session_id],
        |row| row.get(0),
    )?;

    Ok(count > 0)
}

/// Get sessions that use a specific dataset
pub fn get_sessions_for_dataset(
    db: &BioVaultDb,
    dataset_owner: &str,
    dataset_name: &str,
) -> Result<Vec<String>> {
    let mut stmt = db.connection().prepare(
        "SELECT DISTINCT session_id FROM session_datasets
         WHERE dataset_owner = ?1 AND dataset_name = ?2",
    )?;

    let session_ids = stmt
        .query_map(rusqlite::params![dataset_owner, dataset_name], |row| {
            row.get(0)
        })?
        .collect::<Result<Vec<_>, _>>()
        .context("Failed to query sessions for dataset")?;

    Ok(session_ids)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn setup_test_db() -> (TempDir, BioVaultDb) {
        let temp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(temp.path().to_path_buf());
        let db = BioVaultDb::new().unwrap();
        crate::config::clear_test_biovault_home();
        (temp, db)
    }

    #[test]
    fn test_add_and_get_session_datasets() {
        let (_temp, db) = setup_test_db();

        // First create a session
        db.connection()
            .execute(
                "INSERT INTO sessions (session_id, name, session_path, owner, status)
                 VALUES ('test123', 'Test Session', '/tmp/test', 'user@test.com', 'active')",
                [],
            )
            .unwrap();

        // Add a dataset
        let request = AddSessionDatasetRequest {
            session_id: "test123".to_string(),
            dataset_public_url:
                "syft://owner@test.com/public/biovault/datasets/mydata/dataset.yaml".to_string(),
            dataset_owner: "owner@test.com".to_string(),
            dataset_name: "mydata".to_string(),
            role: Some("shared".to_string()),
        };

        let id = add_session_dataset(&db, &request).unwrap();
        assert!(id > 0);

        // Get datasets
        let datasets = get_session_datasets(&db, "test123").unwrap();
        assert_eq!(datasets.len(), 1);
        assert_eq!(datasets[0].dataset_name, "mydata");
        assert_eq!(datasets[0].role, "shared");
    }

    #[test]
    fn test_remove_session_dataset() {
        let (_temp, db) = setup_test_db();

        // Setup
        db.connection()
            .execute(
                "INSERT INTO sessions (session_id, name, session_path, owner, status)
                 VALUES ('test456', 'Test Session', '/tmp/test2', 'user@test.com', 'active')",
                [],
            )
            .unwrap();

        let request = AddSessionDatasetRequest {
            session_id: "test456".to_string(),
            dataset_public_url: "syft://owner@test.com/public/biovault/datasets/data/dataset.yaml"
                .to_string(),
            dataset_owner: "owner@test.com".to_string(),
            dataset_name: "data".to_string(),
            role: None,
        };
        add_session_dataset(&db, &request).unwrap();

        // Remove
        let removed = remove_session_dataset(
            &db,
            "test456",
            "syft://owner@test.com/public/biovault/datasets/data/dataset.yaml",
        )
        .unwrap();
        assert!(removed);

        // Verify removed
        let datasets = get_session_datasets(&db, "test456").unwrap();
        assert!(datasets.is_empty());
    }
}
