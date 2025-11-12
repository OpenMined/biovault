use anyhow::{bail, Result};
use rusqlite::params;
use serde::{Deserialize, Serialize};

use super::{BioVaultDb, FileRecord};

#[derive(Debug, Serialize, Deserialize)]
pub struct CollectionRecord {
    pub id: i64,
    pub name: String,
    pub description: Option<String>,
    pub variable_name: String,
    pub file_count: usize,
    pub created_at: String,
    pub updated_at: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CollectionDetail {
    pub id: i64,
    pub name: String,
    pub description: Option<String>,
    pub variable_name: String,
    pub files: Vec<FileRecord>,
    pub created_at: String,
    pub updated_at: String,
}

/// Generate a snake_case variable name from a display name
pub fn generate_variable_name(name: &str) -> String {
    name.trim()
        .to_lowercase()
        .chars()
        .map(|c| match c {
            'a'..='z' | '0'..='9' => c,
            ' ' | '-' | '_' => '_',
            _ => '\0', // Skip other characters
        })
        .collect::<String>()
        .split('_')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("_")
        .trim_matches('_')
        .to_string()
}

/// Validate that a variable name is valid snake_case
pub fn validate_variable_name(var_name: &str) -> Result<()> {
    if var_name.is_empty() {
        bail!("Variable name cannot be empty");
    }

    if !var_name.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_') {
        bail!("Variable name must be snake_case (lowercase letters, digits, underscores only)");
    }

    if var_name.starts_with('_') || var_name.ends_with('_') {
        bail!("Variable name cannot start or end with underscore");
    }

    if var_name.contains("__") {
        bail!("Variable name cannot contain consecutive underscores");
    }

    if var_name.chars().next().map(|c| c.is_ascii_digit()).unwrap_or(false) {
        bail!("Variable name cannot start with a digit");
    }

    Ok(())
}

/// Create a new collection
pub fn create_collection(
    db: &BioVaultDb,
    name: String,
    description: Option<String>,
    variable_name: Option<String>,
) -> Result<CollectionRecord> {
    let var_name = if let Some(vn) = variable_name {
        validate_variable_name(&vn)?;
        vn
    } else {
        let generated = generate_variable_name(&name);
        if generated.is_empty() {
            bail!("Could not generate a valid variable name from '{}'", name);
        }
        generated
    };

    // Check if variable_name already exists
    let exists: bool = db
        .conn
        .query_row(
            "SELECT COUNT(*) FROM collections WHERE variable_name = ?1",
            params![var_name],
            |row| Ok(row.get::<_, i32>(0)? > 0),
        )
        .unwrap_or(false);

    if exists {
        bail!("Collection with variable_name '{}' already exists", var_name);
    }

    db.conn.execute(
        "INSERT INTO collections (name, description, variable_name) VALUES (?1, ?2, ?3)",
        params![name, description, var_name],
    )?;

    let id = db.conn.last_insert_rowid();

    // Get the created collection with file count
    get_collection_by_id(db, id)
}

/// Get collection by ID
pub fn get_collection_by_id(db: &BioVaultDb, id: i64) -> Result<CollectionRecord> {
    let (name, description, variable_name, created_at, updated_at): (
        String,
        Option<String>,
        String,
        String,
        String,
    ) = db.conn.query_row(
        "SELECT name, description, variable_name, created_at, updated_at FROM collections WHERE id = ?1",
        params![id],
        |row| {
            Ok((
                row.get(0)?,
                row.get(1)?,
                row.get(2)?,
                row.get(3)?,
                row.get(4)?,
            ))
        },
    )?;

    let file_count: i32 = db
        .conn
        .query_row(
            "SELECT COUNT(*) FROM collection_files WHERE collection_id = ?1",
            params![id],
            |row| row.get(0),
        )
        .unwrap_or(0);

    Ok(CollectionRecord {
        id,
        name,
        description,
        variable_name,
        file_count: file_count as usize,
        created_at,
        updated_at,
    })
}

/// Get collection by variable_name or name
pub fn get_collection(db: &BioVaultDb, identifier: &str) -> Result<CollectionRecord> {
    // Try variable_name first (more specific)
    if let Ok(id) = db.conn.query_row(
        "SELECT id FROM collections WHERE variable_name = ?1",
        params![identifier],
        |row| row.get(0),
    ) {
        return get_collection_by_id(db, id);
    }

    // Try name
    if let Ok(id) = db.conn.query_row(
        "SELECT id FROM collections WHERE name = ?1",
        params![identifier],
        |row| row.get(0),
    ) {
        return get_collection_by_id(db, id);
    }

    bail!("Collection '{}' not found", identifier)
}

/// List all collections
pub fn list_collections(db: &BioVaultDb) -> Result<Vec<CollectionRecord>> {
    let mut stmt = db.conn.prepare(
        "SELECT id, name, description, variable_name, created_at, updated_at FROM collections ORDER BY name",
    )?;

    let mut collections = Vec::new();
    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            row.get::<_, String>(1)?,
            row.get::<_, Option<String>>(2)?,
            row.get::<_, String>(3)?,
            row.get::<_, String>(4)?,
            row.get::<_, String>(5)?,
        ))
    })?;

    for row in rows {
        let (id, name, description, variable_name, created_at, updated_at) = row?;

        // Get file count
        let file_count: i32 = db
            .conn
            .query_row(
                "SELECT COUNT(*) FROM collection_files WHERE collection_id = ?1",
                params![id],
                |row| row.get(0),
            )
            .unwrap_or(0);

        collections.push(CollectionRecord {
            id,
            name,
            description,
            variable_name,
            file_count: file_count as usize,
            created_at,
            updated_at,
        });
    }

    Ok(collections)
}

/// Get collection with all files
pub fn get_collection_detail(db: &BioVaultDb, identifier: &str) -> Result<CollectionDetail> {
    let collection = get_collection(db, identifier)?;

    // Get all files in this collection
    let mut stmt = db.conn.prepare(
        "SELECT f.id, f.file_path, f.file_hash, f.file_type, f.file_size, f.data_type,
                f.metadata, f.status, f.processing_error, f.created_at, f.updated_at,
                p.participant_id, p.inferred_sex as participant_inferred_sex
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
         INNER JOIN collection_files cf ON f.id = cf.file_id
         WHERE cf.collection_id = ?1
         ORDER BY f.file_path",
    )?;

    // Collect file IDs first, then fetch genotype metadata separately
    let mut file_rows = Vec::new();
    let files_iter = stmt.query_map(params![collection.id], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            row.get::<_, String>(1)?,
            row.get::<_, String>(2)?,
            row.get::<_, Option<String>>(3)?,
            row.get::<_, Option<i64>>(4)?,
            row.get::<_, Option<String>>(5)?,
            row.get::<_, Option<String>>(6)?,
            row.get::<_, Option<String>>(7)?,
            row.get::<_, Option<String>>(8)?,
            row.get::<_, String>(9)?,
            row.get::<_, String>(10)?,
            row.get::<_, Option<String>>(11)?,
            row.get::<_, Option<String>>(12)?,
        ))
    })?;

    for row in files_iter {
        file_rows.push(row?);
    }

    // Now fetch genotype metadata for genotype files
    let mut files = Vec::new();
    for (id, file_path, file_hash, file_type, file_size, data_type, _metadata, status, processing_error, created_at, updated_at, participant_id, participant_inferred_sex) in file_rows {
        let (source, grch_version, row_count, chromosome_count, inferred_sex) =
            if data_type.as_deref() == Some("Genotype") {
                db.conn
                    .query_row(
                        "SELECT source, grch_version, row_count, chromosome_count, inferred_sex
                         FROM genotype_metadata WHERE file_id = ?1",
                        params![id],
                        |row| {
                            Ok((
                                row.get(0)?,
                                row.get(1)?,
                                row.get(2)?,
                                row.get(3)?,
                                row.get(4)?,
                            ))
                        },
                    )
                    .ok()
                    .unwrap_or((None, None, None, None, None))
            } else {
                (None, None, None, None, None)
            };

        files.push(FileRecord {
            id,
            file_path,
            file_hash,
            file_type,
            file_size: file_size.map(|s| s as u64),
            data_type,
            source,
            grch_version,
            row_count,
            chromosome_count,
            inferred_sex: inferred_sex.or(participant_inferred_sex),
            status,
            processing_error,
            participant_id,
            participant_name: None, // Could be enhanced to fetch participant name
            created_at,
            updated_at,
        });
    }

    Ok(CollectionDetail {
        id: collection.id,
        name: collection.name,
        description: collection.description,
        variable_name: collection.variable_name,
        files,
        created_at: collection.created_at,
        updated_at: collection.updated_at,
    })
}

/// Add files to a collection
pub fn add_files_to_collection(
    db: &BioVaultDb,
    collection_identifier: &str,
    file_ids: Vec<i64>,
) -> Result<usize> {
    let collection = get_collection(db, collection_identifier)?;


    let mut added = 0;
    for file_id in file_ids {
        // Check if file exists
        let exists: bool = db
            .conn
            .query_row(
                "SELECT COUNT(*) FROM files WHERE id = ?1",
                params![file_id],
                |row| Ok(row.get::<_, i32>(0)? > 0),
            )
            .unwrap_or(false);

        if !exists {
            continue; // Skip non-existent files
        }

        // Insert if not already in collection
        match db.conn.execute(
            "INSERT OR IGNORE INTO collection_files (collection_id, file_id) VALUES (?1, ?2)",
            params![collection.id, file_id],
        ) {
            Ok(1) => added += 1,
            Ok(_) => {} // Already exists, ignore
            Err(e) => return Err(e.into()),
        }
    }

    // Update collection updated_at
    db.conn.execute(
        "UPDATE collections SET updated_at = CURRENT_TIMESTAMP WHERE id = ?1",
        params![collection.id],
    )?;

    Ok(added)
}

/// Remove files from a collection
pub fn remove_files_from_collection(
    db: &BioVaultDb,
    collection_identifier: &str,
    file_ids: Vec<i64>,
) -> Result<usize> {
    let collection = get_collection(db, collection_identifier)?;

    let mut removed = 0;
    for file_id in file_ids {
        match db.conn.execute(
            "DELETE FROM collection_files WHERE collection_id = ?1 AND file_id = ?2",
            params![collection.id, file_id],
        ) {
            Ok(count) => removed += count,
            Err(e) => return Err(e.into()),
        }
    }

    // Update collection updated_at
    db.conn.execute(
        "UPDATE collections SET updated_at = CURRENT_TIMESTAMP WHERE id = ?1",
        params![collection.id],
    )?;

    Ok(removed)
}

/// Delete a collection (files are not deleted, only the collection and its associations)
pub fn delete_collection(db: &BioVaultDb, collection_identifier: &str) -> Result<()> {
    let collection = get_collection(db, collection_identifier)?;

    db.conn.execute(
        "DELETE FROM collections WHERE id = ?1",
        params![collection.id],
    )?;

    Ok(())
}

/// Update collection metadata
pub fn update_collection(
    db: &BioVaultDb,
    collection_identifier: &str,
    name: Option<String>,
    description: Option<Option<String>>, // Option<Option> to allow setting to NULL
    variable_name: Option<String>,
) -> Result<CollectionRecord> {
    let collection = get_collection(db, collection_identifier)?;

    let new_name = name.unwrap_or(collection.name);
    let new_description = description.unwrap_or(collection.description);
    let new_var_name = if let Some(vn) = variable_name {
        validate_variable_name(&vn)?;
        // Check if new variable_name conflicts with existing (excluding current collection)
        let exists: bool = db
            .conn
            .query_row(
                "SELECT COUNT(*) FROM collections WHERE variable_name = ?1 AND id != ?2",
                params![vn, collection.id],
                |row| Ok(row.get::<_, i32>(0)? > 0),
            )
            .unwrap_or(false);

        if exists {
            bail!("Collection with variable_name '{}' already exists", vn);
        }
        vn
    } else {
        collection.variable_name
    };

    db.conn.execute(
        "UPDATE collections SET name = ?1, description = ?2, variable_name = ?3, updated_at = CURRENT_TIMESTAMP WHERE id = ?4",
        params![new_name, new_description, new_var_name, collection.id],
    )?;

    get_collection_by_id(db, collection.id)
}

