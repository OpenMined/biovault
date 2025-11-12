use anyhow::{bail, Result};
use rusqlite::params;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

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

    if !var_name
        .chars()
        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
    {
        bail!("Variable name must be snake_case (lowercase letters, digits, underscores only)");
    }

    if var_name.starts_with('_') || var_name.ends_with('_') {
        bail!("Variable name cannot start or end with underscore");
    }

    if var_name.contains("__") {
        bail!("Variable name cannot contain consecutive underscores");
    }

    if var_name
        .chars()
        .next()
        .map(|c| c.is_ascii_digit())
        .unwrap_or(false)
    {
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
        bail!(
            "Collection with variable_name '{}' already exists",
            var_name
        );
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

    // Get file count for this collection
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

        // Get file count for this collection
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
    for (
        id,
        file_path,
        file_hash,
        file_type,
        file_size,
        data_type,
        _metadata,
        status,
        processing_error,
        created_at,
        updated_at,
        participant_id,
        participant_inferred_sex,
    ) in file_rows
    {
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

/// Add files to one or more collections
/// Returns a map of collection identifier -> number of files added
pub fn add_files_to_collections(
    db: &BioVaultDb,
    collection_identifiers: Vec<&str>,
    file_ids: Vec<i64>,
) -> Result<HashMap<String, usize>> {
    if collection_identifiers.is_empty() {
        return Ok(HashMap::new());
    }

    // Resolve all collections first
    let mut collections = Vec::new();
    for identifier in &collection_identifiers {
        let collection = get_collection(db, identifier)?;
        collections.push((identifier.to_string(), collection.id));
    }

    // Validate all files exist (do this once for efficiency)
    let mut valid_file_ids = Vec::new();
    for file_id in &file_ids {
        let exists: bool = db
            .conn
            .query_row(
                "SELECT COUNT(*) FROM files WHERE id = ?1",
                params![file_id],
                |row| Ok(row.get::<_, i32>(0)? > 0),
            )
            .unwrap_or(false);

        if exists {
            valid_file_ids.push(*file_id);
        }
    }

    // Add files to each collection
    let mut results = HashMap::new();
    for (identifier, collection_id) in &collections {
        let mut added = 0;
        for file_id in &valid_file_ids {
            // Insert if not already in collection
            match db.conn.execute(
                "INSERT OR IGNORE INTO collection_files (collection_id, file_id) VALUES (?1, ?2)",
                params![collection_id, file_id],
            ) {
                Ok(1) => added += 1,
                Ok(_) => {} // Already exists, ignore
                Err(e) => return Err(e.into()),
            }
        }

        // Update collection updated_at
        db.conn.execute(
            "UPDATE collections SET updated_at = CURRENT_TIMESTAMP WHERE id = ?1",
            params![collection_id],
        )?;

        results.insert(identifier.clone(), added);
    }

    Ok(results)
}

/// Add files to a single collection (convenience wrapper)
pub fn add_files_to_collection(
    db: &BioVaultDb,
    collection_identifier: &str,
    file_ids: Vec<i64>,
) -> Result<usize> {
    let results = add_files_to_collections(db, vec![collection_identifier], file_ids)?;
    Ok(results.get(collection_identifier).copied().unwrap_or(0))
}

/// Remove files from one or more collections
/// Returns a map of collection identifier -> number of files removed
pub fn remove_files_from_collections(
    db: &BioVaultDb,
    collection_identifiers: Vec<&str>,
    file_ids: Vec<i64>,
) -> Result<HashMap<String, usize>> {
    if collection_identifiers.is_empty() {
        return Ok(HashMap::new());
    }

    // Resolve all collections first
    let mut collections = Vec::new();
    for identifier in &collection_identifiers {
        let collection = get_collection(db, identifier)?;
        collections.push((identifier.to_string(), collection.id));
    }

    // Remove files from each collection
    let mut results = HashMap::new();
    for (identifier, collection_id) in &collections {
        let mut removed = 0;
        for file_id in &file_ids {
            match db.conn.execute(
                "DELETE FROM collection_files WHERE collection_id = ?1 AND file_id = ?2",
                params![collection_id, file_id],
            ) {
                Ok(count) => removed += count,
                Err(e) => return Err(e.into()),
            }
        }

        // Update collection updated_at
        db.conn.execute(
            "UPDATE collections SET updated_at = CURRENT_TIMESTAMP WHERE id = ?1",
            params![collection_id],
        )?;

        results.insert(identifier.clone(), removed);
    }

    Ok(results)
}

/// Remove files from a single collection (convenience wrapper)
pub fn remove_files_from_collection(
    db: &BioVaultDb,
    collection_identifier: &str,
    file_ids: Vec<i64>,
) -> Result<usize> {
    let results = remove_files_from_collections(db, vec![collection_identifier], file_ids)?;
    Ok(results.get(collection_identifier).copied().unwrap_or(0))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use tempfile::TempDir;

    fn setup_test_db() -> (TempDir, BioVaultDb) {
        let temp = TempDir::new().unwrap();
        config::set_test_biovault_home(temp.path());
        let db = BioVaultDb::new().unwrap();
        (temp, db)
    }

    fn cleanup_test_db(temp: TempDir) {
        drop(temp);
        config::clear_test_biovault_home();
    }

    fn create_test_file(db: &BioVaultDb, path: &str) -> i64 {
        use rusqlite::params;
        db.conn
            .execute(
                "INSERT INTO files (file_path, file_hash) VALUES (?1, ?2)",
                params![path, "test_hash"],
            )
            .unwrap();
        db.conn.last_insert_rowid()
    }

    #[test]
    fn test_generate_variable_name() {
        assert_eq!(generate_variable_name("My Collection"), "my_collection");
        assert_eq!(generate_variable_name("Test-Collection"), "test_collection");
        assert_eq!(generate_variable_name("test_collection"), "test_collection");
        assert_eq!(generate_variable_name("  spaced  name  "), "spaced_name");
        assert_eq!(generate_variable_name("Collection123"), "collection123");
    }

    #[test]
    fn test_validate_variable_name() {
        // Valid names
        assert!(validate_variable_name("valid_name").is_ok());
        assert!(validate_variable_name("valid_name_123").is_ok());
        assert!(validate_variable_name("a").is_ok());

        // Invalid names
        assert!(validate_variable_name("").is_err());
        assert!(validate_variable_name("Invalid-Name").is_err());
        assert!(validate_variable_name("invalid name").is_err());
        assert!(validate_variable_name("_invalid").is_err());
        assert!(validate_variable_name("invalid_").is_err());
        assert!(validate_variable_name("invalid__name").is_err());
        assert!(validate_variable_name("123invalid").is_err());
    }

    #[test]
    fn test_create_collection() {
        let (temp, db) = setup_test_db();

        // Create collection without variable name (auto-generated)
        let collection = create_collection(
            &db,
            "Test Collection".to_string(),
            Some("Test description".to_string()),
            None,
        )
        .unwrap();

        assert_eq!(collection.name, "Test Collection");
        assert_eq!(collection.description, Some("Test description".to_string()));
        assert_eq!(collection.variable_name, "test_collection");
        assert_eq!(collection.file_count, 0);

        // Create collection with custom variable name
        let collection2 = create_collection(
            &db,
            "Another Collection".to_string(),
            None,
            Some("custom_var".to_string()),
        )
        .unwrap();

        assert_eq!(collection2.variable_name, "custom_var");

        cleanup_test_db(temp);
    }

    #[test]
    fn test_create_collection_duplicate_variable_name() {
        let (temp, db) = setup_test_db();

        create_collection(
            &db,
            "First".to_string(),
            None,
            Some("duplicate".to_string()),
        )
        .unwrap();

        // Try to create another with same variable name
        let result = create_collection(
            &db,
            "Second".to_string(),
            None,
            Some("duplicate".to_string()),
        );

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already exists"));

        cleanup_test_db(temp);
    }

    #[test]
    fn test_list_collections() {
        let (temp, db) = setup_test_db();

        // Initially empty
        let collections = list_collections(&db).unwrap();
        assert_eq!(collections.len(), 0);

        // Create some collections
        create_collection(&db, "Collection A".to_string(), None, None).unwrap();
        create_collection(&db, "Collection B".to_string(), None, None).unwrap();

        let collections = list_collections(&db).unwrap();
        assert_eq!(collections.len(), 2);
        assert!(collections.iter().any(|c| c.name == "Collection A"));
        assert!(collections.iter().any(|c| c.name == "Collection B"));

        cleanup_test_db(temp);
    }

    #[test]
    fn test_get_collection_by_id() {
        let (temp, db) = setup_test_db();

        let created = create_collection(
            &db,
            "Test Collection".to_string(),
            Some("Description".to_string()),
            None,
        )
        .unwrap();

        let retrieved = get_collection_by_id(&db, created.id).unwrap();
        assert_eq!(retrieved.id, created.id);
        assert_eq!(retrieved.name, "Test Collection");
        assert_eq!(retrieved.description, Some("Description".to_string()));

        cleanup_test_db(temp);
    }

    #[test]
    fn test_get_collection_by_identifier() {
        let (temp, db) = setup_test_db();

        let created = create_collection(
            &db,
            "Test Name".to_string(),
            None,
            Some("test_var".to_string()),
        )
        .unwrap();

        // Get by variable name
        let by_var = get_collection(&db, "test_var").unwrap();
        assert_eq!(by_var.id, created.id);

        // Get by name
        let by_name = get_collection(&db, "Test Name").unwrap();
        assert_eq!(by_name.id, created.id);

        // Non-existent
        assert!(get_collection(&db, "nonexistent").is_err());

        cleanup_test_db(temp);
    }

    #[test]
    fn test_get_collection_detail_empty() {
        let (temp, db) = setup_test_db();

        let collection =
            create_collection(&db, "Empty Collection".to_string(), None, None).unwrap();

        let detail = get_collection_detail(&db, &collection.variable_name).unwrap();
        assert_eq!(detail.id, collection.id);
        assert_eq!(detail.files.len(), 0);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_add_files_to_collection() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(&db, "Test Collection".to_string(), None, None).unwrap();

        // Create test files
        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        let file2_id = create_test_file(&db, "/path/to/file2.txt");
        let file3_id = create_test_file(&db, "/path/to/file3.txt");

        // Add files to collection
        let added = add_files_to_collection(
            &db,
            &collection.variable_name,
            vec![file1_id, file2_id, file3_id],
        )
        .unwrap();

        assert_eq!(added, 3);

        // Verify file count updated
        let updated = get_collection_by_id(&db, collection.id).unwrap();
        assert_eq!(updated.file_count, 3);

        // Add same files again (should not add duplicates)
        let added_again =
            add_files_to_collection(&db, &collection.variable_name, vec![file1_id, file2_id])
                .unwrap();

        assert_eq!(added_again, 0); // Already in collection

        cleanup_test_db(temp);
    }

    #[test]
    fn test_add_files_to_collections_multiple() {
        let (temp, db) = setup_test_db();

        let coll1 = create_collection(&db, "Collection 1".to_string(), None, None).unwrap();
        let coll2 = create_collection(&db, "Collection 2".to_string(), None, None).unwrap();

        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        let file2_id = create_test_file(&db, "/path/to/file2.txt");

        // Add files to multiple collections
        let results = add_files_to_collections(
            &db,
            vec![&coll1.variable_name, &coll2.variable_name],
            vec![file1_id, file2_id],
        )
        .unwrap();

        assert_eq!(results.len(), 2);
        assert_eq!(results.get(&coll1.variable_name), Some(&2));
        assert_eq!(results.get(&coll2.variable_name), Some(&2));

        // Verify both collections have files
        let updated1 = get_collection_by_id(&db, coll1.id).unwrap();
        let updated2 = get_collection_by_id(&db, coll2.id).unwrap();
        assert_eq!(updated1.file_count, 2);
        assert_eq!(updated2.file_count, 2);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_remove_files_from_collection() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(&db, "Test Collection".to_string(), None, None).unwrap();

        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        let file2_id = create_test_file(&db, "/path/to/file2.txt");
        let file3_id = create_test_file(&db, "/path/to/file3.txt");

        // Add files first
        add_files_to_collection(
            &db,
            &collection.variable_name,
            vec![file1_id, file2_id, file3_id],
        )
        .unwrap();

        // Remove some files
        let removed =
            remove_files_from_collection(&db, &collection.variable_name, vec![file1_id, file2_id])
                .unwrap();

        assert_eq!(removed, 2);

        // Verify file count updated
        let updated = get_collection_by_id(&db, collection.id).unwrap();
        assert_eq!(updated.file_count, 1);

        // Verify detail shows correct files
        let detail = get_collection_detail(&db, &collection.variable_name).unwrap();
        assert_eq!(detail.files.len(), 1);
        assert_eq!(detail.files[0].id, file3_id);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_remove_files_from_collections_multiple() {
        let (temp, db) = setup_test_db();

        let coll1 = create_collection(&db, "Collection 1".to_string(), None, None).unwrap();
        let coll2 = create_collection(&db, "Collection 2".to_string(), None, None).unwrap();

        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        let file2_id = create_test_file(&db, "/path/to/file2.txt");

        // Add files to both collections
        add_files_to_collections(
            &db,
            vec![&coll1.variable_name, &coll2.variable_name],
            vec![file1_id, file2_id],
        )
        .unwrap();

        // Remove from both collections
        let results = remove_files_from_collections(
            &db,
            vec![&coll1.variable_name, &coll2.variable_name],
            vec![file1_id],
        )
        .unwrap();

        assert_eq!(results.len(), 2);
        assert_eq!(results.get(&coll1.variable_name), Some(&1));
        assert_eq!(results.get(&coll2.variable_name), Some(&1));

        cleanup_test_db(temp);
    }

    #[test]
    fn test_update_collection() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(
            &db,
            "Original Name".to_string(),
            Some("Original description".to_string()),
            Some("original_var".to_string()),
        )
        .unwrap();

        // Update name only
        let updated = update_collection(
            &db,
            &collection.variable_name,
            Some("New Name".to_string()),
            None,
            None,
        )
        .unwrap();

        assert_eq!(updated.name, "New Name");
        assert_eq!(
            updated.description,
            Some("Original description".to_string())
        );
        assert_eq!(updated.variable_name, "original_var");

        // Update description (set to None)
        let updated =
            update_collection(&db, &updated.variable_name, None, Some(None), None).unwrap();

        assert_eq!(updated.description, None);

        // Update variable name
        let updated = update_collection(
            &db,
            &updated.variable_name,
            None,
            None,
            Some("new_var_name".to_string()),
        )
        .unwrap();

        assert_eq!(updated.variable_name, "new_var_name");

        cleanup_test_db(temp);
    }

    #[test]
    fn test_update_collection_duplicate_variable_name() {
        let (temp, db) = setup_test_db();

        let _coll1 = create_collection(
            &db,
            "Collection 1".to_string(),
            None,
            Some("var1".to_string()),
        )
        .unwrap();
        let coll2 = create_collection(
            &db,
            "Collection 2".to_string(),
            None,
            Some("var2".to_string()),
        )
        .unwrap();

        // Try to update coll2 to use coll1's variable name
        let result = update_collection(
            &db,
            &coll2.variable_name,
            None,
            None,
            Some("var1".to_string()),
        );

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already exists"));

        cleanup_test_db(temp);
    }

    #[test]
    fn test_delete_collection() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(&db, "To Delete".to_string(), None, None).unwrap();

        // Add some files
        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        add_files_to_collection(&db, &collection.variable_name, vec![file1_id]).unwrap();

        // Delete collection
        delete_collection(&db, &collection.variable_name).unwrap();

        // Verify collection is gone
        assert!(get_collection(&db, &collection.variable_name).is_err());

        // Verify files still exist (not deleted)
        let file_exists: bool = db
            .conn
            .query_row(
                "SELECT COUNT(*) FROM files WHERE id = ?1",
                params![file1_id],
                |row| Ok(row.get::<_, i32>(0)? > 0),
            )
            .unwrap_or(false);
        assert!(file_exists);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_add_files_nonexistent_file() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(&db, "Test Collection".to_string(), None, None).unwrap();

        // Try to add non-existent file ID
        let added = add_files_to_collection(
            &db,
            &collection.variable_name,
            vec![99999], // Non-existent file ID
        )
        .unwrap();

        // Should not error, but should not add anything
        assert_eq!(added, 0);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_add_files_to_collections_empty() {
        let (temp, db) = setup_test_db();

        // Empty collection identifiers
        let results = add_files_to_collections(&db, vec![], vec![1, 2, 3]).unwrap();
        assert_eq!(results.len(), 0);

        cleanup_test_db(temp);
    }

    #[test]
    fn test_get_collection_detail_with_files() {
        let (temp, db) = setup_test_db();

        let collection = create_collection(&db, "Test Collection".to_string(), None, None).unwrap();

        let file1_id = create_test_file(&db, "/path/to/file1.txt");
        let file2_id = create_test_file(&db, "/path/to/file2.txt");

        add_files_to_collection(&db, &collection.variable_name, vec![file1_id, file2_id]).unwrap();

        let detail = get_collection_detail(&db, &collection.variable_name).unwrap();
        assert_eq!(detail.files.len(), 2);
        assert!(detail.files.iter().any(|f| f.id == file1_id));
        assert!(detail.files.iter().any(|f| f.id == file2_id));

        cleanup_test_db(temp);
    }
}
