use serial_test::serial;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

// Helper to set up isolated test environment
fn setup_test_env() -> TempDir {
    let temp = TempDir::new().unwrap();
    let home = temp.path().join(".biovault");
    fs::create_dir_all(&home).unwrap();
    biovault::config::set_test_biovault_home(&home);
    temp
}

// Helper to tear down test environment
fn teardown_test_env() {
    biovault::config::clear_test_biovault_home();
}

/// Integration test: Full import → query → export → drop cycle
#[test]
#[serial]
fn test_full_import_export_cycle() {
    // Setup test environment
    let temp = setup_test_env();

    // Initialize database
    let db = biovault::data::BioVaultDb::new().unwrap();

    // Create test CSV file
    let csv_path = temp.path().join("test_data.csv");
    let mut csv_file = fs::File::create(&csv_path).unwrap();
    writeln!(csv_file, "participant_id,measurement,test_date").unwrap();
    writeln!(csv_file, "P001,100,2024-01-01").unwrap();
    writeln!(csv_file, "P002,200,2024-01-02").unwrap();
    writeln!(csv_file, "P003,150,2024-01-03").unwrap();
    csv_file.flush().unwrap();

    // Step 1: Import CSV (creates participants and table)
    let runtime = tokio::runtime::Runtime::new().unwrap();
    runtime
        .block_on(async {
            biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
                file_path: csv_path.to_str().unwrap().to_string(),
                table: "test_measurements".to_string(),
                participant_col: None,
                map: None,
                format: "csv".to_string(),
                on_mismatch: "create".to_string(),
                dry_run: false,
                allow_overwrite: false,
            })
            .await
        })
        .unwrap();

    // Step 2: Verify table exists
    let tables: Vec<String> = db
        .conn
        .prepare(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='user_test_measurements'",
        )
        .unwrap()
        .query_map([], |row| row.get(0))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(tables.len(), 1);
    assert_eq!(tables[0], "user_test_measurements");

    // Step 3: Query the data
    let export_path = temp.path().join("export.csv");
    runtime
        .block_on(async {
            biovault::cli::commands::sql::run(
                "SELECT * FROM user_test_measurements ORDER BY participant_id".to_string(),
                Some(export_path.to_str().unwrap().to_string()),
                "csv".to_string(),
                false,
                false,
            )
            .await
        })
        .unwrap();

    // Step 4: Verify exported content
    let exported = fs::read_to_string(&export_path).unwrap();
    assert!(exported.contains("participant_id"));
    assert!(exported.contains("P001"));
    assert!(exported.contains("P002"));
    assert!(exported.contains("P003"));
    assert!(exported.contains("100"));
    assert!(exported.contains("200"));
    assert!(exported.contains("150"));

    // Count rows (header + 3 data rows = 4 lines)
    let line_count = exported.lines().count();
    assert_eq!(line_count, 4);

    // Step 5: Drop the table
    runtime
        .block_on(async {
            biovault::cli::commands::sql::drop("user_test_measurements".to_string(), true).await
        })
        .unwrap();

    // Step 6: Verify table was dropped
    let tables_after: Vec<String> = db
        .conn
        .prepare(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='user_test_measurements'",
        )
        .unwrap()
        .query_map([], |row| row.get(0))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(tables_after.len(), 0);

    teardown_test_env();
}

/// Integration test: Import with column mapping
#[test]
#[serial]
fn test_import_with_column_mapping() {
    let temp = setup_test_env();

    let _db = biovault::data::BioVaultDb::new().unwrap();

    // Create CSV with non-standard column names
    let csv_path = temp.path().join("mapped_data.csv");
    let mut csv_file = fs::File::create(&csv_path).unwrap();
    writeln!(csv_file, "subject_id,value,test_date").unwrap();
    writeln!(csv_file, "P001,100,2024-01-01").unwrap();
    writeln!(csv_file, "P002,200,2024-01-02").unwrap();
    csv_file.flush().unwrap();

    // Import with column mapping
    let runtime = tokio::runtime::Runtime::new().unwrap();
    runtime
        .block_on(async {
            biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
                file_path: csv_path.to_str().unwrap().to_string(),
                table: "mapped_test".to_string(),
                participant_col: Some("subject_id".to_string()),
                map: Some("subject_id:participant_id,value:measurement".to_string()),
                format: "csv".to_string(),
                on_mismatch: "create".to_string(),
                dry_run: false,
                allow_overwrite: false,
            })
            .await
        })
        .unwrap();

    // Verify table structure has mapped column names
    runtime
        .block_on(async {
            biovault::cli::commands::sql::structure(Some("user_mapped_test".to_string()), false)
                .await
        })
        .unwrap();

    // Query to verify data is accessible with new column names
    let export_path = temp.path().join("mapped_export.csv");
    runtime
        .block_on(async {
            biovault::cli::commands::sql::run(
                "SELECT participant_id, measurement FROM user_mapped_test".to_string(),
                Some(export_path.to_str().unwrap().to_string()),
                "csv".to_string(),
                false,
                false,
            )
            .await
        })
        .unwrap();

    let exported = fs::read_to_string(&export_path).unwrap();
    assert!(exported.contains("participant_id"));
    assert!(exported.contains("measurement"));
    assert!(!exported.contains("subject_id"));
    assert!(!exported.contains("value"));

    // Cleanup
    runtime
        .block_on(async {
            biovault::cli::commands::sql::drop("user_mapped_test".to_string(), true).await
        })
        .unwrap();

    teardown_test_env();
}

/// Integration test: Participant validation with different strategies
#[test]
#[serial]
fn test_participant_validation_strategies() {
    let temp = setup_test_env();

    let db = biovault::data::BioVaultDb::new().unwrap();

    // Create one existing participant
    db.conn
        .execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();

    // Create CSV with mix of valid and invalid participants
    let csv_path = temp.path().join("mixed_data.csv");
    let mut csv_file = fs::File::create(&csv_path).unwrap();
    writeln!(csv_file, "participant_id,measurement").unwrap();
    writeln!(csv_file, "P001,100").unwrap();
    writeln!(csv_file, "P999,200").unwrap();
    csv_file.flush().unwrap();

    let runtime = tokio::runtime::Runtime::new().unwrap();

    // Test 1: Error strategy (should fail)
    let result = runtime.block_on(async {
        biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
            file_path: csv_path.to_str().unwrap().to_string(),
            table: "error_test".to_string(),
            participant_col: None,
            map: None,
            format: "csv".to_string(),
            on_mismatch: "error".to_string(),
            dry_run: false,
            allow_overwrite: false,
        })
        .await
    });
    assert!(result.is_err());

    // Test 2: Skip strategy (should import only P001)
    runtime
        .block_on(async {
            biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
                file_path: csv_path.to_str().unwrap().to_string(),
                table: "skip_test".to_string(),
                participant_col: None,
                map: None,
                format: "csv".to_string(),
                on_mismatch: "skip".to_string(),
                dry_run: false,
                allow_overwrite: false,
            })
            .await
        })
        .unwrap();

    let count: i64 = db
        .conn
        .query_row("SELECT COUNT(*) FROM user_skip_test", [], |row| row.get(0))
        .unwrap();
    assert_eq!(count, 1);

    // Test 3: Create strategy (should create P999 and import both)
    runtime
        .block_on(async {
            biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
                file_path: csv_path.to_str().unwrap().to_string(),
                table: "create_test".to_string(),
                participant_col: None,
                map: None,
                format: "csv".to_string(),
                on_mismatch: "create".to_string(),
                dry_run: false,
                allow_overwrite: false,
            })
            .await
        })
        .unwrap();

    let count: i64 = db
        .conn
        .query_row("SELECT COUNT(*) FROM user_create_test", [], |row| {
            row.get(0)
        })
        .unwrap();
    assert_eq!(count, 2);

    // Verify P999 was created
    let participant_count: i64 = db
        .conn
        .query_row(
            "SELECT COUNT(*) FROM participants WHERE participant_id='P999'",
            [],
            |row| row.get(0),
        )
        .unwrap();
    assert_eq!(participant_count, 1);

    // Cleanup
    runtime.block_on(async {
        let _ = biovault::cli::commands::sql::drop("user_skip_test".to_string(), true).await;
        let _ = biovault::cli::commands::sql::drop("user_create_test".to_string(), true).await;
    });

    teardown_test_env();
}

/// Integration test: Export formats (CSV, TSV, JSON)
#[test]
#[serial]
fn test_export_formats() {
    let temp = setup_test_env();

    let db = biovault::data::BioVaultDb::new().unwrap();

    // Create test data
    db.conn
        .execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();

    db.conn
        .execute(
            "CREATE TABLE user_format_test (
                id INTEGER PRIMARY KEY,
                participant_id TEXT,
                value INTEGER
            )",
            [],
        )
        .unwrap();

    db.conn
        .execute(
            "INSERT INTO user_format_test (participant_id, value) VALUES ('P001', 100)",
            [],
        )
        .unwrap();

    let runtime = tokio::runtime::Runtime::new().unwrap();

    // Test CSV export
    let csv_path = temp.path().join("export.csv");
    runtime
        .block_on(async {
            biovault::cli::commands::sql::run(
                "SELECT * FROM user_format_test".to_string(),
                Some(csv_path.to_str().unwrap().to_string()),
                "csv".to_string(),
                false,
                false,
            )
            .await
        })
        .unwrap();

    let csv_content = fs::read_to_string(&csv_path).unwrap();
    assert!(csv_content.contains(","));
    assert!(csv_content.contains("P001"));

    // Test TSV export
    let tsv_path = temp.path().join("export.tsv");
    runtime
        .block_on(async {
            biovault::cli::commands::sql::run(
                "SELECT * FROM user_format_test".to_string(),
                Some(tsv_path.to_str().unwrap().to_string()),
                "tsv".to_string(),
                false,
                false,
            )
            .await
        })
        .unwrap();

    let tsv_content = fs::read_to_string(&tsv_path).unwrap();
    assert!(tsv_content.contains("\t"));
    assert!(tsv_content.contains("P001"));

    // Test JSON export
    let json_path = temp.path().join("export.json");
    runtime
        .block_on(async {
            biovault::cli::commands::sql::run(
                "SELECT * FROM user_format_test".to_string(),
                Some(json_path.to_str().unwrap().to_string()),
                "json".to_string(),
                false,
                false,
            )
            .await
        })
        .unwrap();

    let json_content = fs::read_to_string(&json_path).unwrap();
    assert!(json_content.contains("\"participant_id\""));
    assert!(json_content.contains("\"value\""));
    assert!(json_content.contains("P001"));

    // Cleanup
    runtime.block_on(async {
        let _ = biovault::cli::commands::sql::drop("user_format_test".to_string(), true).await;
    });

    teardown_test_env();
}

/// Integration test: Protected table validation
#[test]
#[serial]
fn test_protected_tables() {
    let temp = setup_test_env();

    let _db = biovault::data::BioVaultDb::new().unwrap();

    // Create CSV
    let csv_path = temp.path().join("protected_test.csv");
    let mut csv_file = fs::File::create(&csv_path).unwrap();
    writeln!(csv_file, "participant_id,value").unwrap();
    writeln!(csv_file, "P001,100").unwrap();
    csv_file.flush().unwrap();

    let runtime = tokio::runtime::Runtime::new().unwrap();

    // Try to import into protected table (should fail)
    let result = runtime.block_on(async {
        biovault::cli::commands::sql::import(biovault::cli::commands::sql::ImportParams {
            file_path: csv_path.to_str().unwrap().to_string(),
            table: "participants".to_string(), // Protected table
            participant_col: None,
            map: None,
            format: "csv".to_string(),
            on_mismatch: "create".to_string(),
            dry_run: false,
            allow_overwrite: true, // Even with allow_overwrite
        })
        .await
    });

    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("protected"));

    teardown_test_env();
}
