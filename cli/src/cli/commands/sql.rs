use anyhow::Result;
use colored::Colorize;
use rusqlite::Connection;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::data::BioVaultDb;

/// List all tables in the database
pub async fn tables() -> Result<()> {
    let db = BioVaultDb::new()?;
    let table_list = list_tables(&db.conn)?;

    for table in table_list {
        println!("{}", table);
    }

    Ok(())
}

/// Show table structure/schema
pub async fn structure(table: Option<String>, json: bool) -> Result<()> {
    let db = BioVaultDb::new()?;

    if let Some(table_name) = table {
        let structure = get_table_structure(&db.conn, &table_name)?;
        if json {
            println!("{}", serde_json::to_string_pretty(&structure)?);
        } else {
            print_table_structure(&table_name, &structure);
        }
    } else {
        let tables = list_tables(&db.conn)?;
        if json {
            let mut all_structures = HashMap::new();
            for table_name in &tables {
                all_structures.insert(
                    table_name.clone(),
                    get_table_structure(&db.conn, table_name)?,
                );
            }
            println!("{}", serde_json::to_string_pretty(&all_structures)?);
        } else {
            for table_name in &tables {
                let structure = get_table_structure(&db.conn, table_name)?;
                print_table_structure(table_name, &structure);
                println!();
            }
        }
    }

    Ok(())
}

/// Execute SQL query
pub async fn run(
    query: String,
    output: Option<String>,
    format: String,
    allow_write: bool,
    allow_ddl: bool,
) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Detect operation type
    let op_type = detect_sql_operation(&query);

    // Check permissions
    match op_type {
        SqlOperation::Read => {}
        SqlOperation::Write => {
            if !allow_write {
                anyhow::bail!("Write operations require --allow-write flag");
            }
        }
        SqlOperation::Ddl => {
            if !allow_ddl {
                anyhow::bail!("DDL operations require --allow-ddl flag");
            }
        }
        SqlOperation::Dangerous => {
            anyhow::bail!("Dangerous operation detected. This operation is not allowed.");
        }
    }

    // Check for SQL injection patterns
    if has_sql_injection_risk(&query) {
        anyhow::bail!("Potential SQL injection detected. Please review your query.");
    }

    // Execute query
    let results = execute_query(&db.conn, &query)?;

    // Output results
    if let Some(output_path) = output {
        export_results(&results, &output_path, &format)?;
        println!(
            "{}",
            format!("✓ Exported {} rows to {}", results.rows.len(), output_path).green()
        );
    } else {
        print_results(&results, 10);
    }

    Ok(())
}

/// Parameters for CSV/TSV import
pub struct ImportParams {
    pub file_path: String,
    pub table: String,
    pub participant_col: Option<String>,
    pub map: Option<String>,
    pub format: String,
    pub on_mismatch: String,
    pub dry_run: bool,
    pub allow_overwrite: bool,
}

/// Import CSV/TSV data into database
pub async fn import(params: ImportParams) -> Result<()> {
    let ImportParams {
        file_path,
        table,
        participant_col,
        map,
        format,
        on_mismatch,
        dry_run,
        allow_overwrite,
    } = params;
    let db = BioVaultDb::new()?;

    // Check if the base table name (before user_ prefix) is protected
    if is_protected_table(&table) {
        anyhow::bail!("Cannot import into protected table: {}", table);
    }

    // Apply user_ prefix
    let table_name = apply_user_prefix(&table);

    // Check if table exists
    let table_exists = table_exists(&db.conn, &table_name)?;
    if table_exists && !allow_overwrite {
        anyhow::bail!(
            "Table '{}' already exists. Use --allow-overwrite to replace it.",
            table_name
        );
    }

    // Parse column mapping
    let column_mapping = if let Some(map_str) = map {
        parse_column_mapping(&map_str)?
    } else {
        HashMap::new()
    };

    // Parse CSV/TSV
    let csv_data = parse_csv_file(&file_path, &format)?;

    // Determine participant column name (before or after mapping)
    let participant_col_name = participant_col
        .clone()
        .unwrap_or_else(|| "participant_id".to_string());

    // Check if participant column exists in original headers
    let original_participant_col = if csv_data.headers.contains(&participant_col_name) {
        participant_col_name.clone()
    } else {
        // If not found, it might be the mapped name - find the original name
        let reverse_mapping: HashMap<String, String> = column_mapping
            .iter()
            .map(|(k, v)| (v.clone(), k.clone()))
            .collect();

        reverse_mapping
            .get(&participant_col_name)
            .cloned()
            .unwrap_or(participant_col_name.clone())
    };

    // Apply column mapping
    let (mapped_headers, mapped_rows) = apply_column_mapping(csv_data, &column_mapping)?;

    // Find participant column in mapped headers
    let final_participant_col = column_mapping
        .get(&original_participant_col)
        .cloned()
        .unwrap_or(original_participant_col);
    let participant_col_idx = mapped_headers
        .iter()
        .position(|h| h == &final_participant_col)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "Participant column '{}' not found in CSV headers after mapping. Headers: {:?}",
                final_participant_col,
                mapped_headers
            )
        })?;

    // Extract participant IDs
    let participant_ids: Vec<String> = mapped_rows
        .iter()
        .map(|row| row[participant_col_idx].clone())
        .collect();

    // Validate participant IDs
    let validation = validate_participant_ids(&db.conn, &participant_ids)?;

    // Print validation report
    print_validation_report(&validation, &on_mismatch);

    // Handle mismatches
    let mismatch_action = parse_mismatch_action(&on_mismatch)?;
    match mismatch_action {
        MismatchAction::Error => {
            if !validation.invalid.is_empty() {
                anyhow::bail!("Invalid participant IDs found. Use --on-mismatch skip or create.");
            }
        }
        MismatchAction::Skip => {
            if !validation.invalid.is_empty() {
                println!(
                    "{}",
                    format!(
                        "Skipping {} rows with invalid participant IDs",
                        validation.invalid.len()
                    )
                    .yellow()
                );
            }
        }
        MismatchAction::Create => {
            if !validation.invalid.is_empty() {
                println!(
                    "{}",
                    format!("Will create {} new participants", validation.invalid.len()).yellow()
                );
            }
        }
    }

    if dry_run {
        println!("{}", "✓ Dry run complete - no changes made".green());
        return Ok(());
    }

    // Execute import
    let rows_imported = execute_import(ExecuteImportParams {
        conn: &db.conn,
        table: &table_name,
        headers: &mapped_headers,
        rows: &mapped_rows,
        participant_col_idx,
        validation: &validation,
        mismatch_action,
        allow_overwrite,
    })?;

    println!(
        "{}",
        format!("✓ Imported {} rows into {}", rows_imported, table_name)
            .green()
            .bold()
    );

    Ok(())
}

/// Drop a user table
pub async fn drop(table: String, confirm: bool) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Ensure it's a user table
    if !table.starts_with("user_") {
        anyhow::bail!("Can only drop tables with 'user_' prefix. Table: {}", table);
    }

    // Check if table exists
    if !table_exists(&db.conn, &table)? {
        anyhow::bail!("Table '{}' does not exist", table);
    }

    // Get row count
    let row_count: i64 =
        db.conn
            .query_row(&format!("SELECT COUNT(*) FROM {}", table), [], |row| {
                row.get(0)
            })?;

    // Confirm deletion
    if !confirm {
        println!("Table: {}", table.cyan());
        println!("Rows: {}", row_count);
        println!("Type: User table");
        println!();
        println!(
            "{}",
            "WARNING: This will permanently delete the table and all data."
                .yellow()
                .bold()
        );
        print!("Type the table name to confirm: ");
        std::io::stdout().flush()?;

        let mut input = String::new();
        std::io::stdin().read_line(&mut input)?;
        let input = input.trim();

        if input != table {
            anyhow::bail!("Table name did not match. Aborted.");
        }
    }

    // Drop table
    db.conn.execute(&format!("DROP TABLE {}", table), [])?;

    println!("{}", format!("✓ Table {} dropped", table).green().bold());

    Ok(())
}

// ============================================================================
// Helper Functions
// ============================================================================

#[derive(Debug)]
enum SqlOperation {
    Read,
    Write,
    Ddl,
    Dangerous,
}

fn detect_sql_operation(query: &str) -> SqlOperation {
    let query_upper = query.trim().to_uppercase();

    // Check for dangerous operations first
    if query_upper.contains("DROP DATABASE")
        || query_upper.contains("DROP SCHEMA")
        || query_upper.contains("SHUTDOWN")
        || query_upper.contains("ATTACH DATABASE")
    {
        return SqlOperation::Dangerous;
    }

    // Check for DDL
    if query_upper.starts_with("CREATE")
        || query_upper.starts_with("ALTER")
        || query_upper.starts_with("DROP")
        || query_upper.starts_with("TRUNCATE")
    {
        return SqlOperation::Ddl;
    }

    // Check for write operations
    if query_upper.starts_with("INSERT")
        || query_upper.starts_with("UPDATE")
        || query_upper.starts_with("DELETE")
        || query_upper.starts_with("REPLACE")
    {
        return SqlOperation::Write;
    }

    SqlOperation::Read
}

fn has_sql_injection_risk(query: &str) -> bool {
    let patterns = [
        ";",  // Multiple statements
        "--", // SQL comments
        "/*", // Block comments
        "*/",
        "UNION",  // Union-based injection
        "OR 1=1", // Always true conditions
        "OR '1'='1",
        "' OR '",
        "\" OR \"",
    ];

    let query_upper = query.to_uppercase();
    for pattern in &patterns {
        if query_upper.contains(pattern) {
            return true;
        }
    }

    false
}

fn list_tables(conn: &Connection) -> Result<Vec<String>> {
    let mut stmt = conn.prepare(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
    )?;

    let tables = stmt
        .query_map([], |row| row.get(0))?
        .collect::<Result<Vec<String>, _>>()?;

    Ok(tables)
}

#[derive(Debug, serde::Serialize)]
struct ColumnInfo {
    name: String,
    type_name: String,
    nullable: bool,
    default_value: Option<String>,
    primary_key: bool,
}

#[derive(Debug, serde::Serialize)]
struct TableStructure {
    columns: Vec<ColumnInfo>,
    indexes: Vec<String>,
    foreign_keys: Vec<String>,
}

fn get_table_structure(conn: &Connection, table: &str) -> Result<TableStructure> {
    let mut columns = Vec::new();

    // Get column info
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table))?;
    let column_iter = stmt.query_map([], |row| {
        Ok(ColumnInfo {
            name: row.get(1)?,
            type_name: row.get(2)?,
            nullable: row.get::<_, i32>(3)? == 0,
            default_value: row.get(4)?,
            primary_key: row.get::<_, i32>(5)? == 1,
        })
    })?;

    for col in column_iter {
        columns.push(col?);
    }

    // Get indexes
    let mut indexes = Vec::new();
    let mut stmt = conn.prepare(&format!("PRAGMA index_list({})", table))?;
    let index_iter = stmt.query_map([], |row| row.get::<_, String>(1))?;

    for idx in index_iter {
        indexes.push(idx?);
    }

    // Get foreign keys
    let mut foreign_keys = Vec::new();
    let mut stmt = conn.prepare(&format!("PRAGMA foreign_key_list({})", table))?;
    let fk_iter = stmt.query_map([], |row| {
        let table: String = row.get(2)?;
        let from: String = row.get(3)?;
        let to: String = row.get(4)?;
        Ok(format!("{} -> {}({})", from, table, to))
    })?;

    for fk in fk_iter {
        foreign_keys.push(fk?);
    }

    Ok(TableStructure {
        columns,
        indexes,
        foreign_keys,
    })
}

fn print_table_structure(table_name: &str, structure: &TableStructure) {
    println!("{} {}", "TABLE:".bold(), table_name.cyan());

    for col in &structure.columns {
        let mut flags = Vec::new();
        if col.primary_key {
            flags.push("PRIMARY KEY".to_string());
        }
        if !col.nullable {
            flags.push("NOT NULL".to_string());
        }
        if let Some(default) = &col.default_value {
            flags.push(format!("DEFAULT {}", default));
        }

        let flag_str = if flags.is_empty() {
            String::new()
        } else {
            format!(" ({})", flags.join(", "))
        };

        println!(
            "  ├─ {}: {}{}",
            col.name.green(),
            col.type_name.yellow(),
            flag_str.dimmed()
        );
    }

    if !structure.foreign_keys.is_empty() {
        for fk in &structure.foreign_keys {
            println!("  └─ FK: {}", fk.blue());
        }
    }

    if !structure.indexes.is_empty() {
        println!("  Indexes: {}", structure.indexes.join(", ").dimmed());
    }
}

#[derive(Debug)]
struct QueryResults {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn execute_query(conn: &Connection, query: &str) -> Result<QueryResults> {
    let mut stmt = conn.prepare(query)?;
    let column_count = stmt.column_count();

    let headers: Vec<String> = stmt.column_names().iter().map(|s| s.to_string()).collect();

    let rows_iter = stmt.query_map([], |row| {
        let mut values = Vec::new();
        for i in 0..column_count {
            let value: Result<String, _> = row.get(i);
            values.push(value.unwrap_or_else(|_| "NULL".to_string()));
        }
        Ok(values)
    })?;

    let mut rows = Vec::new();
    for row in rows_iter {
        rows.push(row?);
    }

    Ok(QueryResults { headers, rows })
}

fn print_results(results: &QueryResults, limit: usize) {
    if results.rows.is_empty() {
        println!("{}", "No results".yellow());
        return;
    }

    // Print headers
    println!("{}", results.headers.join(" | ").bold());

    // Print rows (limited)
    let display_count = results.rows.len().min(limit);
    for row in &results.rows[..display_count] {
        println!("{}", row.join(" | "));
    }

    if results.rows.len() > limit {
        println!(
            "{}",
            format!("... and {} more rows", results.rows.len() - limit).dimmed()
        );
        println!("{}", "Use --output to export all results".yellow());
    }

    println!();
    println!("{} rows", results.rows.len());
}

fn export_results(results: &QueryResults, path: &str, format: &str) -> Result<()> {
    let mut file = File::create(path)?;

    match format {
        "json" => {
            let json_rows: Vec<HashMap<String, String>> = results
                .rows
                .iter()
                .map(|row| {
                    results
                        .headers
                        .iter()
                        .zip(row.iter())
                        .map(|(h, v)| (h.clone(), v.clone()))
                        .collect()
                })
                .collect();
            writeln!(file, "{}", serde_json::to_string_pretty(&json_rows)?)?;
        }
        "tsv" => {
            writeln!(file, "{}", results.headers.join("\t"))?;
            for row in &results.rows {
                writeln!(file, "{}", row.join("\t"))?;
            }
        }
        _ => {
            // Default to CSV
            let mut wtr = csv::Writer::from_writer(file);
            wtr.write_record(&results.headers)?;
            for row in &results.rows {
                wtr.write_record(row)?;
            }
            wtr.flush()?;
        }
    }

    Ok(())
}

fn apply_user_prefix(table: &str) -> String {
    if table.starts_with("user_") {
        table.to_string()
    } else {
        format!("user_{}", table)
    }
}

fn is_protected_table(table: &str) -> bool {
    let protected = [
        "participants",
        "files",
        "projects",
        "runs",
        "run_participants",
        "messages",
        "schema_version",
        "genotype_metadata",
        "dev_envs",
    ];

    protected.contains(&table)
}

fn table_exists(conn: &Connection, table: &str) -> Result<bool> {
    let count: i64 = conn.query_row(
        "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?",
        [table],
        |row| row.get(0),
    )?;

    Ok(count > 0)
}

fn parse_column_mapping(map_str: &str) -> Result<HashMap<String, String>> {
    let mut mapping = HashMap::new();

    for pair in map_str.split(',') {
        let parts: Vec<&str> = pair.split(':').collect();
        if parts.len() != 2 {
            anyhow::bail!("Invalid mapping format: '{}'. Expected 'old:new'", pair);
        }
        mapping.insert(parts[0].trim().to_string(), parts[1].trim().to_string());
    }

    Ok(mapping)
}

#[derive(Debug)]
struct CsvData {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn parse_csv_file(path: &str, format: &str) -> Result<CsvData> {
    let mut reader = if format == "tsv" {
        csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?
    } else {
        csv::Reader::from_path(path)?
    };

    let headers = reader.headers()?.iter().map(|s| s.to_string()).collect();

    let mut rows = Vec::new();
    for result in reader.records() {
        let record = result?;
        let row: Vec<String> = record.iter().map(|s| s.to_string()).collect();
        rows.push(row);
    }

    Ok(CsvData { headers, rows })
}

fn apply_column_mapping(
    data: CsvData,
    mapping: &HashMap<String, String>,
) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let mapped_headers: Vec<String> = data
        .headers
        .iter()
        .map(|h| mapping.get(h).unwrap_or(h).clone())
        .collect();

    Ok((mapped_headers, data.rows))
}

#[derive(Debug)]
struct ParticipantValidation {
    valid: Vec<String>,
    invalid: Vec<String>,
}

fn validate_participant_ids(conn: &Connection, ids: &[String]) -> Result<ParticipantValidation> {
    let mut valid = Vec::new();
    let mut invalid = Vec::new();

    for id in ids {
        let count: i64 = conn.query_row(
            "SELECT COUNT(*) FROM participants WHERE participant_id = ?",
            [id],
            |row| row.get(0),
        )?;

        if count > 0 {
            if !valid.contains(id) {
                valid.push(id.clone());
            }
        } else if !invalid.contains(id) {
            invalid.push(id.clone());
        }
    }

    Ok(ParticipantValidation { valid, invalid })
}

fn print_validation_report(validation: &ParticipantValidation, on_mismatch: &str) {
    println!("{}", "Validation Report:".bold());
    println!(
        "  {} Valid participant IDs: {}",
        "✓".green(),
        validation.valid.len()
    );

    if !validation.invalid.is_empty() {
        println!(
            "  {} Invalid participant IDs: {} {}",
            "✗".red(),
            validation.invalid.len(),
            format!("[{}]", validation.invalid.join(", ")).dimmed()
        );
        println!();
        println!("  Action: --on-mismatch={}", on_mismatch.yellow());
    }
    println!();
}

#[derive(Debug, Copy, Clone)]
enum MismatchAction {
    Error,
    Skip,
    Create,
}

fn parse_mismatch_action(action: &str) -> Result<MismatchAction> {
    match action {
        "error" => Ok(MismatchAction::Error),
        "skip" => Ok(MismatchAction::Skip),
        "create" => Ok(MismatchAction::Create),
        _ => anyhow::bail!(
            "Invalid on-mismatch action: {}. Use error, skip, or create.",
            action
        ),
    }
}

struct ExecuteImportParams<'a> {
    conn: &'a Connection,
    table: &'a str,
    headers: &'a [String],
    rows: &'a [Vec<String>],
    participant_col_idx: usize,
    validation: &'a ParticipantValidation,
    mismatch_action: MismatchAction,
    allow_overwrite: bool,
}

fn execute_import(params: ExecuteImportParams) -> Result<usize> {
    let ExecuteImportParams {
        conn,
        table,
        headers,
        rows,
        participant_col_idx,
        validation,
        mismatch_action,
        allow_overwrite,
    } = params;
    // Start transaction
    conn.execute("BEGIN TRANSACTION", [])?;

    // Drop existing table if overwrite is allowed
    if allow_overwrite && table_exists(conn, table)? {
        conn.execute(&format!("DROP TABLE {}", table), [])?;
    }

    // Create table
    let column_defs: Vec<String> = headers.iter().map(|h| format!("{} TEXT", h)).collect();

    let create_sql = format!(
        "CREATE TABLE IF NOT EXISTS {} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            {},
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )",
        table,
        column_defs.join(", ")
    );

    conn.execute(&create_sql, [])?;

    // Add foreign key constraint for participant_id
    // Note: SQLite doesn't support adding FK after table creation, so we'd need to recreate the table
    // For now, we'll skip this and just ensure referential integrity in the import logic

    // Create participants if needed
    if matches!(mismatch_action, MismatchAction::Create) {
        for id in &validation.invalid {
            conn.execute(
                "INSERT OR IGNORE INTO participants (participant_id) VALUES (?)",
                [id],
            )?;
        }
    }

    // Insert rows
    let mut imported_count = 0;
    let placeholders = headers.iter().map(|_| "?").collect::<Vec<_>>().join(", ");
    let insert_sql = format!(
        "INSERT INTO {} ({}) VALUES ({})",
        table,
        headers.join(", "),
        placeholders
    );

    for row in rows {
        let participant_id = &row[participant_col_idx];

        // Skip invalid participants if action is Skip
        if matches!(mismatch_action, MismatchAction::Skip)
            && validation.invalid.contains(participant_id)
        {
            continue;
        }

        // Convert row to rusqlite params
        let params: Vec<&dyn rusqlite::ToSql> =
            row.iter().map(|v| v as &dyn rusqlite::ToSql).collect();
        conn.execute(&insert_sql, params.as_slice())?;
        imported_count += 1;
    }

    // Commit transaction
    conn.execute("COMMIT", [])?;

    Ok(imported_count)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rusqlite::Connection;
    use tempfile::TempDir;

    fn setup_test_db() -> (TempDir, Connection) {
        let temp = TempDir::new().unwrap();
        let db_path = temp.path().join("test.db");
        let conn = Connection::open(&db_path).unwrap();

        // Initialize basic schema
        conn.execute(
            "CREATE TABLE participants (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                participant_id TEXT UNIQUE NOT NULL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )
        .unwrap();

        (temp, conn)
    }

    #[test]
    fn test_detect_sql_operation() {
        assert!(matches!(
            detect_sql_operation("SELECT * FROM users"),
            SqlOperation::Read
        ));
        assert!(matches!(
            detect_sql_operation("select id from test"),
            SqlOperation::Read
        ));
        assert!(matches!(
            detect_sql_operation("INSERT INTO users VALUES (1)"),
            SqlOperation::Write
        ));
        assert!(matches!(
            detect_sql_operation("UPDATE users SET name='test'"),
            SqlOperation::Write
        ));
        assert!(matches!(
            detect_sql_operation("DELETE FROM users"),
            SqlOperation::Write
        ));
        assert!(matches!(
            detect_sql_operation("CREATE TABLE test (id INT)"),
            SqlOperation::Ddl
        ));
        assert!(matches!(
            detect_sql_operation("DROP TABLE users"),
            SqlOperation::Ddl
        ));
        assert!(matches!(
            detect_sql_operation("ALTER TABLE users ADD COLUMN test TEXT"),
            SqlOperation::Ddl
        ));
        assert!(matches!(
            detect_sql_operation("DROP DATABASE test"),
            SqlOperation::Dangerous
        ));
    }

    #[test]
    fn test_has_sql_injection_risk() {
        assert!(has_sql_injection_risk(
            "SELECT * FROM users; DROP TABLE users"
        ));
        assert!(has_sql_injection_risk(
            "SELECT * FROM users WHERE id = '1' OR '1'='1'"
        ));
        assert!(has_sql_injection_risk("SELECT * FROM users -- comment"));
        assert!(has_sql_injection_risk("SELECT * FROM users /* comment */"));
        assert!(has_sql_injection_risk(
            "SELECT * FROM users UNION SELECT * FROM passwords"
        ));
        assert!(!has_sql_injection_risk("SELECT * FROM users WHERE id = ?"));
        assert!(!has_sql_injection_risk("SELECT name, email FROM users"));
    }

    #[test]
    fn test_list_tables() {
        let (_temp, conn) = setup_test_db();

        conn.execute("CREATE TABLE test_table (id INTEGER PRIMARY KEY)", [])
            .unwrap();

        let tables = list_tables(&conn).unwrap();

        assert!(tables.contains(&"participants".to_string()));
        assert!(tables.contains(&"test_table".to_string()));
        assert!(!tables.iter().any(|t| t.starts_with("sqlite_")));
    }

    #[test]
    fn test_get_table_structure() {
        let (_temp, conn) = setup_test_db();

        let structure = get_table_structure(&conn, "participants").unwrap();

        assert_eq!(structure.columns.len(), 3);
        assert_eq!(structure.columns[0].name, "id");
        assert!(structure.columns[0].primary_key);
        assert_eq!(structure.columns[1].name, "participant_id");
        assert!(!structure.columns[1].nullable);
    }

    #[test]
    fn test_apply_user_prefix() {
        assert_eq!(apply_user_prefix("test"), "user_test");
        assert_eq!(apply_user_prefix("user_test"), "user_test");
        assert_eq!(apply_user_prefix("measurements"), "user_measurements");
    }

    #[test]
    fn test_is_protected_table() {
        assert!(is_protected_table("participants"));
        assert!(is_protected_table("files"));
        assert!(is_protected_table("messages"));
        assert!(!is_protected_table("user_custom"));
        assert!(!is_protected_table("user_measurements"));
    }

    #[test]
    fn test_table_exists() {
        let (_temp, conn) = setup_test_db();

        assert!(table_exists(&conn, "participants").unwrap());
        assert!(!table_exists(&conn, "nonexistent").unwrap());
    }

    #[test]
    fn test_parse_column_mapping() {
        let mapping = parse_column_mapping("id:participant_id,val:measurement").unwrap();

        assert_eq!(mapping.get("id"), Some(&"participant_id".to_string()));
        assert_eq!(mapping.get("val"), Some(&"measurement".to_string()));

        // Test invalid format
        assert!(parse_column_mapping("invalid_format").is_err());
        assert!(parse_column_mapping("a:b:c").is_err());
    }

    #[test]
    fn test_parse_csv_file() {
        use std::io::Write;

        let temp = TempDir::new().unwrap();
        let csv_path = temp.path().join("test.csv");

        let mut file = std::fs::File::create(&csv_path).unwrap();
        writeln!(file, "participant_id,value").unwrap();
        writeln!(file, "P001,100").unwrap();
        writeln!(file, "P002,200").unwrap();
        file.flush().unwrap();

        let data = parse_csv_file(csv_path.to_str().unwrap(), "csv").unwrap();

        assert_eq!(data.headers, vec!["participant_id", "value"]);
        assert_eq!(data.rows.len(), 2);
        assert_eq!(data.rows[0], vec!["P001", "100"]);
        assert_eq!(data.rows[1], vec!["P002", "200"]);
    }

    #[test]
    fn test_apply_column_mapping() {
        let data = CsvData {
            headers: vec!["subject_id".to_string(), "val".to_string()],
            rows: vec![
                vec!["P001".to_string(), "100".to_string()],
                vec!["P002".to_string(), "200".to_string()],
            ],
        };

        let mut mapping = HashMap::new();
        mapping.insert("subject_id".to_string(), "participant_id".to_string());
        mapping.insert("val".to_string(), "measurement".to_string());

        let (headers, rows) = apply_column_mapping(data, &mapping).unwrap();

        assert_eq!(headers, vec!["participant_id", "measurement"]);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0], vec!["P001", "100"]);
    }

    #[test]
    fn test_validate_participant_ids() {
        let (_temp, conn) = setup_test_db();

        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();
        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P002')",
            [],
        )
        .unwrap();

        let ids = vec!["P001".to_string(), "P002".to_string(), "P003".to_string()];
        let validation = validate_participant_ids(&conn, &ids).unwrap();

        assert_eq!(validation.valid.len(), 2);
        assert_eq!(validation.invalid.len(), 1);
        assert!(validation.valid.contains(&"P001".to_string()));
        assert!(validation.valid.contains(&"P002".to_string()));
        assert!(validation.invalid.contains(&"P003".to_string()));
    }

    #[test]
    fn test_parse_mismatch_action() {
        assert!(matches!(
            parse_mismatch_action("error").unwrap(),
            MismatchAction::Error
        ));
        assert!(matches!(
            parse_mismatch_action("skip").unwrap(),
            MismatchAction::Skip
        ));
        assert!(matches!(
            parse_mismatch_action("create").unwrap(),
            MismatchAction::Create
        ));
        assert!(parse_mismatch_action("invalid").is_err());
    }

    #[test]
    fn test_execute_query() {
        let (_temp, conn) = setup_test_db();

        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();

        let results = execute_query(&conn, "SELECT * FROM participants").unwrap();

        assert_eq!(results.headers.len(), 3);
        assert_eq!(results.rows.len(), 1);
        assert!(results.headers.contains(&"participant_id".to_string()));
    }

    #[test]
    fn test_export_results_csv() {
        use std::fs;

        let temp = TempDir::new().unwrap();
        let export_path = temp.path().join("export.csv");

        let results = QueryResults {
            headers: vec!["id".to_string(), "name".to_string()],
            rows: vec![
                vec!["1".to_string(), "Alice".to_string()],
                vec!["2".to_string(), "Bob".to_string()],
            ],
        };

        export_results(&results, export_path.to_str().unwrap(), "csv").unwrap();

        let content = fs::read_to_string(&export_path).unwrap();
        assert!(content.contains("id,name"));
        assert!(content.contains("Alice"));
        assert!(content.contains("Bob"));
    }

    #[test]
    fn test_export_results_json() {
        use std::fs;

        let temp = TempDir::new().unwrap();
        let export_path = temp.path().join("export.json");

        let results = QueryResults {
            headers: vec!["id".to_string(), "name".to_string()],
            rows: vec![vec!["1".to_string(), "Alice".to_string()]],
        };

        export_results(&results, export_path.to_str().unwrap(), "json").unwrap();

        let content = fs::read_to_string(&export_path).unwrap();
        assert!(content.contains("\"id\""));
        assert!(content.contains("\"name\""));
        assert!(content.contains("Alice"));
    }

    #[test]
    fn test_execute_import() {
        let (_temp, conn) = setup_test_db();

        // Add participants
        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();
        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P002')",
            [],
        )
        .unwrap();

        let headers = vec!["participant_id".to_string(), "measurement".to_string()];
        let rows = vec![
            vec!["P001".to_string(), "100".to_string()],
            vec!["P002".to_string(), "200".to_string()],
        ];

        let validation = ParticipantValidation {
            valid: vec!["P001".to_string(), "P002".to_string()],
            invalid: vec![],
        };

        let count = execute_import(ExecuteImportParams {
            conn: &conn,
            table: "user_test",
            headers: &headers,
            rows: &rows,
            participant_col_idx: 0,
            validation: &validation,
            mismatch_action: MismatchAction::Error,
            allow_overwrite: false,
        })
        .unwrap();

        assert_eq!(count, 2);

        // Verify data was inserted
        let result: i64 = conn
            .query_row("SELECT COUNT(*) FROM user_test", [], |row| row.get(0))
            .unwrap();
        assert_eq!(result, 2);
    }

    #[test]
    fn test_execute_import_with_skip() {
        let (_temp, conn) = setup_test_db();

        // Add only one participant
        conn.execute(
            "INSERT INTO participants (participant_id) VALUES ('P001')",
            [],
        )
        .unwrap();

        let headers = vec!["participant_id".to_string(), "measurement".to_string()];
        let rows = vec![
            vec!["P001".to_string(), "100".to_string()],
            vec!["P999".to_string(), "200".to_string()], // Invalid participant
        ];

        let validation = ParticipantValidation {
            valid: vec!["P001".to_string()],
            invalid: vec!["P999".to_string()],
        };

        let count = execute_import(ExecuteImportParams {
            conn: &conn,
            table: "user_test",
            headers: &headers,
            rows: &rows,
            participant_col_idx: 0,
            validation: &validation,
            mismatch_action: MismatchAction::Skip,
            allow_overwrite: false,
        })
        .unwrap();

        // Should only import 1 row (P001), skip P999
        assert_eq!(count, 1);
    }

    #[test]
    fn test_execute_import_with_create() {
        let (_temp, conn) = setup_test_db();

        let headers = vec!["participant_id".to_string(), "measurement".to_string()];
        let rows = vec![
            vec!["P001".to_string(), "100".to_string()],
            vec!["P002".to_string(), "200".to_string()],
        ];

        let validation = ParticipantValidation {
            valid: vec![],
            invalid: vec!["P001".to_string(), "P002".to_string()],
        };

        let count = execute_import(ExecuteImportParams {
            conn: &conn,
            table: "user_test",
            headers: &headers,
            rows: &rows,
            participant_col_idx: 0,
            validation: &validation,
            mismatch_action: MismatchAction::Create,
            allow_overwrite: false,
        })
        .unwrap();

        assert_eq!(count, 2);

        // Verify participants were created
        let result: i64 = conn
            .query_row("SELECT COUNT(*) FROM participants", [], |row| row.get(0))
            .unwrap();
        assert_eq!(result, 2);
    }
}
