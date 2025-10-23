# BioVault SQL CLI Commands - Design Doc

## Overview
CLI commands for safe SQL operations on BioVault database with participant validation and AI-friendly output.

## Command Structure

### 1. List Tables
```bash
bv sql tables [list]
```
**Output**: Simple list of table names, one per line
**Use case**: Quick reference, AI consumption

---

### 2. Table Structure

#### Single Table
```bash
bv sql structure <table_name>
```

#### All Tables (default)
```bash
bv sql structure [all]
```

**Output Format**:
```
TABLE: participants
├─ id: INTEGER PRIMARY KEY
├─ participant_id: TEXT UNIQUE NOT NULL
├─ created_at: TIMESTAMP
└─ INDEX: idx_participant_id

TABLE: measurements
├─ id: INTEGER PRIMARY KEY
├─ participant_id: INTEGER FOREIGN KEY -> participants(id)
...
```

**Also support JSON flag**:
```bash
bv sql structure --json
```

---

### 3. Run Arbitrary SQL

```bash
bv sql run "<query>" [OPTIONS]
```

**Options**:
- `--output <path>` or `-o <path>`: Export results to CSV
- `--format <csv|tsv|json>`: Output format (default: csv)
- `--no-sanitize`: Skip injection check (dangerous, requires confirmation)

**Features**:
- Read-only by default (SELECT only)
- `--allow-write` flag for INSERT/UPDATE/DELETE (requires confirmation)
- SQL injection detection
- Prepared statement execution
- Result preview in terminal (first 10 rows)
- Full export to file

**Examples**:
```bash
# Query and display
bv sql run "SELECT * FROM participants LIMIT 10"

# Export to CSV
bv sql run "SELECT * FROM measurements WHERE date > '2024-01-01'" -o results.csv

# Export as TSV
bv sql run "SELECT * FROM participants" -o data.tsv --format tsv

# Write operation (with confirmation)
bv sql run "UPDATE participants SET status='active'" --allow-write
```

---

### 4. Import Data

```bash
bv sql import <file_path> [OPTIONS]
```

**Options**:
- `--table <name>`: Target table name (required, auto-prefixed with `user_`)
- `--participant-col <name>`: Column containing participant_id (default: "participant_id")
- `--format <csv|tsv>`: Auto-detect by extension if not specified
- `--map <old:new,...>`: Rename columns during import (e.g., "subject_id:participant_id,val:measurement")
- `--on-mismatch <error|skip|create>`: How to handle participant_id mismatches
  - `error`: Stop import (default)
  - `skip`: Skip rows with invalid participant_id
  - `create`: Create new participant records
- `--dry-run`: Validate without importing
- `--allow-overwrite`: Allow overwriting existing user table (requires confirmation)

**Table Naming & Safety**:
- All imported tables prefixed with `user_` (e.g., `--table measurements` → `user_measurements`)
- Protected tables (core BioVault schema) cannot be overwritten:
  - `participants`, `measurements`, `biomarkers`, etc.
- Existing `user_*` tables require `--allow-overwrite` flag
- Prevents accidental data loss

**Process**:
1. Parse file & apply column mapping
2. Validate table name (add `user_` prefix, check conflicts)
3. Detect participant_id column (after mapping)
4. Validate all participant_ids exist in database
5. Report validation results
6. Confirm before import (unless `--yes`)
7. Create table with participant_id FK constraint
8. Import data with transaction rollback on error

**Examples**:
```bash
# Basic import with validation → creates user_custom_measurements
bv sql import data.csv --table custom_measurements

# Custom participant column + rename
bv sql import survey.tsv --table survey_results --participant-col subject_id \
  --map "subject_id:participant_id,response:answer"

# Create missing participants
bv sql import new_data.csv --table lab_results --on-mismatch create

# Dry run to validate
bv sql import data.csv --table test --dry-run

# Overwrite existing user table
bv sql import updated.csv --table measurements --allow-overwrite
```

**Validation Report**:
```
Parsing: data.csv (150 rows, 5 columns)
Participant column: participant_id
Validation:
  ✓ 145 valid participant_ids
  ✗ 5 invalid participant_ids: [P001, P002, P003, P004, P005]

Action: --on-mismatch=error (default)
Options:
  • --on-mismatch skip: Import 145 rows, skip 5
  • --on-mismatch create: Create 5 new participants, import all 150

Proceed? [y/N]
```

---

## Security Considerations

### SQL Injection Prevention
1. Use prepared statements/parameterized queries
2. Parse and validate SQL AST before execution
3. Whitelist allowed operations by default
4. Flag dangerous operations: DROP, TRUNCATE, ALTER
5. Require explicit `--allow-dangerous` for DDL

### Read/Write Protection
- Default: SELECT only
- `--allow-write`: DML (INSERT, UPDATE, DELETE)
- `--allow-ddl`: DDL (CREATE, ALTER, DROP) - requires double confirmation

### Participant Data Protection
- Validate participant_id references
- Transaction rollback on partial failure
- Backup recommendation before bulk operations
- Audit logging for write operations

---

## Implementation Notes

### Dependencies
- SQL parser (e.g., `sqlparser` crate)
- CSV/TSV parser (`csv` crate)
- Database connection (existing `bv` infrastructure)

### Error Handling
- Clear error messages with suggestions
- Graceful degradation
- Transaction support for atomic operations
- Progress bars for large imports

### Output Formatting
- Terminal: Pretty table with truncation
- CSV/TSV: Standard RFC 4180 format
- JSON: Structured with metadata
- Support for streaming large results

---

## Testing Strategy

### Unit Tests

**SQL Parser & Sanitizer**:
```rust
#[test]
fn test_detect_sql_injection() {
    assert!(is_dangerous("SELECT * FROM users; DROP TABLE users"));
    assert!(is_dangerous("SELECT * FROM users WHERE id = '1' OR '1'='1'"));
    assert!(!is_dangerous("SELECT * FROM users WHERE id = ?"));
}

#[test]
fn test_sql_operation_detection() {
    assert_eq!(detect_operation("SELECT * FROM users"), SqlOp::Read);
    assert_eq!(detect_operation("INSERT INTO users"), SqlOp::Write);
    assert_eq!(detect_operation("DROP TABLE users"), SqlOp::Ddl);
}
```

**CSV/TSV Parser**:
```rust
#[test]
fn test_parse_csv_with_headers() {
    let data = "participant_id,value\nP001,100\nP002,200";
    let result = parse_csv(data).unwrap();
    assert_eq!(result.headers, vec!["participant_id", "value"]);
    assert_eq!(result.rows.len(), 2);
}

#[test]
fn test_detect_participant_column() {
    let headers = vec!["subject_id", "value"];
    let col = find_participant_column(&headers, Some("subject_id"));
    assert_eq!(col, Ok(0));
}
```

**Column Mapping**:
```rust
#[test]
fn test_parse_column_mapping() {
    let mapping = parse_mapping("subject_id:participant_id,val:measurement").unwrap();
    assert_eq!(mapping.get("subject_id"), Some(&"participant_id".to_string()));
    assert_eq!(mapping.get("val"), Some(&"measurement".to_string()));
}

#[test]
fn test_apply_column_mapping() {
    let headers = vec!["subject_id", "val", "date"];
    let mapping = hashmap!{"subject_id" => "participant_id", "val" => "measurement"};
    let result = apply_mapping(headers, &mapping);
    assert_eq!(result, vec!["participant_id", "measurement", "date"]);
}
```

**Table Name Validation**:
```rust
#[test]
fn test_apply_user_prefix() {
    assert_eq!(apply_user_prefix("measurements"), "user_measurements");
    assert_eq!(apply_user_prefix("user_measurements"), "user_measurements");
}

#[test]
fn test_protect_core_tables() {
    assert!(is_protected_table("participants"));
    assert!(is_protected_table("measurements"));
    assert!(!is_protected_table("user_custom"));
}
```

**Participant Validation**:
```rust
#[test]
fn test_validate_participant_ids() {
    let db = setup_test_db();
    insert_participant(&db, "P001");
    insert_participant(&db, "P002");

    let ids = vec!["P001", "P002", "P003"];
    let result = validate_participant_ids(&db, &ids);

    assert_eq!(result.valid, vec!["P001", "P002"]);
    assert_eq!(result.invalid, vec!["P003"]);
}
```

### Integration Tests

**Full Import → Query → Export → Delete Cycle**:
```rust
#[test]
fn test_full_import_export_cycle() {
    let temp_db = setup_test_database();
    let test_csv = create_test_csv();

    // Step 1: Setup participants
    insert_participant(&temp_db, "P001");
    insert_participant(&temp_db, "P002");
    insert_participant(&temp_db, "P003");

    // Step 2: Import CSV
    let import_result = run_import(
        &temp_db,
        test_csv.path(),
        "test_data",
        None, // participant_col
        None, // mapping
        OnMismatch::Error,
    );
    assert!(import_result.is_ok());
    assert_eq!(import_result.unwrap().rows_imported, 3);

    // Step 3: Verify table exists with correct schema
    let tables = list_tables(&temp_db);
    assert!(tables.contains(&"user_test_data".to_string()));

    let structure = get_table_structure(&temp_db, "user_test_data");
    assert!(structure.contains("participant_id"));
    assert!(structure.contains("FOREIGN KEY"));

    // Step 4: Query imported data
    let query_result = run_query(
        &temp_db,
        "SELECT * FROM user_test_data ORDER BY participant_id",
        false, // allow_write
    );
    assert!(query_result.is_ok());
    let rows = query_result.unwrap();
    assert_eq!(rows.len(), 3);

    // Step 5: Export to CSV
    let export_path = temp_dir().join("export.csv");
    let export_result = export_query_results(
        &temp_db,
        "SELECT * FROM user_test_data",
        &export_path,
        ExportFormat::Csv,
    );
    assert!(export_result.is_ok());
    assert!(export_path.exists());

    // Step 6: Verify exported content
    let exported = fs::read_to_string(&export_path).unwrap();
    assert!(exported.contains("participant_id"));
    assert!(exported.contains("P001"));

    // Step 7: Delete table
    let delete_result = run_query(
        &temp_db,
        "DROP TABLE user_test_data",
        true, // allow_ddl
    );
    assert!(delete_result.is_ok());

    let tables_after = list_tables(&temp_db);
    assert!(!tables_after.contains(&"user_test_data".to_string()));
}

#[test]
fn test_import_with_column_mapping() {
    let temp_db = setup_test_database();
    let csv_content = "subject_id,measurement,test_date\nP001,100,2024-01-01\n";
    let test_csv = create_csv_file(csv_content);

    insert_participant(&temp_db, "P001");

    // Import with column mapping
    let mapping = parse_mapping("subject_id:participant_id").unwrap();
    let result = run_import(
        &temp_db,
        test_csv.path(),
        "mapped_data",
        Some("subject_id"),
        Some(mapping),
        OnMismatch::Error,
    );

    assert!(result.is_ok());

    // Verify column was renamed
    let structure = get_table_structure(&temp_db, "user_mapped_data");
    assert!(structure.contains("participant_id"));
    assert!(!structure.contains("subject_id"));
}

#[test]
fn test_import_with_mismatch_strategies() {
    let temp_db = setup_test_database();
    let csv_content = "participant_id,value\nP001,100\nP999,200\n";
    let test_csv = create_csv_file(csv_content);

    insert_participant(&temp_db, "P001");

    // Test error strategy
    let error_result = run_import(
        &temp_db,
        test_csv.path(),
        "error_test",
        None,
        None,
        OnMismatch::Error,
    );
    assert!(error_result.is_err());

    // Test skip strategy
    let skip_result = run_import(
        &temp_db,
        test_csv.path(),
        "skip_test",
        None,
        None,
        OnMismatch::Skip,
    );
    assert!(skip_result.is_ok());
    assert_eq!(skip_result.unwrap().rows_imported, 1);

    // Test create strategy
    let create_result = run_import(
        &temp_db,
        test_csv.path(),
        "create_test",
        None,
        None,
        OnMismatch::Create,
    );
    assert!(create_result.is_ok());
    assert_eq!(create_result.unwrap().rows_imported, 2);

    // Verify new participant was created
    let participants = run_query(
        &temp_db,
        "SELECT participant_id FROM participants WHERE participant_id = 'P999'",
        false,
    );
    assert_eq!(participants.unwrap().len(), 1);
}

#[test]
fn test_protect_core_tables() {
    let temp_db = setup_test_database();
    let csv_content = "participant_id,value\nP001,100\n";
    let test_csv = create_csv_file(csv_content);

    // Attempt to overwrite protected table
    let result = run_import(
        &temp_db,
        test_csv.path(),
        "participants", // Protected table
        None,
        None,
        OnMismatch::Error,
    );

    // Should fail even with allow_overwrite
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("protected"));
}
```

### CLI Integration Tests
```bash
#!/bin/bash
# test_sql_cli.sh

set -e

# Setup test database
TEST_DB=$(mktemp -d)/test.db
bv --db $TEST_DB init

# Add test participants
bv --db $TEST_DB add-participant P001
bv --db $TEST_DB add-participant P002

# Create test CSV
cat > test_data.csv << EOF
participant_id,measurement,date
P001,100,2024-01-01
P002,200,2024-01-02
EOF

# Test 1: Import
bv --db $TEST_DB sql import test_data.csv --table measurements
echo "✓ Import successful"

# Test 2: List tables
TABLES=$(bv --db $TEST_DB sql tables)
echo "$TABLES" | grep -q "user_measurements" || exit 1
echo "✓ Table created"

# Test 3: Query
RESULTS=$(bv --db $TEST_DB sql run "SELECT * FROM user_measurements")
echo "$RESULTS" | grep -q "P001" || exit 1
echo "✓ Query successful"

# Test 4: Export
bv --db $TEST_DB sql run "SELECT * FROM user_measurements" -o export.csv
grep -q "P001,100" export.csv || exit 1
echo "✓ Export successful"

# Test 5: Structure
STRUCTURE=$(bv --db $TEST_DB sql structure user_measurements)
echo "$STRUCTURE" | grep -q "FOREIGN KEY" || exit 1
echo "✓ Structure includes FK constraint"

# Test 6: Delete
bv --db $TEST_DB sql run "DROP TABLE user_measurements" --allow-ddl
TABLES_AFTER=$(bv --db $TEST_DB sql tables)
echo "$TABLES_AFTER" | grep -qv "user_measurements" || exit 1
echo "✓ Table deleted"

echo "All tests passed!"
```

## Open Questions

1. ✅ **Column Mapping**: Implemented with `--map` flag

2. **Batch Size**: For large imports, should we support chunking?
   ```bash
   bv sql import huge.csv --batch-size 1000
   ```

3. **Connection String**: Support for non-default database paths?
   ```bash
   bv sql --db /path/to/db.sqlite run "SELECT ..."
   ```

4. **Export Limits**: Should we limit export sizes or warn on large results?

5. **SQL Dialect**: SQLite-specific features vs standard SQL?

6. **Type Inference**: Auto-detect column types during import or default to TEXT?

### 5. Drop Tables (Cleanup)

```bash
bv sql drop <table_name> [OPTIONS]
```

**Options**:
- `--confirm`: Skip confirmation prompt (dangerous)

**Safety**:
- Only allows dropping `user_*` tables
- Protected tables cannot be dropped
- Requires confirmation with table name input
- Shows row count before deletion

**Example**:
```bash
bv sql drop user_measurements

# Output:
Table: user_measurements
Rows: 150
Type: User table

WARNING: This will permanently delete the table and all data.
Type the table name to confirm: user_measurements
✓ Table dropped

# Shorthand
bv sql drop user_measurements --confirm
```

---

## Future Enhancements

- `bv sql query-builder` - Interactive query builder
- `bv sql backup` - Export entire database
- `bv sql migrate` - Schema migration support
- `bv sql analyze` - Query performance analysis
- Template queries (common operations)
- Saved queries/aliases
- `bv sql run --file query.sql` - Execute SQL from file

---

## Priority Implementation Order

1. ✅ `bv sql tables` - Simple, foundational
2. ✅ `bv sql structure` - AI-friendly schema export
3. ✅ `bv sql run` (read-only) - Core query functionality
4. ⏭️ `bv sql run` (write) - Add safety guards
5. ⏭️ `bv sql import` - Most complex, high value
6. ⏭️ `bv sql drop` - Cleanup utility

---

## Summary

**Core Commands**:
```bash
# Discovery
bv sql tables                    # List all tables
bv sql structure [table]        # Show schema (AI-friendly)

# Query & Export
bv sql run "SELECT ..." [-o file.csv] [--format csv|tsv|json]

# Import with Safety
bv sql import data.csv --table name \
  [--map "old:new"] \
  [--participant-col col] \
  [--on-mismatch error|skip|create]

# Cleanup
bv sql drop user_table_name
```

**Key Safety Features**:
- ✅ `user_` prefix on all imported tables
- ✅ Protected table validation
- ✅ Participant ID validation with FK constraints
- ✅ Column mapping for flexible imports
- ✅ Transaction-based imports with rollback
- ✅ SQL injection prevention
- ✅ Comprehensive unit & integration tests
