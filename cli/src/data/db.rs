use anyhow::{Context, Result};
use rusqlite::{params, Connection};
use std::fs;
use std::path::{Path, PathBuf};
use tracing::info;

pub struct BioVaultDb {
    pub conn: Connection,
}

impl BioVaultDb {
    /// Open or create the BioVault database with automatic migration
    pub fn new() -> Result<Self> {
        let db_path = get_biovault_db_path()?;

        // Check if migration needed
        if needs_migration()? {
            migrate_from_messages_db(&db_path)?;
        }

        let conn = Connection::open(&db_path)
            .with_context(|| format!("Failed to open database at {:?}", db_path))?;

        // Enable WAL mode for better concurrency
        conn.pragma_update(None, "journal_mode", "WAL")?;
        conn.pragma_update(None, "busy_timeout", 5000)?;

        // Initialize/update schema
        Self::init_schema(&conn)?;

        Ok(Self { conn })
    }

    /// Initialize schema from SQL file
    fn init_schema(conn: &Connection) -> Result<()> {
        let schema = include_str!("../schema.sql");

        if let Err(err) = conn.execute_batch(schema) {
            let err_msg = err.to_string();

            // Older databases might be missing newly introduced columns (e.g. runs.pipeline_id).
            // Attempt to migrate the legacy table layout and retry once before bailing out.
            if err_msg.contains("no such column") {
                migrate_runs_table(conn)?;
                conn.execute_batch(schema)?;
            } else {
                return Err(err.into());
            }
        }

        // Run migrations for existing databases
        Self::run_migrations(conn)?;

        // Check/update schema version
        let current_version = get_schema_version(conn)?;
        if current_version.is_none() {
            conn.execute(
                "INSERT INTO schema_version (version) VALUES (?1)",
                ["2.0.0"],
            )?;
            info!("Initialized schema version 2.0.0");
        }

        Ok(())
    }

    /// Run migrations for existing databases
    fn run_migrations(conn: &Connection) -> Result<()> {
        // Check if data_type column exists in files table
        if table_exists(conn, "files")? {
            if !column_exists(conn, "files", "data_type")? {
                info!("Adding data_type column to files table");
                conn.execute(
                    "ALTER TABLE files ADD COLUMN data_type TEXT DEFAULT 'Unknown'",
                    [],
                )?;
                conn.execute(
                    "CREATE INDEX IF NOT EXISTS idx_files_data_type ON files(data_type)",
                    [],
                )?;
                info!("Migration complete: added data_type column and index");
            }

            // Add status column to files if it doesn't exist
            if !column_exists(conn, "files", "status")? {
                info!("Adding status column to files table");
                conn.execute(
                    "ALTER TABLE files ADD COLUMN status TEXT DEFAULT 'complete'",
                    [],
                )?;
                conn.execute(
                    "CREATE INDEX IF NOT EXISTS idx_files_status ON files(status)",
                    [],
                )?;
                info!("Migration complete: added status column and index");
            }

            // Add processing_error column to files if it doesn't exist
            if !column_exists(conn, "files", "processing_error")? {
                info!("Adding processing_error column to files table");
                conn.execute("ALTER TABLE files ADD COLUMN processing_error TEXT", [])?;
                info!("Migration complete: added processing_error column");
            }

            // Add queue_added_at column to files if it doesn't exist
            if !column_exists(conn, "files", "queue_added_at")? {
                info!("Adding queue_added_at column to files table");
                conn.execute("ALTER TABLE files ADD COLUMN queue_added_at DATETIME", [])?;
                info!("Migration complete: added queue_added_at column");
            }
        }

        // Add inferred_sex to participants if it doesn't exist
        if table_exists(conn, "participants")?
            && !column_exists(conn, "participants", "inferred_sex")?
        {
            info!("Adding inferred_sex column to participants table");
            conn.execute("ALTER TABLE participants ADD COLUMN inferred_sex TEXT", [])?;
            info!("Migration complete: added inferred_sex to participants");
        }

        // Drop source and grch_version columns from files table if they exist
        // (these should only be in genotype_metadata table)
        let source_exists = column_exists(conn, "files", "source")?;
        let grch_exists = column_exists(conn, "files", "grch_version")?;

        if table_exists(conn, "files")? && (source_exists || grch_exists) {
            info!("Dropping source and grch_version columns from files table");

            // SQLite doesn't support DROP COLUMN directly, need to recreate table
            conn.execute("BEGIN TRANSACTION", [])?;

            // Create new table without source/grch_version
            conn.execute(
                "CREATE TABLE files_new (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    participant_id INTEGER,
                    file_path TEXT UNIQUE NOT NULL,
                    file_hash TEXT NOT NULL,
                    file_type TEXT,
                    file_size INTEGER,
                    metadata TEXT,
                    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                    data_type TEXT DEFAULT 'Unknown',
                    status TEXT DEFAULT 'complete',
                    processing_error TEXT,
                    queue_added_at DATETIME,
                    FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE SET NULL
                )",
                [],
            )?;

            // Copy data
            conn.execute(
                "INSERT INTO files_new (id, participant_id, file_path, file_hash, file_type, file_size,
                    metadata, created_at, updated_at, data_type, status, processing_error, queue_added_at)
                SELECT id, participant_id, file_path, file_hash, file_type, file_size,
                    metadata, created_at, updated_at, data_type, status, processing_error, queue_added_at
                FROM files",
                [],
            )?;

            // Drop old table
            conn.execute("DROP TABLE files", [])?;

            // Rename new table
            conn.execute("ALTER TABLE files_new RENAME TO files", [])?;

            // Recreate indexes
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_files_participant_id ON files(participant_id)",
                [],
            )?;
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_files_file_type ON files(file_type)",
                [],
            )?;
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_files_hash ON files(file_hash)",
                [],
            )?;
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_files_data_type ON files(data_type)",
                [],
            )?;
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_files_status ON files(status)",
                [],
            )?;

            conn.execute("COMMIT", [])?;

            info!("Migration complete: removed source and grch_version columns from files table");
        }

        // Create genotype_metadata table if it doesn't exist
        conn.execute(
            "CREATE TABLE IF NOT EXISTS genotype_metadata (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                file_id INTEGER UNIQUE NOT NULL,
                source TEXT,
                grch_version TEXT,
                row_count INTEGER,
                chromosome_count INTEGER,
                inferred_sex TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
            )",
            [],
        )?;
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_genotype_file_id ON genotype_metadata(file_id)",
            [],
        )?;

        // Add inferred_sex column to genotype_metadata if it doesn't exist
        if !column_exists(conn, "genotype_metadata", "inferred_sex")? {
            info!("Adding inferred_sex column to genotype_metadata table");
            conn.execute(
                "ALTER TABLE genotype_metadata ADD COLUMN inferred_sex TEXT",
                [],
            )?;
            info!("Migration complete: added inferred_sex column");
        }

        // Add jupyter_port column to dev_envs if it doesn't exist
        if table_exists(conn, "dev_envs")? && !column_exists(conn, "dev_envs", "jupyter_port")? {
            info!("Adding jupyter_port column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_port INTEGER", [])?;
            info!("Migration complete: added jupyter_port column");
        }

        // Add jupyter_pid column to dev_envs if it doesn't exist
        if table_exists(conn, "dev_envs")? && !column_exists(conn, "dev_envs", "jupyter_pid")? {
            info!("Adding jupyter_pid column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_pid INTEGER", [])?;
            info!("Migration complete: added jupyter_pid column");
        }

        // Add jupyter_url column to dev_envs if it doesn't exist
        if table_exists(conn, "dev_envs")? && !column_exists(conn, "dev_envs", "jupyter_url")? {
            info!("Adding jupyter_url column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_url TEXT", [])?;
            info!("Migration complete: added jupyter_url column");
        }

        // Add jupyter_token column to dev_envs if it doesn't exist
        if table_exists(conn, "dev_envs")? && !column_exists(conn, "dev_envs", "jupyter_token")? {
            info!("Adding jupyter_token column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_token TEXT", [])?;
            info!("Migration complete: added jupyter_token column");
        }

        // Add metadata column to runs if it doesn't exist
        migrate_runs_table(conn)?;

        Ok(())
    }

    /// Get the underlying connection (for advanced use)
    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}

fn get_biovault_db_path() -> Result<PathBuf> {
    Ok(crate::config::get_biovault_home()?.join("biovault.db"))
}

fn get_messages_db_path() -> Result<PathBuf> {
    Ok(crate::config::get_biovault_home()?.join("messages.db"))
}

fn needs_migration() -> Result<bool> {
    let biovault_db = get_biovault_db_path()?;
    let messages_db = get_messages_db_path()?;

    // Migrate if messages.db exists but biovault.db doesn't
    Ok(messages_db.exists() && !biovault_db.exists())
}

fn migrate_from_messages_db(target_path: &Path) -> Result<()> {
    let messages_db = get_messages_db_path()?;

    println!("🔄 Migrating messages.db → biovault.db...");
    info!("Starting database migration");

    // 1. Backup old DB
    let backup_path = messages_db.with_extension("db.backup");
    fs::copy(&messages_db, &backup_path).context("Failed to backup messages.db")?;
    println!("✓ Backed up to messages.db.backup");
    info!("Backed up messages.db to {:?}", backup_path);

    // 2. Copy to new location
    fs::copy(&messages_db, target_path).context("Failed to copy to biovault.db")?;
    println!("✓ Copied to biovault.db");
    info!("Copied to biovault.db at {:?}", target_path);

    // 3. Open and upgrade schema
    let conn = Connection::open(target_path)?;

    // Add schema version table
    conn.execute(
        "CREATE TABLE IF NOT EXISTS schema_version (
            version TEXT PRIMARY KEY,
            applied_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )",
        [],
    )?;

    // Add new tables (messages table already exists from copy)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS participants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            participant_id TEXT UNIQUE NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE IF NOT EXISTS files (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            participant_id INTEGER,
            file_path TEXT UNIQUE NOT NULL,
            file_hash TEXT NOT NULL,
            file_type TEXT,
            file_size INTEGER,
            metadata TEXT,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE SET NULL
        )",
        [],
    )?;

    // Add indexes
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_participant_id ON participants(participant_id)",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_files_participant_id ON files(participant_id)",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_files_file_type ON files(file_type)",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_files_hash ON files(file_hash)",
        [],
    )?;

    // Desktop tables (optional - can be added later when desktop migrates)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS projects (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE NOT NULL,
            author TEXT NOT NULL,
            workflow TEXT NOT NULL,
            template TEXT NOT NULL,
            project_path TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE IF NOT EXISTS runs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            project_id INTEGER NOT NULL,
            work_dir TEXT NOT NULL,
            participant_count INTEGER NOT NULL,
            status TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (project_id) REFERENCES projects(id) ON DELETE CASCADE
        )",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_runs_project_id ON runs(project_id)",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status)",
        [],
    )?;

    conn.execute(
        "CREATE TABLE IF NOT EXISTS run_participants (
            run_id INTEGER NOT NULL,
            participant_id INTEGER NOT NULL,
            FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE,
            FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE CASCADE,
            PRIMARY KEY (run_id, participant_id)
        )",
        [],
    )?;

    // Note: Pipeline tables are now in schema.sql (base schema)

    // Mark schema version
    conn.execute(
        "INSERT INTO schema_version (version) VALUES (?1)",
        ["2.0.0"],
    )?;

    println!("✓ Schema upgraded to v2.0.0");
    println!("✅ Migration complete!");
    info!("Migration complete - schema version 2.0.0");

    Ok(())
}

fn migrate_runs_table(conn: &Connection) -> Result<()> {
    if !table_exists(conn, "runs")? {
        return Ok(());
    }

    add_column_if_missing(conn, "runs", "pipeline_id", "INTEGER")?;
    add_column_if_missing(conn, "runs", "step_id", "INTEGER")?;
    add_column_if_missing(conn, "runs", "results_dir", "TEXT")?;
    add_column_if_missing(conn, "runs", "participant_count", "INTEGER")?;
    add_column_if_missing(conn, "runs", "metadata", "TEXT")?;
    add_column_if_missing(conn, "runs", "completed_at", "DATETIME")?;

    if column_exists(conn, "runs", "pipeline_id")? {
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_runs_pipeline_id ON runs(pipeline_id)",
            [],
        )?;
    }

    if column_exists(conn, "runs", "step_id")? {
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_runs_step_id ON runs(step_id)",
            [],
        )?;
    }

    if column_exists(conn, "runs", "status")? {
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status)",
            [],
        )?;
    }

    Ok(())
}

fn table_exists(conn: &Connection, table: &str) -> Result<bool> {
    conn.query_row(
        "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?1",
        params![table],
        |row| row.get::<_, i32>(0),
    )
    .map(|count| count > 0)
    .map_err(Into::into)
}

fn column_exists(conn: &Connection, table: &str, column: &str) -> Result<bool> {
    // Table names are trusted (internal constants) but still escape single quotes defensively.
    let escaped_table = table.replace('\'', "''");
    let sql = format!(
        "SELECT COUNT(*) FROM pragma_table_info('{}') WHERE name=?1",
        escaped_table
    );

    conn.query_row(&sql, params![column], |row| row.get::<_, i32>(0))
        .map(|count| count > 0)
        .map_err(Into::into)
}

fn add_column_if_missing(
    conn: &Connection,
    table: &str,
    column: &str,
    definition: &str,
) -> Result<bool> {
    if !table_exists(conn, table)? || column_exists(conn, table, column)? {
        return Ok(false);
    }

    // Handle legacy bug where the column was added using only the type name
    if let Some(legacy_name) = legacy_column_name(definition) {
        if column_exists(conn, table, legacy_name)? {
            rename_column(conn, table, legacy_name, column)?;
            return Ok(true);
        }
    }

    let sql = format!(
        "ALTER TABLE {} ADD COLUMN {} {}",
        quote_ident(table),
        quote_ident(column),
        definition
    );
    conn.execute(&sql, [])?;
    Ok(true)
}

fn legacy_column_name(definition: &str) -> Option<&str> {
    definition.split_whitespace().next()
}

fn rename_column(conn: &Connection, table: &str, from: &str, to: &str) -> Result<()> {
    if !column_exists(conn, table, from)? {
        return Ok(());
    }

    let sql = format!(
        "ALTER TABLE {} RENAME COLUMN {} TO {}",
        quote_ident(table),
        quote_ident(from),
        quote_ident(to)
    );
    conn.execute(&sql, [])?;
    Ok(())
}

fn quote_ident(name: &str) -> String {
    let mut quoted = String::with_capacity(name.len() + 2);
    quoted.push('"');
    for ch in name.chars() {
        if ch == '"' {
            quoted.push('"');
            quoted.push('"');
        } else {
            quoted.push(ch);
        }
    }
    quoted.push('"');
    quoted
}

fn get_schema_version(conn: &Connection) -> Result<Option<String>> {
    match conn.query_row(
        "SELECT version FROM schema_version ORDER BY applied_at DESC LIMIT 1",
        [],
        |row| row.get(0),
    ) {
        Ok(version) => Ok(Some(version)),
        Err(rusqlite::Error::QueryReturnedNoRows) => Ok(None),
        Err(_) => Ok(None), // Table doesn't exist yet
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use tempfile::TempDir;

    #[test]
    fn test_db_creation() {
        let temp = TempDir::new().unwrap();
        config::set_test_biovault_home(temp.path());

        let db = BioVaultDb::new().unwrap();

        // Check that schema_version table exists
        let version: Option<String> = db
            .conn
            .query_row("SELECT version FROM schema_version LIMIT 1", [], |row| {
                row.get(0)
            })
            .ok();

        assert!(version.is_some());
        assert_eq!(version.unwrap(), "2.0.0");

        config::clear_test_biovault_home();
    }

    #[test]
    fn test_tables_created() {
        let temp = TempDir::new().unwrap();
        config::set_test_biovault_home(temp.path());

        let db = BioVaultDb::new().unwrap();

        // Check all expected tables exist
        let tables: Vec<String> = db
            .conn
            .prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
            .unwrap()
            .query_map([], |row| row.get(0))
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert!(tables.contains(&"schema_version".to_string()));
        assert!(tables.contains(&"messages".to_string()));
        assert!(tables.contains(&"participants".to_string()));
        assert!(tables.contains(&"files".to_string()));
        assert!(tables.contains(&"projects".to_string()));
        assert!(tables.contains(&"runs".to_string()));
        assert!(tables.contains(&"run_participants".to_string()));

        config::clear_test_biovault_home();
    }
}
