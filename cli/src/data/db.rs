use anyhow::{Context, Result};
use rusqlite::Connection;
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
        conn.execute_batch(schema)?;

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
        let column_exists: bool = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='data_type'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !column_exists {
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

        // Add inferred_sex to participants if it doesn't exist
        let sex_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('participants') WHERE name='inferred_sex'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !sex_exists {
            info!("Adding inferred_sex column to participants table");
            conn.execute("ALTER TABLE participants ADD COLUMN inferred_sex TEXT", [])?;
            info!("Migration complete: added inferred_sex to participants");
        }

        // Add status column to files if it doesn't exist
        let status_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='status'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !status_exists {
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
        let error_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='processing_error'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !error_exists {
            info!("Adding processing_error column to files table");
            conn.execute("ALTER TABLE files ADD COLUMN processing_error TEXT", [])?;
            info!("Migration complete: added processing_error column");
        }

        // Add queue_added_at column to files if it doesn't exist
        let queue_added_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='queue_added_at'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !queue_added_exists {
            info!("Adding queue_added_at column to files table");
            conn.execute("ALTER TABLE files ADD COLUMN queue_added_at DATETIME", [])?;
            info!("Migration complete: added queue_added_at column");
        }

        // Drop source and grch_version columns from files table if they exist
        // (these should only be in genotype_metadata table)
        let source_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='source'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        let grch_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('files') WHERE name='grch_version'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if source_exists || grch_exists {
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
        let inferred_sex_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('genotype_metadata') WHERE name='inferred_sex'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !inferred_sex_exists {
            info!("Adding inferred_sex column to genotype_metadata table");
            conn.execute(
                "ALTER TABLE genotype_metadata ADD COLUMN inferred_sex TEXT",
                [],
            )?;
            info!("Migration complete: added inferred_sex column");
        }

        // Add jupyter_port column to dev_envs if it doesn't exist
        let port_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('dev_envs') WHERE name='jupyter_port'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !port_exists {
            info!("Adding jupyter_port column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_port INTEGER", [])?;
            info!("Migration complete: added jupyter_port column");
        }

        // Add jupyter_pid column to dev_envs if it doesn't exist
        let pid_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('dev_envs') WHERE name='jupyter_pid'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !pid_exists {
            info!("Adding jupyter_pid column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_pid INTEGER", [])?;
            info!("Migration complete: added jupyter_pid column");
        }

        // Add jupyter_url column to dev_envs if it doesn't exist
        let url_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('dev_envs') WHERE name='jupyter_url'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !url_exists {
            info!("Adding jupyter_url column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_url TEXT", [])?;
            info!("Migration complete: added jupyter_url column");
        }

        // Add jupyter_token column to dev_envs if it doesn't exist
        let token_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('dev_envs') WHERE name='jupyter_token'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !token_exists {
            info!("Adding jupyter_token column to dev_envs table");
            conn.execute("ALTER TABLE dev_envs ADD COLUMN jupyter_token TEXT", [])?;
            info!("Migration complete: added jupyter_token column");
        }

        // Add metadata column to runs if it doesn't exist
        let metadata_exists = conn
            .query_row(
                "SELECT COUNT(*) FROM pragma_table_info('runs') WHERE name='metadata'",
                [],
                |row| row.get(0),
            )
            .map(|count: i32| count > 0)
            .unwrap_or(false);

        if !metadata_exists {
            info!("Adding metadata column to runs table");
            conn.execute("ALTER TABLE runs ADD COLUMN metadata TEXT", [])?;
            info!("Migration complete: added metadata column to runs");
        }

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
