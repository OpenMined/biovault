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
    use tempfile::TempDir;

    #[test]
    fn test_db_creation() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

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
    }

    #[test]
    fn test_tables_created() {
        let temp = TempDir::new().unwrap();
        std::env::set_var("BIOVAULT_HOME", temp.path());

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
    }
}
