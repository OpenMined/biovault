use anyhow::Result;
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::path::Path;

use super::BioVaultDb;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DevEnv {
    pub id: i64,
    pub project_path: String,
    pub python_version: String,
    pub env_type: String,
    pub jupyter_installed: bool,
    pub jupyter_port: Option<i32>,
    pub jupyter_pid: Option<i32>,
    pub created_at: String,
    pub last_used_at: String,
}

impl BioVaultDb {
    /// Register a new development environment
    pub fn register_dev_env(
        &self,
        project_path: &Path,
        python_version: &str,
        env_type: &str,
        jupyter_installed: bool,
    ) -> Result<i64> {
        let path_str = project_path
            .canonicalize()?
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?
            .to_string();

        // Check if env already exists for this path
        if let Some(existing) = self.get_dev_env(&path_str)? {
            // Update existing
            self.conn.execute(
                "UPDATE dev_envs
                 SET python_version = ?1, env_type = ?2, jupyter_installed = ?3, last_used_at = CURRENT_TIMESTAMP
                 WHERE id = ?4",
                params![python_version, env_type, jupyter_installed as i32, existing.id],
            )?;
            return Ok(existing.id);
        }

        // Insert new
        self.conn.execute(
            "INSERT INTO dev_envs (project_path, python_version, env_type, jupyter_installed)
             VALUES (?1, ?2, ?3, ?4)",
            params![path_str, python_version, env_type, jupyter_installed as i32],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Get dev environment by project path
    pub fn get_dev_env(&self, project_path: &str) -> Result<Option<DevEnv>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, project_path, python_version, env_type, jupyter_installed, jupyter_port, jupyter_pid, created_at, last_used_at
             FROM dev_envs
             WHERE project_path = ?1",
        )?;

        let env = stmt
            .query_row([project_path], |row| {
                Ok(DevEnv {
                    id: row.get(0)?,
                    project_path: row.get(1)?,
                    python_version: row.get(2)?,
                    env_type: row.get(3)?,
                    jupyter_installed: row.get::<_, i32>(4)? != 0,
                    jupyter_port: row.get(5)?,
                    jupyter_pid: row.get(6)?,
                    created_at: row.get(7)?,
                    last_used_at: row.get(8)?,
                })
            })
            .optional()?;

        Ok(env)
    }

    /// List all development environments
    pub fn list_dev_envs(&self) -> Result<Vec<DevEnv>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, project_path, python_version, env_type, jupyter_installed, jupyter_port, jupyter_pid, created_at, last_used_at
             FROM dev_envs
             ORDER BY last_used_at DESC",
        )?;

        let envs = stmt
            .query_map([], |row| {
                Ok(DevEnv {
                    id: row.get(0)?,
                    project_path: row.get(1)?,
                    python_version: row.get(2)?,
                    env_type: row.get(3)?,
                    jupyter_installed: row.get::<_, i32>(4)? != 0,
                    jupyter_port: row.get(5)?,
                    jupyter_pid: row.get(6)?,
                    created_at: row.get(7)?,
                    last_used_at: row.get(8)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(envs)
    }

    /// Update last used timestamp for a dev environment
    pub fn update_dev_env_last_used(&self, project_path: &Path) -> Result<()> {
        let canonical_path = project_path.canonicalize()?;
        let path_str = canonical_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?;

        self.conn.execute(
            "UPDATE dev_envs SET last_used_at = CURRENT_TIMESTAMP WHERE project_path = ?1",
            params![path_str],
        )?;

        Ok(())
    }

    /// Update Jupyter session info (port and PID)
    pub fn update_jupyter_session(
        &self,
        project_path: &Path,
        port: Option<i32>,
        pid: Option<i32>,
    ) -> Result<()> {
        let canonical_path = project_path.canonicalize()?;
        let path_str = canonical_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?;

        self.conn.execute(
            "UPDATE dev_envs SET jupyter_port = ?1, jupyter_pid = ?2, last_used_at = CURRENT_TIMESTAMP WHERE project_path = ?3",
            params![port, pid, path_str],
        )?;

        Ok(())
    }

    /// Delete a dev environment record
    pub fn delete_dev_env(&self, project_path: &str) -> Result<()> {
        self.conn.execute(
            "DELETE FROM dev_envs WHERE project_path = ?1",
            params![project_path],
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    fn setup_test_db() -> (TempDir, BioVaultDb) {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path());
        let db = BioVaultDb::new().unwrap();
        (tmp, db)
    }

    fn teardown_test() {
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_register_and_list_dev_envs() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        let id = db
            .register_dev_env(&project_path, "3.12", "jupyter", true)
            .unwrap();
        assert!(id > 0);

        let envs = db.list_dev_envs().unwrap();
        assert_eq!(envs.len(), 1);
        assert_eq!(envs[0].python_version, "3.12");
        assert!(envs[0].jupyter_installed);

        teardown_test();
    }

    #[test]
    fn test_get_dev_env() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_dev_env(&project_path, "3.12", "jupyter", true)
            .unwrap();

        let canonical_path = project_path.canonicalize().unwrap();
        let env = db.get_dev_env(canonical_path.to_str().unwrap()).unwrap();
        assert!(env.is_some());
        assert_eq!(env.unwrap().python_version, "3.12");

        teardown_test();
    }

    #[test]
    fn test_update_dev_env_last_used() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_dev_env(&project_path, "3.12", "jupyter", true)
            .unwrap();

        // Small delay to ensure timestamp changes
        std::thread::sleep(std::time::Duration::from_millis(10));

        db.update_dev_env_last_used(&project_path).unwrap();

        let canonical_path = project_path.canonicalize().unwrap();
        let env = db.get_dev_env(canonical_path.to_str().unwrap()).unwrap();
        assert!(env.is_some());

        teardown_test();
    }

    #[test]
    fn test_update_existing_dev_env() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        // Register first time
        let id1 = db
            .register_dev_env(&project_path, "3.11", "jupyter", false)
            .unwrap();

        // Register again with different values - should update
        let id2 = db
            .register_dev_env(&project_path, "3.12", "jupyter", true)
            .unwrap();

        assert_eq!(id1, id2); // Same ID means it was updated, not inserted

        let canonical_path = project_path.canonicalize().unwrap();
        let env = db
            .get_dev_env(canonical_path.to_str().unwrap())
            .unwrap()
            .unwrap();
        assert_eq!(env.python_version, "3.12");
        assert!(env.jupyter_installed);

        teardown_test();
    }
}
