use anyhow::{Context, Result};
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::path::Path;

use super::BioVaultDb;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Project {
    pub id: i64,
    pub name: String,
    pub version: String,
    pub author: String,
    pub workflow: String,
    pub template: String,
    pub project_path: String,
    pub created_at: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ProjectYaml {
    pub name: String,
    #[serde(default = "default_version")]
    pub version: String,
    pub author: String,
    pub workflow: String,
    pub template: String,
    pub assets: Vec<String>,
}

fn default_version() -> String {
    "1.0.0".to_string()
}

pub struct UpdateProjectParams<'a> {
    pub name: &'a str,
    pub version: &'a str,
    pub author: &'a str,
    pub workflow: &'a str,
    pub template: &'a str,
    pub project_path: &'a Path,
}

impl BioVaultDb {
    /// List all projects
    pub fn list_projects(&self) -> Result<Vec<Project>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, version, author, workflow, template, project_path, created_at
             FROM projects
             ORDER BY created_at DESC",
        )?;

        let projects = stmt
            .query_map([], |row| {
                Ok(Project {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    version: row.get(2)?,
                    author: row.get(3)?,
                    workflow: row.get(4)?,
                    template: row.get(5)?,
                    project_path: row.get(6)?,
                    created_at: row.get(7)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(projects)
    }

    /// Get project by name, name@version, or ID
    pub fn get_project(&self, identifier: &str) -> Result<Option<Project>> {
        // Try parsing as ID first
        if let Ok(id) = identifier.parse::<i64>() {
            let mut stmt = self.conn.prepare(
                "SELECT id, name, version, author, workflow, template, project_path, created_at
                 FROM projects
                 WHERE id = ?1",
            )?;

            let project = stmt
                .query_row([id], |row| {
                    Ok(Project {
                        id: row.get(0)?,
                        name: row.get(1)?,
                        version: row.get(2)?,
                        author: row.get(3)?,
                        workflow: row.get(4)?,
                        template: row.get(5)?,
                        project_path: row.get(6)?,
                        created_at: row.get(7)?,
                    })
                })
                .optional()?;

            return Ok(project);
        }

        // Check if identifier includes version (name@version)
        if let Some((name, version)) = identifier.split_once('@') {
            let mut stmt = self.conn.prepare(
                "SELECT id, name, version, author, workflow, template, project_path, created_at
                 FROM projects
                 WHERE name = ?1 AND version = ?2",
            )?;

            let project = stmt
                .query_row([name, version], |row| {
                    Ok(Project {
                        id: row.get(0)?,
                        name: row.get(1)?,
                        version: row.get(2)?,
                        author: row.get(3)?,
                        workflow: row.get(4)?,
                        template: row.get(5)?,
                        project_path: row.get(6)?,
                        created_at: row.get(7)?,
                    })
                })
                .optional()?;

            return Ok(project);
        }

        // Otherwise treat as name, get latest version
        let mut stmt = self.conn.prepare(
            "SELECT id, name, version, author, workflow, template, project_path, created_at
             FROM projects
             WHERE name = ?1
             ORDER BY created_at DESC
             LIMIT 1",
        )?;

        let project = stmt
            .query_row([identifier], |row| {
                Ok(Project {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    version: row.get(2)?,
                    author: row.get(3)?,
                    workflow: row.get(4)?,
                    template: row.get(5)?,
                    project_path: row.get(6)?,
                    created_at: row.get(7)?,
                })
            })
            .optional()?;

        Ok(project)
    }

    /// Register a project in the database
    pub fn register_project(
        &self,
        name: &str,
        version: &str,
        author: &str,
        workflow: &str,
        template: &str,
        project_path: &Path,
    ) -> Result<i64> {
        let path_str = project_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?;

        // Check if project with this name and version already exists
        let identifier = format!("{}@{}", name, version);
        if let Some(existing) = self.get_project(&identifier)? {
            anyhow::bail!(
                "Project '{}' version {} already exists (id: {}). Use --overwrite to replace.",
                name,
                version,
                existing.id
            );
        }

        self.conn.execute(
            "INSERT INTO projects (name, version, author, workflow, template, project_path)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![name, version, author, workflow, template, path_str],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update an existing project
    pub fn update_project(
        &self,
        name: &str,
        version: &str,
        author: &str,
        workflow: &str,
        template: &str,
        project_path: &Path,
    ) -> Result<()> {
        let path_str = project_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?;

        let rows_affected = self.conn.execute(
            "UPDATE projects
             SET author = ?1, workflow = ?2, template = ?3, project_path = ?4
             WHERE name = ?5 AND version = ?6",
            params![author, workflow, template, path_str, name, version],
        )?;

        if rows_affected == 0 {
            anyhow::bail!("Project '{}' version {} not found", name, version);
        }

        Ok(())
    }

    /// Update an existing project by ID
    pub fn update_project_by_id(
        &self,
        project_id: i64,
        params: UpdateProjectParams<'_>,
    ) -> Result<()> {
        let name = params.name;
        let version = params.version;
        let author = params.author;
        let workflow = params.workflow;
        let template = params.template;
        let project_path = params.project_path;
        let path_str = project_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid project path"))?;

        // Check if another project has the same (name, version) combination
        let mut stmt = self
            .conn
            .prepare("SELECT id FROM projects WHERE name = ?1 AND version = ?2 AND id != ?3")
            .context("Failed to prepare project uniqueness query")?;

        let conflict = stmt
            .query_row(params![name, version, project_id], |row| {
                row.get::<_, i64>(0)
            })
            .optional()?;

        if conflict.is_some() {
            anyhow::bail!(
                "Project '{}' version {} is already used by a different project",
                name,
                version
            );
        }

        let rows_affected = self
            .conn
            .execute(
                "UPDATE projects
                 SET name = ?1, version = ?2, author = ?3, workflow = ?4, template = ?5, project_path = ?6
                 WHERE id = ?7",
                params![name, version, author, workflow, template, path_str, project_id],
            )
            .context("Failed to update project record")?;

        if rows_affected == 0 {
            anyhow::bail!("Project id {} not found", project_id);
        }

        Ok(())
    }

    fn runs_project_reference_column(&self) -> Result<Option<&'static str>> {
        if self.table_has_column("runs", "project_id")? {
            Ok(Some("project_id"))
        } else if self.table_has_column("runs", "step_id")? {
            Ok(Some("step_id"))
        } else {
            Ok(None)
        }
    }

    fn table_has_column(&self, table: &str, column: &str) -> Result<bool> {
        // table/column names are controlled internally, so simple escaping is sufficient
        let table = table.replace('\'', "''");
        let column = column.replace('\'', "''");
        let sql = format!(
            "SELECT COUNT(*) FROM pragma_table_info('{}') WHERE name='{}'",
            table, column
        );

        let count: i64 = self.conn.query_row(&sql, [], |row| row.get(0))?;
        Ok(count > 0)
    }

    /// Delete a project from the database
    pub fn delete_project(&self, identifier: &str) -> Result<Project> {
        // Get the project first
        let project = self
            .get_project(identifier)?
            .ok_or_else(|| anyhow::anyhow!("Project '{}' not found", identifier))?;

        if let Some(column) = self.runs_project_reference_column()? {
            let query = format!("SELECT COUNT(*) FROM runs WHERE {} = ?1", column);
            let run_count: i64 = self
                .conn
                .query_row(&query, params![project.id], |row| row.get(0))?;

            if run_count > 0 {
                anyhow::bail!(
                    "Cannot delete project '{}': {} associated run(s) exist. Delete runs first.",
                    project.name,
                    run_count
                );
            }
        }

        // Delete the project
        self.conn
            .execute("DELETE FROM projects WHERE id = ?1", params![project.id])?;

        Ok(project)
    }

    /// Count runs for a project
    pub fn count_project_runs(&self, project_id: i64) -> Result<i64> {
        if let Some(column) = self.runs_project_reference_column()? {
            let query = format!("SELECT COUNT(*) FROM runs WHERE {} = ?1", column);
            let count: i64 = self
                .conn
                .query_row(&query, params![project_id], |row| row.get(0))?;
            Ok(count)
        } else {
            Ok(0)
        }
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
    fn test_register_and_list_projects() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        let id = db
            .register_project(
                "test",
                "1.0.0",
                "author@example.com",
                "workflow.nf",
                "default",
                &project_path,
            )
            .unwrap();
        assert!(id > 0);

        let projects = db.list_projects().unwrap();
        assert_eq!(projects.len(), 1);
        assert_eq!(projects[0].name, "test");
        assert_eq!(projects[0].version, "1.0.0");
        assert_eq!(projects[0].author, "author@example.com");

        teardown_test();
    }

    #[test]
    fn test_get_project_by_name() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_project(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &project_path,
        )
        .unwrap();

        let project = db.get_project("test").unwrap();
        assert!(project.is_some());
        let p = project.unwrap();
        assert_eq!(p.name, "test");
        assert_eq!(p.version, "1.0.0");

        teardown_test();
    }

    #[test]
    fn test_get_project_by_id() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        let id = db
            .register_project(
                "test",
                "1.0.0",
                "author@example.com",
                "workflow.nf",
                "default",
                &project_path,
            )
            .unwrap();

        let project = db.get_project(&id.to_string()).unwrap();
        assert!(project.is_some());
        assert_eq!(project.unwrap().id, id);

        teardown_test();
    }

    #[test]
    fn test_duplicate_project_fails() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_project(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &project_path,
        )
        .unwrap();

        let result = db.register_project(
            "test",
            "1.0.0",
            "other@example.com",
            "workflow.nf",
            "default",
            &project_path,
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already exists"));

        teardown_test();
    }

    #[test]
    fn test_update_project() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_project(
            "test",
            "1.0.0",
            "old@example.com",
            "old.nf",
            "old",
            &project_path,
        )
        .unwrap();

        db.update_project(
            "test",
            "1.0.0",
            "new@example.com",
            "new.nf",
            "new",
            &project_path,
        )
        .unwrap();

        let project = db.get_project("test").unwrap().unwrap();
        assert_eq!(project.author, "new@example.com");
        assert_eq!(project.workflow, "new.nf");

        teardown_test();
    }

    #[test]
    fn test_delete_project() {
        let (tmp, db) = setup_test_db();
        let project_path = tmp.path().join("test-project");
        fs::create_dir_all(&project_path).unwrap();

        db.register_project(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &project_path,
        )
        .unwrap();

        let deleted = db.delete_project("test").unwrap();
        assert_eq!(deleted.name, "test");
        assert_eq!(deleted.version, "1.0.0");

        let project = db.get_project("test").unwrap();
        assert!(project.is_none());

        teardown_test();
    }

    #[test]
    fn test_delete_nonexistent_project_fails() {
        let (_tmp, db) = setup_test_db();

        let result = db.delete_project("nonexistent");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));

        teardown_test();
    }
}
