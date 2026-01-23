use anyhow::{Context, Result};
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::path::Path;

use super::BioVaultDb;
use crate::module_spec::ModuleFile;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Module {
    pub id: i64,
    pub name: String,
    pub version: String,
    pub author: String,
    pub workflow: String,
    pub template: String,
    pub module_path: String,
    pub created_at: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ModuleYaml {
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

impl ModuleYaml {
    pub fn from_str(raw: &str) -> Result<Self> {
        let module: ModuleFile = serde_yaml::from_str(raw)?;
        Ok(ModuleYaml::from_module_file(module))
    }

    pub fn from_file(path: &Path) -> Result<Self> {
        let raw = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read module spec at {}", path.display()))?;
        ModuleYaml::from_str(&raw)
    }

    fn from_module_file(module: ModuleFile) -> Self {
        let author = module.metadata.authors.first().cloned().unwrap_or_default();
        let runner = module.spec.runner.clone().unwrap_or_default();
        let workflow = runner
            .entrypoint
            .clone()
            .unwrap_or_else(|| "workflow.nf".to_string());
        let template = runner
            .template
            .clone()
            .unwrap_or_else(|| "dynamic-nextflow".to_string());
        let assets = module
            .spec
            .assets
            .unwrap_or_default()
            .into_iter()
            .map(|asset| asset.path)
            .collect();

        ModuleYaml {
            name: module.metadata.name,
            version: module.metadata.version,
            author,
            workflow,
            template,
            assets,
        }
    }
}

pub struct UpdateModuleParams<'a> {
    pub name: &'a str,
    pub version: &'a str,
    pub author: &'a str,
    pub workflow: &'a str,
    pub template: &'a str,
    pub module_path: &'a Path,
}

impl BioVaultDb {
    /// List all modules
    pub fn list_modules(&self) -> Result<Vec<Module>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, version, author, workflow, template, module_path, created_at
             FROM modules
             ORDER BY created_at DESC",
        )?;

        let modules = stmt
            .query_map([], |row| {
                Ok(Module {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    version: row.get(2)?,
                    author: row.get(3)?,
                    workflow: row.get(4)?,
                    template: row.get(5)?,
                    module_path: row.get(6)?,
                    created_at: row.get(7)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(modules)
    }

    /// Get module by name, name@version, or ID
    pub fn get_module(&self, identifier: &str) -> Result<Option<Module>> {
        // Try parsing as ID first
        if let Ok(id) = identifier.parse::<i64>() {
            let mut stmt = self.conn.prepare(
                "SELECT id, name, version, author, workflow, template, module_path, created_at
                 FROM modules
                 WHERE id = ?1",
            )?;

            let module = stmt
                .query_row([id], |row| {
                    Ok(Module {
                        id: row.get(0)?,
                        name: row.get(1)?,
                        version: row.get(2)?,
                        author: row.get(3)?,
                        workflow: row.get(4)?,
                        template: row.get(5)?,
                        module_path: row.get(6)?,
                        created_at: row.get(7)?,
                    })
                })
                .optional()?;

            return Ok(module);
        }

        // Check if identifier includes version (name@version)
        if let Some((name, version)) = identifier.split_once('@') {
            let mut stmt = self.conn.prepare(
                "SELECT id, name, version, author, workflow, template, module_path, created_at
                 FROM modules
                 WHERE name = ?1 AND version = ?2",
            )?;

            let module = stmt
                .query_row([name, version], |row| {
                    Ok(Module {
                        id: row.get(0)?,
                        name: row.get(1)?,
                        version: row.get(2)?,
                        author: row.get(3)?,
                        workflow: row.get(4)?,
                        template: row.get(5)?,
                        module_path: row.get(6)?,
                        created_at: row.get(7)?,
                    })
                })
                .optional()?;

            return Ok(module);
        }

        // Otherwise treat as name, get latest version
        let mut stmt = self.conn.prepare(
            "SELECT id, name, version, author, workflow, template, module_path, created_at
             FROM modules
             WHERE name = ?1
             ORDER BY created_at DESC
             LIMIT 1",
        )?;

        let module = stmt
            .query_row([identifier], |row| {
                Ok(Module {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    version: row.get(2)?,
                    author: row.get(3)?,
                    workflow: row.get(4)?,
                    template: row.get(5)?,
                    module_path: row.get(6)?,
                    created_at: row.get(7)?,
                })
            })
            .optional()?;

        Ok(module)
    }

    /// Register a module in the database
    pub fn register_module(
        &self,
        name: &str,
        version: &str,
        author: &str,
        workflow: &str,
        template: &str,
        module_path: &Path,
    ) -> Result<i64> {
        let path_str = module_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid module path"))?;

        // Check if module with this name and version already exists
        let identifier = format!("{}@{}", name, version);
        if let Some(existing) = self.get_module(&identifier)? {
            anyhow::bail!(
                "Module '{}' version {} already exists (id: {}). Use --overwrite to replace.",
                name,
                version,
                existing.id
            );
        }

        self.conn.execute(
            "INSERT INTO modules (name, version, author, workflow, template, module_path)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![name, version, author, workflow, template, path_str],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update an existing module
    pub fn update_module(
        &self,
        name: &str,
        version: &str,
        author: &str,
        workflow: &str,
        template: &str,
        module_path: &Path,
    ) -> Result<()> {
        let path_str = module_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid module path"))?;

        let rows_affected = self.conn.execute(
            "UPDATE modules
             SET author = ?1, workflow = ?2, template = ?3, module_path = ?4
             WHERE name = ?5 AND version = ?6",
            params![author, workflow, template, path_str, name, version],
        )?;

        if rows_affected == 0 {
            anyhow::bail!("Module '{}' version {} not found", name, version);
        }

        Ok(())
    }

    /// Update an existing module by ID
    pub fn update_module_by_id(
        &self,
        module_id: i64,
        params: UpdateModuleParams<'_>,
    ) -> Result<()> {
        let name = params.name;
        let version = params.version;
        let author = params.author;
        let workflow = params.workflow;
        let template = params.template;
        let module_path = params.module_path;
        let path_str = module_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid module path"))?;

        // Check if another module has the same (name, version) combination
        let mut stmt = self
            .conn
            .prepare("SELECT id FROM modules WHERE name = ?1 AND version = ?2 AND id != ?3")
            .context("Failed to prepare module uniqueness query")?;

        let conflict = stmt
            .query_row(params![name, version, module_id], |row| {
                row.get::<_, i64>(0)
            })
            .optional()?;

        if conflict.is_some() {
            anyhow::bail!(
                "Module '{}' version {} is already used by a different module",
                name,
                version
            );
        }

        let rows_affected = self
            .conn
            .execute(
                "UPDATE modules
                 SET name = ?1, version = ?2, author = ?3, workflow = ?4, template = ?5, module_path = ?6
                 WHERE id = ?7",
                params![name, version, author, workflow, template, path_str, module_id],
            )
            .context("Failed to update module record")?;

        if rows_affected == 0 {
            anyhow::bail!("Module id {} not found", module_id);
        }

        Ok(())
    }

    fn runs_module_reference_column(&self) -> Result<Option<&'static str>> {
        if self.table_has_column("flow_runs", "module_id")? {
            Ok(Some("module_id"))
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

    /// Delete a module from the database
    pub fn delete_module(&self, identifier: &str) -> Result<Module> {
        // Get the module first
        let module = self
            .get_module(identifier)?
            .ok_or_else(|| anyhow::anyhow!("Module '{}' not found", identifier))?;

        if let Some(column) = self.runs_module_reference_column()? {
            let query = format!("SELECT COUNT(*) FROM flow_runs WHERE {} = ?1", column);
            let run_count: i64 = self
                .conn
                .query_row(&query, params![module.id], |row| row.get(0))?;

            if run_count > 0 {
                anyhow::bail!(
                    "Cannot delete module '{}': {} associated run(s) exist. Delete runs first.",
                    module.name,
                    run_count
                );
            }
        }

        // Delete the module
        self.conn
            .execute("DELETE FROM modules WHERE id = ?1", params![module.id])?;

        Ok(module)
    }

    /// Count runs for a module
    pub fn count_module_runs(&self, module_id: i64) -> Result<i64> {
        if let Some(column) = self.runs_module_reference_column()? {
            let query = format!("SELECT COUNT(*) FROM flow_runs WHERE {} = ?1", column);
            let count: i64 = self
                .conn
                .query_row(&query, params![module_id], |row| row.get(0))?;
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
    fn test_register_and_list_modules() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        let id = db
            .register_module(
                "test",
                "1.0.0",
                "author@example.com",
                "workflow.nf",
                "default",
                &module_path,
            )
            .unwrap();
        assert!(id > 0);

        let modules = db.list_modules().unwrap();
        assert_eq!(modules.len(), 1);
        assert_eq!(modules[0].name, "test");
        assert_eq!(modules[0].version, "1.0.0");
        assert_eq!(modules[0].author, "author@example.com");

        teardown_test();
    }

    #[test]
    fn test_get_module_by_name() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        db.register_module(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &module_path,
        )
        .unwrap();

        let module = db.get_module("test").unwrap();
        assert!(module.is_some());
        let p = module.unwrap();
        assert_eq!(p.name, "test");
        assert_eq!(p.version, "1.0.0");

        teardown_test();
    }

    #[test]
    fn test_get_module_by_id() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        let id = db
            .register_module(
                "test",
                "1.0.0",
                "author@example.com",
                "workflow.nf",
                "default",
                &module_path,
            )
            .unwrap();

        let module = db.get_module(&id.to_string()).unwrap();
        assert!(module.is_some());
        assert_eq!(module.unwrap().id, id);

        teardown_test();
    }

    #[test]
    fn test_duplicate_module_fails() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        db.register_module(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &module_path,
        )
        .unwrap();

        let result = db.register_module(
            "test",
            "1.0.0",
            "other@example.com",
            "workflow.nf",
            "default",
            &module_path,
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already exists"));

        teardown_test();
    }

    #[test]
    fn test_update_module() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        db.register_module(
            "test",
            "1.0.0",
            "old@example.com",
            "old.nf",
            "old",
            &module_path,
        )
        .unwrap();

        db.update_module(
            "test",
            "1.0.0",
            "new@example.com",
            "new.nf",
            "new",
            &module_path,
        )
        .unwrap();

        let module = db.get_module("test").unwrap().unwrap();
        assert_eq!(module.author, "new@example.com");
        assert_eq!(module.workflow, "new.nf");

        teardown_test();
    }

    #[test]
    fn test_delete_module() {
        let (tmp, db) = setup_test_db();
        let module_path = tmp.path().join("test-module");
        fs::create_dir_all(&module_path).unwrap();

        db.register_module(
            "test",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &module_path,
        )
        .unwrap();

        let deleted = db.delete_module("test").unwrap();
        assert_eq!(deleted.name, "test");
        assert_eq!(deleted.version, "1.0.0");

        let module = db.get_module("test").unwrap();
        assert!(module.is_none());

        teardown_test();
    }

    #[test]
    fn test_delete_nonexistent_module_fails() {
        let (_tmp, db) = setup_test_db();

        let result = db.delete_module("nonexistent");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));

        teardown_test();
    }
}
