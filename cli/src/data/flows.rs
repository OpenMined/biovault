use anyhow::Result;
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use super::BioVaultDb;
use crate::flow_spec::FlowSpec;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Flow {
    pub id: i64,
    pub name: String,
    pub flow_path: String,
    pub created_at: String,
    pub updated_at: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spec: Option<FlowSpec>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Run {
    pub id: i64,
    pub flow_id: Option<i64>,   // NULL for standalone module runs
    pub module_id: Option<i64>, // NULL for flow runs
    pub status: String,
    pub work_dir: String,
    pub results_dir: Option<String>,
    pub participant_count: Option<i32>, // Only for standalone module runs
    pub metadata: Option<String>, // JSON: { "input_overrides": {...}, "parameter_overrides": {...} }
    pub created_at: String,
    pub completed_at: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RunConfig {
    pub id: i64,
    pub flow_id: i64,
    pub name: String,
    pub config_data: serde_json::Value, // { "inputs": {...}, "parameters": {...} }
    pub created_at: String,
    pub updated_at: String,
}

impl BioVaultDb {
    /// List all flows
    pub fn list_flows(&self) -> Result<Vec<Flow>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, flow_path, created_at, updated_at
             FROM flows
             ORDER BY created_at DESC",
        )?;

        let mut flows = stmt
            .query_map([], |row| {
                Ok(Flow {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    flow_path: row.get(2)?,
                    created_at: row.get(3)?,
                    updated_at: row.get(4)?,
                    spec: None,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        // Load spec from flow.yaml for each flow
        for flow in &mut flows {
            let yaml_path = PathBuf::from(&flow.flow_path).join("flow.yaml");
            if yaml_path.exists() {
                match FlowSpec::load(&yaml_path) {
                    Ok(spec) => {
                        flow.spec = Some(spec);
                    }
                    Err(e) => {
                        eprintln!(
                            "Warning: Failed to parse flow.yaml for '{}': {}",
                            flow.name, e
                        );
                        eprintln!("  Path: {}", yaml_path.display());
                    }
                }
            }
        }

        Ok(flows)
    }

    /// Get flow by ID
    pub fn get_flow(&self, flow_id: i64) -> Result<Option<Flow>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, flow_path, created_at, updated_at
             FROM flows
             WHERE id = ?1",
        )?;

        let flow = stmt
            .query_row([flow_id], |row| {
                Ok(Flow {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    flow_path: row.get(2)?,
                    created_at: row.get(3)?,
                    updated_at: row.get(4)?,
                    spec: None,
                })
            })
            .optional()?;

        // Load spec if found
        if let Some(mut f) = flow {
            let yaml_path = PathBuf::from(&f.flow_path).join("flow.yaml");
            if yaml_path.exists() {
                match FlowSpec::load(&yaml_path) {
                    Ok(spec) => {
                        f.spec = Some(spec);
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to parse flow.yaml: {}", e);
                        eprintln!("  Path: {}", yaml_path.display());
                    }
                }
            }
            Ok(Some(f))
        } else {
            Ok(None)
        }
    }

    /// Register a flow in the database
    pub fn register_flow(&self, name: &str, flow_path: &str) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO flows (name, flow_path, created_at, updated_at)
             VALUES (?1, ?2, datetime('now'), datetime('now'))",
            params![name, flow_path],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update flow timestamp
    pub fn touch_flow(&self, flow_id: i64) -> Result<()> {
        self.conn.execute(
            "UPDATE flows SET updated_at = datetime('now') WHERE id = ?1",
            params![flow_id],
        )?;
        Ok(())
    }

    /// Delete a flow
    pub fn delete_flow(&self, flow_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM flows WHERE id = ?1", params![flow_id])?;
        Ok(())
    }

    /// List all runs (both flow and standalone module runs)
    pub fn list_runs(&self) -> Result<Vec<Run>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, flow_id, module_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM flow_runs
             ORDER BY created_at DESC",
        )?;

        let runs = stmt
            .query_map([], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    flow_id: row.get(1)?,
                    module_id: row.get(2)?,
                    status: row.get(3)?,
                    work_dir: row.get(4)?,
                    results_dir: row.get(5)?,
                    participant_count: row.get(6)?,
                    metadata: row.get(7)?,
                    created_at: row.get(8)?,
                    completed_at: row.get(9)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(runs)
    }

    /// List only flow runs
    pub fn list_flow_runs(&self) -> Result<Vec<Run>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, flow_id, module_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM flow_runs
             WHERE flow_id IS NOT NULL
             ORDER BY created_at DESC",
        )?;

        let runs = stmt
            .query_map([], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    flow_id: row.get(1)?,
                    module_id: row.get(2)?,
                    status: row.get(3)?,
                    work_dir: row.get(4)?,
                    results_dir: row.get(5)?,
                    participant_count: row.get(6)?,
                    metadata: row.get(7)?,
                    created_at: row.get(8)?,
                    completed_at: row.get(9)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(runs)
    }

    /// Get a specific run
    pub fn get_run(&self, run_id: i64) -> Result<Option<Run>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, flow_id, module_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM flow_runs
             WHERE id = ?1",
        )?;

        let run = stmt
            .query_row([run_id], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    flow_id: row.get(1)?,
                    module_id: row.get(2)?,
                    status: row.get(3)?,
                    work_dir: row.get(4)?,
                    results_dir: row.get(5)?,
                    participant_count: row.get(6)?,
                    metadata: row.get(7)?,
                    created_at: row.get(8)?,
                    completed_at: row.get(9)?,
                })
            })
            .optional()?;

        Ok(run)
    }

    /// Backwards compat alias
    pub fn get_flow_run(&self, run_id: i64) -> Result<Option<Run>> {
        self.get_run(run_id)
    }

    /// Create a flow run record
    pub fn create_flow_run(
        &self,
        flow_id: i64,
        work_dir: &str,
        results_dir: Option<&str>,
    ) -> Result<i64> {
        self.create_flow_run_with_metadata(flow_id, work_dir, results_dir, None)
    }

    /// Create a flow run record with metadata
    pub fn create_flow_run_with_metadata(
        &self,
        flow_id: i64,
        work_dir: &str,
        results_dir: Option<&str>,
        metadata: Option<&str>,
    ) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO flow_runs (flow_id, module_id, status, work_dir, results_dir, metadata, created_at)
             VALUES (?1, NULL, 'running', ?2, ?3, ?4, datetime('now'))",
            params![flow_id, work_dir, results_dir, metadata],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Create a standalone module run record
    pub fn create_module_run(
        &self,
        module_id: i64,
        work_dir: &str,
        participant_count: i32,
    ) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO flow_runs (flow_id, module_id, status, work_dir, participant_count, created_at)
             VALUES (NULL, ?1, 'running', ?2, ?3, datetime('now'))",
            params![module_id, work_dir, participant_count],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update run status (works for both flow and module runs)
    pub fn update_run_status(&self, run_id: i64, status: &str, completed: bool) -> Result<()> {
        if completed {
            self.conn.execute(
                "UPDATE flow_runs SET status = ?1, completed_at = datetime('now') WHERE id = ?2",
                params![status, run_id],
            )?;
        } else {
            self.conn.execute(
                "UPDATE flow_runs SET status = ?1 WHERE id = ?2",
                params![status, run_id],
            )?;
        }
        Ok(())
    }

    /// Backwards compat alias
    pub fn update_flow_run_status(&self, run_id: i64, status: &str, completed: bool) -> Result<()> {
        self.update_run_status(run_id, status, completed)
    }

    /// Delete a run (flow or module)
    pub fn delete_run(&self, run_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM flow_runs WHERE id = ?1", params![run_id])?;
        Ok(())
    }

    /// Backwards compat alias
    pub fn delete_flow_run(&self, run_id: i64) -> Result<()> {
        self.delete_run(run_id)
    }

    // =========================================================================
    // Run Configurations
    // =========================================================================

    /// Save a run configuration for a flow
    pub fn save_flow_run_config(
        &self,
        flow_id: i64,
        name: &str,
        config_data: &serde_json::Value,
    ) -> Result<i64> {
        let config_json = serde_json::to_string(config_data)?;

        self.conn.execute(
            "INSERT INTO flow_run_configs (flow_id, name, config_data, created_at, updated_at)
             VALUES (?1, ?2, ?3, datetime('now'), datetime('now'))",
            params![flow_id, name, config_json],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// List run configurations for a flow
    pub fn list_flow_run_configs(&self, flow_id: i64) -> Result<Vec<RunConfig>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, flow_id, name, config_data, created_at, updated_at
             FROM flow_run_configs
             WHERE flow_id = ?1
             ORDER BY created_at DESC
             LIMIT 10", // Keep last 10
        )?;

        let configs = stmt
            .query_map([flow_id], |row| {
                let config_json: String = row.get(3)?;
                let config_data =
                    serde_json::from_str(&config_json).unwrap_or(serde_json::json!({}));

                Ok(RunConfig {
                    id: row.get(0)?,
                    flow_id: row.get(1)?,
                    name: row.get(2)?,
                    config_data,
                    created_at: row.get(4)?,
                    updated_at: row.get(5)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(configs)
    }

    /// Get a specific run configuration
    pub fn get_flow_run_config(&self, config_id: i64) -> Result<Option<RunConfig>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, flow_id, name, config_data, created_at, updated_at
             FROM flow_run_configs
             WHERE id = ?1",
        )?;

        let config = stmt
            .query_row([config_id], |row| {
                let config_json: String = row.get(3)?;
                let config_data =
                    serde_json::from_str(&config_json).unwrap_or(serde_json::json!({}));

                Ok(RunConfig {
                    id: row.get(0)?,
                    flow_id: row.get(1)?,
                    name: row.get(2)?,
                    config_data,
                    created_at: row.get(4)?,
                    updated_at: row.get(5)?,
                })
            })
            .optional()?;

        Ok(config)
    }

    /// Delete a run configuration
    pub fn delete_flow_run_config(&self, config_id: i64) -> Result<()> {
        self.conn.execute(
            "DELETE FROM flow_run_configs WHERE id = ?1",
            params![config_id],
        )?;
        Ok(())
    }
}
