use anyhow::Result;
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use super::BioVaultDb;
use crate::pipeline_spec::{resolve_pipeline_spec_path, PipelineSpec};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Pipeline {
    pub id: i64,
    pub name: String,
    pub pipeline_path: String,
    pub created_at: String,
    pub updated_at: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spec: Option<PipelineSpec>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Run {
    pub id: i64,
    pub pipeline_id: Option<i64>, // NULL for standalone step runs
    pub step_id: Option<i64>,     // NULL for pipeline runs
    pub status: String,
    pub work_dir: String,
    pub results_dir: Option<String>,
    pub participant_count: Option<i32>, // Only for step runs
    pub metadata: Option<String>, // JSON: { "input_overrides": {...}, "parameter_overrides": {...} }
    pub created_at: String,
    pub completed_at: Option<String>,
}

// Backwards compat alias
pub type PipelineRun = Run;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RunConfig {
    pub id: i64,
    pub pipeline_id: i64,
    pub name: String,
    pub config_data: serde_json::Value, // { "inputs": {...}, "parameters": {...} }
    pub created_at: String,
    pub updated_at: String,
}

impl BioVaultDb {
    /// List all pipelines
    pub fn list_pipelines(&self) -> Result<Vec<Pipeline>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, pipeline_path, created_at, updated_at
             FROM pipelines
             ORDER BY created_at DESC",
        )?;

        let mut pipelines = stmt
            .query_map([], |row| {
                Ok(Pipeline {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    pipeline_path: row.get(2)?,
                    created_at: row.get(3)?,
                    updated_at: row.get(4)?,
                    spec: None,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        // Load spec from flow.yaml/pipeline.yaml for each pipeline
        for pipeline in &mut pipelines {
            let yaml_path =
                resolve_pipeline_spec_path(PathBuf::from(&pipeline.pipeline_path).as_path());
            if yaml_path.exists() {
                match PipelineSpec::load(&yaml_path) {
                    Ok(spec) => {
                        pipeline.spec = Some(spec);
                    }
                    Err(e) => {
                        eprintln!(
                            "Warning: Failed to parse flow spec for '{}': {}",
                            pipeline.name, e
                        );
                        eprintln!("  Path: {}", yaml_path.display());
                    }
                }
            }
        }

        Ok(pipelines)
    }

    /// Get pipeline by ID
    pub fn get_pipeline(&self, pipeline_id: i64) -> Result<Option<Pipeline>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, name, pipeline_path, created_at, updated_at
             FROM pipelines
             WHERE id = ?1",
        )?;

        let pipeline = stmt
            .query_row([pipeline_id], |row| {
                Ok(Pipeline {
                    id: row.get(0)?,
                    name: row.get(1)?,
                    pipeline_path: row.get(2)?,
                    created_at: row.get(3)?,
                    updated_at: row.get(4)?,
                    spec: None,
                })
            })
            .optional()?;

        // Load spec if found
        if let Some(mut p) = pipeline {
            let yaml_path = resolve_pipeline_spec_path(PathBuf::from(&p.pipeline_path).as_path());
            if yaml_path.exists() {
                match PipelineSpec::load(&yaml_path) {
                    Ok(spec) => {
                        p.spec = Some(spec);
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to parse flow spec: {}", e);
                        eprintln!("  Path: {}", yaml_path.display());
                    }
                }
            }
            Ok(Some(p))
        } else {
            Ok(None)
        }
    }

    /// Register a pipeline in the database
    pub fn register_pipeline(&self, name: &str, pipeline_path: &str) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO pipelines (name, pipeline_path, created_at, updated_at)
             VALUES (?1, ?2, datetime('now'), datetime('now'))",
            params![name, pipeline_path],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update pipeline timestamp
    pub fn touch_pipeline(&self, pipeline_id: i64) -> Result<()> {
        self.conn.execute(
            "UPDATE pipelines SET updated_at = datetime('now') WHERE id = ?1",
            params![pipeline_id],
        )?;
        Ok(())
    }

    /// Delete a pipeline
    pub fn delete_pipeline(&self, pipeline_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM pipelines WHERE id = ?1", params![pipeline_id])?;
        Ok(())
    }

    /// List all runs (both pipeline and standalone step runs)
    pub fn list_runs(&self) -> Result<Vec<Run>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, step_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM runs
             ORDER BY created_at DESC",
        )?;

        let runs = stmt
            .query_map([], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
                    step_id: row.get(2)?,
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

    /// List only pipeline runs
    pub fn list_pipeline_runs(&self) -> Result<Vec<Run>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, step_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM runs
             WHERE pipeline_id IS NOT NULL
             ORDER BY created_at DESC",
        )?;

        let runs = stmt
            .query_map([], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
                    step_id: row.get(2)?,
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
            "SELECT id, pipeline_id, step_id, status, work_dir, results_dir, participant_count, metadata, created_at, completed_at
             FROM runs
             WHERE id = ?1",
        )?;

        let run = stmt
            .query_row([run_id], |row| {
                Ok(Run {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
                    step_id: row.get(2)?,
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
    pub fn get_pipeline_run(&self, run_id: i64) -> Result<Option<Run>> {
        self.get_run(run_id)
    }

    /// Create a pipeline run record
    pub fn create_pipeline_run(
        &self,
        pipeline_id: i64,
        work_dir: &str,
        results_dir: Option<&str>,
    ) -> Result<i64> {
        self.create_pipeline_run_with_metadata(pipeline_id, work_dir, results_dir, None)
    }

    /// Create a pipeline run record with metadata
    pub fn create_pipeline_run_with_metadata(
        &self,
        pipeline_id: i64,
        work_dir: &str,
        results_dir: Option<&str>,
        metadata: Option<&str>,
    ) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO runs (pipeline_id, step_id, status, work_dir, results_dir, metadata, created_at)
             VALUES (?1, NULL, 'running', ?2, ?3, ?4, datetime('now'))",
            params![pipeline_id, work_dir, results_dir, metadata],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Create a standalone step run record
    pub fn create_step_run(
        &self,
        step_id: i64,
        work_dir: &str,
        participant_count: i32,
    ) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO runs (pipeline_id, step_id, status, work_dir, participant_count, created_at)
             VALUES (NULL, ?1, 'running', ?2, ?3, datetime('now'))",
            params![step_id, work_dir, participant_count],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update run status (works for both pipeline and step runs)
    pub fn update_run_status(&self, run_id: i64, status: &str, completed: bool) -> Result<()> {
        if completed {
            self.conn.execute(
                "UPDATE runs SET status = ?1, completed_at = datetime('now') WHERE id = ?2",
                params![status, run_id],
            )?;
        } else {
            self.conn.execute(
                "UPDATE runs SET status = ?1 WHERE id = ?2",
                params![status, run_id],
            )?;
        }
        Ok(())
    }

    /// Backwards compat alias
    pub fn update_pipeline_run_status(
        &self,
        run_id: i64,
        status: &str,
        completed: bool,
    ) -> Result<()> {
        self.update_run_status(run_id, status, completed)
    }

    /// Delete a run (pipeline or step)
    pub fn delete_run(&self, run_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM runs WHERE id = ?1", params![run_id])?;
        Ok(())
    }

    /// Backwards compat alias
    pub fn delete_pipeline_run(&self, run_id: i64) -> Result<()> {
        self.delete_run(run_id)
    }

    // =========================================================================
    // Run Configurations
    // =========================================================================

    /// Save a run configuration for a pipeline
    pub fn save_run_config(
        &self,
        pipeline_id: i64,
        name: &str,
        config_data: &serde_json::Value,
    ) -> Result<i64> {
        let config_json = serde_json::to_string(config_data)?;

        self.conn.execute(
            "INSERT INTO run_configs (pipeline_id, name, config_data, created_at, updated_at)
             VALUES (?1, ?2, ?3, datetime('now'), datetime('now'))",
            params![pipeline_id, name, config_json],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// List run configurations for a pipeline
    pub fn list_run_configs(&self, pipeline_id: i64) -> Result<Vec<RunConfig>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, name, config_data, created_at, updated_at
             FROM run_configs
             WHERE pipeline_id = ?1
             ORDER BY created_at DESC
             LIMIT 10", // Keep last 10
        )?;

        let configs = stmt
            .query_map([pipeline_id], |row| {
                let config_json: String = row.get(3)?;
                let config_data =
                    serde_json::from_str(&config_json).unwrap_or(serde_json::json!({}));

                Ok(RunConfig {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
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
    pub fn get_run_config(&self, config_id: i64) -> Result<Option<RunConfig>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, name, config_data, created_at, updated_at
             FROM run_configs
             WHERE id = ?1",
        )?;

        let config = stmt
            .query_row([config_id], |row| {
                let config_json: String = row.get(3)?;
                let config_data =
                    serde_json::from_str(&config_json).unwrap_or(serde_json::json!({}));

                Ok(RunConfig {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
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
    pub fn delete_run_config(&self, config_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM run_configs WHERE id = ?1", params![config_id])?;
        Ok(())
    }
}
