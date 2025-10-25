use anyhow::Result;
use rusqlite::{params, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;

use super::BioVaultDb;
use crate::pipeline_spec::PipelineSpec;

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
    pub pipeline_id: Option<i64>,      // NULL for standalone step runs
    pub step_id: Option<i64>,          // NULL for pipeline runs  
    pub status: String,
    pub work_dir: String,
    pub results_dir: Option<String>,
    pub participant_count: Option<i32>, // Only for step runs
    pub created_at: String,
    pub completed_at: Option<String>,
}

// Backwards compat alias
pub type PipelineRun = Run;

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

        // Load spec from pipeline.yaml for each pipeline
        for pipeline in &mut pipelines {
            let yaml_path = PathBuf::from(&pipeline.pipeline_path).join("pipeline.yaml");
            if yaml_path.exists() {
                match fs::read_to_string(&yaml_path) {
                    Ok(content) => {
                        match serde_yaml::from_str::<PipelineSpec>(&content) {
                            Ok(spec) => {
                                pipeline.spec = Some(spec);
                            }
                            Err(e) => {
                                eprintln!("Warning: Failed to parse pipeline.yaml for '{}': {}", pipeline.name, e);
                                eprintln!("  Path: {}", yaml_path.display());
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to read pipeline.yaml for '{}': {}", pipeline.name, e);
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
            let yaml_path = PathBuf::from(&p.pipeline_path).join("pipeline.yaml");
            if yaml_path.exists() {
                match fs::read_to_string(&yaml_path) {
                    Ok(content) => {
                        match serde_yaml::from_str::<PipelineSpec>(&content) {
                            Ok(spec) => {
                                p.spec = Some(spec);
                            }
                            Err(e) => {
                                eprintln!("Warning: Failed to parse pipeline.yaml: {}", e);
                                eprintln!("  Path: {}", yaml_path.display());
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to read pipeline.yaml: {}", e);
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

    /// List all pipeline runs
    pub fn list_pipeline_runs(&self) -> Result<Vec<PipelineRun>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, status, work_dir, results_dir, created_at, completed_at
             FROM pipeline_runs
             ORDER BY created_at DESC",
        )?;

        let runs = stmt
            .query_map([], |row| {
                Ok(PipelineRun {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
                    status: row.get(2)?,
                    work_dir: row.get(3)?,
                    results_dir: row.get(4)?,
                    created_at: row.get(5)?,
                    completed_at: row.get(6)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(runs)
    }

    /// Get a specific pipeline run
    pub fn get_pipeline_run(&self, run_id: i64) -> Result<Option<PipelineRun>> {
        let mut stmt = self.conn.prepare(
            "SELECT id, pipeline_id, status, work_dir, results_dir, created_at, completed_at
             FROM pipeline_runs
             WHERE id = ?1",
        )?;

        let run = stmt
            .query_row([run_id], |row| {
                Ok(PipelineRun {
                    id: row.get(0)?,
                    pipeline_id: row.get(1)?,
                    status: row.get(2)?,
                    work_dir: row.get(3)?,
                    results_dir: row.get(4)?,
                    created_at: row.get(5)?,
                    completed_at: row.get(6)?,
                })
            })
            .optional()?;

        Ok(run)
    }

    /// Create a pipeline run record
    pub fn create_pipeline_run(
        &self,
        pipeline_id: i64,
        work_dir: &str,
        results_dir: Option<&str>,
    ) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO pipeline_runs (pipeline_id, status, work_dir, results_dir, created_at)
             VALUES (?1, 'running', ?2, ?3, datetime('now'))",
            params![pipeline_id, work_dir, results_dir],
        )?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Update pipeline run status
    pub fn update_pipeline_run_status(
        &self,
        run_id: i64,
        status: &str,
        completed: bool,
    ) -> Result<()> {
        if completed {
            self.conn.execute(
                "UPDATE pipeline_runs SET status = ?1, completed_at = datetime('now') WHERE id = ?2",
                params![status, run_id],
            )?;
        } else {
            self.conn.execute(
                "UPDATE pipeline_runs SET status = ?1 WHERE id = ?2",
                params![status, run_id],
            )?;
        }
        Ok(())
    }

    /// Delete a pipeline run
    pub fn delete_pipeline_run(&self, run_id: i64) -> Result<()> {
        self.conn
            .execute("DELETE FROM pipeline_runs WHERE id = ?1", params![run_id])?;
        Ok(())
    }
}

