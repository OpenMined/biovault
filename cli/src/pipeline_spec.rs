use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

use crate::flow_spec::FlowFile;
use crate::spec_format::{detect_spec_format, SpecFormat};

pub const FLOW_YAML_FILE: &str = "flow.yaml";
pub const PIPELINE_YAML_FILE: &str = "pipeline.yaml";

pub fn resolve_pipeline_spec_path(pipeline_root: &Path) -> std::path::PathBuf {
    let flow_path = pipeline_root.join(FLOW_YAML_FILE);
    if flow_path.exists() {
        return flow_path;
    }
    pipeline_root.join(PIPELINE_YAML_FILE)
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineContextSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub literal: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub from_json: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineStepSpec {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uses: Option<String>,
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub where_exec: Option<String>,
    #[serde(default)]
    pub with: BTreeMap<String, YamlValue>,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub publish: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub store: BTreeMap<String, PipelineStoreSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineSpec {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context: Option<PipelineContextSpec>,
    #[serde(default)]
    pub inputs: BTreeMap<String, PipelineInputSpec>,
    #[serde(default)]
    pub steps: Vec<PipelineStepSpec>,
}

impl PipelineSpec {
    pub fn load(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("Failed to read pipeline spec at {}", path.display()))?;
        match detect_spec_format(path, &raw) {
            SpecFormat::Flow => {
                let flow = FlowFile::from_str(&raw)
                    .with_context(|| format!("Failed to parse flow spec at {}", path.display()))?;
                return flow
                    .to_pipeline_spec()
                    .with_context(|| format!("Failed to convert flow spec at {}", path.display()));
            }
            SpecFormat::Module => {
                return Err(anyhow!(
                    "Detected Module spec at {}. Expected flow.yaml.",
                    path.display()
                ));
            }
            SpecFormat::FlowOverlay => {
                return Err(anyhow!(
                    "Detected FlowOverlay spec at {}. Expected flow.yaml.",
                    path.display()
                ));
            }
            SpecFormat::LegacyProject => {
                return Err(anyhow!(
                    "Detected legacy project spec at {}. Expected flow.yaml.",
                    path.display()
                ));
            }
            SpecFormat::LegacyPipeline | SpecFormat::Unknown => {}
        }
        let spec: PipelineSpec = serde_yaml::from_str(&raw)
            .with_context(|| format!("Failed to parse pipeline spec at {}", path.display()))?;
        Ok(spec)
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).with_context(|| {
                    format!("Failed to create pipeline directory {}", parent.display())
                })?;
            }
        }
        let yaml = serde_yaml::to_string(self).context("Failed to serialize pipeline spec")?;
        fs::write(path, yaml)
            .with_context(|| format!("Failed to write pipeline spec to {}", path.display()))?;
        Ok(())
    }

    pub fn ensure_unique_step_ids(&self) -> Result<()> {
        let mut seen = HashSet::new();
        for step in &self.steps {
            if !seen.insert(&step.id) {
                return Err(anyhow!("Duplicate step id detected: {}", step.id));
            }
        }
        Ok(())
    }
}

pub fn value_to_string(value: &YamlValue) -> Option<String> {
    match value {
        YamlValue::String(s) => Some(s.clone()),
        YamlValue::Number(n) => Some(n.to_string()),
        YamlValue::Bool(b) => Some(b.to_string()),
        YamlValue::Null => None,
        _ => None,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum PipelineInputSpec {
    Simple(String),
    Detailed {
        #[serde(rename = "type")]
        raw_type: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        default: Option<String>,
    },
}

impl PipelineInputSpec {
    pub fn from_type(raw_type: &str) -> Self {
        PipelineInputSpec::Simple(raw_type.to_string())
    }

    pub fn raw_type(&self) -> &str {
        match self {
            PipelineInputSpec::Simple(s) => s,
            PipelineInputSpec::Detailed { raw_type, .. } => raw_type,
        }
    }

    pub fn default_literal(&self) -> Option<&str> {
        match self {
            PipelineInputSpec::Simple(_) => None,
            PipelineInputSpec::Detailed { default, .. } => default.as_deref(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum PipelineStoreSpec {
    #[serde(rename = "sql")]
    Sql(PipelineSqlStoreSpec),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineSqlStoreSpec {
    #[serde(default, alias = "destination")]
    pub target: Option<String>,
    pub source: String,
    #[serde(default, alias = "table_name")]
    pub table: Option<String>,
    #[serde(default, alias = "key_column")]
    #[serde(alias = "participant_id_column")]
    pub participant_column: Option<String>,
    #[serde(default)]
    pub overwrite: Option<bool>,
    #[serde(default)]
    pub format: Option<String>,
}
