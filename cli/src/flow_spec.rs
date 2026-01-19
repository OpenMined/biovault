use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crate::pipeline_spec::{
    PipelineInputSpec, PipelineSpec, PipelineSqlStoreSpec, PipelineStepSpec, PipelineStoreSpec,
};
use crate::spec_format::FLOW_API_VERSION;

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowFile {
    #[serde(rename = "apiVersion")]
    pub api_version: String,
    pub kind: String,
    pub metadata: FlowMetadata,
    pub spec: FlowSpec,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowMetadata {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inputs: Option<BTreeMap<String, FlowInputSpec>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modules: Option<BTreeMap<String, FlowModuleRef>>,
    pub steps: Vec<FlowStep>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowInputSpec {
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<YamlValue>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowStep {
    pub id: String,
    pub uses: FlowUses,
    #[serde(default)]
    pub with: BTreeMap<String, YamlValue>,
    #[serde(default)]
    pub publish: BTreeMap<String, String>,
    #[serde(default)]
    pub store: BTreeMap<String, FlowStoreSpec>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowUses {
    Name(String),
    Ref(FlowModuleRef),
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FlowModuleRef {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<FlowModuleSource>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FlowModuleSource {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "ref")]
    pub reference: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowStoreSpec {
    pub kind: String,
    pub source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub table: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub key_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub overwrite: Option<bool>,
}

impl FlowFile {
    pub fn parse_yaml(raw: &str) -> Result<Self> {
        let flow: FlowFile = serde_yaml::from_str(raw).context("Failed to parse flow spec")?;
        Ok(flow)
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).with_context(|| {
                    format!("Failed to create flow directory {}", parent.display())
                })?;
            }
        }
        let yaml = serde_yaml::to_string(self).context("Failed to serialize flow spec")?;
        fs::write(path, yaml)
            .with_context(|| format!("Failed to write flow spec to {}", path.display()))?;
        Ok(())
    }

    pub fn to_pipeline_spec(&self) -> Result<PipelineSpec> {
        if self.kind != "Flow" {
            return Err(anyhow!("Expected Flow kind, got '{}'", self.kind));
        }

        let inputs = self
            .spec
            .inputs
            .as_ref()
            .map(|inputs| {
                inputs
                    .iter()
                    .map(|(name, input)| {
                        let spec = match &input.default {
                            Some(value) => PipelineInputSpec::Detailed {
                                raw_type: input.raw_type.clone(),
                                default: Some(stringify_default(value)),
                            },
                            None => PipelineInputSpec::Simple(input.raw_type.clone()),
                        };
                        (name.clone(), spec)
                    })
                    .collect()
            })
            .unwrap_or_default();

        let mut steps = Vec::with_capacity(self.spec.steps.len());
        for step in &self.spec.steps {
            let uses = resolve_step_uses(step, self.spec.modules.as_ref());
            let with = rewrite_bindings(&step.with);
            let publish = step.publish.clone();
            let store = step
                .store
                .iter()
                .map(|(name, spec)| {
                    let store_spec = PipelineStoreSpec::Sql(PipelineSqlStoreSpec {
                        target: None,
                        source: spec.source.clone(),
                        table: spec.table.clone(),
                        participant_column: spec.key_column.clone(),
                        overwrite: spec.overwrite,
                        format: None,
                    });
                    (name.clone(), store_spec)
                })
                .collect();

            steps.push(PipelineStepSpec {
                id: step.id.clone(),
                uses,
                where_exec: None,
                with,
                publish,
                store,
            });
        }

        Ok(PipelineSpec {
            name: self.metadata.name.clone(),
            description: self.metadata.description.clone(),
            context: None,
            inputs,
            steps,
        })
    }

    pub fn from_pipeline_spec(spec: &PipelineSpec) -> Result<FlowFile> {
        let inputs = if spec.inputs.is_empty() {
            None
        } else {
            let mapped = spec
                .inputs
                .iter()
                .map(|(name, input)| {
                    let default = match input {
                        PipelineInputSpec::Detailed { default, .. } => {
                            default.as_ref().map(|value| {
                                serde_yaml::from_str::<YamlValue>(value)
                                    .unwrap_or_else(|_| YamlValue::String(value.clone()))
                            })
                        }
                        _ => None,
                    };
                    (
                        name.clone(),
                        FlowInputSpec {
                            raw_type: input.raw_type().to_string(),
                            default,
                        },
                    )
                })
                .collect();
            Some(mapped)
        };

        let mut modules: BTreeMap<String, FlowModuleRef> = BTreeMap::new();
        let mut steps: Vec<FlowStep> = Vec::with_capacity(spec.steps.len());

        for step in &spec.steps {
            let module_name = step.id.clone();
            let uses_ref = step.uses.clone().unwrap_or_else(|| module_name.clone());
            let mut module_ref = FlowModuleRef {
                source: None,
                version: None,
                digest: None,
            };

            if !uses_ref.is_empty() {
                let mut source = FlowModuleSource {
                    kind: None,
                    path: None,
                    url: None,
                    reference: None,
                };
                if uses_ref.starts_with("http://")
                    || uses_ref.starts_with("https://")
                    || uses_ref.starts_with("syft://")
                {
                    source.kind = Some("url".to_string());
                    source.url = Some(uses_ref.clone());
                } else {
                    source.kind = Some("local".to_string());
                    source.path = Some(uses_ref.clone());
                }
                module_ref.source = Some(source);
            }

            modules.insert(module_name.clone(), module_ref);

            let store = step
                .store
                .iter()
                .map(|(name, spec)| {
                    let store_spec = match spec {
                        PipelineStoreSpec::Sql(sql) => FlowStoreSpec {
                            kind: "sql".to_string(),
                            source: sql.source.clone(),
                            table: sql.table.clone(),
                            key_column: sql.participant_column.clone(),
                            overwrite: sql.overwrite,
                        },
                    };
                    (name.clone(), store_spec)
                })
                .collect();

            steps.push(FlowStep {
                id: step.id.clone(),
                uses: FlowUses::Name(module_name),
                with: step.with.clone(),
                publish: step.publish.clone(),
                store,
            });
        }

        Ok(FlowFile {
            api_version: FLOW_API_VERSION.to_string(),
            kind: "Flow".to_string(),
            metadata: FlowMetadata {
                name: spec.name.clone(),
                description: spec.description.clone(),
                version: None,
            },
            spec: FlowSpec {
                inputs,
                modules: if modules.is_empty() {
                    None
                } else {
                    Some(modules)
                },
                steps,
            },
        })
    }
}

fn resolve_step_uses(
    step: &FlowStep,
    modules: Option<&BTreeMap<String, FlowModuleRef>>,
) -> Option<String> {
    match &step.uses {
        FlowUses::Name(name) => {
            if let Some(modules) = modules {
                if let Some(module_ref) = modules.get(name) {
                    if let Some(source) = &module_ref.source {
                        if let Some(path) = &source.path {
                            return Some(path.clone());
                        }
                        if let Some(url) = &source.url {
                            return Some(url.clone());
                        }
                        if let Some(reference) = &source.reference {
                            return Some(reference.clone());
                        }
                    }
                }
            }
            Some(name.clone())
        }
        FlowUses::Ref(module_ref) => {
            if let Some(source) = &module_ref.source {
                if let Some(path) = &source.path {
                    return Some(path.clone());
                }
                if let Some(url) = &source.url {
                    return Some(url.clone());
                }
                if let Some(reference) = &source.reference {
                    return Some(reference.clone());
                }
            }
            None
        }
    }
}

fn rewrite_bindings(bindings: &BTreeMap<String, YamlValue>) -> BTreeMap<String, YamlValue> {
    bindings
        .iter()
        .map(|(key, value)| (key.clone(), rewrite_value(value)))
        .collect()
}

fn rewrite_value(value: &YamlValue) -> YamlValue {
    match value {
        YamlValue::String(s) => {
            if let Some(stripped) = s.strip_prefix("steps.") {
                YamlValue::String(format!("step.{}", stripped))
            } else {
                YamlValue::String(s.clone())
            }
        }
        YamlValue::Sequence(items) => {
            YamlValue::Sequence(items.iter().map(rewrite_value).collect())
        }
        YamlValue::Mapping(map) => {
            let mut rewritten = serde_yaml::Mapping::new();
            for (key, value) in map {
                let next = rewrite_value(value);
                rewritten.insert(key.clone(), next);
            }
            YamlValue::Mapping(rewritten)
        }
        _ => value.clone(),
    }
}

fn stringify_default(value: &YamlValue) -> String {
    match value {
        YamlValue::String(s) => s.clone(),
        YamlValue::Number(n) => n.to_string(),
        YamlValue::Bool(b) => b.to_string(),
        _ => serde_yaml::to_string(value)
            .unwrap_or_default()
            .trim()
            .to_string(),
    }
}
