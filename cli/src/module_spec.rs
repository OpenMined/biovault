use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::HashMap;

use crate::project_spec::{InputSpec, OutputSpec, ParameterSpec, ProjectSpec};
use crate::spec_format::FLOW_API_VERSION;

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleFile {
    #[serde(rename = "apiVersion")]
    pub api_version: String,
    pub kind: String,
    pub metadata: ModuleMetadata,
    pub spec: ModuleSpec,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleMetadata {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub authors: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub author: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runner: Option<ModuleRunner>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inputs: Option<Vec<ModulePort>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub outputs: Option<Vec<ModulePort>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parameters: Option<Vec<ModuleParameter>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assets: Option<Vec<ModuleAsset>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleRunner {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entrypoint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub template: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModulePort {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<ModuleFormat>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub optional: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping: Option<HashMap<String, String>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleFormat {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping: Option<HashMap<String, String>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleParameter {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<YamlValue>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleAsset {
    pub path: String,
}

impl ModuleFile {
    pub fn parse_yaml(raw: &str) -> Result<Self> {
        let module: ModuleFile =
            serde_yaml::from_str(raw).context("Failed to parse module spec")?;
        Ok(module)
    }

    pub fn to_project_spec(&self) -> Result<ProjectSpec> {
        if self.kind != "Module" {
            return Err(anyhow!("Expected Module kind, got '{}'", self.kind));
        }

        let author = self
            .metadata
            .authors
            .as_ref()
            .and_then(|authors| authors.first().cloned())
            .or_else(|| self.metadata.author.clone())
            .unwrap_or_else(|| "unknown".to_string());

        let workflow = self
            .spec
            .runner
            .as_ref()
            .and_then(|runner| runner.entrypoint.clone())
            .unwrap_or_else(|| "workflow.nf".to_string());

        let template = self
            .spec
            .runner
            .as_ref()
            .and_then(|runner| runner.template.clone());

        let assets = self
            .spec
            .assets
            .as_ref()
            .map(|assets| assets.iter().map(|asset| asset.path.clone()).collect())
            .unwrap_or_default();

        let parameters = self
            .spec
            .parameters
            .as_ref()
            .map(|params| {
                params
                    .iter()
                    .map(|param| ParameterSpec {
                        name: param.name.clone(),
                        raw_type: param.raw_type.clone(),
                        description: param.description.clone(),
                        default: param.default.clone(),
                        choices: None,
                        advanced: None,
                    })
                    .collect()
            })
            .unwrap_or_default();

        let inputs = self
            .spec
            .inputs
            .as_ref()
            .map(|inputs| {
                inputs
                    .iter()
                    .map(|input| InputSpec {
                        name: input.name.clone(),
                        raw_type: normalize_optional_type(&input.raw_type, input.optional),
                        description: input.description.clone(),
                        format: input.format.as_ref().and_then(|format| format.kind.clone()),
                        path: input.path.clone(),
                        mapping: input.mapping.clone().or_else(|| {
                            input
                                .format
                                .as_ref()
                                .and_then(|format| format.mapping.clone())
                        }),
                    })
                    .collect()
            })
            .unwrap_or_default();

        let outputs = self
            .spec
            .outputs
            .as_ref()
            .map(|outputs| {
                outputs
                    .iter()
                    .map(|output| OutputSpec {
                        name: output.name.clone(),
                        raw_type: normalize_optional_type(&output.raw_type, output.optional),
                        description: output.description.clone(),
                        format: output
                            .format
                            .as_ref()
                            .and_then(|format| format.kind.clone()),
                        path: output.path.clone(),
                    })
                    .collect()
            })
            .unwrap_or_default();

        Ok(ProjectSpec {
            name: self.metadata.name.clone(),
            author,
            workflow,
            template,
            version: self.metadata.version.clone(),
            assets,
            parameters,
            inputs,
            outputs,
        })
    }

    pub fn from_project_spec(spec: &ProjectSpec) -> ModuleFile {
        let inputs = spec
            .inputs
            .iter()
            .map(module_port_from_input)
            .collect::<Vec<_>>();
        let outputs = spec
            .outputs
            .iter()
            .map(module_port_from_output)
            .collect::<Vec<_>>();
        let parameters = spec
            .parameters
            .iter()
            .map(|param| ModuleParameter {
                name: param.name.clone(),
                raw_type: param.raw_type.clone(),
                description: param.description.clone(),
                default: param.default.clone(),
            })
            .collect::<Vec<_>>();
        let assets = spec
            .assets
            .iter()
            .map(|path| ModuleAsset { path: path.clone() })
            .collect::<Vec<_>>();

        ModuleFile {
            api_version: FLOW_API_VERSION.to_string(),
            kind: "Module".to_string(),
            metadata: ModuleMetadata {
                name: spec.name.clone(),
                version: spec.version.clone().or_else(|| Some("0.1.0".to_string())),
                description: None,
                authors: Some(vec![spec.author.clone()]),
                author: None,
            },
            spec: ModuleSpec {
                runner: Some(ModuleRunner {
                    kind: Some("nextflow".to_string()),
                    entrypoint: Some(spec.workflow.clone()),
                    template: spec.template.clone(),
                }),
                inputs: Some(inputs),
                outputs: Some(outputs),
                parameters: Some(parameters),
                assets: Some(assets),
            },
        }
    }
}

fn normalize_optional_type(raw_type: &str, optional: Option<bool>) -> String {
    if optional.unwrap_or(false) && !raw_type.ends_with('?') {
        format!("{}?", raw_type)
    } else {
        raw_type.to_string()
    }
}

fn module_port_from_input(input: &InputSpec) -> ModulePort {
    let (raw_type, optional) = split_optional_type(&input.raw_type);
    ModulePort {
        name: input.name.clone(),
        raw_type,
        description: input.description.clone(),
        format: input.format.as_ref().map(|kind| ModuleFormat {
            kind: Some(kind.clone()),
            mapping: input.mapping.clone(),
        }),
        optional,
        path: input.path.clone(),
        mapping: None,
    }
}

fn module_port_from_output(output: &OutputSpec) -> ModulePort {
    let (raw_type, optional) = split_optional_type(&output.raw_type);
    ModulePort {
        name: output.name.clone(),
        raw_type,
        description: output.description.clone(),
        format: output.format.as_ref().map(|kind| ModuleFormat {
            kind: Some(kind.clone()),
            mapping: None,
        }),
        optional,
        path: output.path.clone(),
        mapping: None,
    }
}

fn split_optional_type(raw_type: &str) -> (String, Option<bool>) {
    if let Some(stripped) = raw_type.strip_suffix('?') {
        (stripped.to_string(), Some(true))
    } else {
        (raw_type.to_string(), None)
    }
}
