use crate::module_spec::{
    InputSpec, ModuleAsset, ModuleFile, ModuleFileFormatSpec, ModuleFileInputSpec,
    ModuleFileOutputSpec, ModuleFileParameterSpec, ModuleFileSpec, ModuleRunnerSpec, OutputSpec,
    ParameterSpec,
};
use serde::de::Deserializer;
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftPermissions {
    pub rules: Vec<PermissionRule>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PermissionRule {
    pub pattern: String,
    pub access: AccessControl,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AccessControl {
    pub read: Vec<String>,
    pub write: Vec<String>,
    pub admin: Vec<String>,
}

impl SyftPermissions {
    /// Create permissions with global read access for a single datasite
    pub fn new_for_datasite(datasite_email: &str) -> Self {
        SyftPermissions {
            rules: vec![PermissionRule {
                pattern: "**".to_string(),
                access: AccessControl {
                    read: vec![datasite_email.to_string()],
                    write: vec![],
                    admin: vec![],
                },
            }],
        }
    }

    /// Create permissions with global read access for multiple datasites
    pub fn new_for_datasites(datasites: &[String]) -> Self {
        let mut set = BTreeSet::new();
        for datasite in datasites {
            let trimmed = datasite.trim();
            if !trimmed.is_empty() {
                set.insert(trimmed.to_string());
            }
        }
        let read: Vec<String> = set.into_iter().collect();
        SyftPermissions {
            rules: vec![PermissionRule {
                pattern: "**".to_string(),
                access: AccessControl {
                    read,
                    write: vec![],
                    admin: vec![],
                },
            }],
        }
    }

    /// Add a permission rule for a specific pattern with read/write access for given datasites
    pub fn add_rule(&mut self, pattern: &str, read: Vec<String>, write: Vec<String>) {
        self.rules.push(PermissionRule {
            pattern: pattern.to_string(),
            access: AccessControl {
                read,
                write,
                admin: vec![],
            },
        });
    }

pub fn save(&self, path: &PathBuf) -> anyhow::Result<()> {
        let yaml = serde_yaml::to_string(self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct ModuleYaml {
    pub name: String,
    pub version: String,
    pub author: String,
    pub datasites: Option<Vec<String>>,
    pub participants: Option<Vec<String>>,
    pub workflow: String,
    pub runtime: Option<String>,
    pub assets: Option<Vec<String>>,
    pub parameters: Option<Vec<ParameterSpec>>,
    pub inputs: Option<Vec<InputSpec>>,
    pub outputs: Option<Vec<OutputSpec>>,
    pub b3_hashes: Option<HashMap<String, String>>,
    pub extra: HashMap<String, serde_yaml::Value>,
}

impl ModuleYaml {
    pub fn from_file(path: &PathBuf) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        ModuleYaml::from_str(&content)
    }

    pub fn from_str(content: &str) -> anyhow::Result<Self> {
        let value: serde_yaml::Value = serde_yaml::from_str(content)?;
        let module_file: ModuleFile = serde_yaml::from_value(value.clone())?;
        let datasites = extract_vec(&value, "datasites")?;
        let participants = extract_vec(&value, "participants")?;
        let b3_hashes = extract_map(&value, "b3_hashes")?;
        let extra = extract_extra(&value)?;
        Ok(ModuleYaml::from_module_file(
            module_file,
            datasites,
            participants,
            b3_hashes,
            extra,
        ))
    }

    pub fn save(&self, path: &PathBuf) -> anyhow::Result<()> {
        let module_file = self.to_module_file()?;
        let envelope = ModuleYamlEnvelope {
            module: module_file,
            datasites: self.datasites.as_ref(),
            participants: self.participants.as_ref(),
            b3_hashes: self.b3_hashes.as_ref(),
            extra: &self.extra,
        };
        let yaml = serde_yaml::to_string(&envelope)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }

    fn from_module_file(
        module: ModuleFile,
        datasites: Option<Vec<String>>,
        participants: Option<Vec<String>>,
        b3_hashes: Option<HashMap<String, String>>,
        extra: HashMap<String, serde_yaml::Value>,
    ) -> Self {
        let author = module.metadata.authors.first().cloned().unwrap_or_default();
        let workflow = module
            .spec
            .runner
            .as_ref()
            .and_then(|runner| runner.entrypoint.clone())
            .unwrap_or_else(|| "workflow.nf".to_string());
        let runtime = module
            .spec
            .runner
            .as_ref()
            .and_then(|runner| runner.runtime.clone());
        let assets = module
            .spec
            .assets
            .clone()
            .map(|items| items.into_iter().map(|asset| asset.path).collect());
        let inputs = if module.spec.inputs.is_empty() {
            None
        } else {
            Some(module.spec.inputs.iter().map(to_input_spec).collect())
        };
        let outputs = if module.spec.outputs.is_empty() {
            None
        } else {
            Some(module.spec.outputs.iter().map(to_output_spec).collect())
        };
        let parameters = if module.spec.parameters.is_empty() {
            None
        } else {
            Some(
                module
                    .spec
                    .parameters
                    .iter()
                    .map(to_parameter_spec)
                    .collect(),
            )
        };
        let hashes = b3_hashes.or_else(|| {
            module
                .manifest
                .as_ref()
                .map(|m| m.assets.clone().into_iter().collect())
        });

        ModuleYaml {
            name: module.metadata.name,
            version: module.metadata.version,
            author,
            datasites,
            participants,
            workflow,
            runtime,
            assets,
            parameters,
            inputs,
            outputs,
            b3_hashes: hashes,
            extra,
        }
    }

    fn to_module_file(&self) -> anyhow::Result<ModuleFile> {
        let authors = if self.author.trim().is_empty() {
            Vec::new()
        } else {
            vec![self.author.clone()]
        };
        let inputs = self
            .inputs
            .clone()
            .unwrap_or_default()
            .into_iter()
            .map(to_module_input)
            .collect();
        let outputs = self
            .outputs
            .clone()
            .unwrap_or_default()
            .into_iter()
            .map(to_module_output)
            .collect();
        let parameters = self
            .parameters
            .clone()
            .unwrap_or_default()
            .into_iter()
            .map(to_module_parameter)
            .collect();
        let assets = self.assets.as_ref().map(|paths| {
            paths
                .iter()
                .map(|path| ModuleAsset { path: path.clone() })
                .collect()
        });

        Ok(ModuleFile {
            api_version: crate::module_spec::MODULE_API_VERSION.to_string(),
            kind: "Module".to_string(),
            metadata: crate::module_spec::ModuleFileMetadata {
                name: self.name.clone(),
                version: self.version.clone(),
                description: None,
                authors,
                tags: Vec::new(),
                labels: BTreeMap::new(),
                annotations: BTreeMap::new(),
            },
            spec: ModuleFileSpec {
                runner: Some(ModuleRunnerSpec {
                    kind: None,
                    entrypoint: Some(self.workflow.clone()),
                    runtime: self.runtime.clone(),
                    image: None,
                    command: None,
                    env: BTreeMap::new(),
                    syqure: None,
                }),
                inputs,
                outputs,
                parameters,
                assets,
                sandbox: None,
                timeout: None,
            },
            manifest: self
                .b3_hashes
                .clone()
                .map(|assets| crate::module_spec::ModuleFileManifest {
                    digest: None,
                    assets: assets.into_iter().collect(),
                    created_at: None,
                }),
        })
    }
}

impl Serialize for ModuleYaml {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let module_file = self.to_module_file().map_err(serde::ser::Error::custom)?;
        let envelope = ModuleYamlEnvelope {
            module: module_file,
            datasites: self.datasites.as_ref(),
            participants: self.participants.as_ref(),
            b3_hashes: self.b3_hashes.as_ref(),
            extra: &self.extra,
        };
        envelope.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for ModuleYaml {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = serde_yaml::Value::deserialize(deserializer)?;
        let module_file: ModuleFile =
            serde_yaml::from_value(value.clone()).map_err(serde::de::Error::custom)?;
        let datasites = extract_vec(&value, "datasites").map_err(serde::de::Error::custom)?;
        let participants = extract_vec(&value, "participants").map_err(serde::de::Error::custom)?;
        let b3_hashes = extract_map(&value, "b3_hashes").map_err(serde::de::Error::custom)?;
        let extra = extract_extra(&value).map_err(serde::de::Error::custom)?;
        Ok(ModuleYaml::from_module_file(
            module_file,
            datasites,
            participants,
            b3_hashes,
            extra,
        ))
    }
}

#[derive(Serialize)]
struct ModuleYamlEnvelope<'a> {
    #[serde(flatten)]
    module: ModuleFile,
    #[serde(skip_serializing_if = "Option::is_none")]
    datasites: Option<&'a Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    participants: Option<&'a Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "b3_hashes")]
    b3_hashes: Option<&'a HashMap<String, String>>,
    #[serde(flatten)]
    extra: &'a HashMap<String, serde_yaml::Value>,
}

fn extract_vec(value: &serde_yaml::Value, key: &str) -> anyhow::Result<Option<Vec<String>>> {
    match value.get(key) {
        Some(raw) => Ok(Some(serde_yaml::from_value(raw.clone())?)),
        None => Ok(None),
    }
}

fn extract_map(
    value: &serde_yaml::Value,
    key: &str,
) -> anyhow::Result<Option<HashMap<String, String>>> {
    match value.get(key) {
        Some(raw) => Ok(Some(serde_yaml::from_value(raw.clone())?)),
        None => Ok(None),
    }
}

fn extract_extra(value: &serde_yaml::Value) -> anyhow::Result<HashMap<String, serde_yaml::Value>> {
    let mut extra = HashMap::new();
    if let serde_yaml::Value::Mapping(map) = value {
        for (key, val) in map {
            if let Some(key_str) = key.as_str() {
                if matches!(
                    key_str,
                    "apiVersion"
                        | "kind"
                        | "metadata"
                        | "spec"
                        | "manifest"
                        | "datasites"
                        | "participants"
                        | "b3_hashes"
                ) {
                    continue;
                }
                extra.insert(key_str.to_string(), val.clone());
            }
        }
    }
    Ok(extra)
}

fn to_input_spec(input: &ModuleFileInputSpec) -> InputSpec {
    InputSpec {
        name: input.name.clone(),
        raw_type: input.raw_type.clone(),
        description: input.description.clone(),
        format: input.format.as_ref().and_then(|fmt| fmt.kind.clone()),
        path: None,
        mapping: input.format.as_ref().and_then(|fmt| {
            if fmt.mapping.is_empty() {
                None
            } else {
                Some(fmt.mapping.clone().into_iter().collect())
            }
        }),
    }
}

fn to_output_spec(output: &ModuleFileOutputSpec) -> OutputSpec {
    OutputSpec {
        name: output.name.clone(),
        raw_type: output.raw_type.clone(),
        description: output.description.clone(),
        format: output.format.as_ref().and_then(|fmt| fmt.kind.clone()),
        path: output
            .path
            .clone()
            .or_else(|| output.glob.clone())
            .or_else(|| output.regex.clone()),
    }
}

fn to_parameter_spec(param: &ModuleFileParameterSpec) -> ParameterSpec {
    ParameterSpec {
        name: param.name.clone(),
        raw_type: param.raw_type.clone(),
        description: param.description.clone(),
        default: param.default.clone().map(serde_yaml::Value::String),
        choices: None,
        advanced: None,
    }
}

fn to_module_input(input: InputSpec) -> ModuleFileInputSpec {
    ModuleFileInputSpec {
        name: input.name,
        raw_type: input.raw_type,
        description: input.description,
        format: input.format.map(|kind| ModuleFileFormatSpec {
            kind: Some(kind),
            mapping: input.mapping.unwrap_or_default().into_iter().collect(),
            key_path: None,
            delimiter: None,
            header: None,
            schema: BTreeMap::new(),
        }),
        optional: None,
    }
}

fn to_module_output(output: OutputSpec) -> ModuleFileOutputSpec {
    ModuleFileOutputSpec {
        name: output.name,
        raw_type: output.raw_type,
        description: output.description,
        format: output.format.map(|kind| ModuleFileFormatSpec {
            kind: Some(kind),
            mapping: BTreeMap::new(),
            key_path: None,
            delimiter: None,
            header: None,
            schema: BTreeMap::new(),
        }),
        path: output.path,
        glob: None,
        regex: None,
        cardinality: None,
    }
}

fn to_module_parameter(param: ParameterSpec) -> ModuleFileParameterSpec {
    ModuleFileParameterSpec {
        name: param.name,
        raw_type: param.raw_type,
        default: param.default.as_ref().and_then(|value| match value {
            serde_yaml::Value::String(value) => Some(value.clone()),
            serde_yaml::Value::Number(num) => Some(num.to_string()),
            serde_yaml::Value::Bool(val) => Some(val.to_string()),
            _ => None,
        }),
        description: param.description,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InboxSubmission {
    pub name: String,
    pub author: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub datasites: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub participants: Option<Vec<String>>,
    pub syft_url: String,
    pub status: String,
}

impl InboxSubmission {
    pub fn from_file(path: &PathBuf) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let submission: InboxSubmission = serde_yaml::from_str(&content)?;
        Ok(submission)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn syft_permissions_and_save() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("perm.yaml");
        let perms = SyftPermissions::new_for_datasite("user@example.com");
        assert_eq!(perms.rules.len(), 1);
        perms.save(&p).unwrap();
        let read_back: SyftPermissions =
            serde_yaml::from_str(&fs::read_to_string(&p).unwrap()).unwrap();
        assert_eq!(read_back.rules.len(), 1);
    }

    #[test]
    fn module_yaml_round_trip_and_error() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("module.yaml");
        let proj = ModuleYaml {
            name: "N".into(),
            version: "0.1.0".into(),
            author: "A".into(),
            datasites: Some(vec!["x@y".into()]),
            participants: Some(vec!["P1".into()]),
            workflow: "wf".into(),
            runtime: None,
            assets: Some(vec!["a".into()]),
            parameters: None,
            inputs: None,
            outputs: None,
            b3_hashes: None,
            extra: HashMap::new(),
        };
        proj.save(&p).unwrap();
        let loaded = ModuleYaml::from_file(&p).unwrap();
        assert_eq!(loaded.name, "N");
        assert_eq!(loaded.workflow, "wf");
        assert_eq!(loaded.version, "0.1.0");

        let bad = tmp.path().join("bad.yaml");
        fs::write(&bad, "not: [valid").unwrap();
        assert!(ModuleYaml::from_file(&bad).is_err());
    }

    #[test]
    fn inbox_submission_from_file_and_error() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("inbox.yaml");
        let yaml = r#"
name: X
author: A
datasites: [d@example]
participants: [P]
syft_url: syft://u@example/p
status: queued
"#;
        fs::write(&p, yaml).unwrap();
        let sub = InboxSubmission::from_file(&p).unwrap();
        assert_eq!(sub.name, "X");
        assert_eq!(sub.status, "queued");

        let bad = tmp.path().join("bad.yaml");
        fs::write(&bad, "{").unwrap();
        assert!(InboxSubmission::from_file(&bad).is_err());
    }
}
