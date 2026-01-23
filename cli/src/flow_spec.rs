use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

pub const FLOW_YAML_FILE: &str = "flow.yaml";
pub const FLOW_API_VERSION: &str = "syftbox.openmined.org/v1alpha1";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowFile {
    #[serde(rename = "apiVersion")]
    pub api_version: String,
    pub kind: String,
    pub metadata: FlowMetadata,
    pub spec: FlowFileSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest: Option<FlowManifest>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowManifest {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub modules: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub assets: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub created_at: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowMetadata {
    pub name: String,
    pub version: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub authors: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub labels: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub annotations: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFileSpec {
    /// User-defined variables that can be used in URLs and other templates
    /// Variables can reference other variables and built-in placeholders
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub vars: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inputs: BTreeMap<String, FlowFileInputSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub datasites: Option<FlowDatasitesSpec>,
    /// Top-level coordination setup for progress/state sharing
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordination: Option<FlowCoordinationSpec>,
    /// MPC (Multi-Party Computation) setup for secure channels
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mpc: Option<FlowMpcSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub module_paths: Vec<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub modules: BTreeMap<String, FlowModuleDef>,
    #[serde(default)]
    pub steps: Vec<FlowFileStepSpec>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub outputs: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runtime: Option<FlowRuntimeSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub completion: Option<YamlValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowDatasitesSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub groups: BTreeMap<String, FlowSelector>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowModuleDef {
    Ref(FlowModuleRef),
    Inline(FlowModuleInline),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleInline {
    #[serde(rename = "apiVersion")]
    pub api_version: String,
    pub kind: String,
    pub metadata: ModuleMetadata,
    pub spec: ModuleSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest: Option<ModuleManifest>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleRef {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<FlowModuleSource>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allow_dirty: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interface: Option<FlowModuleInterface>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interface_policy: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assets: Option<Vec<FlowModuleAsset>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub policy: Option<FlowModulePolicy>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trust: Option<FlowModuleTrust>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sandbox: Option<FlowModuleSandbox>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleAsset {
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hash: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleSource {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r#ref: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subpath: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allow_fetch: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleInterface {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inputs: Option<Vec<FlowPortSpec>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub outputs: Option<Vec<FlowPortSpec>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub params: Option<Vec<FlowPortSpec>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowPortSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<FlowFormatSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModulePolicy {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allow_local: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allow_unsigned: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allow_unpinned: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleTrust {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub require_signature: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signature_format: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allowed_signers: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub keyring_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowModuleSandbox {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub enabled: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub network: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filesystem: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allowed_paths: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub env_passthrough: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowStepUses {
    Name(String),
    Ref(FlowModuleRef),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowPermissionRuleSpec {
    pub pattern: String,
    #[serde(default)]
    pub access: FlowPermissionAccessSpec,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowPermissionAccessSpec {
    #[serde(default)]
    pub read: Vec<String>,
    #[serde(default)]
    pub write: Vec<String>,
    #[serde(default)]
    pub admin: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowPermissionSpec {
    /// URL where permissions should be created (must be syft:// scheme)
    pub url: String,
    #[serde(default)]
    pub rules: Vec<FlowPermissionRuleSpec>,
}

/// Coordination spec for setting up progress/state sharing between datasites
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowCoordinationSpec {
    /// Which datasites should create their coordination folders
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub targets: Option<FlowRunTargets>,
    /// Execution strategy (parallel or sequential)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strategy: Option<String>,
    /// Who can read the progress/state files
    /// Use "all" for all datasites, or a list like ["{datasites[*]}"]
    pub share_with: FlowCoordinationShareWith,
    /// Full syft:// URL for coordination folder
    /// Supports variables: {datasite.current}, {flow_name}, {run_id}, {run_path}
    /// Default: {run_path}/_progress
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
}

impl FlowCoordinationSpec {
    /// Get the coordination URL template, defaulting to "{run_path}/_progress"
    pub fn url_template(&self) -> &str {
        self.url.as_deref().unwrap_or("{run_path}/_progress")
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowCoordinationShareWith {
    /// Shorthand: "all" means all datasites
    All(String),
    /// Explicit list of datasites
    List(Vec<String>),
}

impl FlowCoordinationShareWith {
    pub fn is_all(&self) -> bool {
        matches!(self, FlowCoordinationShareWith::All(s) if s == "all")
    }
}

/// MPC (Multi-Party Computation) spec for setting up secure communication channels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowMpcSpec {
    /// Base URL for MPC communication folders
    /// Default: {vars.run_path}/_mpc
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,

    /// Communication topology between parties
    /// - "mesh": everyone <-> everyone (default)
    /// - "star": clients <-> aggregator only (TODO)
    /// - "ring": each party <-> next party (TODO)
    #[serde(default = "FlowMpcSpec::default_topology")]
    pub topology: String,

    /// Folder naming pattern for channels
    /// Default: "{from}_to_{to}"
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pattern: Option<String>,
    // TODO: Custom channel definitions for non-standard topologies
    // #[serde(default, skip_serializing_if = "Vec::is_empty")]
    // pub channels: Vec<FlowMpcChannel>,
}

impl Default for FlowMpcSpec {
    fn default() -> Self {
        Self {
            url: None,
            topology: "mesh".to_string(),
            pattern: None,
        }
    }
}

impl FlowMpcSpec {
    fn default_topology() -> String {
        "mesh".to_string()
    }

    /// Get the MPC base URL template
    pub fn url_template(&self) -> &str {
        self.url.as_deref().unwrap_or("{vars.run_path}/_mpc")
    }

    /// Get the channel folder pattern
    pub fn channel_pattern(&self) -> &str {
        self.pattern.as_deref().unwrap_or("{from}_to_{to}")
    }

    /// Generate channel pairs based on topology
    /// Returns Vec<(from_index, to_index)> for each unidirectional channel
    pub fn generate_channels(&self, party_count: usize) -> Vec<(usize, usize)> {
        match self.topology.as_str() {
            "mesh" => {
                // Everyone to everyone (full mesh)
                let mut channels = Vec::new();
                for from in 0..party_count {
                    for to in 0..party_count {
                        if from != to {
                            channels.push((from, to));
                        }
                    }
                }
                channels
            }
            // TODO: Implement other topologies
            // "star" => {
            //     // All clients (1..n) talk only to aggregator (0)
            //     let mut channels = Vec::new();
            //     for client in 1..party_count {
            //         channels.push((0, client));  // aggregator -> client
            //         channels.push((client, 0));  // client -> aggregator
            //     }
            //     channels
            // }
            // "ring" => {
            //     // Each party talks to next (circular)
            //     let mut channels = Vec::new();
            //     for i in 0..party_count {
            //         let next = (i + 1) % party_count;
            //         channels.push((i, next));
            //         channels.push((next, i));
            //     }
            //     channels
            // }
            _ => {
                // Default to mesh for unknown topologies
                self.generate_channels_mesh(party_count)
            }
        }
    }

    fn generate_channels_mesh(&self, party_count: usize) -> Vec<(usize, usize)> {
        let mut channels = Vec::new();
        for from in 0..party_count {
            for to in 0..party_count {
                if from != to {
                    channels.push((from, to));
                }
            }
        }
        channels
    }
}

/// Parsed conditional input value from `with:` section
/// Supports both simple values and conditional values with only/without
#[derive(Debug, Clone)]
pub struct ConditionalInput {
    pub value: String,
    pub only: Option<Vec<String>>,
    pub without: Option<Vec<String>>,
}

impl ConditionalInput {
    /// Parse a YamlValue into a ConditionalInput
    /// Handles both simple strings and objects with value/only/without
    pub fn from_yaml(yaml: &YamlValue) -> Option<Self> {
        match yaml {
            // Simple string value - applies to all targets
            YamlValue::String(s) => Some(ConditionalInput {
                value: s.clone(),
                only: None,
                without: None,
            }),
            // Object with value/only/without
            YamlValue::Mapping(map) => {
                let value = map
                    .get(YamlValue::String("value".to_string()))
                    .and_then(|v| v.as_str())
                    .map(|s| s.to_string())?;

                let only = map
                    .get(YamlValue::String("only".to_string()))
                    .and_then(|v| match v {
                        YamlValue::String(s) => Some(vec![s.clone()]),
                        YamlValue::Sequence(seq) => {
                            let items: Vec<String> = seq
                                .iter()
                                .filter_map(|item| item.as_str().map(|s| s.to_string()))
                                .collect();
                            if items.is_empty() {
                                None
                            } else {
                                Some(items)
                            }
                        }
                        _ => None,
                    });

                let without =
                    map.get(YamlValue::String("without".to_string()))
                        .and_then(|v| match v {
                            YamlValue::String(s) => Some(vec![s.clone()]),
                            YamlValue::Sequence(seq) => {
                                let items: Vec<String> = seq
                                    .iter()
                                    .filter_map(|item| item.as_str().map(|s| s.to_string()))
                                    .collect();
                                if items.is_empty() {
                                    None
                                } else {
                                    Some(items)
                                }
                            }
                            _ => None,
                        });

                Some(ConditionalInput {
                    value,
                    only,
                    without,
                })
            }
            _ => None,
        }
    }

    /// Check if this input applies to the given target
    /// target_name: the datasite name (e.g., "aggregator@sandbox.local")
    /// groups: map of group names to datasite lists (e.g., "clients" -> ["client1@...", "client2@..."])
    pub fn applies_to(
        &self,
        target_name: &str,
        groups: &std::collections::BTreeMap<String, Vec<String>>,
    ) -> bool {
        // Check "only" constraint
        if let Some(ref only_list) = self.only {
            let matches = only_list.iter().any(|pattern| {
                // Direct match
                if pattern == target_name {
                    return true;
                }
                // Group match
                if let Some(group_members) = groups.get(pattern) {
                    if group_members.contains(&target_name.to_string()) {
                        return true;
                    }
                }
                false
            });
            if !matches {
                return false;
            }
        }

        // Check "without" constraint
        if let Some(ref without_list) = self.without {
            let excluded = without_list.iter().any(|pattern| {
                // Direct match
                if pattern == target_name {
                    return true;
                }
                // Group match
                if let Some(group_members) = groups.get(pattern) {
                    if group_members.contains(&target_name.to_string()) {
                        return true;
                    }
                }
                false
            });
            if excluded {
                return false;
            }
        }

        true
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFileStepSpec {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uses: Option<FlowStepUses>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub depends_on: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run: Option<FlowRunSpec>,
    /// Coordination step: sets up progress/state sharing between datasites
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordination: Option<FlowCoordinationSpec>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub with: BTreeMap<String, YamlValue>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub params: BTreeMap<String, YamlValue>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub publish: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub share: BTreeMap<String, FlowFileShareSpec>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub store: BTreeMap<String, FlowFileStoreSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub permissions: Vec<FlowPermissionSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub when: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub retry: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timeout: Option<YamlValue>,
    /// Barrier: wait for specific step to complete on specified targets before proceeding
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub barrier: Option<FlowBarrierSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowRunSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub targets: Option<FlowRunTargets>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strategy: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub round_robin: Option<FlowRoundRobinSpec>,
}

/// Barrier specification: wait for a step to complete on specified targets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowBarrierSpec {
    /// Step ID to wait for completion
    pub wait_for: String,
    /// Targets to wait for (defaults to "all" if not specified)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub targets: Option<FlowRunTargets>,
    /// Timeout in seconds (defaults to 300 if not specified)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timeout: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowRunTargets {
    One(String),
    Many(Vec<String>),
    Selector(FlowSelector),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowRoundRobinSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub offset: Option<i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowSelector {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFileInputSpec {
    /// Type of the input (e.g., "List[String]"). If not specified, inferred from default value.
    #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
    pub raw_type: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<FlowFormatSpec>,
}

impl FlowFileInputSpec {
    /// Get the raw type, inferring from default if not explicitly specified
    pub fn raw_type(&self) -> String {
        if let Some(ref t) = self.raw_type {
            return t.clone();
        }
        // Infer type from default value
        if let Some(ref default) = self.default {
            return infer_type_from_yaml(default);
        }
        "String".to_string()
    }
}

/// Infer a type string from a YAML value
fn infer_type_from_yaml(value: &YamlValue) -> String {
    match value {
        YamlValue::Null => "String?".to_string(),
        YamlValue::Bool(_) => "Bool".to_string(),
        YamlValue::Number(_) => "Number".to_string(),
        YamlValue::String(_) => "String".to_string(),
        YamlValue::Sequence(seq) => {
            if seq.is_empty() {
                "List[String]".to_string()
            } else {
                let inner = infer_type_from_yaml(&seq[0]);
                format!("List[{}]", inner)
            }
        }
        YamlValue::Mapping(_) => "Map".to_string(),
        YamlValue::Tagged(_) => "String".to_string(),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFormatSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub mapping: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub key_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub delimiter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub header: Option<bool>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub schema: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFileShareSpec {
    pub source: String,
    /// URL where to share the file (must be syft:// scheme)
    pub url: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub permissions: Option<FlowSharePermissions>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[serde(rename = "await")]
    pub r#await: Option<YamlValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowSharePermissions {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub read: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub write: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub admin: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowFileStoreSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    pub source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub table: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub key_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub overwrite: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowRuntimeSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub results_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub work_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resume: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cleanup: Option<YamlValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleMetadata {
    pub name: String,
    pub version: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub authors: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub labels: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub annotations: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleManifest {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub assets: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub created_at: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runner: Option<ModuleRunnerSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<ModulePortSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub outputs: Vec<ModulePortSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parameters: Vec<ModuleParameterSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub assets: Vec<ModuleAssetSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sandbox: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timeout: Option<YamlValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleRunnerSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entrypoint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub template: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub image: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub command: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub env: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModulePortSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<FlowFormatSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub glob: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub regex: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub optional: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cardinality: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleParameterSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleAssetSpec {
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub required: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowContextSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub literal: Option<YamlValue>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub from_json: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum RunsOnSpec {
    One(String),
    Many(Vec<String>),
}

impl RunsOnSpec {
    pub fn as_vec(&self) -> Vec<String> {
        match self {
            RunsOnSpec::One(value) => vec![value.clone()],
            RunsOnSpec::Many(values) => values.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowShareSpec {
    /// Source output to share (e.g., "self.outputs.master_list")
    pub source: String,
    /// URL where to share the file (must be syft:// scheme)
    pub url: String,
    #[serde(default)]
    pub read: Vec<String>,
    #[serde(default)]
    pub write: Vec<String>,
    #[serde(default)]
    pub admin: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowStepSpec {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uses: Option<String>,
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub where_exec: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub foreach: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runs_on: Option<RunsOnSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<String>,
    /// Coordination step: sets up progress/state sharing between datasites
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordination: Option<FlowCoordinationSpec>,
    #[serde(default)]
    pub with: BTreeMap<String, YamlValue>,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub publish: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub share: BTreeMap<String, FlowShareSpec>,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub store: BTreeMap<String, FlowStoreSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub permissions: Vec<FlowPermissionSpec>,
    /// Barrier: wait for specific step to complete on specified targets before proceeding
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub barrier: Option<FlowBarrierSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FlowSpec {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context: Option<FlowContextSpec>,
    /// User-defined variables for URL templates
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub vars: BTreeMap<String, String>,
    /// Top-level coordination setup for progress/state sharing
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordination: Option<FlowCoordinationSpec>,
    /// MPC (Multi-Party Computation) setup for secure channels
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mpc: Option<FlowMpcSpec>,
    #[serde(default)]
    pub inputs: BTreeMap<String, FlowInputSpec>,
    #[serde(default)]
    pub steps: Vec<FlowStepSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub datasites: Vec<String>,
}

impl FlowSpec {
    pub fn load(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("Failed to read flow spec at {}", path.display()))?;
        let flow: FlowFile = serde_yaml::from_str(&raw)
            .with_context(|| format!("Failed to parse flow spec at {}", path.display()))?;
        flow.to_flow_spec()
            .with_context(|| format!("Failed to convert flow spec at {}", path.display()))
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).with_context(|| {
                    format!("Failed to create flow directory {}", parent.display())
                })?;
            }
        }
        let flow = FlowFile::from_flow_spec(self)?;
        let yaml = serde_yaml::to_string(&flow).context("Failed to serialize flow spec")?;
        fs::write(path, yaml)
            .with_context(|| format!("Failed to write flow spec to {}", path.display()))?;
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

impl FlowFile {
    pub fn parse_yaml(raw: &str) -> Result<Self> {
        let flow: FlowFile = serde_yaml::from_str(raw).context("Failed to parse flow spec")?;
        Ok(flow)
    }

    pub fn from_flow_spec(spec: &FlowSpec) -> Result<Self> {
        let metadata = FlowMetadata {
            name: spec.name.clone(),
            version: "0.1.0".to_string(),
            description: spec.description.clone(),
            authors: Vec::new(),
            tags: Vec::new(),
            labels: BTreeMap::new(),
            annotations: BTreeMap::new(),
        };

        let inputs = spec
            .inputs
            .iter()
            .map(|(name, input)| {
                let entry = FlowFileInputSpec {
                    raw_type: Some(input.raw_type().to_string()),
                    description: None,
                    default: input
                        .default_literal()
                        .map(|s| YamlValue::String(s.to_string())),
                    format: None,
                };
                (name.clone(), entry)
            })
            .collect();

        let steps = spec
            .steps
            .iter()
            .map(|step| {
                let with_map = step.with.clone();
                let run = if step.runs_on.is_some()
                    || step.foreach.is_some()
                    || step.order.is_some()
                {
                    Some(FlowRunSpec {
                        targets: step
                            .runs_on
                            .as_ref()
                            .map(|runs| match runs {
                                RunsOnSpec::One(value) => FlowRunTargets::One(value.clone()),
                                RunsOnSpec::Many(values) => FlowRunTargets::Many(values.clone()),
                            })
                            .or_else(|| {
                                step.foreach.as_ref().map(|values| {
                                    if values.len() == 1 {
                                        FlowRunTargets::One(values[0].clone())
                                    } else {
                                        FlowRunTargets::Many(values.clone())
                                    }
                                })
                            }),
                        strategy: step.order.clone(),
                        topology: None,
                        round_robin: None,
                    })
                } else {
                    None
                };

                FlowFileStepSpec {
                    id: step.id.clone(),
                    uses: step.uses.clone().map(FlowStepUses::Name),
                    depends_on: None,
                    run,
                    with: with_map,
                    params: BTreeMap::new(),
                    publish: step.publish.clone(),
                    share: step
                        .share
                        .iter()
                        .map(|(name, share)| {
                            let permissions = FlowSharePermissions {
                                read: share.read.clone(),
                                write: share.write.clone(),
                                admin: share.admin.clone(),
                            };
                            (
                                name.clone(),
                                FlowFileShareSpec {
                                    source: share.source.clone(),
                                    url: share.url.clone(),
                                    permissions: Some(permissions),
                                    r#await: None,
                                },
                            )
                        })
                        .collect(),
                    store: step
                        .store
                        .iter()
                        .map(|(name, store)| {
                            let file_store = match store {
                                FlowStoreSpec::Sql(sql) => FlowFileStoreSpec {
                                    kind: Some("sql".to_string()),
                                    source: sql.source.clone(),
                                    table: sql.table.clone(),
                                    key_column: sql.participant_column.clone(),
                                    overwrite: sql.overwrite,
                                },
                            };
                            (name.clone(), file_store)
                        })
                        .collect(),
                    permissions: step.permissions.clone(),
                    coordination: step.coordination.clone(),
                    when: None,
                    retry: None,
                    timeout: None,
                    barrier: step.barrier.clone(),
                }
            })
            .collect();

        Ok(FlowFile {
            api_version: FLOW_API_VERSION.to_string(),
            kind: "Flow".to_string(),
            metadata,
            spec: FlowFileSpec {
                vars: BTreeMap::new(),
                inputs,
                datasites: None,
                coordination: None,
                mpc: None,
                module_paths: Vec::new(),
                modules: BTreeMap::new(),
                steps,
                outputs: BTreeMap::new(),
                runtime: None,
                completion: None,
            },
            manifest: None,
        })
    }

    pub fn to_flow_spec(&self) -> Result<FlowSpec> {
        let mut inputs = BTreeMap::new();
        for (name, input) in &self.spec.inputs {
            let entry = FlowInputSpec::Detailed {
                raw_type: input.raw_type(),
                default: input.default.as_ref().and_then(value_to_string),
            };
            inputs.insert(name.clone(), entry);
        }

        let datasites_all = resolve_datasites_all(self);
        let groups = resolve_datasite_groups(self, &datasites_all);

        let mut steps = Vec::new();
        for step in &self.spec.steps {
            let uses_value = match step.uses.as_ref() {
                Some(FlowStepUses::Name(name)) => {
                    resolve_step_uses_from_modules(step, &self.spec.modules)
                        .or_else(|| Some(name.clone()))
                }
                Some(FlowStepUses::Ref(module_ref)) => resolve_module_ref(module_ref),
                None => None,
            }
            .or_else(|| resolve_step_uses_from_modules(step, &self.spec.modules));

            let mut with_map = step.with.clone();
            for (key, value) in &step.params {
                with_map.entry(key.clone()).or_insert_with(|| value.clone());
            }
            let (runs_on, foreach, order) =
                run_spec_to_legacy(step.run.as_ref(), &datasites_all, &groups);

            let share = step
                .share
                .iter()
                .map(|(name, share)| {
                    let permissions = share.permissions.clone().unwrap_or_default();
                    (
                        name.clone(),
                        FlowShareSpec {
                            source: share.source.clone(),
                            url: share.url.clone(),
                            read: permissions.read,
                            write: permissions.write,
                            admin: permissions.admin,
                        },
                    )
                })
                .collect();

            let store = step
                .store
                .iter()
                .map(|(name, store)| {
                    let kind = store.kind.as_deref().unwrap_or("sql").to_lowercase();
                    let value = if kind == "sql" {
                        FlowStoreSpec::Sql(FlowSqlStoreSpec {
                            target: None,
                            source: store.source.clone(),
                            table: store.table.clone(),
                            participant_column: store.key_column.clone(),
                            overwrite: store.overwrite,
                            format: None,
                        })
                    } else {
                        FlowStoreSpec::Sql(FlowSqlStoreSpec {
                            target: None,
                            source: store.source.clone(),
                            table: store.table.clone(),
                            participant_column: store.key_column.clone(),
                            overwrite: store.overwrite,
                            format: None,
                        })
                    };
                    (name.clone(), value)
                })
                .collect();

            steps.push(FlowStepSpec {
                id: step.id.clone(),
                uses: uses_value,
                where_exec: None,
                foreach,
                runs_on,
                order,
                coordination: step.coordination.clone(),
                with: with_map,
                publish: step.publish.clone(),
                share,
                store,
                permissions: step.permissions.clone(),
                barrier: step.barrier.clone(),
            });
        }

        Ok(FlowSpec {
            name: self.metadata.name.clone(),
            description: self.metadata.description.clone(),
            context: None,
            vars: self.spec.vars.clone(),
            coordination: self.spec.coordination.clone(),
            mpc: self.spec.mpc.clone(),
            inputs,
            steps,
            datasites: datasites_all,
        })
    }
}

fn run_spec_to_legacy(
    run: Option<&FlowRunSpec>,
    datasites_all: &[String],
    groups: &BTreeMap<String, Vec<String>>,
) -> (Option<RunsOnSpec>, Option<Vec<String>>, Option<String>) {
    let mut runs_on = None;
    let mut foreach = None;
    let mut order = None;

    if let Some(run_spec) = run {
        if let Some(targets) = &run_spec.targets {
            if let Some(resolved) = resolve_targets(targets, datasites_all, groups) {
                runs_on = Some(RunsOnSpec::Many(resolved));
            } else {
                match targets {
                    FlowRunTargets::One(value) => {
                        runs_on = Some(RunsOnSpec::One(value.clone()));
                    }
                    FlowRunTargets::Many(values) => {
                        runs_on = Some(RunsOnSpec::Many(values.clone()));
                    }
                    FlowRunTargets::Selector(selector) => {
                        if !selector.include.is_empty() {
                            runs_on = Some(RunsOnSpec::Many(selector.include.clone()));
                        } else if !selector.exclude.is_empty() {
                            foreach = Some(selector.exclude.clone());
                        }
                    }
                }
            }
        }

        if let Some(strategy) = &run_spec.strategy {
            order = Some(strategy.clone());
        }
    }

    (runs_on, foreach, order)
}

fn resolve_datasites_all(flow: &FlowFile) -> Vec<String> {
    if let Some(datasites) = &flow.spec.datasites {
        if let Some(list) = datasites.all.as_ref().and_then(value_to_string_list) {
            return list;
        }
        if let Some(YamlValue::String(binding)) = datasites.all.as_ref() {
            if binding.trim() == "inputs.datasites" {
                if let Some(input) = flow.spec.inputs.get("datasites") {
                    if let Some(list) = input.default.as_ref().and_then(value_to_string_list) {
                        return list;
                    }
                }
            }
        }
    }
    Vec::new()
}

fn resolve_datasite_groups(
    flow: &FlowFile,
    datasites_all: &[String],
) -> BTreeMap<String, Vec<String>> {
    let mut groups = BTreeMap::new();
    let Some(datasites) = &flow.spec.datasites else {
        return groups;
    };
    for (name, selector) in &datasites.groups {
        let resolved = resolve_selector(selector, datasites_all, &BTreeMap::new());
        if !resolved.is_empty() {
            groups.insert(name.clone(), resolved);
        }
    }
    groups
}

fn resolve_targets(
    targets: &FlowRunTargets,
    datasites_all: &[String],
    groups: &BTreeMap<String, Vec<String>>,
) -> Option<Vec<String>> {
    match targets {
        FlowRunTargets::One(value) => {
            if let Some(group) = groups.get(value) {
                return Some(group.clone());
            }
            if value == "all" && !datasites_all.is_empty() {
                return Some(datasites_all.to_vec());
            }
            None
        }
        FlowRunTargets::Many(values) => {
            let mut resolved = Vec::new();
            for value in values {
                if let Some(group) = groups.get(value) {
                    resolved.extend(group.clone());
                } else {
                    resolved.push(value.clone());
                }
            }
            if resolved.is_empty() {
                None
            } else {
                Some(dedupe_preserve(resolved))
            }
        }
        FlowRunTargets::Selector(selector) => {
            let resolved = resolve_selector(selector, datasites_all, groups);
            if resolved.is_empty() {
                None
            } else {
                Some(resolved)
            }
        }
    }
}

fn resolve_selector(
    selector: &FlowSelector,
    datasites_all: &[String],
    groups: &BTreeMap<String, Vec<String>>,
) -> Vec<String> {
    if let Some(role) = &selector.role {
        if let Some(group) = groups.get(role) {
            return group.clone();
        }
    }

    let mut result = Vec::new();
    for include in &selector.include {
        result.extend(expand_datasite_pattern(include, datasites_all));
    }
    if result.is_empty() {
        result = datasites_all.to_vec();
    }
    if !selector.exclude.is_empty() {
        let mut excluded = Vec::new();
        for entry in &selector.exclude {
            excluded.extend(expand_datasite_pattern(entry, datasites_all));
        }
        result.retain(|item| !excluded.contains(item));
    }

    dedupe_preserve(result)
}

fn expand_datasite_pattern(pattern: &str, datasites_all: &[String]) -> Vec<String> {
    let trimmed = pattern.trim();
    if let Some(inner) = trimmed
        .strip_prefix("{datasites[")
        .and_then(|s| s.strip_suffix("]}"))
    {
        if inner == "*" {
            return datasites_all.to_vec();
        }
        if let Some((start, end)) = inner.split_once(':') {
            let start: usize = start.trim().parse().unwrap_or(0);
            let end: usize = end.trim().parse().unwrap_or(datasites_all.len());
            return datasites_all
                .iter()
                .skip(start)
                .take(end.saturating_sub(start))
                .cloned()
                .collect();
        }
        if let Ok(index) = inner.trim().parse::<usize>() {
            return datasites_all.get(index).cloned().into_iter().collect();
        }
    }
    vec![trimmed.to_string()]
}

fn value_to_string_list(value: &YamlValue) -> Option<Vec<String>> {
    match value {
        YamlValue::Sequence(seq) => {
            let mut out = Vec::new();
            for item in seq {
                if let Some(val) = value_to_string(item) {
                    out.push(val);
                }
            }
            Some(out)
        }
        _ => None,
    }
}

fn dedupe_preserve(values: Vec<String>) -> Vec<String> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for value in values {
        if seen.insert(value.clone()) {
            out.push(value);
        }
    }
    out
}

fn resolve_module_ref(module_ref: &FlowModuleRef) -> Option<String> {
    module_ref
        .source
        .as_ref()
        .and_then(|source| source.path.clone().or(source.url.clone()))
}

fn resolve_step_uses_from_modules(
    step: &FlowFileStepSpec,
    modules: &BTreeMap<String, FlowModuleDef>,
) -> Option<String> {
    let name = match step.uses.as_ref()? {
        FlowStepUses::Name(name) => name.clone(),
        FlowStepUses::Ref(_) => return None,
    };

    let module = modules.get(&name)?;
    match module {
        FlowModuleDef::Ref(reference) => resolve_module_ref(reference),
        FlowModuleDef::Inline(_) => None,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FlowInputSpec {
    Simple(String),
    Detailed {
        #[serde(rename = "type")]
        raw_type: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        default: Option<String>,
    },
}

impl FlowInputSpec {
    pub fn from_type(raw_type: &str) -> Self {
        FlowInputSpec::Simple(raw_type.to_string())
    }

    pub fn raw_type(&self) -> &str {
        match self {
            FlowInputSpec::Simple(s) => s,
            FlowInputSpec::Detailed { raw_type, .. } => raw_type,
        }
    }

    pub fn default_literal(&self) -> Option<&str> {
        match self {
            FlowInputSpec::Simple(_) => None,
            FlowInputSpec::Detailed { default, .. } => default.as_deref(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum FlowStoreSpec {
    #[serde(rename = "sql")]
    Sql(FlowSqlStoreSpec),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowSqlStoreSpec {
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
