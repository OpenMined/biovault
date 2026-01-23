<<<<<<< HEAD
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};

use crate::module_spec::ModuleFile;
use crate::spec_format::{detect_spec_format, SpecFormat};

pub const MODULE_YAML_FILE: &str = "module.yaml";
pub const PROJECT_YAML_FILE: &str = "project.yaml";

pub fn resolve_project_spec_path(project_root: &Path) -> std::path::PathBuf {
    let module_path = project_root.join(MODULE_YAML_FILE);
    if module_path.exists() {
        return module_path;
    }
    project_root.join(PROJECT_YAML_FILE)
}

macro_rules! wln {
    ($buf:expr) => {
        writeln!($buf).map_err(|e| anyhow!(e.to_string()))?
    };
    ($buf:expr, $($arg:tt)*) => {
        writeln!($buf, $($arg)*).map_err(|e| anyhow!(e.to_string()))?
    };
}

pub const MODULE_API_VERSION: &str = "syftbox.openmined.org/v1alpha1";

#[derive(Debug, Clone, Serialize, Deserialize)]
=======
use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::Value as YamlValue;
use std::collections::HashMap;

use crate::project_spec::{InputSpec, OutputSpec, ParameterSpec, ProjectSpec};
use crate::spec_format::FLOW_API_VERSION;

#[derive(Debug, Serialize, Deserialize)]
>>>>>>> main
pub struct ModuleFile {
    #[serde(rename = "apiVersion")]
    pub api_version: String,
    pub kind: String,
<<<<<<< HEAD
    pub metadata: ModuleFileMetadata,
    pub spec: ModuleFileSpec,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest: Option<ModuleFileManifest>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileManifest {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub assets: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub created_at: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileMetadata {
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
pub struct ModuleFileSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runner: Option<ModuleRunnerSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<ModuleFileInputSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub outputs: Vec<ModuleFileOutputSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parameters: Vec<ModuleFileParameterSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assets: Option<Vec<ModuleAsset>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sandbox: Option<serde_yaml::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timeout: Option<serde_yaml::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleRunnerSpec {
=======
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
>>>>>>> main
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entrypoint: Option<String>,
<<<<<<< HEAD
    #[serde(default, skip_serializing_if = "Option::is_none", alias = "template")]
    pub runtime: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub image: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub command: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub env: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub syqure: Option<SyqureRunnerOptions>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SyqureRunnerOptions {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub binary: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub docker_image: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub use_docker: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub analyze: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skip_mhe_setup: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transport: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub poll_ms: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub platform: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileInputSpec {
=======
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub template: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModulePort {
>>>>>>> main
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
<<<<<<< HEAD
    pub format: Option<ModuleFileFormatSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub optional: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileOutputSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<ModuleFileFormatSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub glob: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub regex: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cardinality: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileParameterSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModuleFileFormatSpec {
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
pub struct ModuleAsset {
    pub path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModuleSpec {
    pub name: String,
    pub author: String,
    pub workflow: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none", alias = "template")]
    pub runtime: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub datasites: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub env: BTreeMap<String, String>,
    #[serde(default)]
    pub assets: Vec<String>,
    #[serde(default)]
    pub parameters: Vec<ParameterSpec>,
    #[serde(default)]
    pub inputs: Vec<InputSpec>,
    #[serde(default)]
    pub outputs: Vec<OutputSpec>,
    #[serde(default)]
    pub steps: Vec<ModuleStepSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runner: Option<ModuleRunnerSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<serde_yaml::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub choices: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub advanced: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,
=======
    pub format: Option<ModuleFormat>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub optional: Option<bool>,
>>>>>>> main
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping: Option<HashMap<String, String>>,
}

<<<<<<< HEAD
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputSpec {
=======
#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleFormat {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping: Option<HashMap<String, String>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleParameter {
>>>>>>> main
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
<<<<<<< HEAD
    pub format: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModuleStepSpec {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub foreach: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub env: BTreeMap<String, String>,
    #[serde(default)]
    pub inputs: Vec<InputSpec>,
    #[serde(default)]
    pub outputs: Vec<OutputSpec>,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
enum ParameterType {
    String,
    Bool,
    Enum(Vec<String>),
}

#[derive(Debug, Clone)]
struct RecordField {
    name: String,
    ty: TypeExpr,
}

#[derive(Debug, Clone)]
enum TypeExpr {
    String,
    Bool,
    File,
    Directory,
    ParticipantSheet,
    GenotypeRecord,
    BiovaultContext,
    #[allow(dead_code)]
    List(Box<TypeExpr>),
    #[allow(dead_code)]
    Map(Box<TypeExpr>),
    #[allow(dead_code)]
    Record(Vec<RecordField>),
    Optional(Box<TypeExpr>),
}

impl ModuleSpec {
    pub fn load(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
<<<<<<<< HEAD:cli/src/module_spec.rs
            .with_context(|| format!("Failed to read module spec at {}", path.display()))?;
        let module = ModuleFile::parse_yaml(&raw)
            .with_context(|| format!("Failed to parse module spec at {}", path.display()))?;
        module
            .to_module_spec()
            .with_context(|| format!("Failed to convert module spec at {}", path.display()))
========
            .with_context(|| format!("Failed to read project spec at {}", path.display()))?;
        Self::load_from_str(path, &raw)
    }

    /// Load from raw bytes (useful when reading from decrypted storage)
    pub fn load_from_bytes(path: &Path, bytes: &[u8]) -> Result<Self> {
        let raw = std::str::from_utf8(bytes)
            .with_context(|| format!("Failed to read project spec at {}", path.display()))?;
        Self::load_from_str(path, raw)
    }

    fn load_from_str(path: &Path, raw: &str) -> Result<Self> {
        match detect_spec_format(path, raw) {
            SpecFormat::Flow => {
                return Err(anyhow!(
                    "Detected Flow spec at {}. ProjectSpec loader does not support Flow.",
                    path.display()
                ));
            }
            SpecFormat::Module => {
                let module = ModuleFile::parse_yaml(raw).with_context(|| {
                    format!("Failed to parse module spec at {}", path.display())
                })?;
                return module.to_project_spec().with_context(|| {
                    format!("Failed to convert module spec at {}", path.display())
                });
            }
            SpecFormat::FlowOverlay => {
                return Err(anyhow!(
                    "Detected FlowOverlay spec at {}. ProjectSpec loader expects module.yaml.",
                    path.display()
                ));
            }
            SpecFormat::LegacyPipeline => {
                return Err(anyhow!(
                    "Detected legacy pipeline spec at {}. Expected module.yaml.",
                    path.display()
                ));
            }
            SpecFormat::LegacyProject | SpecFormat::Unknown => {}
        }
        let spec: ProjectSpec = serde_yaml::from_str(raw)
            .with_context(|| format!("Failed to parse project spec at {}", path.display()))?;
        Ok(spec)
>>>>>>>> main:cli/src/project_spec.rs
    }
}

pub fn scaffold_from_spec(mut spec: ModuleSpec, target_dir: &Path) -> Result<ModuleSpec> {
    if target_dir.exists() {
        if target_dir.is_file() {
            bail!("Target {} exists and is a file", target_dir.display());
        }
        if target_dir.read_dir()?.next().is_some() {
            bail!("Target directory {} is not empty", target_dir.display());
        }
    } else {
        fs::create_dir_all(target_dir)
            .with_context(|| format!("Failed to create directory {}", target_dir.display()))?;
    }

    if spec.name.trim().is_empty() {
        bail!("Spec 'name' field cannot be empty");
    }
    if spec.author.trim().is_empty() {
        bail!("Spec 'author' field cannot be empty");
    }
    if spec.workflow.trim().is_empty() {
        bail!("Spec 'workflow' field cannot be empty");
    }

    let runtime_name = spec
        .runtime
        .clone()
        .unwrap_or_else(|| "nextflow".to_string());
    spec.runtime = Some(runtime_name);

<<<<<<<< HEAD:cli/src/module_spec.rs
    let module_yaml_path = target_dir.join("module.yaml");
    let workflow_path = target_dir.join(&spec.workflow);
    let assets_dir = target_dir.join("assets");

    let module = ModuleFile::from_module_spec(&spec);
    let yaml = serde_yaml::to_string(&module).context("Failed to serialize module spec")?;
    fs::write(&module_yaml_path, yaml).context("Failed to write module.yaml")?;
========
    let project_yaml_path = target_dir.join(MODULE_YAML_FILE);
    let workflow_path = target_dir.join(&spec.workflow);
    let assets_dir = target_dir.join("assets");

    let module = ModuleFile::from_project_spec(&spec);
    let yaml = serde_yaml::to_string(&module).context("Failed to serialize module spec")?;
    fs::write(&project_yaml_path, yaml).context("Failed to write module.yaml")?;
>>>>>>>> main:cli/src/project_spec.rs

    if let Some(parent) = workflow_path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create directory {}", parent.display()))?;
    }
    let workflow_contents = generate_workflow_stub(&spec)?;
    fs::write(&workflow_path, workflow_contents).context("Failed to write workflow stub")?;

    // Note: template.nf is NOT written to module folder - it's a security boundary
    // Templates are installed in ~/.biovault/env/{template_name}/ and loaded at runtime

    fs::create_dir_all(&assets_dir).context("Failed to create assets directory")?;
    for asset in &spec.assets {
        if asset.trim().is_empty() {
            continue;
        }
        let asset_path = assets_dir.join(asset);
        if let Some(parent) = asset_path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create directory {}", parent.display()))?;
        }
        if !asset_path.exists() {
            fs::write(&asset_path, format!("# TODO: populate {}\n", asset)).with_context(|| {
                format!(
                    "Failed to create asset placeholder {}",
                    asset_path.display()
                )
            })?;
        }
    }

    Ok(spec)
}

#[allow(dead_code)]
fn parse_parameter_type(spec: &ParameterSpec) -> Result<ParameterType> {
    match spec.raw_type.as_str() {
        "String" => Ok(ParameterType::String),
        "Bool" => Ok(ParameterType::Bool),
        other => {
            if let Some(inner) = other.strip_prefix("Enum[") {
                let inner = inner.strip_suffix(']').ok_or_else(|| {
                    anyhow!("Invalid Enum parameter type '{}': missing closing ]", other)
                })?;
                let choices: Vec<String> = inner
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
                if choices.is_empty() {
                    bail!(
                        "Enum parameter '{}' must declare at least one choice",
                        spec.name
                    );
                }
                Ok(ParameterType::Enum(choices))
            } else {
                bail!(
                    "Unsupported parameter type '{}' for parameter '{}'. Supported: String, Bool, Enum[...]",
                    spec.raw_type,
                    spec.name
                )
            }
        }
    }
}

fn split_top_level(input: &str, delimiter: char) -> Vec<String> {
    let mut parts = Vec::new();
    let mut depth = 0usize;
    let mut start = 0usize;

    for (idx, ch) in input.char_indices() {
        match ch {
            '[' | '{' => depth = depth.saturating_add(1),
            ']' | '}' => depth = depth.saturating_sub(1),
            _ => {}
        }

        if ch == delimiter && depth == 0 {
            parts.push(input[start..idx].trim().to_string());
            start = idx + ch.len_utf8();
        }
    }

    parts.push(input[start..].trim().to_string());
    parts
}

fn split_top_level_once(input: &str, delimiter: char) -> Option<(String, String)> {
    let mut depth = 0usize;
    for (idx, ch) in input.char_indices() {
        match ch {
            '[' | '{' => depth = depth.saturating_add(1),
            ']' | '}' => depth = depth.saturating_sub(1),
            _ => {}
        }

        if ch == delimiter && depth == 0 {
            let left = input[..idx].trim().to_string();
            let right = input[idx + ch.len_utf8()..].trim().to_string();
            return Some((left, right));
        }
    }
    None
}

fn strip_wrapped<'a>(raw: &'a str, prefix: &str, suffix: char) -> Option<&'a str> {
    if raw.len() < prefix.len() + 1 {
        return None;
    }
    if !raw[..prefix.len()].eq_ignore_ascii_case(prefix) {
        return None;
    }
    if !raw.ends_with(suffix) {
        return None;
    }
    Some(raw[prefix.len()..raw.len() - 1].trim())
}

fn parse_type_expr(raw: &str) -> Result<TypeExpr> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        bail!("Type expression cannot be empty");
    }

    let (base, optional) = if let Some(stripped) = trimmed.strip_suffix('?') {
        (stripped.trim(), true)
    } else {
        (trimmed, false)
    };

    let parsed = if let Some(inner) = strip_wrapped(base, "List[", ']') {
        TypeExpr::List(Box::new(parse_type_expr(inner)?))
    } else if let Some(inner) = strip_wrapped(base, "Map[", ']') {
        let parts = split_top_level(inner, ',');
        if parts.len() != 2 {
            bail!("Invalid Map type '{}': expected Map[String, T]", raw);
        }
        if !parts[0].eq_ignore_ascii_case("String") {
            bail!("Invalid Map type '{}': key type must be String", raw);
        }
        TypeExpr::Map(Box::new(parse_type_expr(&parts[1])?))
    } else if let Some(inner) =
        strip_wrapped(base, "Record{", '}').or_else(|| strip_wrapped(base, "Dict{", '}'))
    {
        if inner.is_empty() {
            bail!("Record type '{}' must declare at least one field", raw);
        }
        let mut seen = std::collections::HashSet::new();
        let mut fields = Vec::new();
        for field in split_top_level(inner, ',') {
            let (name, ty_raw) = split_top_level_once(&field, ':')
                .ok_or_else(|| anyhow!("Invalid Record field '{}': expected name: Type", field))?;
            if name.is_empty() {
                bail!("Invalid Record field '{}': missing field name", field);
            }
            if !seen.insert(name.clone()) {
                bail!("Duplicate Record field '{}'", name);
            }
            fields.push(RecordField {
                name,
                ty: parse_type_expr(&ty_raw)?,
            });
        }
        TypeExpr::Record(fields)
    } else {
        let base_lower = base.to_ascii_lowercase();
        match base_lower.as_str() {
            "string" => TypeExpr::String,
            "bool" => TypeExpr::Bool,
            "file" => TypeExpr::File,
            "directory" => TypeExpr::Directory,
            "participantsheet" => TypeExpr::ParticipantSheet,
            "genotyperecord" => TypeExpr::GenotypeRecord,
            "biovaultcontext" => TypeExpr::BiovaultContext,
            other => bail!(
                "Unsupported type '{}' (supported primitives: String, Bool, File, Directory, ParticipantSheet, GenotypeRecord, BiovaultContext, List[..], Map[String, ..], Record{{..}})",
                other
            ),
        }
    };

    if optional {
        Ok(TypeExpr::Optional(Box::new(parsed)))
    } else {
        Ok(parsed)
    }
}

pub fn validate_type_expr(raw: &str) -> Result<()> {
    parse_type_expr(raw).map(|_| ())
}

pub fn types_compatible(expected: &str, actual: &str) -> bool {
    match (parse_type_expr(expected), parse_type_expr(actual)) {
        (Ok(expected_ty), Ok(actual_ty)) => type_exprs_compatible(&expected_ty, &actual_ty),
        _ => normalize_type_fallback(expected) == normalize_type_fallback(actual),
    }
}

fn normalize_type_fallback(value: &str) -> String {
    value.trim().trim_end_matches('?').to_ascii_lowercase()
}

fn type_exprs_compatible(expected: &TypeExpr, actual: &TypeExpr) -> bool {
    match (expected, actual) {
        (TypeExpr::Optional(inner), _) => type_exprs_compatible(inner, actual),
        (_, TypeExpr::Optional(inner)) => type_exprs_compatible(expected, inner),
        (TypeExpr::String, TypeExpr::String)
        | (TypeExpr::Bool, TypeExpr::Bool)
        | (TypeExpr::File, TypeExpr::File)
        | (TypeExpr::Directory, TypeExpr::Directory)
        | (TypeExpr::ParticipantSheet, TypeExpr::ParticipantSheet)
        | (TypeExpr::GenotypeRecord, TypeExpr::GenotypeRecord)
        | (TypeExpr::BiovaultContext, TypeExpr::BiovaultContext) => true,
        (TypeExpr::List(expected_inner), TypeExpr::List(actual_inner)) => {
            type_exprs_compatible(expected_inner, actual_inner)
        }
        (TypeExpr::Map(expected_inner), TypeExpr::Map(actual_inner)) => {
            type_exprs_compatible(expected_inner, actual_inner)
        }
        (TypeExpr::Record(expected_fields), TypeExpr::Record(actual_fields)) => {
            if expected_fields.len() != actual_fields.len() {
                return false;
            }
            for expected_field in expected_fields {
                let Some(actual_field) = actual_fields
                    .iter()
                    .find(|field| field.name == expected_field.name)
                else {
                    return false;
                };
                if !type_exprs_compatible(&expected_field.ty, &actual_field.ty) {
                    return false;
                }
            }
            true
        }
        _ => false,
    }
}

pub fn generate_template_nf(spec: &ModuleSpec) -> Result<String> {
    let mut buf = String::new();

    wln!(&mut buf, "nextflow.enable.dsl=2");
    wln!(&mut buf);

    wln!(
        &mut buf,
        "params.results_dir = params.results_dir ?: 'results'"
    );
    wln!(
        &mut buf,
        "params.work_flow_file = params.work_flow_file ?: 'workflow.nf'"
    );
    wln!(&mut buf);

    let mut param_names = Vec::new();
    for parameter in &spec.parameters {
        let kind = parse_parameter_type(parameter)?;
        let default_line = match (&kind, &parameter.default) {
            (ParameterType::String, Some(value)) => value.as_str().map(|s| {
                format!(
                    "params.{name} = params.{name} ?: '{value}'",
                    name = parameter.name,
                    value = s.replace('\'', "\\'")
                )
            }),
            (ParameterType::Bool, Some(value)) => value.as_bool().map(|boolean| {
                format!(
                    "params.{name} = (params.{name} == null ? {val} : params.{name})",
                    name = parameter.name,
                    val = boolean
                )
            }),
            (ParameterType::Enum(choices), Some(value)) => {
                if let Some(s) = value.as_str() {
                    if !choices.is_empty() && !choices.contains(&s.to_string()) {
                        return Err(anyhow!(
                            "Default '{}' for parameter '{}' is not one of {:?}",
                            s,
                            parameter.name,
                            choices
                        ));
                    }
                    Some(format!(
                        "params.{name} = params.{name} ?: '{value}'",
                        name = parameter.name,
                        value = s.replace('\'', "\\'")
                    ))
                } else {
                    None
                }
            }
            _ => None,
        };
        if let Some(line) = default_line {
            wln!(&mut buf, "{}", line);
        }
        param_names.push(parameter.name.clone());
    }

    wln!(&mut buf);
    if param_names.is_empty() {
        wln!(&mut buf, "def __bv_params = [:]");
    } else {
        wln!(&mut buf, "def __bv_params = [");
        for (idx, name) in param_names.iter().enumerate() {
            let suffix = if idx + 1 == param_names.len() {
                ""
            } else {
                ","
            };
            wln!(&mut buf, "    {}: params.{}{}", name, name, suffix);
        }
        wln!(&mut buf, "]");
    }
    wln!(&mut buf);

    wln!(&mut buf, "def __bv_deep_freeze(value) {{");
    wln!(&mut buf, "    if (value instanceof Map) {{");
    wln!(
        &mut buf,
        "        return value.collectEntries {{ k, v -> [k, __bv_deep_freeze(v)] }}.asImmutable()"
    );
    wln!(&mut buf, "    }}");
    wln!(&mut buf, "    if (value instanceof Collection) {{");
    wln!(
        &mut buf,
        "        return value.collect {{ __bv_deep_freeze(it) }}.asImmutable()"
    );
    wln!(&mut buf, "    }}");
    wln!(&mut buf, "    return value");
    wln!(&mut buf, "}}");
    wln!(&mut buf);

    wln!(&mut buf, "def __bv_context_base = [");
    wln!(&mut buf, "    run_id      : params.run_id,");
    wln!(&mut buf, "    datasite    : params.datasite,");
    wln!(&mut buf, "    user        : params.user,");
    wln!(&mut buf, "    run_timestamp: params.run_timestamp,");
    wln!(&mut buf, "    params      : __bv_params");
    wln!(&mut buf, "]");
    wln!(&mut buf, "    .findAll {{ it.value != null }}");
    wln!(
        &mut buf,
        "def context = __bv_deep_freeze(__bv_context_base)"
    );
    wln!(&mut buf);

    let mut argument_names = Vec::new();
    for input in &spec.inputs {
        let ty = parse_type_expr(&input.raw_type)?;
        if is_context_type(&ty) {
            continue;
        }
        wln!(
            &mut buf,
            "if (!params.containsKey('{}')) params.{} = null",
            input.name,
            input.name
        );
        let binding_line = generate_input_binding(&input.name, &ty);
        wln!(&mut buf, "{}", binding_line);
        if let Some(comment) = generate_input_comment(&ty) {
            wln!(&mut buf, "{}", comment);
        }
        wln!(&mut buf);
        argument_names.push(input.name.clone());
    }

    wln!(
        &mut buf,
        "include {{ USER }} from \"${{params.work_flow_file}}\""
    );
    wln!(&mut buf);

    wln!(&mut buf, "workflow {{");
    wln!(&mut buf, "    file(params.results_dir).mkdirs()");
    wln!(&mut buf);
    wln!(&mut buf, "    USER(");
    let mut call_args = vec!["context".to_string()];
    call_args.extend(argument_names.iter().cloned());
    for (idx, arg) in call_args.iter().enumerate() {
        let suffix = if idx + 1 == call_args.len() { "" } else { "," };
        wln!(&mut buf, "        {}{}", arg, suffix);
    }
    wln!(&mut buf, "    )");
    wln!(&mut buf, "}}");

    Ok(buf)
}

pub fn generate_workflow_stub(spec: &ModuleSpec) -> Result<String> {
    let mut buf = String::new();

    wln!(&mut buf, "nextflow.enable.dsl=2");
    wln!(&mut buf);

    wln!(&mut buf, "workflow USER {{");
    wln!(&mut buf, "    take:");
    wln!(&mut buf, "        context");
    let mut preview_inputs = Vec::new();
    for input in &spec.inputs {
        let ty = parse_type_expr(&input.raw_type)?;
        if is_context_type(&ty) {
            continue;
        }
        wln!(&mut buf, "        {}", input.name);
        preview_inputs.push((input.name.clone(), ty));
    }
    wln!(&mut buf);

    wln!(&mut buf, "    main:");
    wln!(
        &mut buf,
        "        println \"[bv] context params: ${{context.params}}\""
    );
    wln!(
        &mut buf,
        "        println \"[bv] context keys: ${{context.keySet()}}\""
    );

    // If there are Python assets, generate example usage
    let python_assets: Vec<&String> = spec.assets.iter().filter(|a| a.ends_with(".py")).collect();
    if !python_assets.is_empty() && !preview_inputs.is_empty() && !spec.outputs.is_empty() {
        let script_name = python_assets[0];
        let first_input = &preview_inputs[0].0;
        let first_output = &spec.outputs[0].name;

        wln!(&mut buf);
        wln!(&mut buf, "        // Example: Use your Python script");
        wln!(
            &mut buf,
            "        def assetsDir = file(context.params.assets_dir)"
        );
        wln!(
            &mut buf,
            "        def script = file(\"${{assetsDir}}/{}\") ",
            script_name
        );
        wln!(
            &mut buf,
            "        def result_ch = process_data({}, Channel.value(script))",
            first_input
        );
        wln!(&mut buf);
        wln!(&mut buf, "    emit:");
        wln!(&mut buf, "        {} = result_ch", first_output);
        for output in spec.outputs.iter().skip(1) {
            wln!(
                &mut buf,
                "        {} = null // TODO: wire this output",
                output.name
            );
        }
        wln!(&mut buf, "}}");
        wln!(&mut buf);
        wln!(&mut buf, "process process_data {{");
        wln!(&mut buf, "    publishDir params.results_dir, mode: 'copy'");
        wln!(&mut buf);
        wln!(&mut buf, "    input:");
        wln!(&mut buf, "        val input_data");
        wln!(&mut buf, "        path script");
        wln!(&mut buf);
        wln!(&mut buf, "    output:");
        wln!(&mut buf, "        path 'output.txt'");
        wln!(&mut buf);
        wln!(&mut buf, "    script:");
        wln!(&mut buf, "    \"\"\"");
        wln!(&mut buf, "    python3 ${{script}} \\");
        wln!(&mut buf, "        --input \"${{input_data}}\" \\");
        wln!(&mut buf, "        --output output.txt");
        wln!(&mut buf, "    \"\"\"");
        wln!(&mut buf, "}}");
    } else {
        // Original preview-only code
        for (name, ty) in &preview_inputs {
            for line in generate_preview_lines(name, ty, "        ") {
                wln!(&mut buf, "{}", line);
            }
        }
        if spec.outputs.is_empty() {
            wln!(
                &mut buf,
                "        // Emit outputs using 'emit' block when ready"
            );
        } else {
            wln!(&mut buf, "        // Emit declared outputs once produced");
            for output in &spec.outputs {
                wln!(
                    &mut buf,
                    "        // - {} ({})",
                    output.name,
                    output.raw_type
                );
            }
        }
        wln!(&mut buf, "}}");
    }

    Ok(buf)
}

fn generate_preview_lines(name: &str, ty: &TypeExpr, indent: &str) -> Vec<String> {
    match ty {
        TypeExpr::String | TypeExpr::Bool => {
            vec![format!(r#"{indent}println "[bv] {name}: ${{{name}}}""#)]
        }
        TypeExpr::File | TypeExpr::Directory => vec![format!(
            r#"{indent}{name}.view {{ println "[bv] {name}: $it" }}"#
        )],
        TypeExpr::List(_) | TypeExpr::GenotypeRecord => vec![
            format!(r#"{indent}if ({name} instanceof Channel) {{"#),
            format!(r#"{indent}    {name}.view {{ println "[bv] {name} item -> $it" }}"#),
            format!(r#"{indent}}} else if ({name} != null) {{"#),
            format!(
                r#"{indent}    ({name} as List).take(5).each {{ println "[bv] {name} item -> $it" }}"#
            ),
            format!(r#"{indent}}} else {{"#),
            format!(r#"{indent}    println "[bv] {name}: (empty)""#),
            format!(r#"{indent}}}"#),
        ],
        TypeExpr::ParticipantSheet => vec![
            format!(r#"{indent}def __preview_{name} = {name}"#),
            format!(r#"{indent}if (__preview_{name}?.rows instanceof Channel) {{"#),
            format!(
                r#"{indent}    (__preview_{name}.rows as Channel).view {{ println "[bv] {name} row -> $it" }}"#
            ),
            format!(r#"{indent}}} else if (__preview_{name}?.rows) {{"#),
            format!(
                r#"{indent}    (__preview_{name}.rows as List).take(5).each {{ println "[bv] {name} row -> $it" }}"#
            ),
            format!(r#"{indent}}} else if (__preview_{name} instanceof Channel) {{"#),
            format!(r#"{indent}    __preview_{name}.view {{ println "[bv] {name} -> $it" }}"#),
            format!(r#"{indent}}} else if (__preview_{name} != null) {{"#),
            format!(r#"{indent}    println "[bv] {name}: ${{__preview_{name}}}""#),
            format!(r#"{indent}}} else {{"#),
            format!(r#"{indent}    println "[bv] {name}: (empty)""#),
            format!(r#"{indent}}}"#),
        ],
        TypeExpr::Map(_) | TypeExpr::Record(_) => {
            vec![format!(r#"{indent}println "[bv] {name}: ${{{name}}}""#)]
        }
        TypeExpr::Optional(inner) => {
            let mut lines = Vec::new();
            lines.push(format!(r#"{indent}if ({name} != null) {{"#));
            for line in generate_preview_lines(name, inner, &(indent.to_owned() + "    ")) {
                lines.push(line);
            }
            lines.push(format!(r#"{indent}}} else {{"#));
            lines.push(format!(r#"{indent}    println "[bv] {name}: (null)""#));
            lines.push(format!(r#"{indent}}}"#));
            lines
        }
        TypeExpr::BiovaultContext => Vec::new(),
    }
}

#[allow(dead_code)]
fn generate_input_binding(name: &str, ty: &TypeExpr) -> String {
    match ty {
        TypeExpr::String => format!("def {name} = params.{name}"),
        TypeExpr::Bool => format!("def {name} = params.{name}"),
        TypeExpr::File => format!("def {name} = Channel.fromPath(params.{name})"),
        TypeExpr::Directory => format!("def {name} = Channel.fromPath(params.{name})"),
        TypeExpr::BiovaultContext => "def context = context".to_string(),
        TypeExpr::ParticipantSheet
        | TypeExpr::GenotypeRecord
        | TypeExpr::List(_)
        | TypeExpr::Map(_)
        | TypeExpr::Record(_) => format!("def {name} = params.{name}"),
        TypeExpr::Optional(inner) => generate_input_binding(name, inner),
    }
}

#[allow(dead_code)]
fn generate_input_comment(ty: &TypeExpr) -> Option<String> {
    match ty {
        TypeExpr::ParticipantSheet => {
            Some("// TODO: Replace with sheet loader for ParticipantSheet".to_string())
        }
        TypeExpr::GenotypeRecord => {
            Some("// TODO: Map params value to GenotypeRecord channel".to_string())
        }
        TypeExpr::List(inner) => Some(format!(
            "// TODO: Populate List[{}] channel",
            describe_type(inner)
        )),
        TypeExpr::Map(inner) => Some(format!(
            "// TODO: Populate Map[String, {}] value",
            describe_type(inner)
        )),
        TypeExpr::Record(fields) => {
            let field_list = fields
                .iter()
                .map(|field| format!("{}: {}", field.name, describe_type(&field.ty)))
                .collect::<Vec<_>>()
                .join(", ");
            Some(format!("// TODO: Populate Record{{{}}} value", field_list))
        }
        TypeExpr::Optional(inner) => generate_input_comment(inner),
        _ => None,
    }
}

#[allow(dead_code)]
fn describe_type(ty: &TypeExpr) -> String {
    match ty {
        TypeExpr::String => "String".to_string(),
        TypeExpr::Bool => "Bool".to_string(),
        TypeExpr::File => "File".to_string(),
        TypeExpr::Directory => "Directory".to_string(),
        TypeExpr::ParticipantSheet => "ParticipantSheet".to_string(),
        TypeExpr::GenotypeRecord => "GenotypeRecord".to_string(),
        TypeExpr::BiovaultContext => "BiovaultContext".to_string(),
        TypeExpr::List(inner) => format!("List[{}]", describe_type(inner)),
        TypeExpr::Map(inner) => format!("Map[String, {}]", describe_type(inner)),
        TypeExpr::Record(fields) => {
            let field_list = fields
                .iter()
                .map(|field| format!("{}: {}", field.name, describe_type(&field.ty)))
                .collect::<Vec<_>>()
                .join(", ");
            format!("Record{{{}}}", field_list)
        }
        TypeExpr::Optional(inner) => format!("{}?", describe_type(inner)),
    }
}

fn is_context_type(ty: &TypeExpr) -> bool {
    match ty {
        TypeExpr::BiovaultContext => true,
        TypeExpr::Optional(inner) => is_context_type(inner),
        _ => false,
    }
}

// Type info for UI
#[derive(Debug, Clone, Serialize)]
pub struct TypeInfo {
    pub base_types: Vec<String>,
    pub common_types: Vec<String>,
}

impl TypeExpr {
    /// Returns all primitive type names
    fn all_primitives() -> Vec<&'static str> {
        vec![
            "String",
            "Bool",
            "File",
            "Directory",
            "ParticipantSheet",
            "GenotypeRecord",
            "BiovaultContext",
        ]
    }

    /// Returns common composite types for inputs
    fn common_input_composites() -> Vec<&'static str> {
        vec![
            "List[File]",
            "List[Directory]",
            "List[GenotypeRecord]",
            "Map[String, File]",
            "Record{bed: File, bim: File, fam: File}",
        ]
    }
}

pub fn get_supported_input_types() -> TypeInfo {
    let mut base_types: Vec<String> = TypeExpr::all_primitives()
        .into_iter()
        .map(|s| s.to_string())
        .collect();

    // Add common composites to base types for inputs
    base_types.extend(
        TypeExpr::common_input_composites()
            .into_iter()
            .map(|s| s.to_string()),
    );

    TypeInfo {
        base_types,
        common_types: vec![
            "File".to_string(),
            "Directory".to_string(),
            "String".to_string(),
            "List[File]".to_string(),
            "Map[String, File]".to_string(),
        ],
    }
}

pub fn get_supported_output_types() -> TypeInfo {
    // Outputs don't support BiovaultContext
    let base_types: Vec<String> = TypeExpr::all_primitives()
        .into_iter()
        .filter(|&t| t != "BiovaultContext")
        .map(|s| s.to_string())
        .collect();

    TypeInfo {
        base_types,
        common_types: vec!["File".to_string(), "Directory".to_string()],
    }
}

pub fn get_supported_parameter_types() -> Vec<String> {
    // Parameters only support String, Bool, and Enum
    vec![
        "String".to_string(),
        "Bool".to_string(),
        "Enum[...]".to_string(),
    ]
}

pub fn get_common_formats() -> Vec<String> {
    vec![
        "csv".to_string(),
        "tsv".to_string(),
        "txt".to_string(),
        "json".to_string(),
        "yaml".to_string(),
        "yml".to_string(),
        "vcf".to_string(),
        "fasta".to_string(),
        "fastq".to_string(),
    ]
}

/// Generate a starter Python script for processing data
pub fn generate_python_script_template(script_name: &str) -> String {
    format!(
        r#"#!/usr/bin/env python3
"""
{script_name}

Starter script for data processing.
Add your custom logic here.
"""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input file path")
    parser.add_argument("--output", required=True, help="Output file path")
    args = parser.parse_args()

    print("🔄 Processing data...")
    print(f"   Input: {{args.input}}")
    print(f"   Output: {{args.output}}")

    # Your processing logic here
    # Example: Read, process, and write data
    with open(args.input) as f_in:
        data = f_in.read()
        # TODO: Add your processing logic
        processed = data  # Passthrough for now

    with open(args.output, 'w') as f_out:
        f_out.write(processed)

    print("✅ Processing complete!")


if __name__ == "__main__":
    main()
"#,
        script_name = script_name
    )
=======
    pub default: Option<YamlValue>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModuleAsset {
    pub path: String,
>>>>>>> main
}

impl ModuleFile {
    pub fn parse_yaml(raw: &str) -> Result<Self> {
        let module: ModuleFile =
            serde_yaml::from_str(raw).context("Failed to parse module spec")?;
        Ok(module)
    }

<<<<<<< HEAD
    pub fn from_module_spec(spec: &ModuleSpec) -> Self {
        let authors = if spec.author.trim().is_empty() {
            Vec::new()
        } else {
            vec![spec.author.clone()]
        };

        let version = spec.version.clone().unwrap_or_else(|| "0.1.0".to_string());

        let runner = ModuleRunnerSpec {
            kind: Some(infer_runner_kind(spec)),
            entrypoint: Some(spec.workflow.clone()),
            runtime: spec.runtime.clone(),
            image: None,
            command: None,
            env: spec.env.clone(),
            syqure: None,
        };

        let inputs = spec
            .inputs
            .iter()
            .map(|input| ModuleFileInputSpec {
                name: input.name.clone(),
                raw_type: input.raw_type.clone(),
                description: input.description.clone(),
                format: format_from_internal(input.format.as_deref(), input.mapping.as_ref()),
                optional: None,
            })
            .collect();

        let outputs = spec
            .outputs
            .iter()
            .map(|output| ModuleFileOutputSpec {
                name: output.name.clone(),
                raw_type: output.raw_type.clone(),
                description: output.description.clone(),
                format: format_from_internal(output.format.as_deref(), None),
                path: output.path.clone(),
                glob: None,
                regex: None,
                cardinality: None,
            })
            .collect();

        let parameters = spec
            .parameters
            .iter()
            .map(|param| ModuleFileParameterSpec {
                name: param.name.clone(),
                raw_type: param.raw_type.clone(),
                default: param.default.as_ref().and_then(value_to_string),
                description: param.description.clone(),
            })
            .collect();

        let assets = if spec.assets.is_empty() {
            None
        } else {
            Some(
                spec.assets
                    .iter()
                    .map(|path| ModuleAsset { path: path.clone() })
                    .collect(),
            )
        };

        ModuleFile {
            api_version: MODULE_API_VERSION.to_string(),
            kind: "Module".to_string(),
            metadata: ModuleFileMetadata {
                name: spec.name.clone(),
                version,
                description: spec.description.clone(),
                authors,
                tags: Vec::new(),
                labels: BTreeMap::new(),
                annotations: BTreeMap::new(),
            },
            spec: ModuleFileSpec {
                runner: Some(runner),
                inputs,
                outputs,
                parameters,
                assets,
                sandbox: None,
                timeout: None,
            },
            manifest: None,
        }
    }

    pub fn to_module_spec(&self) -> Result<ModuleSpec> {
        let runner = self.spec.runner.clone().unwrap_or_default();
        let workflow = runner
            .entrypoint
            .clone()
            .unwrap_or_else(|| "workflow.nf".to_string());
        let author = self.metadata.authors.first().cloned().unwrap_or_default();

        let inputs = self
            .spec
            .inputs
            .iter()
            .map(|input| InputSpec {
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
            })
            .collect();

        let outputs = self
            .spec
            .outputs
            .iter()
            .map(|output| OutputSpec {
                name: output.name.clone(),
                raw_type: output.raw_type.clone(),
                description: output.description.clone(),
                format: output.format.as_ref().and_then(|fmt| fmt.kind.clone()),
                path: output
                    .path
                    .clone()
                    .or_else(|| output.glob.clone())
                    .or_else(|| output.regex.clone()),
            })
            .collect();

        let parameters = self
            .spec
            .parameters
            .iter()
            .map(|param| ParameterSpec {
                name: param.name.clone(),
                raw_type: param.raw_type.clone(),
                description: param.description.clone(),
                default: param.default.clone().map(serde_yaml::Value::String),
                choices: None,
                advanced: None,
            })
            .collect();
=======
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
>>>>>>> main

        let assets = self
            .spec
            .assets
<<<<<<< HEAD
            .clone()
            .unwrap_or_default()
            .into_iter()
            .map(|asset| asset.path)
            .collect();

        Ok(ModuleSpec {
            name: self.metadata.name.clone(),
            author,
            workflow,
            description: self.metadata.description.clone(),
            runtime: runner.runtime.clone(),
            version: Some(self.metadata.version.clone()),
            datasites: None,
            env: runner.env.clone(),
=======
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
>>>>>>> main
            assets,
            parameters,
            inputs,
            outputs,
<<<<<<< HEAD
            steps: Vec::new(),
            runner: Some(runner),
        })
    }
}

fn infer_runner_kind(spec: &ModuleSpec) -> String {
    if let Some(runtime) = &spec.runtime {
        match runtime.as_str() {
            "shell" => return "shell".to_string(),
            "syqure" => return "syqure".to_string(),
            _ => {}
        }
    }
    if spec.workflow.ends_with(".py") {
        return "python".to_string();
    }
    "nextflow".to_string()
}

fn format_from_internal(
    kind: Option<&str>,
    mapping: Option<&HashMap<String, String>>,
) -> Option<ModuleFileFormatSpec> {
    let kind = kind
        .map(|value| value.trim().to_string())
        .filter(|v| !v.is_empty());
    let mapping: BTreeMap<String, String> =
        mapping.cloned().unwrap_or_default().into_iter().collect();
    if kind.is_none() && mapping.is_empty() {
        return None;
    }
    Some(ModuleFileFormatSpec {
        kind,
        mapping,
        key_path: None,
        delimiter: None,
        header: None,
        schema: BTreeMap::new(),
    })
}

fn value_to_string(value: &serde_yaml::Value) -> Option<String> {
    match value {
        serde_yaml::Value::String(s) => Some(s.clone()),
        serde_yaml::Value::Number(n) => Some(n.to_string()),
        serde_yaml::Value::Bool(b) => Some(b.to_string()),
        serde_yaml::Value::Null => None,
        _ => None,
    }
}

/// Scaffold a blank module with optional script generation
pub fn scaffold_blank_module(
    spec: ModuleSpec,
    target_dir: &Path,
    create_python_script: bool,
    script_name: Option<&str>,
) -> Result<ModuleSpec> {
    let mut updated_spec = spec.clone();

    // Scaffold base module structure
    scaffold_from_spec(spec, target_dir)?;

    // Add Python script if requested
    if create_python_script {
        let assets_dir = target_dir.join("assets");
        fs::create_dir_all(&assets_dir).context("Failed to create assets directory")?;

        let script_filename = script_name.unwrap_or("process.py");
        let script_path = assets_dir.join(script_filename);
        let script_content = generate_python_script_template(script_filename);

        fs::write(&script_path, script_content)
            .with_context(|| format!("Failed to write {}", script_path.display()))?;

        // Add to assets list
        updated_spec.assets = vec![script_filename.to_string()];

        // Update module.yaml with assets
<<<<<<<< HEAD:cli/src/module_spec.rs
        let module_yaml_path = target_dir.join("module.yaml");
        let module = ModuleFile::from_module_spec(&updated_spec);
        let yaml =
            serde_yaml::to_string(&module).context("Failed to serialize updated module spec")?;
        fs::write(&module_yaml_path, yaml).context("Failed to update module.yaml with assets")?;
========
        let project_yaml_path = target_dir.join(MODULE_YAML_FILE);
        let module = ModuleFile::from_project_spec(&updated_spec);
        let yaml =
            serde_yaml::to_string(&module).context("Failed to serialize updated module spec")?;
        fs::write(&project_yaml_path, yaml).context("Failed to update module.yaml with assets")?;
>>>>>>>> main:cli/src/project_spec.rs
    }

    Ok(updated_spec)
=======
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
>>>>>>> main
}
