use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};

macro_rules! wln {
    ($buf:expr) => {
        writeln!($buf).map_err(|e| anyhow!(e.to_string()))?
    };
    ($buf:expr, $($arg:tt)*) => {
        writeln!($buf, $($arg)*).map_err(|e| anyhow!(e.to_string()))?
    };
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectSpec {
    pub name: String,
    pub author: String,
    pub workflow: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub template: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default)]
    pub assets: Vec<String>,
    #[serde(default)]
    pub parameters: Vec<ParameterSpec>,
    #[serde(default)]
    pub inputs: Vec<InputSpec>,
    #[serde(default)]
    pub outputs: Vec<OutputSpec>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping: Option<HashMap<String, String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputSpec {
    pub name: String,
    #[serde(rename = "type")]
    pub raw_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
enum ParameterType {
    String,
    Bool,
    Enum(Vec<String>),
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
    Optional(Box<TypeExpr>),
}

impl ProjectSpec {
    pub fn load(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("Failed to read project spec at {}", path.display()))?;
        let spec: ProjectSpec = serde_yaml::from_str(&raw)
            .with_context(|| format!("Failed to parse project spec at {}", path.display()))?;
        Ok(spec)
    }
}

pub fn scaffold_from_spec(mut spec: ProjectSpec, target_dir: &Path) -> Result<ProjectSpec> {
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

    let template_name = spec
        .template
        .clone()
        .unwrap_or_else(|| "dynamic-nextflow".to_string());
    spec.template = Some(template_name);

    let project_yaml_path = target_dir.join("project.yaml");
    let workflow_path = target_dir.join(&spec.workflow);
    let assets_dir = target_dir.join("assets");

    let yaml = serde_yaml::to_string(&spec).context("Failed to serialize project spec")?;
    fs::write(&project_yaml_path, yaml).context("Failed to write project.yaml")?;

    if let Some(parent) = workflow_path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create directory {}", parent.display()))?;
    }
    let workflow_contents = generate_workflow_stub(&spec)?;
    fs::write(&workflow_path, workflow_contents).context("Failed to write workflow stub")?;

    // Note: template.nf is NOT written to project folder - it's a security boundary
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

fn parse_type_expr(raw: &str) -> Result<TypeExpr> {
    let trimmed = raw.trim();
    let (base, optional) = if let Some(stripped) = trimmed.strip_suffix('?') {
        (stripped.trim(), true)
    } else {
        (trimmed, false)
    };

    let parsed = if let Some(inner) = base.strip_prefix("List[") {
        let inner = inner
            .strip_suffix(']')
            .ok_or_else(|| anyhow!("Invalid List type '{}': missing closing ]", raw))?;
        TypeExpr::List(Box::new(parse_type_expr(inner)?))
    } else if let Some(inner) = base.strip_prefix("Map[String,") {
        let inner = inner
            .strip_suffix(']')
            .ok_or_else(|| anyhow!("Invalid Map type '{}': missing closing ]", raw))?;
        let inner = inner.trim_start_matches(',').trim();
        if inner.is_empty() {
            bail!("Map type '{}' is missing value type", raw);
        }
        TypeExpr::Map(Box::new(parse_type_expr(inner)?))
    } else {
        match base {
            "String" => TypeExpr::String,
            "Bool" => TypeExpr::Bool,
            "File" => TypeExpr::File,
            "Directory" => TypeExpr::Directory,
            "ParticipantSheet" => TypeExpr::ParticipantSheet,
            "GenotypeRecord" => TypeExpr::GenotypeRecord,
            "BiovaultContext" => TypeExpr::BiovaultContext,
            other => bail!(
                "Unsupported type '{}' (supported primitives: String, Bool, File, Directory, ParticipantSheet, GenotypeRecord, BiovaultContext, List[..], Map[String, ..])",
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

pub fn generate_template_nf(spec: &ProjectSpec) -> Result<String> {
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

fn generate_workflow_stub(spec: &ProjectSpec) -> Result<String> {
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
        TypeExpr::Map(_) => vec![format!(r#"{indent}println "[bv] {name}: ${{{name}}}""#)],
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
        | TypeExpr::Map(_) => format!("def {name} = params.{name}"),
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
