use anyhow::{Context, Result};
use blake3::Hasher;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

use crate::module_spec::{InputSpec, ModuleSpec, OutputSpec, ParameterSpec};
use crate::types::ModuleYaml;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModuleMetadata {
    pub name: String,
    pub author: String,
    pub workflow: String,
    #[serde(alias = "template")]
    pub runtime: Option<String>,
    #[serde(default)]
    pub version: Option<String>,
    pub assets: Vec<String>,
    #[serde(default)]
    pub parameters: Vec<ParameterSpec>,
    #[serde(default)]
    pub inputs: Vec<InputSpec>,
    #[serde(default)]
    pub outputs: Vec<OutputSpec>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ModuleFileNode {
    pub name: String,
    pub path: String,
    pub is_dir: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub children: Vec<ModuleFileNode>,
}

pub fn load_module_metadata(module_root: &Path) -> Result<Option<ModuleMetadata>> {
    let yaml_path = module_root.join("module.yaml");
    if !yaml_path.exists() {
        return Ok(None);
    }

    let spec = ModuleSpec::load(&yaml_path)
        .with_context(|| format!("Invalid module spec in {}", yaml_path.display()))?;

    Ok(Some(ModuleMetadata {
        name: spec.name,
        author: spec.author,
        workflow: spec.workflow,
        runtime: spec.runtime,
        version: Some(spec.version.unwrap_or_else(|| "1.0.0".to_string())),
        assets: spec.assets,
        parameters: spec.parameters,
        inputs: spec.inputs,
        outputs: spec.outputs,
    }))
}

pub fn save_module_metadata(module_root: &Path, metadata: &ModuleMetadata) -> Result<()> {
    let yaml_path = module_root.join("module.yaml");
    let mut assets: Vec<String> = metadata
        .assets
        .iter()
        .map(|s| s.replace('\\', "/"))
        .collect();
    assets.sort();
    assets.dedup();

    let spec = ModuleSpec {
        name: metadata.name.clone(),
        author: metadata.author.clone(),
        workflow: metadata.workflow.clone(),
        description: None,
        runtime: metadata.runtime.clone(),
        version: Some(
            metadata
                .version
                .clone()
                .unwrap_or_else(|| "1.0.0".to_string()),
        ),
        datasites: None,
        env: Default::default(),
        assets,
        parameters: metadata.parameters.clone(),
        inputs: metadata.inputs.clone(),
        outputs: metadata.outputs.clone(),
        steps: Vec::new(),
        runner: None,
    };

    let module = crate::module_spec::ModuleFile::from_module_spec(&spec);
    let yaml_str = serde_yaml::to_string(&module)
        .with_context(|| format!("Failed to serialize {}", yaml_path.display()))?;
    fs::write(&yaml_path, yaml_str)
        .with_context(|| format!("Failed to write {}", yaml_path.display()))?;

    Ok(())
}

pub fn module_yaml_hash(module_root: &Path) -> Result<Option<String>> {
    let yaml_path = module_root.join("module.yaml");
    if !yaml_path.exists() {
        return Ok(None);
    }

    let bytes =
        fs::read(&yaml_path).with_context(|| format!("Failed to read {}", yaml_path.display()))?;
    if let Ok(module) = serde_yaml::from_slice::<ModuleYaml>(&bytes) {
        if let Some(hashes) = &module.b3_hashes {
            let mut hash_content = String::new();
            hash_content.push_str(&module.name);

            let workflow_key = if module.workflow.trim().is_empty() {
                "workflow.nf"
            } else {
                module.workflow.as_str()
            };
            if let Some(workflow_hash) = hashes
                .get(workflow_key)
                .or_else(|| hashes.get("workflow.nf"))
            {
                hash_content.push_str(workflow_hash);
            }

            let mut sorted_hashes: Vec<_> = hashes.iter().collect();
            sorted_hashes.sort_by_key(|entry| entry.0);
            for (file, hash) in sorted_hashes {
                hash_content.push_str(file);
                hash_content.push_str(hash);
            }

            return Ok(Some(
                blake3::hash(hash_content.as_bytes()).to_hex().to_string(),
            ));
        }
    }

    let mut hasher = Hasher::new();
    hasher.update(&bytes);
    Ok(Some(hasher.finalize().to_hex().to_string()))
}

pub fn build_module_file_tree(module_root: &Path) -> Result<Vec<ModuleFileNode>> {
    fn collect(root: &Path, current: &Path) -> Result<Vec<ModuleFileNode>> {
        let mut entries = Vec::new();
        for entry in fs::read_dir(current)
            .with_context(|| format!("Failed to read directory {}", current.display()))?
        {
            let entry = entry?;
            let path = entry.path();
            let name = entry.file_name().to_string_lossy().to_string();

            if name == "." || name == ".." {
                continue;
            }

            if name == ".venv" {
                continue;
            }

            let relative = path
                .strip_prefix(root)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");

            let file_type = entry.file_type()?;
            if file_type.is_dir() {
                let children = collect(root, &path)?;
                entries.push(ModuleFileNode {
                    name,
                    path: relative,
                    is_dir: true,
                    children,
                });
            } else if name != "module.yaml" {
                entries.push(ModuleFileNode {
                    name,
                    path: relative,
                    is_dir: false,
                    children: Vec::new(),
                });
            }
        }

        entries.sort_by(|a, b| match (a.is_dir, b.is_dir) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => a.name.to_lowercase().cmp(&b.name.to_lowercase()),
        });

        Ok(entries)
    }

    if !module_root.exists() {
        return Ok(Vec::new());
    }

    collect(module_root, module_root)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn load_and_save_metadata_round_trip() {
        let tmp = TempDir::new().unwrap();
        let root = tmp.path();

        let metadata = ModuleMetadata {
            name: "example".into(),
            author: "test@example.com".into(),
            workflow: "workflow.nf".into(),
            runtime: Some("sheet".into()),
            version: Some("1.2.3".into()),
            assets: vec!["schema.yaml".into(), "src/main.py".into()],
            parameters: vec![crate::module_spec::ParameterSpec {
                name: "flag".into(),
                raw_type: "Bool".into(),
                description: None,
                default: None,
                choices: None,
                advanced: None,
            }],
            inputs: vec![crate::module_spec::InputSpec {
                name: "rows".into(),
                raw_type: "List[GenotypeRecord]".into(),
                description: None,
                format: None,
                path: None,
                mapping: None,
            }],
            outputs: vec![crate::module_spec::OutputSpec {
                name: "sheet".into(),
                raw_type: "ParticipantSheet".into(),
                description: None,
                format: None,
                path: Some("results/sheet.csv".into()),
            }],
        };

        save_module_metadata(root, &metadata).unwrap();
        let loaded = load_module_metadata(root).unwrap().unwrap();
        assert_eq!(loaded.name, "example");
        assert_eq!(loaded.author, "test@example.com");
        assert_eq!(loaded.workflow, "workflow.nf");
        assert_eq!(loaded.runtime.as_deref(), Some("sheet"));
        assert_eq!(loaded.assets.len(), 2);
        assert_eq!(loaded.version.as_deref(), Some("1.2.3"));
        assert_eq!(loaded.parameters.len(), 1);
        assert_eq!(loaded.inputs.len(), 1);
        assert_eq!(loaded.outputs.len(), 1);

        let digest = module_yaml_hash(root).unwrap().unwrap();
        assert_eq!(digest.len(), 64);

        // Modify file to ensure hash changes
        let yaml_path = root.join("module.yaml");
        let mut yaml = std::fs::read_to_string(&yaml_path).unwrap();
        yaml.push_str("# comment\n");
        std::fs::write(&yaml_path, yaml).unwrap();

        let new_digest = module_yaml_hash(root).unwrap().unwrap();
        assert_eq!(new_digest.len(), 64);
        assert_ne!(digest, new_digest);
    }

    #[test]
    fn build_tree_skips_module_yaml() {
        let tmp = TempDir::new().unwrap();
        let root = tmp.path();
        fs::write(root.join("module.yaml"), "name: test").unwrap();
        fs::write(root.join("workflow.nf"), "// workflow").unwrap();
        fs::create_dir_all(root.join("bioscript")).unwrap();
        fs::write(root.join("bioscript/lib.rs"), "").unwrap();

        let tree = build_module_file_tree(root).unwrap();
        assert_eq!(tree.len(), 2);
        assert!(tree.iter().any(|node| node.name == "workflow.nf"));
        assert!(tree
            .iter()
            .any(|node| node.name == "bioscript" && node.is_dir));
    }
}
