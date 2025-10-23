use anyhow::{Context, Result};
use blake3::Hasher;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

use crate::project_spec::{InputSpec, OutputSpec, ParameterSpec, ProjectSpec};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectMetadata {
    pub name: String,
    pub author: String,
    pub workflow: String,
    pub template: Option<String>,
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
pub struct ProjectFileNode {
    pub name: String,
    pub path: String,
    pub is_dir: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub children: Vec<ProjectFileNode>,
}

pub fn load_project_metadata(project_root: &Path) -> Result<Option<ProjectMetadata>> {
    let yaml_path = project_root.join("project.yaml");
    if !yaml_path.exists() {
        return Ok(None);
    }

    let spec = ProjectSpec::load(&yaml_path)
        .with_context(|| format!("Invalid project spec in {}", yaml_path.display()))?;

    Ok(Some(ProjectMetadata {
        name: spec.name,
        author: spec.author,
        workflow: spec.workflow,
        template: spec.template,
        version: Some(spec.version.unwrap_or_else(|| "1.0.0".to_string())),
        assets: spec.assets,
        parameters: spec.parameters,
        inputs: spec.inputs,
        outputs: spec.outputs,
    }))
}

pub fn save_project_metadata(project_root: &Path, metadata: &ProjectMetadata) -> Result<()> {
    let yaml_path = project_root.join("project.yaml");
    let mut assets: Vec<String> = metadata
        .assets
        .iter()
        .map(|s| s.replace('\\', "/"))
        .collect();
    assets.sort();
    assets.dedup();

    let spec = ProjectSpec {
        name: metadata.name.clone(),
        author: metadata.author.clone(),
        workflow: metadata.workflow.clone(),
        template: metadata.template.clone(),
        version: Some(
            metadata
                .version
                .clone()
                .unwrap_or_else(|| "1.0.0".to_string()),
        ),
        assets,
        parameters: metadata.parameters.clone(),
        inputs: metadata.inputs.clone(),
        outputs: metadata.outputs.clone(),
    };

    let yaml_str = serde_yaml::to_string(&spec)
        .with_context(|| format!("Failed to serialize {}", yaml_path.display()))?;
    fs::write(&yaml_path, yaml_str)
        .with_context(|| format!("Failed to write {}", yaml_path.display()))?;

    Ok(())
}

pub fn project_yaml_hash(project_root: &Path) -> Result<Option<String>> {
    let yaml_path = project_root.join("project.yaml");
    if !yaml_path.exists() {
        return Ok(None);
    }

    let bytes =
        fs::read(&yaml_path).with_context(|| format!("Failed to read {}", yaml_path.display()))?;
    let mut hasher = Hasher::new();
    hasher.update(&bytes);
    Ok(Some(hasher.finalize().to_hex().to_string()))
}

pub fn build_project_file_tree(project_root: &Path) -> Result<Vec<ProjectFileNode>> {
    fn collect(root: &Path, current: &Path) -> Result<Vec<ProjectFileNode>> {
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
                entries.push(ProjectFileNode {
                    name,
                    path: relative,
                    is_dir: true,
                    children,
                });
            } else if name != "project.yaml" {
                entries.push(ProjectFileNode {
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

    if !project_root.exists() {
        return Ok(Vec::new());
    }

    collect(project_root, project_root)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn load_and_save_metadata_round_trip() {
        let tmp = TempDir::new().unwrap();
        let root = tmp.path();

        let metadata = ProjectMetadata {
            name: "example".into(),
            author: "test@example.com".into(),
            workflow: "workflow.nf".into(),
            template: Some("sheet".into()),
            version: Some("1.2.3".into()),
            assets: vec!["schema.yaml".into(), "src/main.py".into()],
            parameters: vec![crate::project_spec::ParameterSpec {
                name: "flag".into(),
                raw_type: "Bool".into(),
                description: None,
                default: None,
                choices: None,
                advanced: None,
            }],
            inputs: vec![crate::project_spec::InputSpec {
                name: "rows".into(),
                raw_type: "List[GenotypeRecord]".into(),
                description: None,
                format: None,
                path: None,
                mapping: None,
            }],
            outputs: vec![crate::project_spec::OutputSpec {
                name: "sheet".into(),
                raw_type: "ParticipantSheet".into(),
                description: None,
                format: None,
                path: Some("results/sheet.csv".into()),
            }],
        };

        save_project_metadata(root, &metadata).unwrap();
        let loaded = load_project_metadata(root).unwrap().unwrap();
        assert_eq!(loaded.name, "example");
        assert_eq!(loaded.author, "test@example.com");
        assert_eq!(loaded.workflow, "workflow.nf");
        assert_eq!(loaded.template.as_deref(), Some("sheet"));
        assert_eq!(loaded.assets.len(), 2);
        assert_eq!(loaded.version.as_deref(), Some("1.2.3"));
        assert_eq!(loaded.parameters.len(), 1);
        assert_eq!(loaded.inputs.len(), 1);
        assert_eq!(loaded.outputs.len(), 1);

        let digest = project_yaml_hash(root).unwrap().unwrap();
        assert_eq!(digest.len(), 64);

        // Modify file to ensure hash changes
        let yaml_path = root.join("project.yaml");
        let mut yaml = std::fs::read_to_string(&yaml_path).unwrap();
        yaml.push_str("# comment\n");
        std::fs::write(&yaml_path, yaml).unwrap();

        let new_digest = project_yaml_hash(root).unwrap().unwrap();
        assert_eq!(new_digest.len(), 64);
        assert_ne!(digest, new_digest);
    }

    #[test]
    fn build_tree_skips_project_yaml() {
        let tmp = TempDir::new().unwrap();
        let root = tmp.path();
        fs::write(root.join("project.yaml"), "name: test").unwrap();
        fs::write(root.join("workflow.nf"), "// workflow").unwrap();
        fs::create_dir_all(root.join("bioscript")).unwrap();
        fs::write(root.join("bioscript/lib.rs"), "").unwrap();

        let tree = build_project_file_tree(root).unwrap();
        assert_eq!(tree.len(), 2);
        assert!(tree.iter().any(|node| node.name == "workflow.nf"));
        assert!(tree
            .iter()
            .any(|node| node.name == "bioscript" && node.is_dir));
    }
}
