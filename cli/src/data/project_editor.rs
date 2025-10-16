use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use serde_yaml::{Mapping, Value};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectMetadata {
    pub name: String,
    pub author: String,
    pub workflow: String,
    pub template: Option<String>,
    pub assets: Vec<String>,
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

    let content = fs::read_to_string(&yaml_path)
        .with_context(|| format!("Failed to read {}", yaml_path.display()))?;
    let value: Value = serde_yaml::from_str(&content)
        .with_context(|| format!("Invalid YAML in {}", yaml_path.display()))?;

    let mapping = value.as_mapping().cloned().unwrap_or_else(Mapping::new);

    Ok(Some(ProjectMetadata {
        name: mapping
            .get(Value::String("name".into()))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_string(),
        author: mapping
            .get(Value::String("author".into()))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_string(),
        workflow: mapping
            .get(Value::String("workflow".into()))
            .and_then(Value::as_str)
            .unwrap_or("workflow.nf")
            .to_string(),
        template: mapping
            .get(Value::String("template".into()))
            .and_then(Value::as_str)
            .map(|s| s.to_string()),
        assets: mapping
            .get(Value::String("assets".into()))
            .and_then(Value::as_sequence)
            .map(|seq| {
                seq.iter()
                    .filter_map(Value::as_str)
                    .map(|s| s.replace('\\', "/"))
                    .collect()
            })
            .unwrap_or_default(),
    }))
}

pub fn save_project_metadata(project_root: &Path, metadata: &ProjectMetadata) -> Result<()> {
    let yaml_path = project_root.join("project.yaml");
    let mut value = if yaml_path.exists() {
        let content = fs::read_to_string(&yaml_path)
            .with_context(|| format!("Failed to read {}", yaml_path.display()))?;
        serde_yaml::from_str::<Value>(&content)
            .with_context(|| format!("Invalid YAML in {}", yaml_path.display()))?
    } else {
        Value::Mapping(Mapping::new())
    };

    let map = value
        .as_mapping_mut()
        .ok_or_else(|| anyhow::anyhow!("project.yaml must be a mapping"))?;

    map.insert(
        Value::String("name".into()),
        Value::String(metadata.name.clone()),
    );
    map.insert(
        Value::String("author".into()),
        Value::String(metadata.author.clone()),
    );
    map.insert(
        Value::String("workflow".into()),
        Value::String(metadata.workflow.clone()),
    );

    if let Some(template) = metadata.template.as_ref().and_then(|t| {
        let trimmed = t.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }) {
        map.insert(Value::String("template".into()), Value::String(template));
    } else {
        map.remove(Value::String("template".into()));
    }

    let mut assets: Vec<String> = metadata
        .assets
        .iter()
        .map(|s| s.replace('\\', "/"))
        .collect();
    assets.sort();
    assets.dedup();

    map.insert(
        Value::String("assets".into()),
        Value::Sequence(
            assets
                .into_iter()
                .map(Value::String)
                .collect(),
        ),
    );

    let yaml_str = serde_yaml::to_string(&value)
        .with_context(|| format!("Failed to serialize {}", yaml_path.display()))?;
    fs::write(&yaml_path, yaml_str)
        .with_context(|| format!("Failed to write {}", yaml_path.display()))?;

    Ok(())
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
            assets: vec!["schema.yaml".into(), "src/main.py".into()],
        };

        save_project_metadata(root, &metadata).unwrap();
        let loaded = load_project_metadata(root).unwrap().unwrap();
        assert_eq!(loaded.name, "example");
        assert_eq!(loaded.author, "test@example.com");
        assert_eq!(loaded.workflow, "workflow.nf");
        assert_eq!(loaded.template.as_deref(), Some("sheet"));
        assert_eq!(loaded.assets.len(), 2);
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
