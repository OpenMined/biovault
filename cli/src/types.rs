use serde::{Deserialize, Serialize};
use std::collections::HashMap;
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

    pub fn save(&self, path: &PathBuf) -> anyhow::Result<()> {
        let yaml = serde_yaml::to_string(self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectYaml {
    pub name: String,
    pub author: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub datasites: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub participants: Option<Vec<String>>,
    pub workflow: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assets: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub b3_hashes: Option<HashMap<String, String>>,
}

impl ProjectYaml {
    pub fn from_file(path: &PathBuf) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let project: ProjectYaml = serde_yaml::from_str(&content)?;
        Ok(project)
    }

    pub fn save(&self, path: &PathBuf) -> anyhow::Result<()> {
        let yaml = serde_yaml::to_string(self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}
