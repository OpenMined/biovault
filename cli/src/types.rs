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
            rules: vec![
                // Global read for the recipient datasite
                PermissionRule {
                    pattern: "**".to_string(),
                    access: AccessControl {
                        read: vec![datasite_email.to_string()],
                        write: vec![],
                        admin: vec![],
                    },
                },
                // Allow recipient to write results back into this submission
                PermissionRule {
                    pattern: "results/**/*".to_string(),
                    access: AccessControl {
                        read: vec![datasite_email.to_string()],
                        write: vec![datasite_email.to_string()],
                        admin: vec![],
                    },
                },
            ],
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
    pub template: Option<String>,
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
        assert_eq!(perms.rules.len(), 2);
        perms.save(&p).unwrap();
        let read_back: SyftPermissions =
            serde_yaml::from_str(&fs::read_to_string(&p).unwrap()).unwrap();
        assert_eq!(read_back.rules.len(), 2);
    }

    #[test]
    fn project_yaml_round_trip_and_error() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("project.yaml");
        let proj = ProjectYaml {
            name: "N".into(),
            author: "A".into(),
            datasites: Some(vec!["x@y".into()]),
            participants: Some(vec!["P1".into()]),
            workflow: "wf".into(),
            template: None,
            assets: Some(vec!["a".into()]),
            b3_hashes: None,
        };
        proj.save(&p).unwrap();
        let loaded = ProjectYaml::from_file(&p).unwrap();
        assert_eq!(loaded.name, "N");
        assert_eq!(loaded.workflow, "wf");

        let bad = tmp.path().join("bad.yaml");
        fs::write(&bad, "not: [valid").unwrap();
        assert!(ProjectYaml::from_file(&bad).is_err());
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
