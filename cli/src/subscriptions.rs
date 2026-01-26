use anyhow::{Context, Result};
use globset::Glob;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Subscriptions {
    #[serde(default = "default_version")]
    pub version: u32,
    #[serde(default)]
    pub defaults: Defaults,
    #[serde(default)]
    pub rules: Vec<Rule>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Defaults {
    #[serde(default)]
    pub action: Action,
}

impl Default for Defaults {
    fn default() -> Self {
        Self {
            action: Action::Block,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Rule {
    #[serde(default)]
    pub action: Action,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub datasite: Option<String>,
    pub path: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Action {
    Allow,
    Pause,
    Block,
    #[serde(other)]
    Unknown,
}

impl Default for Action {
    fn default() -> Self {
        Action::Block
    }
}

impl Action {
    fn normalize(self) -> Action {
        match self {
            Action::Unknown => Action::Block,
            other => other,
        }
    }
}

fn default_version() -> u32 {
    1
}

pub fn default_rules() -> Vec<Rule> {
    vec![
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "public/crypto/did.json".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "public/biovault/datasets.yaml".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "public/biovault/datasets/*/dataset.yaml".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "app_data/biovault/*.yaml".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "**/syft.pub.yaml".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "**/*.request".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "**/*.response".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "shared/biovault/**".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "shared/flows/**".to_string(),
        },
        Rule {
            action: Action::Allow,
            datasite: Some("*".to_string()),
            path: "shared/syqure/**".to_string(),
        },
    ]
}

pub fn default_config() -> Subscriptions {
    Subscriptions {
        version: 1,
        defaults: Defaults {
            action: Action::Block,
        },
        rules: default_rules(),
    }
}

pub fn load(path: &Path) -> Result<Subscriptions> {
    let raw = match std::fs::read_to_string(path) {
        Ok(raw) => raw,
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            return Ok(default_config())
        }
        Err(err) => return Err(err.into()),
    };
    let mut cfg: Subscriptions = serde_yaml::from_str(&raw)
        .with_context(|| format!("Invalid syft.sub.yaml: {}", path.display()))?;
    normalize_config(&mut cfg);
    Ok(cfg)
}

pub fn save(path: &Path, cfg: &Subscriptions) -> Result<()> {
    let mut normalized = cfg.clone();
    normalize_config(&mut normalized);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let raw = serde_yaml::to_string(&normalized)?;
    let tmp = path.with_extension("tmp");
    std::fs::write(&tmp, raw)?;
    std::fs::rename(tmp, path)?;
    Ok(())
}

pub fn action_for_path(cfg: &Subscriptions, owner: &str, rel_path: &str) -> Action {
    let rel = normalize_path(rel_path);
    let (datasite, rest) = split_datasite(&rel);
    if datasite.is_empty() {
        return cfg.defaults.action.clone().normalize();
    }
    if datasite.eq_ignore_ascii_case(owner) {
        return Action::Allow;
    }

    let mut action = cfg.defaults.action.clone().normalize();
    let full_path = rel.clone();
    for rule in &cfg.rules {
        if !rule_matches(rule, &datasite, &full_path, &rest) {
            continue;
        }
        action = rule.action.clone().normalize();
    }
    action
}

fn normalize_config(cfg: &mut Subscriptions) {
    if cfg.version == 0 {
        cfg.version = 1;
    }
    cfg.defaults.action = cfg.defaults.action.clone().normalize();
    for rule in &mut cfg.rules {
        rule.action = rule.action.clone().normalize();
    }
}

fn rule_matches(rule: &Rule, datasite: &str, full_path: &str, rest: &str) -> bool {
    if rule.path.trim().is_empty() {
        return false;
    }
    if let Some(ds) = &rule.datasite {
        if !matches_glob(ds, datasite) {
            return false;
        }
        return matches_glob(&rule.path, rest);
    }
    matches_glob(&rule.path, full_path)
}

fn matches_glob(pattern: &str, target: &str) -> bool {
    let pattern = normalize_path(pattern);
    let target = normalize_path(target);
    let matcher = match Glob::new(&pattern) {
        Ok(glob) => glob.compile_matcher(),
        Err(_) => return false,
    };
    matcher.is_match(target)
}

fn normalize_path(raw: &str) -> String {
    let mut out = raw.replace('\\', "/");
    while out.starts_with('/') {
        out.remove(0);
    }
    while out.contains("//") {
        out = out.replace("//", "/");
    }
    out
}

fn split_datasite(rel: &str) -> (String, String) {
    if rel.is_empty() {
        return ("".to_string(), "".to_string());
    }
    let mut parts = rel.splitn(2, '/');
    let ds = parts.next().unwrap_or("").to_string();
    let rest = parts.next().unwrap_or("").to_string();
    (ds, rest)
}
