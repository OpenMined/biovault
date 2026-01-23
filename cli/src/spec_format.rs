use serde_yaml::Value as YamlValue;
use std::path::Path;

pub const FLOW_API_VERSION: &str = "syftbox.openmined.org/v1alpha1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecFormat {
    Flow,
    Module,
    FlowOverlay,
    LegacyPipeline,
    LegacyProject,
    Unknown,
}

pub fn detect_spec_format(path: &Path, content: &str) -> SpecFormat {
    if let Ok(yaml) = serde_yaml::from_str::<YamlValue>(content) {
        if let Some(api_version) = yaml.get("apiVersion").and_then(|v| v.as_str()) {
            if api_version == FLOW_API_VERSION {
                if let Some(kind) = yaml.get("kind").and_then(|v| v.as_str()) {
                    return match kind {
                        "Flow" => SpecFormat::Flow,
                        "Module" => SpecFormat::Module,
                        "FlowOverlay" => SpecFormat::FlowOverlay,
                        _ => SpecFormat::Unknown,
                    };
                }
                return SpecFormat::Unknown;
            }
        }
    }

    match path.file_name().and_then(|name| name.to_str()) {
        Some("pipeline.yaml") | Some("pipeline.yml") => SpecFormat::LegacyPipeline,
        Some("project.yaml") | Some("project.yml") => SpecFormat::LegacyProject,
        _ => SpecFormat::Unknown,
    }
}
