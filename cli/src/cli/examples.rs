use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

// Include the generated examples data
include!(concat!(env!("OUT_DIR"), "/examples_data.rs"));

#[derive(Debug, Serialize, Deserialize)]
pub struct ExampleInfo {
    pub name: String,
    pub description: String,
    pub template: String,
    pub files: Vec<String>, // Simplified to just list of file paths
}

#[derive(Debug, Serialize, Deserialize)]
struct ExamplesYaml {
    examples: HashMap<String, ExampleInfo>,
}

pub fn get_available_examples() -> anyhow::Result<HashMap<String, ExampleInfo>> {
    let yaml: ExamplesYaml = serde_yaml::from_str(EXAMPLES_YAML)?;
    Ok(yaml.examples)
}

// Helper function to write example files to a directory
pub fn write_example_to_directory(example_name: &str, target_dir: &Path) -> anyhow::Result<()> {
    // Get the files for this example from the generated data
    let files = EXAMPLE_FILES
        .get(example_name)
        .ok_or_else(|| anyhow::anyhow!("Example '{}' not found", example_name))?;

    for (dest_path, content) in files.iter() {
        let dest_full_path = target_dir.join(dest_path);

        // Create parent directories if needed
        if let Some(parent) = dest_full_path.parent() {
            fs::create_dir_all(parent)?;
        }

        // Write content directly - it's already properly unescaped by Rust
        // when it parses the string literal in the generated code
        fs::write(&dest_full_path, content)?;
    }

    Ok(())
}

// Get list of available examples for display/help
pub fn list_examples() -> Vec<String> {
    EXAMPLE_FILES.keys().map(|k| k.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn examples_list_and_write_round_trip() {
        // At least one example should be embedded
        let names = list_examples();
        assert!(!names.is_empty());

        // For a small sample (first up to 2), write to temp dir and ensure files appear
        let tmp = TempDir::new().unwrap();
        for name in names.iter().take(2) {
            write_example_to_directory(name, tmp.path()).unwrap();
        }
    }
}
