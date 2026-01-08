use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=examples/");
    println!("cargo:rerun-if-changed=examples/examples.yaml");
    println!("cargo:rerun-if-changed=../biovault-beaver/python/src/beaver/__init__.py");
    println!("cargo:rerun-if-changed=../../biovault-beaver/python/src/beaver/__init__.py");
    println!("cargo:rerun-if-env-changed=BIOVAULT_BEAVER_DIR");
    println!("cargo:rerun-if-env-changed=WORKSPACE_ROOT");

    // Expose the biovault-beaver version (from submodule) to the biovault CLI crate at compile time.
    // This is used to pin the PyPI dependency when creating a Jupyter venv.
    println!(
        "cargo:rustc-env=BEAVER_VERSION={}",
        beaver_version_from_submodule()
    );

    // Generate a rust file with embedded examples data
    generate_examples_module();
}

fn beaver_version_from_submodule() -> String {
    let fallback = "0.1.30".to_string();
    let mut candidates = Vec::new();

    if let Ok(path) = env::var("BIOVAULT_BEAVER_DIR") {
        candidates.push(PathBuf::from(path).join("python/src/beaver/__init__.py"));
    }

    if let Ok(root) = env::var("WORKSPACE_ROOT") {
        candidates.push(PathBuf::from(root).join("biovault-beaver/python/src/beaver/__init__.py"));
    }

    if let Ok(manifest_dir) = env::var("CARGO_MANIFEST_DIR") {
        let manifest = PathBuf::from(manifest_dir);
        candidates.push(manifest.join("../biovault-beaver/python/src/beaver/__init__.py"));
        if let Some(root) = manifest.parent().and_then(|p| p.parent()) {
            candidates.push(root.join("biovault-beaver/python/src/beaver/__init__.py"));
        }
    }

    let beaver_init_path = candidates.into_iter().find(|path| path.exists());
    let content = match beaver_init_path.and_then(|path| fs::read_to_string(path).ok()) {
        Some(c) => c,
        None => return fallback,
    };

    for line in content.lines() {
        let trimmed = line.trim();
        if !trimmed.starts_with("__version__") {
            continue;
        }

        if let Some(start) = trimmed.find('"') {
            if let Some(end) = trimmed[start + 1..].find('"') {
                return trimmed[start + 1..start + 1 + end].to_string();
            }
        }
    }

    fallback
}

fn generate_examples_module() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("examples_data.rs");

    let examples_yaml =
        fs::read_to_string("examples/examples.yaml").expect("Failed to read examples.yaml");

    let yaml: serde_yaml::Value =
        serde_yaml::from_str(&examples_yaml).expect("Failed to parse examples.yaml");

    let mut code = String::new();

    // Start with the examples YAML content embedded
    code.push_str(&format!(
        "pub const EXAMPLES_YAML: &str = r#\"{}\"#;\n\n",
        examples_yaml
    ));

    // Create a HashMap of example files
    code.push_str("lazy_static::lazy_static! {\n");
    code.push_str("    pub static ref EXAMPLE_FILES: std::collections::HashMap<&'static str, std::collections::HashMap<&'static str, &'static str>> = {\n");
    code.push_str("        let mut examples = std::collections::HashMap::new();\n");

    if let Some(examples) = yaml.get("examples").and_then(|e| e.as_mapping()) {
        for (key, value) in examples {
            let example_name = key.as_str().unwrap();
            code.push_str("        {\n");
            code.push_str("            let mut files = std::collections::HashMap::new();\n");

            if let Some(files) = value.get("files").and_then(|f| f.as_sequence()) {
                for file_entry in files {
                    // Handle simplified format: just a string with the file path
                    if let Some(file_rel_path) = file_entry.as_str() {
                        // Source path is example_name/file_rel_path
                        let src_path = format!("examples/{}/{}", example_name, file_rel_path);
                        // Destination is the same relative path
                        let dest_path = file_rel_path;

                        if Path::new(&src_path).exists() {
                            let content = fs::read_to_string(&src_path)
                                .unwrap_or_else(|_| panic!("Failed to read {}", src_path));

                            // Escape the content for inclusion in Rust string literal
                            let escaped = content
                                .replace('\\', "\\\\") // Escape backslashes first
                                .replace('"', "\\\"") // Then quotes
                                .replace('\n', "\\n") // Then newlines
                                .replace('\r', "\\r")
                                .replace('\t', "\\t");

                            code.push_str(&format!(
                                "            files.insert(\"{}\", \"{}\");\n",
                                dest_path, escaped
                            ));
                        } else {
                            eprintln!("Warning: File not found: {}", src_path);
                        }
                    }
                }
            }

            code.push_str(&format!(
                "            examples.insert(\"{}\", files);\n",
                example_name
            ));
            code.push_str("        }\n");
        }
    }

    code.push_str("        examples\n");
    code.push_str("    };\n");
    code.push_str("}\n");

    fs::write(&dest_path, code).expect("Failed to write examples_data.rs");
    println!("Generated: {}", dest_path.display());
}
