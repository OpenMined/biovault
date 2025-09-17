use std::env;
use std::fs;
use std::path::Path;

fn main() {
    println!("cargo:rerun-if-changed=examples/");
    println!("cargo:rerun-if-changed=examples/examples.yaml");

    // Generate a rust file with embedded examples data
    generate_examples_module();
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
