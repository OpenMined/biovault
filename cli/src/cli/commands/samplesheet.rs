use anyhow::{Context, Result};
use regex::Regex;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use tracing::{debug, info, warn};

pub struct ExtractedFields {
    participant_id: Option<String>,
    date: Option<String>,
    other_fields: HashMap<String, String>,
}

pub async fn create(
    input_dir: String,
    output_file: String,
    file_filter: Option<String>,
    extract_cols: Option<String>,
    ignore_errors: bool,
) -> Result<()> {
    let input_path = Path::new(&input_dir);
    let output_path = Path::new(&output_file);

    if !input_path.exists() {
        anyhow::bail!("Input directory does not exist: {}", input_dir);
    }

    if !input_path.is_dir() {
        anyhow::bail!("Input path is not a directory: {}", input_dir);
    }

    if output_path.exists() {
        anyhow::bail!(
            "Output file already exists: {}. Please choose a different name.",
            output_file
        );
    }

    let filter_pattern = file_filter.unwrap_or_else(|| "*".to_string());

    let extract_pattern = extract_cols.unwrap_or_default();

    info!("Creating sample sheet from: {}", input_dir);
    info!("Output file: {}", output_file);
    info!("File filter: {}", filter_pattern);
    if !extract_pattern.is_empty() {
        info!("Extract pattern: {}", extract_pattern);
    }

    let entries = fs::read_dir(input_path)?
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            if let Ok(metadata) = entry.metadata() {
                metadata.is_file()
            } else {
                false
            }
        })
        .map(|entry| entry.path())
        .collect::<Vec<_>>();

    let filtered_files = filter_files(&entries, &filter_pattern)?;

    info!("Found {} files matching filter", filtered_files.len());
    for file in &filtered_files {
        println!("  - {}", file.display());
    }
    println!("Total files matching filter: {}", filtered_files.len());

    let mut csv_rows = Vec::new();
    let mut matched_count = 0;
    let mut failed_files = Vec::new();

    if extract_pattern.is_empty() {
        for file_path in &filtered_files {
            if let Some(file_name) = file_path.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    let participant_id = file_name_str
                        .strip_suffix(".txt")
                        .or_else(|| file_name_str.strip_suffix(".csv"))
                        .unwrap_or(file_name_str)
                        .to_string();
                    let absolute_path = file_path
                        .canonicalize()
                        .unwrap_or_else(|_| file_path.clone());
                    csv_rows.push((participant_id, absolute_path));
                    matched_count += 1;
                }
            }
        }
    } else {
        let pattern_regex = create_extraction_regex(&extract_pattern)?;

        for file_path in &filtered_files {
            if let Some(file_name) = file_path.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    match extract_fields(file_name_str, &pattern_regex, &extract_pattern) {
                        Some(fields) => {
                            if let Some(participant_id) = fields.participant_id {
                                let absolute_path = file_path
                                    .canonicalize()
                                    .unwrap_or_else(|_| file_path.clone());
                                csv_rows.push((participant_id, absolute_path));
                                matched_count += 1;
                            }
                        }
                        None => {
                            failed_files.push(file_name_str.to_string());
                        }
                    }
                }
            }
        }
    }

    if !failed_files.is_empty() {
        warn!("Files that didn't match the extraction pattern:");
        for file in &failed_files {
            println!("  ✗ {}", file);
        }
        if !ignore_errors && !failed_files.is_empty() {
            anyhow::bail!(
                "{} files didn't match the extraction pattern. Use --ignore to add them anyway.",
                failed_files.len()
            );
        }
    }

    // Sort by participant_id when extract_cols is used
    if !extract_pattern.is_empty() {
        csv_rows.sort_by(|a, b| a.0.cmp(&b.0));
        info!("Sorted {} participants by participant_id", csv_rows.len());
    }

    let mut csv_content = String::from("participant_id,genotype_file_path\n");
    for (participant_id, file_path) in &csv_rows {
        csv_content.push_str(&format!("{},{}\n", participant_id, file_path.display()));
    }

    fs::write(output_path, csv_content)
        .with_context(|| format!("Failed to write CSV file: {}", output_file))?;

    println!("\n✓ Sample sheet created: {}", output_file);
    println!("  Added {} participants to CSV", matched_count);
    if !failed_files.is_empty() && ignore_errors {
        println!(
            "  Ignored {} files that didn't match pattern",
            failed_files.len()
        );
    }

    Ok(())
}

fn filter_files(files: &[PathBuf], pattern: &str) -> Result<Vec<PathBuf>> {
    if pattern == "*" {
        return Ok(files.to_vec());
    }

    let glob = globset::GlobBuilder::new(pattern)
        .literal_separator(false)
        .build()
        .context("Invalid file filter pattern")?;
    let matcher = glob.compile_matcher();

    Ok(files
        .iter()
        .filter(|path| {
            if let Some(file_name) = path.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    return matcher.is_match(file_name_str);
                }
            }
            false
        })
        .cloned()
        .collect())
}

fn create_extraction_regex(pattern: &str) -> Result<Regex> {
    let mut regex_pattern = pattern.to_string();

    let field_regex = Regex::new(r"\{([^}]+)\}")?;
    let mut replacements = Vec::new();

    for cap in field_regex.captures_iter(pattern) {
        if let Some(field_match) = cap.get(0) {
            if let Some(field_name) = cap.get(1) {
                let field_str = field_match.as_str();
                let name_str = field_name.as_str();
                let capture_group = format!("(?P<{}>[^_]+)", name_str.replace('-', "_"));
                replacements.push((field_str.to_string(), capture_group));
            }
        }
    }

    for (field_str, capture_group) in replacements {
        regex_pattern = regex_pattern.replace(&field_str, &capture_group);
    }

    regex_pattern = regex_pattern.replace(".", "\\.");

    debug!("Generated regex pattern: {}", regex_pattern);

    Regex::new(&regex_pattern).context("Failed to create extraction regex")
}

fn extract_fields(
    filename: &str,
    pattern_regex: &Regex,
    _original_pattern: &str,
) -> Option<ExtractedFields> {
    if let Some(captures) = pattern_regex.captures(filename) {
        let mut fields = ExtractedFields {
            participant_id: None,
            date: None,
            other_fields: HashMap::new(),
        };

        if let Some(participant_id) = captures.name("participant_id") {
            fields.participant_id = Some(participant_id.as_str().to_string());
        }
        if let Some(date) = captures.name("date") {
            fields.date = Some(date.as_str().to_string());
        }

        for (idx, mat) in captures.iter().enumerate() {
            if idx == 0 {
                continue;
            }
            if let Some(m) = mat {
                if let Some(name) = pattern_regex.capture_names().nth(idx).flatten() {
                    if name != "participant_id" && name != "date" {
                        fields
                            .other_fields
                            .insert(name.to_string(), m.as_str().to_string());
                    }
                }
            }
        }

        Some(fields)
    } else {
        None
    }
}
