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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn filter_files_respects_glob_patterns() {
        let base = TempDir::new().unwrap();
        let file_a = base.path().join("sample_A.txt");
        let file_b = base.path().join("control.csv");
        fs::write(&file_a, "data").unwrap();
        fs::write(&file_b, "data").unwrap();

        let inputs = vec![file_a.clone(), file_b.clone()];
        let filtered = filter_files(&inputs, "sample_*.txt").unwrap();

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0], file_a);

        let all = filter_files(&inputs, "*").unwrap();
        assert_eq!(all.len(), 2);
    }

    #[test]
    fn create_extraction_regex_builds_named_groups() {
        let regex = create_extraction_regex("participant_{participant_id}_{date}.txt").unwrap();
        let caps = regex.captures("participant_1234_20240101.txt").unwrap();
        assert_eq!(caps.name("participant_id").unwrap().as_str(), "1234");
        assert_eq!(caps.name("date").unwrap().as_str(), "20240101");
    }

    #[test]
    fn extract_fields_returns_named_and_other_values() {
        let regex =
            create_extraction_regex("participant_{participant_id}_{date}_{batch_id}.txt").unwrap();
        let fields = extract_fields(
            "participant_9876_20231201_batch42.txt",
            &regex,
            "participant_{participant_id}_{date}_{batch_id}.txt",
        )
        .expect("expected captures");

        assert_eq!(fields.participant_id.as_deref(), Some("9876"));
        assert_eq!(fields.date.as_deref(), Some("20231201"));
        assert_eq!(
            fields.other_fields.get("batch_id"),
            Some(&"batch42".to_string())
        );
    }

    #[tokio::test]
    async fn create_generates_csv_without_extract_pattern() {
        let input = TempDir::new().unwrap();
        let output = input.path().join("sheet.csv");

        let file_a = input.path().join("alpha.txt");
        let file_b = input.path().join("beta.csv");
        fs::write(&file_a, "A").unwrap();
        fs::write(&file_b, "B").unwrap();

        create(
            input.path().to_string_lossy().to_string(),
            output.to_string_lossy().to_string(),
            Some("*.txt".to_string()),
            None,
            false,
        )
        .await
        .expect("csv creation succeeds");

        let csv = fs::read_to_string(&output).unwrap();
        assert!(csv.contains("participant_id,genotype_file_path"));
        assert!(csv.contains("alpha"));
        assert!(!csv.contains("beta"));
    }

    #[tokio::test]
    async fn create_errors_when_pattern_does_not_match_and_ignore_off() {
        let input = TempDir::new().unwrap();
        let output = input.path().join("sheet.csv");

        let file = input.path().join("sample_foo.txt");
        fs::write(&file, "data").unwrap();

        let result = create(
            input.path().to_string_lossy().to_string(),
            output.to_string_lossy().to_string(),
            Some("*.txt".to_string()),
            Some("sample_{participant_id}_{date}.txt".to_string()),
            false,
        )
        .await;

        assert!(result.is_err());
        assert!(!output.exists());
    }
}
