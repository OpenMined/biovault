use anyhow::Result;
use chrono::{DateTime, Utc};
use regex::Regex;
use rusqlite::params;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

#[cfg(test)]
use sha2::{Digest, Sha256};
#[cfg(test)]
use std::io::Read;

use super::{BioVaultDb, GenotypeMetadata};

#[derive(Debug, Serialize, Deserialize)]
pub struct ExtensionInfo {
    pub extension: String,
    pub count: usize,
    pub total_size: u64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FileInfo {
    pub path: String,
    pub size: u64,
    pub modified: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ScanResult {
    pub path: String,
    pub scanned_at: String,
    pub extensions: Vec<ExtensionInfo>,
    pub files: Vec<FileInfo>,
    pub total_files: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FileRecord {
    pub id: i64,
    pub file_path: String,
    pub file_hash: String,
    pub file_type: Option<String>,
    pub file_size: Option<u64>,
    pub data_type: Option<String>,
    pub source: Option<String>,
    pub grch_version: Option<String>,
    pub row_count: Option<i64>,
    pub chromosome_count: Option<i64>,
    pub inferred_sex: Option<String>,
    pub status: Option<String>,
    pub processing_error: Option<String>,
    pub participant_id: Option<String>,
    pub participant_name: Option<String>,
    pub created_at: String,
    pub updated_at: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ImportResult {
    pub imported: usize,
    pub skipped: usize,
    pub errors: Vec<String>,
    pub files: Vec<FileRecord>,
}

#[derive(Debug, Clone)]
pub struct CsvFileImport {
    pub file_path: String,
    pub participant_id: Option<String>,
    pub data_type: Option<String>,
    pub source: Option<String>,
    pub grch_version: Option<String>,
    pub row_count: Option<i64>,
    pub chromosome_count: Option<i64>,
    pub inferred_sex: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PatternSuggestion {
    pub pattern: String,
    pub description: String,
    pub example: String,
    pub sample_extractions: Vec<(String, String)>, // (filename, extracted_id)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SuggestPatternsResult {
    pub suggestions: Vec<PatternSuggestion>,
    pub sample_files: Vec<String>,
}

/// Scan directory for files by type
pub fn scan(path: &str, extension: Option<&str>, recursive: bool) -> Result<ScanResult> {
    let path_buf = PathBuf::from(path);
    if !path_buf.exists() {
        anyhow::bail!("Path does not exist: {}", path);
    }

    let mut extension_counts: HashMap<String, (usize, u64)> = HashMap::new();
    let mut files = Vec::new();

    let walker = if recursive {
        WalkDir::new(path).follow_links(false)
    } else {
        WalkDir::new(path).max_depth(1).follow_links(false)
    };

    for entry in walker.into_iter().filter_map(|e| e.ok()) {
        if !entry.file_type().is_file() {
            continue;
        }

        let file_path = entry.path();
        let file_ext = file_path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| format!(".{}", e));

        // Filter by extension if specified
        if let Some(filter) = extension {
            let filters: Vec<&str> = filter.split(',').collect();
            if let Some(ref ext) = file_ext {
                if !filters.iter().any(|f| {
                    let normalized_filter = if f.starts_with('.') {
                        f.to_string()
                    } else {
                        format!(".{}", f)
                    };
                    ext == &normalized_filter
                }) {
                    continue;
                }
            } else {
                continue; // Skip files without extension when filter is set
            }
        }

        // Get file size
        let size = entry.metadata().ok().map(|m| m.len()).unwrap_or(0);

        // Get modified time
        let modified = entry
            .metadata()
            .ok()
            .and_then(|m| m.modified().ok())
            .map(|t| {
                let datetime: DateTime<Utc> = t.into();
                datetime.to_rfc3339()
            });

        // Track extension stats
        if let Some(ext) = &file_ext {
            let entry = extension_counts.entry(ext.clone()).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += size;
        }

        // Add to files list
        files.push(FileInfo {
            path: file_path.to_string_lossy().to_string(),
            size,
            modified,
        });
    }

    // Convert extension counts to sorted vec
    let mut extensions: Vec<ExtensionInfo> = extension_counts
        .into_iter()
        .map(|(extension, (count, total_size))| ExtensionInfo {
            extension,
            count,
            total_size,
        })
        .collect();

    extensions.sort_by(|a, b| b.count.cmp(&a.count));

    Ok(ScanResult {
        path: path.to_string(),
        scanned_at: Utc::now().to_rfc3339(),
        extensions,
        total_files: files.len(),
        files,
    })
}

/// Import files into database
pub fn import(
    db: &BioVaultDb,
    path: &str,
    extension: Option<&str>,
    recursive: bool,
    pattern: Option<&str>,
) -> Result<ImportResult> {
    // First, scan for files
    let scan_result = scan(path, extension, recursive)?;

    let mut imported = 0;
    let mut skipped = 0;
    let mut errors = Vec::new();
    let mut imported_file_ids = Vec::new();

    for file_info in scan_result.files {
        // Extract participant ID if pattern is provided
        let participant_id = if let Some(pat) = pattern {
            // Pass full path to extract_id_from_pattern (needed for {parent} and other path-based patterns)
            extract_id_from_pattern(&file_info.path, pat)
        } else {
            None
        };

        match import_file_with_participant(db, &file_info, participant_id.as_deref()) {
            Ok(Some(file_id)) => {
                imported += 1;
                imported_file_ids.push(file_id);
            }
            Ok(None) => {
                skipped += 1;
            }
            Err(e) => {
                errors.push(format!("{}: {}", file_info.path, e));
            }
        }
    }

    // Fetch imported file records
    let mut files = Vec::new();
    for file_id in imported_file_ids {
        if let Ok(Some(record)) = get_file_by_id(db, file_id) {
            files.push(record);
        }
    }

    Ok(ImportResult {
        imported,
        skipped,
        errors,
        files,
    })
}

/// Import files from CSV with all metadata
pub fn import_from_csv(db: &BioVaultDb, csv_imports: Vec<CsvFileImport>) -> Result<ImportResult> {
    let mut imported = 0;
    let mut skipped = 0;
    let mut errors = Vec::new();
    let mut imported_file_ids = Vec::new();

    for csv_row in csv_imports {
        match import_file_with_metadata(db, &csv_row) {
            Ok(Some(file_id)) => {
                imported += 1;
                imported_file_ids.push(file_id);
            }
            Ok(None) => {
                skipped += 1;
            }
            Err(e) => {
                errors.push(format!("{}: {}", csv_row.file_path, e));
            }
        }
    }

    // Fetch imported file records
    let mut files = Vec::new();
    for file_id in imported_file_ids {
        if let Ok(Some(record)) = get_file_by_id(db, file_id) {
            files.push(record);
        }
    }

    Ok(ImportResult {
        imported,
        skipped,
        errors,
        files,
    })
}

/// Import single file with full metadata from CSV
fn import_file_with_metadata(db: &BioVaultDb, csv_row: &CsvFileImport) -> Result<Option<i64>> {
    use crate::cli::download_cache::calculate_blake3;
    use std::path::Path;

    let file_path = &csv_row.file_path;

    // Calculate file hash using BLAKE3
    let path = Path::new(file_path);
    if !path.exists() {
        anyhow::bail!("File not found: {}", file_path);
    }

    let metadata = std::fs::metadata(path)?;
    let file_size = metadata.len();

    let file_hash = calculate_blake3(path).unwrap_or_else(|e| {
        eprintln!("Warning: Failed to hash {}: {}", file_path, e);
        format!("error_{}", file_size)
    });

    // Extract file type (extension)
    let file_type = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| format!(".{}", e));

    // Get or create participant if ID provided
    let db_participant_id: Option<i64> = if let Some(pid) = &csv_row.participant_id {
        Some(get_or_create_participant(db, pid)?)
    } else {
        None
    };

    // Check if file already exists
    let existing: Option<(String, i64)> = db
        .conn
        .query_row(
            "SELECT file_hash, id FROM files WHERE file_path = ?1",
            params![file_path],
            |row| Ok((row.get(0)?, row.get(1)?)),
        )
        .ok();

    let file_id = if let Some((existing_hash, existing_id)) = existing {
        if existing_hash == file_hash {
            // File already imported - update metadata
            db.conn.execute(
                "UPDATE files SET participant_id = ?1, data_type = ?2, source = ?3, grch_version = ?4, updated_at = CURRENT_TIMESTAMP WHERE id = ?5",
                params![
                    db_participant_id,
                    csv_row.data_type,
                    csv_row.source,
                    csv_row.grch_version,
                    existing_id
                ],
            )?;
            existing_id
        } else {
            // File exists but hash changed - update everything
            db.conn.execute(
                "UPDATE files SET file_hash = ?1, file_size = ?2, participant_id = ?3, data_type = ?4, source = ?5, grch_version = ?6, updated_at = CURRENT_TIMESTAMP WHERE id = ?7",
                params![
                    file_hash,
                    file_size as i64,
                    db_participant_id,
                    csv_row.data_type,
                    csv_row.source,
                    csv_row.grch_version,
                    existing_id
                ],
            )?;
            existing_id
        }
    } else {
        // Insert new file
        db.conn.execute(
            "INSERT INTO files (participant_id, file_path, file_hash, file_type, file_size, data_type, source, grch_version) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                db_participant_id,
                file_path,
                file_hash,
                file_type,
                file_size as i64,
                csv_row.data_type,
                csv_row.source,
                csv_row.grch_version
            ],
        )?;
        db.conn.last_insert_rowid()
    };

    // Save genotype metadata if provided and data_type is Genotype
    if csv_row.data_type.as_deref() == Some("Genotype")
        && (csv_row.row_count.is_some() || csv_row.chromosome_count.is_some())
    {
        // Delete existing metadata if any
        db.conn.execute(
            "DELETE FROM genotype_metadata WHERE file_id = ?1",
            params![file_id],
        )?;

        // Insert new metadata
        db.conn.execute(
                "INSERT INTO genotype_metadata (file_id, source, grch_version, row_count, chromosome_count) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![
                    file_id,
                    csv_row.source,
                    csv_row.grch_version,
                    csv_row.row_count,
                    csv_row.chromosome_count
                ],
            )?;
    }

    // Update participant inferred_sex if provided
    if let (Some(pid), Some(sex)) = (db_participant_id, &csv_row.inferred_sex) {
        db.conn.execute(
            "UPDATE participants SET inferred_sex = ?1 WHERE id = ?2",
            params![sex, pid],
        )?;
    }

    Ok(Some(file_id))
}

/// Import single file into database with optional participant assignment
fn import_file_with_participant(
    db: &BioVaultDb,
    file_info: &FileInfo,
    participant_id: Option<&str>,
) -> Result<Option<i64>> {
    use crate::cli::download_cache::calculate_blake3;
    use std::path::Path;

    let file_path = &file_info.path;

    // Calculate file hash using BLAKE3
    let file_hash = calculate_blake3(Path::new(file_path)).unwrap_or_else(|e| {
        eprintln!("Warning: Failed to hash {}: {}", file_path, e);
        format!("error_{}", file_info.size)
    });

    // Extract file type (extension)
    let file_type = Path::new(file_path)
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| format!(".{}", e));

    // Get or create participant if ID provided
    let db_participant_id: Option<i64> = if let Some(pid) = participant_id {
        Some(get_or_create_participant(db, pid)?)
    } else {
        None
    };

    // Check if file already exists
    let existing: Option<(String, i64)> = db
        .conn
        .query_row(
            "SELECT file_hash, id FROM files WHERE file_path = ?1",
            params![file_path],
            |row| Ok((row.get(0)?, row.get(1)?)),
        )
        .ok();

    if let Some((existing_hash, existing_id)) = existing {
        if existing_hash == file_hash {
            // File already imported - update participant if provided
            if db_participant_id.is_some() {
                db.conn.execute(
                    "UPDATE files SET participant_id = ?1, updated_at = CURRENT_TIMESTAMP WHERE id = ?2",
                    params![db_participant_id, existing_id],
                )?;
            }
            return Ok(None); // Skip, already exists
        } else {
            // File exists but hash changed - update it
            db.conn.execute(
                "UPDATE files SET file_hash = ?1, file_size = ?2, participant_id = ?3, updated_at = CURRENT_TIMESTAMP WHERE id = ?4",
                params![file_hash, file_info.size as i64, db_participant_id, existing_id],
            )?;
            return Ok(Some(existing_id));
        }
    }

    // Insert new file
    db.conn.execute(
        "INSERT INTO files (participant_id, file_path, file_hash, file_type, file_size) VALUES (?1, ?2, ?3, ?4, ?5)",
        params![db_participant_id, file_path, file_hash, file_type, file_info.size as i64],
    )?;

    let file_id = db.conn.last_insert_rowid();
    Ok(Some(file_id))
}

/// Get or create participant by ID
fn get_or_create_participant(db: &BioVaultDb, participant_id: &str) -> Result<i64> {
    // Try to get existing participant
    match db.conn.query_row(
        "SELECT id FROM participants WHERE participant_id = ?1",
        params![participant_id],
        |row| row.get(0),
    ) {
        Ok(id) => Ok(id),
        Err(rusqlite::Error::QueryReturnedNoRows) => {
            // Create new participant
            db.conn.execute(
                "INSERT INTO participants (participant_id) VALUES (?1)",
                params![participant_id],
            )?;
            Ok(db.conn.last_insert_rowid())
        }
        Err(e) => Err(e.into()),
    }
}

/// Get file by ID
pub fn get_file_by_id(db: &BioVaultDb, file_id: i64) -> Result<Option<FileRecord>> {
    let result = db.conn.query_row(
        "SELECT f.id, f.file_path, f.file_hash, f.file_type, f.file_size,
                f.data_type, g.source, g.grch_version,
                g.row_count, g.chromosome_count, g.inferred_sex,
                f.status, f.processing_error,
                p.participant_id, p.participant_id, f.created_at, f.updated_at
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
         LEFT JOIN genotype_metadata g ON f.id = g.file_id
         WHERE f.id = ?1",
        params![file_id],
        |row| {
            Ok(FileRecord {
                id: row.get(0)?,
                file_path: row.get(1)?,
                file_hash: row.get(2)?,
                file_type: row.get(3)?,
                file_size: row.get::<_, Option<i64>>(4)?.map(|s| s as u64),
                data_type: row.get(5)?,
                source: row.get(6)?,
                grch_version: row.get(7)?,
                row_count: row.get(8)?,
                chromosome_count: row.get(9)?,
                inferred_sex: row.get(10)?,
                status: row.get(11)?,
                processing_error: row.get(12)?,
                participant_id: row.get(13)?,
                participant_name: row.get(14)?,
                created_at: row.get(15)?,
                updated_at: row.get(16)?,
            })
        },
    );

    match result {
        Ok(record) => Ok(Some(record)),
        Err(rusqlite::Error::QueryReturnedNoRows) => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// List files with optional filters
pub fn list_files(
    db: &BioVaultDb,
    extension: Option<&str>,
    participant: Option<&str>,
    unassigned: bool,
    limit: Option<usize>,
) -> Result<Vec<FileRecord>> {
    let mut query = String::from(
        "SELECT f.id, f.file_path, f.file_hash, f.file_type, f.file_size,
                f.data_type, g.source, g.grch_version,
                g.row_count, g.chromosome_count, g.inferred_sex,
                f.status, f.processing_error,
                p.participant_id, p.participant_id, f.created_at, f.updated_at
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
         LEFT JOIN genotype_metadata g ON f.id = g.file_id
         WHERE 1=1",
    );

    let mut params_vec: Vec<Box<dyn rusqlite::ToSql>> = Vec::new();

    if let Some(ext) = extension {
        query.push_str(" AND f.file_type = ?");
        params_vec.push(Box::new(ext.to_string()));
    }

    if let Some(pid) = participant {
        query.push_str(" AND p.participant_id = ?");
        params_vec.push(Box::new(pid.to_string()));
    }

    if unassigned {
        query.push_str(" AND f.participant_id IS NULL");
    }

    query.push_str(" ORDER BY f.created_at DESC");

    if let Some(l) = limit {
        query.push_str(&format!(" LIMIT {}", l));
    }

    let params_refs: Vec<&dyn rusqlite::ToSql> = params_vec.iter().map(|p| p.as_ref()).collect();

    let mut stmt = db.conn.prepare(&query)?;
    let rows = stmt.query_map(params_refs.as_slice(), |row| {
        Ok(FileRecord {
            id: row.get(0)?,
            file_path: row.get(1)?,
            file_hash: row.get(2)?,
            file_type: row.get(3)?,
            file_size: row.get::<_, Option<i64>>(4)?.map(|s| s as u64),
            data_type: row.get(5)?,
            source: row.get(6)?,
            grch_version: row.get(7)?,
            row_count: row.get(8)?,
            chromosome_count: row.get(9)?,
            inferred_sex: row.get(10)?,
            status: row.get(11)?,
            processing_error: row.get(12)?,
            participant_id: row.get(13)?,
            participant_name: row.get(14)?,
            created_at: row.get(15)?,
            updated_at: row.get(16)?,
        })
    })?;

    let mut files = Vec::new();
    for row in rows {
        files.push(row?);
    }

    Ok(files)
}

/// Suggest participant ID patterns from filenames
pub fn suggest_patterns(
    path: &str,
    extension: Option<&str>,
    recursive: bool,
) -> Result<SuggestPatternsResult> {
    // First scan for files
    let scan_result = scan(path, extension, recursive)?;

    if scan_result.files.is_empty() {
        return Ok(SuggestPatternsResult {
            suggestions: vec![],
            sample_files: vec![],
        });
    }

    // Extract filenames
    let filenames: Vec<String> = scan_result
        .files
        .iter()
        .filter_map(|f| Path::new(&f.path).file_name())
        .filter_map(|n| n.to_str())
        .map(|s| s.to_string())
        .collect();

    if filenames.is_empty() {
        return Ok(SuggestPatternsResult {
            suggestions: vec![],
            sample_files: vec![],
        });
    }

    let mut suggestions = Vec::new();
    let first = &filenames[0];

    // Pattern 1: case_XXXX_ pattern
    if first.contains("case_") {
        let re = Regex::new(r"case_(\d+)").unwrap();
        if re.is_match(first) {
            let pattern = "case_{id}_*".to_string();
            let mut sample_extractions = Vec::new();

            for filename in filenames.iter().take(3) {
                if let Some(captures) = re.captures(filename) {
                    if let Some(id) = captures.get(1) {
                        sample_extractions.push((filename.clone(), id.as_str().to_string()));
                    }
                }
            }

            suggestions.push(PatternSuggestion {
                pattern: pattern.clone(),
                description: "Case ID pattern (e.g., case_0001_...)".to_string(),
                example: first.clone(),
                sample_extractions,
            });
        }
    }

    // Pattern 2: XXXX_X_X_ pattern (numbers at start)
    let re_start = Regex::new(r"^(\d+)_").unwrap();
    if re_start.is_match(first) {
        let pattern = "{id}_*".to_string();
        let mut sample_extractions = Vec::new();

        for filename in filenames.iter().take(3) {
            if let Some(captures) = re_start.captures(filename) {
                if let Some(id) = captures.get(1) {
                    sample_extractions.push((filename.clone(), id.as_str().to_string()));
                }
            }
        }

        suggestions.push(PatternSuggestion {
            pattern,
            description: "Leading ID pattern (e.g., 001_...)".to_string(),
            example: first.clone(),
            sample_extractions,
        });
    }

    // Pattern 3: Generic number sequences (3+ digits)
    let re_numbers = Regex::new(r"\d{3,}").unwrap();
    if let Some(mat) = re_numbers.find(first) {
        let start = mat.start();
        let end = mat.end();
        let before = &first[..start];
        let after = &first[end..];
        let pattern = format!("{}{{id}}{}", before, after);

        let mut sample_extractions = Vec::new();
        for filename in filenames.iter().take(3) {
            if let Some(mat) = re_numbers.find(filename) {
                sample_extractions.push((filename.clone(), mat.as_str().to_string()));
            }
        }

        suggestions.push(PatternSuggestion {
            pattern,
            description: "Numeric ID pattern".to_string(),
            example: first.clone(),
            sample_extractions,
        });
    }

    // Pattern 4: participant_XXX or sample_XXX
    let re_participant = Regex::new(r"(participant|sample|subject|patient)_(\d+)").unwrap();
    if let Some(captures) = re_participant.captures(first) {
        let prefix = captures.get(1).unwrap().as_str();
        let pattern = format!("{}_{{id}}_*", prefix);

        let mut sample_extractions = Vec::new();
        for filename in filenames.iter().take(3) {
            if let Some(captures) = re_participant.captures(filename) {
                if let Some(id) = captures.get(2) {
                    sample_extractions.push((filename.clone(), id.as_str().to_string()));
                }
            }
        }

        suggestions.push(PatternSuggestion {
            pattern,
            description: format!("{} ID pattern", prefix).to_string(),
            example: first.clone(),
            sample_extractions,
        });
    }

    // Pattern 5: Directory name as ID
    // Check if files are organized in subdirectories where directory name is the ID
    let first_file_path = scan_result.files.first().map(|f| Path::new(&f.path));
    if let Some(first_path) = first_file_path {
        if let Some(parent) = first_path.parent() {
            if let Some(dir_name) = parent.file_name().and_then(|n| n.to_str()) {
                // Check if directory name looks like an ID (contains numbers or is alphanumeric)
                let re_dir_id = Regex::new(r"^\d+$|^[a-zA-Z0-9_-]+$").unwrap();
                if re_dir_id.is_match(dir_name) && dir_name.len() >= 3 {
                    let mut sample_extractions = Vec::new();
                    let mut seen_dirs = std::collections::HashSet::new();

                    for file_scan in scan_result.files.iter().take(10) {
                        if let Some(file_parent) = Path::new(&file_scan.path).parent() {
                            if let Some(parent_name) =
                                file_parent.file_name().and_then(|n| n.to_str())
                            {
                                if seen_dirs.insert(parent_name.to_string()) {
                                    sample_extractions.push((
                                        format!("{}/...", parent_name),
                                        parent_name.to_string(),
                                    ));
                                }
                                if sample_extractions.len() >= 3 {
                                    break;
                                }
                            }
                        }
                    }

                    if !sample_extractions.is_empty() {
                        suggestions.push(PatternSuggestion {
                            pattern: "{parent}".to_string(),
                            description: "Directory name as participant ID".to_string(),
                            example: format!("{}/...", dir_name),
                            sample_extractions,
                        });
                    }
                }
            }
        }
    }

    Ok(SuggestPatternsResult {
        suggestions,
        sample_files: filenames.into_iter().take(5).collect(),
    })
}

/// Extract participant ID from filepath using glob-style pattern
///
/// Pattern Syntax:
/// - `{id}` - Marks the capture group - this is what gets extracted as the participant ID
/// - `*` - Matches any characters (non-greedy, doesn't cross `/` boundaries)
/// - `**` - Matches any characters including `/` (for multi-level paths)
/// - `?` - Matches single character
/// - Literal text matches exactly
///
/// Special Shortcuts (for convenience):
/// - `{parent}` - Shortcut for extracting parent directory name
/// - `{filename}` - Shortcut for extracting filename without extension
/// - `{basename}` - Shortcut for extracting filename with extension
///
/// Examples:
/// - `{id}/*` → `/data/PT001/scan.txt` extracts `PT001` (parent directory)
/// - `hu{id}/*` → `/data/hu627574/file.txt` extracts `627574` (after "hu" in parent)
/// - `{id}_*.txt` → `123456_data.txt` extracts `123456` (before underscore)
/// - `case_{id}_*` → `case_789_report.pdf` extracts `789` (between delimiters)
/// - `*/{id}/*.txt` → `/path/PT001/scan.txt` extracts `PT001` (middle directory)
/// - `participant_{id}.{sample_id}.*` → `participant_PT001.S001.csv` extracts `PT001` (first capture)
///
/// The `{id}` token can be used multiple times, but only the first occurrence will be captured.
pub fn extract_id_from_pattern(filepath: &str, pattern: &str) -> Option<String> {
    let path = Path::new(filepath);

    // Handle special shortcut patterns
    match pattern {
        "{parent}" | "{dirname}" | "{dir}" | "{id}/*" => {
            return path
                .parent()
                .and_then(|p| p.file_name())
                .and_then(|n| n.to_str())
                .map(|s| s.to_string());
        }
        "{filename}" => {
            return path
                .file_stem()
                .and_then(|n| n.to_str())
                .map(|s| s.to_string());
        }
        "{basename}" => {
            return path
                .file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.to_string());
        }
        _ => {}
    }

    // Pattern must contain {id} to capture
    if !pattern.contains("{id}") {
        return None;
    }

    // Build regex from glob-style pattern
    // Determine if pattern should match against full path or just filename
    let match_target = if pattern.contains('/') {
        // Pattern has path separators, match against full path
        filepath
    } else {
        // Pattern is filename-only, match against just the filename
        path.file_name()?.to_str()?
    };

    // Convert pattern to regex
    let mut regex_pattern = String::new();
    let mut chars = pattern.chars().peekable();

    while let Some(ch) = chars.next() {
        match ch {
            '*' => {
                // Check for ** (multi-level wildcard)
                if chars.peek() == Some(&'*') {
                    chars.next(); // consume second *
                    regex_pattern.push_str(".*"); // Match anything including /
                } else {
                    // Single * doesn't cross directory boundaries
                    regex_pattern.push_str("[^/]*");
                }
            }
            '?' => regex_pattern.push('.'), // Single character wildcard
            '.' | '+' | '^' | '$' | '(' | ')' | '[' | ']' | '\\' | '|' => {
                // Escape regex special characters
                regex_pattern.push('\\');
                regex_pattern.push(ch);
            }
            '{' => {
                // Check for {id} token
                let rest: String = chars.clone().collect();
                if rest.starts_with("id}") {
                    // Consume "id}"
                    chars.next(); // i
                    chars.next(); // d
                    chars.next(); // }
                                  // Insert capture group for ID (alphanumeric, underscore, hyphen)
                    regex_pattern.push_str(r"([a-zA-Z0-9_-]+)");
                } else {
                    // Unknown token, treat as literal
                    regex_pattern.push(ch);
                }
            }
            _ => regex_pattern.push(ch),
        }
    }

    // Anchor the pattern to match the full string
    let regex_pattern = format!("^{}$", regex_pattern);

    // Try to match and extract ID
    let re = Regex::new(&regex_pattern).ok()?;
    let captures = re.captures(match_target)?;
    captures.get(1).map(|m| m.as_str().to_string())
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ParticipantRecord {
    pub id: i64,
    pub participant_id: String,
    pub created_at: String,
    pub file_count: i64,
}

/// List all participants with file counts
pub fn list_participants(db: &BioVaultDb) -> Result<Vec<ParticipantRecord>> {
    let mut stmt = db.conn.prepare(
        "SELECT p.id, p.participant_id, p.created_at, COUNT(f.id) as file_count
         FROM participants p
         LEFT JOIN files f ON f.participant_id = p.id
         GROUP BY p.id, p.participant_id, p.created_at
         ORDER BY p.created_at DESC",
    )?;

    let participants = stmt
        .query_map([], |row| {
            Ok(ParticipantRecord {
                id: row.get(0)?,
                participant_id: row.get(1)?,
                created_at: row.get(2)?,
                file_count: row.get(3)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;

    Ok(participants)
}

/// Delete a participant and remove all associated files
pub fn delete_participant(db: &BioVaultDb, id: i64) -> Result<usize> {
    let files_deleted = db
        .conn
        .execute("DELETE FROM files WHERE participant_id = ?1", params![id])?;

    let rows = db
        .conn
        .execute("DELETE FROM participants WHERE id = ?1", params![id])?;

    if rows == 0 {
        anyhow::bail!("Participant with id {} not found", id);
    }

    Ok(files_deleted)
}

/// Delete multiple participants and remove all associated files
pub fn delete_participants_bulk(db: &BioVaultDb, ids: &[i64]) -> Result<usize> {
    if ids.is_empty() {
        return Ok(0);
    }

    let placeholders = ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");

    let delete_files_query = format!(
        "DELETE FROM files WHERE participant_id IN ({})",
        placeholders
    );
    db.conn
        .execute(&delete_files_query, rusqlite::params_from_iter(ids.iter()))?;

    let delete_participants_query =
        format!("DELETE FROM participants WHERE id IN ({})", placeholders);
    let rows = db.conn.execute(
        &delete_participants_query,
        rusqlite::params_from_iter(ids.iter()),
    )?;

    Ok(rows)
}

/// Delete a file record from the catalog
pub fn delete_file(db: &BioVaultDb, file_id: i64) -> Result<()> {
    let rows = db
        .conn
        .execute("DELETE FROM files WHERE id = ?1", params![file_id])?;

    if rows == 0 {
        anyhow::bail!("File with id {} not found", file_id);
    }

    Ok(())
}

/// Delete multiple file records from the catalog
pub fn delete_files_bulk(db: &BioVaultDb, ids: &[i64]) -> Result<usize> {
    if ids.is_empty() {
        return Ok(0);
    }

    let placeholders = ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");
    let delete_query = format!("DELETE FROM files WHERE id IN ({})", placeholders);
    let rows = db
        .conn
        .execute(&delete_query, rusqlite::params_from_iter(ids.iter()))?;

    Ok(rows)
}

/// Link a file to a participant
pub fn link_file_to_participant(
    db: &BioVaultDb,
    file_id: i64,
    participant_id: &str,
) -> Result<FileRecord> {
    // Get or create participant
    let db_participant_id = get_or_create_participant(db, participant_id)?;

    // Update file
    let rows = db.conn.execute(
        "UPDATE files SET participant_id = ?1, updated_at = CURRENT_TIMESTAMP WHERE id = ?2",
        params![db_participant_id, file_id],
    )?;

    if rows == 0 {
        anyhow::bail!("File with id {} not found", file_id);
    }

    // Return updated file record
    get_file_by_id(db, file_id)?.ok_or_else(|| anyhow::anyhow!("File not found after update"))
}

/// Bulk link multiple files to participants
/// Accepts a HashMap of file_path -> participant_id
pub fn link_files_bulk(
    db: &BioVaultDb,
    file_participant_map: &std::collections::HashMap<String, String>,
) -> Result<usize> {
    if file_participant_map.is_empty() {
        return Ok(0);
    }

    let mut updated = 0;

    // Start transaction for better performance
    let tx = db.conn.unchecked_transaction()?;

    for (file_path, participant_id) in file_participant_map {
        // Get or create participant
        let db_participant_id = {
            match tx.query_row(
                "SELECT id FROM participants WHERE participant_id = ?1",
                params![participant_id],
                |row| row.get(0),
            ) {
                Ok(id) => id,
                Err(rusqlite::Error::QueryReturnedNoRows) => {
                    tx.execute(
                        "INSERT INTO participants (participant_id) VALUES (?1)",
                        params![participant_id],
                    )?;
                    tx.last_insert_rowid()
                }
                Err(e) => return Err(e.into()),
            }
        };

        // Update file by path
        let rows = tx.execute(
            "UPDATE files SET participant_id = ?1, updated_at = CURRENT_TIMESTAMP WHERE file_path = ?2",
            params![db_participant_id, file_path],
        )?;

        updated += rows;
    }

    tx.commit()?;
    Ok(updated)
}

/// Unlink a file from its participant
pub fn unlink_file(db: &BioVaultDb, file_id: i64) -> Result<FileRecord> {
    let rows = db.conn.execute(
        "UPDATE files SET participant_id = NULL, updated_at = CURRENT_TIMESTAMP WHERE id = ?1",
        params![file_id],
    )?;

    if rows == 0 {
        anyhow::bail!("File with id {} not found", file_id);
    }

    // Return updated file record
    get_file_by_id(db, file_id)?.ok_or_else(|| anyhow::anyhow!("File not found after update"))
}

// File hashing

pub fn hash_file(path: &str) -> Result<String> {
    let content = fs::read(path)?;
    Ok(blake3::hash(&content).to_hex().to_string())
}

// Fast import - add files to queue for background processing

pub fn import_files_as_pending(db: &BioVaultDb, files: Vec<CsvFileImport>) -> Result<ImportResult> {
    let conn = db.connection();
    let mut imported = 0;
    let mut skipped = 0;
    let mut errors = Vec::new();

    for file_info in files {
        // Check if file already exists
        let existing: Option<i64> = conn
            .query_row(
                "SELECT id FROM files WHERE file_path = ?1",
                [&file_info.file_path],
                |row| row.get(0),
            )
            .ok();

        if existing.is_some() {
            skipped += 1;
            continue;
        }

        // Get or create participant
        let participant_id = if let Some(pid) = &file_info.participant_id {
            match get_or_create_participant(db, pid) {
                Ok(id) => Some(id),
                Err(e) => {
                    errors.push(format!("Failed to create participant {}: {}", pid, e));
                    continue;
                }
            }
        } else {
            None
        };

        // Get file size
        let file_size = match std::fs::metadata(&file_info.file_path) {
            Ok(meta) => Some(meta.len() as i64),
            Err(_) => None,
        };

        // Insert file with status='pending' and placeholder hash
        let result = conn.execute(
            "INSERT INTO files (participant_id, file_path, file_hash, file_type, file_size, data_type, status, queue_added_at, created_at, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, 'pending', CURRENT_TIMESTAMP, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
            rusqlite::params![
                participant_id,
                &file_info.file_path,
                "pending", // Placeholder hash until processed
                None::<String>, // file_type
                file_size,
                file_info.data_type.as_deref().unwrap_or("Unknown"),
            ],
        );

        match result {
            Ok(_) => {
                // If this is a genotype file with metadata, create genotype_metadata row
                let data_type = file_info.data_type.as_deref().unwrap_or("Unknown");
                if data_type == "Genotype"
                    && (file_info.source.is_some() || file_info.grch_version.is_some())
                {
                    let file_id = conn.last_insert_rowid();
                    let meta_result = conn.execute(
                        "INSERT INTO genotype_metadata (file_id, source, grch_version, created_at, updated_at)
                         VALUES (?1, ?2, ?3, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
                        rusqlite::params![
                            file_id,
                            file_info.source.as_deref(),
                            file_info.grch_version.as_deref(),
                        ],
                    );
                    if let Err(e) = meta_result {
                        errors.push(format!(
                            "Failed to create metadata for {}: {}",
                            file_info.file_path, e
                        ));
                    }
                }
                imported += 1;
            }
            Err(e) => errors.push(format!("Failed to import {}: {}", file_info.file_path, e)),
        }
    }

    Ok(ImportResult {
        imported,
        skipped,
        errors,
        files: Vec::new(), // Files are in pending state, not yet fully imported
    })
}

// Queue processing support

#[derive(Debug, Clone)]
pub struct PendingFile {
    pub id: i64,
    pub file_path: String,
    pub data_type: Option<String>,
    pub participant_id: Option<i64>,
}

pub fn get_pending_files(db: &BioVaultDb, limit: usize) -> Result<Vec<PendingFile>> {
    let conn = db.connection();

    let mut stmt = conn.prepare(
        "SELECT id, file_path, data_type, participant_id
         FROM files
         WHERE status = 'pending'
         ORDER BY queue_added_at ASC
         LIMIT ?1",
    )?;

    let files = stmt
        .query_map([limit], |row| {
            Ok(PendingFile {
                id: row.get(0)?,
                file_path: row.get(1)?,
                data_type: row.get(2)?,
                participant_id: row.get(3)?,
            })
        })?
        .collect::<std::result::Result<Vec<_>, _>>()?;

    Ok(files)
}

pub fn update_file_status(
    db: &BioVaultDb,
    file_id: i64,
    status: &str,
    error: Option<&str>,
) -> Result<()> {
    let conn = db.connection();

    conn.execute(
        "UPDATE files
         SET status = ?1,
             processing_error = ?2,
             updated_at = CURRENT_TIMESTAMP
         WHERE id = ?3",
        rusqlite::params![status, error, file_id],
    )?;

    Ok(())
}

pub fn get_genotype_metadata(db: &BioVaultDb, file_id: i64) -> Result<Option<GenotypeMetadata>> {
    let conn = db.connection();

    let result = conn.query_row(
        "SELECT source, grch_version, row_count, chromosome_count, inferred_sex
         FROM genotype_metadata
         WHERE file_id = ?1",
        params![file_id],
        |row| {
            Ok(GenotypeMetadata {
                data_type: "Genotype".to_string(),
                source: row.get(0)?,
                grch_version: row.get(1)?,
                row_count: row.get(2)?,
                chromosome_count: row.get(3)?,
                inferred_sex: row.get(4)?,
            })
        },
    );

    match result {
        Ok(meta) => Ok(Some(meta)),
        Err(rusqlite::Error::QueryReturnedNoRows) => Ok(None),
        Err(e) => Err(e.into()),
    }
}

pub fn update_file_from_queue(
    db: &BioVaultDb,
    file_id: i64,
    hash: &str,
    metadata: Option<&GenotypeMetadata>,
) -> Result<()> {
    let conn = db.connection();

    // Update file hash and data type
    if let Some(meta) = metadata {
        conn.execute(
            "UPDATE files
             SET file_hash = ?1,
                 data_type = ?2,
                 status = 'complete',
                 updated_at = CURRENT_TIMESTAMP
             WHERE id = ?3",
            rusqlite::params![hash, meta.data_type, file_id],
        )?;

        // Insert or update genotype metadata if this is a genotype file
        if meta.data_type == "Genotype" {
            conn.execute(
                "INSERT OR REPLACE INTO genotype_metadata
                 (file_id, source, grch_version, row_count, chromosome_count, inferred_sex, created_at, updated_at)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
                rusqlite::params![
                    file_id,
                    meta.source,
                    meta.grch_version,
                    meta.row_count,
                    meta.chromosome_count,
                    meta.inferred_sex,
                ],
            )?;

            // Update participant's inferred sex if we have it
            if let Some(ref sex) = meta.inferred_sex {
                conn.execute(
                    "UPDATE participants
                     SET inferred_sex = ?1
                     WHERE id = (SELECT participant_id FROM files WHERE id = ?2)
                       AND participant_id IS NOT NULL",
                    rusqlite::params![sex, file_id],
                )?;
            }
        }
    } else {
        // Just update the hash and mark as complete
        conn.execute(
            "UPDATE files
             SET file_hash = ?1,
                 status = 'complete',
                 updated_at = CURRENT_TIMESTAMP
             WHERE id = ?2",
            rusqlite::params![hash, file_id],
        )?;
    }

    Ok(())
}

#[cfg(test)]
/// Calculate SHA256 hash of file (test-only)
fn calculate_file_hash(path: &str) -> Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0; 8192];

    loop {
        let bytes_read = file.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_calculate_file_hash() {
        let temp = TempDir::new().unwrap();
        let file_path = temp.path().join("test.txt");
        fs::write(&file_path, b"hello world").unwrap();

        let hash = calculate_file_hash(file_path.to_str().unwrap()).unwrap();
        // SHA256 of "hello world"
        assert_eq!(
            hash,
            "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        );
    }

    #[test]
    fn test_scan_directory() {
        let temp = TempDir::new().unwrap();

        // Create test files
        fs::write(temp.path().join("file1.txt"), b"content1").unwrap();
        fs::write(temp.path().join("file2.txt"), b"content2").unwrap();
        fs::write(temp.path().join("file3.vcf"), b"content3").unwrap();

        // Scan all files
        let result = scan(temp.path().to_str().unwrap(), None, false).unwrap();
        assert_eq!(result.total_files, 3);
        assert_eq!(result.extensions.len(), 2); // .txt and .vcf

        // Scan only .txt files
        let result = scan(temp.path().to_str().unwrap(), Some(".txt"), false).unwrap();
        assert_eq!(result.total_files, 2);
    }

    #[test]
    fn test_extract_id_from_pattern_parent() {
        // Test {parent} token - extract parent directory name
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{parent}").unwrap();
        assert_eq!(id, "huE922FC");

        // Test {dirname} alias
        let id = extract_id_from_pattern(file_path, "{dirname}").unwrap();
        assert_eq!(id, "huE922FC");

        // Test {dir} alias
        let id = extract_id_from_pattern(file_path, "{dir}").unwrap();
        assert_eq!(id, "huE922FC");
    }

    #[test]
    fn test_extract_id_from_pattern_filename() {
        // Test {filename} token - filename without extension
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{filename}").unwrap();
        assert_eq!(id, "AncestryDNA");

        // Test {basename} token - filename with extension
        let id = extract_id_from_pattern(file_path, "{basename}").unwrap();
        assert_eq!(id, "AncestryDNA.txt");
    }

    #[test]
    fn test_extract_id_from_pattern_legacy() {
        // Test legacy {id}/* pattern
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{id}/*").unwrap();
        assert_eq!(id, "huE922FC");
    }

    #[test]
    fn test_extract_id_from_pattern_filename_pattern() {
        // Test {id} in filename pattern
        let file_path = "/data/files/123456_sample.txt";
        let id = extract_id_from_pattern(file_path, "{id}_*").unwrap();
        assert_eq!(id, "123456");

        // Test alphanumeric IDs
        let file_path = "/data/files/ABC123_sample.txt";
        let id = extract_id_from_pattern(file_path, "{id}_*").unwrap();
        assert_eq!(id, "ABC123");
    }
}
