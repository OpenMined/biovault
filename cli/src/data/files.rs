use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use regex::Regex;
use rusqlite::params;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

use super::BioVaultDb;

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
            let filename = Path::new(&file_info.path)
                .file_name()
                .and_then(|n| n.to_str());

            if let Some(fname) = filename {
                extract_id_from_pattern(fname, pat)
            } else {
                None
            }
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

/// Import single file into database with optional participant assignment
fn import_file_with_participant(
    db: &BioVaultDb,
    file_info: &FileInfo,
    participant_id: Option<&str>,
) -> Result<Option<i64>> {
    let file_path = &file_info.path;

    // Calculate file hash
    let file_hash = calculate_file_hash(file_path)?;

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

/// Calculate SHA256 hash of file
fn calculate_file_hash(path: &str) -> Result<String> {
    let mut file =
        fs::File::open(path).with_context(|| format!("Failed to open file: {}", path))?;

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

/// Get file by ID
pub fn get_file_by_id(db: &BioVaultDb, file_id: i64) -> Result<Option<FileRecord>> {
    let result = db.conn.query_row(
        "SELECT f.id, f.file_path, f.file_hash, f.file_type, f.file_size,
                p.participant_id, p.participant_id, f.created_at, f.updated_at
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
         WHERE f.id = ?1",
        params![file_id],
        |row| {
            Ok(FileRecord {
                id: row.get(0)?,
                file_path: row.get(1)?,
                file_hash: row.get(2)?,
                file_type: row.get(3)?,
                file_size: row.get::<_, Option<i64>>(4)?.map(|s| s as u64),
                participant_id: row.get(5)?,
                participant_name: row.get(6)?,
                created_at: row.get(7)?,
                updated_at: row.get(8)?,
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
                p.participant_id, p.participant_id, f.created_at, f.updated_at
         FROM files f
         LEFT JOIN participants p ON f.participant_id = p.id
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
            participant_id: row.get(5)?,
            participant_name: row.get(6)?,
            created_at: row.get(7)?,
            updated_at: row.get(8)?,
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

    Ok(SuggestPatternsResult {
        suggestions,
        sample_files: filenames.into_iter().take(5).collect(),
    })
}

/// Extract participant ID from filename using pattern (public for CLI use)
pub fn extract_id_from_pattern(filename: &str, pattern: &str) -> Option<String> {
    if !pattern.contains("{id}") {
        return None;
    }

    let regex_pattern = pattern
        .replace(".", "\\.")
        .replace("*", ".*")
        .replace("{id}", r"(\d+)");

    let re = Regex::new(&regex_pattern).ok()?;
    let captures = re.captures(filename)?;
    captures.get(1).map(|m| m.as_str().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
