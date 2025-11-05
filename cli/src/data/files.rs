use anyhow::{anyhow, bail, Result};
use chrono::{DateTime, Utc};
use regex::Regex;
use rusqlite::params;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
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
    pub regex_pattern: String,
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

    let compiled_pattern = pattern.map(compile_pattern).transpose()?;

    for file_info in scan_result.files {
        let extracted_id = compiled_pattern
            .as_ref()
            .and_then(|matcher| matcher.extract(&file_info.path));

        match import_file_with_participant(db, &file_info, extracted_id.as_deref()) {
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
pub fn import_from_csv(
    db: &BioVaultDb,
    csv_imports: Vec<CsvFileImport>,
    run_analysis: bool,
) -> Result<ImportResult> {
    let mut imported = 0;
    let mut skipped = 0;
    let mut errors = Vec::new();
    let mut imported_file_ids = Vec::new();

    for csv_row in csv_imports {
        match import_file_with_metadata(db, &csv_row, run_analysis) {
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
fn import_file_with_metadata(
    db: &BioVaultDb,
    csv_row: &CsvFileImport,
    run_analysis: bool,
) -> Result<Option<i64>> {
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

    // Build metadata snapshot from CSV values
    let mut metadata = GenotypeMetadata {
        data_type: csv_row
            .data_type
            .clone()
            .unwrap_or_else(|| "Unknown".to_string()),
        source: csv_row.source.clone(),
        grch_version: csv_row.grch_version.clone(),
        row_count: csv_row.row_count,
        chromosome_count: csv_row.chromosome_count,
        inferred_sex: csv_row.inferred_sex.clone(),
    };

    // Detect genotype metadata when not provided or still unknown
    if run_analysis {
        match crate::data::detect_genotype_metadata(file_path) {
            Ok(detected) => {
                if metadata.data_type == "Unknown" {
                    metadata.data_type = detected.data_type;
                }
                if metadata.source.is_none() {
                    metadata.source = detected.source;
                }
                if metadata.grch_version.is_none() {
                    metadata.grch_version = detected.grch_version;
                }
            }
            Err(e) => {
                eprintln!(
                    "Warning: Failed to detect metadata for {}: {}",
                    file_path, e
                );
            }
        }

        if metadata.data_type == "Genotype"
            && (metadata.row_count.is_none()
                || metadata.chromosome_count.is_none()
                || metadata.inferred_sex.is_none())
        {
            match crate::data::analyze_genotype_file(file_path) {
                Ok(analysis) => {
                    if metadata.row_count.is_none() {
                        metadata.row_count = analysis.row_count;
                    }
                    if metadata.chromosome_count.is_none() {
                        metadata.chromosome_count = analysis.chromosome_count;
                    }
                    if metadata.inferred_sex.is_none() {
                        metadata.inferred_sex = analysis.inferred_sex;
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to analyze {}: {}", file_path, e);
                }
            }
        }
    }

    // Persist metadata using shared queue update helper so genotype_metadata stays in sync
    update_file_from_queue(db, file_id, &file_hash, Some(&metadata))?;

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
    let scan_result = scan(path, extension, recursive)?;

    if scan_result.files.is_empty() {
        return Ok(SuggestPatternsResult {
            suggestions: vec![],
            sample_files: vec![],
        });
    }

    let entries = collect_file_entries(&scan_result.files);

    if entries.is_empty() {
        return Ok(SuggestPatternsResult {
            suggestions: vec![],
            sample_files: vec![],
        });
    }

    let mut candidates = Vec::new();
    candidates.extend(suggest_directory_candidates(&entries));
    candidates.extend(suggest_leading_numeric_candidates(&entries));
    candidates.extend(suggest_alpha_prefix_numeric_candidates(&entries));
    candidates.extend(suggest_generic_numeric_candidates(&entries));

    let suggestions = coalesce_candidates(candidates)
        .into_iter()
        .map(|candidate| PatternSuggestion {
            pattern: candidate.pattern,
            regex_pattern: candidate.regex_pattern,
            description: candidate.description,
            example: candidate.example,
            sample_extractions: candidate.sample_extractions,
        })
        .collect();

    let sample_files = entries
        .iter()
        .filter_map(|entry| {
            if entry.file_name.is_empty() {
                None
            } else {
                Some(entry.file_name.clone())
            }
        })
        .take(5)
        .collect();

    Ok(SuggestPatternsResult {
        suggestions,
        sample_files,
    })
}

#[derive(Clone, Debug)]
struct FileEntry {
    file_name: String,
    stem: Option<String>,
    parent: Option<String>,
}

#[derive(Clone, Debug)]
struct PatternCandidate {
    pattern: String,
    regex_pattern: String,
    description: String,
    example: String,
    sample_extractions: Vec<(String, String)>,
    coverage: usize,
}

#[derive(Clone, Copy, Debug)]
enum PatternScope {
    Path,
    Parent,
    Basename,
    Stem,
}

#[derive(Clone, Debug)]
pub struct CompiledPattern {
    scope: PatternScope,
    regex: Regex,
}

impl CompiledPattern {
    fn new(scope: PatternScope, regex: Regex) -> Self {
        Self { scope, regex }
    }

    pub fn extract(&self, filepath: &str) -> Option<String> {
        let target = match self.scope {
            PatternScope::Path => Some(filepath.to_string()),
            PatternScope::Parent => Path::new(filepath)
                .parent()
                .and_then(|p| p.file_name())
                .and_then(|n| n.to_str())
                .map(|s| s.to_string()),
            PatternScope::Basename => Path::new(filepath)
                .file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.to_string()),
            PatternScope::Stem => Path::new(filepath)
                .file_stem()
                .and_then(|n| n.to_str())
                .map(|s| s.to_string()),
        }?;

        let caps = self.regex.captures(&target)?;
        if let Some(id) = caps.name("id") {
            return Some(id.as_str().to_string());
        }

        for idx in 1..caps.len() {
            if let Some(mat) = caps.get(idx) {
                return Some(mat.as_str().to_string());
            }
        }

        None
    }
}

/// Compile a participant ID pattern into a reusable matcher.
///
/// Patterns may be expressed as raw regular expressions, structured templates
/// (`{scope:template}`), or legacy glob strings containing `{id}`. Structured
/// templates support scopes such as `path`, `parent`, `filename`, `basename`,
/// and `stem`, with glob-style wildcards (`*`, `**`, `?`).
pub fn compile_pattern(pattern: &str) -> Result<CompiledPattern> {
    compile_pattern_internal(pattern)
}

/// Convenience helper that compiles a pattern and applies it to a single path.
pub fn extract_id_from_pattern(filepath: &str, pattern: &str) -> Result<Option<String>> {
    let compiled = compile_pattern_internal(pattern)?;
    Ok(compiled.extract(filepath))
}

fn compile_pattern_internal(pattern: &str) -> Result<CompiledPattern> {
    let trimmed = pattern.trim();
    if trimmed.is_empty() {
        bail!("Pattern cannot be empty");
    }

    match trimmed {
        "{parent}" | "{dirname}" | "{dir}" | "{id}/*" => {
            let regex = Regex::new(r"^(?P<id>[A-Za-z0-9._-]+)$").unwrap();
            return Ok(CompiledPattern::new(PatternScope::Parent, regex));
        }
        "{filename}" => {
            let regex = Regex::new(r"^(?P<id>.+)$").unwrap();
            return Ok(CompiledPattern::new(PatternScope::Stem, regex));
        }
        "{basename}" => {
            let regex = Regex::new(r"^(?P<id>.+)$").unwrap();
            return Ok(CompiledPattern::new(PatternScope::Basename, regex));
        }
        _ => {}
    }

    if let Some((scope_label, template)) = parse_structured_template(trimmed) {
        let scope = parse_scope(&scope_label)?;
        let regex_source = build_regex_from_template(&template, scope)?;
        let regex = Regex::new(&regex_source).map_err(|err| {
            anyhow!(
                "Invalid regex derived from template '{}': {}",
                template,
                err
            )
        })?;
        ensure_has_capture(&regex)?;
        return Ok(CompiledPattern::new(scope, regex));
    }

    if trimmed.contains("{id}") {
        let scope = if trimmed.contains('/') {
            PatternScope::Path
        } else {
            PatternScope::Basename
        };
        let regex_source = build_regex_from_template(trimmed, scope)?;
        let regex = Regex::new(&regex_source)
            .map_err(|err| anyhow!("Invalid regex derived from pattern '{}': {}", trimmed, err))?;
        ensure_has_capture(&regex)?;
        return Ok(CompiledPattern::new(scope, regex));
    }

    let regex =
        Regex::new(trimmed).map_err(|err| anyhow!("Invalid regex '{}': {}", trimmed, err))?;
    ensure_has_capture(&regex)?;
    Ok(CompiledPattern::new(PatternScope::Path, regex))
}

fn ensure_has_capture(regex: &Regex) -> Result<()> {
    let mut has_capture = false;
    for name in regex.capture_names().flatten() {
        if name == "id" {
            has_capture = true;
            break;
        }
    }
    if !has_capture && regex.captures_len() > 1 {
        has_capture = true;
    }

    if !has_capture {
        bail!("Regex must contain at least one capture group (use (?P<id>...))");
    }

    Ok(())
}

fn parse_structured_template(pattern: &str) -> Option<(String, String)> {
    if !(pattern.starts_with('{') && pattern.ends_with('}')) {
        return None;
    }

    let inner = &pattern[1..pattern.len() - 1];
    let mut parts = inner.splitn(2, ':');
    let scope = parts.next()?.trim();
    let template = parts.next()?.trim();

    if scope.is_empty() || template.is_empty() {
        return None;
    }

    Some((scope.to_string(), template.to_string()))
}

fn parse_scope(label: &str) -> Result<PatternScope> {
    match label.to_lowercase().as_str() {
        "path" | "full" => Ok(PatternScope::Path),
        "parent" | "dir" | "dirname" | "folder" | "directory" => Ok(PatternScope::Parent),
        "filename" | "basename" => Ok(PatternScope::Basename),
        "stem" | "name" => Ok(PatternScope::Stem),
        other => bail!(
            "Unknown pattern scope '{}'. Expected one of path, parent, filename, stem",
            other
        ),
    }
}

fn build_regex_from_template(template: &str, scope: PatternScope) -> Result<String> {
    if !template.contains("{id}") {
        bail!("Template must include {{id}} placeholder");
    }

    let mut regex = String::from("^");
    let mut chars = template.chars().peekable();
    let mut found_id = false;

    while let Some(ch) = chars.next() {
        match ch {
            '*' => {
                if chars.peek() == Some(&'*') {
                    chars.next();
                    regex.push_str(".*");
                } else {
                    match scope {
                        PatternScope::Path => regex.push_str("[^/]*"),
                        _ => regex.push_str(".*"),
                    }
                }
            }
            '?' => regex.push('.'),
            '{' => {
                let mut lookahead = String::new();
                while let Some(&next_ch) = chars.peek() {
                    if next_ch == '}' {
                        chars.next();
                        break;
                    }
                    lookahead.push(next_ch);
                    chars.next();
                }

                if lookahead == "id" {
                    if found_id {
                        bail!("Template may contain only one {{id}} placeholder");
                    }
                    found_id = true;
                    let next_char = chars.peek().copied();
                    let base_class = match scope {
                        PatternScope::Path => "[^/]",
                        _ => match next_char {
                            Some('_') => "[^_]",
                            Some('-') => "[^-]",
                            Some('.') => "[^.]",
                            Some(' ') => "[^ ]",
                            _ => "[A-Za-z0-9._-]",
                        },
                    };

                    let capture = format!("(?P<id>{}+)", base_class);
                    regex.push_str(&capture);
                } else {
                    regex.push_str(&escape_literal(&format!("{{{}}}", lookahead)));
                }
            }
            '.' | '+' | '^' | '$' | '(' | ')' | '[' | ']' | '|' | '\\' => {
                regex.push('\\');
                regex.push(ch);
            }
            _ => regex.push(ch),
        }
    }

    if !found_id {
        bail!("Template must include {{id}} placeholder");
    }

    regex.push('$');
    Ok(regex)
}

fn escape_literal(text: &str) -> String {
    let mut escaped = String::new();
    for ch in text.chars() {
        match ch {
            '.' | '+' | '^' | '$' | '(' | ')' | '[' | ']' | '{' | '}' | '|' | '\\' => {
                escaped.push('\\');
                escaped.push(ch);
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
enum TokenClass {
    Numeric,
    Hex,
    Alpha,
    Alphanumeric,
    Other,
}

fn collect_file_entries(files: &[FileInfo]) -> Vec<FileEntry> {
    files
        .iter()
        .filter_map(|info| {
            let path = Path::new(&info.path);
            let file_name = path.file_name()?.to_str()?.to_string();
            if file_name.starts_with('.') {
                return None;
            }
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let parent = path
                .parent()
                .and_then(|p| p.file_name())
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());

            Some(FileEntry {
                file_name,
                stem,
                parent,
            })
        })
        .collect()
}

fn coalesce_candidates(mut candidates: Vec<PatternCandidate>) -> Vec<PatternCandidate> {
    if candidates.is_empty() {
        return candidates;
    }

    let mut map: HashMap<String, PatternCandidate> = HashMap::new();

    for candidate in candidates.drain(..) {
        map.entry(candidate.pattern.clone())
            .and_modify(|existing| {
                if candidate.coverage > existing.coverage {
                    existing.coverage = candidate.coverage;
                    existing.description = candidate.description.clone();
                    existing.example = candidate.example.clone();
                    existing.regex_pattern = candidate.regex_pattern.clone();
                }

                for sample in &candidate.sample_extractions {
                    if !existing.sample_extractions.contains(sample)
                        && existing.sample_extractions.len() < 5
                    {
                        existing.sample_extractions.push(sample.clone());
                    }
                }
            })
            .or_insert(candidate);
    }

    let mut deduped: Vec<PatternCandidate> = map.into_values().collect();
    deduped.sort_by(|a, b| {
        b.coverage
            .cmp(&a.coverage)
            .then_with(|| a.pattern.cmp(&b.pattern))
    });
    deduped
}

fn suggest_directory_candidates(entries: &[FileEntry]) -> Vec<PatternCandidate> {
    let mut parent_map: HashMap<&str, Vec<&FileEntry>> = HashMap::new();

    for entry in entries {
        if let Some(parent) = entry.parent.as_deref() {
            if parent.starts_with('.') {
                continue;
            }
            parent_map.entry(parent).or_default().push(entry);
        }
    }

    if parent_map.len() < 2 {
        return Vec::new();
    }

    let id_like_parents: Vec<(&str, &Vec<&FileEntry>)> = parent_map
        .iter()
        .filter_map(|(name, files)| {
            if is_id_like_parent(name) {
                Some((*name, files))
            } else {
                None
            }
        })
        .collect();

    if id_like_parents.len() < 2 {
        return Vec::new();
    }

    let id_like_parent_ratio = id_like_parents.len() as f64 / parent_map.len() as f64;
    if parent_map.len() >= 4 && id_like_parent_ratio < 0.6 {
        return Vec::new();
    }

    let id_like_coverage: usize = id_like_parents.iter().map(|(_, files)| files.len()).sum();
    if id_like_coverage == 0 {
        return Vec::new();
    }

    let mut candidates = Vec::new();

    // Generic parent directory suggestion limited to ID-like parent names
    let mut seen = HashSet::new();
    let mut sample_extractions = Vec::new();
    let mut example = String::new();
    for (parent, files) in &id_like_parents {
        if seen.insert(*parent) {
            if let Some(entry) = files.first() {
                let source = format!("{}/{}", parent, entry.file_name);
                if example.is_empty() {
                    example = source.clone();
                }
                sample_extractions.push((source, (*parent).to_string()));
                if sample_extractions.len() >= 3 {
                    break;
                }
            }
        }
    }

    if !sample_extractions.is_empty() {
        candidates.push(PatternCandidate {
            pattern: "{parent:{id}}".to_string(),
            regex_pattern: r".*/(?P<id>[A-Za-z0-9._-]+)/[^/]+$".to_string(),
            description: "Use parent directory names as participant IDs".to_string(),
            example: example.clone(),
            sample_extractions,
            coverage: id_like_coverage,
        });
    }

    // Look for shared prefixes in the ID-like parent directory names
    let parent_names: Vec<&str> = id_like_parents.iter().map(|(name, _)| *name).collect();
    if parent_names.len() < 2 {
        return candidates;
    }

    let prefix = longest_common_prefix(&parent_names);

    if prefix.len() >= 2 {
        let mut remainders: Vec<String> = Vec::new();
        for parent in &parent_names {
            if let Some(rest) = parent.strip_prefix(&prefix) {
                if rest.is_empty() {
                    remainders.clear();
                    break;
                }
                remainders.push(rest.to_string());
            } else {
                remainders.clear();
                break;
            }
        }

        if !remainders.is_empty() {
            let remainder_refs: Vec<&str> = remainders.iter().map(|r| r.as_str()).collect();
            let class = classify_token_slice(&remainder_refs);
            let unique_remainders: HashSet<&str> = remainder_refs.iter().copied().collect();
            let consistent_length = remainders
                .iter()
                .map(|r| r.len())
                .collect::<HashSet<_>>()
                .len()
                == 1;

            if matches!(
                class,
                TokenClass::Numeric | TokenClass::Hex | TokenClass::Alphanumeric
            ) && unique_remainders.len() >= 2
                && consistent_length
            {
                let mut sample_extractions = Vec::new();
                let mut example = String::new();
                let mut seen_dirs = HashSet::new();

                for (parent, files) in &id_like_parents {
                    if parent.starts_with(&prefix) && seen_dirs.insert(*parent) {
                        if let Some(entry) = files.first() {
                            let remainder = parent[prefix.len()..].to_string();
                            let source = format!("{}/{}", parent, entry.file_name);
                            if example.is_empty() {
                                example = source.clone();
                            }
                            sample_extractions.push((source, remainder));
                            if sample_extractions.len() >= 3 {
                                break;
                            }
                        }
                    }
                }

                if !sample_extractions.is_empty() {
                    let coverage = id_like_parents
                        .iter()
                        .filter(|(name, _)| name.starts_with(&prefix))
                        .map(|(_, files)| files.len())
                        .sum();

                    if coverage > 0 {
                        let descriptor = match class {
                            TokenClass::Numeric => "digits",
                            TokenClass::Hex => "hexadecimal IDs",
                            TokenClass::Alpha => "letters",
                            TokenClass::Alphanumeric => "alphanumeric IDs",
                            TokenClass::Other => "IDs",
                        };

                        let template_prefix = sanitize_glob_fragment(&prefix);
                        let remainder_len = remainders.first().map(|r| r.len()).unwrap_or(1);
                        let capture_pattern = match class {
                            TokenClass::Numeric => format!("\\d{{{}}}", remainder_len),
                            TokenClass::Hex => format!("[A-Fa-f0-9]{{{}}}", remainder_len),
                            TokenClass::Alphanumeric => format!("[A-Za-z0-9]{{{}}}", remainder_len),
                            TokenClass::Alpha => format!("[A-Za-z]{{{}}}", remainder_len),
                            TokenClass::Other => format!("[^/]{{{}}}", remainder_len),
                        };

                        let regex_pattern = format!(
                            r".*/(?P<id>{}{})/[^/]+$",
                            escape_literal(&prefix),
                            capture_pattern
                        );

                        candidates.push(PatternCandidate {
                            pattern: format!("{{parent:{}{{id}}}}", template_prefix),
                            regex_pattern,
                            description: format!(
                                "Directories with prefix '{}' followed by {}",
                                prefix, descriptor
                            ),
                            example,
                            sample_extractions,
                            coverage,
                        });
                    }
                }
            }
        }
    }

    candidates
}

fn suggest_leading_numeric_candidates(entries: &[FileEntry]) -> Vec<PatternCandidate> {
    let mut records = Vec::new();
    // Match 1-8 digits at the start of the filename (to catch 01, 02, 001, 000000, 12345678, etc.)
    let re = Regex::new(r"^(\d{1,8})").unwrap();

    for entry in entries {
        if let Some(stem) = &entry.stem {
            if let Some(mat) = re.captures(stem) {
                let digits = mat.get(1).unwrap().as_str().to_string();
                let next_char = stem.chars().nth(digits.len());
                let delimiter = match next_char {
                    Some('_') | Some('-') | Some('.') => Some(next_char.unwrap()),
                    _ => None,
                };

                records.push((entry.clone(), digits, delimiter));
            }
        }
    }

    // Lower threshold to 2 files instead of 3 to be more permissive
    if records.len() < 2 {
        return Vec::new();
    }

    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for (_, id, _) in &records {
        *length_counts.entry(id.len()).or_insert(0) += 1;
    }

    let (dominant_len, dominant_count) = length_counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .unwrap();

    if dominant_count < 2 {
        return Vec::new();
    }

    // Lower coverage threshold to 60% to be more permissive
    if dominant_count * 100 / records.len() < 60 {
        return Vec::new();
    }

    records.retain(|(_, id, _)| id.len() == dominant_len);

    let unique_ids: HashSet<&str> = records.iter().map(|(_, id, _)| id.as_str()).collect();
    if unique_ids.len() < 2 {
        return Vec::new();
    }

    let mut delimiter_counts: HashMap<char, usize> = HashMap::new();
    for (_, _, delimiter) in &records {
        if let Some(d) = delimiter {
            *delimiter_counts.entry(*d).or_insert(0) += 1;
        }
    }

    let chosen_delimiter = delimiter_counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(delim, _)| delim);

    let mut sample_extractions = Vec::new();
    let mut example = String::new();
    for (entry, id, _) in &records {
        sample_extractions.push((entry.file_name.clone(), id.clone()));
        if example.is_empty() {
            example = entry.file_name.clone();
        }
        if sample_extractions.len() >= 3 {
            break;
        }
    }

    let template = match chosen_delimiter {
        Some(d) => format!("{{id}}{}*", d),
        None => "{id}*".to_string(),
    };
    let pattern = format!("{{stem:{}}}", template);

    let delimiter_desc = match chosen_delimiter {
        Some('_') => "underscore",
        Some('-') => "dash",
        Some('.') => "dot",
        None => "", // No delimiter description used
        Some(_) => "delimiter",
    };

    let base_description = format!("Leading numeric ID ({} digits)", dominant_len);
    let description = if delimiter_desc.is_empty() {
        base_description.clone()
    } else {
        format!("{} before {}", base_description, delimiter_desc)
    };

    let regex_pattern = match chosen_delimiter {
        Some(d) => format!(
            "^(?P<id>\\d{{{}}}){}.*$",
            dominant_len,
            escape_literal(&d.to_string())
        ),
        None => format!("^(?P<id>\\d{{{}}}).*$", dominant_len),
    };

    vec![PatternCandidate {
        pattern,
        regex_pattern,
        description,
        example,
        sample_extractions,
        coverage: records.len(),
    }]
}

fn suggest_alpha_prefix_numeric_candidates(entries: &[FileEntry]) -> Vec<PatternCandidate> {
    let re = Regex::new(r"^([A-Za-z]{2,}[_-])(\d{3,})").unwrap();
    let mut records: HashMap<String, Vec<(FileEntry, String)>> = HashMap::new();

    for entry in entries {
        if let Some(stem) = &entry.stem {
            if let Some(caps) = re.captures(stem) {
                let prefix = caps.get(1).unwrap().as_str().to_string();
                let id = caps.get(2).unwrap().as_str().to_string();
                records.entry(prefix).or_default().push((entry.clone(), id));
            }
        }
    }

    let mut candidates = Vec::new();

    for (prefix, matches) in records {
        if matches.len() < 3 {
            continue;
        }

        let unique_ids: HashSet<&str> = matches.iter().map(|(_, id)| id.as_str()).collect();
        if unique_ids.len() < 2 {
            continue;
        }

        let mut sample_extractions = Vec::new();
        let mut example = String::new();
        for (entry, id) in &matches {
            sample_extractions.push((entry.file_name.clone(), id.clone()));
            if example.is_empty() {
                example = entry.file_name.clone();
            }
            if sample_extractions.len() >= 3 {
                break;
            }
        }

        let mut lengths: Vec<usize> = matches.iter().map(|(_, id)| id.len()).collect();
        lengths.sort_unstable();
        let min_len = *lengths.first().unwrap_or(&1);
        let max_len = *lengths.last().unwrap_or(&min_len);
        let capture = if min_len == max_len {
            format!("\\d{{{}}}", min_len)
        } else {
            format!("\\d{{{},{}}}", min_len, max_len)
        };

        let regex_pattern = format!("^{}(?P<id>{}).*$", escape_literal(&prefix), capture);

        let template = format!("{}{{id}}*", prefix);
        let pattern = format!("{{stem:{}}}", template);

        candidates.push(PatternCandidate {
            pattern,
            regex_pattern,
            description: format!(
                "Prefix '{}' followed by digits",
                prefix.trim_end_matches(&['-', '_'][..])
            ),
            example,
            sample_extractions,
            coverage: matches.len(),
        });
    }

    candidates
}

fn suggest_generic_numeric_candidates(entries: &[FileEntry]) -> Vec<PatternCandidate> {
    let re = Regex::new(r"(\d{3,})").unwrap();
    let mut context_map: HashMap<(String, String), Vec<(FileEntry, String)>> = HashMap::new();

    for entry in entries {
        if let Some(stem) = &entry.stem {
            if let Some(mat) = re.find(stem) {
                let id = mat.as_str().to_string();
                let before = stem[..mat.start()].to_string();
                let after = stem[mat.end()..].to_string();
                context_map
                    .entry((before, after))
                    .or_default()
                    .push((entry.clone(), id));
            }
        }
    }

    let mut candidates = Vec::new();

    for ((before, after), matches) in context_map {
        if matches.len() < 3 {
            continue;
        }

        let unique_ids: HashSet<&str> = matches.iter().map(|(_, id)| id.as_str()).collect();
        if unique_ids.len() < 2 {
            continue;
        }

        let mut sample_extractions = Vec::new();
        let mut example = String::new();
        for (entry, id) in &matches {
            sample_extractions.push((entry.file_name.clone(), id.clone()));
            if example.is_empty() {
                example = entry.file_name.clone();
            }
            if sample_extractions.len() >= 3 {
                break;
            }
        }

        let template = format!(
            "{}{{id}}{}",
            sanitize_glob_fragment(&before),
            sanitize_glob_fragment(&after)
        );
        let pattern = format!("{{stem:{}}}", template);

        let mut lengths: Vec<usize> = matches.iter().map(|(_, id)| id.len()).collect();
        lengths.sort_unstable();
        let min_len = *lengths.first().unwrap_or(&3);
        let max_len = *lengths.last().unwrap_or(&min_len);
        let capture = if min_len == max_len {
            format!("\\d{{{}}}", min_len)
        } else {
            format!("\\d{{{},{}}}", min_len, max_len)
        };

        let regex_pattern = format!(
            "^{}(?P<id>{}){}$",
            escape_literal(&before),
            capture,
            escape_literal(&after)
        );

        candidates.push(PatternCandidate {
            pattern,
            regex_pattern,
            description: "Repeated numeric sequence".to_string(),
            example,
            sample_extractions,
            coverage: matches.len(),
        });
    }

    candidates
}

fn is_id_like_parent(name: &str) -> bool {
    if name.len() < 3 {
        return false;
    }

    matches!(
        classify_token(name),
        TokenClass::Numeric | TokenClass::Hex | TokenClass::Alphanumeric
    )
}

fn sanitize_glob_fragment(fragment: &str) -> String {
    if fragment.is_empty() {
        return String::new();
    }

    let mut sanitized = String::new();
    for ch in fragment.chars() {
        match ch {
            '*' | '?' | '{' | '}' | '[' | ']' => {
                sanitized.push('[');
                sanitized.push(ch);
                sanitized.push(']');
            }
            _ => sanitized.push(ch),
        }
    }
    sanitized
}

fn longest_common_prefix(strings: &[&str]) -> String {
    if strings.is_empty() {
        return String::new();
    }

    let first = strings[0];
    let mut end = first.len();

    for s in strings.iter().skip(1) {
        let mut idx = 0;
        let max_len = std::cmp::min(end, s.len());
        while idx < max_len && first.as_bytes()[idx] == s.as_bytes()[idx] {
            idx += 1;
        }
        end = idx;
        if end == 0 {
            break;
        }
    }

    first[..end].to_string()
}

fn classify_token(token: &str) -> TokenClass {
    if token.is_empty() {
        return TokenClass::Other;
    }

    let mut has_alpha = false;
    let mut has_digit = false;
    let mut has_other = false;

    for ch in token.chars() {
        match ch {
            '0'..='9' => has_digit = true,
            'a'..='z' | 'A'..='Z' => has_alpha = true,
            '-' | '_' => {}
            _ => {
                has_other = true;
                break;
            }
        }
    }

    if has_other {
        TokenClass::Other
    } else if has_alpha && has_digit {
        if token.chars().all(|c| c.is_ascii_hexdigit()) {
            TokenClass::Hex
        } else {
            TokenClass::Alphanumeric
        }
    } else if has_digit {
        TokenClass::Numeric
    } else if has_alpha {
        TokenClass::Alpha
    } else {
        TokenClass::Other
    }
}

fn classify_token_slice(tokens: &[&str]) -> TokenClass {
    let mut classes: HashSet<TokenClass> = HashSet::new();
    for token in tokens {
        classes.insert(classify_token(token));
    }

    if classes.len() == 1 {
        *classes.iter().next().unwrap()
    } else if classes.contains(&TokenClass::Other) {
        TokenClass::Other
    } else if classes.contains(&TokenClass::Alphanumeric)
        || classes.contains(&TokenClass::Hex)
        || (classes.contains(&TokenClass::Alpha) && classes.contains(&TokenClass::Numeric))
    {
        TokenClass::Alphanumeric
    } else if classes.contains(&TokenClass::Alpha) {
        TokenClass::Alpha
    } else if classes.contains(&TokenClass::Numeric) {
        TokenClass::Numeric
    } else {
        TokenClass::Other
    }
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
/// Also deletes participants that have no files left after deletion
pub fn delete_file(db: &BioVaultDb, file_id: i64) -> Result<()> {
    // First, get the participant ID for this file before deletion
    let participant_id: Option<i64> = db
        .conn
        .query_row(
            "SELECT participant_id FROM files WHERE id = ?1",
            params![file_id],
            |row| row.get(0),
        )
        .ok();

    let rows = db
        .conn
        .execute("DELETE FROM files WHERE id = ?1", params![file_id])?;

    if rows == 0 {
        anyhow::bail!("File with id {} not found", file_id);
    }

    // Clean up orphaned participant if this file had one and it now has no files
    if let Some(pid) = participant_id {
        let file_count: i64 = db.conn.query_row(
            "SELECT COUNT(*) FROM files WHERE participant_id = ?1",
            params![pid],
            |row| row.get(0),
        )?;

        if file_count == 0 {
            let deleted = db
                .conn
                .execute("DELETE FROM participants WHERE id = ?1", params![pid])?;

            if deleted > 0 {
                eprintln!(
                    "🧹 Cleaned up orphaned participant (id: {}) after deleting file",
                    pid
                );
            }
        }
    }

    Ok(())
}

/// Delete multiple file records from the catalog
/// Also deletes participants that have no files left after deletion
pub fn delete_files_bulk(db: &BioVaultDb, ids: &[i64]) -> Result<usize> {
    if ids.is_empty() {
        return Ok(0);
    }

    // First, get the participant IDs that will be affected before deletion
    let mut affected_participant_ids = Vec::new();
    {
        let placeholders = ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");
        let select_query = format!(
            "SELECT DISTINCT participant_id FROM files WHERE id IN ({}) AND participant_id IS NOT NULL",
            placeholders
        );
        let mut stmt = db.conn.prepare(&select_query)?;
        let participant_iter = stmt.query_map(rusqlite::params_from_iter(ids.iter()), |row| {
            Ok(row.get::<_, i64>(0)?)
        })?;

        for result in participant_iter {
            if let Ok(pid) = result {
                affected_participant_ids.push(pid);
            }
        }
    }

    // Delete the files
    let placeholders = ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");
    let delete_query = format!("DELETE FROM files WHERE id IN ({})", placeholders);
    let rows = db
        .conn
        .execute(&delete_query, rusqlite::params_from_iter(ids.iter()))?;

    // Clean up orphaned participants (participants with no files left)
    // Check each affected participant to see if they still have files
    let mut deleted_participants = 0;
    for participant_id in affected_participant_ids {
        let file_count: i64 = db.conn.query_row(
            "SELECT COUNT(*) FROM files WHERE participant_id = ?1",
            params![participant_id],
            |row| row.get(0),
        )?;

        if file_count == 0 {
            let deleted = db.conn.execute(
                "DELETE FROM participants WHERE id = ?1",
                params![participant_id],
            )?;
            deleted_participants += deleted;
        }
    }

    if deleted_participants > 0 {
        eprintln!(
            "🧹 Cleaned up {} orphaned participant(s) after deleting files",
            deleted_participants
        );
    }

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

        // Extract file type (extension)
        let file_type = Path::new(&file_info.file_path)
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| format!(".{}", e));

        // Insert file with status='pending' and placeholder hash
        let result = conn.execute(
            "INSERT INTO files (participant_id, file_path, file_hash, file_type, file_size, data_type, status, queue_added_at, created_at, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, 'pending', CURRENT_TIMESTAMP, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
            rusqlite::params![
                participant_id,
                &file_info.file_path,
                "pending", // Placeholder hash until processed
                file_type,
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

// Queue information support

#[derive(Debug, Serialize, Deserialize)]
pub struct QueueInfo {
    pub total_pending: usize,
    pub processing_count: usize,
    pub queue_position: Option<usize>, // Position of specific file if file_id provided
    pub currently_processing: Option<QueueFileInfo>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct QueueFileInfo {
    pub id: i64,
    pub file_path: String,
}

/// Get queue information - total pending count, processing count, and optionally queue position for a specific file
pub fn get_queue_info(db: &BioVaultDb, file_id: Option<i64>) -> Result<QueueInfo> {
    let conn = db.connection();

    // Get total pending count
    let total_pending: i64 = conn.query_row(
        "SELECT COUNT(*) FROM files WHERE status = 'pending'",
        [],
        |row| row.get::<_, i64>(0),
    )?;

    // Get processing count
    let processing_count: i64 = conn.query_row(
        "SELECT COUNT(*) FROM files WHERE status = 'processing'",
        [],
        |row| row.get::<_, i64>(0),
    )?;

    // Get currently processing file (if any)
    let currently_processing = conn
        .query_row(
            "SELECT id, file_path FROM files WHERE status = 'processing' ORDER BY updated_at ASC LIMIT 1",
            [],
            |row| {
                Ok(QueueFileInfo {
                    id: row.get(0)?,
                    file_path: row.get(1)?,
                })
            },
        )
        .ok();

    // Get queue position for specific file if provided
    let queue_position = if let Some(fid) = file_id {
        // First check if file is currently processing
        let is_processing: Result<i64, rusqlite::Error> = conn.query_row(
            "SELECT COUNT(*) FROM files WHERE status = 'processing' AND id = ?1",
            [fid],
            |row| row.get::<_, i64>(0),
        );

        if let Ok(count) = is_processing {
            if count > 0 {
                // Currently processing = position 0
                Some(0)
            } else {
                // Count how many files are ahead of this one in the queue
                let position: Result<i64, rusqlite::Error> = conn.query_row(
                    "SELECT COUNT(*) FROM files 
                     WHERE status = 'pending' 
                     AND queue_added_at < (
                         SELECT queue_added_at FROM files WHERE id = ?1 AND status = 'pending'
                     )",
                    [fid],
                    |row| row.get::<_, i64>(0),
                );

                position.ok().map(|p| (p + 1) as usize) // Add 1 since position is 1-indexed
            }
        } else {
            None
        }
    } else {
        None
    };

    Ok(QueueInfo {
        total_pending: total_pending as usize,
        processing_count: processing_count as usize,
        queue_position,
        currently_processing,
    })
}

/// Clear all pending and processing files from the queue (delete them)
/// This stops any ongoing imports and clears the queue
/// Returns the number of files deleted
pub fn clear_pending_queue(db: &BioVaultDb) -> Result<usize> {
    let conn = db.connection();

    // Get all pending and processing file IDs first to track what will be deleted
    let queue_ids: Vec<i64> = conn
        .prepare("SELECT id FROM files WHERE status IN ('pending', 'processing')")?
        .query_map([], |row| row.get::<_, i64>(0))?
        .collect::<std::result::Result<Vec<_>, _>>()?;

    if queue_ids.is_empty() {
        return Ok(0);
    }

    // Get participant IDs that will be affected
    let mut affected_participant_ids = Vec::new();
    {
        let placeholders = queue_ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");
        let select_query = format!(
            "SELECT DISTINCT participant_id FROM files WHERE id IN ({}) AND participant_id IS NOT NULL",
            placeholders
        );
        let mut stmt = conn.prepare(&select_query)?;
        let participant_iter = stmt
            .query_map(rusqlite::params_from_iter(queue_ids.iter()), |row| {
                Ok(row.get::<_, i64>(0)?)
            })?;

        for result in participant_iter {
            if let Ok(pid) = result {
                affected_participant_ids.push(pid);
            }
        }
    }

    // Delete all pending and processing files
    let deleted = conn.execute(
        "DELETE FROM files WHERE status IN ('pending', 'processing')",
        [],
    )?;

    // Clean up orphaned participants (participants with no files left)
    for participant_id in affected_participant_ids {
        let file_count: i64 = conn.query_row(
            "SELECT COUNT(*) FROM files WHERE participant_id = ?1",
            params![participant_id],
            |row| row.get(0),
        )?;

        if file_count == 0 {
            let _ = conn.execute(
                "DELETE FROM participants WHERE id = ?1",
                params![participant_id],
            )?;
        }
    }

    Ok(deleted)
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
    fn test_suggest_patterns_leading_numeric() {
        let temp = TempDir::new().unwrap();
        let filenames = [
            "000000_carika.txt",
            "103704_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt",
            "111442_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt",
            "117292_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt",
            "123364_X_X_GSAv3-DTC_GRCh38-07-01-2025.txt",
            "256789_Combined_Genome.txt",
            "356789_Eric_Uhden_Full_20110718111059.txt",
        ];

        for filename in filenames.iter() {
            fs::write(temp.path().join(filename), b"content").unwrap();
        }

        let result = suggest_patterns(temp.path().to_str().unwrap(), Some(".txt"), false).unwrap();

        assert!(!result
            .suggestions
            .iter()
            .any(|suggestion| suggestion.pattern == "{parent:{id}}"));

        assert!(result
            .suggestions
            .iter()
            .any(|suggestion| suggestion.pattern == "{stem:{id}_*}"));

        let lead_numeric = result
            .suggestions
            .iter()
            .find(|suggestion| suggestion.pattern == "{stem:{id}_*}")
            .unwrap();

        // Ensure we captured canonical six-digit IDs
        assert!(lead_numeric
            .sample_extractions
            .iter()
            .any(|(_, id)| id.len() == 6));
        assert_eq!(lead_numeric.regex_pattern, "^(?P<id>\\d{6})_.*$");
    }

    #[test]
    fn test_suggest_patterns_directory_prefix() {
        let temp = TempDir::new().unwrap();
        let dirs = [
            ("hu17DFDB", vec!["23andMe_Genotyping.txt"]),
            ("hu44DCFF", vec!["JKP001_genotypes.txt"]),
            ("hu2D53F2", vec!["hu2D53F2_20120421013417.txt"]),
            (
                "hu836D0A",
                vec!["genome_Maureen_Markov_Full_20100823192336.txt"],
            ),
            ("huB714CA", vec!["huB714CA_20110726215545.txt"]),
        ];

        for (dir, files) in dirs.iter() {
            let dir_path = temp.path().join(dir);
            fs::create_dir_all(&dir_path).unwrap();
            for file in files {
                fs::write(dir_path.join(file), b"content").unwrap();
            }
        }

        let result = suggest_patterns(temp.path().to_str().unwrap(), None, true).unwrap();

        assert!(result
            .suggestions
            .iter()
            .any(|suggestion| suggestion.pattern == "{parent:{id}}"));

        let parent_pattern = result
            .suggestions
            .iter()
            .find(|suggestion| suggestion.pattern == "{parent:{id}}")
            .unwrap();

        assert_eq!(
            parent_pattern.regex_pattern,
            r".*/(?P<id>[A-Za-z0-9._-]+)/[^/]+$"
        );

        assert!(result
            .suggestions
            .iter()
            .any(|suggestion| suggestion.pattern == "{parent:hu{id}}"));

        let dir_prefix = result
            .suggestions
            .iter()
            .find(|suggestion| suggestion.pattern == "{parent:hu{id}}")
            .unwrap();

        assert!(dir_prefix
            .sample_extractions
            .iter()
            .all(|(_, id)| id.len() == 6));
        assert_eq!(
            dir_prefix.regex_pattern,
            r".*/(?P<id>hu[A-Fa-f0-9]{6})/[^/]+$"
        );
    }

    #[test]
    fn test_extract_id_from_pattern_parent() {
        // Test {parent} token - extract parent directory name
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{parent:{id}}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "huE922FC");

        // Test {dirname} alias
        let id = extract_id_from_pattern(file_path, "{dirname}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "huE922FC");

        // Test {dir} alias
        let id = extract_id_from_pattern(file_path, "{dir}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "huE922FC");
    }

    #[test]
    fn test_extract_id_from_pattern_filename() {
        // Test {filename} token - filename without extension
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{filename}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "AncestryDNA");

        // Test {basename} token - filename with extension
        let id = extract_id_from_pattern(file_path, "{basename}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "AncestryDNA.txt");
    }

    #[test]
    fn test_extract_id_from_pattern_legacy() {
        // Test legacy {id}/* pattern
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, "{id}/*")
            .unwrap()
            .unwrap();
        assert_eq!(id, "huE922FC");
    }

    #[test]
    fn test_extract_id_from_pattern_filename_pattern() {
        // Test {id} in filename pattern
        let file_path = "/data/files/123456_sample.txt";
        let id = extract_id_from_pattern(file_path, "{stem:{id}_*}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "123456");

        // Test alphanumeric IDs
        let file_path = "/data/files/ABC123_sample.txt";
        let id = extract_id_from_pattern(file_path, "{stem:{id}_*}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "ABC123");
    }

    #[test]
    fn test_extract_id_with_structured_template() {
        let file_path = "/data/files/genome_Full_20120106210128.txt";
        let id = extract_id_from_pattern(file_path, "{basename:genome_Full_{id}.txt}")
            .unwrap()
            .unwrap();
        assert_eq!(id, "20120106210128");
    }

    #[test]
    fn test_extract_id_with_regex() {
        let file_path = "/data/genotype_files/huE922FC/AncestryDNA.txt";
        let id = extract_id_from_pattern(file_path, r"(?P<id>hu[0-9A-F]{6})")
            .unwrap()
            .unwrap();
        assert_eq!(id, "huE922FC");
    }
}
