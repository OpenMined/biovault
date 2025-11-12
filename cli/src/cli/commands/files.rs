use anyhow::Result;
use colored::Colorize;

use crate::data::{self, BioVaultDb, CliResponse};

pub async fn scan(
    path: String,
    extension: Option<String>,
    recursive: bool,
    format: String,
) -> Result<()> {
    let result = data::scan(&path, extension.as_deref(), recursive)?;

    if format == "json" {
        let response = CliResponse::new(&result);
        println!("{}", response.to_json()?);
    } else {
        // Table format
        println!("{}", format!("📊 Scan Results: {}", path).bold());
        println!();

        if !result.extensions.is_empty() {
            println!("{}", "Extensions Found:".bold());
            for ext_info in &result.extensions {
                let size_mb = ext_info.total_size as f64 / 1_048_576.0;
                println!(
                    "  {}  {} files    {:.1} MB",
                    ext_info.extension.cyan(),
                    ext_info.count,
                    size_mb
                );
            }
            println!();
        }

        println!(
            "Total: {} files",
            result.total_files.to_string().green().bold()
        );

        if result.total_files > 0 && result.total_files <= 10 {
            println!();
            println!("{}", "Files:".bold());
            for file in &result.files {
                let size_kb = file.size as f64 / 1024.0;
                println!("  {} ({:.1} KB)", file.path, size_kb);
            }
        }
    }

    Ok(())
}

pub async fn import(
    path: String,
    extension: Option<String>,
    recursive: bool,
    pattern: Option<String>,
    collection: Option<String>,
    collection_description: Option<String>,
    collection_var_name: Option<String>,
    dry_run: bool,
    non_interactive: bool,
    format: String,
) -> Result<()> {
    // First, do a scan to show what will be imported
    let scan_result = data::scan(&path, extension.as_deref(), recursive)?;

    if scan_result.files.is_empty() {
        println!("{}", "No files found to import.".yellow());
        return Ok(());
    }

    let compiled_pattern = pattern
        .as_ref()
        .map(|pat| data::compile_pattern(pat))
        .transpose()?;

    // Show preview
    if format != "json" {
        println!("{}", "📋 Import Preview:".bold());
        println!("  Path: {}", path);
        if let Some(ext) = &extension {
            println!("  Extension filter: {}", ext.cyan());
        }
        if let Some(pat) = &pattern {
            println!("  Pattern: {}", pat.green().bold());
        }
        println!(
            "  Files to import: {}",
            scan_result.total_files.to_string().cyan().bold()
        );

        // Show sample extractions if pattern provided
        if let Some(compiled) = &compiled_pattern {
            println!();
            println!("{}", "Sample participant ID extractions:".bold());
            for file_info in scan_result.files.iter().take(5) {
                let filename = std::path::Path::new(&file_info.path)
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(&file_info.path);

                if let Some(id) = compiled.extract(&file_info.path) {
                    println!(
                        "  {} → participant: {}",
                        filename.dimmed(),
                        id.cyan().bold()
                    );
                } else {
                    println!("  {} → {}", filename.dimmed(), "no match".red());
                }
            }
            if scan_result.files.len() > 5 {
                println!("  ... and {} more", scan_result.files.len() - 5);
            }
        }
        println!();
    }

    // Dry run - show what would happen and exit
    if dry_run {
        println!(
            "{}",
            "🔍 Dry run - no files will be imported.".yellow().bold()
        );
        println!(
            "{}",
            "Remove --dry-run flag to proceed with import.".dimmed()
        );
        return Ok(());
    }

    // Interactive confirmation (unless --non-interactive)
    if !non_interactive && format != "json" {
        use dialoguer::Confirm;

        let proceed = Confirm::new()
            .with_prompt("Proceed with import?")
            .default(true)
            .interact()?;

        if !proceed {
            println!("{}", "Import cancelled.".yellow());
            return Ok(());
        }
    }

    // Perform the actual import
    let db = BioVaultDb::new()?;
    let result = data::import(
        &db,
        &path,
        extension.as_deref(),
        recursive,
        pattern.as_deref(),
    )?;

    // Handle collection if specified
    let collection_id = if let Some(collection_name) = &collection {
        // Get or create collection
        let collection_record = match data::get_collection(&db, collection_name) {
            Ok(c) => c,
            Err(_) => {
                // Collection doesn't exist, create it
                let suggested_var_name = data::generate_variable_name(collection_name);
                let var_name = collection_var_name.unwrap_or(suggested_var_name);
                
                if format != "json" {
                    println!(
                        "{}",
                        format!("📦 Creating collection '{}'...", collection_name).dimmed()
                    );
                }
                
                data::create_collection(
                    &db,
                    collection_name.clone(),
                    collection_description.clone(),
                    Some(var_name),
                )?
            }
        };

        // Add imported files to collection
        let collection_var_name = collection_record.variable_name.clone();
        if !result.files.is_empty() {
            let file_ids: Vec<i64> = result.files.iter().map(|f| f.id).collect();
            let added = data::add_files_to_collection(&db, &collection_var_name, file_ids)?;
            
            if format != "json" && added > 0 {
                println!(
                    "{}",
                    format!("  ✓ Added {} file(s) to collection '{}'", added, collection_name)
                        .dimmed()
                );
            }
        }

        Some(collection_record.id)
    } else {
        None
    };

    if format == "json" {
        let mut json_result = serde_json::json!({
            "imported": result.imported,
            "skipped": result.skipped,
            "errors": result.errors,
            "files": result.files
        });
        if let Some(cid) = collection_id {
            json_result["collection_id"] = serde_json::json!(cid);
        }
        let response = CliResponse::new(&json_result);
        println!("{}", response.to_json()?);
    } else {
        // Table format
        println!();
        if result.imported > 0 {
            println!(
                "{}",
                format!("✓ Imported {} files", result.imported)
                    .green()
                    .bold()
            );
        }
        if result.skipped > 0 {
            println!("⊘ Skipped {} files (already imported)", result.skipped);
        }
        if !result.errors.is_empty() {
            println!("{}", format!("✗ {} errors:", result.errors.len()).red());
            for error in &result.errors {
                println!("  {}", error.red());
            }
        }

        if !result.files.is_empty() && result.files.len() <= 10 {
            println!();
            println!("{}", "Files imported:".bold());
            println!(
                "  {}  {}  {}  {}  {}",
                "ID".bold(),
                "Path".bold(),
                "Hash".bold(),
                "Type".bold(),
                "Participant".bold()
            );
            for file in &result.files {
                let hash_short = if file.file_hash.len() > 8 {
                    format!("{}...", &file.file_hash[..8])
                } else {
                    file.file_hash.clone()
                };
                let participant = file.participant_name.as_deref().unwrap_or("-");
                println!(
                    "  {}  {}  {}  {}  {}",
                    file.id,
                    file.file_path,
                    hash_short,
                    file.file_type.as_deref().unwrap_or("-"),
                    participant
                );
            }
        } else if result.imported > 10 {
            println!();
            println!(
                "{}",
                format!("✓ {} files imported successfully", result.imported).green()
            );
            println!(
                "{}",
                "Use 'bv files list' to see all imported files.".dimmed()
            );
        }
    }

    Ok(())
}

pub async fn list(
    extension: Option<String>,
    participant: Option<String>,
    unassigned: bool,
    limit: Option<usize>,
    format: String,
) -> Result<()> {
    let db = BioVaultDb::new()?;
    let files = data::list_files(
        &db,
        extension.as_deref(),
        participant.as_deref(),
        unassigned,
        limit,
    )?;

    if format == "json" {
        #[derive(serde::Serialize)]
        struct ListResponse {
            total: usize,
            files: Vec<data::FileRecord>,
        }

        let response = CliResponse::new(ListResponse {
            total: files.len(),
            files,
        });
        println!("{}", response.to_json()?);
    } else {
        // Table format
        if files.is_empty() {
            println!("{}", "No files found.".yellow());
            println!("Use 'bv files import' to add files.");
            return Ok(());
        }

        println!("{}", format!("Files ({})", files.len()).bold());
        println!();
        println!(
            "  {}  {}  {}  {}  {}  {}",
            "ID".bold(),
            "Path".bold(),
            "Hash".bold(),
            "Type".bold(),
            "Size".bold(),
            "Participant".bold()
        );

        for file in &files {
            let hash_short = if file.file_hash.len() > 8 {
                format!("{}...", &file.file_hash[..8])
            } else {
                file.file_hash.clone()
            };
            let participant = file.participant_name.as_deref().unwrap_or("-");
            let size = if let Some(s) = file.file_size {
                format!("{:.1} KB", s as f64 / 1024.0)
            } else {
                "-".to_string()
            };

            println!(
                "  {}  {}  {}  {}  {}  {}",
                file.id,
                file.file_path,
                hash_short,
                file.file_type.as_deref().unwrap_or("-"),
                size,
                participant
            );
        }
    }

    Ok(())
}

pub async fn info(file: String, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Try to parse as ID first, then as path
    let file_record = if let Ok(file_id) = file.parse::<i64>() {
        data::get_file_by_id(&db, file_id)?
    } else {
        // TODO: Implement get_file_by_path
        anyhow::bail!("File path lookup not yet implemented. Use file ID instead.");
    };

    if let Some(record) = file_record {
        if format == "json" {
            let response = CliResponse::new(&record);
            println!("{}", response.to_json()?);
        } else {
            // Table format
            println!("{}", format!("File: {}", record.id).bold());
            println!("  Path:        {}", record.file_path);
            println!("  Hash:        {}", record.file_hash);
            println!(
                "  Type:        {}",
                record.file_type.as_deref().unwrap_or("-")
            );
            if let Some(size) = record.file_size {
                println!("  Size:        {:.2} MB", size as f64 / 1_048_576.0);
            }
            println!(
                "  Participant: {}",
                record.participant_name.as_deref().unwrap_or("-")
            );
            println!("  Created:     {}", record.created_at);
            println!("  Updated:     {}", record.updated_at);
        }
    } else {
        println!("{}", format!("File not found: {}", file).red());
    }

    Ok(())
}

pub async fn suggest_patterns(
    path: String,
    extension: Option<String>,
    recursive: bool,
    format: String,
) -> Result<()> {
    let result = data::suggest_patterns(&path, extension.as_deref(), recursive)?;

    if format == "json" {
        let response = CliResponse::new(&result);
        println!("{}", response.to_json()?);
    } else {
        // Table format
        if result.suggestions.is_empty() {
            println!("{}", "No patterns detected.".yellow());
            println!("Files may not contain identifiable participant ID patterns.");
            return Ok(());
        }

        println!("{}", "🔍 Detected Patterns:".bold());
        println!();

        for (i, suggestion) in result.suggestions.iter().enumerate() {
            println!(
                "{}. {} - {}",
                (i + 1).to_string().cyan(),
                suggestion.pattern.green().bold(),
                suggestion.description
            );
            println!("   Regex: {}", suggestion.regex_pattern.cyan());
            println!("   Example: {}", suggestion.example.dimmed());

            if !suggestion.sample_extractions.is_empty() {
                println!("   Sample extractions:");
                for (filename, id) in &suggestion.sample_extractions {
                    println!(
                        "     {} → participant ID: {}",
                        filename.dimmed(),
                        id.cyan().bold()
                    );
                }
            }
            println!();
        }

        println!("{}", "Sample files analyzed:".bold());
        for file in &result.sample_files {
            println!("  • {}", file.dimmed());
        }
        println!();

        println!("{}", "💡 To test a pattern:".bold());
        println!("   Run: bv files test-pattern <pattern>");
        println!();
        println!("{}", "💡 To import with a pattern:".bold());
        println!("   (Coming in Phase 2)");
    }

    Ok(())
}

pub async fn delete(id: i64, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get file info before deletion for display
    let file_record = data::get_file_by_id(&db, id)?;

    let file = match file_record {
        Some(f) => f,
        None => {
            if format == "json" {
                let response = CliResponse::new(serde_json::json!({
                    "error": format!("File with id {} not found", id)
                }));
                println!("{}", response.to_json()?);
            } else {
                println!("{}", format!("File not found: {}", id).red());
            }
            return Ok(());
        }
    };

    if format != "json" {
        println!("Deleting file record: {}", file.file_path);
        if let Some(participant) = &file.participant_name {
            println!("  Participant: {}", participant.cyan());
        }
    }

    data::delete_file(&db, id)?;

    if format == "json" {
        #[derive(serde::Serialize)]
        struct DeleteResponse {
            deleted_id: i64,
            file_path: String,
        }

        let response = CliResponse::new(DeleteResponse {
            deleted_id: id,
            file_path: file.file_path,
        });
        println!("{}", response.to_json()?);
    } else {
        println!("{}", format!("✓ Deleted file record {}", id).green().bold());
    }

    Ok(())
}

pub async fn delete_bulk(ids: Vec<i64>, format: String) -> Result<()> {
    if ids.is_empty() {
        if format == "json" {
            let response = CliResponse::new(serde_json::json!({
                "deleted": 0,
                "errors": []
            }));
            println!("{}", response.to_json()?);
        }
        return Ok(());
    }

    let db = BioVaultDb::new()?;

    if format != "json" {
        println!(
            "{}",
            format!("Deleting {} file records...", ids.len()).bold()
        );
    }

    let mut deleted = 0;
    let mut errors = Vec::new();

    for id in &ids {
        match data::delete_file(&db, *id) {
            Ok(_) => {
                deleted += 1;
                if format != "json" {
                    println!("  {} Deleted file ID {}", "✓".green(), id);
                }
            }
            Err(e) => {
                errors.push(format!("Failed to delete file {}: {}", id, e));
                if format != "json" {
                    println!("  {} Failed to delete file ID {}: {}", "✗".red(), id, e);
                }
            }
        }
    }

    if format == "json" {
        let response = CliResponse::new(serde_json::json!({
            "deleted": deleted,
            "errors": errors
        }));
        println!("{}", response.to_json()?);
    } else {
        println!();
        println!(
            "{}",
            format!("✓ Deleted {} of {} file records", deleted, ids.len())
                .green()
                .bold()
        );
        if !errors.is_empty() {
            println!("{}", format!("✗ {} errors", errors.len()).red());
        }
    }

    Ok(())
}

pub async fn link(file_id: i64, participant: String, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;

    let updated_file = data::link_file_to_participant(&db, file_id, &participant)?;

    if format == "json" {
        let response = CliResponse::new(&updated_file);
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!("✓ Linked file {} to participant {}", file_id, participant)
                .green()
                .bold()
        );
        println!("  File: {}", updated_file.file_path);
        println!("  Participant: {}", participant.cyan());
    }

    Ok(())
}

pub async fn link_bulk(file_participant_json: String, format: String) -> Result<()> {
    use std::collections::HashMap;

    // Parse JSON input: {"file_path": "participant_id", ...}
    let file_participant_map: HashMap<String, String> =
        serde_json::from_str(&file_participant_json)
            .map_err(|e| anyhow::anyhow!("Failed to parse JSON: {}", e))?;

    let db = BioVaultDb::new()?;
    let updated = data::link_files_bulk(&db, &file_participant_map)?;

    if format == "json" {
        let result = serde_json::json!({
            "success": true,
            "updated": updated
        });
        println!("{}", serde_json::to_string_pretty(&result)?);
    } else {
        println!(
            "{}",
            format!("✓ Linked {} files to participants", updated)
                .green()
                .bold()
        );
    }

    Ok(())
}

pub async fn detect(files: Vec<String>, format: String) -> Result<()> {
    use std::collections::HashMap;

    let mut results = HashMap::new();

    for file_path in &files {
        match data::detect_genotype_metadata(file_path) {
            Ok(metadata) => {
                results.insert(file_path.clone(), metadata);
            }
            Err(e) => {
                eprintln!("Warning: Failed to detect {}: {}", file_path, e);
            }
        }
    }

    if format == "json" {
        let response = serde_json::json!({
            "success": true,
            "detections": results
        });
        println!("{}", serde_json::to_string_pretty(&response)?);
    } else {
        println!("{}", "File Type Detection Results:".bold());
        println!();

        for (file_path, metadata) in &results {
            println!("  {}", file_path.cyan());
            println!("    Data Type: {}", metadata.data_type);
            if let Some(source) = &metadata.source {
                println!("    Source: {}", source);
            }
            if let Some(version) = &metadata.grch_version {
                println!("    GRCh Version: {}", version);
            }
            println!();
        }

        let genotype_count = results
            .values()
            .filter(|m| m.data_type == "Genotype")
            .count();
        println!(
            "  Detected {} genotype files out of {} total",
            genotype_count,
            results.len()
        );
    }

    Ok(())
}

pub async fn unlink(file_id: i64, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;

    let updated_file = data::unlink_file(&db, file_id)?;

    if format == "json" {
        let response = CliResponse::new(&updated_file);
        println!("{}", response.to_json()?);
    } else {
        println!("{}", format!("✓ Unlinked file {}", file_id).green().bold());
        println!("  File: {}", updated_file.file_path);
        println!("  Participant: {}", "-".dimmed());
    }

    Ok(())
}

pub async fn analyze(files: Vec<String>, format: String) -> Result<()> {
    use std::collections::HashMap;

    let mut results = HashMap::new();

    for file_path in &files {
        match data::analyze_genotype_file(file_path) {
            Ok(metadata) => {
                results.insert(file_path.clone(), metadata);
            }
            Err(e) => {
                eprintln!("Warning: Failed to analyze {}: {}", file_path, e);
            }
        }
    }

    if format == "json" {
        let response = serde_json::json!({
            "success": true,
            "analysis": results
        });
        println!("{}", serde_json::to_string_pretty(&response)?);
    } else {
        println!("{}", "File Analysis Results:".bold());
        println!();

        for (file_path, metadata) in &results {
            println!("{}", file_path.cyan());
            if let Some(rows) = metadata.row_count {
                println!("  Rows: {}", rows.to_string().green());
            }
            if let Some(chroms) = metadata.chromosome_count {
                println!("  Chromosomes: {}", chroms.to_string().green());
            }
            if let Some(sex) = &metadata.inferred_sex {
                println!("  Inferred Sex: {}", sex.yellow());
            }
            println!();
        }

        println!(
            "Analyzed {} files",
            results.len().to_string().green().bold()
        );
    }

    Ok(())
}

pub async fn hash(files: Vec<String>, format: String) -> Result<()> {
    use crate::cli::download_cache::calculate_blake3;
    use std::collections::HashMap;
    use std::path::Path;

    let mut results = HashMap::new();

    for file_path in &files {
        let path = Path::new(file_path);

        if !path.exists() {
            eprintln!("Warning: File not found: {}", file_path);
            continue;
        }

        match calculate_blake3(path) {
            Ok(hash) => {
                results.insert(file_path.clone(), hash);
            }
            Err(e) => {
                eprintln!("Warning: Failed to hash {}: {}", file_path, e);
            }
        }
    }

    if format == "json" {
        let response = serde_json::json!({
            "success": true,
            "hashes": results
        });
        println!("{}", serde_json::to_string_pretty(&response)?);
    } else {
        println!("{}", "File Hashes (BLAKE3):".bold());
        println!();

        for (file_path, hash) in &results {
            println!("{}", file_path.cyan());
            println!("  {}", hash);
            println!();
        }

        println!("Hashed {} files", results.len().to_string().green().bold());
    }

    Ok(())
}

pub async fn export_csv(
    path: String,
    extension: Option<String>,
    recursive: bool,
    pattern: Option<String>,
    output: String,
) -> Result<()> {
    use std::fs::File;
    use std::io::Write;

    // Scan for files
    let scan_result = data::scan(&path, extension.as_deref(), recursive)?;

    if scan_result.files.is_empty() {
        println!("{}", "No files found to export.".yellow());
        return Ok(());
    }

    println!(
        "{}",
        format!("📊 Found {} files", scan_result.total_files).bold()
    );

    let compiled_pattern = pattern
        .as_ref()
        .map(|pat| data::compile_pattern(pat))
        .transpose()?;

    // Create CSV file
    let mut csv_file = File::create(&output)?;

    // Write header
    writeln!(
        csv_file,
        "file_path,participant_id,data_type,source,grch_version,row_count,chromosome_count,inferred_sex"
    )?;

    // Write rows
    let mut count = 0;
    for file_info in scan_result.files {
        // Extract participant ID if pattern provided
        let participant_id = compiled_pattern
            .as_ref()
            .and_then(|compiled| compiled.extract(&file_info.path))
            .unwrap_or_default();

        writeln!(csv_file, "{},{},,,,,,", file_info.path, participant_id)?;
        count += 1;
    }

    println!();
    println!(
        "{}",
        format!("✓ Exported {} files to {}", count, output)
            .green()
            .bold()
    );

    Ok(())
}

pub async fn detect_csv(input_csv: String, output: String) -> Result<()> {
    use csv::{ReaderBuilder, WriterBuilder};
    use std::fs::File;
    use std::io::Write;

    println!(
        "{}",
        format!("🔍 Detecting file types from {}", input_csv).bold()
    );

    // Read input CSV
    let file = File::open(&input_csv)?;
    let mut rdr = ReaderBuilder::new().from_reader(file);

    let mut rows: Vec<std::collections::HashMap<String, String>> = Vec::new();
    for result in rdr.deserialize() {
        let row: std::collections::HashMap<String, String> = result?;
        rows.push(row);
    }

    if rows.is_empty() {
        println!("{}", "No files found in CSV.".yellow());
        return Ok(());
    }

    let total_files = rows.len();
    println!("{}", format!("📋 Processing {} files", total_files).bold());

    // Detect each file
    for (i, row) in rows.iter_mut().enumerate() {
        if let Some(file_path) = row.get("file_path") {
            print!("\r🔍 Detecting... {}/{}", i + 1, total_files);
            std::io::stdout().flush()?;

            match data::detect_genotype_metadata(file_path) {
                Ok(metadata) => {
                    row.insert("data_type".to_string(), metadata.data_type);
                    // Always insert source and grch_version, even if empty
                    row.insert("source".to_string(), metadata.source.unwrap_or_default());
                    row.insert(
                        "grch_version".to_string(),
                        metadata.grch_version.unwrap_or_default(),
                    );
                }
                Err(e) => {
                    eprintln!("\nWarning: Failed to detect {}: {}", file_path, e);
                    // Set defaults on error
                    row.insert("data_type".to_string(), "Unknown".to_string());
                    row.insert("source".to_string(), String::new());
                    row.insert("grch_version".to_string(), String::new());
                }
            }
        }
    }
    println!();

    // Write output CSV
    let out_file = File::create(&output)?;
    let mut wtr = WriterBuilder::new().from_writer(out_file);

    // Write header
    wtr.write_record([
        "file_path",
        "participant_id",
        "data_type",
        "source",
        "grch_version",
        "row_count",
        "chromosome_count",
        "inferred_sex",
    ])?;

    // Write rows
    for row in &rows {
        wtr.write_record([
            row.get("file_path").unwrap_or(&String::new()),
            row.get("participant_id").unwrap_or(&String::new()),
            row.get("data_type").unwrap_or(&String::new()),
            row.get("source").unwrap_or(&String::new()),
            row.get("grch_version").unwrap_or(&String::new()),
            row.get("row_count").unwrap_or(&String::new()),
            row.get("chromosome_count").unwrap_or(&String::new()),
            row.get("inferred_sex").unwrap_or(&String::new()),
        ])?;
    }

    wtr.flush()?;

    println!(
        "{}",
        format!("✓ Updated CSV written to {}", output)
            .green()
            .bold()
    );

    Ok(())
}

pub async fn analyze_csv(input_csv: String, output: String) -> Result<()> {
    use csv::{ReaderBuilder, WriterBuilder};
    use std::fs::File;
    use std::io::Write;

    println!(
        "{}",
        format!("🧬 Analyzing files from {}", input_csv).bold()
    );

    // Read input CSV
    let file = File::open(&input_csv)?;
    let mut rdr = ReaderBuilder::new().from_reader(file);

    let mut rows: Vec<std::collections::HashMap<String, String>> = Vec::new();
    for result in rdr.deserialize() {
        let row: std::collections::HashMap<String, String> = result?;
        rows.push(row);
    }

    if rows.is_empty() {
        println!("{}", "No files found in CSV.".yellow());
        return Ok(());
    }

    let total_files = rows.len();
    println!(
        "{}",
        format!(
            "📋 Processing {} files (this may take a while...)",
            total_files
        )
        .bold()
    );

    // Analyze each file
    for (i, row) in rows.iter_mut().enumerate() {
        if let Some(file_path) = row.get("file_path") {
            print!("\r🧬 Analyzing... {}/{}", i + 1, total_files);
            std::io::stdout().flush()?;

            match data::analyze_genotype_file(file_path) {
                Ok(metadata) => {
                    // Always insert these fields, even if empty
                    row.insert(
                        "row_count".to_string(),
                        metadata
                            .row_count
                            .map(|c| c.to_string())
                            .unwrap_or_default(),
                    );
                    row.insert(
                        "chromosome_count".to_string(),
                        metadata
                            .chromosome_count
                            .map(|c| c.to_string())
                            .unwrap_or_default(),
                    );
                    row.insert(
                        "inferred_sex".to_string(),
                        metadata.inferred_sex.unwrap_or_default(),
                    );
                }
                Err(e) => {
                    eprintln!("\nWarning: Failed to analyze {}: {}", file_path, e);
                    // Set defaults on error
                    row.insert("row_count".to_string(), String::new());
                    row.insert("chromosome_count".to_string(), String::new());
                    row.insert("inferred_sex".to_string(), String::new());
                }
            }
        }
    }
    println!();

    // Write output CSV
    let out_file = File::create(&output)?;
    let mut wtr = WriterBuilder::new().from_writer(out_file);

    // Write header
    wtr.write_record([
        "file_path",
        "participant_id",
        "data_type",
        "source",
        "grch_version",
        "row_count",
        "chromosome_count",
        "inferred_sex",
    ])?;

    // Write rows
    for row in &rows {
        wtr.write_record([
            row.get("file_path").unwrap_or(&String::new()),
            row.get("participant_id").unwrap_or(&String::new()),
            row.get("data_type").unwrap_or(&String::new()),
            row.get("source").unwrap_or(&String::new()),
            row.get("grch_version").unwrap_or(&String::new()),
            row.get("row_count").unwrap_or(&String::new()),
            row.get("chromosome_count").unwrap_or(&String::new()),
            row.get("inferred_sex").unwrap_or(&String::new()),
        ])?;
    }

    wtr.flush()?;

    println!(
        "{}",
        format!("✓ Updated CSV written to {}", output)
            .green()
            .bold()
    );

    Ok(())
}

pub async fn import_csv(
    csv_path: String,
    non_interactive: bool,
    format: String,
    save_skipped: Option<String>,
) -> Result<()> {
    use csv::ReaderBuilder;
    use serde::Deserialize;
    use std::fs::File;

    #[derive(Debug, Deserialize, Clone)]
    struct CsvRow {
        file_path: String,
        #[serde(default)]
        participant_id: Option<String>,
        #[serde(default)]
        data_type: Option<String>,
        #[serde(default)]
        source: Option<String>,
        #[serde(default)]
        grch_version: Option<String>,
        #[serde(default)]
        row_count: Option<i64>,
        #[serde(default)]
        chromosome_count: Option<i64>,
        #[serde(default)]
        inferred_sex: Option<String>,
    }

    // Read CSV file
    let file = File::open(&csv_path)?;
    let mut rdr = ReaderBuilder::new().has_headers(true).from_reader(file);

    let mut rows: Vec<CsvRow> = Vec::new();
    for result in rdr.deserialize() {
        let row: CsvRow = result?;
        rows.push(row);
    }

    if rows.is_empty() {
        println!("{}", "No files found in CSV.".yellow());
        return Ok(());
    }

    if format != "json" {
        println!("{}", format!("📋 CSV Import Preview: {}", csv_path).bold());
        println!(
            "  Files to import: {}",
            rows.len().to_string().cyan().bold()
        );

        use std::collections::BTreeMap;
        let mut ext_counts: BTreeMap<String, usize> = BTreeMap::new();
        for row in &rows {
            let ext = std::path::Path::new(&row.file_path)
                .extension()
                .and_then(|e| e.to_str())
                .map(|s| s.to_lowercase())
                .unwrap_or_else(|| "(none)".to_string());
            *ext_counts.entry(ext).or_default() += 1;
        }
        println!("  File extensions:");
        for (ext, count) in &ext_counts {
            println!("    {:>6}: {}", ext, count);
        }

        let data_type_filled = rows
            .iter()
            .filter(|r| {
                r.data_type
                    .as_ref()
                    .map(|s| !s.trim().is_empty())
                    .unwrap_or(false)
            })
            .count();
        println!("  Rows with data_type: {}/{}", data_type_filled, rows.len());
        println!();
    }

    let mut skipped_rows: Vec<CsvRow> = Vec::new();

    if !non_interactive && format != "json" {
        use dialoguer::Confirm;

        let is_blank_or_unknown = |value: &Option<String>| {
            value
                .as_ref()
                .map(|s| s.trim().is_empty() || s.eq_ignore_ascii_case("unknown"))
                .unwrap_or(true)
        };

        let mut skipped_total = 0usize;

        // Check if row is missing any of the four critical fields: participant_id, data_type, source, grch_version
        let is_incomplete = |row: &CsvRow| -> bool {
            is_blank_or_unknown(&row.participant_id)
                || is_blank_or_unknown(&row.data_type)
                || is_blank_or_unknown(&row.source)
                || is_blank_or_unknown(&row.grch_version)
        };

        let incomplete_count = rows.iter().filter(|row| is_incomplete(row)).count();

        // Step 1: If there are incomplete rows, ask if user wants to skip them
        if incomplete_count > 0 {
            println!(
                "{}",
                format!(
                    "  ⚠️  {} row(s) are missing metadata (participant_id, data_type, source, or grch_version)",
                    incomplete_count
                )
                .yellow()
            );

            if Confirm::new()
                .with_prompt(format!(
                    "Skip all {} incomplete row(s) and import only complete files?",
                    incomplete_count
                ))
                .default(false)
                .interact()?
            {
                // User wants to skip incomplete files - collect them before filtering
                let to_skip: Vec<CsvRow> = rows
                    .iter()
                    .filter(|row| is_incomplete(row))
                    .cloned()
                    .collect();
                skipped_rows.extend(to_skip);
                rows.retain(|row| !is_incomplete(row));
                skipped_total += incomplete_count;

                if rows.is_empty() {
                    println!("{}", "No rows left to import after filtering.".yellow());
                    return Ok(());
                }
            } else {
                // User doesn't want to skip - offer to run detect-csv
                println!();
                if Confirm::new()
                    .with_prompt("Run metadata detection (detect-csv) to fill in missing fields?")
                    .default(true)
                    .interact()?
                {
                    // Run detect-csv equivalent inline
                    use std::io::Write;
                    println!("{}", "🔍 Detecting metadata...".bold());

                    let total_rows = rows.len();
                    for (i, row) in rows.iter_mut().enumerate() {
                        print!("\r  Processing... {}/{}", i + 1, total_rows);
                        std::io::stdout().flush()?;

                        match data::detect_genotype_metadata(&row.file_path) {
                            Ok(metadata) => {
                                // Only update if currently blank/unknown
                                if is_blank_or_unknown(&row.data_type) {
                                    row.data_type = Some(metadata.data_type);
                                }
                                if is_blank_or_unknown(&row.source) {
                                    row.source = metadata.source;
                                }
                                if is_blank_or_unknown(&row.grch_version) {
                                    row.grch_version = metadata.grch_version;
                                }
                            }
                            Err(e) => {
                                eprintln!("\n  Warning: Failed to detect {}: {}", row.file_path, e);
                            }
                        }
                    }
                    println!();
                    println!("{}", "  ✓ Detection complete".green());
                    println!();
                }

                // Step 2: Check if there are still incomplete rows, and ask category by category
                // Recount after each skip to avoid double-counting rows with multiple missing fields

                let missing_participant_id = rows
                    .iter()
                    .filter(|r| is_blank_or_unknown(&r.participant_id))
                    .count();
                if missing_participant_id > 0 {
                    println!("  {} row(s) missing participant_id", missing_participant_id);
                    if Confirm::new()
                        .with_prompt(format!("Skip these {} row(s)?", missing_participant_id))
                        .default(false)
                        .interact()?
                    {
                        let to_skip: Vec<CsvRow> = rows
                            .iter()
                            .filter(|row| is_blank_or_unknown(&row.participant_id))
                            .cloned()
                            .collect();
                        let before = rows.len();
                        rows.retain(|row| !is_blank_or_unknown(&row.participant_id));
                        let actually_removed = before - rows.len();
                        skipped_rows.extend(to_skip);
                        skipped_total += actually_removed;
                    }
                }

                let missing_data_type = rows
                    .iter()
                    .filter(|r| is_blank_or_unknown(&r.data_type))
                    .count();
                if missing_data_type > 0 {
                    println!("  {} row(s) missing data_type", missing_data_type);
                    if Confirm::new()
                        .with_prompt(format!("Skip these {} row(s)?", missing_data_type))
                        .default(false)
                        .interact()?
                    {
                        let to_skip: Vec<CsvRow> = rows
                            .iter()
                            .filter(|row| is_blank_or_unknown(&row.data_type))
                            .cloned()
                            .collect();
                        let before = rows.len();
                        rows.retain(|row| !is_blank_or_unknown(&row.data_type));
                        let actually_removed = before - rows.len();
                        skipped_rows.extend(to_skip);
                        skipped_total += actually_removed;
                    }
                }

                let missing_source = rows
                    .iter()
                    .filter(|r| is_blank_or_unknown(&r.source))
                    .count();
                if missing_source > 0 {
                    println!("  {} row(s) missing source", missing_source);
                    if Confirm::new()
                        .with_prompt(format!("Skip these {} row(s)?", missing_source))
                        .default(false)
                        .interact()?
                    {
                        let to_skip: Vec<CsvRow> = rows
                            .iter()
                            .filter(|row| is_blank_or_unknown(&row.source))
                            .cloned()
                            .collect();
                        let before = rows.len();
                        rows.retain(|row| !is_blank_or_unknown(&row.source));
                        let actually_removed = before - rows.len();
                        skipped_rows.extend(to_skip);
                        skipped_total += actually_removed;
                    }
                }

                let missing_grch = rows
                    .iter()
                    .filter(|r| is_blank_or_unknown(&r.grch_version))
                    .count();
                if missing_grch > 0 {
                    println!("  {} row(s) missing grch_version", missing_grch);
                    if Confirm::new()
                        .with_prompt(format!("Skip these {} row(s)?", missing_grch))
                        .default(false)
                        .interact()?
                    {
                        let to_skip: Vec<CsvRow> = rows
                            .iter()
                            .filter(|row| is_blank_or_unknown(&row.grch_version))
                            .cloned()
                            .collect();
                        let before = rows.len();
                        rows.retain(|row| !is_blank_or_unknown(&row.grch_version));
                        let actually_removed = before - rows.len();
                        skipped_rows.extend(to_skip);
                        skipped_total += actually_removed;
                    }
                }

                if rows.is_empty() {
                    println!("{}", "No rows left to import after filtering.".yellow());
                    return Ok(());
                }
            }
        }

        if skipped_total > 0 {
            println!(
                "{}",
                format!(
                    "  ⏭️  Skipped {} row(s) due to incomplete metadata",
                    skipped_total
                )
                .yellow()
            );
            println!();
        }
    }

    // Perform the actual import as pending (instant - no hashing or metadata detection)
    let db = BioVaultDb::new()?;
    let result = data::import_files_as_pending(
        &db,
        rows.iter()
            .map(|r| data::CsvFileImport {
                file_path: r.file_path.clone(),
                participant_id: r
                    .participant_id
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                data_type: r
                    .data_type
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                source: r.source.as_ref().filter(|s| !s.trim().is_empty()).cloned(),
                grch_version: r
                    .grch_version
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                row_count: r.row_count,
                chromosome_count: r.chromosome_count,
                inferred_sex: r
                    .inferred_sex
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
            })
            .collect(),
        None, // No collection for CSV imports (collections handled separately)
    )?;

    if format == "json" {
        let response = CliResponse::new(&result);
        println!("{}", response.to_json()?);
    } else {
        println!();
        if result.imported > 0 {
            println!(
                "{}",
                format!(
                    "✓ Added {} files to queue (status: pending)",
                    result.imported
                )
                .green()
                .bold()
            );
        }
        if result.skipped > 0 {
            println!("⊘ Skipped {} files (already imported)", result.skipped);
        }
        if !result.errors.is_empty() {
            println!("{}", format!("✗ {} errors:", result.errors.len()).red());
            for error in &result.errors {
                println!("  {}", error.red());
            }
        }
        println!();
        println!(
            "{}",
            "ℹ️  Files added with status='pending' and queue_added_at timestamp.".dimmed()
        );
        println!(
            "{}",
            "   Run 'bv files process-queue' to hash files and detect metadata.".dimmed()
        );
    }

    // Print and/or save skipped rows
    if !skipped_rows.is_empty() && format != "json" {
        println!();
        println!(
            "{}",
            format!("⏭️  {} row(s) were skipped:", skipped_rows.len())
                .yellow()
                .bold()
        );
        for row in &skipped_rows {
            println!("  • {}", row.file_path);
        }

        // Save to CSV if flag provided
        if let Some(output_path) = save_skipped {
            use csv::WriterBuilder;
            use std::fs::File;

            let out_file = File::create(&output_path)?;
            let mut wtr = WriterBuilder::new().from_writer(out_file);

            // Write header
            wtr.write_record([
                "file_path",
                "participant_id",
                "data_type",
                "source",
                "grch_version",
                "row_count",
                "chromosome_count",
                "inferred_sex",
            ])?;

            // Write skipped rows
            for row in &skipped_rows {
                wtr.write_record([
                    row.file_path.as_str(),
                    row.participant_id.as_deref().unwrap_or(""),
                    row.data_type.as_deref().unwrap_or(""),
                    row.source.as_deref().unwrap_or(""),
                    row.grch_version.as_deref().unwrap_or(""),
                    &row.row_count.map(|c| c.to_string()).unwrap_or_default(),
                    &row.chromosome_count
                        .map(|c| c.to_string())
                        .unwrap_or_default(),
                    row.inferred_sex.as_deref().unwrap_or(""),
                ])?;
            }

            wtr.flush()?;

            println!();
            println!(
                "{}",
                format!("✓ Skipped rows saved to {}", output_path)
                    .green()
                    .bold()
            );
        }
    }

    Ok(())
}

pub async fn import_pending(csv_path: String, format: String) -> Result<()> {
    use csv::ReaderBuilder;
    use serde::Deserialize;
    use std::fs::File;

    #[derive(Debug, Deserialize)]
    struct CsvRow {
        file_path: String,
        #[serde(default)]
        participant_id: Option<String>,
        #[serde(default)]
        data_type: Option<String>,
        #[serde(default)]
        source: Option<String>,
        #[serde(default)]
        grch_version: Option<String>,
    }

    // Read CSV file
    let file = File::open(&csv_path)?;
    let mut rdr = ReaderBuilder::new().has_headers(true).from_reader(file);

    let mut rows: Vec<CsvRow> = Vec::new();
    for result in rdr.deserialize() {
        let row: CsvRow = result?;
        rows.push(row);
    }

    if rows.is_empty() {
        println!("{}", "No files found in CSV.".yellow());
        return Ok(());
    }

    if format != "json" {
        println!("{}", format!("📋 Fast Import: {}", csv_path).bold());
        println!(
            "  Files to add to queue: {}",
            rows.len().to_string().cyan().bold()
        );
        println!();
    }

    // Perform the fast import (no hashing, just add as pending)
    let db = BioVaultDb::new()?;
    let result = data::import_files_as_pending(
        &db,
        rows.iter()
            .map(|r| data::CsvFileImport {
                file_path: r.file_path.clone(),
                participant_id: r
                    .participant_id
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                data_type: r
                    .data_type
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                source: r.source.as_ref().filter(|s| !s.trim().is_empty()).cloned(),
                grch_version: r
                    .grch_version
                    .as_ref()
                    .filter(|s| !s.trim().is_empty())
                    .cloned(),
                row_count: None,
                chromosome_count: None,
                inferred_sex: None,
            })
            .collect(),
        None, // No collection for CSV imports (collections handled separately)
    )?;

    if format == "json" {
        let response = CliResponse::new(&result);
        println!("{}", response.to_json()?);
    } else {
        println!();
        if result.imported > 0 {
            println!(
                "{}",
                format!("✓ Added {} files to queue", result.imported)
                    .green()
                    .bold()
            );
        }
        if result.skipped > 0 {
            println!("⊘ Skipped {} files (already imported)", result.skipped);
        }
        if !result.errors.is_empty() {
            println!("{}", format!("✗ {} errors:", result.errors.len()).red());
            for error in &result.errors {
                println!("  {}", error.red());
            }
        }
        println!();
        println!(
            "{}",
            "💡 Run 'bv files process-queue' to process the queue".yellow()
        );
    }

    Ok(())
}

pub async fn process_queue(limit: usize, daemon: bool, format: String) -> Result<()> {
    use std::time::Duration;
    use tokio::time::sleep;

    let db = BioVaultDb::new()?;

    if daemon {
        if format != "json" {
            println!("{}", "🔄 Starting queue processor daemon...".bold().cyan());
            println!("Press Ctrl+C to stop\n");
        }

        loop {
            let processed = process_batch(&db, limit, &format).await?;

            if processed == 0 && format != "json" {
                // No files to process, wait before checking again
                print!(".");
                std::io::Write::flush(&mut std::io::stdout())?;
                sleep(Duration::from_secs(5)).await;
            }
        }
    } else {
        // One-time processing
        process_batch(&db, limit, &format).await?;
    }

    Ok(())
}

async fn process_batch(db: &BioVaultDb, limit: usize, format: &str) -> Result<usize> {
    // Get pending files from database
    let pending_files = data::get_pending_files(db, limit)?;

    if pending_files.is_empty() {
        if format != "json" {
            println!("{}", "No pending files in queue.".yellow());
        }
        return Ok(0);
    }

    if format != "json" {
        println!(
            "{}",
            format!("📦 Processing {} pending files...", pending_files.len()).bold()
        );
        println!();
    }

    let mut processed = 0;
    let mut errors = 0;

    for (index, file) in pending_files.iter().enumerate() {
        if format != "json" {
            println!(
                "  [{}/{}] {}",
                index + 1,
                pending_files.len(),
                file.file_path
            );
        }

        // Mark as processing
        data::update_file_status(db, file.id, "processing", None)?;

        // Process the file: hash + detect metadata
        match process_single_file(db, file).await {
            Ok(_) => {
                data::update_file_status(db, file.id, "complete", None)?;
                processed += 1;
                if format != "json" {
                    println!("    {}", "✓ Complete".green());
                }
            }
            Err(e) => {
                let error_msg = format!("{}", e);
                data::update_file_status(db, file.id, "error", Some(&error_msg))?;
                errors += 1;
                if format != "json" {
                    println!("    {}", format!("✗ Error: {}", error_msg).red());
                }
            }
        }
    }

    if format == "json" {
        use serde_json::json;
        let result = json!({
            "processed": processed,
            "errors": errors,
            "total": pending_files.len()
        });
        println!("{}", serde_json::to_string_pretty(&result)?);
    } else {
        println!();
        println!(
            "{}",
            format!("✓ Processed {} files", processed).green().bold()
        );
        if errors > 0 {
            println!("{}", format!("✗ {} errors", errors).red());
        }
    }

    Ok(processed)
}

async fn process_single_file(db: &BioVaultDb, file: &data::PendingFile) -> Result<()> {
    // 1. Hash the file
    let hash = data::hash_file(&file.file_path)?;

    // 2. Detect genotype metadata if not already set
    let mut metadata = if file.data_type.as_deref() == Some("Unknown") || file.data_type.is_none() {
        data::detect_genotype_metadata(&file.file_path).ok()
    } else if file.data_type.as_deref() == Some("Genotype") {
        // Already detected as Genotype, load existing metadata if available
        match data::get_genotype_metadata(db, file.id) {
            Ok(Some(existing)) => Some(existing),
            _ => {
                // No existing metadata, create placeholder
                Some(data::GenotypeMetadata {
                    data_type: "Genotype".to_string(),
                    source: None,
                    grch_version: None,
                    row_count: None,
                    chromosome_count: None,
                    inferred_sex: None,
                })
            }
        }
    } else {
        None
    };

    // 3. If this is a Genotype file, analyze it for row counts, chromosomes, sex
    if let Some(ref mut meta) = metadata {
        if meta.data_type == "Genotype" {
            match data::analyze_genotype_file(&file.file_path) {
                Ok(analysis) => {
                    // Merge analysis data into metadata (preserve existing source/grch if not in analysis)
                    if meta.row_count.is_none() {
                        meta.row_count = analysis.row_count;
                    }
                    if meta.chromosome_count.is_none() {
                        meta.chromosome_count = analysis.chromosome_count;
                    }
                    if meta.inferred_sex.is_none() {
                        meta.inferred_sex = analysis.inferred_sex.clone();
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to analyze {}: {}", file.file_path, e);
                    // Continue with basic metadata
                }
            }
        }
    }

    // 4. Update the file in database
    data::update_file_from_queue(db, file.id, &hash, metadata.as_ref())?;

    Ok(())
}
