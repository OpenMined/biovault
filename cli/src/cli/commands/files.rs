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
        if let Some(pat) = &pattern {
            println!();
            println!("{}", "Sample participant ID extractions:".bold());
            for file_info in scan_result.files.iter().take(5) {
                let filename = std::path::Path::new(&file_info.path)
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(&file_info.path);

                if let Some(id) = data::extract_id_from_pattern(filename, pat) {
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

    if format == "json" {
        let response = CliResponse::new(&result);
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
