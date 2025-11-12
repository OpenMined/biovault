use anyhow::Result;
use colored::Colorize;

use crate::data::{self, generate_variable_name, BioVaultDb, CliResponse};

pub async fn create(
    name: String,
    description: Option<String>,
    var_name: Option<String>,
    format: String,
) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Show suggested variable name if not provided
    let suggested_var_name = generate_variable_name(&name);

    if var_name.is_none() && format != "json" {
        println!(
            "{}",
            format!("💡 Suggested variable name: {}", suggested_var_name).dimmed()
        );
    }

    let collection = data::create_collection(&db, name.clone(), description, var_name)?;

    if format == "json" {
        let response = CliResponse::new(&collection);
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!("✓ Created collection '{}'", collection.name)
                .green()
                .bold()
        );
        println!("  Variable name: {}", collection.variable_name.cyan());
        if let Some(desc) = &collection.description {
            println!("  Description: {}", desc);
        }
        println!("  ID: {}", collection.id);
    }

    Ok(())
}

pub async fn list(format: String) -> Result<()> {
    let db = BioVaultDb::new()?;
    let collections = data::list_collections(&db)?;

    if format == "json" {
        let response = CliResponse::new(&collections);
        println!("{}", response.to_json()?);
    } else {
        if collections.is_empty() {
            println!("{}", "No collections found.".yellow());
            return Ok(());
        }

        println!("{}", "Collections:".bold());
        println!();
        println!(
            "  {}  {}  {}  {}",
            "ID".bold(),
            "Name".bold(),
            "Variable Name".bold(),
            "Files".bold()
        );

        for collection in &collections {
            println!(
                "  {}  {}  {}  {}",
                collection.id,
                collection.name,
                collection.variable_name.cyan(),
                collection.file_count
            );
        }
    }

    Ok(())
}

pub async fn show(identifier: String, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;
    let collection = data::get_collection_detail(&db, &identifier)?;

    if format == "json" {
        let response = CliResponse::new(&collection);
        println!("{}", response.to_json()?);
    } else {
        println!("{}", format!("Collection: {}", collection.name).bold());
        println!("  ID: {}", collection.id);
        println!("  Variable name: {}", collection.variable_name.cyan());
        if let Some(desc) = &collection.description {
            println!("  Description: {}", desc);
        }
        println!("  Files: {}", collection.files.len());
        println!("  Created: {}", collection.created_at);
        println!("  Updated: {}", collection.updated_at);

        if !collection.files.is_empty() {
            println!();
            println!("{}", "Files:".bold());
            println!(
                "  {}  {}  {}  {}  {}",
                "ID".bold(),
                "Path".bold(),
                "Type".bold(),
                "Size".bold(),
                "Participant".bold()
            );

            for file in &collection.files {
                let size_str = file
                    .file_size
                    .map(|s| {
                        if s < 1024 {
                            format!("{} B", s)
                        } else if s < 1_048_576 {
                            format!("{:.1} KB", s as f64 / 1024.0)
                        } else {
                            format!("{:.1} MB", s as f64 / 1_048_576.0)
                        }
                    })
                    .unwrap_or_else(|| "-".to_string());

                let file_type = file.file_type.as_deref().unwrap_or("-");
                let participant = file.participant_id.as_deref().unwrap_or("-");

                println!(
                    "  {}  {}  {}  {}  {}",
                    file.id, file.file_path, file_type, size_str, participant
                );
            }
        }
    }

    Ok(())
}

pub async fn add_files(collection: String, file_ids: Vec<i64>, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;
    let added = data::add_files_to_collection(&db, &collection, file_ids.clone())?;

    if format == "json" {
        let json_data = serde_json::json!({
            "collection": collection,
            "file_ids": file_ids,
            "added": added
        });
        let response = CliResponse::new(&json_data);
        println!("{}", response.to_json()?);
    } else if added > 0 {
        println!(
            "{}",
            format!("✓ Added {} file(s) to collection '{}'", added, collection)
                .green()
                .bold()
        );
    } else {
        println!(
            "{}",
            "⊘ No files added (all files were already in collection)"
                .to_string()
                .yellow()
        );
    }

    Ok(())
}

pub async fn remove_files(collection: String, file_ids: Vec<i64>, format: String) -> Result<()> {
    // Prevent removing files from "Unsorted Files" from CLI - it's a UI-only virtual collection
    if collection == "unsorted_files" || collection == "Unsorted Files" {
        anyhow::bail!("Cannot remove files from 'Unsorted Files' - it is a virtual collection. Files appear there when not assigned to any collection.");
    }

    let db = BioVaultDb::new()?;
    let removed = data::remove_files_from_collection(&db, &collection, file_ids.clone())?;

    if format == "json" {
        let json_data = serde_json::json!({
            "collection": collection,
            "file_ids": file_ids,
            "removed": removed
        });
        let response = CliResponse::new(&json_data);
        println!("{}", response.to_json()?);
    } else if removed > 0 {
        println!(
            "{}",
            format!(
                "✓ Removed {} file(s) from collection '{}'",
                removed, collection
            )
            .green()
            .bold()
        );
    } else {
        println!("{}", "⊘ No files removed".to_string().yellow());
    }

    Ok(())
}

pub async fn delete(identifier: String, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;
    data::delete_collection(&db, &identifier)?;

    if format == "json" {
        let json_data = serde_json::json!({
            "deleted": identifier
        });
        let response = CliResponse::new(&json_data);
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!("✓ Deleted collection '{}'", identifier)
                .green()
                .bold()
        );
    }

    Ok(())
}

pub async fn update(
    identifier: String,
    name: Option<String>,
    description: Option<String>,
    var_name: Option<String>,
    format: String,
) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Handle empty string as None for description (to clear it)
    let description_opt = description.map(|d| if d.is_empty() { None } else { Some(d) });

    let collection = data::update_collection(&db, &identifier, name, description_opt, var_name)?;

    if format == "json" {
        let response = CliResponse::new(&collection);
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!("✓ Updated collection '{}'", collection.name)
                .green()
                .bold()
        );
        println!("  Variable name: {}", collection.variable_name.cyan());
        if let Some(desc) = &collection.description {
            println!("  Description: {}", desc);
        }
    }

    Ok(())
}
