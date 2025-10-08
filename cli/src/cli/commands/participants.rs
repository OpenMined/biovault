use anyhow::Result;
use colored::Colorize;
use std::collections::HashMap;

use crate::data::{self, BioVaultDb, CliResponse};

pub async fn list(format: String) -> Result<()> {
    let db = BioVaultDb::new()?;
    let participants = data::list_participants(&db)?;

    if format == "json" {
        let response = CliResponse::new(&participants);
        println!("{}", response.to_json()?);
    } else {
        // Table format
        if participants.is_empty() {
            println!("{}", "No participants found.".yellow());
            println!("Files will create participants automatically during import with --pattern.");
            return Ok(());
        }

        println!(
            "{}",
            format!("Participants ({})", participants.len()).bold()
        );
        println!();
        println!(
            "  {}  {}  {}  {}",
            "ID".bold(),
            "Participant ID".bold(),
            "Files".bold(),
            "Created At".bold()
        );

        for p in &participants {
            println!(
                "  {}  {}  {}  {}",
                p.id,
                p.participant_id.cyan(),
                p.file_count.to_string().green(),
                p.created_at.dimmed()
            );
        }
    }

    Ok(())
}

pub async fn delete(id: i64, format: String) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get participant info before deletion for display
    let participants = data::list_participants(&db)?;
    let participant = participants
        .iter()
        .find(|p| p.id == id)
        .ok_or_else(|| anyhow::anyhow!("Participant with id {} not found", id))?;

    if format != "json" {
        println!(
            "Deleting participant: {} (id: {})",
            participant.participant_id.cyan(),
            id
        );
        println!(
            "  {} files will be deleted from the catalog",
            participant.file_count
        );
    }

    let files_deleted = data::delete_participant(&db, id)?;

    if format == "json" {
        #[derive(serde::Serialize)]
        struct DeleteResponse {
            deleted_id: i64,
            participant_id: String,
            files_deleted: usize,
        }

        let response = CliResponse::new(DeleteResponse {
            deleted_id: id,
            participant_id: participant.participant_id.clone(),
            files_deleted,
        });
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!(
                "✓ Deleted participant {} (deleted {} files)",
                participant.participant_id, files_deleted
            )
            .green()
            .bold()
        );
    }

    Ok(())
}

pub async fn delete_bulk(ids: Vec<i64>, format: String) -> Result<()> {
    #[derive(serde::Serialize)]
    struct DeletedParticipant {
        deleted_id: i64,
        participant_id: String,
        files_deleted: usize,
    }

    if ids.is_empty() {
        if format == "json" {
            #[derive(serde::Serialize)]
            struct EmptyResponse {
                deleted: usize,
                files_deleted: usize,
                results: Vec<DeletedParticipant>,
                not_found: Vec<i64>,
            }

            let response = CliResponse::new(EmptyResponse {
                deleted: 0,
                files_deleted: 0,
                results: Vec::new(),
                not_found: Vec::new(),
            });
            println!("{}", response.to_json()?);
        }
        return Ok(());
    }

    let db = BioVaultDb::new()?;
    let participants = data::list_participants(&db)?;
    let mut participant_map: HashMap<i64, _> =
        participants.into_iter().map(|p| (p.id, p)).collect();

    let mut valid = Vec::new();
    let mut not_found = Vec::new();

    for id in ids {
        if let Some(p) = participant_map.remove(&id) {
            valid.push(p);
        } else {
            not_found.push(id);
        }
    }

    if valid.is_empty() {
        if format == "json" {
            #[derive(serde::Serialize)]
            struct Response {
                deleted: usize,
                files_deleted: usize,
                results: Vec<DeletedParticipant>,
                not_found: Vec<i64>,
            }

            let response = CliResponse::new(Response {
                deleted: 0,
                files_deleted: 0,
                results: Vec::new(),
                not_found,
            });
            println!("{}", response.to_json()?);
        } else {
            println!(
                "{}",
                "No matching participants were found for the specified IDs.".yellow()
            );
            if !not_found.is_empty() {
                println!("Missing IDs: {}", format_ids(&not_found));
            }
        }
        return Ok(());
    }

    let valid_ids: Vec<i64> = valid.iter().map(|p| p.id).collect();
    let deleted = data::delete_participants_bulk(&db, &valid_ids)?;

    let total_files_deleted: usize = valid.iter().map(|p| p.file_count as usize).sum();

    if format == "json" {
        #[derive(serde::Serialize)]
        struct Response {
            deleted: usize,
            files_deleted: usize,
            results: Vec<DeletedParticipant>,
            not_found: Vec<i64>,
        }

        let results = valid
            .iter()
            .map(|p| DeletedParticipant {
                deleted_id: p.id,
                participant_id: p.participant_id.clone(),
                files_deleted: p.file_count as usize,
            })
            .collect();

        let response = CliResponse::new(Response {
            deleted,
            files_deleted: total_files_deleted,
            results,
            not_found,
        });
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!("Deleting {} participant(s)...", valid.len()).bold()
        );
        for p in &valid {
            println!(
                "{}",
                format!(
                    "  ✓ Deleted participant {} (deleted {} files)",
                    p.participant_id.cyan(),
                    p.file_count
                )
                .green()
            );
        }

        if !not_found.is_empty() {
            println!(
                "{}",
                format!(
                    "⚠️  The following participant IDs were not found: {}",
                    format_ids(&not_found)
                )
                .yellow()
            );
        }

        println!(
            "{}",
            format!(
                "✓ Deleted {} participant(s), removed {} file(s)",
                deleted, total_files_deleted
            )
            .green()
            .bold()
        );
    }

    Ok(())
}

fn format_ids(ids: &[i64]) -> String {
    ids.iter()
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}
