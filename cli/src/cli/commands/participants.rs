use anyhow::Result;
use colored::Colorize;

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
            "  {} files will be unlinked (participant_id set to NULL)",
            participant.file_count
        );
    }

    data::delete_participant(&db, id)?;

    if format == "json" {
        #[derive(serde::Serialize)]
        struct DeleteResponse {
            deleted_id: i64,
            participant_id: String,
            files_unlinked: i64,
        }

        let response = CliResponse::new(DeleteResponse {
            deleted_id: id,
            participant_id: participant.participant_id.clone(),
            files_unlinked: participant.file_count,
        });
        println!("{}", response.to_json()?);
    } else {
        println!(
            "{}",
            format!(
                "✓ Deleted participant {} (unlinked {} files)",
                participant.participant_id, participant.file_count
            )
            .green()
            .bold()
        );
    }

    Ok(())
}
