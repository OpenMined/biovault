use crate::error::Result;
use crate::types::InboxSubmission;
use colored::Colorize;
use dialoguer::{theme::ColorfulTheme, Select};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

pub async fn list(show_all: bool, show_full: bool) -> Result<()> {
    // Get inbox directory
    let inbox_path = get_inbox_path()?;

    if !inbox_path.exists() {
        println!("📭 No submissions in inbox");
        return Ok(());
    }

    // Collect all inbox submissions
    let mut submissions = Vec::new();

    // Walk through inbox directory structure: ~/.biovault/inbox/{sender_email}/{project_name}-yyyy-mm-dd-hash.yaml
    for entry in WalkDir::new(&inbox_path)
        .min_depth(2)
        .max_depth(2)
        .follow_links(false)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("yaml"))
    {
        let path = entry.path();

        // Extract sender email from path
        let sender = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        // Try to load submission
        match InboxSubmission::from_file(&path.to_path_buf()) {
            Ok(submission) => {
                // Filter out rejected unless --all is specified
                if !show_all && submission.status == "rejected" {
                    continue;
                }

                // Extract date from filename (format: name-yyyy-mm-dd-hash.yaml)
                let filename = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");

                let date_str = extract_date_from_filename(filename);

                submissions.push((sender.to_string(), submission, path.to_path_buf(), date_str));
            }
            Err(e) => {
                eprintln!("Warning: Failed to load submission from {:?}: {}", path, e);
            }
        }
    }

    if submissions.is_empty() {
        if !show_all {
            println!("📭 No active submissions in inbox (use --all to show rejected)");
        } else {
            println!("📭 No submissions in inbox");
        }
        return Ok(());
    }

    // Sort by date (newest first)
    submissions.sort_by(|a, b| b.3.cmp(&a.3));

    // Display submissions
    let filtered_text = if !show_all { " active" } else { "" };
    println!(
        "📬 {}{} submission(s) in inbox:\n",
        submissions.len(),
        filtered_text
    );

    if show_full {
        // Full detailed view
        for (i, (sender, submission, path, _date)) in submissions.iter().enumerate() {
            display_full_submission(i + 1, sender, submission, path);
        }
    } else {
        // Concise table view
        display_concise_list(&submissions);
    }

    Ok(())
}

fn display_concise_list(submissions: &[(String, InboxSubmission, PathBuf, String)]) {
    // Print header
    println!(
        "{:<4} {:<12} {:<20} {:<25} {:<15} {:<30}",
        "#".bold(),
        "Date".bold(),
        "Project".bold(),
        "From".bold(),
        "Participants".bold(),
        "ID".bold()
    );
    println!("{}", "-".repeat(110));

    for (i, (sender, submission, path, date)) in submissions.iter().enumerate() {
        let status_icon = match submission.status.as_str() {
            "pending" => "⏳",
            "approved" => "✅",
            "rejected" => "❌",
            _ => "❓",
        };

        let participants_count = submission
            .participants
            .as_ref()
            .map(|p| p.len().to_string())
            .unwrap_or_else(|| "0".to_string());

        let id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        // Truncate long names for better formatting
        let project_name = if submission.name.len() > 18 {
            format!("{}...", &submission.name[..15])
        } else {
            submission.name.clone()
        };

        let sender_display = if sender.len() > 23 {
            format!("{}...", &sender[..20])
        } else {
            sender.clone()
        };

        let id_display = if id.len() > 28 {
            format!("{}...", &id[..25])
        } else {
            id.to_string()
        };

        println!(
            "{:<4} {:<12} {} {:<18} {:<25} {:<15} {:<30}",
            format!("{}", i + 1),
            date,
            status_icon,
            project_name,
            sender_display.cyan(),
            participants_count,
            id_display.dimmed()
        );
    }
}

fn display_full_submission(index: usize, sender: &str, submission: &InboxSubmission, path: &Path) {
    let status_icon = match submission.status.as_str() {
        "pending" => "⏳",
        "approved" => "✅",
        "rejected" => "❌",
        _ => "❓",
    };

    println!("{}. {} {}", index, status_icon, submission.name.bold());
    println!("   From: {}", sender.cyan());
    println!("   Author: {}", submission.author);
    println!("   Status: {}", format_status(&submission.status));

    if let Some(datasites) = &submission.datasites {
        if !datasites.is_empty() {
            println!("   Datasites: {}", datasites.join(", "));
        }
    }

    if let Some(participants) = &submission.participants {
        if !participants.is_empty() {
            println!("   Participants: {} participant(s)", participants.len());
            if participants.len() <= 3 {
                for participant in participants {
                    println!("     - {}", participant.dimmed());
                }
            }
        }
    }

    println!("   Syft URL: {}", submission.syft_url.dimmed());

    // Extract filename for display
    if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
        println!("   ID: {}", filename.dimmed());
    }

    println!();
}

fn extract_date_from_filename(filename: &str) -> String {
    // Expected format: name-yyyy-mm-dd-hash
    // The hash is 8 characters, so we look for date pattern before it
    if filename.len() < 19 {
        // Minimum: x-yyyy-mm-dd-xxxxxxxx
        return "unknown".to_string();
    }

    // Split by dash and look from the end
    let parts: Vec<&str> = filename.split('-').collect();
    if parts.len() < 4 {
        return "unknown".to_string();
    }

    // The last part is the hash, the 3 before it should be the date
    let len = parts.len();

    // Validate that we have year-month-day pattern
    if let (Ok(year), Ok(month), Ok(day)) = (
        parts[len - 4].parse::<u32>(),
        parts[len - 3].parse::<u32>(),
        parts[len - 2].parse::<u32>(),
    ) {
        // Basic validation
        if (2020..=2100).contains(&year) && (1..=12).contains(&month) && (1..=31).contains(&day) {
            return format!("{:04}-{:02}-{:02}", year, month, day);
        }
    }

    "unknown".to_string()
}

fn get_inbox_path() -> Result<PathBuf> {
    let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
        PathBuf::from(test_home)
    } else {
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
    };
    Ok(home_dir.join(".biovault").join("inbox"))
}

fn format_status(status: &str) -> String {
    match status {
        "pending" => status.yellow().to_string(),
        "approved" => status.green().to_string(),
        "rejected" => status.red().to_string(),
        _ => status.to_string(),
    }
}

pub async fn show(reference: &str, show_all: bool) -> Result<()> {
    let submissions = load_submissions(show_all)?;

    if submissions.is_empty() {
        println!("📭 No submissions found");
        return Ok(());
    }

    // Try to parse as number first
    if let Ok(index) = reference.parse::<usize>() {
        if index > 0 && index <= submissions.len() {
            show_submission_detail(&submissions[index - 1]);
            return Ok(());
        } else {
            println!(
                "❌ Invalid index: {}. Valid range: 1-{}",
                index,
                submissions.len()
            );
            return Ok(());
        }
    }

    // Try to match by partial hash or name
    let reference_lower = reference.to_lowercase();
    let matches: Vec<_> = submissions
        .iter()
        .filter(|(_, submission, path, _)| {
            // Check if it matches the project name
            if submission.name.to_lowercase().contains(&reference_lower) {
                return true;
            }
            // Check if it matches partial ID/hash
            if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                if filename.to_lowercase().contains(&reference_lower) {
                    return true;
                }
            }
            false
        })
        .collect();

    match matches.len() {
        0 => {
            println!("❌ No submission found matching: {}", reference);
            println!("💡 Try using: index number (1,2,3...), partial hash, or project name");
        }
        1 => {
            show_submission_detail(matches[0]);
        }
        _ => {
            println!("⚠️  Multiple submissions match '{}':", reference);
            for (i, (sender, submission, path, _)) in matches.iter().enumerate() {
                let id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                println!(
                    "  {}. {} from {} ({})",
                    i + 1,
                    submission.name.bold(),
                    sender.cyan(),
                    &id[..id.len().min(8)].dimmed()
                );
            }
            println!("\n💡 Use a more specific reference or the full hash");
        }
    }

    Ok(())
}

pub async fn interactive(show_all: bool) -> Result<()> {
    let submissions = load_submissions(show_all)?;

    if submissions.is_empty() {
        println!("📭 No submissions found");
        return Ok(());
    }

    loop {
        // Create display items for selection
        let items: Vec<String> = submissions
            .iter()
            .map(|(sender, submission, path, date)| {
                let status_icon = match submission.status.as_str() {
                    "pending" => "⏳",
                    "approved" => "✅",
                    "rejected" => "❌",
                    _ => "❓",
                };

                let id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                let short_id = if id.len() > 8 {
                    &id[id.len() - 8..]
                } else {
                    id
                };

                format!(
                    "{} {} - {} from {} [{}]",
                    status_icon, submission.name, date, sender, short_id
                )
            })
            .collect();

        // Add exit option
        let mut menu_items = items.clone();
        menu_items.push("📤 Exit".to_string());

        let selection = Select::with_theme(&ColorfulTheme::default())
            .with_prompt("📬 Select submission to view details (↑↓ to navigate, Enter to select, Esc to exit)")
            .items(&menu_items)
            .default(0)
            .interact_opt()
            .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

        match selection {
            Some(index) if index < submissions.len() => {
                println!();
                show_submission_detail(&submissions[index]);

                println!("\n{}", "─".repeat(80));
                println!("Press Enter to continue...");
                let mut _input = String::new();
                std::io::stdin()
                    .read_line(&mut _input)
                    .map_err(|e| anyhow::anyhow!("Failed to read input: {}", e))?;
            }
            _ => {
                // Exit selected or cancelled
                println!("👋 Exiting inbox");
                break;
            }
        }
    }

    Ok(())
}

fn load_submissions(show_all: bool) -> Result<Vec<(String, InboxSubmission, PathBuf, String)>> {
    let inbox_path = get_inbox_path()?;

    if !inbox_path.exists() {
        return Ok(Vec::new());
    }

    let mut submissions = Vec::new();

    for entry in WalkDir::new(&inbox_path)
        .min_depth(2)
        .max_depth(2)
        .follow_links(false)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("yaml"))
    {
        let path = entry.path();

        let sender = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        match InboxSubmission::from_file(&path.to_path_buf()) {
            Ok(submission) => {
                if !show_all && submission.status == "rejected" {
                    continue;
                }

                let filename = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");

                let date_str = extract_date_from_filename(filename);

                submissions.push((sender.to_string(), submission, path.to_path_buf(), date_str));
            }
            Err(e) => {
                eprintln!("Warning: Failed to load submission from {:?}: {}", path, e);
            }
        }
    }

    // Sort by date (newest first)
    submissions.sort_by(|a, b| b.3.cmp(&a.3));

    Ok(submissions)
}

fn show_submission_detail(submission_data: &(String, InboxSubmission, PathBuf, String)) {
    let (sender, submission, path, date) = submission_data;

    let status_icon = match submission.status.as_str() {
        "pending" => "⏳",
        "approved" => "✅",
        "rejected" => "❌",
        _ => "❓",
    };

    println!("{}", "═".repeat(80));
    println!("{} {}", status_icon, submission.name.bold().cyan());
    println!("{}", "─".repeat(80));

    println!("📅 Date:        {}", date.bold());
    println!("📧 From:        {}", sender.cyan());
    println!("✍️  Author:      {}", submission.author);
    println!("📊 Status:      {}", format_status(&submission.status));

    if let Some(datasites) = &submission.datasites {
        if !datasites.is_empty() {
            println!("🏢 Datasites:   {}", datasites.join(", "));
        }
    }

    if let Some(participants) = &submission.participants {
        if !participants.is_empty() {
            println!(
                "👥 Participants: {} participant(s)",
                participants.len().to_string().bold()
            );
            for (i, participant) in participants.iter().enumerate() {
                if i < 5 {
                    println!("     {}. {}", i + 1, participant.dimmed());
                } else if i == 5 {
                    println!("     ... and {} more", participants.len() - 5);
                    break;
                }
            }
        }
    } else {
        println!("👥 Participants: None specified");
    }

    println!("🔗 Syft URL:    {}", submission.syft_url.dimmed());

    if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
        println!("🆔 ID:          {}", filename.dimmed());

        // Show short ID for easy reference
        if filename.len() > 8 {
            let short_id = &filename[filename.len() - 8..];
            println!("🏷️  Short ID:    {}", short_id.yellow());
        }
    }

    println!("{}", "═".repeat(80));
}
