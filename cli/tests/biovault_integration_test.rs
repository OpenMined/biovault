#![cfg(feature = "e2e-tests")]

use anyhow::{anyhow, Context, Result};
use serde::{de::DeserializeOwned, Deserialize};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::thread;
use std::time::Duration;
use walkdir::WalkDir;

const COUNT_LINES_PROJECT: &str = "count-lines";
const COUNT_LINES_PROJECT_NAME: &str = "pipeline-count-lines";
const GENOTYPE_FIXTURES_DIR: &str = "tests/data/genotype_files";
const ALL_SAMPLESHEET_NAME: &str = "all_participants.csv";
const COUNT_LINES_PROJECT_YAML: &str =
    include_str!("../examples/pipeline/count-lines/project.yaml");
const COUNT_LINES_WORKFLOW: &str = include_str!("../examples/pipeline/count-lines/workflow.nf");
const COUNT_LINES_SCRIPT: &str =
    include_str!("../examples/pipeline/count-lines/assets/count_lines.py");
const MESSAGE_WAIT_ATTEMPTS: usize = 10;
const MESSAGE_WAIT_DELAY_SECS: u64 = 2;

fn wait_threshold() -> Option<usize> {
    match env::var("TEST_WAIT_FOR_STEPS") {
        Ok(value) => {
            let normalized = value.trim().to_lowercase();
            match normalized.as_str() {
                "" | "false" | "0" => None,
                "true" | "on" => Some(1),
                _ => normalized.parse::<usize>().ok().filter(|n| *n >= 1),
            }
        }
        Err(_) => None,
    }
}

fn should_pause_between_steps(step_number: usize) -> bool {
    wait_threshold()
        .map(|threshold| step_number >= threshold)
        .unwrap_or(false)
}

fn pause_between_steps(step_number: usize, step_label: &str) {
    if !should_pause_between_steps(step_number) {
        return;
    }

    println!(
        "\n🔎 {} complete. Inspect state as needed, then press Enter to continue...",
        step_label
    );
    print!("Continue> ");
    let _ = io::stdout().flush();
    let mut buffer = String::new();
    let _ = io::stdin().read_line(&mut buffer);
}

fn strip_ansi_codes(input: &str) -> String {
    let mut result = String::with_capacity(input.len());
    let mut chars = input.chars();

    while let Some(c) = chars.next() {
        if c == '\u{1b}' {
            match chars.next() {
                Some('[') => {
                    for next in chars.by_ref() {
                        if ('@'..='~').contains(&next) {
                            break;
                        }
                    }
                }
                Some(']') => {
                    for next in chars.by_ref() {
                        if next == '\u{07}' {
                            break;
                        }
                    }
                }
                _ => {}
            }
        } else {
            result.push(c);
        }
    }

    result
}

fn parse_json_output<T: DeserializeOwned>(stdout: &str, context: &str) -> Result<T> {
    let cleaned = strip_ansi_codes(stdout);
    let trimmed = cleaned.trim_start();
    let start = trimmed
        .find(['{', '['])
        .ok_or_else(|| anyhow!("{}: no JSON payload found", context))?;
    let json_slice = &trimmed[start..];
    let value: serde_json::Value = serde_json::from_str(json_slice)
        .map_err(|e| anyhow!("{}: failed to parse JSON: {}", context, e))?;

    match serde_json::from_value::<T>(value.clone()) {
        Ok(parsed) => Ok(parsed),
        Err(primary_err) => {
            if let Some(data_value) = value.get("data") {
                serde_json::from_value::<T>(data_value.clone()).map_err(|data_err| {
                    anyhow!(
                        "{}: failed to parse response data: {}; original error: {}",
                        context,
                        data_err,
                        primary_err
                    )
                })
            } else {
                Err(anyhow!(
                    "{}: failed to parse JSON: {}",
                    context,
                    primary_err
                ))
            }
        }
    }
}

struct ProcessOutcome {
    run_results: PathBuf,
    submission_folder: String,
}

#[derive(Debug, Deserialize)]
struct FilesListData {
    pub total: usize,
    pub files: Vec<FilesListEntry>,
}

#[derive(Debug, Deserialize)]
struct FilesListEntry {
    pub file_path: String,
    #[serde(default)]
    pub participant_id: Option<String>,
    #[serde(default)]
    pub participant_name: Option<String>,
}

/// End-to-end integration test for BioVault
#[test]
#[ignore]
fn test_biovault_e2e() -> Result<()> {
    println!("\n🧬 Testing BioVault end-to-end workflow!");

    // Get configuration from environment
    let client1_email =
        env::var("SYFTBOX_CLIENT1_EMAIL").unwrap_or_else(|_| "client1@syftbox.net".to_string());
    let client2_email =
        env::var("SYFTBOX_CLIENT2_EMAIL").unwrap_or_else(|_| "client2@syftbox.net".to_string());
    let test_clients_dir =
        env::var("TEST_CLIENTS_DIR").unwrap_or_else(|_| "./test-clients".to_string());
    let test_mode = env::var("TEST_MODE").unwrap_or_else(|_| "docker".to_string());

    println!("  Mode: {}", test_mode);
    println!("  Client1: {}", client1_email);
    println!("  Client2: {}", client2_email);
    println!("  Test dir: {}", test_clients_dir);

    // First, verify we can find the bv binary
    println!("\n📦 Verifying BioVault binary...");
    match get_bv_binary_path() {
        Ok(path) => println!("  Found binary at: {}", path.display()),
        Err(e) => {
            eprintln!("  ERROR: {}", e);
            return Err(e);
        }
    }

    // Build paths for client directories
    let (client1_base, client2_base) = if test_mode == "local" {
        (
            PathBuf::from(&test_clients_dir).join(&client1_email),
            PathBuf::from(&test_clients_dir).join(&client2_email),
        )
    } else {
        (
            PathBuf::from(&test_clients_dir)
                .join(&client1_email)
                .join("SyftBox"),
            PathBuf::from(&test_clients_dir)
                .join(&client2_email)
                .join("SyftBox"),
        )
    };

    // Test 1: Initialize BioVault for both clients
    println!("\n📝 Test 1: Initialize BioVault for both clients");
    test_biovault_init(&client1_base, &client1_email, &test_mode)?;
    test_biovault_init(&client2_base, &client2_email, &test_mode)?;
    pause_between_steps(1, "Test 1: Initialize BioVault for both clients");

    // Test 2: Stage genotype fixtures for processor
    println!("\n📝 Test 2: Stage genotype fixtures for processor");
    let client2_genotype_dir = stage_genotype_fixtures(&client2_base)?;
    let expected_file_count = count_files_in_dir(&client2_genotype_dir)?;
    println!(
        "  Prepared {} genotype files for testing",
        expected_file_count
    );
    pause_between_steps(2, "Test 2: Stage genotype fixtures for processor");

    // Test 3: Import genotype data for processor
    println!("\n📝 Test 3: Import genotype data for processor");
    let imported_count =
        test_import_genotype_data(&client2_base, &client2_genotype_dir, &test_mode)?;
    assert_eq!(
        imported_count, expected_file_count,
        "Imported genotype catalog should match staged fixtures"
    );
    pause_between_steps(3, "Test 3: Import genotype data for processor");

    // Test 4: Prepare dynamic project template for sender
    println!("\n📝 Test 4: Prepare dynamic project template");
    test_prepare_dynamic_project(&client1_base, &test_mode)?;
    pause_between_steps(4, "Test 4: Prepare dynamic project template");

    // Test 5: Submit project from client1 to client2
    println!("\n📝 Test 5: Submit project from client1 to client2");
    test_submit_project(&client1_base, &client2_email, &test_mode)?;
    pause_between_steps(5, "Test 5: Submit project from client1 to client2");

    // Wait for submission to sync
    thread::sleep(Duration::from_secs(10));

    // Test 6: Check client2 received the submission
    println!("\n📝 Test 6: Check client2 received submission");
    test_check_submission(&client2_base, &client1_email, &test_mode)?;
    pause_between_steps(6, "Test 6: Check client2 received submission");

    // Test 7: Client2 processes the request on imported data
    println!("\n📝 Test 7: Client2 processes request on imported data");
    let process_outcome = test_process_request(
        &client2_base,
        &client1_email,
        expected_file_count,
        &test_mode,
    )?;
    pause_between_steps(7, "Test 7: Client2 processes request on imported data");

    println!("\n📝 Test 8: Approve project to return results");
    let approved_csv = approve_project_request(
        &client2_base,
        &client1_base,
        &client1_email,
        &process_outcome.submission_folder,
        &test_mode,
    )?;

    let final_results = approved_csv.unwrap_or(process_outcome.run_results.clone());
    pause_between_steps(8, "Test 8: Approve project to return results");

    // Test 9: Validate processor results
    println!("\n📝 Test 9: Validate processor results");
    test_validate_processor_results(&final_results, expected_file_count)?;
    pause_between_steps(9, "Test 9: Validate processor results");

    // Test 10: Archive and remove permissions
    println!("\n📝 Test 10: Archive and remove permissions");
    test_archive_project(&client1_base, &client2_email, &test_mode)?;
    pause_between_steps(10, "Test 10: Archive and remove permissions");

    println!(
        "
  📬 Dumping RPC request/response JSON files for debugging..."
    );
    dump_rpc_messages(&client1_base.join("datasites"));
    dump_rpc_messages(&client2_base.join("datasites"));
    println!("\n✅ All BioVault tests completed successfully!");
    Ok(())
}

fn test_biovault_init(client_base: &Path, email: &str, test_mode: &str) -> Result<()> {
    // Ensure the client base directory exists
    if !client_base.exists() {
        std::fs::create_dir_all(client_base)?;
    }

    // Clean up any previous BioVault config
    let biovault_home = client_base.join(".biovault");
    if biovault_home.exists() {
        fs::remove_dir_all(&biovault_home)?;
        println!("Cleaned up existing BioVault config");
    }

    // Run bv init with --quiet flag to auto-accept prompts
    let output = run_bv_command(client_base, test_mode, &["init", "--quiet", email])?;
    let stdout_str = String::from_utf8_lossy(&output.stdout);
    println!("bv init output: {}", stdout_str);
    if !output.status.success() {
        eprintln!(
            "bv init stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Err(anyhow::anyhow!("bv init failed"));
    }

    // Run bv config to verify the setup
    let config_output = run_bv_command(client_base, test_mode, &["config"])?;
    let config_str = String::from_utf8_lossy(&config_output.stdout);
    println!("bv config output:\n{}", config_str);

    // Verify config contains expected values
    assert!(
        config_str.contains(email),
        "Config should contain email: {}",
        email
    );
    assert!(
        config_str.contains("BioVault Configuration"),
        "Config should show BioVault configuration"
    );
    assert!(
        config_str.contains("SyftBox data dir:"),
        "Config should show SyftBox data dir"
    );

    // Only check RPC directory if init actually ran (not skipped)
    if !stdout_str.contains("Skipping initialization") {
        // Check that RPC directory was created with correct permissions
        let rpc_dir = client_base
            .join("datasites")
            .join(email)
            .join("app_data")
            .join("biovault")
            .join("rpc");

        // Debug: show what directories exist
        println!("Looking for RPC dir at: {}", rpc_dir.display());
        println!("client_base exists: {}", client_base.exists());
        if client_base.exists() {
            println!("Contents of client_base:");
            if let Ok(entries) = fs::read_dir(client_base) {
                for entry in entries.flatten() {
                    println!("  - {}", entry.path().display());
                }
            }
        }

        assert!(
            rpc_dir.exists(),
            "RPC directory should exist at: {}",
            rpc_dir.display()
        );

        let permission_file = rpc_dir.join("syft.pub.yaml");
        assert!(permission_file.exists(), "RPC permission file should exist");

        // Verify permission file format
        // Note: SyftBox may reformat the file, so we just check that it has the expected pattern
        let content = fs::read_to_string(&permission_file)?;
        assert!(
            content.contains("- pattern:") || content.contains("  - pattern:"),
            "Permission file should have pattern rules"
        );
    }

    println!("✓ BioVault initialized for {}", email);
    Ok(())
}

fn test_import_genotype_data(
    client_base: &Path,
    genotype_dir: &Path,
    test_mode: &str,
) -> Result<usize> {
    if !genotype_dir.exists() {
        return Err(anyhow!(
            "Genotype fixtures folder missing: {}",
            genotype_dir.display()
        ));
    }

    let client_base_abs = canonicalize_path(client_base);
    let csv_path = client_base_abs.join("genotype_catalog.csv");
    if csv_path.exists() {
        fs::remove_file(&csv_path)?;
    }

    // Export fixtures to CSV with participant IDs derived from filenames
    let genotype_dir_abs = canonicalize_path(genotype_dir);
    let export_args = [
        "files".to_string(),
        "export-csv".to_string(),
        genotype_dir_abs.to_string_lossy().to_string(),
        "--pattern".to_string(),
        "{filename}".to_string(),
        "-o".to_string(),
        csv_path.to_string_lossy().to_string(),
    ];
    let export_slices: Vec<&str> = export_args.iter().map(|s| s.as_str()).collect();
    let export_output = run_bv_command(client_base, test_mode, &export_slices)?;
    let export_stdout = String::from_utf8_lossy(&export_output.stdout);
    let export_stderr = String::from_utf8_lossy(&export_output.stderr);
    println!("files export-csv stdout:\n{}", export_stdout);
    if !export_output.status.success() {
        eprintln!("files export-csv stderr:\n{}", export_stderr);
        return Err(anyhow!("bv files export-csv failed"));
    }
    assert!(
        csv_path.exists(),
        "Exported CSV should exist at {}",
        csv_path.display()
    );

    // Detect metadata for exported files
    let detect_args = [
        "files".to_string(),
        "detect-csv".to_string(),
        csv_path.to_string_lossy().to_string(),
        "-o".to_string(),
        csv_path.to_string_lossy().to_string(),
    ];
    let detect_slices: Vec<&str> = detect_args.iter().map(|s| s.as_str()).collect();
    let detect_output = run_bv_command(client_base, test_mode, &detect_slices)?;
    let detect_stdout = String::from_utf8_lossy(&detect_output.stdout);
    let detect_stderr = String::from_utf8_lossy(&detect_output.stderr);
    println!("files detect-csv stdout:\n{}", detect_stdout);
    if !detect_output.status.success() {
        eprintln!("files detect-csv stderr:\n{}", detect_stderr);
        return Err(anyhow!("bv files detect-csv failed"));
    }

    // Import the detected CSV into the database catalog
    let import_args = [
        "files".to_string(),
        "import-csv".to_string(),
        csv_path.to_string_lossy().to_string(),
        "--non-interactive".to_string(),
        "--format".to_string(),
        "json".to_string(),
    ];
    let import_slices: Vec<&str> = import_args.iter().map(|s| s.as_str()).collect();
    let import_output = run_bv_command(client_base, test_mode, &import_slices)?;
    let import_stdout = String::from_utf8_lossy(&import_output.stdout);
    println!("files import-csv stdout:\n{}", import_stdout);
    if !import_output.status.success() {
        eprintln!(
            "files import-csv stderr:\n{}",
            String::from_utf8_lossy(&import_output.stderr)
        );
        return Err(anyhow!("bv files import-csv failed"));
    }

    // List cataloged files to confirm
    let list_output = run_bv_command(
        client_base,
        test_mode,
        &["files", "list", "--format", "json"],
    )?;
    if !list_output.status.success() {
        eprintln!(
            "files list stderr:\n{}",
            String::from_utf8_lossy(&list_output.stderr)
        );
        return Err(anyhow!("bv files list failed"));
    }
    let list_stdout = String::from_utf8_lossy(&list_output.stdout);
    let files_data: FilesListData =
        parse_json_output(&list_stdout, "Failed to parse files list JSON output")?;
    println!("  Catalog now tracks {} genotype files", files_data.total);
    assert!(
        !files_data.files.is_empty(),
        "Catalog should contain imported genotype files"
    );

    Ok(files_data.total)
}

fn test_prepare_dynamic_project(client_base: &Path, _test_mode: &str) -> Result<()> {
    let client_base_abs = canonicalize_path(client_base);
    let project_dir = client_base_abs.join(COUNT_LINES_PROJECT);
    if project_dir.exists() {
        fs::remove_dir_all(&project_dir)?;
    }

    if project_dir.exists() {
        fs::remove_dir_all(&project_dir).with_context(|| {
            format!(
                "failed to remove existing project directory {}",
                project_dir.display()
            )
        })?;
    }

    fs::create_dir_all(project_dir.join("assets"))
        .with_context(|| format!("failed to create assets dir at {}", project_dir.display()))?;
    fs::write(project_dir.join("project.yaml"), COUNT_LINES_PROJECT_YAML)
        .with_context(|| format!("failed to write project.yaml at {}", project_dir.display()))?;
    fs::write(project_dir.join("workflow.nf"), COUNT_LINES_WORKFLOW)
        .with_context(|| format!("failed to write workflow.nf at {}", project_dir.display()))?;
    fs::write(
        project_dir.join("assets/count_lines.py"),
        COUNT_LINES_SCRIPT,
    )
    .with_context(|| {
        format!(
            "failed to write count_lines.py at {}",
            project_dir.display()
        )
    })?;

    assert!(
        project_dir.exists(),
        "Project directory should exist after creation"
    );
    assert!(
        project_dir.join("workflow.nf").exists(),
        "workflow.nf should exist"
    );
    assert!(
        project_dir.join("project.yaml").exists(),
        "project.yaml should exist"
    );

    println!("✓ Project prepared: {}", COUNT_LINES_PROJECT);
    Ok(())
}

fn test_submit_project(client_base: &Path, recipient_email: &str, test_mode: &str) -> Result<()> {
    // Submit the dynamic count-lines project to the recipient
    let project_dir = client_base.join(COUNT_LINES_PROJECT);
    let project_dir_abs = canonicalize_path(&project_dir);

    if !project_dir_abs.exists() {
        println!("⚠ Project directory not found: {}", project_dir.display());
        return Ok(());
    }

    println!(
        "Submitting project from directory: {}",
        client_base.display()
    );
    println!("Project directory exists: {}", project_dir_abs.exists());

    let project_arg = project_dir_abs.to_string_lossy().to_string();
    let submit_args = [
        "submit",
        project_arg.as_str(),
        recipient_email,
        "--non-interactive",
        "--force",
    ];
    let output = run_bv_command(client_base, test_mode, &submit_args)?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    println!("Submit stdout: {}", stdout);
    if !stderr.is_empty() {
        println!("Submit stderr: {}", stderr);
    }

    if !output.status.success() {
        eprintln!(
            "Submit failed with exit code: {}",
            output.status.code().unwrap_or(-1)
        );
        return Err(anyhow::anyhow!("project submit failed"));
    }

    // Check for expected output
    if stdout.contains("Project submitted successfully")
        || stdout.contains("Project message prepared")
    {
        println!("✓ Project submitted to {}", recipient_email);

        // Extract submission location if present
        if let Some(location_line) = stdout.lines().find(|l| l.contains("Location:")) {
            println!("  {}", location_line.trim());
        }

        // Verify message was sent by checking the sent messages
        println!("\n  Verifying message was sent...");
        let msg_output = run_bv_command(client_base, test_mode, &["message", "list", "--sent"])?;

        let msg_stdout = String::from_utf8_lossy(&msg_output.stdout);
        if msg_stdout.contains(recipient_email) && msg_stdout.contains("Project Request") {
            println!("  ✓ Message successfully sent to {}", recipient_email);
        } else {
            println!("  ⚠ Could not verify message was sent");
            println!("  Messages output: {}", msg_stdout);
        }
    } else {
        println!("⚠ Submission may have succeeded but output unclear");
    }

    Ok(())
}

fn test_check_submission(client_base: &Path, sender_email: &str, test_mode: &str) -> Result<()> {
    // Check if submission folder was created in client2's datasite
    let submissions_path = client_base
        .join("datasites")
        .join(sender_email)
        .join("shared")
        .join("biovault")
        .join("submissions");

    if submissions_path.exists() {
        println!(
            "✓ Submission folder exists at: {}",
            submissions_path.display()
        );

        // List submission folders
        let mut found_submission = false;
        for entry in fs::read_dir(&submissions_path)? {
            let entry = entry?;
            if entry.file_type()?.is_dir() {
                let submission_name = entry.file_name().to_string_lossy().to_string();
                if submission_name.starts_with(COUNT_LINES_PROJECT_NAME) {
                    println!("  Found submission: {}", submission_name);

                    // Check for expected files
                    let submission_dir = entry.path();
                    let workflow_file = submission_dir.join("workflow.nf");
                    let project_file = submission_dir.join("project.yaml");
                    let permission_file = submission_dir.join("syft.pub.yaml");

                    if workflow_file.exists() && project_file.exists() && permission_file.exists() {
                        println!("    ✓ All expected files present");
                        found_submission = true;
                    } else {
                        println!("    ⚠ Some expected files missing");
                    }
                }
            }
        }

        if !found_submission {
            println!("  ⚠ No {} submission found yet", COUNT_LINES_PROJECT_NAME);
        }
    } else {
        println!(
            "⚠ Submission folder not found at: {}",
            submissions_path.display()
        );
        println!("  Messages may not have synced yet");
    }

    println!("\n  Checking message state via CLI...");
    let (maybe_message_id, list_stdout) = match wait_for_project_message_id(
        client_base,
        test_mode,
        sender_email,
        true,
        MESSAGE_WAIT_ATTEMPTS,
        Duration::from_secs(MESSAGE_WAIT_DELAY_SECS),
        "check_submission_message",
    ) {
        Ok(result) => result,
        Err(err) => {
            println!("  ⚠ Unable to query unread project messages: {}", err);
            return Ok(());
        }
    };

    if let Some(message_id) = maybe_message_id {
        println!(
            "  Found unread project message [{}] from {}",
            message_id, sender_email
        );

        let read_output = run_bv_command(
            client_base,
            test_mode,
            &["message", "read", &message_id, "--non-interactive"],
        )?;

        if read_output.status.success() {
            println!("  ✓ Marked message as read via CLI");

            let verify_output = run_bv_command(
                client_base,
                test_mode,
                &["message", "list", "--projects", "--unread", "--json"],
            )?;

            if verify_output.status.success() {
                let verify_stdout = String::from_utf8_lossy(&verify_output.stdout);
                if find_project_message_id(&verify_stdout, sender_email).is_some() {
                    println!(
                        "  ⚠ Message still appears unread after read command. CLI output:\n{}",
                        verify_stdout
                    );
                } else {
                    println!("  ✓ Message no longer appears in unread list");
                }
            } else {
                println!(
                    "  ⚠ Unable to verify unread messages after read: {}",
                    String::from_utf8_lossy(&verify_output.stderr)
                );
            }
        } else {
            println!(
                "  ⚠ Failed to mark message as read: {}",
                String::from_utf8_lossy(&read_output.stderr)
            );
        }
    } else {
        println!(
            "  ⚠ Could not locate unread project message for {} in CLI output:\n{}",
            sender_email, list_stdout
        );
    }

    Ok(())
}

fn test_process_request(
    client_base: &Path,
    sender_email: &str,
    expected_count: usize,
    test_mode: &str,
) -> Result<ProcessOutcome> {
    println!("Generating samplesheet from imported catalog...");
    let client_base_abs = canonicalize_path(client_base);
    let samplesheet_path = client_base_abs.join(ALL_SAMPLESHEET_NAME);
    if samplesheet_path.exists() {
        fs::remove_file(&samplesheet_path)?;
    }
    let catalog_count =
        generate_samplesheet_from_db(&client_base_abs, test_mode, &samplesheet_path)?;
    println!("  Catalog returned {} files", catalog_count);
    assert_eq!(
        catalog_count, expected_count,
        "Catalog entries should match staged genotype files"
    );

    let submissions_root = client_base_abs
        .join("datasites")
        .join(sender_email)
        .join("shared")
        .join("biovault")
        .join("submissions");

    if let Err(err) = run_bv_command(client_base, test_mode, &["message", "sync"]) {
        println!("⚠ Failed to sync messages before processing: {}", err);
    }

    if !submissions_root.exists() {
        return Err(anyhow!(
            "Submission folder not found at {}",
            submissions_root.display()
        ));
    }

    let mut submissions: Vec<(String, PathBuf)> = Vec::new();
    for entry in fs::read_dir(&submissions_root)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            let name = entry.file_name().to_string_lossy().to_string();
            if name.starts_with(COUNT_LINES_PROJECT_NAME) {
                submissions.push((name, entry.path()));
            }
        }
    }

    if submissions.is_empty() {
        println!(
            "⚠ No {} submission found yet, waiting for sync...",
            COUNT_LINES_PROJECT_NAME
        );
        let awaited_path = wait_for_submission_directory(
            client_base,
            test_mode,
            &submissions_root,
            COUNT_LINES_PROJECT_NAME,
            30,
            Duration::from_secs(2),
        )?;
        let submission_name = awaited_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| COUNT_LINES_PROJECT_NAME.to_string());
        println!("Using awaited submission directory: {}", submission_name);
        submissions.push((submission_name, awaited_path));
    }

    submissions.sort_by(|a, b| a.0.cmp(&b.0));
    let (submission_name, _submission_dir) = submissions.pop().unwrap();
    println!("Using submission directory: {}", submission_name);
    let (maybe_message_id, list_stdout) = wait_for_project_message_id(
        client_base,
        test_mode,
        sender_email,
        false,
        MESSAGE_WAIT_ATTEMPTS,
        Duration::from_secs(MESSAGE_WAIT_DELAY_SECS),
        "process_request",
    )?;
    let message_id = maybe_message_id.ok_or_else(|| {
        anyhow!(
            "Could not locate project message for {} in CLI output:\n{}",
            sender_email,
            list_stdout
        )
    })?;

    let process_args = [
        "message",
        "process",
        &message_id,
        "--test",
        "--participant",
        "ALL",
        "--non-interactive",
    ];
    let process_output = run_bv_command(&client_base_abs, test_mode, &process_args)?;
    let process_stdout = String::from_utf8_lossy(&process_output.stdout);
    println!("processor CLI stdout:\n{}", process_stdout);
    if !process_output.status.success() {
        return Err(anyhow!(
            "Processing project via CLI failed: {}",
            String::from_utf8_lossy(&process_output.stderr)
        ));
    }

    let results_dir = parse_results_dir_from_stdout(&process_stdout)
        .or_else(|| {
            let fallback = private_results_dir(&client_base_abs, &submission_name, "ALL", true);
            if fallback.exists() {
                Some(fallback)
            } else {
                None
            }
        })
        .ok_or_else(|| anyhow!("Unable to determine processor results directory"))?;

    let run_results = results_dir.join("line_counts.csv");
    if !run_results.exists() {
        return Err(anyhow!(
            "Processor results file not found at {}",
            run_results.display()
        ));
    }

    Ok(ProcessOutcome {
        run_results,
        submission_folder: submission_name,
    })
}

fn wait_for_project_message_id(
    client_base: &Path,
    test_mode: &str,
    sender_email: &str,
    unread_only: bool,
    attempts: usize,
    delay: Duration,
    context: &str,
) -> Result<(Option<String>, String)> {
    let mut last_stdout = String::new();
    for attempt in 1..=attempts {
        println!(
            "  [message-wait:{context}] Attempt {}/{} (unread_only={})",
            attempt, attempts, unread_only
        );

        let sync_output = run_bv_command(client_base, test_mode, &["message", "sync"])?;
        if !sync_output.status.success() {
            println!(
                "  [message-wait:{context}] 'bv message sync' exited with {}",
                sync_output.status
            );
            let sync_stderr = String::from_utf8_lossy(&sync_output.stderr);
            if !sync_stderr.trim().is_empty() {
                println!("  [message-wait:{context}] sync stderr:\n{}", sync_stderr);
            }
        }

        let mut args = vec!["message", "list", "--projects", "--json"];
        if unread_only {
            args.insert(3, "--unread");
        }
        let list_output = run_bv_command(client_base, test_mode, &args)?;
        let list_stdout = String::from_utf8_lossy(&list_output.stdout).to_string();
        last_stdout = list_stdout.clone();
        let list_stderr = String::from_utf8_lossy(&list_output.stderr);
        println!(
            "  [message-wait:{context}] 'bv {}' exited with {}",
            args.join(" "),
            list_output.status
        );
        if !list_stderr.trim().is_empty() {
            println!("  [message-wait:{context}] list stderr:\n{}", list_stderr);
        }

        if !list_output.status.success() {
            return Err(anyhow!(
                "'bv {}' failed with status {}. Stdout:\n{}\nStderr:\n{}",
                args.join(" "),
                list_output.status,
                list_stdout,
                list_stderr
            ));
        }

        let entry_count = serde_json::from_str::<serde_json::Value>(list_stdout.trim())
            .ok()
            .and_then(|value| value.as_array().map(|arr| arr.len()))
            .unwrap_or(0);
        println!(
            "  [message-wait:{context}] Parsed {} project entries ({} bytes)",
            entry_count,
            list_stdout.trim().len()
        );

        if let Some(message_id) = find_project_message_id(&list_stdout, sender_email) {
            println!(
                "  [message-wait:{context}] Found message {} from {}",
                message_id, sender_email
            );
            return Ok((Some(message_id), list_stdout));
        }

        if attempt < attempts {
            println!(
                "  [message-wait:{context}] Message not visible yet, sleeping {:?} before retry",
                delay
            );
            thread::sleep(delay);
        }
    }

    println!(
        "  [message-wait:{context}] Exhausted {} attempts without locating message from {}",
        attempts, sender_email
    );
    Ok((None, last_stdout))
}

fn wait_for_submission_directory(
    client_base: &Path,
    test_mode: &str,
    submissions_root: &Path,
    prefix: &str,
    attempts: usize,
    delay: Duration,
) -> Result<PathBuf> {
    for attempt in 0..attempts {
        if submissions_root.exists() {
            let mut matches: Vec<(String, PathBuf)> = Vec::new();
            for entry in fs::read_dir(submissions_root)? {
                let entry = entry?;
                if entry.file_type()?.is_dir() {
                    let name = entry.file_name().to_string_lossy().to_string();
                    if name.starts_with(prefix) {
                        matches.push((name, entry.path()));
                    }
                }
            }

            if !matches.is_empty() {
                matches.sort_by(|a, b| a.0.cmp(&b.0));
                return Ok(matches.pop().unwrap().1);
            }
        }

        if attempt + 1 < attempts {
            if let Err(err) = run_bv_command(client_base, test_mode, &["message", "sync"]) {
                println!("⚠ Failed to sync messages while waiting: {}", err);
            }
            thread::sleep(delay);
        }
    }

    Err(anyhow!(
        "Timed out waiting for {} submission directory under {}",
        prefix,
        submissions_root.display()
    ))
}

fn test_validate_processor_results(results_csv: &Path, expected_count: usize) -> Result<()> {
    println!("Examining processor results at {}", results_csv.display());
    if !results_csv.exists() {
        return Err(anyhow!(
            "Results file not found at {}",
            results_csv.display()
        ));
    }

    let contents = fs::read_to_string(results_csv)?;
    let mut lines = contents.lines();
    let header = lines.next().unwrap_or("");
    assert!(
        header.contains("participant_id") && header.contains("line_count"),
        "Results header should include participant_id and line_count"
    );

    let mut rows = 0usize;
    let mut non_zero = 0usize;
    for (idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        rows += 1;
        if idx < 5 {
            println!("  row{}: {}", idx + 1, line);
        }
        if let Some(count_field) = line.split(',').next_back() {
            if count_field.trim().parse::<usize>().unwrap_or(0) > 0 {
                non_zero += 1;
            }
        }
    }

    println!("  Results contain {} data rows", rows);
    assert_eq!(
        rows, expected_count,
        "Results should include one row per genotype file"
    );
    assert!(non_zero > 0, "At least one line_count should be non-zero");

    Ok(())
}

fn test_archive_project(client_base: &Path, _other_email: &str, test_mode: &str) -> Result<()> {
    // Archive the project message to revoke write permissions

    // Find the project message to archive
    println!("Looking for project message to archive...");
    let list_output = run_bv_command(
        client_base,
        test_mode,
        &["message", "list", "--sent", "--projects"],
    );

    if let Ok(output) = list_output {
        let stdout = String::from_utf8_lossy(&output.stdout);

        // Find the dynamic project message
        if stdout.contains(COUNT_LINES_PROJECT_NAME) {
            // Extract message ID
            let msg_id =
                if let Some(line) = stdout.lines().find(|l| l.contains("[") && l.contains("]")) {
                    if let Some(start) = line.find('[') {
                        line.find(']').map(|end| line[start + 1..end].to_string())
                    } else {
                        None
                    }
                } else {
                    None
                };

            if let Some(id) = msg_id {
                println!("Found project message to archive: {}", id);

                // Archive the message
                let archive_output =
                    run_bv_command(client_base, test_mode, &["message", "archive", &id]);

                if let Ok(output) = archive_output {
                    if output.status.success() {
                        println!("✓ Project archived - write permissions revoked");

                        // Verify permissions were updated
                        let submissions_path = client_base
                            .join("datasites")
                            .join(
                                client_base
                                    .file_name()
                                    .unwrap()
                                    .to_string_lossy()
                                    .to_string(),
                            )
                            .join("shared")
                            .join("biovault")
                            .join("submissions");

                        if submissions_path.exists() {
                            for entry in fs::read_dir(&submissions_path)?.flatten() {
                                if entry.file_type()?.is_dir() {
                                    let name = entry.file_name().to_string_lossy().to_string();
                                    if name.starts_with(COUNT_LINES_PROJECT_NAME) {
                                        let perm_file = entry.path().join("syft.pub.yaml");
                                        if perm_file.exists() {
                                            let content = fs::read_to_string(&perm_file)?;
                                            if !content.contains("results/**/*") {
                                                println!("  ✓ Verified: results write rule removed from permissions");
                                            } else {
                                                println!(
                                                    "  ⚠ Warning: results write rule still present"
                                                );
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        let stderr = String::from_utf8_lossy(&output.stderr);
                        println!("⚠ Archive failed: {}", stderr);
                    }
                }
            } else {
                println!("⚠ Could not find message ID to archive");
            }
        } else {
            println!("⚠ No project message found to archive");
        }
    }

    Ok(())
}

fn stage_genotype_fixtures(client_base: &Path) -> Result<PathBuf> {
    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(GENOTYPE_FIXTURES_DIR);
    if !source.exists() {
        return Err(anyhow!("Fixture directory not found: {}", source.display()));
    }

    let client_base_abs = canonicalize_path(client_base);
    let target = client_base_abs.join("genotype_fixtures");
    if target.exists() {
        fs::remove_dir_all(&target)?;
    }
    copy_dir_recursive(&source, &target)?;
    Ok(canonicalize_path(&target))
}

fn count_files_in_dir(dir: &Path) -> Result<usize> {
    if !dir.exists() {
        return Ok(0);
    }

    let mut count = 0usize;
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if entry.file_type()?.is_file() {
            count += 1;
        }
    }
    Ok(count)
}

fn copy_dir_recursive(src: &Path, dest: &Path) -> Result<()> {
    if dest.exists() {
        fs::remove_dir_all(dest)?;
    }
    fs::create_dir_all(dest)?;

    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let path = entry.path();
        let target_path = dest.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            copy_dir_recursive(&path, &target_path)?;
        } else {
            if let Some(parent) = target_path.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::copy(&path, &target_path)?;
        }
    }
    Ok(())
}

fn derive_participant_from_path(path: &Path) -> String {
    path.file_stem()
        .and_then(|s| s.to_str())
        .or_else(|| path.file_name().and_then(|s| s.to_str()))
        .unwrap_or("unknown")
        .to_string()
}

fn generate_samplesheet_from_db(
    client_base: &Path,
    test_mode: &str,
    output_path: &Path,
) -> Result<usize> {
    let output = run_bv_command(
        client_base,
        test_mode,
        &["files", "list", "--format", "json"],
    )?;
    if !output.status.success() {
        eprintln!(
            "files list stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Err(anyhow!("bv files list failed"));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let response: FilesListData =
        parse_json_output(&stdout, "Failed to parse files list JSON output")?;

    if response.files.is_empty() {
        println!("⚠ Catalog returned no files");
    }

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file = fs::File::create(output_path)?;
    writeln!(file, "participant_id,genotype_file_path")?;

    let mut rows = Vec::new();
    for entry in &response.files {
        let participant = entry
            .participant_name
            .as_ref()
            .or(entry.participant_id.as_ref())
            .cloned()
            .unwrap_or_else(|| derive_participant_from_path(Path::new(&entry.file_path)));

        let canonical = Path::new(&entry.file_path)
            .canonicalize()
            .unwrap_or_else(|_| PathBuf::from(&entry.file_path));

        rows.push((participant, canonical.to_string_lossy().to_string()));
    }

    rows.sort_by(|a, b| a.0.cmp(&b.0));
    for (participant, path) in &rows {
        writeln!(file, "{},{}", participant, path)?;
    }

    println!(
        "  Samplesheet generated at {} with {} rows",
        output_path.display(),
        rows.len()
    );

    Ok(rows.len())
}

fn canonicalize_path(path: &Path) -> PathBuf {
    path.canonicalize().unwrap_or_else(|_| path.to_path_buf())
}

fn get_bv_binary_path() -> Result<PathBuf> {
    // Try multiple possible locations for the binary
    // In CI, tests might run from different directories
    let possible_paths = vec![
        PathBuf::from("cli/target/release/bv"),
        PathBuf::from("target/release/bv"),
        PathBuf::from("../target/release/bv"),
        PathBuf::from("../../target/release/bv"),
        PathBuf::from("./bv"), // Sometimes the binary is copied to current dir
    ];

    for path in &possible_paths {
        if path.exists() {
            return Ok(path.canonicalize()?);
        }
    }

    // If not found, print debugging info
    eprintln!("Debug: Current directory is {:?}", std::env::current_dir());
    eprintln!("Debug: Searched for binary in:");
    for path in &possible_paths {
        eprintln!("  - {:?} (exists: {})", path, path.exists());
    }

    // Also check if we can find it via cargo metadata
    if let Ok(output) = Command::new("cargo")
        .args(["metadata", "--format-version", "1", "--no-deps"])
        .output()
    {
        if output.status.success() {
            // Try to parse target directory from cargo metadata
            let metadata = String::from_utf8_lossy(&output.stdout);
            if let Some(target_start) = metadata.find("\"target_directory\":") {
                let target_rest = &metadata[target_start + 20..];
                if let Some(end_quote) = target_rest.find('"') {
                    let target_dir = &target_rest[1..end_quote];
                    let binary_path = PathBuf::from(target_dir).join("release/bv");
                    if binary_path.exists() {
                        return Ok(binary_path.canonicalize()?);
                    }
                }
            }
        }
    }

    Err(anyhow::anyhow!(
        "BioVault binary not found. Run 'cargo build --release' first"
    ))
}

fn run_bv_command(
    client_base: &Path,
    test_mode: &str,
    args: &[&str],
) -> Result<std::process::Output> {
    let client_base_abs = canonicalize_path(client_base);
    // Copy bv binary to temp location to avoid path issues
    let bv_source = get_bv_binary_path()?;
    let bv_temp = PathBuf::from("/tmp/bv-test");
    fs::copy(&bv_source, &bv_temp)?;

    let mut cmd = if test_mode == "local" {
        // For local mode, run bv through sbenv which sets all the right env vars
        // Use the sbenv shell shim that will auto-rebuild if needed
        let sbenv_path = client_base_abs
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("sbenv/sbenv")
            .canonicalize()
            .unwrap_or_else(|_| {
                // Try to find sbenv in the parent directories
                let mut current = client_base_abs.clone();
                loop {
                    let candidate = current.join("sbenv/sbenv");
                    if candidate.exists() {
                        return candidate;
                    }
                    if !current.pop() {
                        break;
                    }
                }
                // Last resort: assume sbenv is in PATH
                PathBuf::from("sbenv")
            });

        // Create a script that activates sbenv and runs bv with args
        // Don't set BIOVAULT_HOME - let bv use SYFTBOX_DATA_DIR/.biovault naturally
        let args_str = args
            .iter()
            .map(|a| format!("'{}'", a))
            .collect::<Vec<_>>()
            .join(" ");
        let script = format!(
            r#"set -e
cd '{}'
eval "$('{}' activate --quiet)"
echo "Debug: SYFTBOX_CONFIG_PATH=$SYFTBOX_CONFIG_PATH" >&2
echo "Debug: SYFTBOX_DATA_DIR=$SYFTBOX_DATA_DIR" >&2
echo "Debug: Running bv with args: {}" >&2
'{}' {} < /dev/null"#,
            client_base_abs.display(),
            sbenv_path.display(),
            args_str,
            bv_temp.display(),
            args_str
        );

        let mut cmd = Command::new("bash");
        cmd.arg("-c").arg(&script);
        cmd
    } else {
        // For Docker mode, the SyftBox directory structure is:
        // test-clients-docker/client1@syftbox.net/SyftBox/
        // We want .biovault to be created at test-clients-docker/client1@syftbox.net/SyftBox/.biovault
        // So we set SYFTBOX_DATA_DIR to the SyftBox directory itself
        let mut cmd = Command::new(&bv_temp);
        cmd.current_dir(&client_base_abs);
        // Use absolute path to avoid any relative path confusion
        let abs_client_base = client_base_abs.clone();
        cmd.env("SYFTBOX_DATA_DIR", &abs_client_base);
        cmd.args(args);
        cmd
    };

    Ok(cmd.output()?)
}

fn find_project_message_id(output: &str, sender_email: &str) -> Option<String> {
    let trimmed = output.trim();
    if trimmed.is_empty() {
        return None;
    }

    let json_value: serde_json::Value = serde_json::from_str(trimmed).ok()?;
    let messages = json_value.as_array()?;

    for msg in messages {
        let from = msg.get("from").and_then(|v| v.as_str()).unwrap_or("");
        if !from.eq_ignore_ascii_case(sender_email) {
            continue;
        }

        let msg_type = msg
            .get("message_type")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .to_lowercase();
        if msg_type != "project" {
            continue;
        }

        if let Some(id) = msg.get("id").and_then(|v| v.as_str()) {
            return Some(id.to_string());
        }
    }

    None
}

fn parse_results_dir_from_stdout(stdout: &str) -> Option<PathBuf> {
    stdout.lines().find_map(|line| {
        let trimmed = line.trim();
        trimmed
            .strip_prefix("Results saved to: ")
            .map(|path| PathBuf::from(path.trim()))
    })
}

fn private_results_dir(
    client_base: &Path,
    submission_folder: &str,
    participant: &str,
    test_run: bool,
) -> PathBuf {
    let results_root = client_base
        .join("private")
        .join("app_data")
        .join("biovault")
        .join("submissions")
        .join(submission_folder);
    let results_dir_name = if test_run {
        "results-test"
    } else {
        "results-real"
    };
    results_root.join(results_dir_name).join(participant)
}

fn dump_rpc_messages(datasites_root: &Path) {
    let app_path = datasites_root
        .join("app_data")
        .join("biovault")
        .join("rpc")
        .join("message");
    if !app_path.exists() {
        return;
    }
    println!("  Inspecting RPC envelope dir: {}", app_path.display());
    for entry in WalkDir::new(&app_path).into_iter().filter_map(|e| e.ok()) {
        if entry.file_type().is_file() {
            let path = entry.path();
            let ext = path.extension().and_then(|e| e.to_str());
            if matches!(ext, Some("json") | Some("request") | Some("response")) {
                println!("    File: {}", path.display());
                match fs::read_to_string(path) {
                    Ok(content) => match serde_json::from_str::<serde_json::Value>(&content) {
                        Ok(json) => {
                            println!("{}", serde_json::to_string_pretty(&json).unwrap_or(content))
                        }
                        Err(_) => println!("{}", content),
                    },
                    Err(err) => println!("    ⚠  Could not read {}: {}", path.display(), err),
                }
            }
        }
    }
}

fn approve_project_request(
    client_base: &Path,
    sender_base: &Path,
    sender_email: &str,
    submission_folder: &str,
    test_mode: &str,
) -> Result<Option<PathBuf>> {
    // Ensure latest message state before inspecting
    if let Err(err) = run_bv_command(client_base, test_mode, &["message", "sync"]) {
        println!("⚠️  Failed to sync messages before approval: {}", err);
    }

    let inbox_output = run_bv_command(
        client_base,
        test_mode,
        &["message", "list", "--projects", "--json"],
    )?;
    if !inbox_output.status.success() {
        println!(
            "⚠️  Unable to list project messages: {}",
            String::from_utf8_lossy(&inbox_output.stderr)
        );
        return Ok(None);
    }
    let stdout = String::from_utf8_lossy(&inbox_output.stdout);

    let message_id = find_project_message_id(&stdout, sender_email);

    if let Some(id) = message_id {
        let approve_args = [
            "message",
            "process",
            &id,
            "--real",
            "--participant",
            "ALL",
            "--approve",
            "--non-interactive",
        ];
        let output = run_bv_command(client_base, test_mode, &approve_args)?;
        if output.status.success() {
            println!("  ✓ Approved project {}", sender_email);
        } else {
            println!(
                "  ⚠  Approval may have failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            return Ok(None);
        }
    } else {
        println!(
            "⚠️  Could not find message ID to approve for {}",
            sender_email
        );
        return Ok(None);
    }

    let sender_results = sender_base
        .join("datasites")
        .join(sender_email)
        .join("shared")
        .join("biovault")
        .join("submissions")
        .join(submission_folder)
        .join("results")
        .join("line_counts.csv");

    if sender_results.exists() {
        println!("  ✓ Results available at {}", sender_results.display());
        Ok(Some(sender_results))
    } else {
        println!(
            "  ⚠  Results file not yet found at {}",
            sender_results.display()
        );
        Ok(None)
    }
}
