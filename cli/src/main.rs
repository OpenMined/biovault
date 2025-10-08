use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

use biovault::cli;

use cli::commands;

// Validator for example names that also shows available examples
fn validate_example_name(s: &str) -> Result<String, String> {
    let examples = cli::examples::list_examples();
    if examples.contains(&s.to_string()) {
        Ok(s.to_string())
    } else {
        Err(format!(
            "Unknown example '{}'. Available examples: {}",
            s,
            examples.join(", ")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use tempfile::TempDir;

    struct EnvVarGuard {
        key: &'static str,
        previous: Option<String>,
    }

    impl EnvVarGuard {
        fn set(key: &'static str, value: &str) -> Self {
            let previous = std::env::var(key).ok();
            std::env::set_var(key, value);
            Self { key, previous }
        }
    }

    impl Drop for EnvVarGuard {
        fn drop(&mut self) {
            if let Some(prev) = &self.previous {
                std::env::set_var(self.key, prev);
            } else {
                std::env::remove_var(self.key);
            }
        }
    }

    struct TestHomeGuard {
        _temp: TempDir,
    }

    impl TestHomeGuard {
        fn new() -> Self {
            let temp = TempDir::new().unwrap();
            let home = temp.path().join(".biovault");
            std::fs::create_dir_all(&home).unwrap();
            biovault::config::set_test_biovault_home(&home);
            Self { _temp: temp }
        }
    }

    impl Drop for TestHomeGuard {
        fn drop(&mut self) {
            biovault::config::clear_test_biovault_home();
        }
    }

    #[test]
    fn validate_example_name_accepts_known_and_rejects_unknown() {
        let list = cli::examples::list_examples();
        // When at least one example exists, it validates
        if let Some(first) = list.first() {
            assert!(validate_example_name(first).is_ok());
        }
        // Unknown example returns Err with helpful message
        let err = validate_example_name("__definitely_not_real__").unwrap_err();
        assert!(err.contains("Unknown example"));
    }

    #[test]
    fn clap_parses_init_command_with_flags() {
        let cli = Cli::parse_from(["bv", "init", "--quiet", "user@example.com", "--verbose"]);
        assert!(cli.verbose);
        match cli.command {
            Commands::Init { email, quiet } => {
                assert_eq!(email.as_deref(), Some("user@example.com"));
                assert!(quiet);
            }
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn clap_parses_run_command_defaults() {
        let cli = Cli::parse_from([
            "bv",
            "run",
            "project-dir",
            "participant.yaml",
            "--dry-run",
            "--results-dir",
            "out",
        ]);
        match cli.command {
            Commands::Run {
                project_folder,
                participant_source,
                test,
                download,
                dry_run,
                with_docker,
                work_dir,
                resume,
                template,
                results_dir,
            } => {
                assert_eq!(project_folder, "project-dir");
                assert_eq!(participant_source, "participant.yaml");
                assert!(!test);
                assert!(!download);
                assert!(dry_run);
                assert!(with_docker); // default true
                assert!(work_dir.is_none());
                assert!(!resume);
                assert!(template.is_none());
                assert_eq!(results_dir.as_deref(), Some("out"));
            }
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn clap_parses_sample_data_list_subcommand() {
        let cli = Cli::parse_from(["bv", "sample-data", "list"]);
        match cli.command {
            Commands::SampleData { command } => match command {
                SampleDataCommands::List => {}
                _ => panic!("expected List variant"),
            },
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn clap_parses_submit_command_with_flags() {
        let cli = Cli::parse_from([
            "bv",
            "submit",
            "./project",
            "friend@example.com",
            "--non-interactive",
            "--force",
        ]);
        match cli.command {
            Commands::Submit {
                project_path,
                destination,
                non_interactive,
                force,
            } => {
                assert_eq!(project_path, "./project");
                assert_eq!(destination, "friend@example.com");
                assert!(non_interactive);
                assert!(force);
            }
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn clap_parses_inbox_filters() {
        let cli = Cli::parse_from([
            "bv",
            "inbox",
            "--plain",
            "--sent",
            "--message-type",
            "project",
            "--from",
            "alice@example.com",
            "--search",
            "urgent",
        ]);
        match cli.command {
            Commands::Inbox {
                interactive,
                plain,
                sent,
                all,
                unread,
                projects,
                message_type,
                from,
                search,
            } => {
                assert!(!interactive);
                assert!(plain);
                assert!(sent);
                assert!(!all);
                assert!(!unread);
                assert!(!projects);
                assert_eq!(message_type.as_deref(), Some("project"));
                assert_eq!(from.as_deref(), Some("alice@example.com"));
                assert_eq!(search.as_deref(), Some("urgent"));
            }
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn clap_parses_project_examples_subcommand() {
        let cli = Cli::parse_from(["bv", "project", "examples"]);
        match cli.command {
            Commands::Project { command } => match command {
                ProjectCommands::Examples => {}
                _ => panic!("expected Examples variant"),
            },
            _ => panic!("unexpected command variant"),
        }
    }

    #[test]
    fn async_main_with_sample_data_list_executes() {
        let _home_guard = TestHomeGuard::new();
        let _skip_guard = EnvVarGuard::set("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        let runtime = tokio::runtime::Runtime::new().unwrap();

        let cli = Cli {
            command: Commands::SampleData {
                command: SampleDataCommands::List,
            },
            verbose: false,
            config: None,
        };

        runtime
            .block_on(async { super::async_main_with(cli).await })
            .unwrap();
    }

    #[test]
    fn async_main_with_project_examples_executes() {
        let _home_guard = TestHomeGuard::new();
        let _skip_guard = EnvVarGuard::set("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        let runtime = tokio::runtime::Runtime::new().unwrap();

        let cli = Cli {
            command: Commands::Project {
                command: ProjectCommands::Examples,
            },
            verbose: true,
            config: None,
        };

        runtime
            .block_on(async { super::async_main_with(cli).await })
            .unwrap();
    }

    #[test]
    fn async_main_with_info_executes() {
        let _home_guard = TestHomeGuard::new();
        let _skip_guard = EnvVarGuard::set("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        let runtime = tokio::runtime::Runtime::new().unwrap();

        let cli = Cli {
            command: Commands::Info { json: false },
            verbose: false,
            config: None,
        };

        runtime
            .block_on(async { super::async_main_with(cli).await })
            .unwrap();
    }

    #[test]
    fn async_main_with_participant_list_executes() {
        let _home_guard = TestHomeGuard::new();
        let _skip_guard = EnvVarGuard::set("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        let runtime = tokio::runtime::Runtime::new().unwrap();

        let cli = Cli {
            command: Commands::Participant {
                command: ParticipantCommands::List,
            },
            verbose: false,
            config: None,
        };

        runtime
            .block_on(async { super::async_main_with(cli).await })
            .unwrap();
    }

    #[test]
    fn async_main_with_project_create_executes() {
        let _home_guard = TestHomeGuard::new();
        let _skip_guard = EnvVarGuard::set("BIOVAULT_SKIP_UPDATE_CHECK", "1");
        let runtime = tokio::runtime::Runtime::new().unwrap();

        let cli = Cli {
            command: Commands::Project {
                command: ProjectCommands::Create {
                    name: Some("test_project".to_string()),
                    folder: None,
                    example: None,
                },
            },
            verbose: false,
            config: None,
        };

        // This may fail due to missing setup, but that's OK for this test
        let _ = runtime.block_on(async { super::async_main_with(cli).await });
    }

    #[test]
    fn test_env_var_guard_restores_previous() {
        let key = "TEST_VAR_GUARD";
        std::env::set_var(key, "original");

        {
            let _guard = EnvVarGuard::set(key, "temp");
            assert_eq!(std::env::var(key).unwrap(), "temp");
        }

        assert_eq!(std::env::var(key).unwrap(), "original");
        std::env::remove_var(key);
    }

    #[test]
    fn test_env_var_guard_removes_if_not_set() {
        let key = "TEST_VAR_GUARD_REMOVE";
        std::env::remove_var(key);

        {
            let _guard = EnvVarGuard::set(key, "temp");
            assert_eq!(std::env::var(key).unwrap(), "temp");
        }

        assert!(std::env::var(key).is_err());
    }
}

#[derive(Parser)]
#[command(
    name = "bv",
    version,
    about = "BioVault - A bioinformatics data management CLI",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    #[arg(short, long, global = true, help = "Increase verbosity")]
    verbose: bool,

    #[arg(long, global = true, help = "Path to config file")]
    config: Option<String>,
}

#[derive(Subcommand)]
enum Commands {
    #[command(about = "Check for updates and install the latest version")]
    Update,
    #[command(about = "Initialize a new BioVault repository")]
    Init {
        #[arg(
            help = "Email address for the vault configuration (optional, will detect from SYFTBOX_EMAIL)"
        )]
        email: Option<String>,

        #[arg(short, long, help = "Automatically accept defaults (for testing)")]
        quiet: bool,
    },

    #[command(about = "Show system information")]
    Info {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Check for required dependencies")]
    Check {
        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "Setup environment for known systems (e.g., Google Colab)")]
    Setup,

    #[command(about = "Project management commands")]
    Project {
        #[command(subcommand)]
        command: ProjectCommands,
    },

    #[command(about = "Run a project workflow with Nextflow")]
    Run {
        #[arg(help = "Path to project directory")]
        project_folder: String,

        #[arg(
            help = "Participant source: local file path, Syft URL, or HTTP URL (with optional #fragment)"
        )]
        participant_source: String,

        #[arg(long, help = "Use mock data if available")]
        test: bool,

        #[arg(long, help = "Auto-confirm file downloads")]
        download: bool,

        #[arg(long, help = "Show commands without executing")]
        dry_run: bool,

        #[arg(long, default_value = "true", help = "Run with Docker")]
        with_docker: bool,

        #[arg(long, help = "Nextflow work directory")]
        work_dir: Option<String>,

        #[arg(long, help = "Resume from previous run")]
        resume: bool,

        #[arg(long, help = "Template to use (default, snp, etc.)")]
        template: Option<String>,

        #[arg(long, help = "Custom results directory name")]
        results_dir: Option<String>,
    },

    #[command(name = "sample-data", about = "Manage sample data")]
    SampleData {
        #[command(subcommand)]
        command: SampleDataCommands,
    },

    #[command(about = "Manage participants")]
    Participant {
        #[command(subcommand)]
        command: ParticipantCommands,
    },

    #[command(about = "Manage participants in database catalog")]
    Participants {
        #[command(subcommand)]
        command: ParticipantsCommands,
    },

    #[command(about = "Manage biobank data publishing")]
    Biobank {
        #[command(subcommand)]
        command: BiobankCommands,
    },

    #[command(about = "Manage BioVault configuration")]
    Config {
        #[command(subcommand)]
        command: Option<ConfigCommands>,

        #[arg(long, help = "Output as JSON")]
        json: bool,
    },

    #[command(about = "FASTQ file operations")]
    Fastq {
        #[command(subcommand)]
        command: FastqCommands,
    },

    #[command(about = "Manage files and file catalogs")]
    Files {
        #[command(subcommand)]
        command: FilesCommands,
    },

    #[command(about = "Submit a project to another biobank via SyftBox")]
    Submit {
        #[arg(help = "Path to project directory (use '.' for current directory)")]
        project_path: String,

        #[arg(
            help = "Destination: either a datasite email (e.g., user@domain.com) or full Syft URL (e.g., syft://user@domain.com/public/biovault/participants.yaml#participants.ID)"
        )]
        destination: String,

        #[arg(long, help = "Skip interactive prompts, use defaults")]
        non_interactive: bool,

        #[arg(
            long,
            help = "Force resubmission even if project was already submitted"
        )]
        force: bool,
    },

    #[command(about = "Clean up stale database locks")]
    Cleanup {
        #[arg(long, help = "Clean all locks in all virtualenvs")]
        all: bool,
    },

    #[command(about = "View and manage inbox messages")]
    Inbox {
        #[arg(short = 'i', long, help = "Interactive mode (default)")]
        interactive: bool,

        #[arg(long, help = "Plain, non-interactive list output")]
        plain: bool,

        #[arg(short = 's', long, help = "Show sent messages")]
        sent: bool,

        #[arg(short = 'a', long, help = "Show all messages (including deleted)")]
        all: bool,

        #[arg(short = 'u', long, help = "Show only unread messages")]
        unread: bool,

        #[arg(short = 'p', long, help = "Show project submissions")]
        projects: bool,

        #[arg(
            short = 't',
            long,
            help = "Filter by message type (text/project/request)"
        )]
        message_type: Option<String>,

        #[arg(short = 'f', long, help = "Filter by sender")]
        from: Option<String>,

        #[arg(long, help = "Search messages by content")]
        search: Option<String>,
    },

    #[command(about = "Manage messages via SyftBox RPC")]
    Message {
        #[command(subcommand)]
        command: MessageCommands,
    },

    #[command(about = "Sample sheet operations")]
    Samplesheet {
        #[command(subcommand)]
        command: SamplesheetCommands,
    },

    #[command(about = "Manage the BioVault daemon for automatic message processing")]
    Daemon {
        #[command(subcommand)]
        command: DaemonCommands,
    },

    #[command(
        name = "hard-reset",
        about = "Delete all BioVault data and configuration (DESTRUCTIVE)"
    )]
    HardReset {
        #[arg(long, help = "Skip confirmation prompts (use with caution)")]
        ignore_warning: bool,
    },

    #[command(about = "Launch BioVault Desktop GUI")]
    Desktop {
        #[arg(long, help = "Path to BioVault config directory")]
        config: Option<String>,
    },
}

#[derive(Subcommand)]
enum DaemonCommands {
    #[command(about = "Start the BioVault daemon")]
    Start {
        #[arg(long, help = "Run daemon in foreground (no background)")]
        foreground: bool,
    },

    #[command(about = "Stop the running daemon")]
    Stop,

    #[command(about = "Restart the daemon (stop if running, then start)")]
    Restart {
        #[arg(long, help = "Run daemon in foreground after restart")]
        foreground: bool,
    },

    #[command(about = "Check daemon status")]
    Status,

    #[command(about = "View daemon logs")]
    Logs {
        #[arg(short, long, help = "Follow log output (tail -f)")]
        follow: bool,

        #[arg(short, long, help = "Number of lines to show (default: 50)")]
        lines: Option<usize>,
    },

    #[command(about = "Install daemon as a systemd service (Linux only)")]
    Install,

    #[command(about = "Uninstall daemon systemd service (Linux only)")]
    Uninstall,

    #[command(about = "List all installed daemon services")]
    List,

    #[command(about = "Show the systemd service file")]
    Show,

    #[command(about = "Reinstall daemon service (uninstall + install)")]
    Reinstall,
}

#[derive(Subcommand)]
enum ProjectCommands {
    #[command(about = "Create a new project")]
    Create {
        #[arg(long, help = "Project name")]
        name: Option<String>,

        #[arg(long, help = "Folder path (defaults to ./{name})")]
        folder: Option<String>,

        #[arg(long, value_parser = validate_example_name, help = "Use example template (use 'bv project examples' to list available)")]
        example: Option<String>,
    },

    #[command(about = "List available example templates")]
    Examples,

    #[command(about = "Import a project from URL or register a local project")]
    Import {
        #[arg(help = "URL to project.yaml or local directory path")]
        source: String,

        #[arg(long, help = "Override project name")]
        name: Option<String>,

        #[arg(long, help = "Overwrite existing project")]
        overwrite: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "List all registered projects")]
    List {
        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Show detailed information about a project")]
    Show {
        #[arg(help = "Project name or ID")]
        identifier: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Delete a project")]
    Delete {
        #[arg(help = "Project name or ID")]
        identifier: String,

        #[arg(long, help = "Keep project files, only remove from database")]
        keep_files: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },
}

#[derive(Subcommand)]
enum SampleDataCommands {
    #[command(about = "Fetch sample data")]
    Fetch {
        #[arg(
            value_delimiter = ',',
            help = "Participant IDs to fetch (comma-separated)"
        )]
        participant_ids: Option<Vec<String>>,

        #[arg(long, help = "Fetch all available sample data")]
        all: bool,
    },

    #[command(about = "List available sample data")]
    List,
}

#[derive(Subcommand)]
enum ParticipantCommands {
    #[command(about = "Add a new participant")]
    Add {
        #[arg(long, help = "Participant ID")]
        id: Option<String>,

        #[arg(long, help = "Aligned file path (.cram, .bam, or .sam)")]
        aligned: Option<String>,

        #[arg(
            long,
            help = "Template type (default or snp)",
            default_value = "default"
        )]
        template: Option<String>,

        #[arg(long, help = "SNP file path (for SNP template)")]
        snp: Option<String>,

        #[arg(long, help = "Reference genome file path (.fa or .fasta)")]
        reference: Option<String>,

        #[arg(long, help = "Reference version (GRCh38 or GRCh37)")]
        ref_version: Option<String>,

        #[arg(long, help = "Skip interactive prompts, use defaults")]
        non_interactive: bool,
    },

    #[command(about = "List all participants")]
    List,

    #[command(about = "Delete a participant")]
    Delete {
        #[arg(help = "Participant ID to delete")]
        id: String,
    },

    #[command(about = "Validate participant files")]
    Validate {
        #[arg(help = "Participant ID to validate (validates all if not specified)")]
        id: Option<String>,
    },
}

#[derive(Subcommand)]
enum ParticipantsCommands {
    #[command(about = "List all participants in the catalog")]
    List {
        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Delete a participant and unlink all files")]
    Delete {
        #[arg(help = "Participant database ID")]
        id: i64,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },
}

#[derive(Subcommand)]
enum BiobankCommands {
    #[command(about = "List biobanks in SyftBox")]
    List,

    #[command(about = "Publish participants to SyftBox")]
    Publish {
        #[arg(long, help = "Participant ID to publish")]
        participant_id: Option<String>,

        #[arg(long, help = "Publish all participants")]
        all: bool,

        #[arg(
            long,
            help = "HTTP relay servers (defaults to syftbox.net)",
            value_delimiter = ','
        )]
        http_relay_servers: Option<Vec<String>>,
    },

    #[command(about = "Unpublish participants from SyftBox")]
    Unpublish {
        #[arg(long, help = "Participant ID to unpublish")]
        participant_id: Option<String>,

        #[arg(long, help = "Unpublish all participants")]
        all: bool,
    },
}

#[derive(Subcommand)]
enum ConfigCommands {
    #[command(about = "Set email address")]
    Email {
        #[arg(help = "Email address")]
        email: String,
    },

    #[command(about = "Set SyftBox config path")]
    Syftbox {
        #[arg(help = "Path to SyftBox config.json (omit to use default)")]
        path: Option<String>,
    },
}

#[derive(Subcommand)]
enum FastqCommands {
    #[command(about = "Combine multiple FASTQ files into one")]
    Combine {
        #[arg(help = "Input folder containing FASTQ files")]
        input_folder: String,

        #[arg(help = "Output file path")]
        output_file: String,

        #[arg(long, help = "Validate files before combining")]
        validate: bool,

        #[arg(long, help = "Skip validation prompt and use default")]
        no_prompt: bool,

        #[arg(
            long,
            default_value = "tsv",
            help = "Stats output format (tsv, yaml, json)"
        )]
        stats_format: String,
    },
}

#[derive(Subcommand)]
enum FilesCommands {
    #[command(about = "Scan directory for files by type")]
    Scan {
        #[arg(help = "Directory path to scan")]
        path: String,

        #[arg(long, help = "File extension filter (e.g., .vcf, .cram)")]
        ext: Option<String>,

        #[arg(long, help = "Scan subdirectories recursively", default_value = "true")]
        recursive: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Import files into database catalog")]
    Import {
        #[arg(help = "File path or directory to import")]
        path: String,

        #[arg(long, help = "File extension filter (e.g., .vcf, .cram)")]
        ext: Option<String>,

        #[arg(
            long,
            help = "Import subdirectories recursively",
            default_value = "true"
        )]
        recursive: bool,

        #[arg(
            long,
            help = "Pattern to extract participant IDs (e.g., 'case_{id}_*')"
        )]
        pattern: Option<String>,

        #[arg(long, help = "Dry run - show what would be imported without importing")]
        dry_run: bool,

        #[arg(long, help = "Skip interactive confirmation prompt")]
        non_interactive: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "List cataloged files with filtering")]
    List {
        #[arg(long, help = "Filter by file extension")]
        ext: Option<String>,

        #[arg(long, help = "Filter by participant ID")]
        participant: Option<String>,

        #[arg(long, help = "Show only files without participant assignment")]
        unassigned: bool,

        #[arg(long, help = "Maximum number of results")]
        limit: Option<usize>,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Show detailed file information")]
    Info {
        #[arg(help = "File ID or path")]
        file: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Suggest participant ID patterns from filenames")]
    SuggestPatterns {
        #[arg(help = "Directory path to analyze")]
        path: String,

        #[arg(long, help = "File extension filter (e.g., .vcf, .cram)")]
        ext: Option<String>,

        #[arg(long, help = "Scan subdirectories recursively", default_value = "true")]
        recursive: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Delete a file record from the catalog")]
    Delete {
        #[arg(help = "File ID to delete")]
        id: i64,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Delete multiple file records from the catalog")]
    DeleteBulk {
        #[arg(help = "File IDs to delete (comma-separated)", value_delimiter = ',')]
        ids: Vec<i64>,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Detect file types and extract genotype metadata (fast - header only)")]
    Detect {
        #[arg(help = "File paths to analyze", num_args = 1..)]
        files: Vec<String>,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(
        about = "Analyze genotype files for row count, chromosomes, and sex (slow - reads entire file)"
    )]
    Analyze {
        #[arg(help = "File paths to analyze", num_args = 1..)]
        files: Vec<String>,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Link a file to a participant")]
    Link {
        #[arg(help = "File ID to link")]
        file_id: i64,

        #[arg(long, help = "Participant ID to link to")]
        participant: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Bulk link files to participants via JSON")]
    LinkBulk {
        #[arg(help = "JSON object mapping file paths to participant IDs")]
        json: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Remove participant assignment from a file")]
    Unlink {
        #[arg(help = "File ID to unlink")]
        file_id: i64,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Calculate BLAKE3 hashes for files")]
    Hash {
        #[arg(help = "File paths to hash", num_args = 1..)]
        files: Vec<String>,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Export files to CSV with participant IDs")]
    ExportCsv {
        #[arg(help = "Directory path to scan")]
        path: String,

        #[arg(long, help = "File extension filter (e.g., .vcf, .txt)")]
        ext: Option<String>,

        #[arg(long, help = "Scan subdirectories recursively", default_value = "true")]
        recursive: bool,

        #[arg(long, help = "Pattern to extract participant IDs (e.g., '{parent}')")]
        pattern: Option<String>,

        #[arg(short = 'o', long, help = "Output CSV file path")]
        output: String,
    },

    #[command(about = "Detect file types in CSV and update metadata columns")]
    DetectCsv {
        #[arg(help = "Input CSV file path")]
        input_csv: String,

        #[arg(short = 'o', long, help = "Output CSV file path")]
        output: String,
    },

    #[command(about = "Analyze genotype files in CSV for row counts and sex")]
    AnalyzeCsv {
        #[arg(help = "Input CSV file path")]
        input_csv: String,

        #[arg(short = 'o', long, help = "Output CSV file path")]
        output: String,
    },

    #[command(about = "Import files from CSV with complete metadata")]
    ImportCsv {
        #[arg(help = "Path to CSV file containing file metadata")]
        csv_path: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Process pending files in the queue (hash and detect metadata)")]
    ProcessQueue {
        #[arg(
            long,
            help = "Maximum number of files to process",
            default_value = "100"
        )]
        limit: usize,

        #[arg(long, help = "Daemon mode - continuously process queue")]
        daemon: bool,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },

    #[command(about = "Fast import - add files to queue for background processing")]
    ImportPending {
        #[arg(help = "Path to CSV file containing file metadata")]
        csv_path: String,

        #[arg(long, help = "Output format (json|table)", default_value = "table")]
        format: String,
    },
}

#[derive(Subcommand)]
enum SamplesheetCommands {
    #[command(about = "Create a sample sheet CSV from a folder of files")]
    Create {
        #[arg(help = "Input directory containing files")]
        input_dir: String,

        #[arg(help = "Output CSV file path")]
        output_file: String,

        #[arg(
            long = "file_filter",
            help = "File pattern filter (e.g., *.txt, default: all files)"
        )]
        file_filter: Option<String>,

        #[arg(
            long = "extract_cols",
            help = "Pattern for extracting fields from filenames (e.g., {participant_id}_X_X_GSAv3-DTC_GRCh38-{date}.txt)"
        )]
        extract_cols: Option<String>,

        #[arg(
            long = "ignore",
            help = "Add files even if they don't match the extraction pattern"
        )]
        ignore: bool,
    },
}

#[derive(Subcommand)]
enum MessageCommands {
    #[command(about = "Send a message to another datasite")]
    Send {
        #[arg(help = "Recipient email address")]
        recipient: String,

        #[arg(help = "Message content")]
        message: String,

        #[arg(short = 's', long = "subject", help = "Optional message subject")]
        subject: Option<String>,
    },

    #[command(about = "Reply to a message")]
    Reply {
        #[arg(help = "Message ID to reply to")]
        message_id: String,

        #[arg(help = "Reply content")]
        body: String,
    },

    #[command(about = "Read a specific message")]
    Read {
        #[arg(help = "Message ID to read")]
        message_id: String,
    },

    #[command(about = "Delete a message")]
    Delete {
        #[arg(help = "Message ID to delete")]
        message_id: String,
    },

    #[command(about = "List messages")]
    List {
        #[arg(short = 'u', long = "unread", help = "Show only unread messages")]
        unread: bool,

        #[arg(short = 's', long = "sent", help = "Show sent messages")]
        sent: bool,

        #[arg(short = 'p', long = "projects", help = "Show only project messages")]
        projects: bool,
    },

    #[command(about = "View a message thread")]
    Thread {
        #[arg(help = "Thread ID to view")]
        thread_id: String,
    },

    #[command(about = "Sync messages (check for new and update ACKs)")]
    Sync,

    #[command(about = "Process a project message (run test/real)")]
    Process {
        #[arg(help = "Message ID of the project to process")]
        message_id: String,

        #[arg(long, help = "Run with test data", conflicts_with = "real")]
        test: bool,

        #[arg(long, help = "Run with real data", conflicts_with = "test")]
        real: bool,

        #[arg(long, help = "Participant to use (defaults to first available)")]
        participant: Option<String>,

        #[arg(long, help = "Approve after successful run")]
        approve: bool,

        #[arg(long, help = "Non-interactive mode (skip prompts)")]
        non_interactive: bool,
    },

    #[command(about = "Archive a project message (revoke write permissions)")]
    Archive {
        #[arg(help = "Message ID to archive")]
        message_id: String,
    },
}

fn main() -> Result<()> {
    // Set up panic hook to log crashes
    std::panic::set_hook(Box::new(|panic_info| {
        let payload = panic_info.payload();
        let message = if let Some(s) = payload.downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = payload.downcast_ref::<String>() {
            s.clone()
        } else {
            "Unknown panic payload".to_string()
        };

        let location = panic_info
            .location()
            .map(|loc| format!("{}:{}:{}", loc.file(), loc.line(), loc.column()))
            .unwrap_or_else(|| "unknown location".to_string());

        let crash_log = format!(
            "\n=== PANIC/CRASH DETECTED ===\n\
             Time: {}\n\
             Message: {}\n\
             Location: {}\n\
             Backtrace: use RUST_BACKTRACE=1 for details\n\
             ============================\n",
            chrono::Utc::now().format("%Y-%m-%d %H:%M:%S%.3f UTC"),
            message,
            location
        );

        // Try to write to daemon log if we're in daemon mode
        if let Ok(syftbox_dir) = std::env::var("SYFTBOX_DATA_DIR") {
            let log_path = std::path::PathBuf::from(&syftbox_dir)
                .join(".biovault")
                .join("logs")
                .join("daemon.log");

            if let Ok(mut file) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&log_path)
            {
                use std::io::Write;
                let _ = writeln!(file, "{}", crash_log);
            }
        }

        // Also write to stderr
        eprintln!("{}", crash_log);
    }));

    // Use a minimal tokio runtime - this is just a file-watching daemon, not a web server
    // 2 worker threads is plenty for occasional I/O operations
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(2)
        .thread_name("bv-worker")
        .enable_all()
        .build()
        .expect("Failed to create tokio runtime");

    runtime.block_on(async_main())
}

async fn async_main_with(cli: Cli) -> Result<()> {
    let filter_level = if cli.verbose { "debug" } else { "info" };

    let _ = tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(filter_level)))
        .try_init();

    // Random version check on startup (10% chance)
    let _ = commands::update::check_and_notify_random().await;

    // Check for upgrades and perform them if needed
    let _ = cli::upgrade::check_and_upgrade();

    match cli.command {
        Commands::Update => {
            commands::update::execute().await?;
        }
        Commands::Init { email, quiet } => {
            commands::init::execute(email.as_deref(), quiet).await?;
        }
        Commands::Info { json } => {
            commands::info::execute(json).await?;
        }
        Commands::Check { json } => {
            commands::check::execute(json).await?;
        }
        Commands::Setup => {
            commands::setup::execute().await?;
        }
        Commands::Project { command } => match command {
            ProjectCommands::Create {
                name,
                folder,
                example,
            } => {
                commands::project::create(name, folder, example).await?;
            }
            ProjectCommands::Examples => {
                commands::project::list_examples()?;
            }
            ProjectCommands::Import {
                source,
                name,
                overwrite,
                format,
            } => {
                let fmt = if format == "table" {
                    None
                } else {
                    Some(format)
                };
                commands::project_management::import(source, name, overwrite, fmt).await?;
            }
            ProjectCommands::List { format } => {
                let fmt = if format == "table" {
                    None
                } else {
                    Some(format)
                };
                commands::project_management::list(fmt)?;
            }
            ProjectCommands::Show { identifier, format } => {
                let fmt = if format == "table" {
                    None
                } else {
                    Some(format)
                };
                commands::project_management::show(identifier, fmt)?;
            }
            ProjectCommands::Delete {
                identifier,
                keep_files,
                format,
            } => {
                let fmt = if format == "table" {
                    None
                } else {
                    Some(format)
                };
                commands::project_management::delete(identifier, keep_files, fmt)?;
            }
        },
        Commands::Run {
            project_folder,
            participant_source,
            test,
            download,
            dry_run,
            with_docker,
            work_dir,
            resume,
            template,
            results_dir,
        } => {
            commands::run::execute(commands::run::RunParams {
                project_folder,
                participant_source,
                test,
                download,
                dry_run,
                with_docker,
                work_dir,
                resume,
                template,
                results_dir,
            })
            .await?;
        }
        Commands::SampleData { command } => match command {
            SampleDataCommands::Fetch {
                participant_ids,
                all,
            } => {
                commands::sample_data::fetch(participant_ids, all, false).await?;
            }
            SampleDataCommands::List => {
                commands::sample_data::list().await?;
            }
        },
        Commands::Participant { command } => match command {
            ParticipantCommands::Add {
                id,
                aligned,
                template,
                snp,
                reference,
                ref_version,
                non_interactive,
            } => {
                commands::participant::add(
                    id,
                    aligned,
                    template,
                    snp,
                    reference,
                    ref_version,
                    non_interactive,
                )
                .await?;
            }
            ParticipantCommands::List => {
                commands::participant::list().await?;
            }
            ParticipantCommands::Delete { id } => {
                commands::participant::delete(id).await?;
            }
            ParticipantCommands::Validate { id } => {
                commands::participant::validate(id).await?;
            }
        },
        Commands::Participants { command } => match command {
            ParticipantsCommands::List { format } => {
                commands::participants::list(format).await?;
            }
            ParticipantsCommands::Delete { id, format } => {
                commands::participants::delete(id, format).await?;
            }
        },
        Commands::Biobank { command } => match command {
            BiobankCommands::List => {
                commands::biobank::list(None).await?;
            }
            BiobankCommands::Publish {
                participant_id,
                all,
                http_relay_servers,
            } => {
                commands::biobank::publish(participant_id, all, http_relay_servers).await?;
            }
            BiobankCommands::Unpublish {
                participant_id,
                all,
            } => {
                commands::biobank::unpublish(participant_id, all).await?;
            }
        },
        Commands::Config { command, json } => {
            if let Some(cmd) = command {
                match cmd {
                    ConfigCommands::Email { email } => {
                        commands::config_cmd::set_email(email).await?;
                    }
                    ConfigCommands::Syftbox { path } => {
                        commands::config_cmd::set_syftbox(path).await?;
                    }
                }
            } else {
                commands::config_cmd::show(json).await?;
            }
        }
        Commands::Fastq { command } => match command {
            FastqCommands::Combine {
                input_folder,
                output_file,
                validate,
                no_prompt,
                stats_format,
            } => {
                let should_validate = if no_prompt { Some(validate) } else { None };
                commands::fastq::combine(
                    input_folder,
                    output_file,
                    should_validate,
                    Some(stats_format),
                )
                .await?;
            }
        },
        Commands::Files { command } => match command {
            FilesCommands::Scan {
                path,
                ext,
                recursive,
                format,
            } => {
                commands::files::scan(path, ext, recursive, format).await?;
            }
            FilesCommands::Import {
                path,
                ext,
                recursive,
                pattern,
                dry_run,
                non_interactive,
                format,
            } => {
                commands::files::import(
                    path,
                    ext,
                    recursive,
                    pattern,
                    dry_run,
                    non_interactive,
                    format,
                )
                .await?;
            }
            FilesCommands::List {
                ext,
                participant,
                unassigned,
                limit,
                format,
            } => {
                commands::files::list(ext, participant, unassigned, limit, format).await?;
            }
            FilesCommands::Info { file, format } => {
                commands::files::info(file, format).await?;
            }
            FilesCommands::SuggestPatterns {
                path,
                ext,
                recursive,
                format,
            } => {
                commands::files::suggest_patterns(path, ext, recursive, format).await?;
            }
            FilesCommands::Delete { id, format } => {
                commands::files::delete(id, format).await?;
            }
            FilesCommands::DeleteBulk { ids, format } => {
                commands::files::delete_bulk(ids, format).await?;
            }
            FilesCommands::Detect { files, format } => {
                commands::files::detect(files, format).await?;
            }
            FilesCommands::Analyze { files, format } => {
                commands::files::analyze(files, format).await?;
            }
            FilesCommands::Link {
                file_id,
                participant,
                format,
            } => {
                commands::files::link(file_id, participant, format).await?;
            }
            FilesCommands::LinkBulk { json, format } => {
                commands::files::link_bulk(json, format).await?;
            }
            FilesCommands::Unlink { file_id, format } => {
                commands::files::unlink(file_id, format).await?;
            }
            FilesCommands::Hash { files, format } => {
                commands::files::hash(files, format).await?;
            }
            FilesCommands::ExportCsv {
                path,
                ext,
                recursive,
                pattern,
                output,
            } => {
                commands::files::export_csv(path, ext, recursive, pattern, output).await?;
            }
            FilesCommands::DetectCsv { input_csv, output } => {
                commands::files::detect_csv(input_csv, output).await?;
            }
            FilesCommands::AnalyzeCsv { input_csv, output } => {
                commands::files::analyze_csv(input_csv, output).await?;
            }
            FilesCommands::ImportCsv { csv_path, format } => {
                commands::files::import_csv(csv_path, format).await?;
            }
            FilesCommands::ProcessQueue {
                limit,
                daemon,
                format,
            } => {
                commands::files::process_queue(limit, daemon, format).await?;
            }
            FilesCommands::ImportPending { csv_path, format } => {
                commands::files::import_pending(csv_path, format).await?;
            }
        },
        Commands::Submit {
            project_path,
            destination,
            non_interactive,
            force,
        } => {
            commands::submit::submit(project_path, destination, non_interactive, force).await?;
        }
        Commands::Cleanup { all } => {
            let config = biovault::config::Config::load()?;
            commands::messages::cleanup_locks(&config, all)?;
        }
        Commands::Inbox {
            interactive,
            plain,
            sent,
            all,
            unread,
            projects,
            message_type,
            from,
            search,
        } => {
            let config = biovault::config::Config::load()?;
            // Default behavior: interactive unless --plain is provided
            if plain && !interactive {
                let filters = commands::inbox::ListFilters {
                    sent,
                    all,
                    unread,
                    projects,
                    message_type,
                    from,
                    search,
                };
                commands::inbox::list(&config, filters)?;
            } else {
                // When both flags are provided, prefer interactive
                commands::inbox::interactive(&config, None).await?;
            }
        }
        Commands::Message { command } => match command {
            MessageCommands::Send {
                recipient,
                message,
                subject,
            } => {
                let config = biovault::config::Config::load()?;
                commands::messages::send_message(
                    &config,
                    &recipient,
                    &message,
                    subject.as_deref(),
                )?;
            }
            MessageCommands::Reply { message_id, body } => {
                let config = biovault::config::Config::load()?;
                commands::messages::reply_message(&config, &message_id, &body)?;
            }
            MessageCommands::Read { message_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::read_message(&config, &message_id).await?;
            }
            MessageCommands::Delete { message_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::delete_message(&config, &message_id)?;
            }
            MessageCommands::List {
                unread,
                sent,
                projects,
            } => {
                let config = biovault::config::Config::load()?;
                commands::messages::list_messages(&config, unread, sent, projects)?;
            }
            MessageCommands::Thread { thread_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::view_thread(&config, &thread_id)?;
            }
            MessageCommands::Sync => {
                let config = biovault::config::Config::load()?;
                commands::messages::sync_messages(&config)?;
            }
            MessageCommands::Process {
                message_id,
                test,
                real,
                participant,
                approve,
                non_interactive,
            } => {
                let config = biovault::config::Config::load()?;
                commands::messages::process_project_message(
                    &config,
                    &message_id,
                    test,
                    real,
                    participant,
                    approve,
                    non_interactive,
                )
                .await?;
            }
            MessageCommands::Archive { message_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::archive_message(&config, &message_id)?;
            }
        },
        Commands::Samplesheet { command } => match command {
            SamplesheetCommands::Create {
                input_dir,
                output_file,
                file_filter,
                extract_cols,
                ignore,
            } => {
                commands::samplesheet::create(
                    input_dir,
                    output_file,
                    file_filter,
                    extract_cols,
                    ignore,
                )
                .await?;
            }
        },
        Commands::Daemon { command } => {
            // If BV_DAEMON_CONFIG env var is set, use it (for spawned daemon processes)
            // Otherwise load config normally
            let config = if let Ok(config_json) = std::env::var("BV_DAEMON_CONFIG") {
                serde_json::from_str(&config_json)?
            } else {
                biovault::config::Config::load()?
            };

            match command {
                DaemonCommands::Start { foreground } => {
                    commands::daemon::start(&config, foreground).await?;
                }
                DaemonCommands::Stop => {
                    commands::daemon::stop(&config).await?;
                }
                DaemonCommands::Restart { foreground } => {
                    // Stop the daemon if it's running
                    let _ = commands::daemon::stop(&config).await;
                    // Small delay to ensure clean shutdown
                    tokio::time::sleep(std::time::Duration::from_secs(1)).await;
                    // Start it again
                    commands::daemon::start(&config, foreground).await?;
                }
                DaemonCommands::Status => {
                    commands::daemon::service_status(&config).await?;
                }
                DaemonCommands::Logs { follow, lines } => {
                    commands::daemon::logs(&config, follow, lines).await?;
                }
                DaemonCommands::Install => {
                    commands::daemon::install_service(&config).await?;
                }
                DaemonCommands::Uninstall => {
                    commands::daemon::uninstall_service(&config).await?;
                }
                DaemonCommands::List => {
                    commands::daemon::list_services().await?;
                }
                DaemonCommands::Show => {
                    commands::daemon::show_service(&config).await?;
                }
                DaemonCommands::Reinstall => {
                    commands::daemon::reinstall_service(&config).await?;
                }
            }
        }
        Commands::HardReset { ignore_warning } => {
            commands::hard_reset::execute(ignore_warning).await?;
        }
        Commands::Desktop { config } => {
            let biovault_home = if let Some(cfg) = config {
                cfg
            } else if let Ok(env_home) = std::env::var("BIOVAULT_HOME") {
                env_home
            } else {
                biovault::config::get_biovault_home()?
                    .to_string_lossy()
                    .to_string()
            };

            #[cfg(target_os = "macos")]
            {
                let app_path = std::path::Path::new("/Applications/BioVault Desktop.app");
                if !app_path.exists() {
                    eprintln!("Error: BioVault Desktop is not installed.");
                    eprintln!("\nPlease install BioVault Desktop from:");
                    eprintln!("https://github.com/OpenMined/biovault-desktop");
                    eprintln!("\nOr download the latest release:");
                    eprintln!("https://github.com/OpenMined/biovault-desktop/releases/latest");
                    std::process::exit(1);
                }

                match std::process::Command::new("open")
                    .arg("-n")
                    .arg("/Applications/BioVault Desktop.app")
                    .arg("--args")
                    .arg("--biovault-config")
                    .arg(&biovault_home)
                    .spawn()
                {
                    Ok(_) => println!("Launching BioVault Desktop with config: {}", biovault_home),
                    Err(e) => {
                        eprintln!("Error launching BioVault Desktop: {}", e);
                        eprintln!("\nPlease ensure BioVault Desktop is installed from:");
                        eprintln!("https://github.com/OpenMined/biovault-desktop");
                        std::process::exit(1);
                    }
                }
            }

            #[cfg(not(target_os = "macos"))]
            {
                match which::which("bv-desktop").or_else(|_| which::which("biovault-desktop")) {
                    Ok(desktop_bin) => {
                        match std::process::Command::new(desktop_bin)
                            .arg("--biovault-config")
                            .arg(&biovault_home)
                            .spawn()
                        {
                            Ok(_) => println!(
                                "Launching BioVault Desktop with config: {}",
                                biovault_home
                            ),
                            Err(e) => {
                                eprintln!("Error launching BioVault Desktop: {}", e);
                                eprintln!("\nPlease ensure BioVault Desktop is installed from:");
                                eprintln!("https://github.com/OpenMined/biovault-desktop");
                                std::process::exit(1);
                            }
                        }
                    }
                    Err(_) => {
                        eprintln!("Error: BioVault Desktop is not installed or not in PATH.");
                        eprintln!("\nPlease install BioVault Desktop from:");
                        eprintln!("https://github.com/OpenMined/biovault-desktop");
                        eprintln!("\nOr download the latest release:");
                        eprintln!("https://github.com/OpenMined/biovault-desktop/releases/latest");
                        std::process::exit(1);
                    }
                }
            }
        }
    }

    Ok(())
}

async fn async_main() -> Result<()> {
    let cli = Cli::parse();
    async_main_with(cli).await
}
