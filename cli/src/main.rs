use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

use biovault::cli;

use cli::commands;

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
    },

    #[command(about = "Show system information")]
    Info,

    #[command(about = "Check for required dependencies")]
    Check,

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

    #[command(about = "Manage biobank data publishing")]
    Biobank {
        #[command(subcommand)]
        command: BiobankCommands,
    },

    #[command(about = "Manage BioVault configuration")]
    Config {
        #[command(subcommand)]
        command: Option<ConfigCommands>,
    },

    #[command(about = "FASTQ file operations")]
    Fastq {
        #[command(subcommand)]
        command: FastqCommands,
    },

    #[command(about = "Submit a project to another biobank via SyftBox")]
    Submit {
        #[arg(help = "Path to project directory (use '.' for current directory)")]
        project_path: String,

        #[arg(
            help = "Destination: either a datasite email (e.g., user@domain.com) or full Syft URL (e.g., syft://user@domain.com/public/biovault/participants.yaml#participants.ID)"
        )]
        destination: String,
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
}

#[derive(Subcommand)]
enum ProjectCommands {
    #[command(about = "Create a new project")]
    Create {
        #[arg(long, help = "Project name")]
        name: Option<String>,

        #[arg(long, help = "Folder path (defaults to ./{name})")]
        folder: Option<String>,
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
    },

    #[command(about = "View a message thread")]
    Thread {
        #[arg(help = "Thread ID to view")]
        thread_id: String,
    },

    #[command(about = "Sync messages (check for new and update ACKs)")]
    Sync,
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    let filter_level = if cli.verbose { "debug" } else { "info" };

    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(filter_level)))
        .init();

    // Random version check on startup (10% chance)
    let _ = commands::update::check_and_notify_random().await;

    match cli.command {
        Commands::Update => {
            commands::update::execute().await?;
        }
        Commands::Init { email } => {
            commands::init::execute(email.as_deref()).await?;
        }
        Commands::Info => {
            commands::info::execute().await?;
        }
        Commands::Check => {
            commands::check::execute().await?;
        }
        Commands::Setup => {
            commands::setup::execute().await?;
        }
        Commands::Project { command } => match command {
            ProjectCommands::Create { name, folder } => {
                commands::project::create(name, folder).await?;
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
            })
            .await?;
        }
        Commands::SampleData { command } => match command {
            SampleDataCommands::Fetch {
                participant_ids,
                all,
            } => {
                commands::sample_data::fetch(participant_ids, all).await?;
            }
            SampleDataCommands::List => {
                commands::sample_data::list().await?;
            }
        },
        Commands::Participant { command } => match command {
            ParticipantCommands::Add { id, aligned } => {
                commands::participant::add(id, aligned).await?;
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
        Commands::Config { command } => {
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
                commands::config_cmd::show().await?;
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
        Commands::Submit {
            project_path,
            destination,
        } => {
            commands::submit::submit(project_path, destination).await?;
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
                commands::inbox::interactive(&config, None)?;
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
                commands::messages::read_message(&config, &message_id)?;
            }
            MessageCommands::Delete { message_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::delete_message(&config, &message_id)?;
            }
            MessageCommands::List { unread } => {
                let config = biovault::config::Config::load()?;
                commands::messages::list_messages(&config, unread)?;
            }
            MessageCommands::Thread { thread_id } => {
                let config = biovault::config::Config::load()?;
                commands::messages::view_thread(&config, &thread_id)?;
            }
            MessageCommands::Sync => {
                let config = biovault::config::Config::load()?;
                commands::messages::sync_messages(&config)?;
            }
        },
    }

    Ok(())
}
