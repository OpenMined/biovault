use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

mod cli;
mod config;
mod error;

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
    #[command(about = "Initialize a new BioVault repository")]
    Init {
        #[arg(help = "Email address for the vault configuration")]
        email: String,
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

        #[arg(help = "Path to patient file (YAML)")]
        patient_file: String,

        #[arg(
            long,
            value_delimiter = ',',
            help = "Comma-separated list of patient IDs"
        )]
        patients: Option<Vec<String>>,

        #[arg(long, help = "Process single patient")]
        patient: Option<String>,

        #[arg(long, help = "Process all patients in file")]
        all: bool,

        #[arg(long, help = "Run TEST patient only")]
        test: bool,

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
        #[arg(value_delimiter = ',', help = "Patient IDs to fetch (comma-separated)")]
        patient_ids: Option<Vec<String>>,

        #[arg(long, help = "Fetch all available sample data")]
        all: bool,
    },

    #[command(about = "List available sample data")]
    List,
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    let filter_level = if cli.verbose { "debug" } else { "info" };

    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(filter_level)))
        .init();

    match cli.command {
        Commands::Init { email } => {
            info!("Initializing BioVault with email: {}", email);
            commands::init::execute(&email).await?;
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
            patient_file,
            patients,
            patient,
            all,
            test,
            dry_run,
            with_docker,
            work_dir,
            resume,
        } => {
            commands::run::execute(commands::run::RunParams {
                project_folder,
                patient_file,
                patients,
                patient,
                all,
                test,
                dry_run,
                with_docker,
                work_dir,
                resume,
            })
            .await?;
        }
        Commands::SampleData { command } => match command {
            SampleDataCommands::Fetch { patient_ids, all } => {
                commands::sample_data::fetch(patient_ids, all).await?;
            }
            SampleDataCommands::List => {
                commands::sample_data::list().await?;
            }
        },
    }

    Ok(())
}
