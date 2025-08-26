use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

mod cli;
mod config;

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
    }

    Ok(())
}
