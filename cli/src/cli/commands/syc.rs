use crate::config::{get_biovault_home, Config};
use crate::syftbox::syc as syc_utils;
use crate::Result;
use anyhow::Context;
use clap::Subcommand;
use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub enum SycCommands {
    /// Import another participant's public crypto bundle
    Import {
        #[arg(help = "Path to the bundle JSON (e.g., datasites/<email>/public/crypto/did.json)")]
        bundle: PathBuf,

        #[arg(long, help = "Expected identity email inside the bundle")]
        expected_identity: Option<String>,
    },
    /// Read a file through the storage layer (creates shadow if encrypted)
    Read {
        #[arg(help = "Path to the file to read (relative to SyftBox data directory)")]
        file: PathBuf,

        #[arg(long, help = "Output the file contents to stdout")]
        output: bool,
    },
}

pub async fn handle(command: SycCommands, config: &Config) -> Result<()> {
    match command {
        SycCommands::Import {
            bundle,
            expected_identity,
        } => {
            let data_root = config.get_syftbox_data_dir()?;
            let vault_home = get_biovault_home()?;
            let vault_path = syc_utils::vault_path_for_home(&vault_home);
            let parsed = syc_utils::import_public_bundle(
                &bundle,
                expected_identity.as_deref(),
                &vault_path,
                Some(&data_root),
                Some(config.email.as_str()),
            )
            .with_context(|| format!("failed to import bundle from {}", bundle.display()))?;
            println!("✓ Imported Syft Crypto bundle for {}", parsed.identity);
            println!("  fingerprint: {}", parsed.fingerprint);
        }
        SycCommands::Read { file, output } => {
            use crate::syftbox::storage::SyftBoxStorage;

            let data_dir = config.get_syftbox_data_dir()?;
            let storage = SyftBoxStorage::new(&data_dir);

            // Make path absolute if needed
            let file_path = if file.is_absolute() {
                file.clone()
            } else {
                data_dir.join(&file)
            };

            // Read through storage layer (creates shadow if encrypted)
            let contents = storage
                .read_with_shadow(&file_path)
                .with_context(|| format!("failed to read file: {}", file_path.display()))?;

            if output {
                // Output to stdout (useful for piping)
                use std::io::Write;
                std::io::stdout().write_all(&contents)?;
            } else {
                // Just confirm read without outputting
                println!("✓ Read {} bytes from {}", contents.len(), file.display());
                println!("  (shadow created if file was encrypted)");
            }
        }
    }
    Ok(())
}
