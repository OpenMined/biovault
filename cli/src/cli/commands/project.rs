use crate::{config::Config, error::Result};
use anyhow::Context;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;

const PROJECT_YAML_TEMPLATE: &str = include_str!("../../templates/project.yaml");
const WORKFLOW_NF_TEMPLATE: &str = include_str!("../../templates/workflow.nf");

pub async fn create(name: Option<String>, folder: Option<String>) -> Result<()> {
    // Get or prompt for project name
    let project_name = match name {
        Some(n) => n,
        None => {
            print!("Project name: ");
            io::stdout().flush()?;
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            input.trim().to_string()
        }
    };

    if project_name.is_empty() {
        return Err(anyhow::anyhow!("Project name cannot be empty").into());
    }

    // Get or determine folder path
    let project_dir = match folder {
        Some(f) => PathBuf::from(f),
        None => {
            // Default to current working directory + project name
            let cwd = std::env::current_dir().context("Failed to get current working directory")?;
            let default_path = cwd.join(&project_name);

            // Ask for confirmation
            println!("Project will be created at: {}", default_path.display());
            print!("Continue? [Y/n]: ");
            io::stdout().flush()?;

            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            let input = input.trim().to_lowercase();

            if input.is_empty() || input == "y" || input == "yes" {
                default_path
            } else if input == "n" || input == "no" {
                print!("Enter custom path: ");
                io::stdout().flush()?;
                let mut custom_path = String::new();
                io::stdin().read_line(&mut custom_path)?;
                PathBuf::from(custom_path.trim())
            } else {
                return Err(
                    anyhow::anyhow!("Invalid response. Please run the command again.").into(),
                );
            }
        }
    };

    // Check if folder already exists
    if project_dir.exists() {
        return Err(anyhow::anyhow!(
            "Folder '{}' already exists. Please choose a different name or location.",
            project_dir.display()
        )
        .into());
    }

    // Get email from config
    // For testing: allow overriding home directory via environment variable
    let home_dir = if let Ok(test_home) = std::env::var("BIOVAULT_TEST_HOME") {
        std::path::PathBuf::from(test_home)
    } else {
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?
    };
    let config_file = home_dir.join(".biovault").join("config.yaml");

    let email = if config_file.exists() {
        let config = Config::from_file(&config_file)?;
        config.email
    } else {
        return Err(anyhow::anyhow!(
            "BioVault not initialized. Please run 'bv init <email>' first"
        )
        .into());
    };

    // Create project structure
    println!("Creating project '{}'...", project_name);

    // Create directories
    fs::create_dir_all(&project_dir).with_context(|| {
        format!(
            "Failed to create project directory: {}",
            project_dir.display()
        )
    })?;

    let assets_dir = project_dir.join("assets");
    fs::create_dir_all(&assets_dir).with_context(|| {
        format!(
            "Failed to create assets directory: {}",
            assets_dir.display()
        )
    })?;

    // Create project.yaml
    let project_yaml_content = PROJECT_YAML_TEMPLATE
        .replace("{project_name}", &project_name)
        .replace("{email}", &email);

    let project_yaml_path = project_dir.join("project.yaml");
    fs::write(&project_yaml_path, project_yaml_content).with_context(|| {
        format!(
            "Failed to write project.yaml: {}",
            project_yaml_path.display()
        )
    })?;

    // Create workflow.nf
    let workflow_path = project_dir.join("workflow.nf");
    fs::write(&workflow_path, WORKFLOW_NF_TEMPLATE)
        .with_context(|| format!("Failed to write workflow.nf: {}", workflow_path.display()))?;

    println!(
        "✓ Project created successfully at: {}",
        project_dir.display()
    );

    // Get the actual folder name for display (might be different from project name)
    let folder_name = project_dir
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&project_name);

    println!("\nProject structure:");
    println!("  {}/", folder_name);
    println!("    ├── project.yaml");
    println!("    ├── assets/");
    println!("    └── workflow.nf");

    Ok(())
}
