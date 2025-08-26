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

    // Get or default folder path
    let folder_path = folder.unwrap_or_else(|| format!("./{}", project_name));
    let project_dir = PathBuf::from(&folder_path);

    // Check if folder already exists
    if project_dir.exists() {
        return Err(anyhow::anyhow!(
            "Folder '{}' already exists. Please choose a different name or location.",
            project_dir.display()
        )
        .into());
    }

    // Get email from config
    let home_dir =
        dirs::home_dir().ok_or_else(|| anyhow::anyhow!("Could not determine home directory"))?;
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
    println!("\nProject structure:");
    println!("  {}/", project_name);
    println!("    ├── project.yaml");
    println!("    ├── assets/");
    println!("    └── workflow.nf");

    Ok(())
}
