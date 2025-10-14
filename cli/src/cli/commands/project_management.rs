use crate::config;
use crate::data::{BioVaultDb, Project, ProjectYaml};
use crate::error::Result;
use anyhow::Context;
use colored::Colorize;
use std::fs;
use std::path::{Path, PathBuf};

/// Import a project from URL or register a local project directory
pub async fn import(
    source: String,
    name: Option<String>,
    overwrite: bool,
    format: Option<String>,
) -> Result<()> {
    let db = BioVaultDb::new()?;
    let quiet = format.as_deref() == Some("json");

    // Determine if source is a URL or local path
    let is_url = source.starts_with("http://") || source.starts_with("https://");

    let project = if is_url {
        import_from_url(&db, &source, name.clone(), overwrite, quiet).await?
    } else {
        import_from_local(&db, &source, name.clone(), overwrite, quiet)?
    };

    if quiet {
        let response = crate::data::CliResponse::new(project);
        println!("{}", serde_json::to_string_pretty(&response)?);
    } else if is_url {
        println!(
            "\n✅ Project '{}' imported successfully!",
            project.name.green()
        );
        println!("   Location: {}", project.project_path);
    } else {
        println!(
            "\n✅ Project '{}' registered successfully!",
            project.name.green()
        );
        println!("   Location: {}", project.project_path);
        println!("\n💡 This project is referenced in-place.");
        println!("   Any changes you make will be reflected when you run it.");
    }

    Ok(())
}

pub async fn import_project_record(
    source: String,
    name: Option<String>,
    overwrite: bool,
) -> Result<Project> {
    let db = BioVaultDb::new()?;
    let is_url = source.starts_with("http://") || source.starts_with("https://");

    if is_url {
        import_from_url(&db, &source, name, overwrite, true).await
    } else {
        import_from_local(&db, &source, name, overwrite, true)
    }
}

/// Create a new project and register it in the database (for desktop app)
pub async fn create_project_record(name: String, example: Option<String>) -> Result<Project> {
    let db = BioVaultDb::new()?;

    // Check if project already exists
    if let Some(existing) = db.get_project(&name)? {
        return Err(anyhow::anyhow!(
            "Project '{}' already exists (id: {}). Please choose a different name.",
            name,
            existing.id
        )
        .into());
    }

    // Get BIOVAULT_HOME and create projects directory
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;

    let project_folder = projects_dir.join(&name);
    let project_folder_str = project_folder.to_string_lossy().to_string();

    // Use the existing project::create function to create files
    crate::cli::commands::project::create(
        Some(name.clone()),
        Some(project_folder_str.clone()),
        example,
    )
    .await?;

    // Load the created project.yaml to get details
    let yaml_path = project_folder.join("project.yaml");
    let yaml_content =
        fs::read_to_string(&yaml_path).context("Failed to read project.yaml after creation")?;

    let project_yaml: ProjectYaml =
        serde_yaml::from_str(&yaml_content).context("Failed to parse project.yaml")?;

    // Register the project in the database
    db.register_project(
        &name,
        &project_yaml.author,
        &project_yaml.workflow,
        &project_yaml.template,
        &project_folder,
    )?;

    // Get the registered project from the database and return it
    let project = db
        .get_project(&name)?
        .ok_or_else(|| anyhow::anyhow!("Project was created but not found in database"))?;

    Ok(project)
}

async fn import_from_url(
    db: &BioVaultDb,
    url: &str,
    name_override: Option<String>,
    overwrite: bool,
    quiet: bool,
) -> Result<Project> {
    if !quiet {
        println!("📥 Importing project from: {}", url.cyan());
    }

    // Convert github.com URLs to raw.githubusercontent.com
    let raw_url = url
        .replace("github.com", "raw.githubusercontent.com")
        .replace("/blob/", "/");

    if !quiet {
        println!("📄 Downloading project.yaml...");
    }

    // Download project.yaml
    let yaml_content = download_file(&raw_url).await?;
    let yaml_str = String::from_utf8(yaml_content).context("Invalid UTF-8 in project.yaml")?;

    // Parse project.yaml
    let project_yaml: ProjectYaml =
        serde_yaml::from_str(&yaml_str).context("Failed to parse project.yaml")?;

    let project_name = name_override.unwrap_or(project_yaml.name.clone());

    // Check if project exists
    if !overwrite {
        if let Some(existing) = db.get_project(&project_name)? {
            return Err(anyhow::anyhow!(
                "Project '{}' already exists (id: {}). Use --overwrite to replace.",
                project_name,
                existing.id
            )
            .into());
        }
    }

    // Create project directory in BIOVAULT_HOME/projects/
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;

    let project_dir = projects_dir.join(&project_name);
    if project_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&project_dir)?;
        } else {
            return Err(anyhow::anyhow!(
                "Project directory already exists: {}",
                project_dir.display()
            )
            .into());
        }
    }

    fs::create_dir_all(&project_dir)?;

    // Save project.yaml
    let yaml_path = project_dir.join("project.yaml");
    fs::write(&yaml_path, yaml_str)?;

    if !quiet {
        println!("✓ Saved project.yaml");
    }

    // Download workflow file
    let base_url = raw_url
        .rsplit_once('/')
        .map(|(base, _)| base)
        .ok_or_else(|| anyhow::anyhow!("Invalid URL format"))?;

    let workflow_url = format!("{}/{}", base_url, project_yaml.workflow);

    if !quiet {
        println!("📄 Downloading {}...", project_yaml.workflow);
    }

    let workflow_content = download_file(&workflow_url).await?;
    let workflow_path = project_dir.join(&project_yaml.workflow);
    fs::write(&workflow_path, workflow_content)?;

    if !quiet {
        println!("✓ Saved {}", project_yaml.workflow);
    }

    // Download assets
    if !project_yaml.assets.is_empty() {
        let assets_dir = project_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Downloading {} assets...", project_yaml.assets.len());
        }

        let assets_url = format!("{}/assets", base_url);

        for asset in &project_yaml.assets {
            let asset_url = format!("{}/{}", assets_url, asset);

            if !quiet {
                println!("  - {}", asset);
            }

            let asset_content = download_file(&asset_url).await?;
            let asset_path = assets_dir.join(asset);

            // Create parent directories if asset has a path
            if let Some(parent) = asset_path.parent() {
                fs::create_dir_all(parent)?;
            }

            fs::write(&asset_path, asset_content)?;
        }

        if !quiet {
            println!("✓ Downloaded all assets");
        }
    }

    // Register in database
    if overwrite {
        if db.get_project(&project_name)?.is_some() {
            db.update_project(
                &project_name,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        } else {
            db.register_project(
                &project_name,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &project_yaml.author,
            &project_yaml.workflow,
            &project_yaml.template,
            &project_dir,
        )?;
    }

    let project = db
        .get_project(&project_name)?
        .ok_or_else(|| anyhow::anyhow!("Project '{}' not found after import", project_name))?;

    Ok(project)
}

fn import_from_local(
    db: &BioVaultDb,
    path: &str,
    name_override: Option<String>,
    overwrite: bool,
    quiet: bool,
) -> Result<Project> {
    let project_path = PathBuf::from(path);
    if !project_path.exists() {
        return Err(anyhow::anyhow!("Path does not exist: {}", path).into());
    }

    if !project_path.is_dir() {
        return Err(anyhow::anyhow!("Path is not a directory: {}", path).into());
    }

    // Look for project.yaml
    let yaml_path = project_path.join("project.yaml");
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No project.yaml found in directory: {}",
            project_path.display()
        )
        .into());
    }

    if !quiet {
        println!("📁 Registering local project: {}", path.cyan());
    }

    // Parse project.yaml
    let yaml_str = fs::read_to_string(&yaml_path)?;
    let project_yaml: ProjectYaml = serde_yaml::from_str(&yaml_str)?;

    let project_name = name_override.unwrap_or(project_yaml.name.clone());

    // Check if project exists
    if !overwrite {
        if let Some(existing) = db.get_project(&project_name)? {
            return Err(anyhow::anyhow!(
                "Project '{}' already exists (id: {}). Use --overwrite to replace.",
                project_name,
                existing.id
            )
            .into());
        }
    }

    // Register in database (using original path, not copied)
    if overwrite {
        if db.get_project(&project_name)?.is_some() {
            db.update_project(
                &project_name,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_path,
            )?;
        } else {
            db.register_project(
                &project_name,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_path,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &project_yaml.author,
            &project_yaml.workflow,
            &project_yaml.template,
            &project_path,
        )?;
    }

    let project = db
        .get_project(&project_name)?
        .ok_or_else(|| anyhow::anyhow!("Project '{}' not found after import", project_name))?;

    Ok(project)
}

/// List all projects
pub fn list(format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;
    let projects = db.list_projects()?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(projects);
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    if projects.is_empty() {
        println!("📭 No projects found");
        println!("\n💡 Import a project:");
        println!("   bv project import <url>");
        println!("   bv project import /path/to/project");
        return Ok(());
    }

    println!("\n📦 {} project(s):\n", projects.len());
    println!(
        "{:<4} {:<20} {:<25} {:<15} {:<20}",
        "ID".bold(),
        "Name".bold(),
        "Author".bold(),
        "Workflow".bold(),
        "Created".bold()
    );
    println!("{}", "─".repeat(90));

    for project in projects {
        println!(
            "{:<4} {:<20} {:<25} {:<15} {:<20}",
            project.id,
            truncate(&project.name, 18),
            truncate(&project.author, 23),
            truncate(&project.workflow, 13),
            &project.created_at[..19] // Trim timestamp to date+time
        );
    }

    println!();

    Ok(())
}

/// Show detailed information about a project
pub fn show(identifier: String, format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;
    let project = db
        .get_project(&identifier)?
        .ok_or_else(|| anyhow::anyhow!("Project '{}' not found", identifier))?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(project);
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    let run_count = db.count_project_runs(project.id)?;

    println!("\n{}", "═".repeat(80));
    println!("📦 {}", project.name.bold().cyan());
    println!("{}", "─".repeat(80));
    println!("ID:           {}", project.id);
    println!("Author:       {}", project.author);
    println!("Workflow:     {}", project.workflow);
    println!("Template:     {}", project.template);
    println!("Location:     {}", project.project_path);
    println!("Created:      {}", project.created_at);
    println!("Runs:         {}", run_count);
    println!("{}", "═".repeat(80));
    println!();

    // Show project.yaml contents
    let yaml_path = Path::new(&project.project_path).join("project.yaml");
    if yaml_path.exists() {
        println!("📄 project.yaml:");
        println!("{}", "─".repeat(80));
        let yaml_content = fs::read_to_string(&yaml_path)?;
        for line in yaml_content.lines().take(20) {
            println!("   {}", line.dimmed());
        }
        println!("{}", "─".repeat(80));
        println!();
    }

    Ok(())
}

/// Delete a project
pub fn delete(identifier: String, keep_files: bool, format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get project before deleting
    let project = db.delete_project(&identifier)?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(serde_json::json!({
            "deleted": project,
            "files_kept": keep_files
        }));
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    println!(
        "✅ Deleted project '{}' from database",
        project.name.green()
    );

    // Delete files if requested
    if !keep_files {
        let project_path = Path::new(&project.project_path);
        if project_path.exists() {
            fs::remove_dir_all(project_path)?;
            println!("✅ Deleted project files: {}", project_path.display());
        }
    } else {
        println!("📁 Project files kept: {}", project.project_path);
    }

    Ok(())
}

// Helper function to download files
async fn download_file(url: &str) -> Result<Vec<u8>> {
    let response = reqwest::get(url)
        .await
        .map_err(|e| anyhow::anyhow!("Failed to download {}: {}", url, e))?;

    if !response.status().is_success() {
        return Err(
            anyhow::anyhow!("Failed to download {}: HTTP {}", url, response.status()).into(),
        );
    }

    let bytes = response
        .bytes()
        .await
        .map_err(|e| anyhow::anyhow!("Failed to read response bytes: {}", e))?
        .to_vec();
    Ok(bytes)
}

// Helper function to truncate strings
fn truncate(s: &str, max_len: usize) -> String {
    if s.len() > max_len {
        format!("{}...", &s[..max_len - 3])
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    fn setup_test() -> TempDir {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path());
        tmp
    }

    fn teardown_test() {
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn test_import_from_local() {
        let tmp = setup_test();

        // Create a test project directory
        let project_dir = tmp.path().join("test-project");
        fs::create_dir_all(&project_dir).unwrap();

        let yaml_content = r#"
name: test-project
author: test@example.com
workflow: workflow.nf
template: default
assets: []
"#;
        fs::write(project_dir.join("project.yaml"), yaml_content).unwrap();
        fs::write(project_dir.join("workflow.nf"), "// workflow").unwrap();

        // Import the project
        let project = import_from_local(
            &BioVaultDb::new().unwrap(),
            project_dir.to_str().unwrap(),
            None,
            false,
            true,
        )
        .expect("local project import should succeed");
        assert_eq!(project.name, "test-project");

        // Verify it's in the database
        let db = BioVaultDb::new().unwrap();
        let project = db.get_project("test-project").unwrap();
        assert!(project.is_some());

        teardown_test();
    }

    #[test]
    fn test_list_projects() {
        let tmp = setup_test();
        let db = BioVaultDb::new().unwrap();

        // Create some test projects
        let project_path = tmp.path().join("project1");
        fs::create_dir_all(&project_path).unwrap();
        db.register_project(
            "project1",
            "author@example.com",
            "workflow.nf",
            "default",
            &project_path,
        )
        .unwrap();

        // List should work
        let result = list(None);
        assert!(result.is_ok());

        teardown_test();
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("short", 10), "short");
        assert_eq!(truncate("verylongstring", 10), "verylon...");
    }
}
