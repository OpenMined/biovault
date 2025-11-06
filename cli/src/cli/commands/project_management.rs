use crate::cli::examples;
use crate::config;
use crate::data::{BioVaultDb, Project, ProjectYaml};
use crate::error::Result;
use crate::pipeline_spec::PipelineSpec;
use anyhow::Context;
use colored::Colorize;
use std::fs;
use std::path::{Path, PathBuf};

// Standard filenames - can be made configurable in the future
const PROJECT_YAML_FILE: &str = "project.yaml";
const PIPELINE_YAML_FILE: &str = "pipeline.yaml";

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

pub fn create_project_record(
    name: String,
    example: Option<String>,
    target_dir: Option<PathBuf>,
) -> Result<Project> {
    let project_name = name.trim();
    if project_name.is_empty() {
        return Err(anyhow::anyhow!("Project name cannot be empty").into());
    }
    if project_name == "."
        || project_name == ".."
        || project_name.contains('/')
        || project_name.contains('\\')
    {
        return Err(anyhow::anyhow!("Project name contains invalid characters").into());
    }

    let example = example.and_then(|raw| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    });

    let email_config = crate::config::Config::load()
        .map(|cfg| cfg.email)
        .unwrap_or_default();
    let email_trimmed = email_config.trim().to_string();

    let db = BioVaultDb::new()?;
    if let Some(existing) = db.get_project(project_name)? {
        return Err(anyhow::anyhow!(
            "Project '{}' already exists (id: {}). Please choose a different name.",
            project_name,
            existing.id
        )
        .into());
    }

    let (project_dir, cleanup_on_error) = if let Some(dir) = target_dir {
        let project_dir = dir;

        if project_dir.exists() {
            if !project_dir.is_dir() {
                return Err(anyhow::anyhow!(
                    "Target path '{}' exists but is not a directory",
                    project_dir.display()
                )
                .into());
            }

            let yaml_path = project_dir.join(PROJECT_YAML_FILE);
            if yaml_path.exists() {
                return Err(anyhow::anyhow!(
                    "A project.yaml already exists at {}",
                    yaml_path.display()
                )
                .into());
            }

            (project_dir, false)
        } else {
            fs::create_dir_all(&project_dir)?;
            (project_dir, true)
        }
    } else {
        let biovault_home = config::get_biovault_home()?;
        let projects_dir = biovault_home.join("projects");
        fs::create_dir_all(&projects_dir)?;

        let project_dir = projects_dir.join(project_name);
        if project_dir.exists() {
            return Err(anyhow::anyhow!(
                "A project named '{}' already exists at {}",
                project_name,
                project_dir.display()
            )
            .into());
        }

        fs::create_dir_all(&project_dir)?;
        (project_dir, true)
    };

    let project_dir_for_cleanup = project_dir.clone();
    let project_name_owned = project_name.to_string();

    let build_result: Result<Project> = (|| {
        if let Some(ref example_name) = example {
            examples::write_example_to_directory(example_name, &project_dir).with_context(
                || {
                    format!(
                        "Failed to scaffold project from example '{}': {}",
                        example_name,
                        project_dir.display()
                    )
                },
            )?;
        } else {
            use crate::project_spec::{self, ProjectSpec};

            let author_placeholder = if email_trimmed.is_empty() {
                "user@example.com".to_string()
            } else {
                email_trimmed.clone()
            };

            // Remove the directory we just created so scaffold_from_spec can create it
            if cleanup_on_error {
                fs::remove_dir_all(&project_dir)?;
            }

            let minimal_spec = ProjectSpec {
                name: project_name_owned.clone(),
                author: author_placeholder,
                workflow: "workflow.nf".to_string(),
                template: Some("dynamic-nextflow".to_string()),
                version: Some("1.0.0".to_string()),
                assets: vec![],
                parameters: vec![],
                inputs: vec![],
                outputs: vec![],
            };

            project_spec::scaffold_from_spec(minimal_spec, &project_dir)
                .context("Failed to scaffold dynamic-nextflow project")?;
        }

        let yaml_path = project_dir.join(PROJECT_YAML_FILE);
        let yaml_str = fs::read_to_string(&yaml_path).context("Failed to read project.yaml")?;

        let mut yaml_value: serde_yaml::Value =
            serde_yaml::from_str(&yaml_str).context("Invalid project.yaml format")?;

        let mut author_field = email_trimmed.clone();
        let mut workflow_field = "workflow.nf".to_string();
        let mut template_field = example
            .clone()
            .unwrap_or_else(|| "dynamic-nextflow".to_string());
        let mut version_field = "1.0.0".to_string();

        if let serde_yaml::Value::Mapping(ref mut map) = yaml_value {
            let name_key = serde_yaml::Value::String("name".to_string());
            map.insert(
                name_key,
                serde_yaml::Value::String(project_name_owned.clone()),
            );

            let author_key = serde_yaml::Value::String("author".to_string());
            if !email_trimmed.is_empty() {
                map.insert(
                    author_key.clone(),
                    serde_yaml::Value::String(email_trimmed.clone()),
                );
            }
            if let Some(serde_yaml::Value::String(existing)) = map.get(&author_key) {
                if !existing.trim().is_empty() {
                    author_field = existing.clone();
                }
            }

            let workflow_key = serde_yaml::Value::String("workflow".to_string());
            if let Some(serde_yaml::Value::String(existing)) = map.get(&workflow_key) {
                if !existing.trim().is_empty() {
                    workflow_field = existing.clone();
                }
            }

            let template_key = serde_yaml::Value::String("template".to_string());
            if !map.contains_key(&template_key) {
                map.insert(
                    template_key.clone(),
                    serde_yaml::Value::String(template_field.clone()),
                );
            }
            if let Some(serde_yaml::Value::String(existing)) = map.get(&template_key) {
                if !existing.trim().is_empty() {
                    template_field = existing.clone();
                }
            }

            let version_key = serde_yaml::Value::String("version".to_string());
            if let Some(serde_yaml::Value::String(existing)) = map.get(&version_key) {
                if !existing.trim().is_empty() {
                    version_field = existing.clone();
                }
            }
        } else {
            return Err(anyhow::anyhow!("project.yaml has unexpected structure").into());
        }

        let updated_yaml =
            serde_yaml::to_string(&yaml_value).context("Failed to serialize project.yaml")?;
        fs::write(&yaml_path, updated_yaml).context("Failed to update project.yaml")?;

        db.register_project(
            &project_name_owned,
            &version_field,
            &author_field,
            &workflow_field,
            &template_field,
            &project_dir,
        )?;

        let created = db
            .get_project(&project_name_owned)?
            .ok_or_else(|| anyhow::anyhow!("Project registration completed but not found"))?;

        Ok(created)
    })();

    match build_result {
        Ok(project) => Ok(project),
        Err(err) => {
            if cleanup_on_error {
                let _ = fs::remove_dir_all(project_dir_for_cleanup);
            }
            Err(err)
        }
    }
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

    // Create project directory path (including version to allow multiple versions)
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;
    let dir_name = format!("{}-{}", project_name, project_yaml.version);
    let project_dir = projects_dir.join(&dir_name);

    // Check if project already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", project_name, project_yaml.version);
        if let Some(mut existing) = db.get_project(&identifier)? {
            // Project with this version exists in DB - check if files are valid
            let existing_path = PathBuf::from(&existing.project_path);
            let existing_yaml = existing_path.join(PROJECT_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid project exists - check if path needs updating to new convention
                let expected_dir_name = format!("{}-{}", project_name, project_yaml.version);
                let expected_path = projects_dir.join(&expected_dir_name);

                if existing_path != expected_path && !existing_path.ends_with(&expected_dir_name) {
                    // Path is in old format (without version) - migrate to new format
                    if !quiet {
                        println!(
                            "   🔄 Migrating project '{}' to versioned path format",
                            project_name.dimmed()
                        );
                    }

                    // Move directory to new location if needed
                    if !expected_path.exists() {
                        fs::rename(&existing_path, &expected_path).with_context(|| {
                            format!(
                                "Failed to migrate project directory from {} to {}",
                                existing_path.display(),
                                expected_path.display()
                            )
                        })?;
                    } else {
                        // New path already exists - remove old one
                        fs::remove_dir_all(&existing_path).ok();
                    }

                    // Update DB with new path
                    db.update_project(
                        &project_name,
                        &project_yaml.version,
                        &project_yaml.author,
                        &project_yaml.workflow,
                        &project_yaml.template,
                        &expected_path,
                    )?;

                    // Update existing struct
                    existing.project_path = expected_path.to_string_lossy().to_string();
                }

                // Reuse existing valid project (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Project '{}' version {} already registered, reusing (id: {})",
                        project_name.dimmed(),
                        project_yaml.version.dimmed(),
                        existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Project DB record exists but files are missing/invalid - clean it up
                if !quiet {
                    println!(
                        "   🧹 Cleaning up orphaned project record for '{}@{}'",
                        project_name.dimmed(),
                        project_yaml.version.dimmed()
                    );
                }
                db.delete_project(&identifier)?;
                // Continue with fresh import below
            }
        }
    }

    // Handle directory conflicts
    if project_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&project_dir)?;
        } else {
            // Directory exists but not in DB (orphaned) - clean it up
            if !quiet {
                println!(
                    "   🧹 Removing orphaned project directory: {}",
                    project_dir.display().to_string().dimmed()
                );
            }
            fs::remove_dir_all(&project_dir)?;
        }
    }

    fs::create_dir_all(&project_dir)?;

    // Save project.yaml
    let yaml_path = project_dir.join(PROJECT_YAML_FILE);
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
        // Check if this exact version exists
        let identifier = format!("{}@{}", project_name, project_yaml.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        } else {
            db.register_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &project_yaml.version,
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

/// Copy a local project to the managed directory (~/.biovault/projects/)
/// This mirrors the behavior of import_from_url but copies from local filesystem
async fn copy_local_project_to_managed(
    db: &BioVaultDb,
    source_path: &Path,
    overwrite: bool,
    quiet: bool,
) -> Result<Project> {
    if !source_path.exists() {
        return Err(
            anyhow::anyhow!("Source path does not exist: {}", source_path.display()).into(),
        );
    }

    if !source_path.is_dir() {
        return Err(
            anyhow::anyhow!("Source path is not a directory: {}", source_path.display()).into(),
        );
    }

    // Look for project.yaml
    let yaml_path = source_path.join(PROJECT_YAML_FILE);
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No project.yaml found in directory: {}",
            source_path.display()
        )
        .into());
    }

    // Read and parse project.yaml
    let yaml_str = fs::read_to_string(&yaml_path)?;
    let project_yaml: ProjectYaml = serde_yaml::from_str(&yaml_str)?;

    let project_name = project_yaml.name.clone();

    // Create project directory in managed location (like import_from_url does)
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;
    let dir_name = format!("{}-{}", project_name, project_yaml.version);
    let project_dir = projects_dir.join(&dir_name);

    // Check if project already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", project_name, project_yaml.version);
        if let Some(existing) = db.get_project(&identifier)? {
            let existing_path = PathBuf::from(&existing.project_path);
            let existing_yaml = existing_path.join(PROJECT_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid project exists - reuse it (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Project '{}' version {} already imported, reusing (id: {})",
                        project_name, project_yaml.version, existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Project DB record exists but files are missing/invalid - clean it up
                if !quiet {
                    println!(
                        "   🧹 Cleaning up orphaned project record for '{}@{}'",
                        project_name, project_yaml.version
                    );
                }
                db.delete_project(&identifier)?;
            }
        }
    }

    // Handle directory conflicts
    if project_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&project_dir)?;
        } else {
            // Directory exists but not in DB (orphaned) - clean it up
            if !quiet {
                println!(
                    "   🧹 Removing orphaned project directory: {}",
                    project_dir.display()
                );
            }
            fs::remove_dir_all(&project_dir)?;
        }
    }

    fs::create_dir_all(&project_dir)?;

    // Copy project.yaml
    let dest_yaml_path = project_dir.join(PROJECT_YAML_FILE);
    fs::copy(&yaml_path, &dest_yaml_path)?;

    if !quiet {
        println!("✓ Copied project.yaml");
    }

    // Copy workflow file
    let workflow_source = source_path.join(&project_yaml.workflow);
    if workflow_source.exists() {
        let workflow_dest = project_dir.join(&project_yaml.workflow);
        if let Some(parent) = workflow_dest.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::copy(&workflow_source, &workflow_dest)?;
        if !quiet {
            println!("✓ Copied {}", project_yaml.workflow);
        }
    } else if !quiet {
        println!(
            "⚠️  Warning: workflow file '{}' not found in source",
            project_yaml.workflow
        );
    }

    // Copy assets
    if !project_yaml.assets.is_empty() {
        let assets_dir = project_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Copying {} assets...", project_yaml.assets.len());
        }

        for asset in &project_yaml.assets {
            let asset_source = source_path.join("assets").join(asset);
            if asset_source.exists() {
                let asset_dest = assets_dir.join(asset);

                // Create parent directories if asset has a path
                if let Some(parent) = asset_dest.parent() {
                    fs::create_dir_all(parent)?;
                }

                fs::copy(&asset_source, &asset_dest)?;

                if !quiet {
                    println!("  ✓ {}", asset);
                }
            } else if !quiet {
                println!("  ⚠️  Warning: asset '{}' not found", asset);
            }
        }

        if !quiet {
            println!("✓ Copied all assets");
        }
    }

    // Register in database
    if overwrite {
        let identifier = format!("{}@{}", project_name, project_yaml.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        } else {
            db.register_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_dir,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &project_yaml.version,
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
    let yaml_path = project_path.join(PROJECT_YAML_FILE);
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

    // Check if project already exists in DB
    if !overwrite {
        if let Some(existing) = db.get_project(&project_name)? {
            // Check if existing project points to the same path
            let existing_path = PathBuf::from(&existing.project_path);
            let canonical_existing = existing_path.canonicalize().ok();
            let canonical_new = project_path.canonicalize().ok();

            if canonical_existing == canonical_new && canonical_existing.is_some() {
                // Same project, already registered - reuse it (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Project '{}' already registered at this path (id: {})",
                        project_name.dimmed(),
                        existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Different path - this is a conflict
                return Err(anyhow::anyhow!(
                    "Project '{}' already exists (id: {}) at different path: {}. Use --overwrite to replace.",
                    project_name,
                    existing.id,
                    existing.project_path
                )
                .into());
            }
        }
    }

    // Register in database (using original path, not copied)
    if overwrite {
        let identifier = format!("{}@{}", project_name, project_yaml.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_path,
            )?;
        } else {
            db.register_project(
                &project_name,
                &project_yaml.version,
                &project_yaml.author,
                &project_yaml.workflow,
                &project_yaml.template,
                &project_path,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &project_yaml.version,
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
    let yaml_path = Path::new(&project.project_path).join(PROJECT_YAML_FILE);
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

/// Context for resolving pipeline step dependencies
#[derive(Debug, Clone)]
pub enum DependencyContext {
    /// Dependencies are resolved relative to a GitHub base URL
    GitHub { base_url: String },
    /// Dependencies are resolved relative to a local filesystem path
    Local { base_path: PathBuf },
}

/// Resolve and import all pipeline step dependencies
///
/// This function handles importing project dependencies for pipeline steps,
/// rewriting the pipeline spec to use registered project names.
///
/// Returns `true` if any dependencies were imported and the spec was updated.
pub async fn resolve_pipeline_dependencies(
    spec: &mut PipelineSpec,
    dependency_context: &DependencyContext,
    pipeline_yaml_path: &Path,
    overwrite: bool,
    quiet: bool,
) -> Result<bool> {
    use colored::Colorize;

    let db = BioVaultDb::new()?;
    let mut any_rewritten = false;

    if spec.steps.is_empty() {
        return Ok(false);
    }

    if !quiet {
        println!(
            "\n{} Importing dependencies ({} steps):",
            "📦".cyan(),
            spec.steps.len()
        );
    }

    // Collect step updates to avoid borrow checker issues
    let mut step_updates: Vec<(usize, String)> = Vec::new();

    for (index, step) in spec.steps.iter().enumerate() {
        if let Some(uses) = &step.uses {
            // Skip absolute local paths
            if uses.starts_with('/') {
                if !quiet {
                    println!(
                        "   {} Step '{}' uses absolute path (keeping as-is)",
                        "ℹ️ ".cyan(),
                        step.id.dimmed()
                    );
                }
                continue;
            }

            // Skip if already a simple name (already registered)
            // BUT: in Local context, ALWAYS try to resolve as path first
            let should_skip_as_registered = match dependency_context {
                DependencyContext::GitHub { .. } => {
                    // In GitHub context, no slashes + not http = registered name
                    !uses.contains('/') && !uses.starts_with("http")
                }
                DependencyContext::Local { .. } => {
                    // In Local context, never skip - always try to resolve as path first
                    // This handles cases like "apol1-classifier" which could be a relative path
                    false // Don't skip - let path resolution logic below handle it
                }
            };

            if should_skip_as_registered {
                if !quiet {
                    println!(
                        "   {} Step '{}' already uses name (keeping as-is)",
                        "ℹ️ ".cyan(),
                        step.id.dimmed()
                    );
                }
                continue;
            }

            // Resolve project reference based on context
            let (should_use_local, local_path_opt, url_opt) =
                if uses.starts_with("http://") || uses.starts_with("https://") {
                    // Absolute URL - use as-is
                    (false, None, Some(format!("{}/project.yaml", uses)))
                } else {
                    // Relative path - resolve based on context
                    match dependency_context {
                        DependencyContext::GitHub { base_url } => (
                            false,
                            None,
                            Some(format!("{}/{}/project.yaml", base_url, uses)),
                        ),
                        DependencyContext::Local { base_path } => {
                            // Resolve exactly like GitHub: base_path + uses = project location
                            // GitHub does: base_url + "/" + uses + "/project.yaml"
                            // Local does: base_path.join(uses) -> check for project.yaml there
                            let project_path = base_path.join(uses);

                            // Try to find project.yaml (same logic as GitHub URL resolution)
                            let project_yaml_path = if project_path.is_dir() {
                                project_path.join(PROJECT_YAML_FILE)
                            } else if project_path
                                .extension()
                                .map(|e| e == "yaml" || e == "yml")
                                .unwrap_or(false)
                            {
                                // Direct path to yaml file
                                project_path.clone()
                            } else {
                                // Assume it's a directory path, look for project.yaml inside
                                project_path.join(PROJECT_YAML_FILE)
                            };

                            // Check if the resolved path exists (mirrors GitHub's existence check)
                            if project_yaml_path.exists() {
                                // Path exists - import it (like GitHub downloads it)
                                let project_dir = if project_path.is_dir() {
                                    project_path
                                } else {
                                    project_yaml_path
                                        .parent()
                                        .unwrap_or(base_path)
                                        .to_path_buf()
                                };
                                (true, Some(project_dir), None)
                            } else {
                                // Path doesn't exist - check if it's already a registered project name
                                // If registered, skip this step (like GitHub does for registered names)
                                let db_check = db.get_project(uses).ok().flatten();
                                if db_check.is_some() {
                                    // It's registered - skip this dependency (already using registered name)
                                    if !quiet {
                                        println!(
                                        "   {} Step '{}' uses registered name '{}' (keeping as-is)",
                                        "ℹ️ ".cyan(),
                                        step.id.dimmed(),
                                        uses
                                    );
                                    }
                                    // Skip this step entirely - it's already using a registered project name
                                    continue;
                                }
                                // Not found as path or registered - will error later when we try to import
                                (false, None, None)
                            }
                        }
                    }
                };

            if !quiet {
                print!("   {} {} ", "•".cyan(), step.id.bold());
            }

            // Import the project (registers in DB automatically)
            // Smart resolution: check if already registered first, then copy if needed
            let import_result = if should_use_local {
                if let Some(local_path) = local_path_opt {
                    // Check if project at this path is already registered
                    let existing_projects =
                        db.list_projects().context("Failed to list projects")?;
                    let already_registered = existing_projects.iter().any(|p| {
                        PathBuf::from(&p.project_path).canonicalize().ok()
                            == local_path.canonicalize().ok()
                    });

                    if already_registered {
                        // Project already registered - parse local path to get project name and reuse
                        let yaml_path = local_path.join(PROJECT_YAML_FILE);
                        if yaml_path.exists() {
                            let yaml_str = fs::read_to_string(&yaml_path)
                                .context("Failed to read project.yaml")?;
                            let project_yaml: ProjectYaml = serde_yaml::from_str(&yaml_str)
                                .context("Failed to parse project.yaml")?;
                            let identifier =
                                format!("{}@{}", project_yaml.name, project_yaml.version);

                            match db.get_project(&identifier) {
                                Ok(Some(project)) => {
                                    // Reuse existing registered project
                                    Ok((project, true))
                                }
                                _ => {
                                    // Not found by identifier, copy to managed for consistency
                                    copy_local_project_to_managed(&db, &local_path, overwrite, true)
                                        .await
                                        .map(|p| (p, true))
                                }
                            }
                        } else {
                            return Err(anyhow::anyhow!(
                                "project.yaml not found at {}",
                                local_path.display()
                            )
                            .into());
                        }
                    } else {
                        // Not registered - copy to managed directory (like GitHub imports)
                        copy_local_project_to_managed(&db, &local_path, overwrite, true)
                            .await
                            .map(|p| (p, true))
                    }
                } else {
                    return Err(anyhow::anyhow!(
                        "Local path not available for dependency '{}'",
                        uses
                    )
                    .into());
                }
            } else if let Some(url) = url_opt {
                import_from_url(&db, &url, None, overwrite, true)
                    .await
                    .map(|p| (p, false))
            } else {
                // Fallback: try URL import as last resort for local context
                if let DependencyContext::Local { base_path: _ } = dependency_context {
                    // Could try to convert local path to URL, but this is unusual
                    return Err(anyhow::anyhow!("Cannot resolve dependency '{}': local path doesn't exist and no URL available", uses).into());
                } else {
                    return Err(anyhow::anyhow!(
                        "Cannot resolve dependency '{}': no valid source",
                        uses
                    )
                    .into());
                }
            };

            match import_result {
                Ok((project, _is_local)) => {
                    if !quiet {
                        println!("{} → {}", "✓ imported".green(), project.name.green());
                    }
                    // Store update for later (avoid borrow checker issue)
                    step_updates.push((index, project.name.clone()));
                    any_rewritten = true;
                }
                Err(e) => {
                    if !quiet {
                        println!("{}: {}", "⚠️  failed".yellow(), e.to_string().dimmed());
                    }
                    // Continue to next step - don't fail entire import
                }
            }
        }
    }

    // Apply step updates now that we're done iterating
    for (index, project_name) in step_updates {
        if let Some(updated_step) = spec.steps.get_mut(index) {
            updated_step.uses = Some(project_name);
        }
    }

    // Always save the spec (with description preserved) after dependency resolution
    // This ensures description is preserved even if dependencies weren't rewritten
    spec.save(pipeline_yaml_path)
        .context("Failed to save pipeline.yaml with description")?;

    if any_rewritten
        && !quiet {
            println!(
                "\n{} Updated pipeline.yaml to use registered names...",
                "🔧".cyan()
            );
        }

    Ok(any_rewritten)
}

/// Import a pipeline from URL with all its step dependencies
pub async fn import_pipeline_with_deps(
    url: &str,
    name_override: Option<String>,
    overwrite: bool,
) -> Result<String> {
    use crate::pipeline_spec::PipelineSpec;
    use colored::Colorize;

    // Convert GitHub URL to raw URL
    let raw_url = url
        .replace("github.com", "raw.githubusercontent.com")
        .replace("/blob/", "/");

    println!("{} Downloading pipeline from {}", "📥".cyan(), url.cyan());

    // Download pipeline YAML
    let yaml_content = download_file(&raw_url).await?;
    let yaml_str = String::from_utf8(yaml_content).context("Invalid UTF-8 in pipeline.yaml")?;

    // Parse pipeline spec
    let mut spec: PipelineSpec =
        serde_yaml::from_str(&yaml_str).context("Failed to parse pipeline.yaml")?;

    let pipeline_name = name_override.unwrap_or_else(|| spec.name.clone());

    println!("{} Pipeline: {}", "📦".cyan(), pipeline_name.bold());
    println!("   Steps: {}", spec.steps.len());

    // Create pipeline directory
    let biovault_home = crate::config::get_biovault_home()?;
    let pipelines_dir = biovault_home.join("pipelines");
    fs::create_dir_all(&pipelines_dir)?;

    let pipeline_dir = pipelines_dir.join(&pipeline_name);

    if pipeline_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&pipeline_dir)?;
        } else {
            return Err(anyhow::anyhow!(
                "Pipeline directory already exists: {}. Use --overwrite to replace.",
                pipeline_dir.display()
            )
            .into());
        }
    }

    fs::create_dir_all(&pipeline_dir)?;

    let pipeline_yaml_path = pipeline_dir.join(PIPELINE_YAML_FILE);

    // Extract base URL for resolving relative project paths
    let base_url = if let Some(idx) = url.rfind('/') {
        url[..idx].to_string()
    } else {
        url.to_string()
    };

    // Import each step's project and rewrite YAML to use registered names
    resolve_pipeline_dependencies(
        &mut spec,
        &DependencyContext::GitHub { base_url },
        &pipeline_yaml_path,
        overwrite,
        false, // quiet = false for CLI output
    )
    .await?;

    // Register pipeline in database (check for existing if overwrite)
    let db = BioVaultDb::new()?;
    let pipeline_id = if overwrite {
        // Check if pipeline with this name already exists
        let existing_pipelines = db.list_pipelines()?;
        if let Some(existing) = existing_pipelines.iter().find(|p| p.name == pipeline_name) {
            // Delete existing pipeline from DB
            db.delete_pipeline(existing.id)?;
        }
        // Register the new pipeline
        db.register_pipeline(&pipeline_name, &pipeline_dir.to_string_lossy())?
    } else {
        db.register_pipeline(&pipeline_name, &pipeline_dir.to_string_lossy())?
    };

    println!(
        "\n{} Pipeline '{}' imported successfully!",
        "✅".green().bold(),
        pipeline_name.bold()
    );
    println!(
        "   Location: {}",
        pipeline_dir.display().to_string().dimmed()
    );
    println!("   ID: {}", pipeline_id);

    Ok(pipeline_dir.to_string_lossy().to_string())
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
        fs::write(project_dir.join(PROJECT_YAML_FILE), yaml_content).unwrap();
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
            "1.0.0",
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
