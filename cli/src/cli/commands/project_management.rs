use crate::cli::examples;
use crate::config;
use crate::data::{BioVaultDb, Project};
use crate::error::Result;
use crate::flow_spec::{FlowFile, FlowModuleSource, FlowUses};
use crate::module_spec::{ModuleFile, ModuleRunner};
use anyhow::Context;
use colored::Colorize;
use std::fs;
use std::path::{Path, PathBuf};

// Standard filenames - can be made configurable in the future
const PROJECT_YAML_FILE: &str = "module.yaml";
const PIPELINE_YAML_FILE: &str = "flow.yaml";

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
                    "A module.yaml already exists at {}",
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
        let yaml_str = fs::read_to_string(&yaml_path).context("Failed to read module.yaml")?;

        let mut module = ModuleFile::parse_yaml(&yaml_str).context("Invalid module.yaml format")?;
        module.metadata.name = project_name_owned.clone();

        if !email_trimmed.is_empty() {
            module.metadata.authors = Some(vec![email_trimmed.clone()]);
            module.metadata.author = None;
        }

        if module
            .metadata
            .version
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .is_none()
        {
            module.metadata.version = Some("0.1.0".to_string());
        }

        let template_default = example
            .clone()
            .unwrap_or_else(|| "dynamic-nextflow".to_string());

        let runner = module.spec.runner.get_or_insert(ModuleRunner {
            kind: Some("nextflow".to_string()),
            entrypoint: Some("workflow.nf".to_string()),
            template: Some(template_default.clone()),
        });

        if runner
            .entrypoint
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .is_none()
        {
            runner.entrypoint = Some("workflow.nf".to_string());
        }

        if runner
            .template
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .is_none()
        {
            runner.template = Some(template_default.clone());
        }

        let author_field = module
            .metadata
            .authors
            .as_ref()
            .and_then(|authors| authors.first().cloned())
            .or_else(|| module.metadata.author.clone())
            .unwrap_or_else(|| "unknown".to_string());

        let workflow_field = runner
            .entrypoint
            .clone()
            .unwrap_or_else(|| "workflow.nf".to_string());
        let template_field = runner
            .template
            .clone()
            .unwrap_or_else(|| "default".to_string());
        let version_field = module
            .metadata
            .version
            .clone()
            .unwrap_or_else(|| "0.1.0".to_string());

        let updated_yaml =
            serde_yaml::to_string(&module).context("Failed to serialize module.yaml")?;
        fs::write(&yaml_path, updated_yaml).context("Failed to update module.yaml")?;

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
        println!("📄 Downloading module.yaml...");
    }

    // Download module.yaml
    let yaml_content = download_file(&raw_url).await?;
    let yaml_str = String::from_utf8(yaml_content).context("Invalid UTF-8 in module.yaml")?;

    let module_info = module_info_from_str(&yaml_str).context("Failed to parse module.yaml")?;

    let project_name = name_override.unwrap_or_else(|| module_info.name.clone());

    // Create project directory path (including version to allow multiple versions)
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;
    let dir_name = format!("{}-{}", project_name, module_info.version);
    let project_dir = projects_dir.join(&dir_name);

    // Check if project already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", project_name, module_info.version);
        if let Some(mut existing) = db.get_project(&identifier)? {
            // Project with this version exists in DB - check if files are valid
            let existing_path = PathBuf::from(&existing.project_path);
            let existing_yaml = existing_path.join(PROJECT_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid project exists - check if path needs updating to new convention
                let expected_dir_name = format!("{}-{}", project_name, module_info.version);
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
                        &module_info.version,
                        &module_info.author,
                        &module_info.workflow,
                        &module_info.template,
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
                        module_info.version.dimmed(),
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
                        module_info.version.dimmed()
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

    // Save module.yaml
    let yaml_path = project_dir.join(PROJECT_YAML_FILE);
    fs::write(&yaml_path, yaml_str)?;

    if !quiet {
        println!("✓ Saved module.yaml");
    }

    // Download workflow file
    let base_url = raw_url
        .rsplit_once('/')
        .map(|(base, _)| base)
        .ok_or_else(|| anyhow::anyhow!("Invalid URL format"))?;

    let workflow_url = format!("{}/{}", base_url, module_info.workflow);

    if !quiet {
        println!("📄 Downloading {}...", module_info.workflow);
    }

    let workflow_content = download_file(&workflow_url).await?;
    let workflow_path = project_dir.join(&module_info.workflow);
    if let Some(parent) = workflow_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&workflow_path, workflow_content)?;

    if !quiet {
        println!("✓ Saved {}", module_info.workflow);
    }

    // Download assets
    if !module_info.assets.is_empty() {
        let assets_dir = project_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Downloading {} assets...", module_info.assets.len());
        }

        for asset in &module_info.assets {
            let asset_rel = asset.strip_prefix("assets/").unwrap_or(asset);
            let asset_url = format!("{}/assets/{}", base_url, asset_rel);

            if !quiet {
                println!("  - {}", asset);
            }

            let asset_content = download_file(&asset_url).await?;
            let asset_path = assets_dir.join(asset_rel);

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
        let identifier = format!("{}@{}", project_name, module_info.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_dir,
            )?;
        } else {
            db.register_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_dir,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &module_info.version,
            &module_info.author,
            &module_info.workflow,
            &module_info.template,
            &project_dir,
        )?;
    }

    let project = db
        .get_project(&project_name)?
        .ok_or_else(|| anyhow::anyhow!("Project '{}' not found after import", project_name))?;

    Ok(project)
}

struct ModuleInfo {
    name: String,
    version: String,
    author: String,
    workflow: String,
    template: String,
    assets: Vec<String>,
}

fn module_info_from_str(raw: &str) -> Result<ModuleInfo> {
    let module = ModuleFile::parse_yaml(raw)?;
    let spec = module.to_project_spec()?;
    Ok(ModuleInfo {
        name: spec.name,
        version: spec.version.unwrap_or_else(|| "0.1.0".to_string()),
        author: spec.author,
        workflow: spec.workflow,
        template: spec.template.unwrap_or_else(|| "default".to_string()),
        assets: spec.assets,
    })
}

enum ModuleImportSource {
    Local { module_dir: PathBuf },
    Url { url: String },
}

impl ModuleImportSource {
    fn cache_key(&self) -> String {
        match self {
            ModuleImportSource::Local { module_dir } => module_dir
                .canonicalize()
                .unwrap_or_else(|_| module_dir.clone())
                .to_string_lossy()
                .to_string(),
            ModuleImportSource::Url { url } => url.clone(),
        }
    }
}

fn resolve_module_source(
    source: &FlowModuleSource,
    dependency_context: &DependencyContext,
    module_hint: &str,
    step_id: &str,
    db: &BioVaultDb,
    quiet: bool,
) -> Result<Option<ModuleImportSource>> {
    use colored::Colorize;

    if let Some(url) = source.url.as_deref() {
        return Ok(Some(ModuleImportSource::Url {
            url: ensure_module_yaml_url(url),
        }));
    }

    if let Some(path) = source.path.as_deref() {
        return match dependency_context {
            DependencyContext::GitHub { base_url } => {
                let trimmed = path.trim_start_matches('/');
                let url = ensure_module_yaml_url(&format!("{}/{}", base_url, trimmed));
                Ok(Some(ModuleImportSource::Url { url }))
            }
            DependencyContext::Local { base_path } => {
                let path_buf = if Path::new(path).is_absolute() {
                    if !quiet {
                        println!(
                            "   {} Step '{}' uses absolute path (keeping as-is)",
                            "ℹ️ ".cyan(),
                            step_id.dimmed()
                        );
                    }
                    return Ok(None);
                } else {
                    base_path.join(path)
                };

                let module_yaml_path = if path_buf.is_dir() {
                    path_buf.join(PROJECT_YAML_FILE)
                } else if path_buf
                    .extension()
                    .map(|e| e == "yaml" || e == "yml")
                    .unwrap_or(false)
                {
                    path_buf.clone()
                } else {
                    path_buf.join(PROJECT_YAML_FILE)
                };

                if module_yaml_path.exists() {
                    let module_dir = if path_buf.is_dir() {
                        path_buf
                    } else {
                        module_yaml_path.parent().unwrap_or(base_path).to_path_buf()
                    };
                    return Ok(Some(ModuleImportSource::Local { module_dir }));
                }

                if db.get_project(module_hint).ok().flatten().is_some() {
                    if !quiet {
                        println!(
                            "   {} Step '{}' uses registered name '{}' (keeping as-is)",
                            "ℹ️ ".cyan(),
                            step_id.dimmed(),
                            module_hint
                        );
                    }
                    return Ok(None);
                }

                Err(anyhow::anyhow!("Cannot resolve module source for '{}'", module_hint).into())
            }
        };
    }

    if let Some(reference) = source.reference.as_deref() {
        if reference.starts_with("http://") || reference.starts_with("https://") {
            return Ok(Some(ModuleImportSource::Url {
                url: ensure_module_yaml_url(reference),
            }));
        }
    }

    Ok(None)
}

fn ensure_module_yaml_url(raw: &str) -> String {
    let trimmed = raw.trim_end_matches('/');
    if trimmed.ends_with(".yaml") || trimmed.ends_with(".yml") {
        trimmed.to_string()
    } else {
        format!("{}/module.yaml", trimmed)
    }
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

    // Look for module.yaml
    let yaml_path = source_path.join(PROJECT_YAML_FILE);
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No module.yaml found in directory: {}",
            source_path.display()
        )
        .into());
    }

    let yaml_str = fs::read_to_string(&yaml_path)?;
    let module_info = module_info_from_str(&yaml_str)?;

    let project_name = module_info.name.clone();

    // Create project directory in managed location (like import_from_url does)
    let biovault_home = config::get_biovault_home()?;
    let projects_dir = biovault_home.join("projects");
    fs::create_dir_all(&projects_dir)?;
    let dir_name = format!("{}-{}", project_name, module_info.version);
    let project_dir = projects_dir.join(&dir_name);

    // Check if project already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", project_name, module_info.version);
        if let Some(existing) = db.get_project(&identifier)? {
            let existing_path = PathBuf::from(&existing.project_path);
            let existing_yaml = existing_path.join(PROJECT_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid project exists - reuse it (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Project '{}' version {} already imported, reusing (id: {})",
                        project_name, module_info.version, existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Project DB record exists but files are missing/invalid - clean it up
                if !quiet {
                    println!(
                        "   🧹 Cleaning up orphaned project record for '{}@{}'",
                        project_name, module_info.version
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

    // Copy module.yaml
    let dest_yaml_path = project_dir.join(PROJECT_YAML_FILE);
    fs::copy(&yaml_path, &dest_yaml_path)?;

    if !quiet {
        println!("✓ Copied module.yaml");
    }

    // Copy workflow file
    let workflow_source = source_path.join(&module_info.workflow);
    if workflow_source.exists() {
        let workflow_dest = project_dir.join(&module_info.workflow);
        if let Some(parent) = workflow_dest.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::copy(&workflow_source, &workflow_dest)?;
        if !quiet {
            println!("✓ Copied {}", module_info.workflow);
        }
    } else if !quiet {
        println!(
            "⚠️  Warning: workflow file '{}' not found in source",
            module_info.workflow
        );
    }

    // Copy assets
    if !module_info.assets.is_empty() {
        let assets_dir = project_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Copying {} assets...", module_info.assets.len());
        }

        for asset in &module_info.assets {
            let asset_rel = asset.strip_prefix("assets/").unwrap_or(asset);
            let asset_source = source_path.join("assets").join(asset_rel);
            if asset_source.exists() {
                let asset_dest = assets_dir.join(asset_rel);

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
        let identifier = format!("{}@{}", project_name, module_info.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_dir,
            )?;
        } else {
            db.register_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_dir,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &module_info.version,
            &module_info.author,
            &module_info.workflow,
            &module_info.template,
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

    // Look for module.yaml
    let yaml_path = project_path.join(PROJECT_YAML_FILE);
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No module.yaml found in directory: {}",
            project_path.display()
        )
        .into());
    }

    if !quiet {
        println!("📁 Registering local project: {}", path.cyan());
    }

    let yaml_str = fs::read_to_string(&yaml_path)?;
    let module_info = module_info_from_str(&yaml_str)?;

    let project_name = name_override.unwrap_or_else(|| module_info.name.clone());

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
        let identifier = format!("{}@{}", project_name, module_info.version);
        if db.get_project(&identifier)?.is_some() {
            db.update_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_path,
            )?;
        } else {
            db.register_project(
                &project_name,
                &module_info.version,
                &module_info.author,
                &module_info.workflow,
                &module_info.template,
                &project_path,
            )?;
        }
    } else {
        db.register_project(
            &project_name,
            &module_info.version,
            &module_info.author,
            &module_info.workflow,
            &module_info.template,
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

    // Show module.yaml contents
    let yaml_path = Path::new(&project.project_path).join(PROJECT_YAML_FILE);
    if yaml_path.exists() {
        println!("📄 module.yaml:");
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

/// Context for resolving flow step dependencies
#[derive(Debug, Clone)]
pub enum DependencyContext {
    /// Dependencies are resolved relative to a GitHub base URL
    GitHub { base_url: String },
    /// Dependencies are resolved relative to a local filesystem path
    Local { base_path: PathBuf },
}

/// Resolve and import all flow step dependencies
///
/// This function handles importing module dependencies for flow steps,
/// rewriting the flow spec to use registered module names.
///
/// Returns `true` if any dependencies were imported and the spec was updated.
pub async fn resolve_pipeline_dependencies(
    flow: &mut FlowFile,
    dependency_context: &DependencyContext,
    flow_yaml_path: &Path,
    overwrite: bool,
    quiet: bool,
) -> Result<bool> {
    use colored::Colorize;

    let db = BioVaultDb::new()?;
    let mut any_rewritten = false;

    if flow.spec.steps.is_empty() {
        return Ok(false);
    }

    if !quiet {
        println!(
            "\n{} Importing dependencies ({} steps):",
            "📦".cyan(),
            flow.spec.steps.len()
        );
    }

    let mut step_updates: Vec<(usize, String, Option<String>)> = Vec::new();
    let mut imported_by_key: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    let mut imported_by_source: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();

    for (index, step) in flow.spec.steps.iter().enumerate() {
        let (module_key, module_ref) = match &step.uses {
            FlowUses::Name(name) => {
                let module_ref = flow
                    .spec
                    .modules
                    .as_ref()
                    .and_then(|modules| modules.get(name))
                    .cloned();
                if module_ref.is_none() {
                    if !quiet {
                        println!(
                            "   {} Step '{}' uses registered name '{}' (keeping as-is)",
                            "ℹ️ ".cyan(),
                            step.id.dimmed(),
                            name
                        );
                    }
                    continue;
                }
                (Some(name.clone()), module_ref)
            }
            FlowUses::Ref(module_ref) => (None, Some(module_ref.clone())),
        };

        let module_ref = match module_ref {
            Some(module_ref) => module_ref,
            None => continue,
        };

        let source = match module_ref.source.as_ref() {
            Some(source) => source,
            None => {
                if !quiet {
                    println!(
                        "   {} Step '{}' uses registered name (keeping as-is)",
                        "ℹ️ ".cyan(),
                        step.id.dimmed()
                    );
                }
                continue;
            }
        };

        if let Some(key) = module_key.as_ref() {
            if let Some(imported_name) = imported_by_key.get(key) {
                step_updates.push((index, imported_name.clone(), module_key.clone()));
                any_rewritten = true;
                continue;
            }
        }

        let resolved_source = match resolve_module_source(
            source,
            dependency_context,
            module_key.as_deref().unwrap_or(step.id.as_str()),
            &step.id,
            &db,
            quiet,
        ) {
            Ok(source) => source,
            Err(e) => {
                if !quiet {
                    println!("{}: {}", "⚠️  failed".yellow(), e.to_string().dimmed());
                }
                continue;
            }
        };

        let resolved_source = match resolved_source {
            Some(source) => source,
            None => continue,
        };

        let source_key = resolved_source.cache_key();
        if let Some(imported_name) = imported_by_source.get(&source_key) {
            step_updates.push((index, imported_name.clone(), module_key.clone()));
            any_rewritten = true;
            if let Some(key) = module_key.as_ref() {
                imported_by_key.insert(key.clone(), imported_name.clone());
            }
            continue;
        }

        if !quiet {
            print!("   {} {} ", "•".cyan(), step.id.bold());
        }

        let import_result = match resolved_source {
            ModuleImportSource::Local { module_dir } => {
                copy_local_project_to_managed(&db, &module_dir, overwrite, true).await
            }
            ModuleImportSource::Url { url } => {
                import_from_url(&db, &url, None, overwrite, true).await
            }
        };

        match import_result {
            Ok(project) => {
                if !quiet {
                    println!("{} → {}", "✓ imported".green(), project.name.green());
                }
                step_updates.push((index, project.name.clone(), module_key.clone()));
                any_rewritten = true;
                imported_by_source.insert(source_key, project.name.clone());
                if let Some(key) = module_key.as_ref() {
                    imported_by_key.insert(key.clone(), project.name.clone());
                }
            }
            Err(e) => {
                if !quiet {
                    println!("{}: {}", "⚠️  failed".yellow(), e.to_string().dimmed());
                }
            }
        }
    }

    for (index, project_name, _) in &step_updates {
        if let Some(updated_step) = flow.spec.steps.get_mut(*index) {
            updated_step.uses = FlowUses::Name(project_name.clone());
        }
    }

    if let Some(modules) = flow.spec.modules.as_mut() {
        for (_, _, module_key) in step_updates {
            if let Some(key) = module_key {
                if let Some(module_ref) = modules.get_mut(&key) {
                    module_ref.source = None;
                }
            }
        }
    }

    flow.save(flow_yaml_path)
        .context("Failed to save flow.yaml with description")?;

    if any_rewritten && !quiet {
        println!(
            "\n{} Updated flow.yaml to use registered names...",
            "🔧".cyan()
        );
    }

    Ok(any_rewritten)
}

/// Import a flow from URL with all its step dependencies
pub async fn import_pipeline_with_deps(
    url: &str,
    name_override: Option<String>,
    overwrite: bool,
) -> Result<String> {
    use colored::Colorize;

    // Convert GitHub URL to raw URL
    let raw_url = url
        .replace("github.com", "raw.githubusercontent.com")
        .replace("/blob/", "/");

    println!("{} Downloading pipeline from {}", "📥".cyan(), url.cyan());

    // Download pipeline YAML
    let yaml_content = download_file(&raw_url).await?;
    let yaml_str = String::from_utf8(yaml_content).context("Invalid UTF-8 in flow.yaml")?;

    let mut flow = FlowFile::parse_yaml(&yaml_str).context("Failed to parse flow.yaml")?;
    if flow.kind != "Flow" {
        return Err(anyhow::anyhow!("Expected Flow kind but found '{}'", flow.kind).into());
    }

    let pipeline_name = name_override.unwrap_or_else(|| flow.metadata.name.clone());

    println!("{} Pipeline: {}", "📦".cyan(), pipeline_name.bold());
    println!("   Steps: {}", flow.spec.steps.len());

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
    let base_url = if let Some((base, _)) = raw_url.rsplit_once('/') {
        base.to_string()
    } else {
        raw_url.clone()
    };

    // Import each step's project and rewrite YAML to use registered names
    resolve_pipeline_dependencies(
        &mut flow,
        &DependencyContext::GitHub { base_url },
        &pipeline_yaml_path,
        overwrite,
        false,
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
apiVersion: syftbox.openmined.org/v1alpha1
kind: Module
metadata:
  name: test-project
  version: 0.1.0
  authors:
    - test@example.com
spec:
  runner:
    kind: nextflow
    entrypoint: workflow.nf
    template: default
  assets: []
  inputs: []
  outputs: []
  parameters: []
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
