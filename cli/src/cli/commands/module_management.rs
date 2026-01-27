use crate::cli::examples;
use crate::config;
use crate::data::{BioVaultDb, Module, ModuleYaml};
use crate::error::Result;
use crate::flow_spec::FlowSpec;
use anyhow::Context;
use colored::Colorize;
use std::fs;
use std::path::{Path, PathBuf};

// Standard filenames - can be made configurable in the future
const MODULE_YAML_FILE: &str = "module.yaml";
const FLOW_YAML_FILE: &str = "flow.yaml";

/// Import a module from URL or register a local module directory
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

    let module = if is_url {
        import_from_url(&db, &source, name.clone(), overwrite, quiet).await?
    } else {
        import_from_local(&db, &source, name.clone(), overwrite, quiet)?
    };

    if quiet {
        let response = crate::data::CliResponse::new(module);
        println!("{}", serde_json::to_string_pretty(&response)?);
    } else if is_url {
        println!(
            "\n✅ Module '{}' imported successfully!",
            module.name.green()
        );
        println!("   Location: {}", module.module_path);
    } else {
        println!(
            "\n✅ Module '{}' registered successfully!",
            module.name.green()
        );
        println!("   Location: {}", module.module_path);
        println!("\n💡 This module is referenced in-place.");
        println!("   Any changes you make will be reflected when you run it.");
    }

    Ok(())
}

pub async fn import_module_record(
    source: String,
    name: Option<String>,
    overwrite: bool,
) -> Result<Module> {
    let db = BioVaultDb::new()?;
    let is_url = source.starts_with("http://") || source.starts_with("https://");

    if is_url {
        import_from_url(&db, &source, name, overwrite, true).await
    } else {
        import_from_local(&db, &source, name, overwrite, true)
    }
}

pub fn create_module_record(
    name: String,
    example: Option<String>,
    target_dir: Option<PathBuf>,
) -> Result<Module> {
    let module_name = name.trim();
    if module_name.is_empty() {
        return Err(anyhow::anyhow!("Module name cannot be empty").into());
    }
    if module_name == "."
        || module_name == ".."
        || module_name.contains('/')
        || module_name.contains('\\')
    {
        return Err(anyhow::anyhow!("Module name contains invalid characters").into());
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
    if let Some(existing) = db.get_module(module_name)? {
        return Err(anyhow::anyhow!(
            "Module '{}' already exists (id: {}). Please choose a different name.",
            module_name,
            existing.id
        )
        .into());
    }

    let (module_dir, cleanup_on_error) = if let Some(dir) = target_dir {
        let module_dir = dir;

        if module_dir.exists() {
            if !module_dir.is_dir() {
                return Err(anyhow::anyhow!(
                    "Target path '{}' exists but is not a directory",
                    module_dir.display()
                )
                .into());
            }

            let yaml_path = module_dir.join(MODULE_YAML_FILE);
            if yaml_path.exists() {
                return Err(anyhow::anyhow!(
                    "A module.yaml already exists at {}",
                    yaml_path.display()
                )
                .into());
            }

            (module_dir, false)
        } else {
            fs::create_dir_all(&module_dir)?;
            (module_dir, true)
        }
    } else {
        let biovault_home = config::get_biovault_home()?;
        let modules_dir = biovault_home.join("modules");
        fs::create_dir_all(&modules_dir)?;

        let module_dir = modules_dir.join(module_name);
        if module_dir.exists() {
            return Err(anyhow::anyhow!(
                "A module named '{}' already exists at {}",
                module_name,
                module_dir.display()
            )
            .into());
        }

        fs::create_dir_all(&module_dir)?;
        (module_dir, true)
    };

    let module_dir_for_cleanup = module_dir.clone();
    let module_name_owned = module_name.to_string();

    let build_result: Result<Module> = (|| {
        if let Some(ref example_name) = example {
            examples::write_example_to_directory(example_name, &module_dir).with_context(|| {
                format!(
                    "Failed to scaffold module from example '{}': {}",
                    example_name,
                    module_dir.display()
                )
            })?;
        } else {
            use crate::module_spec::{self, ModuleSpec};

            let author_placeholder = if email_trimmed.is_empty() {
                "user@example.com".to_string()
            } else {
                email_trimmed.clone()
            };

            // Remove the directory we just created so scaffold_from_spec can create it
            if cleanup_on_error {
                fs::remove_dir_all(&module_dir)?;
            }

            let minimal_spec = ModuleSpec {
                name: module_name_owned.clone(),
                author: author_placeholder,
                workflow: "workflow.nf".to_string(),
                description: None,
                runtime: Some("nextflow".to_string()),
                version: Some("1.0.0".to_string()),
                datasites: None,
                env: Default::default(),
                assets: vec![],
                parameters: vec![],
                inputs: vec![],
                outputs: vec![],
                steps: Vec::new(),
                runner: None,
            };

            module_spec::scaffold_from_spec(minimal_spec, &module_dir)
                .context("Failed to scaffold dynamic-nextflow module")?;
        }

        let yaml_path = module_dir.join(MODULE_YAML_FILE);
        let yaml_str = fs::read_to_string(&yaml_path).context("Failed to read module.yaml")?;

        let mut module_file = crate::module_spec::ModuleFile::parse_yaml(&yaml_str)
            .context("Invalid module.yaml format")?;

        let mut author_field = email_trimmed.clone();
        let mut workflow_field = "workflow.nf".to_string();
        let mut template_field = example
            .clone()
            .unwrap_or_else(|| "dynamic-nextflow".to_string());
        let mut version_field = "1.0.0".to_string();

        module_file.metadata.name = module_name_owned.clone();

        if !email_trimmed.is_empty() {
            module_file.metadata.authors = vec![email_trimmed.clone()];
        }
        if let Some(existing) = module_file.metadata.authors.first() {
            if !existing.trim().is_empty() {
                author_field = existing.clone();
            }
        }

        if module_file.metadata.version.trim().is_empty() {
            module_file.metadata.version = version_field.clone();
        } else {
            version_field = module_file.metadata.version.clone();
        }

        let runner = module_file.spec.runner.get_or_insert_with(Default::default);
        match runner.entrypoint.as_deref().map(str::trim) {
            Some(existing) if !existing.is_empty() => {
                workflow_field = existing.to_string();
            }
            _ => {
                runner.entrypoint = Some(workflow_field.clone());
            }
        }
        match runner.runtime.as_deref().map(str::trim) {
            Some(existing) if !existing.is_empty() => {
                template_field = existing.to_string();
            }
            _ => {
                runner.runtime = Some(template_field.clone());
            }
        }

        let updated_yaml =
            serde_yaml::to_string(&module_file).context("Failed to serialize module.yaml")?;
        fs::write(&yaml_path, updated_yaml).context("Failed to update module.yaml")?;

        db.register_module(
            &module_name_owned,
            &version_field,
            &author_field,
            &workflow_field,
            &template_field,
            &module_dir,
        )?;

        let created = db
            .get_module(&module_name_owned)?
            .ok_or_else(|| anyhow::anyhow!("Module registration completed but not found"))?;

        Ok(created)
    })();

    match build_result {
        Ok(module) => Ok(module),
        Err(err) => {
            if cleanup_on_error {
                let _ = fs::remove_dir_all(module_dir_for_cleanup);
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
) -> Result<Module> {
    if !quiet {
        println!("📥 Importing module from: {}", url.cyan());
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

    // Parse module.yaml
    let module_yaml = ModuleYaml::from_str(&yaml_str).context("Failed to parse module.yaml")?;

    let module_name = name_override.unwrap_or(module_yaml.name.clone());

    // Create module directory path (including version to allow multiple versions)
    let biovault_home = config::get_biovault_home()?;
    let modules_dir = biovault_home.join("modules");
    fs::create_dir_all(&modules_dir)?;
    let dir_name = format!("{}-{}", module_name, module_yaml.version);
    let module_dir = modules_dir.join(&dir_name);

    // Check if module already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", module_name, module_yaml.version);
        if let Some(mut existing) = db.get_module(&identifier)? {
            // Module with this version exists in DB - check if files are valid
            let existing_path = PathBuf::from(&existing.module_path);
            let existing_yaml = existing_path.join(MODULE_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid module exists - check if path needs updating to new convention
                let expected_dir_name = format!("{}-{}", module_name, module_yaml.version);
                let expected_path = modules_dir.join(&expected_dir_name);

                if existing_path != expected_path && !existing_path.ends_with(&expected_dir_name) {
                    // Path is in old format (without version) - migrate to new format
                    if !quiet {
                        println!(
                            "   🔄 Migrating module '{}' to versioned path format",
                            module_name.dimmed()
                        );
                    }

                    // Move directory to new location if needed
                    if !expected_path.exists() {
                        fs::rename(&existing_path, &expected_path).with_context(|| {
                            format!(
                                "Failed to migrate module directory from {} to {}",
                                existing_path.display(),
                                expected_path.display()
                            )
                        })?;
                    } else {
                        // New path already exists - remove old one
                        fs::remove_dir_all(&existing_path).ok();
                    }

                    // Update DB with new path
                    db.update_module(
                        &module_name,
                        &module_yaml.version,
                        &module_yaml.author,
                        &module_yaml.workflow,
                        &module_yaml.template,
                        &expected_path,
                    )?;

                    // Update existing struct
                    existing.module_path = expected_path.to_string_lossy().to_string();
                }

                // Reuse existing valid module (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Module '{}' version {} already registered, reusing (id: {})",
                        module_name.dimmed(),
                        module_yaml.version.dimmed(),
                        existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Module DB record exists but files are missing/invalid - clean it up
                if !quiet {
                    println!(
                        "   🧹 Cleaning up orphaned module record for '{}@{}'",
                        module_name.dimmed(),
                        module_yaml.version.dimmed()
                    );
                }
                db.delete_module(&identifier)?;
                // Continue with fresh import below
            }
        }
    }

    // Handle directory conflicts
    if module_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&module_dir)?;
        } else {
            // Directory exists but not in DB (orphaned) - clean it up
            if !quiet {
                println!(
                    "   🧹 Removing orphaned module directory: {}",
                    module_dir.display().to_string().dimmed()
                );
            }
            fs::remove_dir_all(&module_dir)?;
        }
    }

    fs::create_dir_all(&module_dir)?;

    // Save module.yaml
    let yaml_path = module_dir.join(MODULE_YAML_FILE);
    fs::write(&yaml_path, yaml_str)?;

    if !quiet {
        println!("✓ Saved module.yaml");
    }

    // Download workflow file
    let base_url = raw_url
        .rsplit_once('/')
        .map(|(base, _)| base)
        .ok_or_else(|| anyhow::anyhow!("Invalid URL format"))?;

    let workflow_url = format!("{}/{}", base_url, module_yaml.workflow);

    if !quiet {
        println!("📄 Downloading {}...", module_yaml.workflow);
    }

    let workflow_content = download_file(&workflow_url).await?;
    let workflow_path = module_dir.join(&module_yaml.workflow);
    fs::write(&workflow_path, workflow_content)?;

    if !quiet {
        println!("✓ Saved {}", module_yaml.workflow);
    }

    // Download assets
    if !module_yaml.assets.is_empty() {
        let assets_dir = module_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Downloading {} assets...", module_yaml.assets.len());
        }

        let assets_url = format!("{}/assets", base_url);

        for asset in &module_yaml.assets {
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
        let identifier = format!("{}@{}", module_name, module_yaml.version);
        if db.get_module(&identifier)?.is_some() {
            db.update_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_dir,
            )?;
        } else {
            db.register_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_dir,
            )?;
        }
    } else {
        db.register_module(
            &module_name,
            &module_yaml.version,
            &module_yaml.author,
            &module_yaml.workflow,
            &module_yaml.template,
            &module_dir,
        )?;
    }

    let module = db
        .get_module(&module_name)?
        .ok_or_else(|| anyhow::anyhow!("Module '{}' not found after import", module_name))?;

    Ok(module)
}

/// Copy a local module to the managed directory (~/.biovault/modules/)
/// This mirrors the behavior of import_from_url but copies from local filesystem
async fn copy_local_module_to_managed(
    db: &BioVaultDb,
    source_path: &Path,
    overwrite: bool,
    quiet: bool,
) -> Result<Module> {
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
    let yaml_path = source_path.join(MODULE_YAML_FILE);
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No module.yaml found in directory: {}",
            source_path.display()
        )
        .into());
    }

    // Read and parse module.yaml
    let yaml_str = fs::read_to_string(&yaml_path)?;
    let module_yaml = ModuleYaml::from_str(&yaml_str)?;

    let module_name = module_yaml.name.clone();

    // Create module directory in managed location (like import_from_url does)
    let biovault_home = config::get_biovault_home()?;
    let modules_dir = biovault_home.join("modules");
    fs::create_dir_all(&modules_dir)?;
    let dir_name = format!("{}-{}", module_name, module_yaml.version);
    let module_dir = modules_dir.join(&dir_name);

    // Check if module already exists in DB (with this exact version)
    if !overwrite {
        let identifier = format!("{}@{}", module_name, module_yaml.version);
        if let Some(existing) = db.get_module(&identifier)? {
            let existing_path = PathBuf::from(&existing.module_path);
            let existing_yaml = existing_path.join(MODULE_YAML_FILE);

            if existing_yaml.exists() && existing_path.is_dir() {
                // Valid module exists - reuse it (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Module '{}' version {} already imported, reusing (id: {})",
                        module_name, module_yaml.version, existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Module DB record exists but files are missing/invalid - clean it up
                if !quiet {
                    println!(
                        "   🧹 Cleaning up orphaned module record for '{}@{}'",
                        module_name, module_yaml.version
                    );
                }
                db.delete_module(&identifier)?;
            }
        }
    }

    // Handle directory conflicts
    if module_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&module_dir)?;
        } else {
            // Directory exists but not in DB (orphaned) - clean it up
            if !quiet {
                println!(
                    "   🧹 Removing orphaned module directory: {}",
                    module_dir.display()
                );
            }
            fs::remove_dir_all(&module_dir)?;
        }
    }

    fs::create_dir_all(&module_dir)?;

    // Copy module.yaml
    let dest_yaml_path = module_dir.join(MODULE_YAML_FILE);
    fs::copy(&yaml_path, &dest_yaml_path)?;

    if !quiet {
        println!("✓ Copied module.yaml");
    }

    // Copy workflow file
    let workflow_source = source_path.join(&module_yaml.workflow);
    if workflow_source.exists() {
        let workflow_dest = module_dir.join(&module_yaml.workflow);
        if let Some(parent) = workflow_dest.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::copy(&workflow_source, &workflow_dest)?;
        if !quiet {
            println!("✓ Copied {}", module_yaml.workflow);
        }
    } else if !quiet {
        println!(
            "⚠️  Warning: workflow file '{}' not found in source",
            module_yaml.workflow
        );
    }

    // Copy assets
    if !module_yaml.assets.is_empty() {
        let assets_dir = module_dir.join("assets");
        fs::create_dir_all(&assets_dir)?;

        if !quiet {
            println!("📦 Copying {} assets...", module_yaml.assets.len());
        }

        for asset in &module_yaml.assets {
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
        let identifier = format!("{}@{}", module_name, module_yaml.version);
        if db.get_module(&identifier)?.is_some() {
            db.update_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_dir,
            )?;
        } else {
            db.register_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_dir,
            )?;
        }
    } else {
        db.register_module(
            &module_name,
            &module_yaml.version,
            &module_yaml.author,
            &module_yaml.workflow,
            &module_yaml.template,
            &module_dir,
        )?;
    }

    let module = db
        .get_module(&module_name)?
        .ok_or_else(|| anyhow::anyhow!("Module '{}' not found after import", module_name))?;

    Ok(module)
}

fn import_from_local(
    db: &BioVaultDb,
    path: &str,
    name_override: Option<String>,
    overwrite: bool,
    quiet: bool,
) -> Result<Module> {
    let module_path = PathBuf::from(path);
    if !module_path.exists() {
        return Err(anyhow::anyhow!("Path does not exist: {}", path).into());
    }

    if !module_path.is_dir() {
        return Err(anyhow::anyhow!("Path is not a directory: {}", path).into());
    }

    // Look for module.yaml
    let yaml_path = module_path.join(MODULE_YAML_FILE);
    if !yaml_path.exists() {
        return Err(anyhow::anyhow!(
            "No module.yaml found in directory: {}",
            module_path.display()
        )
        .into());
    }

    if !quiet {
        println!("📁 Registering local module: {}", path.cyan());
    }

    // Parse module.yaml
    let yaml_str = fs::read_to_string(&yaml_path)?;
    let module_yaml = ModuleYaml::from_str(&yaml_str)?;

    let module_name = name_override.unwrap_or(module_yaml.name.clone());

    // Check if module already exists in DB
    if !overwrite {
        if let Some(existing) = db.get_module(&module_name)? {
            // Check if existing module points to the same path
            let existing_path = PathBuf::from(&existing.module_path);
            let canonical_existing = existing_path.canonicalize().ok();
            let canonical_new = module_path.canonicalize().ok();

            if canonical_existing == canonical_new && canonical_existing.is_some() {
                // Same module, already registered - reuse it (idempotent)
                if !quiet {
                    println!(
                        "   ℹ️  Module '{}' already registered at this path (id: {})",
                        module_name.dimmed(),
                        existing.id
                    );
                }
                return Ok(existing);
            } else {
                // Different path - this is a conflict
                return Err(anyhow::anyhow!(
                    "Module '{}' already exists (id: {}) at different path: {}. Use --overwrite to replace.",
                    module_name,
                    existing.id,
                    existing.module_path
                )
                .into());
            }
        }
    }

    // Register in database (using original path, not copied)
    if overwrite {
        let identifier = format!("{}@{}", module_name, module_yaml.version);
        if db.get_module(&identifier)?.is_some() {
            db.update_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_path,
            )?;
        } else {
            db.register_module(
                &module_name,
                &module_yaml.version,
                &module_yaml.author,
                &module_yaml.workflow,
                &module_yaml.template,
                &module_path,
            )?;
        }
    } else {
        db.register_module(
            &module_name,
            &module_yaml.version,
            &module_yaml.author,
            &module_yaml.workflow,
            &module_yaml.template,
            &module_path,
        )?;
    }

    let module = db
        .get_module(&module_name)?
        .ok_or_else(|| anyhow::anyhow!("Module '{}' not found after import", module_name))?;

    Ok(module)
}

/// List all modules
pub fn list(format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;
    let modules = db.list_modules()?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(modules);
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    if modules.is_empty() {
        println!("📭 No modules found");
        println!("\n💡 Import a module:");
        println!("   bv module import <url>");
        println!("   bv module import /path/to/module");
        return Ok(());
    }

    println!("\n📦 {} module(s):\n", modules.len());
    println!(
        "{:<4} {:<20} {:<25} {:<15} {:<20}",
        "ID".bold(),
        "Name".bold(),
        "Author".bold(),
        "Workflow".bold(),
        "Created".bold()
    );
    println!("{}", "─".repeat(90));

    for module in modules {
        println!(
            "{:<4} {:<20} {:<25} {:<15} {:<20}",
            module.id,
            truncate(&module.name, 18),
            truncate(&module.author, 23),
            truncate(&module.workflow, 13),
            &module.created_at[..19] // Trim timestamp to date+time
        );
    }

    println!();

    Ok(())
}

/// Show detailed information about a module
pub fn show(identifier: String, format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;
    let module = db
        .get_module(&identifier)?
        .ok_or_else(|| anyhow::anyhow!("Module '{}' not found", identifier))?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(module);
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    let run_count = db.count_module_runs(module.id)?;

    println!("\n{}", "═".repeat(80));
    println!("📦 {}", module.name.bold().cyan());
    println!("{}", "─".repeat(80));
    println!("ID:           {}", module.id);
    println!("Author:       {}", module.author);
    println!("Workflow:     {}", module.workflow);
    println!("Template:     {}", module.template);
    println!("Location:     {}", module.module_path);
    println!("Created:      {}", module.created_at);
    println!("Runs:         {}", run_count);
    println!("{}", "═".repeat(80));
    println!();

    // Show module.yaml contents
    let yaml_path = Path::new(&module.module_path).join(MODULE_YAML_FILE);
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

/// Delete a module
pub fn delete(identifier: String, keep_files: bool, format: Option<String>) -> Result<()> {
    let db = BioVaultDb::new()?;

    // Get module before deleting
    let module = db.delete_module(&identifier)?;

    if format.as_deref() == Some("json") {
        let response = crate::data::CliResponse::new(serde_json::json!({
            "deleted": module,
            "files_kept": keep_files
        }));
        println!("{}", serde_json::to_string_pretty(&response)?);
        return Ok(());
    }

    println!("✅ Deleted module '{}' from database", module.name.green());

    // Delete files if requested
    if !keep_files {
        let module_path = Path::new(&module.module_path);
        if module_path.exists() {
            fs::remove_dir_all(module_path)?;
            println!("✅ Deleted module files: {}", module_path.display());
        }
    } else {
        println!("📁 Module files kept: {}", module.module_path);
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
pub async fn resolve_flow_dependencies(
    spec: &mut FlowSpec,
    dependency_context: &DependencyContext,
    flow_yaml_path: &Path,
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

            // Resolve module reference based on context
            let (should_use_local, local_path_opt, url_opt) =
                if uses.starts_with("http://") || uses.starts_with("https://") {
                    // Absolute URL - use as-is
                    (false, None, Some(format!("{}/module.yaml", uses)))
                } else {
                    // Relative path - resolve based on context
                    match dependency_context {
                        DependencyContext::GitHub { base_url } => (
                            false,
                            None,
                            Some(format!("{}/{}/module.yaml", base_url, uses)),
                        ),
                        DependencyContext::Local { base_path } => {
                            // Resolve exactly like GitHub: base_path + uses = module location
                            // GitHub does: base_url + "/" + uses + "/module.yaml"
                            // Local does: base_path.join(uses) -> check for module.yaml there
                            let module_path = base_path.join(uses);

                            // Try to find module.yaml (same logic as GitHub URL resolution)
                            let module_yaml_path = if module_path.is_dir() {
                                module_path.join(MODULE_YAML_FILE)
                            } else if module_path
                                .extension()
                                .map(|e| e == "yaml" || e == "yml")
                                .unwrap_or(false)
                            {
                                // Direct path to yaml file
                                module_path.clone()
                            } else {
                                // Assume it's a directory path, look for module.yaml inside
                                module_path.join(MODULE_YAML_FILE)
                            };

                            // Check if the resolved path exists (mirrors GitHub's existence check)
                            if module_yaml_path.exists() {
                                // Path exists - import it (like GitHub downloads it)
                                let module_dir = if module_path.is_dir() {
                                    module_path
                                } else {
                                    module_yaml_path.parent().unwrap_or(base_path).to_path_buf()
                                };
                                (true, Some(module_dir), None)
                            } else {
                                // Path doesn't exist - check if it's already a registered module name
                                // If registered, skip this step (like GitHub does for registered names)
                                let db_check = db.get_module(uses).ok().flatten();
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
                                    // Skip this step entirely - it's already using a registered module name
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

            // Import the module (registers in DB automatically)
            // Smart resolution: check if already registered first, then copy if needed
            let import_result = if should_use_local {
                if let Some(local_path) = local_path_opt {
                    // Check if module at this path is already registered
                    let existing_modules = db.list_modules().context("Failed to list modules")?;
                    let already_registered = existing_modules.iter().any(|p| {
                        PathBuf::from(&p.module_path).canonicalize().ok()
                            == local_path.canonicalize().ok()
                    });

                    if already_registered {
                        // Module already registered - parse local path to get module name and reuse
                        let yaml_path = local_path.join(MODULE_YAML_FILE);
                        if yaml_path.exists() {
                            let yaml_str = fs::read_to_string(&yaml_path)
                                .context("Failed to read module.yaml")?;
                            let module_yaml = ModuleYaml::from_str(&yaml_str)
                                .context("Failed to parse module.yaml")?;
                            let identifier =
                                format!("{}@{}", module_yaml.name, module_yaml.version);

                            match db.get_module(&identifier) {
                                Ok(Some(module)) => {
                                    // Reuse existing registered module
                                    Ok((module, true))
                                }
                                _ => {
                                    // Not found by identifier, copy to managed for consistency
                                    copy_local_module_to_managed(&db, &local_path, overwrite, true)
                                        .await
                                        .map(|p| (p, true))
                                }
                            }
                        } else {
                            return Err(anyhow::anyhow!(
                                "module.yaml not found at {}",
                                local_path.display()
                            )
                            .into());
                        }
                    } else {
                        // Not registered - copy to managed directory (like GitHub imports)
                        copy_local_module_to_managed(&db, &local_path, overwrite, true)
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
                Ok((module, _is_local)) => {
                    if !quiet {
                        println!("{} → {}", "✓ imported".green(), module.name.green());
                    }
                    // Store update for later (avoid borrow checker issue)
                    step_updates.push((index, module.name.clone()));
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
    for (index, module_name) in step_updates {
        if let Some(updated_step) = spec.steps.get_mut(index) {
            updated_step.uses = Some(module_name);
        }
    }

    // Always save the spec (with description preserved) after dependency resolution
    // This ensures description is preserved even if dependencies weren't rewritten
    spec.save(flow_yaml_path)
        .context("Failed to save flow.yaml with description")?;

    if any_rewritten && !quiet {
        println!(
            "\n{} Updated flow.yaml to use registered names...",
            "🔧".cyan()
        );
    }

    Ok(any_rewritten)
}

/// Import a flow from URL or local path with all its step dependencies
pub async fn import_flow_with_deps(
    url: &str,
    name_override: Option<String>,
    overwrite: bool,
) -> Result<String> {
    use crate::flow_spec::{FlowFile, FlowSpec};
    use colored::Colorize;

    // Check if this is a local path (starts with / or file://)
    let is_local = url.starts_with('/') || url.starts_with("file://");
    let local_path = if url.starts_with("file://") {
        url.strip_prefix("file://").unwrap_or(url)
    } else {
        url
    };

    let (yaml_str, dependency_context) = if is_local {
        // Local file path
        let path = PathBuf::from(local_path);
        println!(
            "{} Loading flow from {}",
            "📥".cyan(),
            path.display().to_string().cyan()
        );

        if !path.exists() {
            return Err(anyhow::anyhow!("Local flow file not found: {}", path.display()).into());
        }

        let yaml_str = fs::read_to_string(&path).context("Failed to read flow.yaml")?;

        // Extract base path for resolving relative module paths
        let base_path = path
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."));

        (yaml_str, DependencyContext::Local { base_path })
    } else {
        // Remote URL
        // Convert GitHub URL to raw URL
        let raw_url = url
            .replace("github.com", "raw.githubusercontent.com")
            .replace("/blob/", "/");

        println!("{} Downloading flow from {}", "📥".cyan(), url.cyan());

        // Download flow YAML
        let yaml_content = download_file(&raw_url).await?;
        let yaml_str = String::from_utf8(yaml_content).context("Invalid UTF-8 in flow.yaml")?;

        // Extract base URL for resolving relative module paths
        let base_url = if let Some(idx) = url.rfind('/') {
            url[..idx].to_string()
        } else {
            url.to_string()
        };

        (yaml_str, DependencyContext::GitHub { base_url })
    };

    // Parse flow spec via FlowFile to support Flow schema shapes
    let flow_file = FlowFile::parse_yaml(&yaml_str)
        .inspect_err(|e| {
            eprintln!(
                "Failed to parse flow.yaml from {}. First 200 chars:\n{}",
                url,
                truncate(&yaml_str, 200)
            );
        })
        .context("Failed to parse flow.yaml")?;
    let mut spec: FlowSpec = flow_file
        .to_flow_spec()
        .context("Failed to convert flow spec")?;

    let flow_name = name_override.unwrap_or_else(|| spec.name.clone());

    println!("{} Flow: {}", "📦".cyan(), flow_name.bold());
    println!("   Steps: {}", spec.steps.len());

    // Create flow directory
    let biovault_home = crate::config::get_biovault_home()?;
    let flows_dir = biovault_home.join("flows");
    fs::create_dir_all(&flows_dir)?;

    let flow_dir = flows_dir.join(&flow_name);

    if flow_dir.exists() {
        if overwrite {
            fs::remove_dir_all(&flow_dir)?;
        } else {
            return Err(anyhow::anyhow!(
                "Flow directory already exists: {}. Use --overwrite to replace.",
                flow_dir.display()
            )
            .into());
        }
    }

    fs::create_dir_all(&flow_dir)?;

    let flow_yaml_path = flow_dir.join(FLOW_YAML_FILE);

    // Import each step's module and rewrite YAML to use registered names
    resolve_flow_dependencies(
        &mut spec,
        &dependency_context,
        &flow_yaml_path,
        overwrite,
        false, // quiet = false for CLI output
    )
    .await?;

    // Register flow in database (check for existing if overwrite)
    let db = BioVaultDb::new()?;
    let flow_id = if overwrite {
        // Check if flow with this name already exists
        let existing_flows = db.list_flows()?;
        if let Some(existing) = existing_flows.iter().find(|p| p.name == flow_name) {
            // Delete existing flow from DB
            db.delete_flow(existing.id)?;
        }
        // Register the new flow
        db.register_flow(&flow_name, &flow_dir.to_string_lossy())?
    } else {
        db.register_flow(&flow_name, &flow_dir.to_string_lossy())?
    };

    println!(
        "\n{} Flow '{}' imported successfully!",
        "✅".green().bold(),
        flow_name.bold()
    );
    println!("   Location: {}", flow_dir.display().to_string().dimmed());
    println!("   ID: {}", flow_id);

    Ok(flow_dir.to_string_lossy().to_string())
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

        // Create a test module directory
        let module_dir = tmp.path().join("test-module");
        fs::create_dir_all(&module_dir).unwrap();

        let yaml_content = r#"
apiVersion: syftbox.openmined.org/v1alpha1
kind: Module
metadata:
  name: test-module
  version: 0.1.0
  authors:
    - test@example.com
spec:
  runner:
    kind: nextflow
    template: default
    entrypoint: workflow.nf
  inputs: []
  outputs: []
  parameters: []
  assets: []
"#;
        fs::write(module_dir.join(MODULE_YAML_FILE), yaml_content).unwrap();
        fs::write(module_dir.join("workflow.nf"), "// workflow").unwrap();

        // Import the module
        let module = import_from_local(
            &BioVaultDb::new().unwrap(),
            module_dir.to_str().unwrap(),
            None,
            false,
            true,
        )
        .expect("local module import should succeed");
        assert_eq!(module.name, "test-module");

        // Verify it's in the database
        let db = BioVaultDb::new().unwrap();
        let module = db.get_module("test-module").unwrap();
        assert!(module.is_some());

        teardown_test();
    }

    #[test]
    fn test_list_modules() {
        let tmp = setup_test();
        let db = BioVaultDb::new().unwrap();

        // Create some test modules
        let module_path = tmp.path().join("module1");
        fs::create_dir_all(&module_path).unwrap();
        db.register_module(
            "module1",
            "1.0.0",
            "author@example.com",
            "workflow.nf",
            "default",
            &module_path,
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
