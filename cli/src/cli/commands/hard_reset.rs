use crate::config::{get_biovault_home, Config};
use crate::Result;
use anyhow::Context;
use colored::*;
use dialoguer::{theme::ColorfulTheme, Confirm};
use std::fs;
use std::path::PathBuf;
use tracing::{error, info, warn};

#[derive(Debug)]
pub struct CleanupPath {
    pub description: String,
    pub path: PathBuf,
    pub exists: bool,
}

fn get_cleanup_paths(config: &Config) -> Vec<CleanupPath> {
    let mut paths = Vec::new();

    // 1. .biovault folder
    let biovault_home = get_biovault_home().unwrap_or_else(|_| PathBuf::from("~/.biovault"));
    paths.push(CleanupPath {
        description: "BioVault home directory".to_string(),
        path: biovault_home.clone(),
        exists: biovault_home.exists(),
    });

    // Get SyftBox data directory
    let data_dir = if let Ok(data_dir) = config.get_syftbox_data_dir() {
        Some(data_dir)
    } else if let Ok(data_dir) = std::env::var("SYFTBOX_DATA_DIR") {
        Some(PathBuf::from(data_dir))
    } else {
        None
    };

    if let Some(data_dir) = data_dir {
        // Determine the actual root data_dir and datasite path
        let (real_data_dir, datasite_path) =
            if data_dir.components().any(|c| c.as_os_str() == "datasites")
                && data_dir.to_string_lossy().contains(&config.email)
            {
                // SYFTBOX_DATA_DIR is pointing to the datasite itself
                // We need to find the parent that doesn't contain "datasites"
                let mut parent = data_dir.clone();
                while parent.components().any(|c| c.as_os_str() == "datasites") {
                    if let Some(p) = parent.parent() {
                        parent = p.to_path_buf();
                    } else {
                        break;
                    }
                }
                (parent, data_dir.clone())
            } else {
                // Normal case: data_dir is the root
                let datasite = data_dir.join("datasites").join(&config.email);
                (data_dir.clone(), datasite)
            };

        // Show the datasite base path for context
        info!("Using datasite path: {}", datasite_path.display());

        // 2. Public biovault folder (under datasite)
        let public_biovault = datasite_path.join("public").join("biovault");
        paths.push(CleanupPath {
            description: "└─ public/biovault".to_string(),
            path: public_biovault.clone(),
            exists: public_biovault.exists(),
        });

        // 3. Shared submissions folder (under datasite)
        let shared_submissions = datasite_path
            .join("shared")
            .join("biovault")
            .join("submissions");
        paths.push(CleanupPath {
            description: "└─ shared/biovault/submissions".to_string(),
            path: shared_submissions.clone(),
            exists: shared_submissions.exists(),
        });

        // 4. App data biovault folder (includes RPC) (under datasite)
        let app_data_biovault = datasite_path.join("app_data").join("biovault");
        paths.push(CleanupPath {
            description: "└─ app_data/biovault (includes RPC)".to_string(),
            path: app_data_biovault.clone(),
            exists: app_data_biovault.exists(),
        });

        // 5. Private app_data/biovault (at real DATA_DIR root, not under datasite!)
        let private_biovault = real_data_dir
            .join("private")
            .join("app_data")
            .join("biovault");
        paths.push(CleanupPath {
            description: "Private app_data/biovault".to_string(),
            path: private_biovault.clone(),
            exists: private_biovault.exists(),
        });
    } else {
        warn!("Could not determine SyftBox data directory - some paths may not be cleaned");
    }

    paths
}

fn delete_path(path: &PathBuf) -> anyhow::Result<()> {
    if path.exists() {
        if path.is_dir() {
            fs::remove_dir_all(path)
                .with_context(|| format!("Failed to remove directory: {}", path.display()))?;
        } else {
            fs::remove_file(path)
                .with_context(|| format!("Failed to remove file: {}", path.display()))?;
        }
    }
    Ok(())
}

pub async fn execute(ignore_warning: bool) -> Result<()> {
    println!("\n{}", "⚠️  BioVault Hard Reset".red().bold());
    println!(
        "{}",
        "This will DELETE all BioVault data and configuration!".red()
    );
    println!();

    // Try to load config - if it doesn't exist, create a minimal one
    let config = match Config::get_config_path().and_then(Config::from_file) {
        Ok(c) => {
            info!("Loaded existing config for email: {}", c.email);
            c
        }
        Err(e) => {
            // Create a minimal config for path resolution
            warn!(
                "No existing config found ({}), using minimal config for cleanup",
                e
            );
            Config {
                email: std::env::var("SYFTBOX_EMAIL")
                    .unwrap_or_else(|_| "unknown@email.com".to_string()),
                syftbox_config: None,
                version: None,
                binary_paths: None,
            }
        }
    };

    let paths = get_cleanup_paths(&config);

    // Show what will be deleted
    println!("{}", "The following paths will be deleted:".yellow().bold());
    println!();

    let mut any_exists = false;

    for (i, path_info) in paths.iter().enumerate() {
        // Add spacing and headers for different sections
        if i == 1 && paths.len() > 1 {
            // Display the datasite header - extract it from the first datasite path
            if let Some(first_datasite_path) =
                paths.iter().find(|p| p.description.starts_with("└─"))
            {
                // Get parent of parent to show the datasite root
                if let Some(parent) = first_datasite_path.path.parent() {
                    if let Some(parent2) = parent.parent() {
                        println!("\nDatasite: {}", parent2.display().to_string().cyan());
                    }
                }
            }
        }

        // Add header for private folder (separate from datasite)
        if path_info.description == "Private app_data/biovault" {
            // Extract the data_dir from the private path
            if let Some(parent) = path_info.path.parent() {
                if let Some(parent2) = parent.parent() {
                    if let Some(parent3) = parent2.parent() {
                        println!(
                            "\nSyftBox Data Dir: {}",
                            parent3.display().to_string().cyan()
                        );
                    }
                }
            }
        }

        if path_info.exists {
            any_exists = true;
            if path_info.description.starts_with("└─") {
                println!("  {} {}", "✓".green(), path_info.description);
            } else if path_info.description == "Private app_data/biovault" {
                println!("  {} └─ private/app_data/biovault", "✓".green());
            } else {
                println!(
                    "  {} {} ({})",
                    "✓".green(),
                    path_info.description,
                    path_info.path.display().to_string().dimmed()
                );
            }
        } else if path_info.description.starts_with("└─") {
            println!(
                "  {} {} {}",
                "○".dimmed(),
                path_info.description.dimmed(),
                "(does not exist)".dimmed()
            );
        } else if path_info.description == "Private app_data/biovault" {
            println!(
                "  {} {} {}",
                "○".dimmed(),
                "└─ private/app_data/biovault".dimmed(),
                "(does not exist)".dimmed()
            );
        } else {
            println!(
                "  {} {} ({}) {}",
                "○".dimmed(),
                path_info.description.dimmed(),
                path_info.path.display().to_string().dimmed(),
                "(does not exist)".dimmed()
            );
        }
    }

    if !any_exists {
        println!();
        println!("{}", "No BioVault data found to delete.".green());
        return Ok(());
    }

    println!();

    // Confirmation
    if !ignore_warning {
        let theme = ColorfulTheme::default();

        // First confirmation
        let confirmed = Confirm::with_theme(&theme)
            .with_prompt("Are you sure you want to delete all BioVault data?")
            .default(false)
            .interact()
            .map_err(|e| anyhow::anyhow!("Failed to get confirmation: {}", e))?;

        if !confirmed {
            println!("{}", "Operation cancelled.".yellow());
            return Ok(());
        }

        // Second confirmation for safety
        let really_confirmed = Confirm::with_theme(&theme)
            .with_prompt("This action CANNOT be undone. Continue?")
            .default(false)
            .interact()
            .map_err(|e| anyhow::anyhow!("Failed to get confirmation: {}", e))?;

        if !really_confirmed {
            println!("{}", "Operation cancelled.".yellow());
            return Ok(());
        }
    } else {
        info!("Running in non-interactive mode (--ignore-warning)");
    }

    // Perform deletion
    println!();
    println!("{}", "Deleting BioVault data...".yellow());

    let mut errors = Vec::new();
    let mut successes = 0;

    for path_info in paths {
        if path_info.exists {
            print!("  Deleting {}... ", path_info.description);
            match delete_path(&path_info.path) {
                Ok(_) => {
                    println!("{}", "✓".green());
                    successes += 1;
                }
                Err(e) => {
                    println!("{}", "✗".red());
                    error!("Failed to delete {}: {}", path_info.path.display(), e);
                    errors.push((path_info.description.clone(), e));
                }
            }
        }
    }

    println!();

    if !errors.is_empty() {
        println!("{}", "⚠️  Some paths could not be deleted:".yellow().bold());
        for (desc, err) in errors {
            println!("  {} {}: {}", "✗".red(), desc, err.to_string().red());
        }
        println!();
        println!(
            "{}",
            "You may need to manually delete these paths with appropriate permissions.".yellow()
        );
    }

    if successes > 0 {
        println!(
            "{} {} paths deleted successfully",
            "✓".green().bold(),
            successes
        );
        println!();
        println!("{}", "BioVault has been reset.".green().bold());
        println!("Run {} to set up a fresh installation.", "`bv init`".cyan());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{
        clear_test_biovault_home, clear_test_syftbox_data_dir, set_test_biovault_home,
        set_test_syftbox_data_dir,
    };
    use std::fs;
    use tempfile::TempDir;

    fn test_config() -> Config {
        Config {
            email: "user@example.com".into(),
            syftbox_config: None,
            version: None,
        }
    }

    #[test]
    fn cleanup_paths_cover_root_data_dir() {
        let tmp = TempDir::new().unwrap();
        let bv_home = tmp.path().join("biovault_home");
        fs::create_dir_all(&bv_home).unwrap();
        set_test_biovault_home(&bv_home);

        let root = tmp.path().join("syftbox");
        let datasite = root.join("datasites").join("user@example.com");
        let public = datasite.join("public").join("biovault");
        let shared = datasite.join("shared").join("biovault").join("submissions");
        let app_data = datasite.join("app_data").join("biovault");
        let private = root.join("private").join("app_data").join("biovault");

        fs::create_dir_all(&root).unwrap();
        fs::create_dir_all(&public).unwrap();
        fs::create_dir_all(&shared).unwrap();
        fs::create_dir_all(&app_data).unwrap();
        fs::create_dir_all(&private).unwrap();
        set_test_syftbox_data_dir(&root);

        let config = test_config();
        let paths = get_cleanup_paths(&config);

        assert_eq!(paths.len(), 5);
        assert!(paths.iter().any(|p| p.path == bv_home && p.exists));
        assert!(paths.iter().any(|p| p.path == public && p.exists));
        assert!(paths.iter().any(|p| p.path == shared && p.exists));
        assert!(paths.iter().any(|p| p.path == app_data && p.exists));
        assert!(paths.iter().any(|p| p.path == private && p.exists));

        clear_test_biovault_home();
        clear_test_syftbox_data_dir();
    }

    #[test]
    fn cleanup_paths_handles_datasite_data_dir() {
        let tmp = TempDir::new().unwrap();
        let bv_home = tmp.path().join("biovault_home");
        fs::create_dir_all(&bv_home).unwrap();
        set_test_biovault_home(&bv_home);

        let root = tmp.path().join("syftbox");
        let datasite = root.join("datasites").join("user@example.com");
        let public = datasite.join("public").join("biovault");
        let shared = datasite.join("shared").join("biovault").join("submissions");
        let app_data = datasite.join("app_data").join("biovault");
        let private = root.join("private").join("app_data").join("biovault");

        fs::create_dir_all(&datasite).unwrap();
        fs::create_dir_all(&public).unwrap();
        fs::create_dir_all(&shared).unwrap();
        fs::create_dir_all(&app_data).unwrap();
        fs::create_dir_all(&private).unwrap();
        set_test_syftbox_data_dir(&datasite);

        let config = test_config();
        let paths = get_cleanup_paths(&config);

        assert_eq!(paths.len(), 5);
        assert!(paths.iter().any(|p| p.path == bv_home && p.exists));
        assert!(paths.iter().any(|p| p.path == public && p.exists));
        assert!(paths.iter().any(|p| p.path == shared && p.exists));
        assert!(paths.iter().any(|p| p.path == app_data && p.exists));
        assert!(paths.iter().any(|p| p.path == private && p.exists));

        clear_test_biovault_home();
        clear_test_syftbox_data_dir();
    }

    #[test]
    fn delete_path_removes_files_and_directories() {
        let tmp = TempDir::new().unwrap();
        let file_path = tmp.path().join("file.txt");
        fs::write(&file_path, "data").unwrap();
        delete_path(&file_path).unwrap();
        assert!(!file_path.exists());

        let dir_path = tmp.path().join("dir");
        fs::create_dir_all(dir_path.join("nested")).unwrap();
        delete_path(&dir_path).unwrap();
        assert!(!dir_path.exists());
    }

    #[test]
    fn execute_removes_all_known_paths_in_ignore_mode() {
        let tmp = TempDir::new().unwrap();
        let bv_home = tmp.path().join("biovault_home");
        fs::create_dir_all(&bv_home).unwrap();
        set_test_biovault_home(&bv_home);

        let config_path = bv_home.join("config.yaml");
        test_config().save(&config_path).unwrap();

        let root = tmp.path().join("syftbox");
        let datasite = root.join("datasites").join("user@example.com");
        let public = datasite.join("public").join("biovault");
        let shared = datasite.join("shared").join("biovault").join("submissions");
        let app_data = datasite.join("app_data").join("biovault");
        let private = root.join("private").join("app_data").join("biovault");

        fs::create_dir_all(&public).unwrap();
        fs::create_dir_all(&shared).unwrap();
        fs::create_dir_all(&app_data).unwrap();
        fs::create_dir_all(&private).unwrap();
        fs::write(public.join("file"), "contents").unwrap();
        fs::write(private.join("file"), "contents").unwrap();
        set_test_syftbox_data_dir(&root);

        let runtime = tokio::runtime::Builder::new_current_thread()
            .enable_all()
            .build()
            .unwrap();
        runtime.block_on(async {
            execute(true).await.unwrap();
        });

        assert!(!bv_home.exists());
        assert!(!public.exists());
        assert!(!shared.exists());
        assert!(!app_data.exists());
        assert!(!private.exists());

        clear_test_biovault_home();
        clear_test_syftbox_data_dir();
    }
}
