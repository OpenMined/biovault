use crate::Result;
use serde::{Deserialize, Serialize};
use std::process::Command;

#[derive(Debug, Serialize, Deserialize)]
struct DependencyConfig {
    dependencies: Vec<Dependency>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Dependency {
    name: String,
    check_running: bool,
    install_instructions: String,
    description: String,
}

pub async fn execute() -> Result<()> {
    // Load the deps.yaml file embedded in the binary
    let deps_yaml = include_str!("../../deps.yaml");
    let config: DependencyConfig = serde_yaml::from_str(deps_yaml)?;
    
    println!("BioVault Dependency Check");
    println!("=========================\n");
    
    let mut all_found = true;
    let mut all_running = true;
    
    for dep in &config.dependencies {
        print!("Checking {}... ", dep.name);
        
        // Check if the binary exists in PATH
        let exists = which::which(&dep.name).is_ok();
        
        if !exists {
            all_found = false;
            println!("❌ NOT FOUND");
            println!("  Description: {}", dep.description);
            println!("  Installation instructions:");
            for line in dep.install_instructions.lines() {
                if !line.trim().is_empty() {
                    println!("    {}", line);
                }
            }
            println!();
        } else {
            print!("✓ Found");
            
            // Check if it needs to be running and if it is
            if dep.check_running {
                let is_running = check_if_running(&dep.name);
                if is_running {
                    println!(" (running)");
                } else {
                    all_running = false;
                    println!(" (NOT RUNNING)");
                    println!("  To start {}, run: {}", dep.name, get_start_command(&dep.name));
                }
            } else {
                println!();
            }
        }
    }
    
    println!("\n=========================");
    if all_found && all_running {
        println!("✓ All dependencies satisfied!");
    } else if !all_found {
        println!("⚠️  Some dependencies are missing. Please install them using the instructions above.");
    } else if !all_running {
        println!("⚠️  Some services are not running. Please start them using the commands above.");
    }
    
    Ok(())
}

fn check_if_running(service: &str) -> bool {
    match service {
        "docker" => {
            // Check if Docker daemon is running
            Command::new("docker")
                .arg("info")
                .output()
                .map(|output| output.status.success())
                .unwrap_or(false)
        }
        _ => false,
    }
}

fn get_start_command(service: &str) -> String {
    match service {
        "docker" => "Open Docker Desktop or run 'sudo dockerd' (Linux)".to_string(),
        _ => format!("Start {}", service),
    }
}