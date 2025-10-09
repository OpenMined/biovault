use crate::error::Result;
use anyhow::anyhow;
use std::process::Command;

pub async fn install(version: &str) -> Result<()> {
    println!("📦 Installing Python {} via UV...", version);

    let output = Command::new("uv")
        .args(["python", "install", version])
        .output()?;

    if output.status.success() {
        println!("✅ Python {} installed successfully", version);
        let stdout = String::from_utf8_lossy(&output.stdout);
        if !stdout.trim().is_empty() {
            println!("   {}", stdout.trim());
        }
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow!("Failed to install Python {}: {}", version, stderr).into());
    }

    Ok(())
}

pub async fn list(installed_only: bool) -> Result<()> {
    let mut args = vec!["python", "list"];
    if installed_only {
        args.push("--only-installed");
    }

    let output = Command::new("uv").args(&args).output()?;

    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        println!("{}", stdout);
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow!("Failed to list Python versions: {}", stderr).into());
    }

    Ok(())
}

pub async fn show() -> Result<()> {
    let output = Command::new("uv").args(["python", "find"]).output()?;

    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
        println!("Default Python: {}", stdout);
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        println!("⚠️  Could not find default Python: {}", stderr.trim());
        println!("Try: bv python install 3.12");
    }

    Ok(())
}

pub async fn uninstall(version: &str) -> Result<()> {
    println!("🗑️  Uninstalling Python {}...", version);

    let output = Command::new("uv")
        .args(["python", "uninstall", version])
        .output()?;

    if output.status.success() {
        println!("✅ Python {} uninstalled", version);
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow!("Failed to uninstall Python {}: {}", version, stderr).into());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    #[ignore = "requires UV to be installed"]
    async fn test_show_python() {
        let result = show().await;
        assert!(result.is_ok());
    }

    #[tokio::test]
    #[ignore = "requires UV to be installed"]
    async fn test_list_python() {
        let result = list(false).await;
        assert!(result.is_ok());
    }
}
