use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};

const PERMISSION_FILE_NAME: &str = "syft.pub.yaml";
pub const DEFAULT_RPC_PERMISSION_CONTENT: &str = r#"rules:
- pattern: '**/*.request'
  access:
    admin: []
    read:
    - '*'
    write:
    - '*'
- pattern: '**/*.response'
  access:
    admin: []
    read:
    - '*'
    write:
    - '*'
"#;

pub const DEFAULT_APP_PERMISSION_CONTENT: &str = r#"rules:
- pattern: "**"
  access:
    admin: []
    read:
    - "*"
    write: []
"#;

/// Represents a SyftBox application
#[derive(Debug, Clone)]
pub struct SyftBoxApp {
    pub app_name: String,
    pub email: String,
    pub data_dir: PathBuf,
    pub app_data_dir: PathBuf,
    pub rpc_dir: PathBuf,
}

impl SyftBoxApp {
    /// Create a new SyftBox app, ensuring all necessary directories and files exist
    pub fn new(data_dir: &Path, email: &str, app_name: &str) -> Result<Self> {
        let app_data_dir = data_dir
            .join("datasites")
            .join(email)
            .join("app_data")
            .join(app_name);

        let rpc_dir = app_data_dir.join("rpc");

        let app = Self {
            app_name: app_name.to_string(),
            email: email.to_string(),
            data_dir: data_dir.to_path_buf(),
            app_data_dir,
            rpc_dir,
        };

        app.ensure_initialized()?;
        Ok(app)
    }

    /// Ensure the app directory structure and permission files exist
    fn ensure_initialized(&self) -> Result<()> {
        // Create app_data directory if it doesn't exist
        if !self.app_data_dir.exists() {
            fs::create_dir_all(&self.app_data_dir).with_context(|| {
                format!(
                    "Failed to create app_data directory: {:?}",
                    self.app_data_dir
                )
            })?;
            // quiet: avoid noisy output in normal operations
        }

        // Create app-level permission file if it doesn't exist
        let app_permission_file = self.app_data_dir.join(PERMISSION_FILE_NAME);
        if !app_permission_file.exists() {
            fs::write(&app_permission_file, DEFAULT_APP_PERMISSION_CONTENT).with_context(|| {
                format!(
                    "Failed to create app permission file: {:?}",
                    app_permission_file
                )
            })?;
            // quiet
        }

        // Create RPC directory if it doesn't exist
        if !self.rpc_dir.exists() {
            fs::create_dir_all(&self.rpc_dir)
                .with_context(|| format!("Failed to create RPC directory: {:?}", self.rpc_dir))?;
            // quiet
        }

        // Create RPC permission file if it doesn't exist
        let rpc_permission_file = self.rpc_dir.join(PERMISSION_FILE_NAME);
        if !rpc_permission_file.exists() {
            fs::write(&rpc_permission_file, DEFAULT_RPC_PERMISSION_CONTENT).with_context(|| {
                format!(
                    "Failed to create RPC permission file: {:?}",
                    rpc_permission_file
                )
            })?;
            // quiet
        }

        Ok(())
    }

    /// Get the path for a specific endpoint
    pub fn endpoint_path(&self, endpoint_name: &str) -> PathBuf {
        // Remove leading slash if present
        let clean_name = endpoint_name.trim_start_matches('/');
        self.rpc_dir.join(clean_name)
    }

    /// Register a new endpoint by creating its directory
    pub fn register_endpoint(&self, endpoint_name: &str) -> Result<PathBuf> {
        let endpoint_dir = self.endpoint_path(endpoint_name);

        if !endpoint_dir.exists() {
            fs::create_dir_all(&endpoint_dir).with_context(|| {
                format!("Failed to create endpoint directory: {:?}", endpoint_dir)
            })?;
            // quiet
        }

        Ok(endpoint_dir)
    }

    /// Check if an endpoint exists
    pub fn endpoint_exists(&self, endpoint_name: &str) -> bool {
        self.endpoint_path(endpoint_name).exists()
    }

    /// List all registered endpoints
    pub fn list_endpoints(&self) -> Result<Vec<String>> {
        let mut endpoints = Vec::new();

        if self.rpc_dir.exists() {
            for entry in fs::read_dir(&self.rpc_dir)? {
                let entry = entry?;
                let path = entry.path();

                // Skip the permission file
                if path.is_dir() {
                    if let Some(name) = path.file_name() {
                        if let Some(name_str) = name.to_str() {
                            endpoints.push(format!("/{}", name_str));
                        }
                    }
                }
            }
        }

        Ok(endpoints)
    }

    /// Build a syft:// URL for an endpoint
    pub fn build_syft_url(&self, endpoint_name: &str) -> String {
        let clean_endpoint = endpoint_name.trim_start_matches('/');
        format!(
            "syft://{}/app_data/{}/rpc/{}",
            self.email, self.app_name, clean_endpoint
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_app_initialization() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;

        // Check that directories were created
        assert!(app.rpc_dir.exists());

        // Check that permission file was created
        let permission_file = app.rpc_dir.join(PERMISSION_FILE_NAME);
        assert!(permission_file.exists());

        Ok(())
    }

    #[test]
    fn test_endpoint_registration() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app")?;

        // Register an endpoint
        let endpoint_path = app.register_endpoint("/message")?;
        assert!(endpoint_path.exists());
        assert!(app.endpoint_exists("/message"));

        // List endpoints
        let endpoints = app.list_endpoints()?;
        assert!(endpoints.contains(&"/message".to_string()));

        Ok(())
    }

    #[test]
    fn test_syft_url_building() {
        let temp_dir = TempDir::new().unwrap();
        let app = SyftBoxApp::new(temp_dir.path(), "test@example.com", "test_app").unwrap();

        let url = app.build_syft_url("/message");
        assert_eq!(url, "syft://test@example.com/app_data/test_app/rpc/message");
    }
}
