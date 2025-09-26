#[cfg(test)]
mod hard_reset_tests {
    use anyhow::Result;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    // Helper struct to manage test environment
    struct TestEnvironment {
        _temp_dir: TempDir, // Keep alive for RAII cleanup
        biovault_home: PathBuf,
        syftbox_data_dir: PathBuf,
        email: String,
    }

    impl TestEnvironment {
        fn new() -> Result<Self> {
            let temp_dir = TempDir::new()?;
            let biovault_home = temp_dir.path().join(".biovault");
            let syftbox_data_dir = temp_dir.path().join("syftbox_data");
            let email = "test@example.com".to_string();

            Ok(TestEnvironment {
                _temp_dir: temp_dir,
                biovault_home,
                syftbox_data_dir,
                email,
            })
        }

        fn create_biovault_structure(&self) -> Result<()> {
            // Create .biovault directory structure
            fs::create_dir_all(&self.biovault_home)?;
            fs::write(
                self.biovault_home.join("config.yaml"),
                format!("email: {}\nversion: 0.1.0", self.email),
            )?;

            // Create datasite structure
            let datasite_path = self.syftbox_data_dir.join("datasites").join(&self.email);

            // Public biovault
            let public_biovault = datasite_path.join("public").join("biovault");
            fs::create_dir_all(&public_biovault)?;
            fs::write(
                public_biovault.join("participants.yaml"),
                "participants: []",
            )?;

            // Shared submissions
            let shared_submissions = datasite_path
                .join("shared")
                .join("biovault")
                .join("submissions");
            fs::create_dir_all(&shared_submissions)?;
            fs::write(
                shared_submissions.join("test_submission.yaml"),
                "submission: test",
            )?;

            // App data biovault (includes RPC)
            let app_data_biovault = datasite_path.join("app_data").join("biovault");
            fs::create_dir_all(app_data_biovault.join("rpc"))?;
            fs::write(app_data_biovault.join("rpc").join("messages.json"), "[]")?;

            // Private app_data/biovault (at DATA_DIR root, not under datasite!)
            let private_biovault = self
                .syftbox_data_dir
                .join("private")
                .join("app_data")
                .join("biovault");
            fs::create_dir_all(&private_biovault)?;
            fs::write(private_biovault.join("private_data.yaml"), "data: private")?;

            Ok(())
        }

        fn verify_all_deleted(&self) -> bool {
            // Check that all BioVault directories are gone
            let datasite_path = self.syftbox_data_dir.join("datasites").join(&self.email);

            !self.biovault_home.exists()
                && !datasite_path.join("public").join("biovault").exists()
                && !datasite_path
                    .join("shared")
                    .join("biovault")
                    .join("submissions")
                    .exists()
                && !datasite_path.join("app_data").join("biovault").exists()
                && !self
                    .syftbox_data_dir
                    .join("private")
                    .join("app_data")
                    .join("biovault")
                    .exists()
        }

        fn verify_some_exist(&self) -> bool {
            // Check that at least some paths still exist
            let datasite_path = self.syftbox_data_dir.join("datasites").join(&self.email);

            self.biovault_home.exists()
                || datasite_path.join("public").join("biovault").exists()
                || datasite_path.join("shared").join("biovault").exists()
                || datasite_path.join("app_data").join("biovault").exists()
                || self
                    .syftbox_data_dir
                    .join("private")
                    .join("app_data")
                    .join("biovault")
                    .exists()
        }

        fn get_cleanup_paths(&self) -> Vec<PathBuf> {
            let datasite_path = self.syftbox_data_dir.join("datasites").join(&self.email);

            vec![
                self.biovault_home.clone(),
                datasite_path.join("public").join("biovault"),
                datasite_path
                    .join("shared")
                    .join("biovault")
                    .join("submissions"),
                datasite_path.join("app_data").join("biovault"),
                self.syftbox_data_dir
                    .join("private")
                    .join("app_data")
                    .join("biovault"),
            ]
        }
    }

    // Safety check to prevent accidental deletion of real ~/.biovault
    fn is_safe_test_path(path: &Path) -> bool {
        // Path must be in a temp directory or contain "test" or "tmp"
        let path_str = path.to_string_lossy().to_lowercase();

        // Handle both forward and backslashes (Windows uses backslashes)
        let normalized = path_str.replace('\\', "/");

        normalized.contains("/tmp/")
            || normalized.contains("/temp/")
            || normalized.contains("test")
            || normalized.contains("tempfile")
            || normalized.contains("/var/folders/") // macOS temp dirs
            || normalized.contains("appdata/local/temp") // Windows temp dirs
            || path_str.contains(".tmp") // Rust tempfile crate pattern
    }

    fn delete_path_safely(path: &Path) -> Result<()> {
        // CRITICAL SAFETY CHECK: Never delete real user data
        if !is_safe_test_path(path) {
            panic!(
                "SAFETY: Refusing to delete path that doesn't appear to be in a test directory: {}",
                path.display()
            );
        }

        if path.exists() {
            if path.is_dir() {
                fs::remove_dir_all(path)?;
            } else {
                fs::remove_file(path)?;
            }
        }
        Ok(())
    }

    #[test]
    fn test_hard_reset_deletes_all_paths() -> Result<()> {
        let env = TestEnvironment::new()?;

        // Set environment variables for the test
        std::env::set_var("BIOVAULT_HOME", &env.biovault_home);
        std::env::set_var("SYFTBOX_DATA_DIR", &env.syftbox_data_dir);
        std::env::set_var("SYFTBOX_EMAIL", &env.email);

        // Create the full BioVault structure
        env.create_biovault_structure()?;

        // Verify all paths exist before deletion
        assert!(env.biovault_home.exists(), "BioVault home should exist");
        assert!(
            env.verify_some_exist(),
            "At least some BioVault paths should exist"
        );

        // Get paths and delete them
        let paths = env.get_cleanup_paths();
        for path in &paths {
            if path.exists() {
                delete_path_safely(path)?;
            }
        }

        // Verify all paths are deleted
        assert!(
            env.verify_all_deleted(),
            "All BioVault paths should be deleted"
        );

        // Clean up environment variables
        std::env::remove_var("BIOVAULT_HOME");
        std::env::remove_var("SYFTBOX_DATA_DIR");
        std::env::remove_var("SYFTBOX_EMAIL");

        Ok(())
    }

    #[test]
    fn test_hard_reset_handles_non_existent_paths() -> Result<()> {
        let env = TestEnvironment::new()?;

        // Set environment variables for the test
        std::env::set_var("BIOVAULT_HOME", &env.biovault_home);
        std::env::set_var("SYFTBOX_DATA_DIR", &env.syftbox_data_dir);
        std::env::set_var("SYFTBOX_EMAIL", &env.email);

        // Don't create any structure - paths don't exist

        // Get paths and attempt to delete them (should not error)
        let paths = env.get_cleanup_paths();
        for path in &paths {
            // This should not fail even if paths don't exist
            delete_path_safely(path)?;
        }

        // Verify nothing exists (as expected)
        assert!(
            env.verify_all_deleted(),
            "All paths should remain non-existent"
        );

        // Clean up environment variables
        std::env::remove_var("BIOVAULT_HOME");
        std::env::remove_var("SYFTBOX_DATA_DIR");
        std::env::remove_var("SYFTBOX_EMAIL");

        Ok(())
    }

    #[test]
    fn test_hard_reset_partial_structure() -> Result<()> {
        let env = TestEnvironment::new()?;

        // Set environment variables for the test
        std::env::set_var("BIOVAULT_HOME", &env.biovault_home);
        std::env::set_var("SYFTBOX_DATA_DIR", &env.syftbox_data_dir);
        std::env::set_var("SYFTBOX_EMAIL", &env.email);

        // Create only partial structure
        fs::create_dir_all(&env.biovault_home)?;
        fs::write(
            env.biovault_home.join("config.yaml"),
            format!("email: {}", env.email),
        )?;

        let datasite_path = env.syftbox_data_dir.join("datasites").join(&env.email);
        let public_biovault = datasite_path.join("public").join("biovault");
        fs::create_dir_all(&public_biovault)?;

        // Verify partial structure exists
        assert!(env.biovault_home.exists());
        assert!(public_biovault.exists());

        // Get paths and delete them
        let paths = env.get_cleanup_paths();
        for path in &paths {
            if path.exists() {
                delete_path_safely(path)?;
            }
        }

        // Verify all paths are deleted
        assert!(
            env.verify_all_deleted(),
            "All existing paths should be deleted"
        );

        // Clean up environment variables
        std::env::remove_var("BIOVAULT_HOME");
        std::env::remove_var("SYFTBOX_DATA_DIR");
        std::env::remove_var("SYFTBOX_EMAIL");

        Ok(())
    }

    #[test]
    fn test_safety_check_prevents_real_deletion() {
        // Test that safety check prevents deletion of real paths
        let home_dir = dirs::home_dir().unwrap();
        let real_biovault = home_dir.join(".biovault");

        assert!(
            !is_safe_test_path(&real_biovault),
            "Should reject real ~/.biovault path"
        );

        assert!(
            !is_safe_test_path(&PathBuf::from("/usr/local/bin")),
            "Should reject system paths"
        );

        // Test that temp paths are accepted
        assert!(
            is_safe_test_path(&PathBuf::from("/tmp/test_dir")),
            "Should accept /tmp paths"
        );

        assert!(
            is_safe_test_path(&PathBuf::from("/var/folders/xx/test")),
            "Should accept macOS temp paths"
        );

        assert!(
            is_safe_test_path(&PathBuf::from("/home/user/test_biovault")),
            "Should accept paths with 'test' in them"
        );
    }

    #[test]
    fn test_hard_reset_with_edge_case_syftbox_dir() -> Result<()> {
        let env = TestEnvironment::new()?;

        // Test edge case: SYFTBOX_DATA_DIR points directly to a datasite
        let datasite_dir = env.syftbox_data_dir.join("datasites").join(&env.email);

        std::env::set_var("BIOVAULT_HOME", &env.biovault_home);
        std::env::set_var("SYFTBOX_DATA_DIR", &datasite_dir);
        std::env::set_var("SYFTBOX_EMAIL", &env.email);

        // Create structure with this edge-case configuration
        fs::create_dir_all(&env.biovault_home)?;

        // Create public biovault under the datasite
        let public_biovault = datasite_dir.join("public").join("biovault");
        fs::create_dir_all(&public_biovault)?;

        // Private should still be at the real data_dir root
        let private_biovault = env
            .syftbox_data_dir
            .join("private")
            .join("app_data")
            .join("biovault");
        fs::create_dir_all(&private_biovault)?;

        // Delete paths
        let paths_to_delete = vec![
            env.biovault_home.clone(),
            public_biovault.clone(),
            private_biovault.clone(),
        ];

        for path in &paths_to_delete {
            if path.exists() {
                delete_path_safely(path)?;
            }
        }

        // Verify deletion
        assert!(!env.biovault_home.exists());
        assert!(!public_biovault.exists());
        assert!(!private_biovault.exists());

        // Clean up environment variables
        std::env::remove_var("BIOVAULT_HOME");
        std::env::remove_var("SYFTBOX_DATA_DIR");
        std::env::remove_var("SYFTBOX_EMAIL");

        Ok(())
    }

    #[test]
    #[should_panic(expected = "SAFETY: Refusing to delete")]
    fn test_panic_on_unsafe_path() {
        let home_dir = dirs::home_dir().unwrap();
        let real_biovault = home_dir.join(".biovault");

        // This should panic due to safety check
        let _ = delete_path_safely(&real_biovault);
    }

    #[test]
    fn test_read_only_file_deletion() -> Result<()> {
        let env = TestEnvironment::new()?;

        std::env::set_var("BIOVAULT_HOME", &env.biovault_home);

        // Create a read-only file
        fs::create_dir_all(&env.biovault_home)?;
        let readonly_file = env.biovault_home.join("readonly.txt");
        fs::write(&readonly_file, "read only content")?;

        // Make file read-only on Unix systems
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = fs::metadata(&readonly_file)?.permissions();
            perms.set_mode(0o444); // read-only
            fs::set_permissions(&readonly_file, perms)?;
        }

        // Try to delete - should handle gracefully
        let result = delete_path_safely(&env.biovault_home);

        // On some systems this might fail, on others it might succeed
        // The important thing is it doesn't panic
        if let Err(err) = result {
            // If it failed, make sure it's a permission error
            let err_string = err.to_string();
            assert!(
                err_string.contains("permission") || err_string.contains("Permission"),
                "Expected permission error, got: {}",
                err_string
            );
        }

        // Clean up - restore permissions and delete
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            if readonly_file.exists() {
                let mut perms = fs::metadata(&readonly_file)?.permissions();
                perms.set_mode(0o644); // restore write permission
                fs::set_permissions(&readonly_file, perms)?;
            }
        }

        if env.biovault_home.exists() {
            fs::remove_dir_all(&env.biovault_home)?;
        }

        std::env::remove_var("BIOVAULT_HOME");

        Ok(())
    }
}
