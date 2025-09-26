#![cfg(feature = "e2e-tests")]

use std::env;
use std::fs;
use std::path::PathBuf;
use std::thread;
use std::time::Duration;

#[test]
#[ignore] // Run with: cargo nextest run --ignored
fn test_syftbox_file_sync() {
    // Get configuration from environment variables
    let client1_email =
        env::var("SYFTBOX_CLIENT1_EMAIL").unwrap_or_else(|_| "client1@syftbox.net".to_string());
    let client2_email =
        env::var("SYFTBOX_CLIENT2_EMAIL").unwrap_or_else(|_| "client2@syftbox.net".to_string());
    let test_clients_dir =
        env::var("TEST_CLIENTS_DIR").unwrap_or_else(|_| "./test-clients-docker".to_string());
    let server_url =
        env::var("SYFTBOX_SERVER_URL").unwrap_or_else(|_| "http://localhost:8080".to_string());
    let test_mode = env::var("TEST_MODE").unwrap_or_else(|_| "docker".to_string());

    println!("Testing file sync between clients:");
    println!("  Mode: {}", test_mode);
    println!("  Client1: {}", client1_email);
    println!("  Client2: {}", client2_email);
    println!("  Server: {}", server_url);
    println!("  Test dir: {}", test_clients_dir);

    // Step 1: Verify server is reachable
    let server_check = reqwest::blocking::get(&server_url);
    assert!(
        server_check.is_ok(),
        "Server is not reachable at {}",
        server_url
    );
    println!("✓ Server is reachable");

    // Step 2: Wait for clients to initialize (give them time to create directories)
    println!("⏳ Waiting for clients to initialize...");
    thread::sleep(Duration::from_secs(10));

    // Step 3: Create test file in client1's public directory
    // Docker mode: test-clients-docker/email/SyftBox/datasites/email/public
    // Local mode: test-clients-local/email/datasites/email/public
    let client1_public_dir = if test_mode == "local" {
        PathBuf::from(&test_clients_dir)
            .join(&client1_email)
            .join("datasites")
            .join(&client1_email)
            .join("public")
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client1_email)
            .join("SyftBox")
            .join("datasites")
            .join(&client1_email)
            .join("public")
    };

    // Create directory if it doesn't exist
    fs::create_dir_all(&client1_public_dir).expect("Failed to create client1 public directory");

    // Use a unique filename for each test run to avoid conflicts
    let test_id = format!(
        "{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs()
    );
    let test_file_name = format!("test_sync_{}.txt", test_id);
    let test_file_content = format!("Integration test file created at timestamp {}", test_id);
    let test_file_path = client1_public_dir.join(&test_file_name);

    fs::write(&test_file_path, &test_file_content).expect("Failed to write test file");
    println!("✓ Test file written to: {:?}", test_file_path);

    // Step 4: Wait for file to sync to client2
    let expected_sync_path = if test_mode == "local" {
        PathBuf::from(&test_clients_dir)
            .join(&client2_email)
            .join("datasites")
            .join(&client1_email)
            .join("public")
            .join(&test_file_name)
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client2_email)
            .join("SyftBox")
            .join("datasites")
            .join(&client1_email)
            .join("public")
            .join(&test_file_name)
    };

    println!("⏳ Waiting for file to sync to: {:?}", expected_sync_path);

    let mut synced = false;
    for i in 0..60 {
        if expected_sync_path.exists() {
            synced = true;
            break;
        }
        thread::sleep(Duration::from_secs(1));
        if i % 10 == 0 && i > 0 {
            println!("   Still waiting... ({}s)", i);
        }
    }

    assert!(synced, "File did not sync within 60 seconds");
    println!("✓ File synced successfully!");

    // Step 5: Verify file content
    let synced_content =
        fs::read_to_string(&expected_sync_path).expect("Failed to read synced file");
    assert_eq!(synced_content, test_file_content, "File content mismatch");
    println!("✓ File content verified");

    println!("\n✅ File sync test completed successfully!");
}

#[test]
#[ignore]
fn test_syftbox_multiple_files_sync() {
    let client1_email =
        env::var("SYFTBOX_CLIENT1_EMAIL").unwrap_or_else(|_| "client1@syftbox.net".to_string());
    let client2_email =
        env::var("SYFTBOX_CLIENT2_EMAIL").unwrap_or_else(|_| "client2@syftbox.net".to_string());
    let test_clients_dir =
        env::var("TEST_CLIENTS_DIR").unwrap_or_else(|_| "./test-clients-docker".to_string());
    let test_mode = env::var("TEST_MODE").unwrap_or_else(|_| "docker".to_string());

    // Wait for initialization
    thread::sleep(Duration::from_secs(5));

    // Create multiple test files
    let client1_public_dir = if test_mode == "local" {
        PathBuf::from(&test_clients_dir)
            .join(&client1_email)
            .join("datasites")
            .join(&client1_email)
            .join("public")
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client1_email)
            .join("SyftBox")
            .join("datasites")
            .join(&client1_email)
            .join("public")
    };

    fs::create_dir_all(&client1_public_dir).expect("Failed to create directory");

    let test_files = vec![
        ("file1.txt", "Content of file 1"),
        ("file2.txt", "Content of file 2"),
        ("data.json", r#"{"test": "data", "number": 42}"#),
    ];

    // Write all test files
    for (filename, content) in &test_files {
        let file_path = client1_public_dir.join(filename);
        fs::write(&file_path, content).unwrap_or_else(|_| panic!("Failed to write {}", filename));
        println!("✓ Created {}", filename);
    }

    // Wait for all files to sync
    let client2_sync_dir = if test_mode == "local" {
        PathBuf::from(&test_clients_dir)
            .join(&client2_email)
            .join("datasites")
            .join(&client1_email)
            .join("public")
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client2_email)
            .join("SyftBox")
            .join("datasites")
            .join(&client1_email)
            .join("public")
    };

    let mut all_synced = true;
    for (filename, expected_content) in &test_files {
        let sync_path = client2_sync_dir.join(filename);

        // Wait for file to appear
        let mut file_synced = false;
        for _ in 0..30 {
            if sync_path.exists() {
                file_synced = true;
                break;
            }
            thread::sleep(Duration::from_secs(1));
        }

        if !file_synced {
            println!("✗ {} did not sync", filename);
            all_synced = false;
            continue;
        }

        // Verify content
        let actual_content = fs::read_to_string(&sync_path)
            .unwrap_or_else(|_| panic!("Failed to read {}", filename));
        if actual_content != *expected_content {
            println!("✗ {} content mismatch", filename);
            all_synced = false;
        } else {
            println!("✓ {} synced correctly", filename);
        }
    }

    assert!(all_synced, "Not all files synced correctly");
    println!("\n✅ Multiple file sync test completed successfully!");
}
