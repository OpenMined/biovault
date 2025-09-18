use anyhow::Result;
use reqwest::blocking::Client;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::thread;
use std::time::Duration;

#[test]
#[ignore]
fn test_syftbox_permission_system() -> Result<()> {
    println!("\n🔐 Testing SyftBox permission system:");

    // Get configuration from environment
    let server_url =
        env::var("SYFTBOX_SERVER_URL").unwrap_or_else(|_| "http://localhost:8080".to_string());
    let client1_email =
        env::var("SYFTBOX_CLIENT1_EMAIL").unwrap_or_else(|_| "client1@syftbox.net".to_string());
    let client2_email =
        env::var("SYFTBOX_CLIENT2_EMAIL").unwrap_or_else(|_| "client2@syftbox.net".to_string());
    let bad_client_email =
        env::var("SYFTBOX_BAD_CLIENT_EMAIL").unwrap_or_else(|_| "bad@syftbox.net".to_string());
    let test_clients_dir =
        env::var("TEST_CLIENTS_DIR").unwrap_or_else(|_| "./test-clients".to_string());
    let test_mode = env::var("TEST_MODE").unwrap_or_else(|_| "docker".to_string());

    println!("  Mode: {}", test_mode);
    println!("  Client1: {}", client1_email);
    println!("  Client2: {}", client2_email);
    println!("  Bad Client: {}", bad_client_email);
    println!("  Server: {}", server_url);
    println!("  Test dir: {}", test_clients_dir);

    // Check server is reachable
    let client = Client::new();
    let response = client.get(&server_url).send()?;
    assert!(response.status().is_success(), "Server should be reachable");
    println!("✓ Server is reachable");

    // Wait for clients to initialize
    println!("⏳ Waiting for clients to initialize...");
    thread::sleep(Duration::from_secs(10));

    // Build paths based on test mode
    let (client1_base, client2_base, bad_client_base) = if test_mode == "local" {
        (
            PathBuf::from(&test_clients_dir).join(&client1_email),
            PathBuf::from(&test_clients_dir).join(&client2_email),
            PathBuf::from(&test_clients_dir).join(&bad_client_email),
        )
    } else {
        (
            PathBuf::from(&test_clients_dir)
                .join(&client1_email)
                .join("SyftBox"),
            PathBuf::from(&test_clients_dir)
                .join(&client2_email)
                .join("SyftBox"),
            PathBuf::from(&test_clients_dir)
                .join(&bad_client_email)
                .join("SyftBox"),
        )
    };

    // Test 1: User-only read/write permissions
    println!("\n📝 Test 1: User-only read/write permissions");
    test_user_only_permissions(
        &client1_base,
        &client2_base,
        &bad_client_base,
        &client1_email,
    )?;

    // Test 2: Everyone read, user-only write
    println!("\n📝 Test 2: Everyone read, user-only write");
    test_everyone_read_user_write(
        &client1_base,
        &client2_base,
        &bad_client_base,
        &client1_email,
    )?;

    // Test 3: Specific user permissions
    println!("\n📝 Test 3: Specific user permissions");
    test_specific_user_permissions(
        &client1_base,
        &client2_base,
        &bad_client_base,
        &client1_email,
        &client2_email,
    )?;

    // Test 4 removed for now - template permissions need more clarity

    println!("\n✅ All permission tests completed successfully!");
    Ok(())
}

fn test_user_only_permissions(
    client1_base: &PathBuf,
    client2_base: &PathBuf,
    bad_client_base: &PathBuf,
    client1_email: &str,
) -> Result<()> {
    // Create a private folder in client1's datasite
    let private_dir = client1_base
        .join("datasites")
        .join(client1_email)
        .join("private");
    fs::create_dir_all(&private_dir)?;

    // Create user-specific folder that will match the template
    let user_dir = private_dir.join(client1_email);
    fs::create_dir_all(&user_dir)?;

    // Create permission file with UserEmail template - only users can access their own folder
    let permission_content = r#"rules:
  - pattern: '{{.UserEmail}}/*'
    access:
      admin: []
      read:
        - 'USER'
      write:
        - 'USER'
"#;

    let permission_file = private_dir.join("syft.pub.yaml");
    fs::write(&permission_file, permission_content)?;
    println!("✓ Created user-only permission file with {{.UserEmail}} template");

    // Create a test file in the user-specific folder
    let test_file = user_dir.join("secret.txt");
    fs::write(&test_file, "This is a secret only for the owner")?;
    println!("✓ Created secret file in user-specific folder");

    // Wait for sync
    thread::sleep(Duration::from_secs(5));

    // Check that client2 and bad client don't have access to client1's user folder
    let client2_private = client2_base
        .join("datasites")
        .join(client1_email)
        .join("private")
        .join(client1_email)
        .join("secret.txt");
    let bad_private = bad_client_base
        .join("datasites")
        .join(client1_email)
        .join("private")
        .join(client1_email)
        .join("secret.txt");

    assert!(
        !client2_private.exists(),
        "Client2 should not have access to client1's user-specific file"
    );
    assert!(
        !bad_private.exists(),
        "Bad client should not have access to client1's user-specific file"
    );
    println!("✓ Other clients cannot access user-specific files");

    Ok(())
}

fn test_everyone_read_user_write(
    client1_base: &PathBuf,
    client2_base: &PathBuf,
    bad_client_base: &PathBuf,
    client1_email: &str,
) -> Result<()> {
    // Create a public folder in client1's datasite
    let public_dir = client1_base
        .join("datasites")
        .join(client1_email)
        .join("public");
    fs::create_dir_all(&public_dir)?;

    // Create permission file for everyone read, owner write is implicit
    let permission_content = r#"rules:
  - pattern: '**'
    access:
      admin: []
      read:
        - '*'
      write: []
"#;

    let permission_file = public_dir.join("syft.pub.yaml");
    fs::write(&permission_file, permission_content)?;
    println!("✓ Created everyone-read permission file");

    // Create a test file
    let test_file = public_dir.join("readme.txt");
    let content = "Everyone can read this, but only owner can modify";
    fs::write(&test_file, content)?;
    println!("✓ Created readme file");

    // Wait longer for sync to propagate to all clients
    thread::sleep(Duration::from_secs(15));

    // Check that client2 and bad client have the file (read access)
    let client2_public = client2_base
        .join("datasites")
        .join(client1_email)
        .join("public")
        .join("readme.txt");
    let bad_public = bad_client_base
        .join("datasites")
        .join(client1_email)
        .join("public")
        .join("readme.txt");

    // Wait for file to sync with retries
    let mut synced = false;
    for i in 0..20 {
        if client2_public.exists() && bad_public.exists() {
            synced = true;
            break;
        }
        println!("Waiting for public file sync... attempt {}/20", i + 1);
        thread::sleep(Duration::from_secs(5));
    }

    assert!(synced, "Files should be readable by all clients");

    // Verify content
    let client2_content = fs::read_to_string(&client2_public)?;
    assert_eq!(
        client2_content, content,
        "Client2 should have correct content"
    );

    let bad_content = fs::read_to_string(&bad_public)?;
    assert_eq!(
        bad_content, content,
        "Bad client should have correct content"
    );

    println!("✓ All clients can read public files");

    // Try to write from bad client (should fail or not sync back)
    let bad_write_attempt = bad_public.parent().unwrap().join("hack_attempt.txt");
    fs::write(&bad_write_attempt, "Trying to write where I shouldn't")?;

    // Wait for potential sync
    thread::sleep(Duration::from_secs(5));

    // Check that the file didn't sync to client1
    let client1_hack = client1_base
        .join("datasites")
        .join(client1_email)
        .join("public")
        .join("hack_attempt.txt");
    assert!(
        !client1_hack.exists(),
        "Bad client's write should not sync to owner"
    );
    println!("✓ Non-owner writes are blocked");

    Ok(())
}

fn test_specific_user_permissions(
    client1_base: &PathBuf,
    client2_base: &PathBuf,
    bad_client_base: &PathBuf,
    client1_email: &str,
    client2_email: &str,
) -> Result<()> {
    // Create a shared folder in client1's datasite
    let shared_dir = client1_base
        .join("datasites")
        .join(client1_email)
        .join("shared");
    fs::create_dir_all(&shared_dir)?;

    // Create permission file allowing specific user (client2) to read, owner access is implicit
    let permission_content = format!(
        r#"rules:
  - pattern: '**'
    access:
      admin: []
      read:
        - '{}'
      write: []
"#,
        client2_email
    );

    let permission_file = shared_dir.join("syft.pub.yaml");
    fs::write(&permission_file, permission_content)?;
    println!(
        "✓ Created shared permission file allowing {} to read",
        client2_email
    );

    // Create a test file
    let test_file = shared_dir.join("collaboration.txt");
    fs::write(&test_file, "Client1 owns this, Client2 can read")?;
    println!("✓ Created collaboration file");

    // Wait for sync
    thread::sleep(Duration::from_secs(10));

    // Check that client2 has the file but bad client doesn't
    let client2_shared = client2_base
        .join("datasites")
        .join(client1_email)
        .join("shared")
        .join("collaboration.txt");
    let bad_shared = bad_client_base
        .join("datasites")
        .join(client1_email)
        .join("shared")
        .join("collaboration.txt");

    // Wait for file to sync to client2
    let mut synced = false;
    for _ in 0..12 {
        if client2_shared.exists() {
            synced = true;
            break;
        }
        thread::sleep(Duration::from_secs(5));
    }

    assert!(synced, "Client2 should have access to shared file");
    assert!(
        !bad_shared.exists(),
        "Bad client should not have access to shared file"
    );
    println!("✓ Specific user permissions working correctly");

    Ok(())
}

#[test]
#[ignore]
fn test_permission_updates() -> Result<()> {
    println!("\n🔄 Testing permission updates:");

    let server_url =
        env::var("SYFTBOX_SERVER_URL").unwrap_or_else(|_| "http://localhost:8080".to_string());
    let client1_email =
        env::var("SYFTBOX_CLIENT1_EMAIL").unwrap_or_else(|_| "client1@syftbox.net".to_string());
    let client2_email =
        env::var("SYFTBOX_CLIENT2_EMAIL").unwrap_or_else(|_| "client2@syftbox.net".to_string());
    let test_clients_dir =
        env::var("TEST_CLIENTS_DIR").unwrap_or_else(|_| "./test-clients".to_string());
    let test_mode = env::var("TEST_MODE").unwrap_or_else(|_| "docker".to_string());

    // Build paths
    let client1_base = if test_mode == "local" {
        PathBuf::from(&test_clients_dir).join(&client1_email)
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client1_email)
            .join("SyftBox")
    };

    let client2_base = if test_mode == "local" {
        PathBuf::from(&test_clients_dir).join(&client2_email)
    } else {
        PathBuf::from(&test_clients_dir)
            .join(&client2_email)
            .join("SyftBox")
    };

    // Create a folder that starts private
    let dynamic_dir = client1_base
        .join("datasites")
        .join(&client1_email)
        .join("dynamic");
    fs::create_dir_all(&dynamic_dir)?;

    // Start with private permissions using UserEmail template
    let private_perms = r#"rules:
  - pattern: '{{.UserEmail}}/*'
    access:
      admin: []
      read:
        - 'USER'
      write:
        - 'USER'
"#;

    // Create user-specific folder
    let user_dir = dynamic_dir.join(&client1_email);
    fs::create_dir_all(&user_dir)?;

    fs::write(dynamic_dir.join("syft.pub.yaml"), private_perms)?;
    fs::write(user_dir.join("secret.txt"), "Initially private")?;
    println!("✓ Created private folder with UserEmail template");

    // Wait for sync
    thread::sleep(Duration::from_secs(5));

    // Verify client2 doesn't have access to client1's user folder
    let client2_dynamic = client2_base
        .join("datasites")
        .join(&client1_email)
        .join("dynamic")
        .join(&client1_email)
        .join("secret.txt");
    assert!(
        !client2_dynamic.exists(),
        "File should be private initially"
    );
    println!("✓ File is private");

    // Update permissions to allow everyone to read all files
    let public_perms = r#"rules:
  - pattern: '**/*'
    access:
      admin: []
      read:
        - '*'
      write: []
"#;

    fs::write(dynamic_dir.join("syft.pub.yaml"), public_perms)?;
    println!("✓ Updated permissions to allow everyone to read");

    // Wait for permission update to propagate
    thread::sleep(Duration::from_secs(10));

    // Check if client2 now has access
    let mut updated = false;
    for _ in 0..12 {
        if client2_dynamic.exists() {
            updated = true;
            break;
        }
        thread::sleep(Duration::from_secs(5));
    }

    assert!(
        updated,
        "Client2 should have access after permission update"
    );
    println!("✓ Permission updates propagate correctly");

    Ok(())
}
