use biovault::syftbox::app::DEFAULT_RPC_PERMISSION_CONTENT;

#[test]
fn test_rpc_constant_formatting() {
    println!("Constant content:");
    println!("{}", DEFAULT_RPC_PERMISSION_CONTENT);

    // Check that the constant has proper indentation
    assert!(
        DEFAULT_RPC_PERMISSION_CONTENT.contains("\n  - pattern:"),
        "Constant should have 2-space indentation before dash"
    );

    // Write to temp file and check
    let temp_dir = tempfile::tempdir().unwrap();
    let temp_file = temp_dir.path().join("test.yaml");
    std::fs::write(&temp_file, DEFAULT_RPC_PERMISSION_CONTENT).unwrap();

    let content_back = std::fs::read_to_string(&temp_file).unwrap();
    println!("File content after write:");
    println!("{}", content_back);

    assert!(
        content_back.contains("\n  - pattern:"),
        "File should have 2-space indentation after write"
    );
}
