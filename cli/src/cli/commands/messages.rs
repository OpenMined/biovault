use crate::config::Config;
use crate::syftbox::rpc::{check_requests, process_request, send_response};
use crate::syftbox::{RpcRequest, RpcResponse, SyftBoxApp};
use anyhow::Result;
use serde::{Deserialize, Serialize};

const MESSAGE_ENDPOINT: &str = "/message";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MessagePayload {
    pub message: String,
    pub from: Option<String>,
    pub timestamp: Option<String>,
}

/// Initialize the message endpoint for BioVault
pub fn init_message_endpoint(config: &Config) -> Result<SyftBoxApp> {
    let data_dir = config.get_syftbox_data_dir()?;
    let app = SyftBoxApp::new(&data_dir, &config.email, "biovault")?;

    // Register the message endpoint
    app.register_endpoint(MESSAGE_ENDPOINT)?;

    println!("BioVault RPC initialized for {}", config.email);
    println!(
        "Message endpoint registered at: {}",
        app.build_syft_url(MESSAGE_ENDPOINT)
    );

    Ok(app)
}

/// Check for incoming messages
pub fn check_messages(config: &Config) -> Result<Vec<(String, MessagePayload)>> {
    let app = init_message_endpoint(config)?;
    let requests = check_requests(&app, MESSAGE_ENDPOINT)?;

    let mut messages = Vec::new();

    for (request_path, request) in requests {
        // Extract sender from request
        let sender = request.sender.clone();

        // Try to parse the message payload
        match request.body_as_json::<MessagePayload>() {
            Ok(payload) => {
                messages.push((sender.clone(), payload));

                // Send acknowledgment response
                let response = RpcResponse::ok_json(
                    &request,
                    config.email.clone(),
                    &serde_json::json!({
                        "status": "received",
                        "message": "Message received successfully"
                    }),
                )?;

                send_response(&app, MESSAGE_ENDPOINT, &request_path, &response)?;
            }
            Err(e) => {
                eprintln!("Failed to parse message from {}: {}", sender, e);

                // Send error response
                let error_response = RpcResponse::error(
                    &request,
                    config.email.clone(),
                    400,
                    "Invalid message format",
                );

                send_response(&app, MESSAGE_ENDPOINT, &request_path, &error_response)?;
            }
        }
    }

    if !messages.is_empty() {
        println!("Received {} new message(s)", messages.len());
    }

    Ok(messages)
}

/// Send a message to another datasite
pub fn send_message(config: &Config, recipient_email: &str, message: &str) -> Result<()> {
    let data_dir = config.get_syftbox_data_dir()?;

    // Build path to recipient's BioVault RPC directory
    let recipient_app_dir = data_dir
        .join("datasites")
        .join(recipient_email)
        .join("app_data")
        .join("biovault");

    let recipient_rpc_dir = recipient_app_dir.join("rpc");

    // Check if recipient has BioVault app directory
    if !recipient_app_dir.exists() {
        return Err(anyhow::anyhow!(
            "Recipient {} has not installed BioVault. Directory not found: {:?}",
            recipient_email,
            recipient_app_dir
        ));
    }

    // Create app-level permission file if it doesn't exist
    let app_permission_file = recipient_app_dir.join("syft.pub.yaml");
    if !app_permission_file.exists() {
        std::fs::write(
            &app_permission_file,
            crate::syftbox::app::DEFAULT_APP_PERMISSION_CONTENT,
        )?;
        println!("Created app permission file: {:?}", app_permission_file);
    }

    // Create RPC directory if it doesn't exist
    if !recipient_rpc_dir.exists() {
        std::fs::create_dir_all(&recipient_rpc_dir)?;
        println!(
            "Created RPC directory for recipient: {:?}",
            recipient_rpc_dir
        );

        // Create permission file for RPC directory
        let rpc_permission_file = recipient_rpc_dir.join("syft.pub.yaml");
        if !rpc_permission_file.exists() {
            std::fs::write(
                &rpc_permission_file,
                crate::syftbox::app::DEFAULT_RPC_PERMISSION_CONTENT,
            )?;
            println!("Created RPC permission file: {:?}", rpc_permission_file);
        }
    }

    // Create the message endpoint directory if it doesn't exist
    let recipient_endpoint_dir = recipient_rpc_dir.join("message");
    if !recipient_endpoint_dir.exists() {
        std::fs::create_dir_all(&recipient_endpoint_dir)?;
        println!(
            "Created message endpoint for recipient: {:?}",
            recipient_endpoint_dir
        );
    }

    // Create the message payload
    let payload = MessagePayload {
        message: message.to_string(),
        from: Some(config.email.clone()),
        timestamp: Some(chrono::Utc::now().to_rfc3339()),
    };

    let payload_json = serde_json::to_vec(&payload)?;

    // Build the recipient's syft URL
    let recipient_url = format!("syft://{}/app_data/biovault/rpc/message", recipient_email);

    // Create the request
    let request = RpcRequest::new(
        config.email.clone(),
        recipient_url,
        "POST".to_string(),
        payload_json,
    );

    // Write the request directly to the recipient's endpoint
    let request_filename = format!("{}.request", request.id);
    let request_path = recipient_endpoint_dir.join(request_filename);

    let request_json = serde_json::to_string_pretty(&request)?;
    std::fs::write(&request_path, request_json)?;

    println!(
        "Message sent to {} (request ID: {})",
        recipient_email, request.id
    );
    println!("Request file: {:?}", request_path);

    Ok(())
}

/// Process incoming messages with a custom handler
pub fn process_messages<F>(config: &Config, handler: F) -> Result<()>
where
    F: Fn(&str, &MessagePayload) -> Result<String>,
{
    let app = init_message_endpoint(config)?;
    let requests = check_requests(&app, MESSAGE_ENDPOINT)?;

    for (request_path, request) in requests {
        let sender = request.sender.clone();

        process_request(
            &app,
            MESSAGE_ENDPOINT,
            &request_path,
            &request,
            |req| match req.body_as_json::<MessagePayload>() {
                Ok(payload) => match handler(&sender, &payload) {
                    Ok(response_message) => RpcResponse::ok_json(
                        req,
                        config.email.clone(),
                        &serde_json::json!({
                            "status": "processed",
                            "response": response_message
                        }),
                    ),
                    Err(e) => Ok(RpcResponse::error(
                        req,
                        config.email.clone(),
                        500,
                        &format!("Error processing message: {}", e),
                    )),
                },
                Err(e) => Ok(RpcResponse::error(
                    req,
                    config.email.clone(),
                    400,
                    &format!("Invalid message format: {}", e),
                )),
            },
        )?;
    }

    Ok(())
}

/// List all messages (requests and responses) in the message endpoint
pub fn list_messages(config: &Config) -> Result<()> {
    let app = init_message_endpoint(config)?;

    // Check for incoming requests
    let requests = check_requests(&app, MESSAGE_ENDPOINT)?;

    if !requests.is_empty() {
        println!("\n📥 Incoming Messages (Requests):");
        for (_, request) in requests {
            println!("  ID: {}", request.id);
            println!("  From: {}", request.sender);
            println!("  Created: {}", request.created);

            if let Ok(payload) = request.body_as_json::<MessagePayload>() {
                println!("  Message: {}", payload.message);
            }
            println!();
        }
    } else {
        println!("No incoming messages");
    }

    // Check for responses to messages we sent
    let responses = crate::syftbox::rpc::check_responses(&app, MESSAGE_ENDPOINT)?;

    if !responses.is_empty() {
        println!("\n📤 Message Responses:");
        for (_, response) in responses {
            println!("  ID: {}", response.id);
            println!("  From: {}", response.sender);
            println!("  Status: {}", response.status_code);
            println!("  Created: {}", response.created);

            if let Ok(body) = response.body_as_string() {
                println!("  Response: {}", body);
            }
            println!();
        }
    } else {
        println!("No message responses");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn create_test_config() -> Config {
        Config {
            email: "test@example.com".to_string(),
            syftbox_config: None,
        }
    }

    #[test]
    fn test_init_message_endpoint() -> Result<()> {
        let temp_dir = TempDir::new()?;
        std::env::set_var("SYFTBOX_DATA_DIR", temp_dir.path());
        let config = create_test_config();

        let app = init_message_endpoint(&config)?;
        assert!(app.endpoint_exists(MESSAGE_ENDPOINT));

        Ok(())
    }

    #[test]
    fn test_send_and_check_messages() -> Result<()> {
        let temp_dir = TempDir::new()?;
        std::env::set_var("SYFTBOX_DATA_DIR", temp_dir.path());
        let config = create_test_config();

        // Create recipient's endpoint directory first
        let recipient_rpc_dir = temp_dir
            .path()
            .join("datasites")
            .join("recipient@example.com")
            .join("app_data")
            .join("biovault")
            .join("rpc")
            .join("message");
        std::fs::create_dir_all(&recipient_rpc_dir)?;

        // Send a message
        send_message(&config, "recipient@example.com", "Hello from test")?;

        // The message should appear in the recipient's endpoint
        let request_files: Vec<_> = std::fs::read_dir(&recipient_rpc_dir)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry
                    .path()
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .map(|ext| ext == "request")
                    .unwrap_or(false)
            })
            .collect();

        assert_eq!(request_files.len(), 1);

        Ok(())
    }
}
