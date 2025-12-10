use crate::config::Config;
use crate::defaults::SYFTBOX_DEFAULT_SERVER_URL;
use crate::error::Error as ConfigError;
use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::env;
use std::path::{Path, PathBuf};

const DEFAULT_SERVER_URL: &str = SYFTBOX_DEFAULT_SERVER_URL;
const DEFAULT_CLIENT_URL: &str = "http://localhost:7938";

#[derive(Debug, Serialize)]
struct RequestOtpPayload {
    email: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    name: Option<String>,
}

#[derive(Debug, Serialize)]
struct VerifyOtpPayload {
    email: String,
    code: String,
}

#[derive(Debug, Deserialize)]
struct VerifyResponse {
    #[serde(rename = "accessToken")]
    access_token: String,
    #[serde(rename = "refreshToken")]
    refresh_token: String,
}

#[derive(Debug, Deserialize)]
struct ApiErrorMessage {
    #[serde(default)]
    message: Option<String>,
    #[serde(default)]
    error: Option<String>,
    #[serde(default)]
    detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
struct SyftboxConfigFile {
    #[serde(default)]
    data_dir: String,
    #[serde(default)]
    email: String,
    #[serde(default)]
    server_url: String,
    #[serde(default)]
    client_url: Option<String>,
    #[serde(default)]
    client_token: Option<String>,
    #[serde(default)]
    refresh_token: Option<String>,
    #[serde(flatten)]
    other: HashMap<String, Value>,
}

impl SyftboxConfigFile {
    fn load(path: &Path) -> Result<Self> {
        if path.exists() {
            let content = std::fs::read_to_string(path)
                .with_context(|| format!("Failed to read SyftBox config: {}", path.display()))?;
            let mut parsed: SyftboxConfigFile =
                serde_json::from_str(&content).with_context(|| {
                    format!("Failed to parse SyftBox config JSON: {}", path.display())
                })?;
            // Normalise base fields to avoid stray whitespace.
            parsed.email = parsed.email.trim().to_string();
            parsed.server_url = parsed.server_url.trim().to_string();
            if parsed.data_dir.is_empty() {
                if let Some(dir) = env::var_os("SYFTBOX_DATA_DIR") {
                    parsed.data_dir = dir.to_string_lossy().to_string();
                }
            }
            Ok(parsed)
        } else {
            Ok(Self::default())
        }
    }

    fn save(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Failed to create SyftBox config directory: {}",
                    parent.display()
                )
            })?;
        }
        let json = serde_json::to_string_pretty(self).with_context(|| {
            format!("Failed to serialise SyftBox config for {}", path.display())
        })?;
        std::fs::write(path, json)
            .with_context(|| format!("Failed to write SyftBox config: {}", path.display()))?;
        Ok(())
    }
}

pub async fn request_otp(
    email: Option<String>,
    name: Option<String>,
    server: Option<String>,
) -> Result<()> {
    let config_state = maybe_load_config()?;
    let config_ref = config_state.as_ref().map(|(cfg, _)| cfg);

    let email = resolve_email(email, config_ref)?;
    let server_url = resolve_server(server.as_deref(), config_ref);
    let request_url = join_url(&server_url, "/auth/otp/request");

    let payload = RequestOtpPayload {
        email: email.clone(),
        name: name.clone(),
    };

    let client = reqwest::Client::new();
    let response = client
        .post(&request_url)
        .header("User-Agent", "biovault-cli")
        .json(&payload)
        .send()
        .await
        .with_context(|| format!("Failed to contact {}", request_url))?;

    if !response.status().is_success() {
        return Err(build_api_error("OTP request", response).await);
    }

    println!("Requested OTP for {} via {}", email, server_url);

    if let Some((mut config, path)) = config_state {
        let mut creds = config.syftbox_credentials.unwrap_or_default();
        let mut changed = false;

        if creds.email.as_deref() != Some(&email) {
            creds.email = Some(email.clone());
            changed = true;
        }
        if creds.server_url.as_deref() != Some(&server_url) {
            creds.server_url = Some(server_url.clone());
            changed = true;
        }

        if changed {
            config.syftbox_credentials = Some(creds);
            config.save(&path)?;
        }
    }

    Ok(())
}

pub async fn submit_otp(
    code: &str,
    email: Option<String>,
    server: Option<String>,
    config_path_override: Option<String>,
    data_dir_override: Option<String>,
    client_url_override: Option<String>,
) -> Result<()> {
    let trimmed_code = code.trim();
    if trimmed_code.is_empty() {
        return Err(anyhow!("OTP code cannot be empty"));
    }

    let mut config_state = maybe_load_config()?;
    let config_ref = config_state.as_ref().map(|(cfg, _)| cfg);

    let email = resolve_email(email, config_ref)?;
    let server_url = resolve_server(server.as_deref(), config_ref);
    let verify_url = join_url(&server_url, "/auth/otp/verify");

    let payload = VerifyOtpPayload {
        email: email.clone(),
        code: trimmed_code.to_string(),
    };

    let client = reqwest::Client::new();
    let response = client
        .post(&verify_url)
        .header("User-Agent", "biovault-cli")
        .json(&payload)
        .send()
        .await
        .with_context(|| format!("Failed to contact {}", verify_url))?;

    if !response.status().is_success() {
        return Err(build_api_error("OTP verification", response).await);
    }

    let tokens: VerifyResponse = response
        .json()
        .await
        .context("Failed to parse OTP verification response")?;

    let (mut config, config_path) = match config_state.take() {
        Some(state) => state,
        None => load_or_create_config(&email)?,
    };

    if config.email.is_empty() {
        config.email = email.clone();
    }

    let syftbox_config_path = if let Some(path) = config_path_override {
        config.syftbox_config = Some(path.clone());
        PathBuf::from(path)
    } else {
        config.get_syftbox_config_path()?
    };

    let mut syftbox_config = SyftboxConfigFile::load(&syftbox_config_path)?;

    let data_dir = resolve_data_dir(data_dir_override.as_deref(), &config, &syftbox_config)?;
    std::fs::create_dir_all(&data_dir).with_context(|| {
        format!(
            "Failed to ensure SyftBox data directory exists: {}",
            data_dir.display()
        )
    })?;
    let data_dir_string = data_dir.to_string_lossy().to_string();
    syftbox_config.data_dir = data_dir_string.clone();
    syftbox_config.email = email.clone();
    syftbox_config.server_url = server_url.clone();

    let client_url = resolve_client_url(client_url_override.as_deref(), &config, &syftbox_config);
    syftbox_config.client_url = Some(client_url.clone());
    syftbox_config.client_token = Some(tokens.access_token.clone());
    syftbox_config.refresh_token = Some(tokens.refresh_token.clone());

    syftbox_config.save(&syftbox_config_path)?;

    let mut creds = config.syftbox_credentials.unwrap_or_default();
    creds.email = Some(email.clone());
    creds.server_url = Some(server_url.clone());
    creds.client_url = Some(client_url.clone());
    creds.data_dir = Some(data_dir_string);
    creds.refresh_token = Some(tokens.refresh_token.clone());
    creds.access_token = Some(tokens.access_token.clone());
    config.syftbox_credentials = Some(creds);

    config.save(&config_path)?;

    println!(
        "Stored SyftBox credentials in {} and updated config {}",
        syftbox_config_path.display(),
        config_path.display()
    );

    Ok(())
}

fn join_url(base: &str, path: &str) -> String {
    let clean_base = base.trim_end_matches('/');
    let clean_path = path.trim_start_matches('/');
    format!("{}/{}", clean_base, clean_path)
}

fn resolve_email(email: Option<String>, config: Option<&Config>) -> Result<String> {
    if let Some(email) = email.filter(|s| !s.trim().is_empty()) {
        return Ok(email.trim().to_string());
    }

    if let Some(cfg) = config {
        if !cfg.email.trim().is_empty() {
            return Ok(cfg.email.trim().to_string());
        }
        if let Some(creds) = cfg.syftbox_credentials.as_ref() {
            if let Some(email) = creds.email.as_ref() {
                if !email.trim().is_empty() {
                    return Ok(email.trim().to_string());
                }
            }
        }
    }

    if let Ok(env_email) = env::var("SYFTBOX_EMAIL") {
        if !env_email.trim().is_empty() {
            return Ok(env_email.trim().to_string());
        }
    }

    Err(anyhow!(
        "Email address is required. Provide --email or run 'bv init <email>'."
    ))
}

fn resolve_server(server: Option<&str>, config: Option<&Config>) -> String {
    if let Some(server) = server {
        if !server.trim().is_empty() {
            return server.trim().trim_end_matches('/').to_string();
        }
    }

    if let Some(cfg) = config {
        if let Some(creds) = cfg.syftbox_credentials.as_ref() {
            if let Some(server_url) = creds.server_url.as_ref() {
                if !server_url.trim().is_empty() {
                    return server_url.trim().trim_end_matches('/').to_string();
                }
            }
        }
    }

    if let Ok(env_server) = env::var("SYFTBOX_SERVER_URL") {
        if !env_server.trim().is_empty() {
            return env_server.trim().trim_end_matches('/').to_string();
        }
    }

    DEFAULT_SERVER_URL.to_string()
}

fn resolve_client_url(
    client_url_override: Option<&str>,
    config: &Config,
    syftbox_config: &SyftboxConfigFile,
) -> String {
    if let Some(url) = client_url_override {
        if !url.trim().is_empty() {
            return url.trim().to_string();
        }
    }

    if let Some(url) = syftbox_config.client_url.as_ref() {
        if !url.trim().is_empty() {
            return url.trim().to_string();
        }
    }

    if let Some(creds) = config.syftbox_credentials.as_ref() {
        if let Some(url) = creds.client_url.as_ref() {
            if !url.trim().is_empty() {
                return url.trim().to_string();
            }
        }
    }

    if let Ok(env_url) = env::var("SYFTBOX_CLIENT_URL") {
        if !env_url.trim().is_empty() {
            return env_url.trim().to_string();
        }
    }

    DEFAULT_CLIENT_URL.to_string()
}

fn resolve_data_dir(
    override_dir: Option<&str>,
    config: &Config,
    syftbox_config: &SyftboxConfigFile,
) -> Result<PathBuf> {
    let default_data_dir = Config::default_syftbox_data_dir()?;
    let legacy_default = dirs::home_dir().map(|h| h.join("SyftBox"));

    if let Some(dir) = override_dir {
        if !dir.trim().is_empty() {
            return Ok(PathBuf::from(dir));
        }
    }

    if !syftbox_config.data_dir.trim().is_empty() {
        let candidate = PathBuf::from(syftbox_config.data_dir.trim());
        if let Some(legacy) = &legacy_default {
            if &candidate == legacy {
                return Ok(default_data_dir);
            }
        }
        return Ok(candidate);
    }

    if let Some(creds) = config.syftbox_credentials.as_ref() {
        if let Some(dir) = creds.data_dir.as_ref() {
            if !dir.trim().is_empty() {
                let candidate = PathBuf::from(dir.trim());
                if let Some(legacy) = &legacy_default {
                    if &candidate == legacy {
                        return Ok(default_data_dir);
                    }
                }
                return Ok(candidate);
            }
        }
    }

    if let Ok(env_dir) = env::var("SYFTBOX_DATA_DIR") {
        if !env_dir.trim().is_empty() {
            return Ok(PathBuf::from(env_dir.trim()));
        }
    }

    Ok(default_data_dir)
}

fn maybe_load_config() -> Result<Option<(Config, PathBuf)>> {
    let config_path = Config::get_config_path()?;
    match Config::load() {
        Ok(config) => Ok(Some((config, config_path))),
        Err(ConfigError::NotInitialized) => Ok(None),
        Err(e) => Err(e.into()),
    }
}

fn load_or_create_config(email: &str) -> Result<(Config, PathBuf)> {
    let config_path = Config::get_config_path()?;
    match Config::load() {
        Ok(config) => Ok((config, config_path)),
        Err(ConfigError::NotInitialized) => {
            let config = Config {
                email: email.to_string(),
                syftbox_config: None,
                version: None,
                binary_paths: None,
                syftbox_credentials: None,
            };
            config.save(&config_path)?;
            Ok((config, config_path))
        }
        Err(e) => Err(e.into()),
    }
}

async fn build_api_error(action: &str, response: reqwest::Response) -> anyhow::Error {
    let status = response.status();
    let body = response.text().await.unwrap_or_default();

    if let Ok(parsed) = serde_json::from_str::<ApiErrorMessage>(&body) {
        let message = parsed
            .message
            .or(parsed.error)
            .or(parsed.detail)
            .unwrap_or(body);
        anyhow!("{} failed ({}): {}", action, status, message)
    } else if body.is_empty() {
        anyhow!("{} failed ({})", action, status)
    } else {
        anyhow!("{} failed ({}): {}", action, status, body)
    }
}
