use biovault::cli::commands::key::{handle, KeyCommands};
use biovault::config::Config;
use biovault::syftbox::syc::{parse_public_bundle_file, provision_local_identity};
use std::sync::Mutex;
use syftbox_sdk::{identity_material_from_recovery_key, SyftRecoveryKey};
use tempfile::TempDir;

static ENV_MUTEX: Mutex<()> = Mutex::new(());

struct EnvGuard {
    key: &'static str,
    previous: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let previous = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, previous }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        if let Some(prev) = &self.previous {
            std::env::set_var(self.key, prev);
        } else {
            std::env::remove_var(self.key);
        }
    }
}

fn setup_config(email: &str) -> (TempDir, Config, EnvGuard, EnvGuard) {
    let temp = TempDir::new().unwrap();
    let data_dir = temp.path().join("data");
    std::fs::create_dir_all(&data_dir).unwrap();
    let syft_cfg = temp.path().join("syftbox-config.json");
    std::fs::write(
        &syft_cfg,
        format!(
            r#"{{"data_dir":"{}","email":"{}","server_url":""}}"#,
            data_dir.display(),
            email
        ),
    )
    .unwrap();
    let env_guard = EnvGuard::set("SYFTBOX_CONFIG_PATH", syft_cfg.to_str().unwrap());
    let vault_guard = EnvGuard::set(
        "SYC_VAULT",
        temp.path().join("vault").to_string_lossy().as_ref(),
    );

    let mut config = Config::new(email.to_string());
    config.syftbox_config = Some(syft_cfg.to_string_lossy().to_string());
    (temp, config, env_guard, vault_guard)
}

#[test]
fn generate_wipe_restore_round_trip() {
    let _lock = ENV_MUTEX.lock().unwrap();
    let (temp, config, _env, _vault) = setup_config("alice@example.com");
    let data_dir = temp.path().join("data");
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    // Generate via provisioning to capture mnemonic
    let outcome = provision_local_identity(&config.email, &data_dir, None).unwrap();
    let mnemonic = outcome
        .recovery_mnemonic
        .clone()
        .expect("mnemonic should be present");
    let encrypted_root = syftbox_sdk::syftbox::syc::resolve_encrypted_root(&data_dir);
    let export_path = encrypted_root
        .join(&config.email)
        .join("public")
        .join("crypto")
        .join("did.json");
    let initial = parse_public_bundle_file(&export_path).unwrap();

    // Wipe via CLI
    rt.block_on(async {
        handle(
            KeyCommands::Wipe {
                email: None,
                vault: Some(temp.path().join("vault")),
                data_dir: Some(data_dir.clone()),
                json: false,
            },
            &config,
        )
        .await
        .unwrap();
    });
    assert!(!export_path.exists(), "export should be removed after wipe");

    // Restore via CLI
    rt.block_on(async {
        handle(
            KeyCommands::Restore {
                email: config.email.clone(),
                mnemonic: mnemonic.clone(),
                vault: Some(temp.path().join("vault")),
                data_dir: Some(data_dir.clone()),
                json: false,
            },
            &config,
        )
        .await
        .unwrap();
    });

    let restored = parse_public_bundle_file(&export_path).unwrap();
    assert_eq!(initial.fingerprint, restored.fingerprint);
}

#[test]
fn tofu_conflict_requires_override() {
    let _lock = ENV_MUTEX.lock().unwrap();
    let (temp, config, _env, _vault) = setup_config("bob@example.com");
    let data_dir = temp.path().join("data");
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .unwrap();

    // Baseline identity
    let baseline = provision_local_identity(&config.email, &data_dir, None).unwrap();

    // Create conflicting bundle for same identity
    let other_material =
        identity_material_from_recovery_key(&config.email, &SyftRecoveryKey::generate()).unwrap();
    let conflict_bundle_path = temp.path().join("conflict.json");
    std::fs::write(
        &conflict_bundle_path,
        serde_json::to_vec_pretty(&other_material.public_bundle).unwrap(),
    )
    .unwrap();

    // Import without override should fail
    let err = rt
        .block_on(async {
            handle(
                KeyCommands::Import {
                    bundle: conflict_bundle_path.clone(),
                    email: Some(config.email.clone()),
                    ignore_tofu: false,
                    vault: Some(temp.path().join("vault")),
                    data_dir: Some(data_dir.clone()),
                    json: false,
                },
                &config,
            )
            .await
        })
        .unwrap_err();
    assert!(
        err.to_string().contains("TOFU"),
        "expected TOFU warning, got {err}"
    );

    // Import with --ignore-tofu should succeed and replace fingerprint
    rt.block_on(async {
        handle(
            KeyCommands::Import {
                bundle: conflict_bundle_path.clone(),
                email: Some(config.email.clone()),
                ignore_tofu: true,
                vault: Some(temp.path().join("vault")),
                data_dir: Some(data_dir.clone()),
                json: false,
            },
            &config,
        )
        .await
        .unwrap();
    });
    let replaced = parse_public_bundle_file(&baseline.public_bundle_path).unwrap();
    let new_fp = parse_public_bundle_file(&conflict_bundle_path)
        .unwrap()
        .fingerprint;
    assert_eq!(replaced.fingerprint, new_fp);
}
