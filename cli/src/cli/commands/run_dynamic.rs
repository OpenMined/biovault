use super::run::execute_with_logging;
use crate::data::module_yaml_hash;
use crate::error::Result;
use crate::module_spec::{ModuleSpec, ModuleStepSpec};
use anyhow::Context;
use chrono::Local;
use colored::Colorize;
use regex::Regex;
use serde_json::{json, Value as JsonValue};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::ffi::OsStr;
use std::fs::{self, OpenOptions};
use std::future::Future;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::net::{Ipv4Addr, SocketAddrV4, TcpListener};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

const SEQURE_COMMUNICATION_PORT_STRIDE: usize = 1000;
const SEQURE_DATA_SHARING_PORT_OFFSET: usize = 10_000;
const SEQURE_PORT_BASE_MIN: usize = 20_000;

#[derive(Debug, Clone, Default)]
pub struct DynamicExecutionContext {
    pub current_datasite: Option<String>,
    pub datasites_override: Option<Vec<String>>,
    pub syftbox_data_dir: Option<String>,
    pub run_id: Option<String>,
    pub flow_name: Option<String>,
    pub syqure_port_base: Option<usize>,
    pub tauri_context: bool,
}

tokio::task_local! {
    static EXECUTION_CONTEXT: DynamicExecutionContext;
}

pub async fn with_execution_context<F, T>(ctx: DynamicExecutionContext, fut: F) -> T
where
    F: Future<Output = T>,
{
    // Keep run-specific values task-local so concurrent parties do not race
    // through process-global env vars (a prior desktop-only failure mode).
    EXECUTION_CONTEXT.scope(ctx, fut).await
}

fn execution_context() -> Option<DynamicExecutionContext> {
    EXECUTION_CONTEXT.try_with(|ctx| ctx.clone()).ok()
}

pub fn mpc_comm_port_with_base(
    base: usize,
    local_pid: usize,
    remote_pid: usize,
    parties: usize,
) -> usize {
    let min_pid = std::cmp::min(local_pid, remote_pid);
    let max_pid = std::cmp::max(local_pid, remote_pid);
    let offset_major = min_pid * parties - min_pid * (min_pid + 1) / 2;
    let offset_minor = max_pid - min_pid;
    base + offset_major + offset_minor
}

/// Write a unified syqure diagnostics file to /tmp so both CLI and Tauri
/// produce identical output for comparison.
/// File: /tmp/biovault-syqure-diag-{run_id}-p{party_id}.log
fn write_syqure_diag_file(
    run_id: &str,
    party_id: usize,
    party_count: usize,
    email: &str,
    port_base: usize,
    tcp_proxy: bool,
    env_map: &BTreeMap<String, String>,
    module_path: &Path,
    binary: &str,
    entrypoint: &Path,
) {
    let safe_run_id: String = run_id
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    let path = PathBuf::from(format!(
        "/tmp/biovault-syqure-diag-{}-p{}.log",
        safe_run_id, party_id
    ));
    let caller = if execution_context()
        .map(|ctx| ctx.tauri_context)
        .unwrap_or(false)
    {
        "tauri"
    } else if env::var("TAURI_ENV").is_ok()
        || env::var("TAURI_PLATFORM").is_ok()
        || std::env::current_exe()
            .ok()
            .and_then(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().contains("biovault-desktop"))
            })
            .unwrap_or(false)
    {
        "tauri"
    } else {
        "cli"
    };
    let ts = chrono::Local::now().to_rfc3339();
    let mut lines: Vec<String> = Vec::new();

    lines.push(format!("timestamp={}", ts));
    lines.push(format!("caller={}", caller));
    lines.push(format!("run_id={}", run_id));
    lines.push(format!("party_id={}", party_id));
    lines.push(format!("party_count={}", party_count));
    lines.push(format!("email={}", email));
    lines.push(format!("port_base={}", port_base));
    lines.push(format!("tcp_proxy={}", tcp_proxy));
    lines.push(format!("binary={}", binary));
    lines.push(format!("entrypoint={}", entrypoint.display()));
    lines.push(format!("module_path={}", module_path.display()));
    lines.push(format!(
        "BV_SYFTBOX_BACKEND={}",
        env::var("BV_SYFTBOX_BACKEND").unwrap_or_else(|_| "unset".to_string())
    ));
    lines.push(format!(
        "SYFTBOX_HOTLINK_TCP_PROXY={}",
        env::var("SYFTBOX_HOTLINK_TCP_PROXY").unwrap_or_else(|_| "unset".to_string())
    ));

    let hint_path = syqure_port_base_hint_path(run_id, party_count);
    lines.push(format!(
        "port_hint_file={} exists={}",
        hint_path.display(),
        hint_path.exists()
    ));

    // Env vars (SEQURE_* and BV_*)
    for (k, v) in env_map
        .iter()
        .filter(|(k, _)| k.starts_with("SEQURE_") || k.starts_with("BV_"))
    {
        lines.push(format!("env.{}={}", k, v));
    }

    // TCP proxy port checks
    if tcp_proxy {
        for peer_id in 0..party_count {
            if peer_id == party_id {
                continue;
            }
            let comm_base = port_base + party_id * SEQURE_COMMUNICATION_PORT_STRIDE;
            let port = mpc_comm_port_with_base(comm_base, party_id, peer_id, party_count);
            let listening = TcpListener::bind(SocketAddrV4::new(Ipv4Addr::LOCALHOST, port as u16))
                .map(|_| false)
                .unwrap_or(true);
            lines.push(format!(
                "tcp_port_check: CP{}<->CP{} port={} status={}",
                party_id,
                peer_id,
                port,
                if listening {
                    "LISTENING"
                } else {
                    "NOT_LISTENING"
                }
            ));
        }
    }

    let content = lines.join("\n") + "\n";
    if let Err(e) = fs::write(&path, &content) {
        eprintln!(
            "warn: failed to write syqure diag file {}: {}",
            path.display(),
            e
        );
    } else {
        println!("  Diag file: {}", path.display());
    }
}

fn max_syqure_base_port(party_count: usize) -> usize {
    let parties = party_count.max(2);
    let max_party_base_delta = (parties - 1) * SEQURE_COMMUNICATION_PORT_STRIDE;
    let max_pair_offset = parties * (parties - 1) / 2;
    let reserve =
        max_party_base_delta + max_pair_offset + SEQURE_DATA_SHARING_PORT_OFFSET + parties;
    (u16::MAX as usize).saturating_sub(reserve)
}

fn required_syqure_ports_for_party(
    global_base: usize,
    party_id: usize,
    party_count: usize,
) -> Option<Vec<u16>> {
    let parties = party_count.max(2);
    if party_id >= parties {
        return None;
    }
    let party_base = global_base + party_id * SEQURE_COMMUNICATION_PORT_STRIDE;
    let mut ports: BTreeSet<u16> = BTreeSet::new();
    for remote_pid in 0..parties {
        if party_id == remote_pid {
            continue;
        }
        let port = mpc_comm_port_with_base(party_base, party_id, remote_pid, parties);
        ports.insert(u16::try_from(port).ok()?);
    }
    ports.insert(u16::try_from(party_base + SEQURE_DATA_SHARING_PORT_OFFSET).ok()?);
    Some(ports.into_iter().collect())
}

fn required_syqure_ports_for_all(global_base: usize, party_count: usize) -> Option<Vec<u16>> {
    let parties = party_count.max(2);
    let mut all_ports: BTreeSet<u16> = BTreeSet::new();
    for party_id in 0..parties {
        for port in required_syqure_ports_for_party(global_base, party_id, parties)? {
            all_ports.insert(port);
        }
    }
    Some(all_ports.into_iter().collect())
}

fn tcp_port_is_free(port: u16) -> bool {
    TcpListener::bind(SocketAddrV4::new(Ipv4Addr::LOCALHOST, port)).is_ok()
}

fn syqure_port_base_is_available(global_base: usize, party_count: usize) -> bool {
    let Some(ports) = required_syqure_ports_for_all(global_base, party_count) else {
        return false;
    };
    ports.into_iter().all(tcp_port_is_free)
}

fn syqure_port_base_hint_path(run_id: &str, party_count: usize) -> PathBuf {
    let safe_run_id: String = run_id
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    env::temp_dir().join(format!(
        "biovault-syqure-port-base-{}-p{}.txt",
        safe_run_id, party_count
    ))
}

fn read_syqure_port_base_hint(path: &Path) -> Option<usize> {
    let raw = fs::read_to_string(path).ok()?;
    raw.trim().parse::<usize>().ok()
}

fn validate_syqure_port_base(name: &str, base: usize, max_base: usize) -> Result<()> {
    if base < SEQURE_PORT_BASE_MIN || base > max_base {
        return Err(anyhow::anyhow!(
            "{} must be in range {}..={}, got {}",
            name,
            SEQURE_PORT_BASE_MIN,
            max_base,
            base
        )
        .into());
    }
    Ok(())
}

fn parse_port_base_env(name: &str) -> Result<Option<usize>> {
    match env::var(name) {
        Ok(value) => {
            let trimmed = value.trim();
            if trimmed.is_empty() {
                return Ok(None);
            }
            let parsed = trimmed.parse::<usize>().map_err(|_| {
                anyhow::anyhow!("{} must be a positive integer, got '{}'", name, trimmed)
            })?;
            if parsed == 0 || parsed > u16::MAX as usize {
                return Err(anyhow::anyhow!(
                    "{} must be in range 1..={}, got {}",
                    name,
                    u16::MAX,
                    parsed
                )
                .into());
            }
            Ok(Some(parsed))
        }
        Err(_) => Ok(None),
    }
}

pub fn prepare_syqure_port_base_for_run(
    run_id: &str,
    party_count: usize,
    requested_base_override: Option<usize>,
) -> Result<usize> {
    let max_base = max_syqure_base_port(party_count);
    if max_base <= SEQURE_PORT_BASE_MIN {
        return Err(anyhow::anyhow!(
            "Unable to allocate Syqure base port for {} parties (invalid port range)",
            party_count
        )
        .into());
    }

    if let Some(base) = requested_base_override {
        // In desktop multiparty, all parties in one session must share one base.
        // Do not derive a per-party base, or channels will silently diverge.
        validate_syqure_port_base("context syqure port base", base, max_base)?;
        env::set_var("BV_SYQURE_PORT_BASE", base.to_string());
        return Ok(base);
    }

    if let Some(base) = parse_port_base_env("BV_SYQURE_PORT_BASE")? {
        validate_syqure_port_base("BV_SYQURE_PORT_BASE", base, max_base)?;
        return Ok(base);
    }
    let requested_base = parse_port_base_env("SEQURE_COMMUNICATION_PORT")?;
    if let Some(base) = requested_base {
        validate_syqure_port_base("SEQURE_COMMUNICATION_PORT", base, max_base)?;
        env::set_var("BV_SYQURE_PORT_BASE", base.to_string());
        return Ok(base);
    }

    let hint_path = syqure_port_base_hint_path(run_id, party_count);
    if let Some(base) = read_syqure_port_base_hint(&hint_path) {
        validate_syqure_port_base("cached syqure port base", base, max_base)?;
        env::set_var("BV_SYQURE_PORT_BASE", base.to_string());
        return Ok(base);
    }

    let span = max_base - SEQURE_PORT_BASE_MIN + 1;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    run_id.hash(&mut hasher);
    party_count.hash(&mut hasher);
    let seed_candidate = SEQURE_PORT_BASE_MIN + (hasher.finish() as usize % span);
    let mut selected: Option<usize> = None;
    for offset in 0..span {
        let candidate =
            SEQURE_PORT_BASE_MIN + ((seed_candidate - SEQURE_PORT_BASE_MIN + offset) % span);
        if syqure_port_base_is_available(candidate, party_count) {
            selected = Some(candidate);
            break;
        }
    }
    let selected = selected.ok_or_else(|| {
        anyhow::anyhow!(
            "No free Syqure TCP proxy base port found for {} parties. \
             Tried {} candidate bases in [{}..={}].",
            party_count,
            span,
            SEQURE_PORT_BASE_MIN,
            max_base
        )
    })?;

    match OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&hint_path)
    {
        Ok(mut file) => {
            let _ = writeln!(file, "{}", selected);
        }
        Err(err) if err.kind() == io::ErrorKind::AlreadyExists => {
            if let Some(existing) = read_syqure_port_base_hint(&hint_path) {
                validate_syqure_port_base("cached syqure port base", existing, max_base)?;
                env::set_var("BV_SYQURE_PORT_BASE", existing.to_string());
                return Ok(existing);
            }
            let _ = fs::write(&hint_path, format!("{selected}\n"));
        }
        Err(_) => {}
    }
    env::set_var("BV_SYQURE_PORT_BASE", selected.to_string());
    Ok(selected)
}

fn is_truthy(value: &str) -> bool {
    matches!(
        value.trim().to_ascii_lowercase().as_str(),
        "1" | "true" | "yes" | "on"
    )
}

fn env_flag(names: &[&str]) -> Option<bool> {
    for name in names {
        if let Ok(value) = env::var(name) {
            return Some(is_truthy(&value));
        }
    }
    None
}

fn env_value(name: &str) -> String {
    env::var(name).unwrap_or_else(|_| "<unset>".to_string())
}

fn canonical_hotlink_transport(value: &str) -> String {
    let normalized = value.trim().to_ascii_lowercase();
    if normalized == "webrtc" {
        "hotlink".to_string()
    } else {
        normalized
    }
}

fn syqure_container_proxy_host(container_runtime: &str) -> String {
    if let Ok(value) = env::var("BV_SYQURE_CP_HOST") {
        let trimmed = value.trim();
        if !trimmed.is_empty() {
            return trimmed.to_string();
        }
    }
    if container_runtime.eq_ignore_ascii_case("podman") {
        "host.containers.internal".to_string()
    } else {
        // Prefer host alias (reaches the macOS/Windows host; gateway alias is only the VM router).
        "host.docker.internal".to_string()
    }
}

fn describe_syqure_transport_mode(env_map: &BTreeMap<String, String>) -> String {
    let sequre_transport = env_map
        .get("SEQURE_TRANSPORT")
        .cloned()
        .unwrap_or_else(|| "unknown".to_string())
        .trim()
        .to_ascii_lowercase();
    let bv_transport = canonical_hotlink_transport(&env::var("BV_SYQURE_TRANSPORT").unwrap_or_default());
    let hotlink_enabled = env_flag(&["BV_SYFTBOX_HOTLINK", "SYFTBOX_HOTLINK"]).unwrap_or(false);
    let p2p_enabled =
        env_flag(&["BV_SYFTBOX_HOTLINK_QUIC", "SYFTBOX_HOTLINK_QUIC"]).unwrap_or(true);
    let p2p_only =
        env_flag(&["BV_SYFTBOX_HOTLINK_QUIC_ONLY", "SYFTBOX_HOTLINK_QUIC_ONLY"]).unwrap_or(false);

    if sequre_transport == "tcp" && (bv_transport == "hotlink" || hotlink_enabled) {
        if p2p_only {
            return "hotlink over TCP proxy (p2p-only, websocket fallback disabled)".to_string();
        }
        if p2p_enabled {
            return "hotlink over TCP proxy (p2p preferred, websocket fallback enabled)"
                .to_string();
        }
        return "hotlink over TCP proxy (websocket only, p2p disabled)".to_string();
    }

    if sequre_transport == "file" || bv_transport == "file" {
        return "file transport".to_string();
    }

    format!("{} transport", sequre_transport)
}

fn is_hotlink_transport_mode(env_map: &BTreeMap<String, String>) -> bool {
    let sequre_transport = env_map
        .get("SEQURE_TRANSPORT")
        .map(|v| v.trim().to_ascii_lowercase())
        .unwrap_or_default();
    let bv_transport = canonical_hotlink_transport(&env::var("BV_SYQURE_TRANSPORT").unwrap_or_default());
    let hotlink_enabled = env_flag(&["BV_SYFTBOX_HOTLINK", "SYFTBOX_HOTLINK"]).unwrap_or(false);
    sequre_transport == "tcp" && (bv_transport == "hotlink" || hotlink_enabled)
}

fn count_mpc_files(dir: &Path) -> (usize, Option<u64>) {
    let mut count = 0usize;
    let mut latest_modified: Option<u64> = None;

    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                count += 1;
                if let Ok(metadata) = path.metadata() {
                    if let Ok(modified) = metadata.modified() {
                        let secs = modified
                            .duration_since(UNIX_EPOCH)
                            .unwrap_or_default()
                            .as_secs();
                        latest_modified = Some(latest_modified.map_or(secs, |lm| lm.max(secs)));
                    }
                }
            } else if path.is_dir() {
                let (sub_count, sub_latest) = count_mpc_files(&path);
                count += sub_count;
                if let Some(sl) = sub_latest {
                    latest_modified = Some(latest_modified.map_or(sl, |lm| lm.max(sl)));
                }
            }
        }
    }

    (count, latest_modified)
}

struct HotlinkTelemetrySummary {
    tx_packets: u64,
    tx_bytes: u64,
    tx_p2p_packets: u64,
    rx_packets: u64,
    rx_bytes: u64,
    webrtc_connected: u64,
}

fn read_hotlink_telemetry(path: &Path) -> Option<HotlinkTelemetrySummary> {
    let raw = fs::read_to_string(path).ok()?;
    let v: JsonValue = serde_json::from_str(&raw).ok()?;
    Some(HotlinkTelemetrySummary {
        tx_packets: v.get("tx_packets").and_then(|x| x.as_u64()).unwrap_or(0),
        tx_bytes: v.get("tx_bytes").and_then(|x| x.as_u64()).unwrap_or(0),
        tx_p2p_packets: v
            .get("tx_p2p_packets")
            .or_else(|| v.get("tx_quic_packets"))
            .and_then(|x| x.as_u64())
            .unwrap_or(0),
        rx_packets: v.get("rx_packets").and_then(|x| x.as_u64()).unwrap_or(0),
        rx_bytes: v.get("rx_bytes").and_then(|x| x.as_u64()).unwrap_or(0),
        webrtc_connected: v.get("webrtc_connected").and_then(|x| x.as_u64()).unwrap_or(0),
    })
}

fn fmt_hotlink_bytes(b: u64) -> String {
    if b >= 1_048_576 {
        format!("{:.1}MB", b as f64 / 1_048_576.0)
    } else if b >= 1024 {
        format!("{:.1}KB", b as f64 / 1024.0)
    } else {
        format!("{}B", b)
    }
}

fn multiparty_progress_dirs(
    datasites_root: &str,
    email: &str,
    flow_name: &str,
    run_id: &str,
) -> Vec<PathBuf> {
    let root = PathBuf::from(datasites_root);
    vec![
        root.join(email)
            .join("shared")
            .join("flows")
            .join(flow_name)
            .join(run_id)
            .join("_progress"),
        root.join(email)
            .join("datasites")
            .join(email)
            .join("shared")
            .join("flows")
            .join(flow_name)
            .join(run_id)
            .join("_progress"),
    ]
}

fn peer_failed_step_in_progress(progress_dir: &Path, step_id: &str) -> Option<String> {
    let entries = fs::read_dir(progress_dir).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if !name.ends_with(".json") || name == "state.json" || name == "progress.json" {
            continue;
        }
        let raw = match fs::read_to_string(&path) {
            Ok(raw) => raw,
            Err(_) => continue,
        };
        let parsed: JsonValue = match serde_json::from_str(&raw) {
            Ok(parsed) => parsed,
            Err(_) => continue,
        };
        let parsed_step = parsed
            .get("step_id")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        if parsed_step != step_id {
            continue;
        }
        let status = parsed
            .get("status")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        if !status.eq_ignore_ascii_case("Failed") {
            continue;
        }
        let role = parsed.get("role").and_then(|v| v.as_str()).unwrap_or("");
        if role.is_empty() {
            return Some("unknown".to_string());
        }
        return Some(role.to_string());
    }
    None
}

fn detect_peer_step_failure(
    datasites_root: &str,
    flow_name: &str,
    run_id: &str,
    step_id: &str,
    local_email: &str,
    party_emails: &[String],
) -> Option<String> {
    for peer_email in party_emails {
        if peer_email.eq_ignore_ascii_case(local_email) {
            continue;
        }
        for progress_dir in multiparty_progress_dirs(datasites_root, peer_email, flow_name, run_id)
        {
            if !progress_dir.exists() {
                continue;
            }
            if let Some(role) = peer_failed_step_in_progress(&progress_dir, step_id) {
                return Some(format!("{} ({})", peer_email, role));
            }
        }
    }
    None
}

fn append_desktop_log(message: &str) {
    if let Ok(path) = std::env::var("BIOVAULT_DESKTOP_LOG_FILE") {
        if path.is_empty() {
            return;
        }
        let path = PathBuf::from(path);
        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(&path) {
            let timestamp = Local::now().format("%Y-%m-%dT%H:%M:%S%:z");
            let _ = writeln!(file, "[{}][INFO] {}", timestamp, message);
        }
    }
}

fn shell_quote(value: &OsStr) -> String {
    let s = value.to_string_lossy();
    if s.is_empty() {
        return "''".to_string();
    }
    if s.chars()
        .all(|c| c.is_ascii_alphanumeric() || "-_./:@".contains(c))
    {
        return s.into_owned();
    }
    let escaped = s.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

fn format_command(cmd: &Command) -> String {
    let mut parts = Vec::new();
    parts.push(shell_quote(cmd.get_program()));
    for arg in cmd.get_args() {
        parts.push(shell_quote(arg));
    }
    parts.join(" ")
}

fn spawn_child_stream_capture<R: std::io::Read + Send + 'static>(
    reader: R,
    stream_name: &'static str,
    log_path: Option<PathBuf>,
    echo_stderr: bool,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let mut writer = log_path
            .as_ref()
            .and_then(|path| OpenOptions::new().create(true).append(true).open(path).ok());
        let buffered = BufReader::new(reader);
        for line in buffered.lines().map_while(|res| res.ok()) {
            if echo_stderr && stream_name == "stderr" {
                eprintln!("  [syqure:{}] {}", stream_name, line);
            }
            if let Some(file) = writer.as_mut() {
                let _ = writeln!(file, "[{}] {}", stream_name, line);
            }
        }
    })
}

fn bundled_env_var(name: &str) -> Option<&'static str> {
    match name {
        "java" => Some("BIOVAULT_BUNDLED_JAVA"),
        "nextflow" => Some("BIOVAULT_BUNDLED_NEXTFLOW"),
        "uv" => Some("BIOVAULT_BUNDLED_UV"),
        "syftbox" => Some("SYFTBOX_BINARY"),
        _ => None,
    }
}

fn resolve_binary_path(cfg: Option<&crate::config::Config>, name: &str) -> Option<String> {
    if let Some(cfg) = cfg {
        if let Some(path) = cfg.get_binary_path(name) {
            if !path.is_empty() {
                return Some(path);
            }
        }
    }

    if let Some(env_key) = bundled_env_var(name) {
        if let Ok(env_path) = std::env::var(env_key) {
            let trimmed = env_path.trim();
            if !trimmed.is_empty() {
                return Some(trimmed.to_string());
            }
        }
    }

    None
}

fn build_augmented_path(cfg: Option<&crate::config::Config>) -> Option<String> {
    let mut entries = BTreeSet::new();
    for key in ["nextflow", "java", "docker"] {
        if let Some(bin_path) = resolve_binary_path(cfg, key) {
            if bin_path.is_empty() {
                continue;
            }
            if let Some(parent) = Path::new(&bin_path).parent() {
                entries.insert(parent.to_path_buf());
            }
        }
    }

    if entries.is_empty() {
        return None;
    }

    let mut paths: Vec<PathBuf> = entries.into_iter().collect();
    if let Some(existing) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&existing));
    }

    std::env::join_paths(paths)
        .ok()
        .and_then(|joined| joined.into_string().ok())
}

#[derive(Debug, Clone, Copy)]
pub struct RunSettings {
    /// Whether Docker must be running before launching Nextflow.
    pub require_docker: bool,
}

impl Default for RunSettings {
    fn default() -> Self {
        Self {
            // Current templates assume Docker; keep default strict.
            require_docker: true,
        }
    }
}

/// Check if we should use Docker to run Nextflow
/// - Windows: Always (no native Nextflow)
/// - Linux/macOS ARM64: When forced via env var (native Nextflow works)
/// - Linux/macOS x86: Never by default (native Nextflow works)
fn should_use_docker_for_nextflow() -> bool {
    // Allow forcing container mode for testing on non-Windows platforms
    if std::env::var("BIOVAULT_FORCE_CONTAINER_NEXTFLOW")
        .map(|v| !v.is_empty() && v != "0" && v.to_lowercase() != "false")
        .unwrap_or(false)
    {
        return true;
    }
    cfg!(target_os = "windows")
}

/// Detect if running on ARM64 architecture
fn is_arm64() -> bool {
    std::env::consts::ARCH == "aarch64"
}

/// Check if x86 container emulation should be forced for all containers
fn should_force_x86_containers() -> bool {
    std::env::var("BIOVAULT_FORCE_X86_CONTAINERS")
        .map(|v| !v.is_empty() && v != "0" && v.to_lowercase() != "false")
        .unwrap_or(false)
}

/// Get the appropriate platform flag for Docker containers
/// Returns None if Docker should auto-select (preferred for multi-arch images)
/// Returns Some("linux/amd64") if x86 emulation is needed
/// Returns Some("linux/arm64") only when explicitly required
#[allow(dead_code)]
fn get_container_platform(force_x86: bool) -> Option<&'static str> {
    if force_x86 || should_force_x86_containers() {
        return Some("linux/amd64");
    }
    // Let Docker auto-select - it will use native if available, emulate if not
    None
}

/// Extract container image names from a workflow.nf file
fn extract_containers_from_workflow(workflow_path: &Path) -> Vec<String> {
    let mut containers = Vec::new();

    let content = match fs::read_to_string(workflow_path) {
        Ok(c) => c,
        Err(_) => return containers,
    };

    // Match patterns like: container 'image:tag' or container "image:tag"
    let re = Regex::new(r#"container\s+['"]([^'"]+)['"]"#).ok();
    if let Some(regex) = re {
        for cap in regex.captures_iter(&content) {
            if let Some(image) = cap.get(1) {
                containers.push(image.as_str().to_string());
            }
        }
    }

    // Also match container params.xxx_container patterns
    let param_re = Regex::new(r#"container\s+params\.(\w+)"#).ok();
    if let Some(regex) = param_re {
        for cap in regex.captures_iter(&content) {
            if let Some(_param_name) = cap.get(1) {
                // These are dynamic, we can't check them statically
                // They'll be handled at runtime
            }
        }
    }

    containers
}

/// Check if a container image has arm64 support
/// First checks locally pulled images, then falls back to registry manifest
fn check_container_has_arm64(image: &str) -> bool {
    // First, try to check the local image (works offline)
    let local_output = Command::new("docker")
        .args(["image", "inspect", "--format", "{{.Architecture}}", image])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output();

    if let Ok(out) = local_output {
        if out.status.success() {
            let arch = String::from_utf8_lossy(&out.stdout).trim().to_lowercase();
            if arch == "arm64" || arch == "aarch64" {
                append_desktop_log(&format!(
                    "[Platform] Local image {} is arm64 (native)",
                    image
                ));
                return true;
            } else if arch == "amd64" || arch == "x86_64" {
                append_desktop_log(&format!(
                    "[Platform] Local image {} is {} (will emulate)",
                    image, arch
                ));
                return false;
            }
            // Unknown arch, continue to manifest check
        }
    }

    // Image not pulled locally, try remote manifest inspect (requires network)
    let manifest_output = Command::new("docker")
        .args(["manifest", "inspect", image])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output();

    match manifest_output {
        Ok(out) => {
            if !out.status.success() {
                // If manifest inspect fails (offline or image doesn't exist),
                // assume native support to avoid unnecessary emulation
                // The actual docker pull will fail if image doesn't exist
                append_desktop_log(&format!(
                    "[Platform] Cannot inspect manifest for {} (offline?), assuming native support",
                    image
                ));
                return true;
            }
            let stdout = String::from_utf8_lossy(&out.stdout);
            // Check for arm64/aarch64 architecture in manifest
            let has_arm64 = stdout.contains("\"architecture\": \"arm64\"")
                || stdout.contains("\"architecture\":\"arm64\"")
                || stdout.contains("\"aarch64\"")
                || stdout.contains("/arm64");

            if has_arm64 {
                append_desktop_log(&format!(
                    "[Platform] Registry manifest for {} includes arm64",
                    image
                ));
            } else {
                append_desktop_log(&format!(
                    "[Platform] Registry manifest for {} is x86-only",
                    image
                ));
            }
            has_arm64
        }
        Err(_) => {
            // docker not available, assume native support
            append_desktop_log("[Platform] Docker not available for manifest check");
            true
        }
    }
}

/// Check all containers in a workflow and return true if any lack arm64 support
fn any_container_lacks_arm64(workflow_path: &Path) -> bool {
    if !is_arm64() {
        // Not on arm64, no need to check
        return false;
    }

    let containers = extract_containers_from_workflow(workflow_path);
    if containers.is_empty() {
        return false;
    }

    append_desktop_log(&format!(
        "[Platform] Checking {} container(s) for arm64 support on arm64 host",
        containers.len()
    ));

    for image in &containers {
        if !check_container_has_arm64(image) {
            append_desktop_log(&format!(
                "[Platform] Container {} lacks arm64 support, will use x86 emulation",
                image
            ));
            return true;
        } else {
            append_desktop_log(&format!("[Platform] Container {} has arm64 support", image));
        }
    }

    false
}

/// Find Git Bash on Windows, falling back to "bash" if not found.
#[cfg(target_os = "windows")]
fn resolve_git_bash() -> String {
    let candidates = [
        r"C:\Program Files\Git\bin\bash.exe",
        r"C:\Program Files\Git\usr\bin\bash.exe",
        r"C:\Program Files (x86)\Git\bin\bash.exe",
        r"C:\Program Files (x86)\Git\usr\bin\bash.exe",
    ];
    for candidate in candidates {
        if Path::new(candidate).exists() {
            return candidate.to_string();
        }
    }
    // Fallback to PATH lookup
    "bash".to_string()
}

#[cfg(not(target_os = "windows"))]
fn resolve_git_bash() -> String {
    "bash".to_string()
}

/// Strip the Windows extended path prefix (\\?\) for use with Docker volume mounts.
/// Docker Desktop on Windows expects standard Windows paths like C:\Users\foo,
/// not extended paths like \\?\C:\Users\foo.
#[cfg(target_os = "windows")]
fn strip_extended_path_prefix(path: &str) -> String {
    let stripped = path.trim_start_matches(r"\\?\").trim_start_matches(r"//?/");
    stripped.to_string()
}

#[cfg(not(target_os = "windows"))]
fn strip_extended_path_prefix(path: &str) -> String {
    path.to_string()
}

/// Detect when "docker" is actually a Podman shim (common on Windows).
#[cfg(target_os = "windows")]
fn is_podman_shim(binary: &str) -> bool {
    let mut cmd = Command::new(binary);
    super::configure_child_process(&mut cmd);
    let output = cmd.arg("--version").output();
    if let Ok(output) = output {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let combined = format!("{}{}", stdout, stderr).to_lowercase();
        return combined.contains("podman");
    }
    false
}

/// Convert a Windows path to a container-compatible path for use *inside* the container.
/// - Docker Desktop: C:\Users\foo -> /c/Users/foo
/// - Podman on WSL: C:\Users\foo -> /mnt/c/Users/foo
#[cfg(target_os = "windows")]
fn windows_path_to_container(path: &Path, use_podman_format: bool) -> String {
    let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    let path_str = canonical.to_string_lossy();
    // Convert backslashes to forward slashes
    let unix_path = path_str.replace('\\', "/");
    // Prefix for drive letter mount
    let prefix = if use_podman_format { "/mnt" } else { "" };
    // Convert drive letter: C:/... -> /c/... or /mnt/c/...
    if unix_path.len() >= 2 && unix_path.chars().nth(1) == Some(':') {
        let drive = unix_path
            .chars()
            .next()
            .unwrap()
            .to_lowercase()
            .next()
            .unwrap();
        format!("{}/{}{}", prefix, drive, &unix_path[2..])
    } else if unix_path.starts_with("\\\\?\\") || unix_path.starts_with("//?/") {
        // Handle extended path prefix \\?\C:\... or //?/C:/...
        let stripped = unix_path
            .trim_start_matches("\\\\?\\")
            .trim_start_matches("//?/");
        let stripped = stripped.replace('\\', "/");
        if stripped.len() >= 2 && stripped.chars().nth(1) == Some(':') {
            let drive = stripped
                .chars()
                .next()
                .unwrap()
                .to_lowercase()
                .next()
                .unwrap();
            format!("{}/{}{}", prefix, drive, &stripped[2..])
        } else {
            format!("{}/{}", prefix, stripped)
        }
    } else {
        unix_path
    }
}

/// Convert a Windows path to a host-side mount source.
/// For Windows runtimes, the host side should be a Windows path.
#[cfg(target_os = "windows")]
fn windows_path_to_mount_source(path: &Path, _use_podman_format: bool) -> String {
    let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    strip_extended_path_prefix(canonical.to_string_lossy().as_ref())
}

#[cfg(not(target_os = "windows"))]
fn windows_path_to_container(_path: &Path, _use_podman_format: bool) -> String {
    _path.to_string_lossy().to_string()
}

#[cfg(not(target_os = "windows"))]
fn windows_path_to_mount_source(path: &Path, _use_podman_format: bool) -> String {
    path.to_string_lossy().to_string()
}

/// Recursively remap Windows paths in JSON values to container-compatible paths
fn remap_json_paths_for_container(value: &JsonValue, use_podman_format: bool) -> JsonValue {
    match value {
        JsonValue::String(s) => {
            // Check if it looks like a Windows path (contains backslash or drive letter)
            if s.contains('\\') || (s.len() >= 2 && s.chars().nth(1) == Some(':')) {
                JsonValue::String(windows_path_to_container(Path::new(s), use_podman_format))
            } else {
                value.clone()
            }
        }
        JsonValue::Object(map) => {
            let mut new_map = serde_json::Map::new();
            for (k, v) in map {
                new_map.insert(
                    k.clone(),
                    remap_json_paths_for_container(v, use_podman_format),
                );
            }
            JsonValue::Object(new_map)
        }
        JsonValue::Array(arr) => JsonValue::Array(
            arr.iter()
                .map(|v| remap_json_paths_for_container(v, use_podman_format))
                .collect(),
        ),
        _ => value.clone(),
    }
}

/// Remap paths in JSON to VM-local paths (for Hyper-V mode)
/// Windows paths are converted to point to the vm_data_dir by filename
#[cfg(target_os = "windows")]
fn remap_json_paths_for_vm(value: &JsonValue, vm_data_dir: &str) -> JsonValue {
    match value {
        JsonValue::String(s) => {
            // Check if it looks like a Windows path
            if s.contains('\\') || (s.len() >= 2 && s.chars().nth(1) == Some(':')) {
                let path = Path::new(s);
                let file_name = path.file_name().unwrap_or_default().to_string_lossy();
                JsonValue::String(format!("{}/{}", vm_data_dir, file_name))
            } else {
                value.clone()
            }
        }
        JsonValue::Object(map) => {
            let mut new_map = serde_json::Map::new();
            for (k, v) in map {
                new_map.insert(k.clone(), remap_json_paths_for_vm(v, vm_data_dir));
            }
            JsonValue::Object(new_map)
        }
        JsonValue::Array(arr) => JsonValue::Array(
            arr.iter()
                .map(|v| remap_json_paths_for_vm(v, vm_data_dir))
                .collect(),
        ),
        _ => value.clone(),
    }
}

#[cfg(target_os = "windows")]
/// Remap paths in JSON to VM-local paths using an explicit staging map.
fn remap_json_paths_for_vm_with_map(
    value: &JsonValue,
    vm_data_dir: &str,
    path_map: &BTreeMap<String, String>,
) -> JsonValue {
    match value {
        JsonValue::String(s) => {
            if looks_like_windows_absolute_path(s) {
                let normalized = normalize_windows_path_str(s).to_ascii_lowercase();
                if let Some(mapped) = path_map.get(&normalized) {
                    JsonValue::String(mapped.replace('\\', "/"))
                } else {
                    let path = Path::new(s);
                    let file_name = path.file_name().unwrap_or_default().to_string_lossy();
                    JsonValue::String(format!("{}/{}", vm_data_dir, file_name))
                }
            } else {
                value.clone()
            }
        }
        JsonValue::Object(map) => {
            let mut new_map = serde_json::Map::new();
            for (k, v) in map {
                new_map.insert(
                    k.clone(),
                    remap_json_paths_for_vm_with_map(v, vm_data_dir, path_map),
                );
            }
            JsonValue::Object(new_map)
        }
        JsonValue::Array(arr) => JsonValue::Array(
            arr.iter()
                .map(|v| remap_json_paths_for_vm_with_map(v, vm_data_dir, path_map))
                .collect(),
        ),
        _ => value.clone(),
    }
}

/// Normalize a Windows path string (strip extended prefix, convert to backslashes)
#[cfg(target_os = "windows")]
fn normalize_windows_path_str(s: &str) -> String {
    // Strip extended-length path prefix if present
    let stripped = s
        .strip_prefix("\\\\?\\")
        .or_else(|| s.strip_prefix("//?/"))
        .unwrap_or(s);
    // Handle MSYS/Git Bash style paths like /c/Users/... or /mnt/c/Users/...
    if let Some(rest) = stripped.strip_prefix("/mnt/") {
        if rest.len() >= 3 && rest.as_bytes()[1] == b'/' {
            let drive = rest.chars().next().unwrap_or('c').to_ascii_uppercase();
            let remainder = &rest[2..];
            return format!(
                "{}:\\{}",
                drive,
                remainder.replace('/', "\\").trim_start_matches("\\")
            );
        }
    }
    if stripped.len() >= 3 && stripped.starts_with('/') && stripped.as_bytes()[2] == b'/' {
        let drive = stripped.chars().nth(1).unwrap_or('c').to_ascii_uppercase();
        let remainder = &stripped[3..];
        return format!(
            "{}:\\{}",
            drive,
            remainder.replace('/', "\\").trim_start_matches("\\")
        );
    }
    // Convert forward slashes to backslashes
    stripped.replace('/', "\\")
}

#[cfg(target_os = "windows")]
/// Remap paths in JSON to a flat staging dir using an explicit staging map.
fn remap_json_paths_for_flat_dir_with_map(
    value: &JsonValue,
    flat_data_dir: &str,
    path_map: &BTreeMap<String, String>,
) -> JsonValue {
    match value {
        JsonValue::String(s) => {
            if looks_like_windows_absolute_path(s) {
                let normalized = normalize_windows_path_str(s).to_ascii_lowercase();
                if let Some(mapped) = path_map.get(&normalized) {
                    JsonValue::String(mapped.replace('\\', "/"))
                } else {
                    let path = Path::new(s);
                    let file_name = path.file_name().unwrap_or_default().to_string_lossy();
                    JsonValue::String(format!("{}/{}", flat_data_dir, file_name))
                }
            } else {
                value.clone()
            }
        }
        JsonValue::Object(map) => {
            let mut new_map = serde_json::Map::new();
            for (k, v) in map {
                new_map.insert(
                    k.clone(),
                    remap_json_paths_for_flat_dir_with_map(v, flat_data_dir, path_map),
                );
            }
            JsonValue::Object(new_map)
        }
        JsonValue::Array(arr) => JsonValue::Array(
            arr.iter()
                .map(|v| remap_json_paths_for_flat_dir_with_map(v, flat_data_dir, path_map))
                .collect(),
        ),
        _ => value.clone(),
    }
}

#[cfg(target_os = "windows")]
fn normalize_windows_path_key(path: &Path) -> String {
    normalize_windows_path_str(&path.to_string_lossy()).to_ascii_lowercase()
}

#[cfg(target_os = "windows")]
fn sanitize_filename(name: &str) -> String {
    let mut out = String::new();
    for c in name.chars() {
        if c.is_ascii_alphanumeric() || c == '.' || c == '-' || c == '_' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "file".to_string()
    } else {
        out
    }
}

#[cfg(target_os = "windows")]
fn stage_name_for_path(path: &Path) -> String {
    let normalized = normalize_windows_path_str(&path.to_string_lossy());
    let hash = blake3::hash(normalized.as_bytes()).to_hex().to_string();
    let base = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let base = sanitize_filename(&base);
    format!("{}_{}", &hash[..10], base)
}

#[cfg(target_os = "windows")]
fn env_var_truthy(name: &str) -> bool {
    match std::env::var(name) {
        Ok(value) => {
            let normalized = value.trim().to_ascii_lowercase();
            !normalized.is_empty()
                && normalized != "0"
                && normalized != "false"
                && normalized != "no"
        }
        Err(_) => false,
    }
}

#[cfg(not(target_os = "windows"))]
fn env_var_truthy(name: &str) -> bool {
    match std::env::var(name) {
        Ok(value) => {
            let normalized = value.trim().to_ascii_lowercase();
            !normalized.is_empty()
                && normalized != "0"
                && normalized != "false"
                && normalized != "no"
        }
        Err(_) => false,
    }
}

/// Normalize and lowercase a Windows path for case-insensitive comparisons.
#[cfg(target_os = "windows")]
fn normalize_windows_path_for_compare(path: &Path) -> PathBuf {
    let normalized = normalize_windows_path_str(&path.to_string_lossy());
    PathBuf::from(normalized.to_ascii_lowercase())
}

/// Return true if `target` is within `base` after normalizing Windows paths.
#[cfg(target_os = "windows")]
fn is_path_within(base: &Path, target: &Path) -> bool {
    let base_norm = normalize_windows_path_for_compare(base);
    let target_norm = normalize_windows_path_for_compare(target);
    target_norm.starts_with(&base_norm)
}

#[cfg(not(target_os = "windows"))]
fn is_path_within(base: &Path, target: &Path) -> bool {
    target.starts_with(base)
}

/// Compute a relative path from base to target after normalizing Windows paths.
#[cfg(target_os = "windows")]
fn relative_to_base(base: &Path, target: &Path) -> PathBuf {
    let base_norm = normalize_windows_path_for_compare(base);
    let target_norm = normalize_windows_path_for_compare(target);
    target_norm
        .strip_prefix(&base_norm)
        .unwrap_or(&target_norm)
        .to_path_buf()
}

#[cfg(not(target_os = "windows"))]
fn relative_to_base(base: &Path, target: &Path) -> PathBuf {
    target.strip_prefix(base).unwrap_or(target).to_path_buf()
}

/// Extract all file/directory paths from a JSON value (for Docker volume mounting)
#[cfg(target_os = "windows")]
fn extract_paths_from_json(value: &JsonValue, paths: &mut Vec<PathBuf>) {
    match value {
        JsonValue::String(s) => {
            // Check if it looks like a Windows path (with drive letter)
            if looks_like_windows_absolute_path(s) {
                // Normalize: strip extended prefix and convert forward slashes to backslashes
                let normalized = normalize_windows_path_str(s);
                let path = Path::new(&normalized);
                append_desktop_log(&format!(
                    "[JSON Extract] Found path: {} (normalized: {}, exists: {})",
                    s,
                    normalized,
                    path.exists()
                ));
                if path.exists() {
                    // Get the parent directory for files, or the path itself for directories
                    if path.is_file() {
                        if let Some(parent) = path.parent() {
                            append_desktop_log(&format!(
                                "[JSON Extract] Adding parent dir: {}",
                                parent.display()
                            ));
                            paths.push(parent.to_path_buf());
                        }
                        // If it's a CSV file, also extract paths from inside it
                        if s.to_lowercase().ends_with(".csv") {
                            append_desktop_log(&format!(
                                "[JSON Extract] Reading CSV for embedded paths: {}",
                                s
                            ));
                            if let Ok(content) = fs::read_to_string(path) {
                                extract_paths_from_csv(&content, paths);
                            } else {
                                append_desktop_log(&format!(
                                    "[JSON Extract] ERROR: Failed to read CSV: {}",
                                    s
                                ));
                            }
                        }
                    } else {
                        append_desktop_log(&format!(
                            "[JSON Extract] Adding directory: {}",
                            path.display()
                        ));
                        paths.push(path.to_path_buf());
                    }
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                extract_paths_from_json(v, paths);
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                extract_paths_from_json(v, paths);
            }
        }
        _ => {}
    }
}

/// Check if a string looks like a Windows absolute path
#[cfg(target_os = "windows")]
fn looks_like_windows_absolute_path(s: &str) -> bool {
    // Handle extended-length path prefix: \\?\C:\... or //?/C:/...
    let stripped = s
        .strip_prefix("\\\\?\\")
        .or_else(|| s.strip_prefix("//?/"))
        .unwrap_or(s);

    // Accept MSYS/Git Bash style: /c/... or /mnt/c/...
    if let Some(rest) = stripped.strip_prefix("/mnt/") {
        if rest.len() >= 3 && rest.as_bytes()[1] == b'/' {
            let drive = rest.chars().next().unwrap_or('c');
            return drive.is_ascii_alphabetic();
        }
    }
    if stripped.len() >= 3 && stripped.starts_with('/') && stripped.as_bytes()[2] == b'/' {
        let drive = stripped.chars().nth(1).unwrap_or('c');
        return drive.is_ascii_alphabetic();
    }

    if stripped.len() < 3 {
        return false;
    }

    let chars: Vec<char> = stripped.chars().take(3).collect();
    // Check for drive letter pattern: X:\ or X:/
    if chars.len() >= 3 && chars[1] == ':' {
        let first = chars[0];
        let third = chars[2];
        return first.is_ascii_alphabetic() && (third == '\\' || third == '/');
    }
    false
}

/// Extract Windows paths from CSV content
#[cfg(target_os = "windows")]
fn extract_paths_from_csv(content: &str, paths: &mut Vec<PathBuf>) {
    append_desktop_log(&format!(
        "[CSV Extract] Processing CSV content ({} bytes, {} lines)",
        content.len(),
        content.lines().count()
    ));

    let mut found_count = 0;
    let mut added_count = 0;

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_reader(content.as_bytes());
    for result in reader.records() {
        if let Ok(record) = result {
            for field in record.iter() {
                let field = field.trim();
                // Accept both backslash (C:\) and forward slash (C:/) Windows paths
                if looks_like_windows_absolute_path(field) {
                    found_count += 1;
                    // Normalize: strip extended prefix and convert forward slashes to backslashes
                    let normalized = normalize_windows_path_str(field);
                    let path = Path::new(&normalized);
                    append_desktop_log(&format!(
                        "[CSV Extract] Found path: {} (normalized: {}, exists: {})",
                        field,
                        normalized,
                        path.exists()
                    ));
                    if path.exists() {
                        if path.is_file() {
                            if let Some(parent) = path.parent() {
                                append_desktop_log(&format!(
                                    "[CSV Extract] Adding parent dir: {}",
                                    parent.display()
                                ));
                                paths.push(parent.to_path_buf());
                                added_count += 1;
                            }
                        } else {
                            append_desktop_log(&format!(
                                "[CSV Extract] Adding directory: {}",
                                path.display()
                            ));
                            paths.push(path.to_path_buf());
                            added_count += 1;
                        }
                    }
                }
            }
        }
    }

    append_desktop_log(&format!(
        "[CSV Extract] Found {} paths, added {} mount candidates",
        found_count, added_count
    ));
}

/// Rewrite a CSV file converting Windows paths to container-compatible paths
#[cfg(target_os = "windows")]
fn rewrite_csv_with_container_paths(csv_path: &Path, use_podman_format: bool) -> Result<()> {
    append_desktop_log(&format!("[CSV Rewrite] Processing: {}", csv_path.display()));
    let mut converted_count = 0;

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(csv_path)
        .with_context(|| format!("Failed to open CSV: {}", csv_path.display()))?;
    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .from_writer(vec![]);

    for result in reader.records() {
        let record = result
            .with_context(|| format!("Failed to parse CSV record: {}", csv_path.display()))?;
        let mut out: Vec<String> = Vec::with_capacity(record.len());
        for field in record.iter() {
            if looks_like_windows_absolute_path(field) {
                converted_count += 1;
                out.push(windows_path_to_container(
                    Path::new(field),
                    use_podman_format,
                ));
            } else {
                out.push(field.to_string());
            }
        }
        writer
            .write_record(out)
            .with_context(|| format!("Failed to write CSV record: {}", csv_path.display()))?;
    }

    append_desktop_log(&format!(
        "[CSV Rewrite] Converted {} paths in {}",
        converted_count,
        csv_path.display()
    ));

    let new_content = String::from_utf8(
        writer
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to finalize CSV writer: {}", e))?,
    )
    .context("Failed to encode rewritten CSV")?;
    fs::write(csv_path, new_content)
        .with_context(|| format!("Failed to write converted CSV: {}", csv_path.display()))?;

    Ok(())
}

/// Check if a string looks like a Windows path
#[cfg(target_os = "windows")]
fn is_windows_path(s: &str) -> bool {
    // Regular path: C:\...
    if s.len() >= 2 && s.chars().nth(1) == Some(':') {
        return true;
    }
    // Extended path: \\?\C:\...
    if s.starts_with("\\\\?\\") || s.starts_with("//?/") {
        return true;
    }
    false
}

/// Find and rewrite all CSV files referenced in inputs JSON
#[cfg(target_os = "windows")]
fn rewrite_input_csvs_for_container(value: &JsonValue, use_podman_format: bool) -> Result<()> {
    match value {
        JsonValue::String(s) => {
            if is_windows_path(s) && s.to_lowercase().ends_with(".csv") {
                let path = Path::new(s);
                if path.exists() && path.is_file() {
                    append_desktop_log(&format!("[Pipeline] Rewriting CSV: {}", s));
                    rewrite_csv_with_container_paths(path, use_podman_format)?;
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                rewrite_input_csvs_for_container(v, use_podman_format)?;
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                rewrite_input_csvs_for_container(v, use_podman_format)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Find and rewrite all CSV files referenced in inputs JSON to use VM-local paths.
/// This is used for Hyper-V mode where Windows paths can't be mounted directly.
#[cfg(target_os = "windows")]
fn rewrite_input_csvs_for_vm(
    value: &JsonValue,
    vm_data_dir: &str,
    path_map: &BTreeMap<String, String>,
) -> Result<()> {
    match value {
        JsonValue::String(s) => {
            if is_windows_path(s) && s.to_lowercase().ends_with(".csv") {
                let path = Path::new(s);
                if path.exists() && path.is_file() {
                    append_desktop_log(&format!("[Pipeline] Rewriting CSV for VM: {}", s));
                    rewrite_csv_for_vm(path, vm_data_dir, path_map)?;
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                rewrite_input_csvs_for_vm(v, vm_data_dir, path_map)?;
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                rewrite_input_csvs_for_vm(v, vm_data_dir, path_map)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Rewrite a CSV file to use VM-local paths instead of Windows paths.
/// Windows paths like C:/Users/foo/data/file.txt become /tmp/biovault-xxx/data/file.txt
#[cfg(target_os = "windows")]
fn rewrite_csv_for_vm(
    csv_path: &Path,
    vm_data_dir: &str,
    path_map: &BTreeMap<String, String>,
) -> Result<()> {
    append_desktop_log(&format!(
        "[CSV Rewrite VM] Processing: {}",
        csv_path.display()
    ));
    let mut converted_count = 0;
    let mut missing_count = 0;

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(csv_path)
        .with_context(|| format!("Failed to open CSV: {}", csv_path.display()))?;
    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .from_writer(vec![]);

    for result in reader.records() {
        let record = result
            .with_context(|| format!("Failed to parse CSV record: {}", csv_path.display()))?;
        let mut out: Vec<String> = Vec::with_capacity(record.len());
        for field in record.iter() {
            if looks_like_windows_absolute_path(field) {
                let normalized = normalize_windows_path_str(field);
                if let Some(mapped) = path_map.get(&normalized) {
                    converted_count += 1;
                    out.push(mapped.replace('\\', "/"));
                } else {
                    missing_count += 1;
                    let file_name = Path::new(&normalized)
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy();
                    out.push(format!("{}/{}", vm_data_dir, file_name));
                }
            } else {
                out.push(field.to_string());
            }
        }
        writer
            .write_record(out)
            .with_context(|| format!("Failed to write CSV record: {}", csv_path.display()))?;
    }

    let new_content = String::from_utf8(
        writer
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to finalize CSV writer: {}", e))?,
    )
    .context("Failed to encode rewritten CSV")?;

    if missing_count > 0 {
        append_desktop_log(&format!(
            "[CSV Rewrite VM] Warning: {} paths not in staging map for {}",
            missing_count,
            csv_path.display()
        ));
    }

    append_desktop_log(&format!(
        "[CSV Rewrite VM] Converted {} paths to VM-local in {}",
        converted_count,
        csv_path.display()
    ));
    fs::write(csv_path, new_content)
        .with_context(|| format!("Failed to write converted CSV: {}", csv_path.display()))?;

    Ok(())
}

/// Get unique root directories that need to be mounted (deduplicates nested paths)
#[cfg(target_os = "windows")]
fn get_unique_mount_roots(paths: Vec<PathBuf>) -> Vec<PathBuf> {
    use std::collections::HashSet;

    let mut roots: Vec<PathBuf> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();

    for path in paths {
        // Canonicalize to resolve symlinks and normalize
        let canonical = path.canonicalize().unwrap_or(path);

        // Check if this path is already covered by an existing root
        let mut is_covered = false;
        for root in &roots {
            if canonical.starts_with(root) {
                is_covered = true;
                break;
            }
        }

        if !is_covered && !seen.contains(&canonical) {
            // Remove any existing roots that are children of this new path
            roots.retain(|r| !r.starts_with(&canonical));
            roots.push(canonical.clone());
            seen.insert(canonical);
        }
    }

    roots
}

/// Build a PATH that includes Docker's bin directory (needed for credential helpers on Windows)
#[cfg(target_os = "windows")]
fn build_docker_path(docker_bin: &str) -> Option<String> {
    let docker_path = Path::new(docker_bin);
    let docker_dir = docker_path.parent()?;

    let mut paths = vec![docker_dir.to_path_buf()];
    if let Some(existing) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&existing));
    }

    std::env::join_paths(paths)
        .ok()
        .and_then(|joined| joined.into_string().ok())
}

#[cfg(not(target_os = "windows"))]
fn build_docker_path(_docker_bin: &str) -> Option<String> {
    None
}

fn resolve_docker_config_path() -> Option<String> {
    if let Ok(config) = std::env::var("BIOVAULT_DOCKER_CONFIG") {
        let trimmed = config.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    if let Ok(config) = std::env::var("DOCKER_CONFIG") {
        let trimmed = config.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    None
}

fn apply_docker_config_arg(cmd: &mut Command) {
    if let Some(config) = resolve_docker_config_path() {
        append_desktop_log(&format!("[Pipeline] Using Docker config: {}", config));
        cmd.arg("--config").arg(config);
    }
}

/// Pull a Docker image if not already present (needed on Windows for credential helper PATH issues)
#[cfg(target_os = "windows")]
fn pull_docker_image_if_needed(docker_bin: &str, image: &str) -> Result<()> {
    append_desktop_log(&format!(
        "[Pipeline] Checking if image {} is available...",
        image
    ));

    // Check if image exists locally
    let mut check_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut check_cmd);
    apply_docker_config_arg(&mut check_cmd);
    check_cmd
        .arg("image")
        .arg("inspect")
        .arg(image)
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        check_cmd.env("PATH", &docker_path);
    }

    if check_cmd.status().map(|s| s.success()).unwrap_or(false) {
        append_desktop_log(&format!(
            "[Pipeline] Image {} is already available locally",
            image
        ));
        return Ok(());
    }

    // Pull the image
    append_desktop_log(&format!("[Pipeline] Pulling image {}...", image));
    println!("📦 Pulling Docker image: {}", image);

    let mut pull_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut pull_cmd);
    apply_docker_config_arg(&mut pull_cmd);
    pull_cmd.arg("pull");

    // On Windows, specify platform to avoid manifest resolution issues
    // when Docker daemon might be misconfigured or in Windows container mode.
    // On ARM64, we let Docker auto-select to get native images when available.
    #[cfg(target_os = "windows")]
    if !is_arm64() || should_force_x86_containers() {
        pull_cmd.arg("--platform=linux/amd64");
    }

    pull_cmd.arg(image);

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        pull_cmd.env("PATH", &docker_path);
    }

    let status = pull_cmd.status().context("Failed to execute docker pull")?;

    if !status.success() {
        append_desktop_log(&format!("[Pipeline] Failed to pull image {}", image));
        return Err(anyhow::anyhow!("Failed to pull Docker image: {}", image).into());
    }

    append_desktop_log(&format!("[Pipeline] Successfully pulled image {}", image));
    Ok(())
}

/// Select the appropriate Nextflow runner image based on architecture
/// - ARM64: Use ARM64-native image (avoids Go runtime crashes under QEMU)
/// - x86_64: Use standard image with modern Docker CLI
fn get_nextflow_runner_image() -> &'static str {
    // ARM64 needs native image to avoid Go runtime crashes in Docker CLI
    if is_arm64() {
        // Multi-arch image that includes arm64
        "ghcr.io/openmined/nextflow-runner:25.10.2"
    } else {
        // x86_64 image with modern Docker CLI
        "ghcr.io/openmined/nextflow-runner:25.10.2"
    }
}

#[cfg(target_os = "windows")]
fn ensure_nextflow_runner_image(docker_bin: &str) -> Result<&'static str> {
    // The stock `nextflow/nextflow` image includes an outdated Docker CLI (17.x), which cannot
    // talk to modern Docker daemons. When the workflow uses `container` directives, each task
    // fails with exit code 125 and messages like:
    //   "client version 1.32 is too old. Minimum supported API version is 1.44"
    //
    // We use a pre-built image from GitHub Container Registry that includes a modern docker client.
    // On ARM64, we need native images to avoid Go runtime crashes under QEMU emulation.
    let nextflow_runner_image: &'static str = get_nextflow_runner_image();

    // Check if runner image exists locally
    let mut check_cmd = Command::new(docker_bin);
    super::configure_child_process(&mut check_cmd);
    apply_docker_config_arg(&mut check_cmd);
    check_cmd
        .arg("image")
        .arg("inspect")
        .arg(nextflow_runner_image)
        .stdout(Stdio::null())
        .stderr(Stdio::null());

    if let Some(docker_path) = build_docker_path(docker_bin) {
        check_cmd.env("PATH", &docker_path);
    }

    if check_cmd.status().map(|s| s.success()).unwrap_or(false) {
        append_desktop_log(&format!(
            "[Pipeline] Using Nextflow runner image: {} (arch: {})",
            nextflow_runner_image,
            std::env::consts::ARCH
        ));
        return Ok(nextflow_runner_image);
    }

    // Pull the pre-built image from GitHub Container Registry
    append_desktop_log(&format!(
        "[Pipeline] Pulling Nextflow runner image from {} (arch: {})",
        nextflow_runner_image,
        std::env::consts::ARCH
    ));

    pull_docker_image_if_needed(docker_bin, nextflow_runner_image)?;

    append_desktop_log(&format!(
        "[Pipeline] Nextflow runner image {} ready",
        nextflow_runner_image
    ));
    Ok(nextflow_runner_image)
}

#[cfg(not(target_os = "windows"))]
fn ensure_nextflow_runner_image(docker_bin: &str) -> Result<&'static str> {
    // On non-Windows, if we're using container mode (forced), we need the runner image
    // Otherwise this function shouldn't be called, but return a sensible default
    if should_use_docker_for_nextflow() {
        let nextflow_runner_image: &'static str = get_nextflow_runner_image();

        // Check if image exists locally
        let mut check_cmd = Command::new(docker_bin);
        super::configure_child_process(&mut check_cmd);
        check_cmd
            .arg("image")
            .arg("inspect")
            .arg(nextflow_runner_image)
            .stdout(Stdio::null())
            .stderr(Stdio::null());

        if !check_cmd.status().map(|s| s.success()).unwrap_or(false) {
            // Pull the image
            append_desktop_log(&format!(
                "[Pipeline] Pulling Nextflow runner image: {} (arch: {})",
                nextflow_runner_image,
                std::env::consts::ARCH
            ));
            println!(
                "📦 Pulling Nextflow runner image: {}",
                nextflow_runner_image
            );

            let mut pull_cmd = Command::new(docker_bin);
            super::configure_child_process(&mut pull_cmd);
            let pull_status = pull_cmd
                .args(["pull", nextflow_runner_image])
                .status()
                .context("Failed to pull Nextflow runner image")?;

            if !pull_status.success() {
                return Err(anyhow::anyhow!(
                    "Failed to pull Nextflow runner image: {}",
                    nextflow_runner_image
                )
                .into());
            }
        }

        append_desktop_log(&format!(
            "[Pipeline] Using Nextflow runner image: {} (arch: {})",
            nextflow_runner_image,
            std::env::consts::ARCH
        ));
        return Ok(nextflow_runner_image);
    }

    // Fallback for native Nextflow execution (image not used)
    Ok("nextflow/nextflow:25.10.2")
}

#[cfg(target_os = "linux")]
fn check_docker_running(docker_bin: &str) -> Result<()> {
    let mut cmd = Command::new(docker_bin);
    apply_docker_config_arg(&mut cmd);
    cmd.arg("info")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        cmd.env("PATH", &docker_path);
    }

    let output = cmd
        .output()
        .with_context(|| format!("Failed to execute '{}'", docker_bin))?;

    if output.status.success() {
        return Ok(());
    }

    let stderr = String::from_utf8_lossy(&output.stderr);
    let mut base_msg = format!(
        "Docker daemon is not running ({} exited with {:?}). Please start Docker Desktop or the Docker service and retry.",
        docker_bin,
        output.status.code()
    );

    // Add a more helpful hint when the socket is reachable but permission is denied.
    let lowered = stderr.to_ascii_lowercase();
    if lowered.contains("permission denied") || lowered.contains("connect: permission denied") {
        if std::env::var_os("APPIMAGE").is_some() {
            base_msg.push_str(
                " It looks like this session cannot access /var/run/docker.sock. \
If you recently added your user to the docker group, log out and back in (do not rely on 'newgrp docker' inside the AppImage), then relaunch the app.",
            );
        } else {
            base_msg.push_str(
                " It looks like this session cannot access /var/run/docker.sock. \
If you recently added your user to the docker group, log out and back in, then retry.",
            );
        }
    }

    if !stderr.trim().is_empty() {
        base_msg.push_str(&format!(" Stderr: {}", stderr.trim()));
    }

    Err(anyhow::anyhow!(base_msg).into())
}

#[cfg(not(target_os = "linux"))]
fn check_docker_running(docker_bin: &str) -> Result<()> {
    let mut cmd = Command::new(docker_bin);
    super::configure_child_process(&mut cmd);
    apply_docker_config_arg(&mut cmd);
    cmd.arg("info").stdout(Stdio::null()).stderr(Stdio::null());

    // Add Docker PATH for credential helpers on Windows
    if let Some(docker_path) = build_docker_path(docker_bin) {
        cmd.env("PATH", &docker_path);
    }

    let status = cmd
        .status()
        .with_context(|| format!("Failed to execute '{}'", docker_bin))?;

    if status.success() {
        return Ok(());
    }

    Err(anyhow::anyhow!(
        "Docker daemon is not running ({} exited with {:?}). Please start Docker Desktop or the Docker service and retry.",
        docker_bin,
        status.code()
    )
    .into())
}

/// Check if a container runtime binary is available and working (supports Linux containers).
/// Returns true if `binary info` succeeds.
#[cfg(target_os = "windows")]
fn is_container_runtime_working(binary: &str) -> bool {
    let mut cmd = Command::new(binary);
    super::configure_child_process(&mut cmd);
    cmd.arg("info")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    match cmd.output() {
        Ok(output) => {
            if !output.status.success() {
                return false;
            }
            // Check if the output indicates Linux container support
            // Docker Desktop in Windows container mode shows "OSType: windows"
            // We need "OSType: linux" or Podman (which is always Linux)
            let stdout = String::from_utf8_lossy(&output.stdout);
            // Podman always runs Linux containers, Docker needs to be in Linux mode
            if binary.contains("podman") {
                return true;
            }
            // For Docker, check if it's in Linux container mode
            !stdout.contains("OSType: windows")
        }
        Err(_) => false,
    }
}

/// Get the Podman socket path for mounting into containers.
/// On Windows with WSL, this is typically /run/user/<uid>/podman/podman.sock
#[cfg(target_os = "windows")]
fn get_podman_socket_path() -> Option<String> {
    // Try to get socket path from `podman info`
    let mut cmd = Command::new("podman");
    super::configure_child_process(&mut cmd);
    cmd.args(["info", "--format", "{{.Host.RemoteSocket.Path}}"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    if let Ok(output) = cmd.output() {
        if output.status.success() {
            let socket = String::from_utf8_lossy(&output.stdout).trim().to_string();
            // Strip unix:// prefix if present
            let socket = socket.strip_prefix("unix://").unwrap_or(&socket);
            if !socket.is_empty() {
                append_desktop_log(&format!("[Runtime] Podman socket path: {}", socket));
                return Some(socket.to_string());
            }
        }
    }

    // Fallback to common default path
    append_desktop_log("[Runtime] Using default Podman socket path");
    Some("/run/user/1000/podman/podman.sock".to_string())
}

#[cfg(not(target_os = "windows"))]
fn get_podman_socket_path() -> Option<String> {
    None
}

/// Detect the best available container runtime on Windows.
/// Checks config for explicit preference, then tries Docker, then Podman.
/// Returns the binary path/name for the working runtime.
#[cfg(target_os = "windows")]
fn detect_container_runtime(cfg: Option<&crate::config::Config>) -> Result<String> {
    // Check if user has explicitly configured a container runtime preference
    if let Some(runtime_pref) = cfg.and_then(|c| c.get_binary_path("container_runtime")) {
        let pref = runtime_pref.to_lowercase();
        append_desktop_log(&format!(
            "[Runtime] Config specifies container_runtime: {}",
            pref
        ));

        match pref.as_str() {
            "podman" => {
                // User explicitly wants Podman
                if is_container_runtime_working("podman") {
                    append_desktop_log("[Runtime] Using configured Podman");
                    return Ok("podman".to_string());
                }
                return Err(anyhow::anyhow!(
                    "Container runtime 'podman' is configured but not working. Ensure Podman machine is running."
                ).into());
            }
            "docker" => {
                // User explicitly wants Docker
                if let Some(docker_path) = resolve_binary_path(cfg, "docker") {
                    if is_container_runtime_working(&docker_path) {
                        append_desktop_log("[Runtime] Using configured Docker path");
                        return Ok(docker_path);
                    }
                }
                if is_container_runtime_working("docker") {
                    append_desktop_log("[Runtime] Using system Docker");
                    return Ok("docker".to_string());
                }
                return Err(anyhow::anyhow!(
                    "Container runtime 'docker' is configured but not working. Ensure Docker Desktop is running in Linux container mode."
                ).into());
            }
            "auto" | "" => {
                // Fall through to auto-detection
                append_desktop_log("[Runtime] Auto-detecting container runtime...");
            }
            _ => {
                append_desktop_log(&format!(
                    "[Runtime] Unknown container_runtime '{}', using auto-detection",
                    pref
                ));
            }
        }
    }

    // Auto-detection: try configured Docker path first
    if let Some(docker_path) = resolve_binary_path(cfg, "docker") {
        append_desktop_log(&format!(
            "[Runtime] Checking configured Docker: {}",
            docker_path
        ));
        if is_container_runtime_working(&docker_path) {
            append_desktop_log("[Runtime] Configured Docker is working with Linux containers");
            return Ok(docker_path);
        }
        append_desktop_log("[Runtime] Configured Docker not working or in Windows container mode");
    }

    // Try system Docker
    append_desktop_log("[Runtime] Checking system docker...");
    if is_container_runtime_working("docker") {
        append_desktop_log("[Runtime] System docker is working with Linux containers");
        return Ok("docker".to_string());
    }
    append_desktop_log("[Runtime] System docker not available or not working");

    // Try Podman as fallback
    append_desktop_log("[Runtime] Checking podman...");
    if is_container_runtime_working("podman") {
        append_desktop_log("[Runtime] Podman is working");
        return Ok("podman".to_string());
    }
    append_desktop_log("[Runtime] Podman not available or not working");

    Err(anyhow::anyhow!(
        "No working container runtime found. Please install Docker Desktop (in Linux container mode) or Podman, and ensure the daemon/machine is running."
    ).into())
}

#[cfg(not(target_os = "windows"))]
fn detect_container_runtime(cfg: Option<&crate::config::Config>) -> Result<String> {
    // On non-Windows, just use Docker as before
    Ok(resolve_binary_path(cfg, "docker").unwrap_or_else(|| "docker".to_string()))
}

pub async fn execute_dynamic(
    project_folder: &str,
    args: Vec<String>,
    dry_run: bool,
    resume: bool,
    results_dir: Option<String>,
    run_settings: RunSettings,
) -> Result<()> {
    let project_path = Path::new(project_folder);
    if !project_path.exists() {
        return Err(anyhow::anyhow!("Project folder does not exist: {}", project_folder).into());
    }
    let project_abs = project_path
        .canonicalize()
        .context("Failed to resolve project path")?;

    let nextflow_log_path = project_path.join(".nextflow.log");
    fs::remove_file(&nextflow_log_path).ok();

    let spec_path = crate::project_spec::resolve_project_spec_path(project_path);
    if !spec_path.exists() {
        return Err(anyhow::anyhow!(
            "module.yaml not found in {}. Use 'bv module create' first.",
            project_folder
        )
        .into());
    }

    let spec = ModuleSpec::load(&spec_path)?;
    let runtime = spec.runtime.as_deref().unwrap_or("nextflow");

    match runtime {
        "shell" => {
            return execute_shell(
                &spec,
                project_path,
                args,
                dry_run,
                results_dir,
                run_settings,
            )
            .await;
        }
        "syqure" => {
            return execute_syqure(
                &spec,
                project_path,
                args,
                dry_run,
                results_dir,
                run_settings,
            )
            .await;
        }
        "nextflow" | "dynamic-nextflow" => {
            // Fall through to nextflow execution below
        }
        _ => {
            return Err(anyhow::anyhow!(
                "Unknown runtime '{}'. Supported: nextflow, shell, syqure",
                runtime
            )
            .into());
        }
    }

    let current_datasite = resolve_current_datasite();

    println!("🚀 Running project: {}", spec.name.bold());

    let parsed_args = parse_cli_args(&args)?;
    let nextflow_args = parsed_args.passthrough.clone();
    let nextflow_max_forks = parsed_args
        .nextflow_max_forks
        .or(parse_nextflow_max_forks_env()?);
    if let Some(value) = nextflow_max_forks {
        append_desktop_log(&format!("[Pipeline] Nextflow maxForks set to {}", value));
    }

    validate_no_clashes(&spec, &parsed_args)?;

    let inputs_json = build_inputs_json(&spec, &parsed_args, project_path)?;
    let mut params_json = build_params_json(&spec, &parsed_args)?;

    let assets_dir_path = project_path.join("assets");
    let assets_dir_abs = assets_dir_path
        .canonicalize()
        .unwrap_or_else(|_| assets_dir_path.clone());

    params_json
        .entry("assets_dir".to_string())
        .or_insert_with(|| json!(assets_dir_abs.to_string_lossy().to_string()));

    let results_path = results_dir.as_deref().unwrap_or("results");
    let results_path_buf = if Path::new(results_path).is_absolute() {
        PathBuf::from(results_path)
    } else {
        project_path.join(results_path)
    };
    if !dry_run {
        fs::create_dir_all(&results_path_buf).context("Failed to create results directory")?;
    }
    let results_path_str = results_path_buf.to_string_lossy().to_string();
    #[cfg(target_os = "windows")]
    let requested_results_path_buf = results_path_buf.clone();
    #[cfg(target_os = "windows")]
    let mut results_path_for_run = results_path_buf.clone();
    #[cfg(not(target_os = "windows"))]
    let results_path_for_run = results_path_buf.clone();

    let datasites_override = resolve_datasites_override();
    let template_datasites = if let Some(override_list) = datasites_override {
        override_list
    } else if spec.datasites.as_ref().is_some_and(|d| !d.is_empty()) {
        spec.datasites.clone().unwrap_or_default()
    } else {
        Vec::new()
    };
    let datasite_index = resolve_datasite_index(&template_datasites, current_datasite.as_deref());
    let context = DatasiteContext {
        datasites: template_datasites,
        current: current_datasite.clone(),
        index: datasite_index,
    };
    let env_map = build_shell_env(
        &spec.env,
        &BTreeMap::new(),
        &context,
        project_path,
        &results_path_buf,
        "workflow",
    );

    // Check user workflow exists
    let workflow_path = project_path.join(&spec.workflow);
    if !workflow_path.exists() {
        return Err(anyhow::anyhow!(
            "Workflow file not found: {}. Expected at: {}",
            spec.workflow,
            workflow_path.display()
        )
        .into());
    }

    // Load template from .biovault/env/{template_name}/ (security boundary)
    let biovault_home = crate::config::get_biovault_home()?;
    let biovault_home_abs = biovault_home
        .canonicalize()
        .unwrap_or_else(|_| biovault_home.clone());
    let template_name = match spec.runtime.as_deref() {
        Some("nextflow") | Some("dynamic-nextflow") | None => "dynamic-nextflow",
        Some(other) => other,
    };
    let env_dir = biovault_home_abs.join("env").join(template_name);
    let template_path = env_dir.join("template.nf");

    if matches!(template_name, "dynamic-nextflow" | "nextflow") {
        install_dynamic_template(&biovault_home_abs)?;
    }

    if !template_path.exists() {
        return Err(anyhow::anyhow!(
            "Template not found: {}. Run 'bv init' to install templates.",
            template_path.display()
        )
        .into());
    }

    // Canonicalize paths for Nextflow
    let use_docker = should_use_docker_for_nextflow();
    let template_abs = if use_docker {
        template_path.clone()
    } else {
        template_path
            .canonicalize()
            .context("Failed to resolve template path")?
    };

    let workflow_abs = if use_docker {
        workflow_path.clone()
    } else {
        workflow_path
            .canonicalize()
            .context("Failed to resolve workflow path")?
    };

    let project_spec_abs = if use_docker {
        spec_path.clone()
    } else {
        spec_path
            .canonicalize()
            .context("Failed to resolve project spec path")?
    };

    let inputs_json_str =
        serde_json::to_string(&inputs_json).context("Failed to encode inputs metadata to JSON")?;
    let params_json_str = serde_json::to_string(&params_json)
        .context("Failed to encode parameters metadata to JSON")?;

    let config = crate::config::get_config().ok();
    let nextflow_bin =
        resolve_binary_path(config.as_ref(), "nextflow").unwrap_or_else(|| "nextflow".to_string());

    // Resolve container runtime: check env var first, then auto-detect
    // Supports: BIOVAULT_CONTAINER_RUNTIME=docker|podman or configured docker path
    let docker_bin = if let Ok(runtime) = std::env::var("BIOVAULT_CONTAINER_RUNTIME") {
        let trimmed = runtime.trim();
        if !trimmed.is_empty() {
            append_desktop_log(&format!(
                "[Pipeline] Using container runtime from BIOVAULT_CONTAINER_RUNTIME: {}",
                trimmed
            ));
            trimmed.to_string()
        } else {
            detect_container_runtime(config.as_ref())?
        }
    } else {
        detect_container_runtime(config.as_ref())?
    };

    append_desktop_log(&format!(
        "[Pipeline] Selected container runtime: {}",
        docker_bin
    ));
    println!("🐳 Using container runtime: {}", docker_bin.bold());

    // Check container runtime availability (required for Windows Docker execution and workflow containers)
    // Skip checks in dry-run mode
    if !dry_run && (run_settings.require_docker || use_docker) {
        append_desktop_log("[Pipeline] Checking container runtime availability...");
        if let Err(err) = check_docker_running(&docker_bin) {
            append_desktop_log(&format!(
                "[Pipeline] Container runtime check failed: {}",
                err
            ));
            return Err(err);
        }
        append_desktop_log(&format!(
            "[Pipeline] {} is running (info succeeded)",
            docker_bin
        ));
        println!("✅ {} daemon is running", docker_bin);
    }

    #[cfg(target_os = "windows")]
    let docker_is_podman = is_podman_shim(&docker_bin);
    #[cfg(not(target_os = "windows"))]
    let docker_is_podman = false;

    if docker_is_podman {
        append_desktop_log(
            "[Runtime] Detected Docker binary is a Podman shim; enabling Podman mode",
        );
    }

    // Track Hyper-V mode for post-run cleanup (need to capture before if-else block)
    #[cfg(target_os = "windows")]
    let using_hyperv_mode = use_docker
        && (docker_bin.contains("podman") || docker_is_podman)
        && is_podman_hyperv(
            &docker_bin,
            docker_bin.contains("podman") || docker_is_podman,
        );
    #[cfg(not(target_os = "windows"))]
    let using_hyperv_mode = false;

    #[cfg(target_os = "windows")]
    let hyperv_host_mount = using_hyperv_mode && env_var_truthy("BIOVAULT_HYPERV_MOUNT");
    #[cfg(not(target_os = "windows"))]
    let hyperv_host_mount = false;

    // For Hyper-V mode: capture VM temp dir and results path for post-run copy-back
    #[cfg(target_os = "windows")]
    let mut vm_temp_dir_for_cleanup: Option<String> = None;
    #[cfg(not(target_os = "windows"))]
    let vm_temp_dir_for_cleanup: Option<String> = None;

    // Store docker_results path for post-run copy-back (will be set inside if-else)
    #[cfg(target_os = "windows")]
    let mut vm_results_dir: Option<String> = None;
    #[cfg(not(target_os = "windows"))]
    let vm_results_dir: Option<String> = None;

    #[cfg(target_os = "windows")]
    let mut vm_project_path: Option<String> = None;

    // Track Hyper-V host mount temp dir for optional results copy-back
    #[cfg(target_os = "windows")]
    let mut hyperv_flat_dir: Option<PathBuf> = None;
    #[cfg(target_os = "windows")]
    let mut hyperv_flat_results_dir: Option<PathBuf> = None;
    #[cfg(target_os = "windows")]
    let mut vm_path_map: Option<BTreeMap<String, String>> = None;

    // Build command - use Docker on Windows, native Nextflow elsewhere
    let mut cmd = if use_docker {
        append_desktop_log("[Pipeline] Using Docker to run Nextflow (Windows mode)");

        // Build/get Nextflow runner image with modern Docker CLI
        // Skip image pulling in dry-run mode - use placeholder
        let nextflow_image = if dry_run {
            "ghcr.io/openmined/nextflow-runner:latest"
        } else {
            ensure_nextflow_runner_image(&docker_bin)?
        };

        // Detect if we're using Podman early so we can use correct path format
        // Podman on WSL uses /mnt/c/... while Docker Desktop uses /c/...
        let using_podman = docker_bin.contains("podman") || docker_is_podman;

        // Check if using Podman with Hyper-V backend (which has broken 9P mounts)
        let using_hyperv = using_podman && is_podman_hyperv(&docker_bin, using_podman);
        let hyperv_mount_mode = using_hyperv && hyperv_host_mount;
        if using_hyperv {
            if hyperv_mount_mode {
                append_desktop_log(
                    "[Pipeline] Detected Podman with Hyper-V backend - using host mount mode",
                );
                println!("Using Hyper-V host mount mode (copying inputs to flat dir)");
            } else {
                append_desktop_log(
                    "[Pipeline] Detected Podman with Hyper-V backend - using VM copy mode",
                );
                println!("📋 Using Hyper-V mode (copying files to VM instead of mounting)");
            }
        }

        // For Hyper-V mode: create a temp directory in the VM and track paths
        #[cfg(target_os = "windows")]
        let vm_temp_dir: Option<String> = if using_hyperv && !hyperv_mount_mode && !dry_run {
            let dir = create_vm_temp_dir(&docker_bin, "biovault")?;
            // Capture for post-run cleanup
            vm_temp_dir_for_cleanup = Some(dir.clone());
            Some(dir)
        } else {
            None
        };

        #[cfg(not(target_os = "windows"))]
        let vm_temp_dir: Option<String> = None;

        // Convert all paths to container-compatible format
        // For Hyper-V mode, paths will be in the VM temp directory
        let (
            docker_biovault_home,
            docker_project_path,
            docker_template,
            docker_workflow,
            docker_project_spec,
            docker_results,
            docker_log_path,
        ) = if let Some(ref vm_dir) = vm_temp_dir {
            // Hyper-V mode: use VM-local paths
            let vm_project = format!("{}/project", vm_dir);
            let vm_biovault_home = format!("{}/biovault_home", vm_dir);
            let vm_results = format!("{}/results", vm_dir);

            // Capture for post-run copy-back
            #[cfg(target_os = "windows")]
            {
                vm_results_dir = Some(vm_results.clone());
            }

            // Template/workflow/spec are relative to project when possible
            let template_in_project = is_path_within(&project_abs, &template_abs);
            let template_rel = relative_to_base(&project_abs, &template_abs);
            let workflow_rel = relative_to_base(&project_abs, &workflow_abs);
            let spec_rel = relative_to_base(&project_abs, &project_spec_abs);

            let vm_template = if template_in_project {
                format!("{}/{}", vm_project, template_rel.display()).replace('\\', "/")
            } else {
                format!("{}/env/{}/template.nf", vm_biovault_home, template_name).replace('\\', "/")
            };
            let vm_workflow =
                format!("{}/{}", vm_project, workflow_rel.display()).replace('\\', "/");
            let vm_spec = format!("{}/{}", vm_project, spec_rel.display()).replace('\\', "/");

            #[cfg(target_os = "windows")]
            {
                vm_project_path = Some(vm_project.clone());
            }

            (
                vm_biovault_home,
                vm_project,
                vm_template,
                vm_workflow,
                vm_spec,
                vm_results,
                "/tmp/.nextflow.log".to_string(),
            )
        } else {
            // Normal mode: convert Windows paths to container format
            let docker_biovault_home = windows_path_to_container(&biovault_home_abs, using_podman);
            let docker_project_path = windows_path_to_container(&project_abs, using_podman);
            let docker_template = windows_path_to_container(&template_abs, using_podman);
            let docker_workflow = windows_path_to_container(&workflow_abs, using_podman);
            let docker_project_spec = windows_path_to_container(&project_spec_abs, using_podman);
            let docker_log_path = if using_podman {
                "/tmp/.nextflow.log".to_string()
            } else {
                windows_path_to_container(&nextflow_log_path, using_podman)
            };
            let results_abs = results_path_for_run
                .canonicalize()
                .unwrap_or_else(|_| results_path_for_run.clone());
            let docker_results = windows_path_to_container(&results_abs, using_podman);

            (
                docker_biovault_home,
                docker_project_path,
                docker_template,
                docker_workflow,
                docker_project_spec,
                docker_results,
                docker_log_path,
            )
        };

        let _results_abs = results_path_buf
            .canonicalize()
            .unwrap_or_else(|_| results_path_buf.clone());

        // Extract paths from inputs that need to be mounted (must do before rewriting CSVs)
        // This is Windows-specific: extract paths from CSV files and rewrite them for Docker
        #[cfg(target_os = "windows")]
        let mut inputs_json_value: JsonValue =
            serde_json::to_value(&inputs_json).context("Failed to convert inputs to JSON value")?;
        #[cfg(not(target_os = "windows"))]
        let inputs_json_value: JsonValue =
            serde_json::to_value(&inputs_json).context("Failed to convert inputs to JSON value")?;

        #[cfg(target_os = "windows")]
        if using_hyperv_mode && hyperv_host_mount && !dry_run {
            let host_root = std::env::var("BIOVAULT_HYPERV_HOST_DIR").unwrap_or_else(|_| {
                let drive = std::env::var("SystemDrive").unwrap_or_else(|_| "C:".to_string());
                format!("{}\\bvtemp", drive)
            });
            let unique_id = format!(
                "biovault-{}-{}",
                std::process::id(),
                Local::now().format("%Y%m%d-%H%M%S")
            );
            let flat_root = PathBuf::from(host_root).join(unique_id);
            let flat_data_dir = flat_root.join("data");
            let flat_results_dir = flat_root.join("results");

            fs::create_dir_all(&flat_data_dir).with_context(|| {
                format!(
                    "Failed to create Hyper-V flat data dir: {}",
                    flat_data_dir.display()
                )
            })?;
            fs::create_dir_all(&flat_results_dir).with_context(|| {
                format!(
                    "Failed to create Hyper-V flat results dir: {}",
                    flat_results_dir.display()
                )
            })?;

            append_desktop_log(&format!(
                "[Pipeline] Hyper-V host mount staging dir: {}",
                flat_root.display()
            ));
            println!("Staging inputs in: {}", flat_root.display());

            let mut data_files: Vec<PathBuf> = Vec::new();
            extract_files_from_json(&inputs_json_value, &mut data_files);

            let flat_data_dir_str = flat_data_dir.to_string_lossy().replace('\\', "/");
            let mut path_map: BTreeMap<String, String> = BTreeMap::new();
            let mut stage_pairs: Vec<(PathBuf, PathBuf)> = Vec::new();
            let mut missing_paths: Vec<PathBuf> = Vec::new();

            for data_path in &data_files {
                if !data_path.exists() {
                    missing_paths.push(data_path.clone());
                    continue;
                }
                let key = normalize_windows_path_key(data_path);
                if path_map.contains_key(&key) {
                    append_desktop_log(&format!(
                        "[Pipeline] Duplicate input path (already staged): {}",
                        data_path.display()
                    ));
                    continue;
                }
                let staged_name = stage_name_for_path(data_path);
                let staged_path = flat_data_dir.join(&staged_name);
                let staged_path_str = staged_path.to_string_lossy().replace('\\', "/");
                path_map.insert(key, staged_path_str);
                stage_pairs.push((data_path.clone(), staged_path));
            }

            append_desktop_log(&format!(
                "[Pipeline] Hyper-V staging (flat) inputs: total={}, unique={}, missing={}",
                data_files.len(),
                stage_pairs.len(),
                missing_paths.len()
            ));
            if !path_map.is_empty() {
                let sample = path_map
                    .iter()
                    .take(5)
                    .map(|(k, v)| format!("{} -> {}", k, v))
                    .collect::<Vec<_>>()
                    .join(" | ");
                append_desktop_log(&format!(
                    "[Pipeline] Hyper-V staging map (flat) sample: {}",
                    sample
                ));
            }
            if !missing_paths.is_empty() {
                let sample = missing_paths
                    .iter()
                    .take(5)
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(" | ");
                append_desktop_log(&format!(
                    "[Pipeline] Warning: missing {} input paths (sample): {}",
                    missing_paths.len(),
                    sample
                ));
                if env_var_truthy("BIOVAULT_HYPERV_REQUIRE_ALL_INPUTS") {
                    return Err(anyhow::anyhow!(
                        "Missing {} input paths in Hyper-V host mount staging",
                        missing_paths.len()
                    )
                    .into());
                }
            }

            for (src, dst) in &stage_pairs {
                copy_path_to_path(src, dst)?;
                if dst
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.eq_ignore_ascii_case("csv"))
                    .unwrap_or(false)
                {
                    rewrite_csv_for_flat_dir(dst, &flat_data_dir_str, &path_map)?;
                }
            }

            // Remap inputs JSON to the flat data dir so mounts stay junction-free
            inputs_json_value = remap_json_paths_for_flat_dir_with_map(
                &inputs_json_value,
                &flat_data_dir_str,
                &path_map,
            );

            hyperv_flat_dir = Some(flat_root);
            hyperv_flat_results_dir = Some(flat_results_dir.clone());
            let _ = path_map;
            results_path_for_run = flat_results_dir;
        }

        #[cfg(target_os = "windows")]
        let (mut mount_roots, vm_data_dir) = if let Some(ref vm_dir) = vm_temp_dir {
            // Hyper-V mode: extract data paths and copy them to VM
            let mut data_paths: Vec<PathBuf> = Vec::new();
            extract_files_from_json(&inputs_json_value, &mut data_paths);

            let vm_data = format!("{}/data", vm_dir);
            let mut path_map: BTreeMap<String, String> = BTreeMap::new();
            let mut missing_paths: Vec<PathBuf> = Vec::new();

            for data_path in &data_paths {
                if !data_path.exists() {
                    missing_paths.push(data_path.clone());
                    continue;
                }
                let key = normalize_windows_path_key(data_path);
                if path_map.contains_key(&key) {
                    continue;
                }
                let staged_name = stage_name_for_path(data_path);
                let vm_file_path = format!("{}/{}", vm_data, staged_name);
                path_map.insert(key, vm_file_path);
            }

            append_desktop_log(&format!(
                "[Pipeline] Hyper-V staging (VM) inputs: total={}, unique={}, missing={}",
                data_paths.len(),
                path_map.len(),
                missing_paths.len()
            ));
            if !path_map.is_empty() {
                let sample = path_map
                    .iter()
                    .take(5)
                    .map(|(k, v)| format!("{} -> {}", k, v))
                    .collect::<Vec<_>>()
                    .join(" | ");
                append_desktop_log(&format!(
                    "[Pipeline] Hyper-V staging map (VM) sample: {}",
                    sample
                ));
            }
            if !missing_paths.is_empty() {
                let sample = missing_paths
                    .iter()
                    .take(5)
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(" | ");
                append_desktop_log(&format!(
                    "[Pipeline] Warning: missing {} input paths (sample): {}",
                    missing_paths.len(),
                    sample
                ));
                if env_var_truthy("BIOVAULT_HYPERV_REQUIRE_ALL_INPUTS") {
                    return Err(anyhow::anyhow!(
                        "Missing {} input paths in Hyper-V VM staging",
                        missing_paths.len()
                    )
                    .into());
                }
            }

            // Rewrite CSV files to use VM-local paths before copying them into the VM
            if !dry_run {
                append_desktop_log("[Pipeline] Rewriting CSV files with VM-local paths...");
                rewrite_input_csvs_for_vm(&inputs_json_value, &vm_data, &path_map)?;
            }

            if !dry_run && !data_paths.is_empty() {
                append_desktop_log(&format!(
                    "[Pipeline] Copying {} data files to VM...",
                    data_paths.len()
                ));
                println!("📁 Copying {} data files to VM...", data_paths.len());

                // Create data directory in VM
                let _ = podman_machine_ssh_cmd(&docker_bin)
                    .arg(format!("mkdir -p '{}' && chmod 777 '{}'", vm_data, vm_data))
                    .status();

                // Copy each data file
                for data_path in &data_paths {
                    if data_path.exists() {
                        let key = normalize_windows_path_key(data_path);
                        let vm_file_path = match path_map.get(&key) {
                            Some(path) => path.clone(),
                            None => {
                                append_desktop_log(&format!(
                                    "[Pipeline] Warning: no VM staging path for {}",
                                    data_path.display()
                                ));
                                continue;
                            }
                        };
                        let result = if data_path.is_dir() {
                            copy_dir_to_vm(&docker_bin, data_path, &vm_file_path)
                        } else {
                            copy_file_to_vm(&docker_bin, data_path, &vm_file_path)
                        };
                        if let Err(e) = result {
                            append_desktop_log(&format!(
                                "[Pipeline] Warning: Failed to copy {}: {}",
                                data_path.display(),
                                e
                            ));
                        }
                    }
                }

                // Make all data files readable
                let _ = podman_machine_ssh_cmd(&docker_bin)
                    .arg(format!("chmod -R 777 '{}'", vm_data))
                    .status();
            }

            vm_path_map = Some(path_map);
            (Vec::new(), Some(vm_data))
        } else {
            // Normal mode: extract and mount Windows paths
            let mut data_paths: Vec<PathBuf> = Vec::new();
            extract_paths_from_json(&inputs_json_value, &mut data_paths);
            let roots = get_unique_mount_roots(data_paths);

            // Rewrite CSV files to convert Windows paths to container-compatible paths
            // Skip in dry-run mode since this modifies files
            if !dry_run {
                append_desktop_log(
                    "[Pipeline] Rewriting CSV files with container-compatible paths...",
                );
                append_desktop_log(&format!(
                    "[Pipeline] inputs_json: {}",
                    serde_json::to_string(&inputs_json_value).unwrap_or_default()
                ));
                rewrite_input_csvs_for_container(&inputs_json_value, using_podman)?;
            }

            (roots, None)
        };

        #[cfg(not(target_os = "windows"))]
        let mount_roots: Vec<PathBuf> = Vec::new();

        #[cfg(target_os = "windows")]
        if hyperv_host_mount {
            if let Some(ref flat_results) = hyperv_flat_results_dir {
                mount_roots.push(flat_results.clone());
            }
        }

        #[cfg(target_os = "windows")]
        let docker_results = if using_hyperv_mode && hyperv_host_mount {
            let results_abs = results_path_for_run
                .canonicalize()
                .unwrap_or_else(|_| results_path_for_run.clone());
            windows_path_to_container(&results_abs, using_podman)
        } else {
            docker_results
        };

        append_desktop_log("[Pipeline] Docker path mappings:");
        append_desktop_log(&format!(
            "  biovault_home: {} -> {}",
            biovault_home_abs.display(),
            docker_biovault_home
        ));
        append_desktop_log(&format!(
            "  project_path: {} -> {}",
            project_abs.display(),
            docker_project_path
        ));
        append_desktop_log(&format!(
            "  template: {} -> {}",
            template_abs.display(),
            docker_template
        ));
        append_desktop_log(&format!(
            "[Pipeline] Additional data mounts: {:?}",
            mount_roots
        ));

        // For Hyper-V mode: copy project files to VM before building command
        #[cfg(target_os = "windows")]
        if let Some(ref vm_dir) = vm_temp_dir {
            let vm_project = format!("{}/project", vm_dir);
            let vm_biovault_home = format!("{}/biovault_home", vm_dir);
            append_desktop_log(&format!(
                "[Pipeline] Copying project to VM: {} -> {}",
                project_abs.display(),
                vm_project
            ));
            println!("📁 Copying project files to VM...");
            copy_dir_to_vm(&docker_bin, &project_abs, &vm_project)?;

            // If the template lives outside the project, copy it into the VM
            let template_in_project = is_path_within(&project_abs, &template_abs);
            if !template_in_project {
                let template_src = biovault_home_abs.join("env").join(template_name);
                if template_src.exists() {
                    let vm_template_dir = format!("{}/env/{}", vm_biovault_home, template_name);
                    append_desktop_log(&format!(
                        "[Pipeline] Copying template env to VM: {} -> {}",
                        template_src.display(),
                        vm_template_dir
                    ));
                    copy_dir_to_vm(&docker_bin, &template_src, &vm_template_dir)?;
                }
            }

            // Create results directory in VM
            let vm_results = format!("{}/results", vm_dir);
            let _ = podman_machine_ssh_cmd(&docker_bin)
                .arg(format!(
                    "mkdir -p '{}' && chmod 777 '{}'",
                    vm_results, vm_results
                ))
                .status();

            append_desktop_log(&format!(
                "[Pipeline] VM directories created: project={}, results={}",
                vm_project, vm_results
            ));
        }

        let mut docker_cmd = Command::new(&docker_bin);
        super::configure_child_process(&mut docker_cmd);
        apply_docker_config_arg(&mut docker_cmd);

        // Add Docker PATH for credential helpers on Windows
        if let Some(docker_path) = build_docker_path(&docker_bin) {
            docker_cmd.env("PATH", docker_path);
        }

        docker_cmd.arg("run").arg("--rm");

        // Add container name for graceful stop support
        // Use BIOVAULT_FLOW_RUN_ID if available, otherwise generate a timestamp
        let container_name = std::env::var("BIOVAULT_FLOW_RUN_ID")
            .unwrap_or_else(|_| chrono::Utc::now().format("%Y%m%d%H%M%S").to_string());
        let container_name = format!("bv-nf-{}", container_name);
        docker_cmd.args(["--name", &container_name]);

        // Save container name to file for pause/stop to use
        if let Some(ref results) = results_dir {
            let container_file = Path::new(results).join("flow.container");
            let _ = std::fs::write(&container_file, &container_name);
            append_desktop_log(&format!(
                "[Pipeline] Nextflow container name: {} (saved to {})",
                container_name,
                container_file.display()
            ));
        }

        // On Windows with Docker (not Podman), specify platform.
        // On ARM64, we don't force amd64 since we use ARM64-native images for the runner.
        #[cfg(target_os = "windows")]
        if !using_podman && (!is_arm64() || should_force_x86_containers()) {
            docker_cmd.args(["--platform", "linux/amd64"]);
        }

        // When using Podman, add flags to fix permission issues with mounted volumes
        if using_podman {
            // Handle .nextflow directory permissions
            // If switching from Docker to Podman, the old .nextflow dir may have wrong ownership
            let nextflow_local_dir = project_abs.join(".nextflow");
            if nextflow_local_dir.exists() {
                // Check if we can write to it by trying to create a test file
                let test_file = nextflow_local_dir.join(".podman-write-test");
                match fs::write(&test_file, "test") {
                    Ok(_) => {
                        let _ = fs::remove_file(&test_file);
                    }
                    Err(_) => {
                        // Can't write - likely created by Docker with different permissions
                        // Remove and recreate with correct ownership
                        append_desktop_log(&format!(
                            "[Pipeline] Removing .nextflow directory with incompatible permissions: {}",
                            nextflow_local_dir.display()
                        ));
                        println!(
                            "⚠️  Cleaning .nextflow directory (was created by different container runtime)"
                        );
                        let _ = fs::remove_dir_all(&nextflow_local_dir);
                    }
                }
            }

            // Create .nextflow directory if it doesn't exist
            if !nextflow_local_dir.exists() {
                let _ = fs::create_dir_all(&nextflow_local_dir);
                append_desktop_log(&format!(
                    "[Pipeline] Created .nextflow directory: {}",
                    nextflow_local_dir.display()
                ));
            }

            docker_cmd
                .arg("--userns=keep-id")
                // Disable SELinux labeling which can cause permission issues on mounted volumes
                .arg("--security-opt")
                .arg("label=disable")
                // Set NXF_HOME to a writable location inside the container
                .arg("-e")
                .arg("NXF_HOME=/tmp/.nextflow")
                // Set NXF_TEMP to avoid writing temp files to Windows-mounted paths
                // (Windows mounts don't support POSIX operations like mv with attribute preservation)
                .arg("-e")
                .arg("NXF_TEMP=/tmp");

            // Mount Podman socket and set CONTAINER_HOST for nested container execution
            if let Some(podman_socket) = get_podman_socket_path() {
                append_desktop_log(&format!(
                    "[Pipeline] Mounting Podman socket: {} -> /run/podman/podman.sock",
                    podman_socket
                ));
                docker_cmd
                    .arg("-v")
                    .arg(format!("{}:/run/podman/podman.sock", podman_socket))
                    .arg("-e")
                    .arg("CONTAINER_HOST=unix:///run/podman/podman.sock");
            }
        } else {
            // Mount Docker Desktop engine socket so Nextflow-in-Docker can launch workflow containers.
            //
            // With Docker Desktop (Linux containers), the daemon runs in a Linux VM and provides a Unix socket.
            // Binding it into the Nextflow container allows Nextflow to launch workflow containers.
            docker_cmd
                .arg("-v")
                .arg("/var/run/docker.sock:/var/run/docker.sock");
        }

        // For Hyper-V mode: we don't mount Windows paths (they're copied to VM)
        // For normal mode: mount Windows paths into container
        if let Some(ref vm_dir) = vm_temp_dir {
            // Hyper-V mode: mount the VM temp directory into the container
            // Files were copied to the VM's native filesystem, now mount them into container
            append_desktop_log(&format!(
                "[Pipeline] Hyper-V mode: mounting VM directory {} into container",
                vm_dir
            ));
            docker_cmd.arg("-v").arg(format!("{}:{}", vm_dir, vm_dir));
        } else {
            // Normal mode: mount Windows paths
            docker_cmd
                // Mount the BioVault home directory
                .arg("-v")
                .arg(format!(
                    "{}:{}",
                    windows_path_to_mount_source(&biovault_home, using_podman),
                    docker_biovault_home
                ))
                // Mount the project path (may be same as above, Docker handles duplicates)
                .arg("-v")
                .arg(format!(
                    "{}:{}",
                    windows_path_to_mount_source(project_path, using_podman),
                    docker_project_path
                ));

            // Mount additional data directories discovered from inputs
            for mount_path in &mount_roots {
                let container_mount = windows_path_to_container(mount_path, using_podman);
                let host_mount = windows_path_to_mount_source(mount_path, using_podman);
                append_desktop_log(&format!(
                    "[Pipeline] Adding mount: {} -> {}",
                    mount_path.display(),
                    container_mount
                ));
                docker_cmd
                    .arg("-v")
                    .arg(format!("{}:{}", host_mount, container_mount));
            }
        }

        // Generate runtime-specific Nextflow config (Docker vs Podman)
        let runtime_config_path = if !dry_run {
            // On Windows with Podman, prefer copying inputs to avoid broken symlinks on /mnt/c mounts.
            let stage_in_copy = hyperv_host_mount || (cfg!(target_os = "windows") && using_podman);
            let config_path = generate_runtime_config(
                &project_abs,
                using_podman,
                stage_in_copy,
                nextflow_max_forks,
            )?;
            // For Hyper-V mode, the config file was copied to VM along with project
            // Reference it via the VM project path
            #[cfg(target_os = "windows")]
            if let Some(ref vm_dir) = vm_temp_dir {
                // Copy the config file to VM if it wasn't part of the project copy
                let config_name = config_path
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy();
                let vm_config_path = format!("{}/project/{}", vm_dir, config_name);
                if let Err(e) = copy_file_to_vm(&docker_bin, &config_path, &vm_config_path) {
                    append_desktop_log(&format!(
                        "[Pipeline] Warning: Failed to copy config to VM: {}",
                        e
                    ));
                }
                Some(vm_config_path)
            } else {
                Some(windows_path_to_container(&config_path, using_podman))
            }
            #[cfg(not(target_os = "windows"))]
            Some(config_path.to_string_lossy().to_string())
        } else {
            None
        };

        // Use /tmp as working directory for Podman to avoid Windows mount I/O issues
        // All paths are absolute so this is safe. /tmp always exists in containers.
        let container_work_dir = if using_podman {
            "/tmp".to_string()
        } else {
            docker_project_path.clone()
        };

        docker_cmd
            // Set working directory
            .arg("-w")
            .arg(&container_work_dir)
            // Use Nextflow runner image with modern Docker CLI
            .arg(nextflow_image)
            // Nextflow command (container entrypoint is bash, not nextflow)
            .arg("nextflow")
            // Nextflow arguments
            .arg("-log")
            .arg(&docker_log_path)
            .arg("run");

        // Pass runtime config file to Nextflow
        if let Some(ref config_docker_path) = runtime_config_path {
            docker_cmd.arg("-c").arg(config_docker_path);
        }

        docker_cmd.arg(&docker_template);

        // For Podman: work directory must be on Windows mount so nested containers can access it
        // We use /tmp as working directory (-w) but explicit -work-dir on the mount
        if using_podman {
            let work_dir_path = format!("{}/work", docker_project_path);
            docker_cmd.arg("-work-dir").arg(&work_dir_path);
        }

        if resume {
            docker_cmd.arg("-resume");
        }

        for extra in &nextflow_args {
            docker_cmd.arg(extra);
        }

        // Re-encode JSON with container paths (inputs_json_value already created above)
        let params_json_value: JsonValue =
            serde_json::to_value(&params_json).context("Failed to convert params to JSON value")?;

        // For Hyper-V mode: use VM-local paths. For normal mode: use container mount paths.
        #[cfg(target_os = "windows")]
        let (docker_inputs_json, mut docker_params_json) = if let Some(ref vm_data) = vm_data_dir {
            if let Some(ref map) = vm_path_map {
                (
                    remap_json_paths_for_vm_with_map(&inputs_json_value, vm_data, map),
                    remap_json_paths_for_vm_with_map(&params_json_value, vm_data, map),
                )
            } else {
                (
                    remap_json_paths_for_vm(&inputs_json_value, vm_data),
                    remap_json_paths_for_vm(&params_json_value, vm_data),
                )
            }
        } else {
            (
                remap_json_paths_for_container(&inputs_json_value, using_podman),
                remap_json_paths_for_container(&params_json_value, using_podman),
            )
        };

        #[cfg(not(target_os = "windows"))]
        let (docker_inputs_json, docker_params_json) = (
            remap_json_paths_for_container(&inputs_json_value, using_podman),
            remap_json_paths_for_container(&params_json_value, using_podman),
        );

        #[cfg(target_os = "windows")]
        if vm_data_dir.is_some() {
            if let Some(ref vm_project) = vm_project_path {
                if let JsonValue::Object(ref mut map) = docker_params_json {
                    map.insert(
                        "assets_dir".to_string(),
                        JsonValue::String(format!("{}/assets", vm_project)),
                    );
                }
            }
        }

        let docker_inputs_json_str = serde_json::to_string(&docker_inputs_json)
            .context("Failed to encode Docker inputs metadata to JSON")?;
        let docker_params_json_str = serde_json::to_string(&docker_params_json)
            .context("Failed to encode Docker parameters metadata to JSON")?;

        docker_cmd
            .arg("--work_flow_file")
            .arg(&docker_workflow)
            .arg("--module_spec")
            .arg(&docker_project_spec)
            .arg("--inputs_json")
            .arg(docker_inputs_json_str)
            .arg("--params_json")
            .arg(docker_params_json_str)
            .arg("--results_dir")
            .arg(&docker_results);

        docker_cmd
    } else {
        // Native Nextflow execution (macOS, Linux)
        append_desktop_log(&format!(
            "[Pipeline] Using native nextflow binary: {}",
            nextflow_bin
        ));

        // Generate runtime config for native Docker execution
        // This config sets docker.enabled = true so Nextflow uses Docker for process containers
        let native_runtime_config = if !dry_run && run_settings.require_docker {
            Some(generate_runtime_config(
                &project_abs,
                false,
                false,
                nextflow_max_forks,
            )?)
        } else {
            None
        };

        // Log original PATH from environment
        if let Some(original_path) = std::env::var_os("PATH") {
            append_desktop_log(&format!(
                "[Pipeline] Original PATH from environment: {}",
                original_path.to_string_lossy()
            ));
        } else {
            append_desktop_log("[Pipeline] WARNING: No PATH environment variable found!");
        }

        let mut native_cmd = Command::new(&nextflow_bin);
        super::configure_child_process(&mut native_cmd);

        append_desktop_log("[Pipeline] Preferred binary paths:");
        for binary in ["nextflow", "java", "docker"] {
            if let Some(path) = resolve_binary_path(config.as_ref(), binary) {
                append_desktop_log(&format!("  {} = {}", binary, path));
            } else {
                append_desktop_log(&format!("  {} = <not configured>", binary));
            }
        }

        if let Some(path_env) = build_augmented_path(config.as_ref()) {
            append_desktop_log(&format!(
                "[Pipeline] Final augmented PATH for nextflow: {}",
                path_env
            ));
            native_cmd.env("PATH", path_env);
        } else {
            append_desktop_log(
                "[Pipeline] WARNING: Could not build augmented PATH, using system PATH",
            );
        }

        native_cmd.arg("-log").arg(&nextflow_log_path);

        native_cmd.arg("run");

        // Pass runtime config file to Nextflow (enables docker.enabled = true)
        if let Some(ref config_path) = native_runtime_config {
            native_cmd.arg("-c").arg(config_path);
        }

        native_cmd.arg(&template_abs);

        if resume {
            native_cmd.arg("-resume");
        }

        for extra in &nextflow_args {
            native_cmd.arg(extra);
        }

        native_cmd
            .arg("--work_flow_file")
            .arg(&workflow_abs)
            .arg("--module_spec")
            .arg(&project_spec_abs)
            .arg("--inputs_json")
            .arg(inputs_json_str)
            .arg("--params_json")
            .arg(params_json_str)
            .arg("--results_dir")
            .arg(&results_path_str);

        native_cmd
    };

    for (key, value) in &env_map {
        cmd.env(key, value);
    }

    let display_cmd = format_command(&cmd);

    if dry_run {
        println!("\n🔍 Dry run - would execute:");
        println!("  {}", display_cmd.dimmed());
        append_desktop_log(&format!(
            "[Pipeline] (dry-run) Nextflow command: {}",
            display_cmd
        ));
        return Ok(());
    }

    println!("\n▶️  Executing Nextflow...\n");
    println!("  {}", display_cmd.dimmed());
    append_desktop_log(&format!("[Pipeline] Nextflow command: {}", display_cmd));

    cmd.current_dir(project_path);
    let work_dir = project_path.join("work");
    let status = execute_with_logging(
        cmd,
        Some(nextflow_log_path),
        Some(work_dir),
        Some(project_path.to_path_buf()),
    )
    .context("Failed to execute nextflow")?;

    // For Hyper-V mode: copy results back from VM and cleanup
    #[cfg(target_os = "windows")]
    if using_hyperv_mode {
        if let (Some(ref vm_results), Some(ref vm_dir)) =
            (&vm_results_dir, &vm_temp_dir_for_cleanup)
        {
            append_desktop_log(&format!(
                "[Pipeline] Copying results back from VM: {} -> {}",
                vm_results,
                results_path_buf.display()
            ));
            println!("📁 Copying results back from VM...");

            if let Err(e) = copy_from_vm(&docker_bin, vm_results, &results_path_buf) {
                append_desktop_log(&format!(
                    "[Pipeline] Warning: Failed to copy results from VM: {}",
                    e
                ));
                eprintln!("⚠️  Warning: Failed to copy results from VM: {}", e);
            }

            // Cleanup VM temp directory
            append_desktop_log(&format!(
                "[Pipeline] Cleaning up VM temp directory: {}",
                vm_dir
            ));
            cleanup_vm_dir(&docker_bin, vm_dir);
        }
    }

    // For Hyper-V host mount mode: copy results back from flat dir
    #[cfg(target_os = "windows")]
    if using_hyperv_mode && hyperv_host_mount {
        if let Some(ref flat_results) = hyperv_flat_results_dir {
            if flat_results != &requested_results_path_buf {
                append_desktop_log(&format!(
                    "[Pipeline] Copying results back from flat dir: {} -> {}",
                    flat_results.display(),
                    requested_results_path_buf.display()
                ));
                println!("Copying results back from flat dir...");
                if let Err(e) = copy_dir_recursive(flat_results, &requested_results_path_buf) {
                    append_desktop_log(&format!(
                        "[Pipeline] Warning: Failed to copy results from flat dir: {}",
                        e
                    ));
                    eprintln!("Warning: Failed to copy results from flat dir: {}", e);
                }
            }
        }

        if !env_var_truthy("BIOVAULT_KEEP_HYPERV_HOST_DIR") {
            if let Some(ref flat_root) = hyperv_flat_dir {
                let _ = fs::remove_dir_all(flat_root);
            }
        } else if let Some(ref flat_root) = hyperv_flat_dir {
            println!("Keeping Hyper-V host mount dir: {}", flat_root.display());
        }
    }

    // Suppress unused variable warnings on non-Windows
    let _ = &vm_temp_dir_for_cleanup;
    let _ = &vm_results_dir;
    let _ = &using_hyperv_mode;

    if !status.success() {
        append_desktop_log(&format!(
            "[Pipeline] Nextflow exited with status: {:?}",
            status.code()
        ));
        return Err(
            anyhow::anyhow!("Nextflow execution failed with code: {:?}", status.code()).into(),
        );
    }

    println!("\n✅ Workflow completed successfully!");
    append_desktop_log("[Pipeline] Workflow completed successfully!");
    Ok(())
}

#[derive(Debug, Clone)]
struct DatasiteContext {
    datasites: Vec<String>,
    current: Option<String>,
    index: Option<usize>,
}

fn resolve_current_datasite() -> Option<String> {
    if let Some(current) = execution_context().and_then(|ctx| ctx.current_datasite) {
        if !current.trim().is_empty() {
            return Some(current);
        }
    }
    if let Ok(env_override) = env::var("BIOVAULT_DATASITE_OVERRIDE") {
        if !env_override.trim().is_empty() {
            return Some(env_override.trim().to_string());
        }
    }
    if let Ok(cfg) = crate::config::Config::load() {
        if !cfg.email.trim().is_empty() {
            return Some(cfg.email);
        }
    }
    if let Ok(env_email) = env::var("SYFTBOX_EMAIL") {
        if !env_email.trim().is_empty() {
            return Some(env_email.trim().to_string());
        }
    }
    if let Ok(env_email) = env::var("BIOVAULT_DATASITE") {
        if !env_email.trim().is_empty() {
            return Some(env_email.trim().to_string());
        }
    }
    None
}

fn resolve_datasites_override() -> Option<Vec<String>> {
    if let Some(list) = execution_context().and_then(|ctx| ctx.datasites_override) {
        if !list.is_empty() {
            return Some(list);
        }
    }
    let raw = env::var("BIOVAULT_DATASITES_OVERRIDE").ok()?;
    let list: Vec<String> = raw
        .split(',')
        .map(|entry| entry.trim())
        .filter(|entry| !entry.is_empty())
        .map(|entry| entry.to_string())
        .collect();
    if list.is_empty() {
        None
    } else {
        Some(list)
    }
}

fn resolve_datasite_index(datasites: &[String], current: Option<&str>) -> Option<usize> {
    let current = current?;
    datasites.iter().position(|site| site == current)
}

fn render_template(value: &str, ctx: &DatasiteContext) -> String {
    let mut rendered = value.to_string();
    if let Some(current) = &ctx.current {
        rendered = rendered.replace("{current_datasite}", current);
    }
    if let Some(index) = ctx.index {
        rendered = rendered.replace("{datasites.index}", &index.to_string());
        rendered = rendered.replace("{datasite.index}", &index.to_string());
    }
    rendered.replace("{datasites}", &ctx.datasites.join(","))
}

fn env_key_suffix(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_uppercase()
            } else {
                '_'
            }
        })
        .collect()
}

fn build_shell_env(
    spec_env: &BTreeMap<String, String>,
    step_env: &BTreeMap<String, String>,
    ctx: &DatasiteContext,
    project_path: &Path,
    results_path: &Path,
    step_id: &str,
) -> BTreeMap<String, String> {
    let mut env_map = spec_env.clone();
    for (key, value) in step_env {
        env_map.insert(key.clone(), value.clone());
    }

    let mut rendered = BTreeMap::new();
    for (key, value) in env_map {
        rendered.insert(key, render_template(&value, ctx));
    }

    if let Ok(bv_bin) = env::var("BV_BIN") {
        if !bv_bin.trim().is_empty() {
            rendered.insert("BV_BIN".to_string(), bv_bin);
        }
    }

    rendered.insert(
        "BV_PROJECT_DIR".to_string(),
        project_path.to_string_lossy().to_string(),
    );
    rendered.insert(
        "BV_RESULTS_DIR".to_string(),
        results_path.to_string_lossy().to_string(),
    );
    rendered.insert(
        "BV_ASSETS_DIR".to_string(),
        project_path.join("assets").to_string_lossy().to_string(),
    );
    rendered.insert("BV_STEP_ID".to_string(), step_id.to_string());
    rendered.insert("BV_DATASITES".to_string(), ctx.datasites.join(","));

    if let Some(current) = &ctx.current {
        rendered.insert("BV_CURRENT_DATASITE".to_string(), current.clone());
    }
    if let Some(index) = ctx.index {
        rendered.insert("BV_DATASITE_INDEX".to_string(), index.to_string());
    }

    rendered
}

fn shell_passthrough_args(args: &[String]) -> Vec<String> {
    args.iter()
        .position(|arg| arg == "--")
        .map(|idx| args[idx + 1..].to_vec())
        .unwrap_or_default()
}

fn default_shell_step(spec: &ModuleSpec) -> ModuleStepSpec {
    ModuleStepSpec {
        id: "run".to_string(),
        foreach: spec.datasites.clone(),
        order: None,
        env: BTreeMap::new(),
        inputs: Vec::new(),
        outputs: Vec::new(),
    }
}

fn shell_input_uses_path_resolution(raw_type: &str) -> bool {
    let trimmed = raw_type.trim();
    let normalized = trimmed.trim_end_matches('?');
    matches!(normalized, "File" | "Directory" | "Dir" | "Path")
}

async fn execute_shell(
    spec: &ModuleSpec,
    project_path: &Path,
    args: Vec<String>,
    dry_run: bool,
    results_dir: Option<String>,
    _run_settings: RunSettings,
) -> Result<()> {
    let workflow_path = project_path.join(&spec.workflow);
    if !workflow_path.exists() {
        return Err(anyhow::anyhow!("Workflow file not found: {}", workflow_path.display()).into());
    }

    let results_path = if let Some(dir) = results_dir {
        let path = PathBuf::from(&dir);
        if path.is_absolute() {
            path
        } else {
            project_path.join(dir)
        }
    } else {
        project_path.join("results")
    };
    if !dry_run {
        fs::create_dir_all(&results_path)?;
    }

    let current_datasite = resolve_current_datasite();
    if spec.datasites.is_some() && current_datasite.is_none() {
        return Err(anyhow::anyhow!(
            "datasites specified in project.yaml but current datasite could not be determined"
        )
        .into());
    }

    let parsed_args = parse_cli_args(&args)?;
    let steps = if spec.steps.is_empty() {
        vec![default_shell_step(spec)]
    } else {
        spec.steps.clone()
    };

    for step in steps {
        let passthrough_args = shell_passthrough_args(&args);
        let step_targets = step
            .foreach
            .clone()
            .or_else(|| spec.datasites.clone())
            .unwrap_or_default();

        let should_run = if step_targets.is_empty() {
            true
        } else if let Some(current) = current_datasite.as_ref() {
            step_targets.iter().any(|site| site == current)
        } else {
            false
        };

        if !should_run {
            println!(
                "⏭️  Skipping step '{}' (current datasite not in foreach list)",
                step.id
            );
            continue;
        }

        let template_datasites = if let Some(override_list) = resolve_datasites_override() {
            override_list
        } else if let Some(spec_datasites) = &spec.datasites {
            if spec_datasites.is_empty() {
                step_targets.clone()
            } else {
                spec_datasites.clone()
            }
        } else {
            step_targets.clone()
        };
        let index = resolve_datasite_index(&template_datasites, current_datasite.as_deref());

        let ctx = DatasiteContext {
            datasites: template_datasites,
            current: current_datasite.clone(),
            index,
        };

        let mut env_map = build_shell_env(
            &spec.env,
            &step.env,
            &ctx,
            project_path,
            &results_path,
            &step.id,
        );
        if !env_map.contains_key("BV_RUN_ID") {
            if let Ok(Some(hash)) = module_yaml_hash(project_path) {
                env_map.insert("BV_RUN_ID".to_string(), hash);
            }
        }
        if let Ok(cfg) = crate::config::Config::load() {
            if let Ok(data_dir) = cfg.get_syftbox_data_dir() {
                env_map
                    .entry("BV_SYFTBOX_DATA_DIR".to_string())
                    .or_insert_with(|| data_dir.to_string_lossy().to_string());
                let datasites_root = data_dir.join("datasites");
                env_map
                    .entry("BV_DATASITES_ROOT".to_string())
                    .or_insert_with(|| datasites_root.to_string_lossy().to_string());
            }
        }
        let input_specs = if step.inputs.is_empty() {
            &spec.inputs
        } else {
            &step.inputs
        };
        for input in input_specs {
            let env_key = format!("BV_INPUT_{}", env_key_suffix(&input.name));
            if let Some(arg) = parsed_args.inputs.get(&input.name) {
                let rendered_value = render_template(&arg.value, &ctx);
                let final_value = if shell_input_uses_path_resolution(&input.raw_type) {
                    let input_path = if Path::new(&rendered_value).is_absolute() {
                        PathBuf::from(rendered_value)
                    } else {
                        project_path.join(rendered_value)
                    };
                    input_path.to_string_lossy().to_string()
                } else {
                    rendered_value
                };
                env_map.insert(env_key, final_value);
                continue;
            }
            if let Some(path_template) = input.path.as_deref() {
                let rendered_path = render_template(path_template, &ctx);
                let input_path = if Path::new(&rendered_path).is_absolute() {
                    PathBuf::from(rendered_path)
                } else {
                    project_path.join(rendered_path)
                };
                env_map
                    .entry(env_key)
                    .or_insert_with(|| input_path.to_string_lossy().to_string());
            }
        }
        let output_specs = if step.outputs.is_empty() {
            &spec.outputs
        } else {
            &step.outputs
        };
        for output in output_specs {
            let raw_path = output.path.as_deref().unwrap_or(&output.name);
            let rendered_path = render_template(raw_path, &ctx);
            let output_path = if Path::new(&rendered_path).is_absolute() {
                PathBuf::from(rendered_path)
            } else {
                results_path.join(rendered_path)
            };
            let env_key = format!("BV_OUTPUT_{}", env_key_suffix(&output.name));
            env_map.insert(env_key, output_path.to_string_lossy().to_string());
        }

        let mut cmd = Command::new(resolve_git_bash());
        super::configure_child_process(&mut cmd);
        // Convert Windows path to Unix-style for Git Bash
        #[cfg(target_os = "windows")]
        let bash_workflow_path = windows_path_to_container(&workflow_path, false);
        #[cfg(not(target_os = "windows"))]
        let bash_workflow_path = workflow_path.to_string_lossy().to_string();
        cmd.arg(&bash_workflow_path);
        if !passthrough_args.is_empty() {
            cmd.args(&passthrough_args);
        }
        cmd.current_dir(project_path);
        for (key, value) in &env_map {
            cmd.env(key, value);
        }

        let display_cmd = format_command(&cmd);
        if dry_run {
            println!("\n🔍 Dry run - would execute:");
            println!("  {}", display_cmd.dimmed());
            continue;
        }

        println!("\n▶️  Executing shell workflow for step '{}'...", step.id);
        println!("  {}", display_cmd.dimmed());
        let status = cmd.status().context("Failed to execute shell workflow")?;
        if !status.success() {
            return Err(
                anyhow::anyhow!("Shell workflow exited with code: {:?}", status.code()).into(),
            );
        }
    }

    println!("\n✅ Shell workflow completed successfully!");
    Ok(())
}

async fn execute_syqure(
    spec: &ModuleSpec,
    module_path: &Path,
    args: Vec<String>,
    dry_run: bool,
    results_dir: Option<String>,
    _run_settings: RunSettings,
) -> Result<()> {
    println!("🔐 Running syqure module: {}", spec.name.bold());

    let current_datasite = resolve_current_datasite();

    let results_path = if let Some(dir) = results_dir {
        let path = PathBuf::from(&dir);
        if path.is_absolute() {
            path
        } else {
            module_path.join(dir)
        }
    } else {
        module_path.join("results")
    };
    if !dry_run {
        fs::create_dir_all(&results_path)?;
    }

    let (syqure_binary, use_docker) = resolve_syqure_backend(spec)?;

    let exec_ctx = execution_context();
    let run_id = exec_ctx
        .as_ref()
        .and_then(|c| c.run_id.clone())
        .or_else(|| env::var("BV_RUN_ID").ok())
        .unwrap_or_else(|| chrono::Local::now().format("%Y%m%d%H%M%S").to_string());

    let datasites_root = exec_ctx
        .as_ref()
        .and_then(|c| c.syftbox_data_dir.as_ref())
        .map(|dir| {
            PathBuf::from(dir)
                .join("datasites")
                .to_string_lossy()
                .to_string()
        })
        .or_else(|| env::var("BV_DATASITES_ROOT").ok())
        .or_else(|| {
            env::var("SYFTBOX_DATA_DIR").ok().map(|dir| {
                // Use PathBuf::join to properly handle path separators on all platforms
                PathBuf::from(&dir)
                    .join("datasites")
                    .to_string_lossy()
                    .to_string()
            })
        })
        .ok_or_else(|| {
            anyhow::anyhow!(
                "Syqure runtime requires BV_DATASITES_ROOT or SYFTBOX_DATA_DIR environment variable"
            )
        })?;

    let entrypoint = spec
        .runner
        .as_ref()
        .and_then(|r| r.entrypoint.as_ref())
        .cloned()
        .unwrap_or_else(|| spec.workflow.clone());

    let entrypoint_path = module_path.join(&entrypoint);
    if !entrypoint_path.exists() {
        return Err(
            anyhow::anyhow!("Syqure entrypoint not found: {}", entrypoint_path.display()).into(),
        );
    }

    let poll_ms = spec
        .runner
        .as_ref()
        .and_then(|r| r.syqure.as_ref())
        .and_then(|s| s.poll_ms)
        .unwrap_or(50);

    let raw_transport = spec
        .runner
        .as_ref()
        .and_then(|r| r.syqure.as_ref())
        .and_then(|s| s.transport.as_ref())
        .cloned()
        .unwrap_or_else(|| "file".to_string());
    let transport = canonical_hotlink_transport(&raw_transport);

    // Get platform from spec, defaulting to linux/amd64
    // Syqure containers are x86-only for now, so we always force amd64
    // This will use QEMU emulation on ARM64 systems
    let docker_platform = spec
        .runner
        .as_ref()
        .and_then(|r| r.syqure.as_ref())
        .and_then(|s| s.platform.as_ref())
        .cloned()
        .unwrap_or_else(|| "linux/amd64".to_string());

    let parsed_args = parse_cli_args(&args)?;

    // Build datasite context - same pattern as execute_shell
    let datasites_list = resolve_datasites_override()
        .or_else(|| spec.datasites.clone())
        .unwrap_or_default();
    let datasite_index_ctx = resolve_datasite_index(&datasites_list, current_datasite.as_deref());
    let ctx = DatasiteContext {
        datasites: datasites_list.clone(),
        current: current_datasite.clone(),
        index: datasite_index_ctx,
    };

    // MPC party info - derived from context, not env vars
    let party_emails = if datasites_list.is_empty() {
        env::var("BV_DATASITES").unwrap_or_default()
    } else {
        datasites_list.join(",")
    };

    if party_emails.is_empty() {
        return Err(anyhow::anyhow!(
            "Syqure runtime requires datasites in module spec or BV_DATASITES environment variable"
        )
        .into());
    }

    let party_count = party_emails.split(',').count();
    let current_email = current_datasite.clone().unwrap_or_default();
    let party_id = datasite_index_ctx.unwrap_or(0);
    let syqure_port_base = prepare_syqure_port_base_for_run(
        &run_id,
        party_count,
        exec_ctx.as_ref().and_then(|c| c.syqure_port_base),
    )?;

    let mut env_map = build_shell_env(
        &spec.env,
        &BTreeMap::new(),
        &ctx,
        module_path,
        &results_path,
        "syqure",
    );

    // MPC file dir must match flow.yaml mpc.url pattern: shared/flows/{flow_name}/{run_id}/_mpc
    let flow_name = exec_ctx
        .as_ref()
        .and_then(|c| c.flow_name.clone())
        .or_else(|| env::var("BV_FLOW_NAME").ok())
        .unwrap_or_else(|| "syqure-flow".to_string());
    let file_dir = format!("shared/flows/{}/{}/_mpc", flow_name, run_id);

    env_map.insert("SEQURE_TRANSPORT".to_string(), transport.clone());
    env_map.insert("SEQURE_FILE_DIR".to_string(), file_dir);
    env_map.insert("SEQURE_FILE_POLL_MS".to_string(), poll_ms.to_string());
    env_map.insert("SEQURE_CP_COUNT".to_string(), party_count.to_string());
    env_map.insert("SEQURE_PARTY_EMAILS".to_string(), party_emails.clone());
    env_map.insert("SEQURE_DATASITES_ROOT".to_string(), datasites_root.clone());
    env_map.insert("SEQURE_LOCAL_EMAIL".to_string(), current_email.clone());
    env_map.insert("SEQURE_FILE_KEEP".to_string(), "1".to_string());
    env_map.insert("SEQURE_FILE_DEBUG".to_string(), "1".to_string());
    let mut docker_network_name: Option<String> = None;
    let mut docker_network_subnet: Option<String> = None;
    let mut docker_party_ip: Option<String> = None;

    // Docker+hotlink defaults to direct container TCP over a shared Docker network.
    // Set BV_SYQURE_DOCKER_DIRECT=0 to force legacy proxy mode for debugging.
    let docker_direct_tcp = use_docker
        && transport == "hotlink"
        && env::var("BV_SYQURE_DOCKER_DIRECT")
            .map(|v| is_truthy(&v))
            .unwrap_or(true);

    if docker_direct_tcp {
        if party_count > 253 {
            return Err(anyhow::anyhow!(
                "Syqure Docker direct mode supports up to 253 parties, got {}",
                party_count
            )
            .into());
        }
        let cp_ips: Vec<String> = if let Ok(raw) = env::var("BV_SYQURE_DOCKER_CP_IPS") {
            let values: Vec<String> = raw
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            if values.len() != party_count {
                return Err(anyhow::anyhow!(
                    "BV_SYQURE_DOCKER_CP_IPS must have {} entries (got {})",
                    party_count,
                    values.len()
                )
                .into());
            }
            values
        } else {
            (0..party_count)
                .map(|idx| format!("172.29.0.{}", idx + 2))
                .collect()
        };
        if party_id >= cp_ips.len() {
            return Err(anyhow::anyhow!(
                "Party id {} out of bounds for CP IP list of size {}",
                party_id,
                cp_ips.len()
            )
            .into());
        }
        let network_name =
            env::var("BV_SYQURE_DOCKER_NETWORK").unwrap_or_else(|_| format!("syqure-net-{}", run_id));
        let network_subnet =
            env::var("BV_SYQURE_DOCKER_SUBNET").unwrap_or_else(|_| "172.29.0.0/24".to_string());
        let base_port =
            SEQURE_COMMUNICATION_PORT + (party_id as usize) * SEQURE_COMMUNICATION_PORT_STRIDE;

        env_map.insert("SEQURE_TRANSPORT".to_string(), "tcp".to_string());
        env_map.insert("SEQURE_TCP_PROXY".to_string(), "0".to_string());
        env_map.insert("SEQURE_CP_IPS".to_string(), cp_ips.join(","));
        env_map.insert(
            "SEQURE_COMMUNICATION_PORT".to_string(),
            base_port.to_string(),
        );
        env_map.insert(
            "SEQURE_DATA_SHARING_PORT".to_string(),
            (base_port + SEQURE_DATA_SHARING_PORT_OFFSET).to_string(),
        );
        docker_party_ip = Some(cp_ips[party_id].clone());
        docker_network_name = Some(network_name);
        docker_network_subnet = Some(network_subnet);
        println!(
            "  Syqure Docker direct TCP enabled: network={} ip={}",
            docker_network_name.as_deref().unwrap_or(""),
            docker_party_ip.as_deref().unwrap_or("")
        );
    } else {
        // Legacy behavior: hotlink over host TCP proxy.
        let tcp_proxy = env_map
            .get("SEQURE_TCP_PROXY")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(transport == "hotlink");
        env_map.insert(
            "SEQURE_TCP_PROXY".to_string(),
            if tcp_proxy {
                "1".to_string()
            } else {
                "0".to_string()
            },
        );
        if tcp_proxy && transport == "hotlink" {
            env_map.insert("SEQURE_TRANSPORT".to_string(), "tcp".to_string());
        }
        if tcp_proxy {
            let proxy_ip = if use_docker {
                let container_runtime =
                    env::var("BIOVAULT_CONTAINER_RUNTIME").unwrap_or_else(|_| "docker".to_string());
                syqure_container_proxy_host(&container_runtime)
            } else {
                "127.0.0.1".to_string()
            };
            let proxy_ips = std::iter::repeat_n(proxy_ip.as_str(), party_count)
                .collect::<Vec<_>>()
                .join(",");
            env_map.insert("SEQURE_CP_IPS".to_string(), proxy_ips);
            let base_port =
                SEQURE_COMMUNICATION_PORT + (party_id as usize) * SEQURE_COMMUNICATION_PORT_STRIDE;
            env_map.insert(
                "SEQURE_COMMUNICATION_PORT".to_string(),
                base_port.to_string(),
            );
            env_map.insert(
                "SEQURE_DATA_SHARING_PORT".to_string(),
                (base_port + 10_000).to_string(),
            );
        }
    }
    for pass_through in &["SYQURE_SKIP_BUNDLE", "SYQURE_DEBUG"] {
        if let Ok(val) = env::var(pass_through) {
            env_map.insert(pass_through.to_string(), val);
        }
    }
    if env::var("SYQURE_BUNDLE_CACHE").is_err() {
        let cache_dir = PathBuf::from(&datasites_root)
            .join(".syqure-cache")
            .join(&run_id)
            .join(format!("party-{}", party_id));
        env_map.insert(
            "SYQURE_BUNDLE_CACHE".to_string(),
            cache_dir.to_string_lossy().to_string(),
        );
    }
    let codon_path_missing = env_map
        .get("CODON_PATH")
        .map(|v| v.trim().is_empty())
        .unwrap_or(true);
    if codon_path_missing {
        if let Some(codon_path) = resolve_codon_path_from_syqure_binary(&syqure_binary) {
            env_map.insert(
                "CODON_PATH".to_string(),
                codon_path.to_string_lossy().to_string(),
            );
        }
    }

    // Validate CODON_PATH/lib/codon and Sequre plugin — fail fast with actionable error.
    if let Some(codon) = env_map.get("CODON_PATH") {
        let codon_lib = PathBuf::from(codon).join("lib").join("codon");
        if !codon_lib.exists() {
            return Err(anyhow::anyhow!(
                "ERROR: CODON_PATH is set to '{}' but lib/codon subdirectory does not exist.\n\
                 Expected: {}\n\
                 The syqure binary will fail with exit code 127 on Linux.\n\
                 Check your syqure build or set CODON_PATH to a valid Codon installation.",
                codon,
                codon_lib.display()
            )
            .into());
        }

        #[cfg(target_os = "linux")]
        if codon.contains("/bin/macos-") {
            return Err(anyhow::anyhow!(
                "ERROR: CODON_PATH appears to be macOS on Linux: '{}'.\n\
                 Set CODON_PATH to a Linux codon bundle (e.g. .../bin/linux-x86/codon).",
                codon
            )
            .into());
        }

        #[cfg(target_os = "macos")]
        if codon.contains("/bin/linux-") {
            return Err(anyhow::anyhow!(
                "ERROR: CODON_PATH appears to be Linux on macOS: '{}'.\n\
                 Set CODON_PATH to a macOS codon bundle (e.g. .../bin/macos-arm64/codon).",
                codon
            )
            .into());
        }

        #[cfg(target_os = "linux")]
        let sequre_plugin = codon_lib
            .join("plugins")
            .join("sequre")
            .join("build")
            .join("libsequre.so");
        #[cfg(target_os = "macos")]
        let sequre_plugin = codon_lib
            .join("plugins")
            .join("sequre")
            .join("build")
            .join("libsequre.dylib");
        #[cfg(not(any(target_os = "linux", target_os = "macos")))]
        let sequre_plugin = codon_lib.join("plugins").join("sequre").join("build");

        if !sequre_plugin.exists() {
            return Err(anyhow::anyhow!(
                "ERROR: Sequre plugin library not found for CODON_PATH='{}'.\n\
                 Expected plugin path: {}\n\
                 This usually means CODON_PATH points at the wrong platform bundle.",
                codon,
                sequre_plugin.display()
            )
            .into());
        }

        println!("  Codon lib/codon verified: {}", codon_lib.display());
        println!("  Sequre plugin verified: {}", sequre_plugin.display());
    }

    println!(
        "  Syqure env: SEQURE_TRANSPORT={} SEQURE_TCP_PROXY={} SEQURE_CP_IPS={} SEQURE_COMMUNICATION_PORT={} SYQURE_SKIP_BUNDLE={} CODON_PATH={} SYQURE_DEBUG={}",
        env_map
            .get("SEQURE_TRANSPORT")
            .map(|s| s.as_str())
            .unwrap_or(""),
        env_map
            .get("SEQURE_TCP_PROXY")
            .map(|s| s.as_str())
            .unwrap_or(""),
        env_map.get("SEQURE_CP_IPS").map(|s| s.as_str()).unwrap_or(""),
        env_map
            .get("BV_SYQURE_PORT_BASE")
            .map(|s| s.as_str())
            .unwrap_or(""),
        env_map
            .get("SEQURE_COMMUNICATION_PORT")
            .map(|s| s.as_str())
            .unwrap_or(""),
        env_map
            .get("SYQURE_SKIP_BUNDLE")
            .map(|s| s.as_str())
            .unwrap_or(""),
        env_map.get("CODON_PATH").map(|s| s.as_str()).unwrap_or(""),
        env_map
            .get("SYQURE_DEBUG")
            .map(|s| s.as_str())
            .unwrap_or(""),
        &syqure_binary
    );
    println!(
        "  Syqure transport mode: {}",
        describe_syqure_transport_mode(&env_map)
    );
    println!(
        "  Hotlink flags: BV_SYFTBOX_HOTLINK={} BV_SYFTBOX_HOTLINK_QUIC(legacy p2p)={} BV_SYFTBOX_HOTLINK_QUIC_ONLY(legacy p2p-only)={}",
        env_value("BV_SYFTBOX_HOTLINK"),
        env_value("BV_SYFTBOX_HOTLINK_QUIC"),
        env_value("BV_SYFTBOX_HOTLINK_QUIC_ONLY")
    );
    println!(
        "  Syqure coordination: run_id={} party_id={}/{} email={} BV_SYFTBOX_BACKEND={} SYFTBOX_HOTLINK_TCP_PROXY={}",
        run_id,
        party_id,
        party_count,
        current_email,
        env_value("BV_SYFTBOX_BACKEND"),
        env_value("SYFTBOX_HOTLINK_TCP_PROXY"),
    );
    let hint_path = syqure_port_base_hint_path(&run_id, party_count);
    println!(
        "  Port hint file: {} (exists={})",
        hint_path.display(),
        hint_path.exists()
    );
    if tcp_proxy {
        let effective_comm_port = env_map
            .get("SEQURE_COMMUNICATION_PORT")
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(syqure_port_base);
        let expected_comm_port = syqure_port_base + party_id * SEQURE_COMMUNICATION_PORT_STRIDE;
        let mismatch = effective_comm_port != expected_comm_port;
        println!(
            "  TCP proxy ports: party_id={} party_count={} global_base={} expected_comm_port={} effective_comm_port={} stride={}{}",
            party_id,
            party_count,
            syqure_port_base,
            expected_comm_port,
            effective_comm_port,
            SEQURE_COMMUNICATION_PORT_STRIDE,
            if mismatch { " [MISMATCH]" } else { "" }
        );

        let mut peer_ports = Vec::new();
        for peer_id in 0..party_count {
            if peer_id == party_id {
                continue;
            }
            let port = mpc_comm_port_with_base(effective_comm_port, party_id, peer_id, party_count);
            peer_ports.push(format!("CP{}<->CP{}={}", party_id, peer_id, port));
        }
        if !peer_ports.is_empty() {
            println!(
                "  TCP proxy expected peer ports (from CP{} view): {}",
                party_id,
                peer_ports.join(", ")
            );
        }
    }

    for input in &spec.inputs {
        let env_key = format!("BV_INPUT_{}", env_key_suffix(&input.name));
        let sequre_env_key = format!("SEQURE_INPUT_{}", env_key_suffix(&input.name));
        if let Some(arg) = parsed_args.inputs.get(&input.name) {
            let rendered_value = render_template(&arg.value, &ctx);
            let final_value = if input.raw_type == "File" {
                let input_path = if Path::new(&rendered_value).is_absolute() {
                    PathBuf::from(rendered_value)
                } else {
                    module_path.join(rendered_value)
                };
                input_path.to_string_lossy().to_string()
            } else {
                rendered_value
            };
            env_map.insert(env_key.clone(), final_value.clone());
            env_map.insert(sequre_env_key, final_value);
        } else if let Some(path_template) = input.path.as_deref() {
            let rendered_path = render_template(path_template, &ctx);
            let input_path = if Path::new(&rendered_path).is_absolute() {
                PathBuf::from(rendered_path)
            } else {
                module_path.join(rendered_path)
            };
            let path_str = input_path.to_string_lossy().to_string();
            env_map.insert(env_key, path_str.clone());
            env_map.insert(sequre_env_key, path_str);
        }
    }
    if let Some(array_length_path) = env_map.get("SEQURE_INPUT_ARRAY_LENGTH").cloned() {
        if let Ok(raw) = fs::read_to_string(&array_length_path) {
            let trimmed = raw.trim();
            if !trimmed.is_empty() {
                env_map.insert("SEQURE_ARRAY_LENGTH".to_string(), trimmed.to_string());
            }
        }
    }

    for output in &spec.outputs {
        let env_key = format!("BV_OUTPUT_{}", env_key_suffix(&output.name));
        let sequre_env_key = format!("SEQURE_OUTPUT_{}", env_key_suffix(&output.name));
        if let Some(path_template) = output.path.as_deref() {
            let rendered_path = render_template(path_template, &ctx);
            let output_path = results_path.join(&rendered_path);
            let path_str = output_path.to_string_lossy().to_string();
            env_map.insert(env_key, path_str.clone());
            env_map.insert(sequre_env_key, path_str);
        }
    }

    if !dry_run {
        let skip_done = env::var("BV_FLOW_SKIP_DONE").unwrap_or_default() == "1";
        let force_run = env::var("BV_FLOW_FORCE").unwrap_or_default() == "1";
        if skip_done && !force_run && !spec.outputs.is_empty() {
            let mut all_outputs_present = true;
            for output in &spec.outputs {
                let env_key = format!("BV_OUTPUT_{}", env_key_suffix(&output.name));
                if let Some(path) = env_map.get(&env_key) {
                    if fs::metadata(path).is_err() {
                        all_outputs_present = false;
                        break;
                    }
                } else {
                    all_outputs_present = false;
                    break;
                }
            }
            if all_outputs_present {
                println!("⏭️  Skipping syqure step (outputs present; BV_FLOW_SKIP_DONE=1)");
                return Ok(());
            }
        }
    }

    if dry_run {
        println!("\n📋 [DRY RUN] Syqure execution plan:");
        println!(
            "  Backend: {}",
            if use_docker { "docker" } else { "native" }
        );
        println!("  Binary/Image: {}", syqure_binary);
        println!("  Entrypoint: {}", entrypoint_path.display());
        println!("  Party ID: {}", party_id);
        println!("  Party count: {}", party_count);
        println!("  Run ID: {}", run_id);
        println!("\n  Environment:");
        for (k, v) in env_map
            .iter()
            .filter(|(k, _)| k.starts_with("SEQURE_") || k.starts_with("BV_"))
        {
            println!("    {}={}", k, v);
        }
        return Ok(());
    }

    fs::create_dir_all(&datasites_root)?;

    // Write unified diagnostic file (same format for both CLI and Tauri)
    write_syqure_diag_file(
        &run_id,
        party_id,
        party_count,
        &current_email,
        syqure_port_base,
        tcp_proxy,
        &env_map,
        module_path,
        &syqure_binary,
        &entrypoint_path,
    );

    if tcp_proxy {
        println!("  Pre-launch TCP proxy port check:");
        for peer_id in 0..party_count {
            if peer_id == party_id {
                continue;
            }
            let comm_base = syqure_port_base + party_id * SEQURE_COMMUNICATION_PORT_STRIDE;
            let port = mpc_comm_port_with_base(comm_base, party_id, peer_id, party_count);
            let listening = TcpListener::bind(SocketAddrV4::new(Ipv4Addr::LOCALHOST, port as u16))
                .map(|_| false)
                .unwrap_or(true);
            println!(
                "    CP{}<->CP{} port {} {}",
                party_id,
                peer_id,
                port,
                if listening {
                    "LISTENING (proxy ready)"
                } else {
                    "FREE (proxy NOT yet bound)"
                }
            );
        }
    }

    if use_docker {
        execute_syqure_docker(
            &syqure_binary,
            &entrypoint_path,
            party_id,
            &env_map,
            module_path,
            &datasites_root,
            &run_id,
            &docker_platform,
            docker_network_name.as_deref(),
            docker_network_subnet.as_deref(),
            docker_party_ip.as_deref(),
        )?;
    } else {
        println!(
            "  [trace] spawn_blocking execute_syqure_native party_id={} thread={:?}",
            party_id,
            std::thread::current().id()
        );
        let binary_owned = syqure_binary.to_string();
        let entrypoint_owned = entrypoint_path.to_path_buf();
        let env_map_owned = env_map.clone();
        let sb_result = tokio::task::spawn_blocking(move || {
            execute_syqure_native(&binary_owned, &entrypoint_owned, party_id, &env_map_owned)
        })
        .await;
        println!("  [trace] spawn_blocking returned party_id={}", party_id);
        sb_result.map_err(|e| anyhow::anyhow!("Syqure task join error: {}", e))??;
    }

    println!("\n✅ Syqure execution completed successfully!");
    Ok(())
}

fn resolve_syqure_backend(spec: &ModuleSpec) -> Result<(String, bool)> {
    let syqure_opts = spec.runner.as_ref().and_then(|r| r.syqure.as_ref());

    let force_docker_env = env::var("BV_SYQURE_USE_DOCKER")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false);

    let force_docker_spec = syqure_opts.and_then(|s| s.use_docker).unwrap_or(false);

    let is_windows = cfg!(target_os = "windows");
    let use_docker = force_docker_env || force_docker_spec || is_windows;

    if use_docker {
        let docker_image = syqure_opts
            .and_then(|s| s.docker_image.as_ref())
            .cloned()
            .or_else(|| {
                crate::config::Config::load()
                    .ok()
                    .and_then(|c| c.syqure.and_then(|s| s.docker_image))
            })
            .unwrap_or_else(|| "ghcr.io/madhavajay/syqure-cli:latest".to_string());
        return Ok((docker_image, true));
    }

    if let Some(binary) = syqure_opts.and_then(|s| s.binary.as_ref()) {
        if Path::new(binary).exists() {
            return Ok((binary.clone(), false));
        }
    }

    if let Ok(bin) = env::var("SEQURE_NATIVE_BIN") {
        if Path::new(&bin).exists() {
            return Ok((bin, false));
        }
    }

    if let Ok(config) = crate::config::Config::load() {
        if let Some(syqure_cfg) = config.syqure {
            if let Some(use_docker_cfg) = syqure_cfg.use_docker {
                if use_docker_cfg {
                    let image = syqure_cfg
                        .docker_image
                        .unwrap_or_else(|| "ghcr.io/madhavajay/syqure-cli:latest".to_string());
                    return Ok((image, true));
                }
            }
            if let Some(binary) = syqure_cfg.binary {
                if Path::new(&binary).exists() {
                    return Ok((binary, false));
                }
            }
        }
        if let Some(paths) = config.binary_paths {
            if let Some(syqure_path) = paths.syqure {
                if Path::new(&syqure_path).exists() {
                    return Ok((syqure_path, false));
                }
            }
        }
    }

    let mut which_cmd = Command::new("which");
    super::configure_child_process(&mut which_cmd);
    if let Ok(which_output) = which_cmd.arg("syqure").output() {
        if which_output.status.success() {
            let path = String::from_utf8_lossy(&which_output.stdout)
                .trim()
                .to_string();
            if !path.is_empty() && Path::new(&path).exists() {
                return Ok((path, false));
            }
        }
    }

    Err(anyhow::anyhow!(
        "Syqure binary not found. Set SEQURE_NATIVE_BIN environment variable, \
         configure binary path in biovault config, or use Docker mode (use_docker: true)"
    )
    .into())
}

fn resolve_codon_path_from_syqure_binary(syqure_binary: &str) -> Option<PathBuf> {
    let bin_path = PathBuf::from(syqure_binary);
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    let platforms = ["linux-x86", "linux-x86_64", "linux-amd64", "macos-arm64"];
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    let platforms = ["macos-arm64", "linux-x86", "linux-x86_64", "linux-amd64"];
    #[cfg(not(any(
        all(target_os = "linux", target_arch = "x86_64"),
        all(target_os = "macos", target_arch = "aarch64")
    )))]
    let platforms = ["linux-x86", "linux-x86_64", "linux-amd64", "macos-arm64"];


    for ancestor in bin_path.ancestors() {
        let mut candidates = Vec::new();
        for platform in platforms {
            candidates.push(ancestor.join("bin").join(platform).join("codon"));
        }
        candidates.push(ancestor.join("bin/codon"));
        candidates.push(ancestor.join("codon/install"));
        candidates.push(ancestor.join("codon"));
        for candidate in candidates {
            if candidate.join("lib/codon/stdlib").exists() {
                return Some(candidate);
            }
        }
    }
    None
}

fn execute_syqure_native(
    binary: &str,
    entrypoint: &Path,
    party_id: usize,
    env_map: &BTreeMap<String, String>,
) -> Result<()> {
    println!(
        "  [trace] execute_syqure_native START party_id={} thread={:?} pid={}",
        party_id,
        std::thread::current().id(),
        std::process::id()
    );

    // ── Pre-flight: validate binary ──────────────────────────────────
    // Fail loudly before spawning so CI logs immediately show what's wrong.
    let binary_path = Path::new(binary);
    if !binary_path.exists() {
        return Err(anyhow::anyhow!(
            "FATAL: syqure binary does not exist: {}\n\
             Set SEQURE_NATIVE_BIN to a valid path or build syqure first.",
            binary
        )
        .into());
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = fs::metadata(binary_path)
            .map(|m| m.permissions().mode())
            .unwrap_or(0);
        if mode & 0o111 == 0 {
            return Err(anyhow::anyhow!(
                "FATAL: syqure binary is not executable (mode {:o}): {}\n\
                 Run: chmod +x {}",
                mode,
                binary,
                binary
            )
            .into());
        }
    }
    println!("  Using native syqure: {}", binary);

    // ── Pre-flight: validate CODON_PATH and Codon shared libraries ───
    // The syqure binary dynamically links libcodonrt and libcodonc.
    // On macOS, @loader_path rpaths handle this automatically.
    // On Linux, the dynamic linker needs LD_LIBRARY_PATH to point at
    // the Codon lib directory, otherwise the binary exits with code 127
    // ("command not found" / shared library load failure).
    let codon_path = env_map.get("CODON_PATH").cloned().unwrap_or_default();
    let codon_lib_dir = if !codon_path.is_empty() {
        let lib_dir = PathBuf::from(&codon_path).join("lib").join("codon");
        if !lib_dir.exists() {
            return Err(anyhow::anyhow!(
                "FATAL: CODON_PATH/lib/codon does not exist: {}\n\
                 Refusing to start syqure with invalid CODON_PATH='{}'.",
                lib_dir.display(),
                codon_path
            )
            .into());
        } else {
            println!("  Codon lib dir: {}", lib_dir.display());
            Some(lib_dir)
        }
    } else {
        // Try to find codon libs relative to the syqure binary as a fallback.
        // Dev builds: syqure/target/{release,debug}/syqure → syqure/bin/<platform>/codon/lib/codon
        // Bundled:    resources/syqure/syqure → resources/syqure/lib/codon
        let mut found = None;
        if let Some(bin_parent) = binary_path.parent() {
            let candidates = [
                bin_parent.join("lib").join("codon"),
                bin_parent.join("..").join("lib").join("codon"),
            ];
            for c in &candidates {
                if c.exists() {
                    found = Some(c.clone());
                    break;
                }
            }
        }
        if let Some(ref lib_dir) = found {
            println!(
                "  Codon lib dir (auto-detected near binary): {}",
                lib_dir.display()
            );
        } else {
            return Err(anyhow::anyhow!(
                "FATAL: CODON_PATH not set and no codon libs found near binary.\n\
                 Refusing to start syqure.\n\
                 Binary: {}",
                binary
            )
            .into());
        }
        found
    };

    if let Some(ref lib_dir) = codon_lib_dir {
        #[cfg(target_os = "linux")]
        let sequre_plugin = lib_dir
            .join("plugins")
            .join("sequre")
            .join("build")
            .join("libsequre.so");
        #[cfg(target_os = "macos")]
        let sequre_plugin = lib_dir
            .join("plugins")
            .join("sequre")
            .join("build")
            .join("libsequre.dylib");
        #[cfg(not(any(target_os = "linux", target_os = "macos")))]
        let sequre_plugin = lib_dir.join("plugins").join("sequre").join("build");

        if !sequre_plugin.exists() {
            return Err(anyhow::anyhow!(
                "FATAL: Sequre plugin library not found before syqure start.\n\
                 Expected: {}\n\
                 CODON_PATH='{}'",
                sequre_plugin.display(),
                codon_path
            )
            .into());
        }
        println!("  Sequre plugin verified: {}", sequre_plugin.display());
    }

    let mut cmd = Command::new(binary);
    super::configure_child_process(&mut cmd);
    cmd.arg(entrypoint);
    cmd.arg("--");
    cmd.arg(party_id.to_string());

    if !env_map.contains_key("CODON_PATH") {
        cmd.env_remove("CODON_PATH");
    }
    if !env_map.contains_key("CODON_PLUGIN_PATH") {
        cmd.env_remove("CODON_PLUGIN_PATH");
    }

    for (k, v) in env_map {
        cmd.env(k, v);
    }

    // ── LD_LIBRARY_PATH: ensure Codon shared libs are discoverable ───
    // On Linux, the syqure binary links libcodonrt.so and libcodonc.so
    // dynamically. The compiled rpath ($ORIGIN/lib/codon) only works when
    // libs are co-located with the binary. In dev/CI the libs live under
    // CODON_PATH/lib/codon or the bundle cache, so we must add them to
    // LD_LIBRARY_PATH explicitly. On macOS this is harmless (dyld ignores it).
    if let Some(ref lib_dir) = codon_lib_dir {
        let lib_dir_str = lib_dir.to_string_lossy().to_string();
        // Also include the parent (lib/) for any top-level .so files.
        let lib_parent = lib_dir.parent().map(|p| p.to_string_lossy().to_string());
        let existing_ld = env::var("LD_LIBRARY_PATH").unwrap_or_default();
        let mut parts: Vec<String> = vec![lib_dir_str.clone()];
        if let Some(ref parent) = lib_parent {
            parts.push(parent.clone());
        }
        if !existing_ld.is_empty() {
            parts.push(existing_ld);
        }
        let new_ld = parts.join(":");
        cmd.env("LD_LIBRARY_PATH", &new_ld);
        println!("  LD_LIBRARY_PATH={}", new_ld);
    }
    if env::var("BV_SYQURE_BACKTRACE")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
    {
        cmd.env("RUST_BACKTRACE", "1");
    }

    let display_cmd = format_command(&cmd);
    println!("\n▶️  Executing syqure...");
    println!("  {}", display_cmd.dimmed());

    // Write syqure log to private cache dir (not shared results dir) to avoid
    // leaking keys through shared step outputs.
    let syqure_log_path = env_map
        .get("SYQURE_BUNDLE_CACHE")
        .or_else(|| env_map.get("BV_RESULTS_DIR"))
        .map(PathBuf::from)
        .map(|dir| dir.join("syqure-native.log"));

    // Syqure startup currently removes stale `sock.*` files from the process CWD.
    // In multiparty desktop runs, parties execute concurrently and can otherwise
    // delete each other's active sockets if they share one working directory.
    let syqure_work_dir = env_map
        .get("SYQURE_BUNDLE_CACHE")
        .map(PathBuf::from)
        .or_else(|| env_map.get("BV_RESULTS_DIR").map(PathBuf::from))
        .or_else(|| env::current_dir().ok())
        .unwrap_or_else(|| PathBuf::from("."));
    if let Err(err) = fs::create_dir_all(&syqure_work_dir) {
        return Err(anyhow::anyhow!(
            "Failed to prepare Syqure working directory {}: {}",
            syqure_work_dir.display(),
            err
        )
        .into());
    }
    cmd.current_dir(&syqure_work_dir);
    if let Some(path) = syqure_log_path.as_ref() {
        if let Some(parent) = path.parent() {
            let _ = fs::create_dir_all(parent);
        }
        if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(path) {
            let _ = writeln!(
                f,
                "[{}] starting syqure command: {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S%:z"),
                display_cmd
            );
            for (k, v) in env_map.iter().filter(|(k, _)| {
                k.starts_with("SEQURE_")
                    || k.starts_with("BV_")
                    || k.starts_with("SYFTBOX_")
                    || k.starts_with("SYQURE_")
                    || k.starts_with("CODON_")
                    || k.starts_with("BIOVAULT_")
            }) {
                let _ = writeln!(f, "[env] {}={}", k, v);
            }
            let _ = writeln!(f, "[env] SYQURE_WORK_DIR={}", syqure_work_dir.display());
        }
        println!("  Syqure output log: {}", path.display());
    }

    cmd.stdout(Stdio::piped());
    cmd.stderr(Stdio::piped());

    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        // SAFETY: pre_exec is unsafe because it runs in the child after fork.
        // We only call async-signal-safe libc here.
        unsafe {
            cmd.pre_exec(|| {
                // Start a new process group so we can terminate syqure and its children together.
                let rc = libc::setpgid(0, 0);
                if rc != 0 {
                    return Err(std::io::Error::last_os_error());
                }
                Ok(())
            });
        }
    }

    let datasites_root = env_map
        .get("SEQURE_DATASITES_ROOT")
        .cloned()
        .unwrap_or_default();
    let file_dir = env_map.get("SEQURE_FILE_DIR").cloned().unwrap_or_default();
    let party_emails: Vec<String> = env_map
        .get("SEQURE_PARTY_EMAILS")
        .map(|s| s.split(',').map(|e| e.to_string()).collect())
        .unwrap_or_default();
    let current_email = env_map
        .get("SEQURE_LOCAL_EMAIL")
        .cloned()
        .unwrap_or_default();
    let tcp_transport = env_map
        .get("SEQURE_TRANSPORT")
        .map(|v| v == "tcp")
        .unwrap_or(false);
    let hotlink_transport = is_hotlink_transport_mode(env_map);
    let status_label = if hotlink_transport {
        "  📡 Hotlink (↑tx ↓rx data p2p-status)"
    } else if tcp_transport {
        "  📡 MPC dir activity (files/age)"
    } else {
        "  📡 Syqure channels"
    };

    let stop_flag = Arc::new(AtomicBool::new(false));
    let stop_flag_clone = stop_flag.clone();

    let datasites_root_clone = datasites_root.clone();
    let file_dir_clone = file_dir.clone();
    let party_emails_clone = party_emails.clone();
    let status_label_clone = status_label.to_string();
    let hotlink_transport_clone = hotlink_transport;

    if hotlink_transport {
        println!(
            "  ℹ️  Hotlink stats are per-party packet counters from .syftbox/hotlink_telemetry.json."
        );
    } else if tcp_transport {
        println!(
            "  ℹ️  TCP/hotlink mode: channel counters below are directory activity, not packet latency."
        );
    }

    let monitor_handle = thread::spawn(move || {
        let mut last_status = String::new();
        while !stop_flag_clone.load(Ordering::Relaxed) {
            let now = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs();

            let mut status_parts = Vec::new();
            for party in &party_emails_clone {
                let short_party = party.split('@').next().unwrap_or(party);
                if hotlink_transport_clone {
                    let telemetry_paths = [
                        PathBuf::from(&datasites_root_clone)
                            .join(party)
                            .join(".syftbox")
                            .join("hotlink_telemetry.json"),
                        PathBuf::from(&datasites_root_clone)
                            .join(party)
                            .join("datasites")
                            .join(party)
                            .join(".syftbox")
                            .join("hotlink_telemetry.json"),
                    ];
                    let mut telemetry_found = false;
                    for telemetry_path in telemetry_paths {
                        if let Some(t) = read_hotlink_telemetry(&telemetry_path) {
                            let p2p_status = if t.webrtc_connected > 0 {
                                format!("p2p:UP({})", t.webrtc_connected)
                            } else if t.tx_p2p_packets > 0 {
                                "p2p:used".to_string()
                            } else {
                                "p2p:off".to_string()
                            };
                            let total = t.tx_packets + t.rx_packets;
                            if total == 0 {
                                status_parts.push(format!(
                                    "{}: idle {}",
                                    short_party, p2p_status
                                ));
                            } else {
                                status_parts.push(format!(
                                    "{}: ↑{} ↓{} {} {}",
                                    short_party,
                                    t.tx_packets,
                                    t.rx_packets,
                                    fmt_hotlink_bytes(t.tx_bytes + t.rx_bytes),
                                    p2p_status,
                                ));
                            }
                            telemetry_found = true;
                            break;
                        }
                    }
                    if telemetry_found {
                        continue;
                    }
                    status_parts.push(format!("{}: waiting", short_party));
                    continue;
                }
                let channel_path = PathBuf::from(&datasites_root_clone)
                    .join(party)
                    .join(&file_dir_clone);
                let (count, last_mod) = count_mpc_files(&channel_path);
                let ago = last_mod.map(|lm| now.saturating_sub(lm)).unwrap_or(999);
                status_parts.push(format!("{}:{}f/{}s", short_party, count, ago));
            }

            let status = status_parts.join(" | ");
            if status != last_status {
                println!("{}: {}", status_label_clone, status);
                last_status = status;
            }

            thread::sleep(Duration::from_secs(5));
        }
    });

    let start_time = std::time::Instant::now();
    let syqure_timeout_s = env_map
        .get("BV_SYQURE_NATIVE_TIMEOUT_S")
        .and_then(|v| v.trim().parse::<u64>().ok())
        .or_else(|| {
            env::var("BV_SYQURE_NATIVE_TIMEOUT_S")
                .ok()
                .and_then(|v| v.trim().parse::<u64>().ok())
        })
        .filter(|v| *v > 0)
        .unwrap_or(600);
    let syqure_timeout = Duration::from_secs(syqure_timeout_s);
    println!(
        "  [trace] syqure timeout={}s party_id={}",
        syqure_timeout_s, party_id
    );
    println!("  [trace] spawning syqure child process...");
    let mut child = cmd.spawn().context("Failed to spawn syqure")?;
    println!(
        "  [trace] syqure child spawned pid={} party_id={}",
        child.id(),
        party_id
    );
    let mut io_threads = Vec::new();
    if let Some(stdout) = child.stdout.take() {
        io_threads.push(spawn_child_stream_capture(
            stdout,
            "stdout",
            syqure_log_path.clone(),
            false,
        ));
    }
    if let Some(stderr) = child.stderr.take() {
        io_threads.push(spawn_child_stream_capture(
            stderr,
            "stderr",
            syqure_log_path.clone(),
            true,
        ));
    }
    let flow_name = env_map
        .get("BV_FLOW_NAME")
        .cloned()
        .or_else(|| env::var("BV_FLOW_NAME").ok())
        .unwrap_or_default();
    let run_id = env_map
        .get("BV_RUN_ID")
        .cloned()
        .or_else(|| env::var("BV_RUN_ID").ok())
        .unwrap_or_default();
    let step_id = env::var("BV_MULTIPARTY_STEP_ID").unwrap_or_default();
    let mut peer_failure_reason: Option<String> = None;

    #[cfg(unix)]
    fn kill_process_group(pid: u32) {
        let pgid = -(pid as i32);
        unsafe {
            libc::kill(pgid, libc::SIGTERM);
        }
        // Give it a moment to exit gracefully, then force kill.
        std::thread::sleep(Duration::from_millis(250));
        unsafe {
            libc::kill(pgid, libc::SIGKILL);
        }
    }

    #[cfg(unix)]
    struct SyqureKillGuard {
        pid: u32,
        child_ptr: *mut std::process::Child,
    }

    #[cfg(unix)]
    impl Drop for SyqureKillGuard {
        fn drop(&mut self) {
            unsafe {
                if let Ok(None) = (*self.child_ptr).try_wait() {
                    kill_process_group(self.pid);
                }
            }
        }
    }

    #[cfg(unix)]
    static SYQURE_SHUTDOWN_REQUESTED: AtomicBool = AtomicBool::new(false);

    #[cfg(unix)]
    unsafe extern "C" fn syqure_signal_handler(_: libc::c_int) {
        SYQURE_SHUTDOWN_REQUESTED.store(true, Ordering::SeqCst);
    }

    #[cfg(unix)]
    {
        // This runner can execute multiple Syqure steps in a long-lived desktop process.
        // Reset the global shutdown latch before each execution so stale state from a
        // previous run does not terminate the new process immediately.
        SYQURE_SHUTDOWN_REQUESTED.store(false, Ordering::SeqCst);
        println!(
            "  [trace] installing signal handlers for syqure pid={} party_id={}",
            std::process::id(),
            party_id
        );
        unsafe {
            libc::signal(libc::SIGINT, syqure_signal_handler as libc::sighandler_t);
            libc::signal(libc::SIGTERM, syqure_signal_handler as libc::sighandler_t);
        }
    }

    #[cfg(unix)]
    let _kill_guard = SyqureKillGuard {
        pid: child.id(),
        child_ptr: &mut child as *mut std::process::Child,
    };

    let mut poll_count: u64 = 0;
    let status = loop {
        match child.try_wait().context("Failed to poll syqure")? {
            Some(status) => {
                println!(
                    "  [trace] syqure child exited status={:?} party_id={} after {:.1}s polls={}",
                    status.code(),
                    party_id,
                    start_time.elapsed().as_secs_f64(),
                    poll_count
                );
                break status;
            }
            None => {
                poll_count += 1;
                if peer_failure_reason.is_none() && start_time.elapsed() >= syqure_timeout {
                    let reason = format!(
                        "Syqure timed out after {}s for party_id={}; terminating process group",
                        syqure_timeout_s, party_id
                    );
                    println!("  ⚠️  {}", reason);
                    #[cfg(unix)]
                    kill_process_group(child.id());
                    #[cfg(not(unix))]
                    {
                        let _ = child.kill();
                    }
                    peer_failure_reason = Some(reason);
                }
                // Log every 10s (50 polls at 200ms)
                if poll_count.is_multiple_of(50) {
                    println!(
                        "  [trace] syqure still running party_id={} pid={} elapsed={:.1}s polls={}",
                        party_id,
                        child.id(),
                        start_time.elapsed().as_secs_f64(),
                        poll_count
                    );
                }
                #[cfg(unix)]
                if SYQURE_SHUTDOWN_REQUESTED.load(Ordering::SeqCst) {
                    if peer_failure_reason.is_none() {
                        let reason = "Shutdown requested; stopping local Syqure".to_string();
                        println!("  ⚠️  {}", reason);
                        peer_failure_reason = Some(reason);
                    }
                    kill_process_group(child.id());
                }
                if peer_failure_reason.is_none()
                    && !datasites_root.is_empty()
                    && !flow_name.is_empty()
                    && !run_id.is_empty()
                    && !step_id.is_empty()
                    && !party_emails.is_empty()
                    && !current_email.is_empty()
                {
                    if let Some(peer) = detect_peer_step_failure(
                        &datasites_root,
                        &flow_name,
                        &run_id,
                        &step_id,
                        &current_email,
                        &party_emails,
                    ) {
                        let reason = format!(
                            "Peer failure detected for step '{}' from {}; stopping local Syqure",
                            step_id, peer
                        );
                        println!("  ⚠️  {}", reason);
                        #[cfg(unix)]
                        kill_process_group(child.id());
                        #[cfg(not(unix))]
                        {
                            let _ = child.kill();
                        }
                        peer_failure_reason = Some(reason);
                    }
                }
                thread::sleep(Duration::from_millis(200));
            }
        }
    };

    stop_flag.store(true, Ordering::Relaxed);
    let _ = monitor_handle.join();
    for handle in io_threads {
        let _ = handle.join();
    }

    let elapsed = start_time.elapsed();
    println!(
        "  ⏱️  Syqure step duration: {}.{:03}s",
        elapsed.as_secs(),
        elapsed.subsec_millis()
    );

    if let Some(reason) = peer_failure_reason {
        return Err(anyhow::anyhow!(reason).into());
    }

    if !status.success() {
        let syqure_log_tail = syqure_log_path
            .as_ref()
            .and_then(|path| fs::read_to_string(path).ok())
            .and_then(|contents| {
                let tail = contents
                    .lines()
                    .rev()
                    .take(80)
                    .collect::<Vec<_>>()
                    .into_iter()
                    .rev()
                    .collect::<Vec<_>>()
                    .join("\n");
                if tail.trim().is_empty() {
                    None
                } else {
                    Some(tail)
                }
            });

        #[cfg(unix)]
        {
            use std::os::unix::process::ExitStatusExt;
            if let Some(signal) = status.signal() {
                let mut msg = format!("Syqure exited due to signal: {}", signal);
                if let Some(tail) = syqure_log_tail.as_ref() {
                    msg.push_str("\n--- syqure-native.log tail ---\n");
                    msg.push_str(tail);
                }
                return Err(anyhow::anyhow!(msg).into());
            }
        }
        let mut msg = format!("Syqure exited with code: {:?}", status.code());
        if let Some(tail) = syqure_log_tail {
            msg.push_str("\n--- syqure-native.log tail ---\n");
            msg.push_str(&tail);
        }
        return Err(anyhow::anyhow!(msg).into());
    }

    println!(
        "  [trace] execute_syqure_native END party_id={} ok",
        party_id
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn execute_syqure_docker(
    image: &str,
    entrypoint: &Path,
    party_id: usize,
    env_map: &BTreeMap<String, String>,
    module_path: &Path,
    datasites_root: &str,
    run_id: &str,
    platform: &str,
    docker_network_name: Option<&str>,
    docker_network_subnet: Option<&str>,
    docker_party_ip: Option<&str>,
) -> Result<()> {
    fn to_container_path(path: &Path) -> String {
        path.to_string_lossy().replace('\\', "/")
    }

    fn map_env_path(
        value: &str,
        module_root: &Path,
        results_root: &Path,
        datasites_root: &Path,
        datasites_in_container: &str,
    ) -> Option<String> {
        let stripped = strip_extended_path_prefix(value);
        let value_path = PathBuf::from(&stripped);
        if let Ok(rel_path) = value_path.strip_prefix(results_root) {
            return Some(format!("/results/{}", to_container_path(rel_path)));
        }
        if let Ok(rel_path) = value_path.strip_prefix(module_root) {
            return Some(format!(
                "/workspace/project/{}",
                to_container_path(rel_path)
            ));
        }
        if let Ok(rel_path) = value_path.strip_prefix(datasites_root) {
            return Some(format!(
                "{}/{}",
                datasites_in_container.trim_end_matches('/'),
                to_container_path(rel_path)
            ));
        }
        None
    }

    // Use podman if BIOVAULT_CONTAINER_RUNTIME is set to "podman"
    let container_runtime =
        env::var("BIOVAULT_CONTAINER_RUNTIME").unwrap_or_else(|_| "docker".to_string());

    println!(
        "  Using {} image: {}",
        if container_runtime == "podman" {
            "Podman"
        } else {
            "Docker"
        },
        image
    );

    let container_name = format!("syqure-{}-pid{}", run_id, party_id);
    let keep_containers = env_var_truthy("BIOVAULT_SYQURE_KEEP_CONTAINERS")
        || env_var_truthy("BV_SYQURE_KEEP_CONTAINERS");

    let mut rm_cmd = Command::new(&container_runtime);
    super::configure_child_process(&mut rm_cmd);
    let _ = rm_cmd
        .args(["rm", "-f", &container_name])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();

    if let Some(network_name) = docker_network_name {
        if container_runtime != "docker" {
            return Err(anyhow::anyhow!(
                "Syqure Docker direct mode requires docker runtime, got {}",
                container_runtime
            )
            .into());
        }
        let mut inspect = Command::new(&container_runtime);
        super::configure_child_process(&mut inspect);
        let exists = inspect
            .args(["network", "inspect", network_name])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !exists {
            let mut create = Command::new(&container_runtime);
            super::configure_child_process(&mut create);
            create.args(["network", "create", "--driver", "bridge"]);
            if let Some(subnet) = docker_network_subnet {
                create.args(["--subnet", subnet]);
            }
            create.arg(network_name);
            let created = create.output().with_context(|| {
                format!("Failed to create docker network {}", network_name)
            })?;
            if !created.status.success() {
                let mut recheck = Command::new(&container_runtime);
                super::configure_child_process(&mut recheck);
                let already_exists = recheck
                    .args(["network", "inspect", network_name])
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .map(|s| s.success())
                    .unwrap_or(false);
                if !already_exists {
                    return Err(anyhow::anyhow!(
                        "Failed creating docker network {}: {}",
                        network_name,
                        String::from_utf8_lossy(&created.stderr)
                    )
                    .into());
                }
            }
        }
        println!(
            "  Docker direct network: {} ({})",
            network_name,
            docker_network_subnet.unwrap_or("subnet unspecified")
        );
    }

    let module_path_abs = module_path
        .canonicalize()
        .unwrap_or_else(|_| module_path.to_path_buf());

    let mut datasites_in_container = "/datasites".to_string();
    let results_in_container = "/results";

    let shared_datasites_root = env::var("BV_SHARED_DATASITES_ROOT").ok();
    let effective_datasites_mount = shared_datasites_root.as_deref().unwrap_or(datasites_root);

    let entrypoint_rel = entrypoint.strip_prefix(module_path).unwrap_or(entrypoint);
    let container_entrypoint = format!("/workspace/project/{}", to_container_path(entrypoint_rel));

    let results_root = env_map
        .get("BV_RESULTS_DIR")
        .and_then(|p| {
            let path = PathBuf::from(p);
            path.parent()
                .and_then(|p| p.parent())
                .map(|p| p.to_path_buf())
        })
        .unwrap_or_else(|| module_path.join("results"));

    let module_root_for_match = PathBuf::from(strip_extended_path_prefix(
        module_path_abs.to_string_lossy().as_ref(),
    ));
    let datasites_mount_path = strip_extended_path_prefix(effective_datasites_mount);
    let results_mount_path = strip_extended_path_prefix(results_root.to_string_lossy().as_ref());
    fs::create_dir_all(&datasites_mount_path).with_context(|| {
        format!(
            "Failed to create datasites mount source: {}",
            datasites_mount_path
        )
    })?;
    fs::create_dir_all(&results_mount_path).with_context(|| {
        format!(
            "Failed to create results mount source: {}",
            results_mount_path
        )
    })?;

    let results_root_for_match = PathBuf::from(&results_mount_path);
    let datasites_root_for_match = PathBuf::from(&datasites_mount_path);

    // Docker Desktop can fail mounting nested paths with special characters in parent directories
    // (e.g. ".../client1@sandbox.local/datasites"). Mount the parent and remap root in-container.
    let mut datasites_mount_source = datasites_mount_path.clone();
    let mut datasites_volume_target = datasites_in_container.clone();
    if container_runtime == "docker" && datasites_mount_path.ends_with("/datasites") {
        if let Some(parent) = Path::new(&datasites_mount_path).parent() {
            if let Some(parent_str) = parent.to_str() {
                datasites_mount_source = parent_str.to_string();
                datasites_volume_target = "/datasite-home".to_string();
                datasites_in_container = "/datasite-home/datasites".to_string();
                // Docker Desktop can reject a direct bind mount to paths containing '@'
                // (e.g. ".../sandbox/client1@sandbox.local"). In that case, mount the
                // parent sandbox root and remap datasites to the current datasite path.
                let parent_has_at = parent
                    .file_name()
                    .and_then(|name| name.to_str())
                    .map(|name| name.contains('@'))
                    .unwrap_or(false);
                if parent_has_at {
                    if let Some(sandbox_root) = parent.parent().and_then(|p| p.to_str()) {
                        if let Some(current_ds) = env_map.get("BV_CURRENT_DATASITE") {
                            datasites_mount_source = sandbox_root.to_string();
                            datasites_volume_target = "/datasite-root".to_string();
                            datasites_in_container =
                                format!("/datasite-root/{}/datasites", current_ds);
                            println!(
                                "  Docker mount fallback: using sandbox root {} (datasite {})",
                                datasites_mount_source, current_ds
                            );
                        }
                    }
                }
            }
        }
    }
    fs::create_dir_all(&datasites_mount_source).with_context(|| {
        format!(
            "Failed to create resolved datasites mount source: {}",
            datasites_mount_source
        )
    })?;

    // Fail-fast: verify container→host proxy connectivity before launching the real container.
    let tcp_proxy_enabled = env_map
        .get("SEQURE_TCP_PROXY")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    if tcp_proxy_enabled {
        if let Some(cp_ips) = env_map.get("SEQURE_CP_IPS") {
        if let Some(proxy_host) = cp_ips.split(',').next() {
            let comm_port = env_map
                .get("SEQURE_COMMUNICATION_PORT")
                .map(|s| s.as_str())
                .unwrap_or("9001");
            println!(
                "  Preflight: checking container→host connectivity {}:{}...",
                proxy_host, comm_port
            );
            let probe = Command::new(&container_runtime)
                .args([
                    "run",
                    "--rm",
                    "alpine:3.19",
                    "sh",
                    "-c",
                    &format!("nc -z -w3 '{}' '{}' 2>&1", proxy_host, comm_port),
                ])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .output();
            match probe {
                Ok(output) if output.status.success() => {
                    println!("  Preflight: {}:{} reachable ✓", proxy_host, comm_port);
                }
                Ok(output) => {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    let stdout = String::from_utf8_lossy(&output.stdout);
                    let combined = format!("{}\n{}", stdout, stderr);
                    if combined.to_ascii_lowercase().contains("network is unreachable") {
                        return Err(anyhow::anyhow!(
                            "Container→host proxy route unreachable for {}:{} ({}). Set BV_SYQURE_CP_HOST to a reachable host IP.",
                            proxy_host,
                            comm_port,
                            combined.trim()
                        )
                        .into());
                    }
                    eprintln!(
                        "  Preflight WARNING: {}:{} not reachable from container (proxy may not be listening yet)",
                        proxy_host, comm_port
                    );
                }
                Err(err) => {
                    eprintln!(
                        "  Preflight WARNING: failed to probe {}:{} from container: {}",
                        proxy_host, comm_port, err
                    );
                }
            }
        }
    }
    }

    let mut cmd = Command::new(&container_runtime);
    super::configure_child_process(&mut cmd);
    cmd.args(["run", "--name", &container_name]);
    if !keep_containers {
        cmd.arg("--rm");
    }
    if let Some(network_name) = docker_network_name {
        cmd.args(["--network", network_name]);
    }
    if let Some(party_ip) = docker_party_ip {
        cmd.args(["--ip", party_ip]);
    }

    // Add --platform for docker when specified (podman handles this differently)
    // Syqure containers are x86-only, so platform will typically be "linux/amd64"
    // This uses QEMU emulation on ARM64 systems
    if container_runtime == "docker" && !platform.is_empty() {
        cmd.args(["--platform", platform]);
    }
    if container_runtime == "docker" && cfg!(target_os = "linux") {
        // Ensure Docker host aliases are available on Linux too.
        cmd.args(["--add-host", "host.docker.internal:host-gateway"]);
        cmd.args(["--add-host", "gateway.docker.internal:host-gateway"]);
    }

    for (k, v) in env_map {
        if k == "SEQURE_DATASITES_ROOT" {
            cmd.args(["-e", &format!("{}={}", k, datasites_in_container)]);
        } else if let Some(mapped) = map_env_path(
            v,
            &module_root_for_match,
            &results_root_for_match,
            &datasites_root_for_match,
            &datasites_in_container,
        ) {
            cmd.args(["-e", &format!("{}={}", k, mapped)]);
        } else {
            cmd.args(["-e", &format!("{}={}", k, v)]);
        }
    }

    // Strip Windows extended path prefix for Docker volume mounts
    let module_mount_path = strip_extended_path_prefix(module_path_abs.to_string_lossy().as_ref());
    cmd.args(["-v", &format!("{}:/workspace/project", module_mount_path)]);
    cmd.args([
        "-v",
        &format!("{}:{}", datasites_mount_source, datasites_volume_target),
    ]);
    cmd.args([
        "-v",
        &format!("{}:{}", results_mount_path, results_in_container),
    ]);

    cmd.arg(image);
    cmd.arg("syqure");
    cmd.arg(&container_entrypoint);
    cmd.args(["--", &party_id.to_string()]);

    let display_cmd = format_command(&cmd);
    println!("\n▶️  Executing syqure in Docker...");
    println!("  {}", display_cmd.dimmed());

    let status = cmd.status().context("Failed to execute syqure in Docker")?;
    if !status.success() {
        if keep_containers {
            eprintln!("Syqure container {} failed; dumping logs:", container_name);
            let mut logs_cmd = Command::new(&container_runtime);
            super::configure_child_process(&mut logs_cmd);
            let _ = logs_cmd
                .args(["logs", "--tail", "200", &container_name])
                .status();
        } else {
            eprintln!(
                "Syqure container exited with code {:?}. Re-run with BIOVAULT_SYQURE_KEEP_CONTAINERS=1 to keep containers for logs.",
                status.code()
            );
        }
        return Err(anyhow::anyhow!(
            "Syqure Docker container exited with code: {:?}",
            status.code()
        )
        .into());
    }

    Ok(())
}

#[derive(Debug)]
struct ParsedArgs {
    inputs: HashMap<String, InputArg>,
    params: HashMap<String, String>,
    passthrough: Vec<String>,
    nextflow_max_forks: Option<u32>,
}

#[derive(Debug)]
struct InputArg {
    value: String,
    format_override: Option<String>,
}

fn parse_cli_args(args: &[String]) -> Result<ParsedArgs> {
    let mut inputs = HashMap::new();
    let mut params = HashMap::new();
    let mut format_overrides = HashMap::new();
    let mut passthrough = Vec::new();
    let mut nextflow_max_forks: Option<u32> = None;

    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];

        if arg == "--" {
            passthrough.extend(args[i + 1..].iter().cloned());
            break;
        }

        if !arg.starts_with("--") {
            passthrough.push(arg.clone());
            i += 1;
            continue;
        }

        let key = arg.strip_prefix("--").unwrap();

        if key == "set" {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];

            let (target, val) = value.split_once('=').ok_or_else(|| {
                anyhow::anyhow!(
                    "Invalid --set assignment '{}'. Use inputs.name=value or params.name=value.",
                    value
                )
            })?;

            if let Some(input_name) = target.strip_prefix("inputs.") {
                inputs.insert(
                    input_name.to_string(),
                    InputArg {
                        value: val.to_string(),
                        format_override: None,
                    },
                );
            } else if let Some(param_name) = target.strip_prefix("params.") {
                params.insert(param_name.to_string(), val.to_string());
            } else if let Some(param_name) = target.strip_prefix("param.") {
                params.insert(param_name.to_string(), val.to_string());
            } else {
                return Err(anyhow::anyhow!(
                    "Unsupported --set target '{}'. Expected inputs.<name> or params.<name>.",
                    target
                )
                .into());
            }
            i += 2;
            continue;
        }

        if key == "nxf-max-forks" || key == "nextflow-max-forks" {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = args[i + 1].trim();
            let parsed = value.parse::<u32>().ok().filter(|v| *v > 0);
            if parsed.is_none() {
                return Err(anyhow::anyhow!(
                    "Invalid value for {} (expected positive integer): {}",
                    arg,
                    value
                )
                .into());
            }
            nextflow_max_forks = parsed;
            i += 2;
            continue;
        }

        if key.starts_with("nxf-max-forks=") || key.starts_with("nextflow-max-forks=") {
            let value = key.split_once('=').map(|x| x.1).unwrap_or("").trim();
            let parsed = value.parse::<u32>().ok().filter(|v| *v > 0);
            if parsed.is_none() {
                return Err(anyhow::anyhow!(
                    "Invalid value for {} (expected positive integer): {}",
                    arg,
                    value
                )
                .into());
            }
            nextflow_max_forks = parsed;
            i += 1;
            continue;
        }

        if key.starts_with("param.") {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];
            let param_name = key.strip_prefix("param.").unwrap();
            params.insert(param_name.to_string(), value.clone());
            i += 2;
            continue;
        }

        if key.contains(".format") {
            if i + 1 >= args.len() {
                return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
            }
            let value = &args[i + 1];
            let input_name = key.strip_suffix(".format").unwrap();
            format_overrides.insert(input_name.to_string(), value.clone());
            i += 2;
            continue;
        }

        if key.contains(".mapping.") {
            // Future: support inline mapping overrides
            return Err(
                anyhow::anyhow!("Inline mapping overrides not yet supported: {}", key).into(),
            );
        }

        match key {
            "results-dir" | "results_dir" => {
                i += 2;
                continue;
            }
            _ => {}
        }

        if i + 1 >= args.len() {
            return Err(anyhow::anyhow!("Missing value for argument: {}", arg).into());
        }

        let value = &args[i + 1];
        inputs.insert(
            key.to_string(),
            InputArg {
                value: value.clone(),
                format_override: None,
            },
        );

        i += 2;
    }

    for (input_name, format) in &format_overrides {
        if let Some(input) = inputs.get_mut(input_name) {
            input.format_override = Some(format.clone());
        }
    }

    Ok(ParsedArgs {
        inputs,
        params,
        passthrough,
        nextflow_max_forks,
    })
}

fn parse_nextflow_max_forks_env() -> Result<Option<u32>> {
    let raw = match env::var("BIOVAULT_NXF_MAX_FORKS") {
        Ok(value) => value,
        Err(_) => return Ok(None),
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let parsed = trimmed.parse::<u32>().ok().filter(|v| *v > 0);
    if parsed.is_none() {
        return Err(anyhow::anyhow!(
            "Invalid BIOVAULT_NXF_MAX_FORKS (expected positive integer): {}",
            trimmed
        )
        .into());
    }
    Ok(parsed)
}

fn validate_no_clashes(spec: &ModuleSpec, parsed: &ParsedArgs) -> Result<()> {
    let input_names: Vec<&str> = spec.inputs.iter().map(|i| i.name.as_str()).collect();
    let output_names: Vec<&str> = spec.outputs.iter().map(|o| o.name.as_str()).collect();

    for param_name in parsed.params.keys() {
        if input_names.contains(&param_name.as_str()) {
            return Err(anyhow::anyhow!(
                "Parameter '{}' clashes with input name. Use --param.{} instead.",
                param_name,
                param_name
            )
            .into());
        }
        if output_names.contains(&param_name.as_str()) {
            return Err(anyhow::anyhow!(
                "Parameter '{}' clashes with output name. Use --param.{} instead.",
                param_name,
                param_name
            )
            .into());
        }
    }

    for input_name in parsed.inputs.keys() {
        if !input_names.contains(&input_name.as_str()) {
            println!(
                "⚠️  Warning: Unknown input '{}'. Expected inputs: {}",
                input_name.yellow(),
                input_names.join(", ")
            );
        }
    }

    Ok(())
}

fn build_inputs_json(
    spec: &ModuleSpec,
    parsed: &ParsedArgs,
    _project_path: &Path,
) -> Result<HashMap<String, JsonValue>> {
    let mut inputs_json = HashMap::new();

    for input_spec in &spec.inputs {
        if let Some(input_arg) = parsed.inputs.get(&input_spec.name) {
            let path_str = &input_arg.value;
            let path = Path::new(path_str);

            if !path.exists() {
                return Err(anyhow::anyhow!("Input file not found: {}", path_str).into());
            }

            let format = input_arg
                .format_override
                .as_deref()
                .or(input_spec.format.as_deref())
                .or_else(|| detect_format(path))
                .unwrap_or("unknown");

            inputs_json.insert(
                input_spec.name.clone(),
                json!({
                    "path": path.canonicalize()?.to_string_lossy().to_string(),
                    "type": input_spec.raw_type,
                    "format": format,
                    "mapping": input_spec.mapping,
                }),
            );
        } else if !input_spec.raw_type.ends_with('?') {
            return Err(
                anyhow::anyhow!("Required input '{}' not provided", input_spec.name).into(),
            );
        }
    }

    Ok(inputs_json)
}

fn build_params_json(spec: &ModuleSpec, parsed: &ParsedArgs) -> Result<HashMap<String, JsonValue>> {
    let mut params_json = HashMap::new();

    for param_spec in &spec.parameters {
        let value = if let Some(v) = parsed.params.get(&param_spec.name) {
            match param_spec.raw_type.as_str() {
                "Bool" => {
                    let bool_val = v.parse::<bool>().context(format!(
                        "Parameter '{}' expects Bool, got '{}'",
                        param_spec.name, v
                    ))?;
                    json!(bool_val)
                }
                "String" => json!(v),
                ty if ty.starts_with("Enum") => json!(v),
                unsupported => {
                    return Err(
                        anyhow::anyhow!("Unsupported parameter type: {}", unsupported).into(),
                    );
                }
            }
        } else if let Some(default) = &param_spec.default {
            serde_json::to_value(default)
                .context("Failed to convert default param value to JSON")?
        } else {
            continue;
        };

        params_json.insert(param_spec.name.clone(), value);
    }

    Ok(params_json)
}

fn detect_format(path: &Path) -> Option<&'static str> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .and_then(|ext| match ext.to_lowercase().as_str() {
            "json" => Some("json"),
            "yaml" | "yml" => Some("yaml"),
            "csv" => Some("csv"),
            "tsv" => Some("tsv"),
            "vcf" | "vcf.gz" => Some("vcf"),
            _ => None,
        })
}

/// Check if Podman is using the Hyper-V backend on Windows.
/// This is detected by checking the CONTAINERS_MACHINE_PROVIDER env var
/// or by running `podman machine inspect` to check the VM type.
#[cfg(target_os = "windows")]
fn is_podman_hyperv(docker_bin: &str, using_podman: bool) -> bool {
    if !using_podman {
        return false;
    }

    // Check environment variable first
    if let Ok(provider) = std::env::var("CONTAINERS_MACHINE_PROVIDER") {
        return provider.to_lowercase() == "hyperv";
    }

    // Try to detect from podman machine info
    let mut cmd = if docker_bin.contains("podman") {
        Command::new(docker_bin)
    } else {
        Command::new("podman")
    };
    super::configure_child_process(&mut cmd);
    let output = cmd.arg("machine").arg("inspect").output();

    if let Ok(output) = output {
        let stdout = String::from_utf8_lossy(&output.stdout);
        // Hyper-V machines have ConfigDir containing "hyperv"
        return stdout.contains("hyperv");
    }

    false
}

#[cfg(not(target_os = "windows"))]
fn is_podman_hyperv(_docker_bin: &str, _using_podman: bool) -> bool {
    false
}

#[cfg(target_os = "windows")]
fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<()> {
    fs::create_dir_all(dst)
        .with_context(|| format!("Failed to create directory: {}", dst.display()))?;
    for entry in
        fs::read_dir(src).with_context(|| format!("Failed to read directory: {}", src.display()))?
    {
        let entry = entry.with_context(|| "Failed to read directory entry")?;
        let path = entry.path();
        let file_name = entry.file_name();
        let dest_path = dst.join(file_name);
        if path.is_dir() {
            copy_dir_recursive(&path, &dest_path)?;
        } else {
            fs::copy(&path, &dest_path).with_context(|| {
                format!(
                    "Failed to copy file: {} -> {}",
                    path.display(),
                    dest_path.display()
                )
            })?;
        }
    }
    Ok(())
}

#[cfg(target_os = "windows")]
fn copy_path_to_path(src: &Path, dest_path: &Path) -> Result<()> {
    if let Some(parent) = dest_path.parent() {
        fs::create_dir_all(parent).with_context(|| {
            format!(
                "Failed to create destination dir for {}",
                dest_path.display()
            )
        })?;
    }
    if src.is_dir() {
        copy_dir_recursive(src, dest_path)?;
    } else {
        fs::copy(src, dest_path).with_context(|| {
            format!(
                "Failed to copy file: {} -> {}",
                src.display(),
                dest_path.display()
            )
        })?;
    }
    Ok(())
}

/// Rewrite a CSV file to use host-local flat paths (e.g., C:/bvtemp/...).
#[cfg(target_os = "windows")]
fn rewrite_csv_for_flat_dir(
    csv_path: &Path,
    flat_data_dir: &str,
    path_map: &BTreeMap<String, String>,
) -> Result<()> {
    append_desktop_log(&format!(
        "[CSV Rewrite Flat] Processing: {}",
        csv_path.display()
    ));
    let mut converted_count = 0;
    let mut missing_count = 0;

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(csv_path)
        .with_context(|| format!("Failed to open CSV: {}", csv_path.display()))?;
    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .from_writer(vec![]);

    for result in reader.records() {
        let record = result
            .with_context(|| format!("Failed to parse CSV record: {}", csv_path.display()))?;
        let mut out: Vec<String> = Vec::with_capacity(record.len());
        for field in record.iter() {
            if looks_like_windows_absolute_path(field) {
                let normalized = normalize_windows_path_str(field);
                if let Some(mapped) = path_map.get(&normalized) {
                    converted_count += 1;
                    out.push(mapped.replace('\\', "/"));
                } else {
                    missing_count += 1;
                    let file_name = Path::new(&normalized)
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy();
                    out.push(format!("{}/{}", flat_data_dir, file_name));
                }
            } else {
                out.push(field.to_string());
            }
        }
        writer
            .write_record(out)
            .with_context(|| format!("Failed to write CSV record: {}", csv_path.display()))?;
    }

    let new_content = String::from_utf8(
        writer
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to finalize CSV writer: {}", e))?,
    )
    .context("Failed to encode rewritten CSV")?;

    if missing_count > 0 {
        append_desktop_log(&format!(
            "[CSV Rewrite Flat] Warning: {} paths not in staging map for {}",
            missing_count,
            csv_path.display()
        ));
    }

    append_desktop_log(&format!(
        "[CSV Rewrite Flat] Converted {} paths to flat dir in {}",
        converted_count,
        csv_path.display()
    ));
    fs::write(csv_path, new_content)
        .with_context(|| format!("Failed to write converted CSV: {}", csv_path.display()))?;

    Ok(())
}

/// Extract Windows file paths from CSV content (keeps file paths, not parent dirs).
#[cfg(target_os = "windows")]
fn extract_files_from_csv(content: &str, files: &mut Vec<PathBuf>) {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_reader(content.as_bytes());
    for record in reader.records().flatten() {
        for field in record.iter() {
            let field = field.trim();
            if looks_like_windows_absolute_path(field) {
                let normalized = normalize_windows_path_str(field);
                let path = Path::new(&normalized);
                if path.exists() {
                    files.push(path.to_path_buf());
                }
            }
        }
    }
}

/// Extract Windows file paths from a JSON value (includes CSV files and entries inside them).
#[cfg(target_os = "windows")]
fn extract_files_from_json(value: &JsonValue, files: &mut Vec<PathBuf>) {
    match value {
        JsonValue::String(s) => {
            if looks_like_windows_absolute_path(s) {
                let normalized = normalize_windows_path_str(s);
                let path = Path::new(&normalized);
                if path.exists() {
                    files.push(path.to_path_buf());
                    if s.to_lowercase().ends_with(".csv") {
                        if let Ok(content) = fs::read_to_string(path) {
                            extract_files_from_csv(&content, files);
                        }
                    }
                }
            }
        }
        JsonValue::Object(map) => {
            for v in map.values() {
                extract_files_from_json(v, files);
            }
        }
        JsonValue::Array(arr) => {
            for v in arr {
                extract_files_from_json(v, files);
            }
        }
        _ => {}
    }
}

/// Get the default Podman machine name from `podman system connection list`.
#[cfg(target_os = "windows")]
fn get_default_podman_machine(docker_bin: &str) -> Option<String> {
    let mut cmd = Command::new(docker_bin);
    super::configure_child_process(&mut cmd);
    cmd.args([
        "system",
        "connection",
        "list",
        "--format",
        "{{.Name}}\t{{.Default}}",
    ]);

    let output = cmd.output().ok()?;
    if !output.status.success() {
        return None;
    }

    for line in String::from_utf8_lossy(&output.stdout).lines() {
        let mut parts = line.split('\t');
        let name = parts.next()?.trim();
        let is_default = parts.next().unwrap_or("").trim();
        if is_default == "true" {
            let name = name.trim_end_matches("-root");
            if !name.is_empty() {
                return Some(name.to_string());
            }
        }
    }

    None
}

#[cfg(target_os = "windows")]
fn podman_machine_ssh_cmd(docker_bin: &str) -> Command {
    let mut cmd = Command::new(docker_bin);
    super::configure_child_process(&mut cmd);
    cmd.arg("machine").arg("ssh");
    if let Some(name) = get_default_podman_machine(docker_bin) {
        cmd.arg(name);
    }
    cmd.arg("--");
    cmd
}

/// Copy a file into the Podman VM using stdin pipe.
/// This works around the 9P mount issues on Hyper-V.
#[cfg(target_os = "windows")]
fn copy_file_to_vm(docker_bin: &str, local_path: &Path, vm_path: &str) -> Result<()> {
    use std::io::Read;

    let mut file = fs::File::open(local_path).map_err(|e| {
        anyhow::anyhow!(
            "Failed to open file for VM copy: {}: {}",
            local_path.display(),
            e
        )
    })?;

    let mut contents = Vec::new();
    file.read_to_end(&mut contents)
        .map_err(|e| anyhow::anyhow!("Failed to read file: {}: {}", local_path.display(), e))?;

    let mut cmd = podman_machine_ssh_cmd(docker_bin);
    let mut child = cmd
        .arg(format!("cat > '{}'", vm_path))
        .stdin(Stdio::piped())
        .spawn()
        .map_err(|e| anyhow::anyhow!("Failed to spawn podman machine ssh: {}", e))?;

    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(&contents)
            .map_err(|e| anyhow::anyhow!("Failed to write to VM: {}", e))?;
    }

    let status = child
        .wait()
        .map_err(|e| anyhow::anyhow!("Failed to wait for VM copy: {}", e))?;
    if !status.success() {
        return Err(anyhow::anyhow!("Failed to copy file to VM: {}", local_path.display()).into());
    }

    Ok(())
}

/// Copy a directory recursively into the Podman VM.
#[cfg(target_os = "windows")]
fn copy_dir_to_vm(docker_bin: &str, local_dir: &Path, vm_dir: &str) -> Result<()> {
    // Create the directory in the VM
    let status = podman_machine_ssh_cmd(docker_bin)
        .arg(format!("mkdir -p '{}' && chmod 777 '{}'", vm_dir, vm_dir))
        .status()
        .map_err(|e| anyhow::anyhow!("Failed to create directory in VM: {}", e))?;

    if !status.success() {
        return Err(anyhow::anyhow!("Failed to create directory in VM: {}", vm_dir).into());
    }

    // Copy each file/subdir
    for entry in fs::read_dir(local_dir)
        .map_err(|e| anyhow::anyhow!("Failed to read directory: {}: {}", local_dir.display(), e))?
    {
        let entry = entry.map_err(|e| anyhow::anyhow!("Failed to read directory entry: {}", e))?;
        let path = entry.path();
        let name = entry.file_name();
        let vm_path = format!("{}/{}", vm_dir, name.to_string_lossy());

        if path.is_dir() {
            copy_dir_to_vm(docker_bin, &path, &vm_path)?;
        } else {
            copy_file_to_vm(docker_bin, &path, &vm_path)?;
        }
    }

    // Make files world-readable for nested containers
    let _ = podman_machine_ssh_cmd(docker_bin)
        .arg(format!("chmod -R 777 '{}'", vm_dir))
        .status();

    Ok(())
}

/// Copy a file or directory from the Podman VM back to Windows.
#[cfg(target_os = "windows")]
fn copy_from_vm(docker_bin: &str, vm_path: &str, local_path: &Path) -> Result<()> {
    // Create parent directory if needed
    if let Some(parent) = local_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| anyhow::anyhow!("Failed to create parent directory: {}", e))?;
    }

    // Use tar to copy directory contents
    let output = podman_machine_ssh_cmd(docker_bin)
        .arg(format!(
            "if [ -d '{}' ]; then cd '{}' && tar -cf - .; else cat '{}'; fi",
            vm_path, vm_path, vm_path
        ))
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to read from VM: {}", e))?;

    if !output.status.success() {
        // Directory might be empty or not exist, that's ok for results
        return Ok(());
    }

    // Check if it's a tar archive (starts with tar magic) or plain file
    if output.stdout.len() >= 262 && &output.stdout[257..262] == b"ustar" {
        // It's a tar archive, extract it
        use std::io::Cursor;
        let cursor = Cursor::new(output.stdout);
        let mut archive = tar::Archive::new(cursor);
        archive
            .unpack(local_path)
            .map_err(|e| anyhow::anyhow!("Failed to extract tar archive: {}", e))?;
    } else {
        // Plain file
        fs::write(local_path, &output.stdout).map_err(|e| {
            anyhow::anyhow!("Failed to write file: {}: {}", local_path.display(), e)
        })?;
    }

    Ok(())
}

/// Create a temporary directory in the Podman VM and return its path.
#[cfg(target_os = "windows")]
fn create_vm_temp_dir(docker_bin: &str, prefix: &str) -> Result<String> {
    let output = podman_machine_ssh_cmd(docker_bin)
        .arg(format!(
            "mktemp -d /tmp/{}-XXXXXX && chmod 777 /tmp/{}-*",
            prefix, prefix
        ))
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to create temp directory in VM: {}", e))?;

    if !output.status.success() {
        return Err(anyhow::anyhow!(
            "Failed to create temp directory in VM: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    let vm_path = String::from_utf8_lossy(&output.stdout).trim().to_string();
    Ok(vm_path)
}

/// Clean up a directory in the Podman VM.
#[cfg(target_os = "windows")]
fn cleanup_vm_dir(docker_bin: &str, vm_dir: &str) {
    let _ = podman_machine_ssh_cmd(docker_bin)
        .arg(format!("sudo rm -rf '{}'", vm_dir))
        .status();
}

/// Generate a Nextflow config file for the specified container runtime.
/// Returns the path to the generated config file.
///
/// Platform handling:
/// - ARM64 with all arm64-capable containers: No --platform (native execution)
/// - ARM64 with any x86-only container: --platform linux/amd64 (explicit emulation)
/// - x86_64 or force flag set: --platform linux/amd64 for compatibility
fn generate_runtime_config(
    project_path: &Path,
    using_podman: bool,
    stage_in_copy: bool,
    max_forks: Option<u32>,
) -> Result<PathBuf> {
    let config_path = project_path.join(".biovault-runtime.config");

    // Check for workflow.nf to inspect containers
    let workflow_path = project_path.join("workflow.nf");

    // Determine if we should force x86 platform for task containers
    // Force x86 if:
    // 1. Environment variable is set, OR
    // 2. Not on ARM64 (x86 host), OR
    // 3. On ARM64 but any container lacks ARM64 support
    let containers_need_x86 = if is_arm64() && workflow_path.exists() && !using_podman {
        any_container_lacks_arm64(&workflow_path)
    } else {
        false
    };
    let force_x86 = should_force_x86_containers() || !is_arm64() || containers_need_x86;

    let config_contents = if using_podman {
        // Podman configuration
        // - Use podman.enabled instead of docker.enabled
        // - Use /bin/sh for alpine-based containers (no bash)
        // - Add security-opt to disable SELinux labeling for nested containers
        let mut config = r#"process.executor = 'local'
podman.enabled = true
process.shell = ['/bin/sh', '-ue']
podman.runOptions = '--security-opt label=disable'
"#
        .to_string();
        // Note: we intentionally do not add --platform for podman yet because it is harder to
        // validate across environments. Add it once we can test reliably.
        if stage_in_copy {
            config.push_str("\nprocess.stageInMode = 'copy'\n");
        }
        if let Some(value) = max_forks {
            config.push_str(&format!("\nprocess.maxForks = {}\n", value));
        }
        config
    } else {
        // Docker configuration
        // On ARM64: Don't force platform - let Docker auto-select (native if available, emulate if not)
        // On x86_64 or when forced: Use --platform linux/amd64 for compatibility
        let mut config = if force_x86 {
            r#"process.executor = 'local'
docker.enabled = true
docker.runOptions = '--platform linux/amd64 -u $(id -u):$(id -g)'
"#
            .to_string()
        } else {
            // ARM64: Let Docker auto-select the best platform
            // Multi-arch images will use arm64, x86-only images will be emulated
            r#"process.executor = 'local'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
"#
            .to_string()
        };
        if let Some(value) = max_forks {
            config.push_str(&format!("\nprocess.maxForks = {}\n", value));
        }
        config
    };

    fs::write(&config_path, config_contents).context("Failed to write runtime nextflow.config")?;

    append_desktop_log(&format!(
        "[Pipeline] Generated runtime config at {} (podman={}, force_x86={}, arch={})",
        config_path.display(),
        using_podman,
        force_x86,
        std::env::consts::ARCH
    ));

    Ok(config_path)
}

fn install_dynamic_template(biovault_home: &Path) -> Result<()> {
    let env_dir = biovault_home.join("env").join("dynamic-nextflow");
    if !env_dir.exists() {
        fs::create_dir_all(&env_dir).context("Failed to create dynamic template directory")?;
    }

    let template_path = env_dir.join("template.nf");
    let template_contents = include_str!("../../templates/dynamic/template.nf");
    fs::write(&template_path, template_contents)
        .context("Failed to install dynamic template.nf")?;

    println!("📦 Dynamic template ready at {}", template_path.display());

    // Note: We no longer install a static nextflow.config here.
    // Instead, we generate a runtime-specific config in generate_runtime_config()
    // based on whether Docker or Podman is being used.

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::module_spec::{InputSpec, ModuleSpec};
    use tempfile::TempDir;

    fn sample_project_spec() -> ModuleSpec {
        ModuleSpec {
            name: "test".to_string(),
            author: "author".to_string(),
            workflow: "workflow.nf".to_string(),
            description: None,
            runtime: Some("nextflow".to_string()),
            version: None,
            datasites: None,
            env: Default::default(),
            assets: vec![],
            parameters: vec![],
            inputs: vec![
                InputSpec {
                    name: "samplesheet".to_string(),
                    raw_type: "File".to_string(),
                    description: None,
                    format: Some("csv".to_string()),
                    path: None,
                    mapping: None,
                },
                InputSpec {
                    name: "data_dir".to_string(),
                    raw_type: "Directory".to_string(),
                    description: None,
                    format: None,
                    path: None,
                    mapping: None,
                },
            ],
            outputs: vec![],
            steps: Vec::new(),
            runner: None,
        }
    }

    #[test]
    fn build_inputs_json_handles_file_and_directory() {
        let tmp = TempDir::new().unwrap();
        let file_path = tmp.path().join("participants.csv");
        std::fs::write(&file_path, "id,path\n1,a.txt\n").unwrap();
        let dir_path = tmp.path().join("data");
        std::fs::create_dir_all(&dir_path).unwrap();

        let parsed = ParsedArgs {
            inputs: HashMap::from([
                (
                    "samplesheet".to_string(),
                    InputArg {
                        value: file_path.to_string_lossy().to_string(),
                        format_override: None,
                    },
                ),
                (
                    "data_dir".to_string(),
                    InputArg {
                        value: dir_path.to_string_lossy().to_string(),
                        format_override: None,
                    },
                ),
            ]),
            params: HashMap::new(),
            passthrough: Vec::new(),
            nextflow_max_forks: None,
        };

        let project_spec = sample_project_spec();
        let inputs = build_inputs_json(&project_spec, &parsed, tmp.path()).unwrap();

        let sheet_entry = inputs.get("samplesheet").expect("samplesheet entry");
        assert_eq!(sheet_entry["type"], json!("File"));
        assert_eq!(sheet_entry["format"], json!("csv"));
        let sheet_path = sheet_entry["path"].as_str().unwrap();
        assert_eq!(
            sheet_path,
            file_path.canonicalize().unwrap().to_string_lossy()
        );

        let dir_entry = inputs.get("data_dir").expect("data_dir entry");
        assert_eq!(dir_entry["type"], json!("Directory"));
        let dir_json_path = dir_entry["path"].as_str().unwrap();
        assert_eq!(
            dir_json_path,
            dir_path.canonicalize().unwrap().to_string_lossy()
        );
    }

    #[test]
    fn parse_cli_args_supports_set_inputs_and_params() {
        let args = vec![
            "--set".to_string(),
            "inputs.samplesheet=/tmp/sheet.csv".to_string(),
            "--set".to_string(),
            "params.threshold=0.5".to_string(),
        ];

        let parsed = parse_cli_args(&args).expect("parse --set inputs");

        let sheet = parsed
            .inputs
            .get("samplesheet")
            .expect("samplesheet input parsed");
        assert_eq!(sheet.value, "/tmp/sheet.csv");

        let threshold = parsed
            .params
            .get("threshold")
            .expect("param threshold parsed");
        assert_eq!(threshold, "0.5");
        assert!(parsed.passthrough.is_empty());
    }

    #[test]
    fn parse_cli_args_ignores_results_dir() {
        let args = vec![
            "--results-dir".to_string(),
            "custom_results".to_string(),
            "--samplesheet".to_string(),
            "/tmp/sheet.csv".to_string(),
        ];

        let parsed = parse_cli_args(&args).unwrap();
        assert!(parsed.inputs.contains_key("samplesheet"));
        assert!(!parsed.inputs.contains_key("results-dir"));
        assert!(parsed.passthrough.is_empty());
    }

    #[test]
    fn parse_cli_args_captures_nextflow_flags() {
        let args = vec![
            "--samplesheet".to_string(),
            "/tmp/sheet.csv".to_string(),
            "-with-singularity".to_string(),
            "-profile".to_string(),
            "docker".to_string(),
        ];

        let parsed = parse_cli_args(&args).expect("parse passthrough flags");
        assert_eq!(
            parsed.passthrough,
            vec![
                "-with-singularity".to_string(),
                "-profile".to_string(),
                "docker".to_string()
            ]
        );
    }
}
