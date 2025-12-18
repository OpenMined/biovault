#!/usr/bin/env bash
set -euo pipefail

# GitHub Actions Windows runners often provide `python` but not `python3` on PATH.
# Normalize so the rest of the script can keep using `python3`.
if ! command -v python3 >/dev/null 2>&1; then
  if command -v python >/dev/null 2>&1; then
    python3() { python "$@"; }
  else
    echo "Missing required tool: python3 (or python)" >&2
    exit 1
  fi
fi

usage() {
  cat <<'EOF'
Usage: ./devstack.sh [--clients email1,email2] [--sandbox DIR] [--reset] [--stop] [--status]

Starts or stops the SyftBox devstack (server + clients) without Docker using the
sbdev tool from the syftbox submodule. Defaults to the sandbox clients already
used in the scenario tests.

Options:
  --clients list   Comma-separated client emails (default: client1@sandbox.local,client2@sandbox.local)
  --sandbox DIR    Sandbox root path (default: ./sandbox)
  --client-mode MODE  Client implementation: go|rust|mixed|embedded (default: go)
  --rust-client-bin PATH  Path to Rust syftbox client binary (default: syftbox/rust/target/release/syftbox-rs)
  --skip-rust-build Skip building the Rust client (assumes binary exists)
  --skip-client-daemons Do not launch syftbox client daemons (server+minio only; still writes per-client config.json)
  --reset          Remove any existing devstack state before starting (also removes sandbox on stop)
  --skip-sync-check Skip the sbdev sync probe after boot (faster, less safe)
  --skip-keys      Skip key generation (for manual key management testing)
  --stop           Stop the devstack instead of starting it
  --status         Print the current devstack state (relay/state.json) and exit
  -h, --help       Show this message

Environment (optional defaults when flags not provided):
  BV_DEVSTACK_CLIENT_MODE      go|rust|mixed|embedded
  BV_DEVSTACK_RUST_CLIENT_BIN  Path to Rust client binary
  BV_DEVSTACK_SKIP_RUST_BUILD  Set to 1 to skip building Rust client
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SYFTBOX_DIR="$ROOT_DIR/syftbox"
SANDBOX_DIR="${SANDBOX_DIR:-$ROOT_DIR/sandbox}"
GO_CACHE_DIR="$SYFTBOX_DIR/.gocache"
ACTION="start"
RESET_FLAG=0
SKIP_SYNC_CHECK=0
SKIP_KEYS=0
SKIP_CLIENT_DAEMONS=0
EMBEDDED_MODE=0
RAW_CLIENTS=()
CLIENT_MODE=""
CLIENT_MODE_EXPLICIT=0
RUST_CLIENT_BIN=""
SKIP_RUST_BUILD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clients)
      [[ $# -lt 2 ]] && { echo "Missing value for --clients" >&2; usage >&2; exit 1; }
      RAW_CLIENTS+=("$2")
      shift
      ;;
    --sandbox)
      [[ $# -lt 2 ]] && { echo "Missing value for --sandbox" >&2; usage >&2; exit 1; }
      SANDBOX_DIR="$2"
      shift
      ;;
    --client-mode)
      [[ $# -lt 2 ]] && { echo "Missing value for --client-mode" >&2; usage >&2; exit 1; }
      CLIENT_MODE="$2"
      CLIENT_MODE_EXPLICIT=1
      shift
      ;;
    --rust-client-bin)
      [[ $# -lt 2 ]] && { echo "Missing value for --rust-client-bin" >&2; usage >&2; exit 1; }
      RUST_CLIENT_BIN="$2"
      CLIENT_MODE_EXPLICIT=1
      if [[ -z "$CLIENT_MODE" ]]; then
        CLIENT_MODE="rust"
      fi
      shift
      ;;
    --skip-rust-build)
      SKIP_RUST_BUILD=1
      ;;
    --rust)
      CLIENT_MODE="rust"
      CLIENT_MODE_EXPLICIT=1
      ;;
    --mixed)
      CLIENT_MODE="mixed"
      CLIENT_MODE_EXPLICIT=1
      ;;
    --embedded)
      CLIENT_MODE="embedded"
      CLIENT_MODE_EXPLICIT=1
      ;;
    --reset)
      RESET_FLAG=1
      ;;
    --skip-sync-check)
      SKIP_SYNC_CHECK=1
      ;;
    --skip-keys)
      SKIP_KEYS=1
      ;;
    --skip-client-daemons)
      SKIP_CLIENT_DAEMONS=1
      ;;
    --stop)
      ACTION="stop"
      ;;
    --status)
      ACTION="status"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

if (( ! CLIENT_MODE_EXPLICIT )) && [[ -n "${BV_DEVSTACK_CLIENT_MODE:-}" ]]; then
  CLIENT_MODE="${BV_DEVSTACK_CLIENT_MODE}"
  CLIENT_MODE_EXPLICIT=1
fi

if [[ -z "$RUST_CLIENT_BIN" ]] && [[ -n "${BV_DEVSTACK_RUST_CLIENT_BIN:-}" ]]; then
  RUST_CLIENT_BIN="${BV_DEVSTACK_RUST_CLIENT_BIN}"
  CLIENT_MODE_EXPLICIT=1
  if [[ -z "$CLIENT_MODE" ]]; then
    CLIENT_MODE="rust"
  fi
fi

if (( ! SKIP_RUST_BUILD )) && [[ "${BV_DEVSTACK_SKIP_RUST_BUILD:-0}" == "1" ]]; then
  SKIP_RUST_BUILD=1
fi

abs_path() {
  python3 - <<'PY' "$1"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

SANDBOX_DIR="$(abs_path "$SANDBOX_DIR")"

require_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing required tool: $1" >&2; exit 1; }
}

require_file() {
  [[ -f "$1" ]] || { echo "Missing required file: $1" >&2; exit 1; }
}

is_windows() {
  case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) return 0 ;;
    *) return 1 ;;
  esac
}

declare -a CLIENTS=()
add_client() {
  local raw="$1"
  [[ -z "$raw" ]] && return
  raw="${raw#"${raw%%[![:space:]]*}"}"
  raw="${raw%"${raw##*[![:space:]]}"}"
  [[ -z "$raw" ]] && return
  raw="$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]')"
  for existing in "${CLIENTS[@]:-}"; do
    [[ "$existing" == "$raw" ]] && return
  done
  CLIENTS+=("$raw")
}

if ((${#RAW_CLIENTS[@]})); then
  for block in "${RAW_CLIENTS[@]}"; do
    IFS=',' read -r -a parts <<< "$block"
    for part in "${parts[@]}"; do
      add_client "$part"
    done
  done
fi

if ((${#CLIENTS[@]} == 0)); then
  add_client "client1@sandbox.local"
  add_client "client2@sandbox.local"
fi

BV_BIN="${BV_BIN:-$ROOT_DIR/cli/target/release/bv}"

ensure_bv_binary() {
  if [[ -x "$BV_BIN" ]]; then
    return
  fi
  require_bin cargo
  echo "Building BioVault CLI (release)..."
  (cd "$ROOT_DIR/cli" && cargo build --release)
  [[ -x "$BV_BIN" ]] || { echo "Failed to build BioVault CLI at $BV_BIN" >&2; exit 1; }
}

stop_stack() {
  echo "Stopping SyftBox devstack at $SANDBOX_DIR..."
  # Stop any Jupyter processes in the sandbox first
  pkill -f "jupyter.*$SANDBOX_DIR" 2>/dev/null || true

  # Stop embedded syftbox hosts (best effort).
  if [[ -x "$BV_BIN" ]]; then
    for email in "${CLIENTS[@]}"; do
      local client_dir="$SANDBOX_DIR/$email"
      if [[ -d "$client_dir" ]]; then
        HOME="$client_dir" "$BV_BIN" syftboxd stop >/dev/null 2>&1 || true
      fi
    done
  fi

  if [[ -d "$SYFTBOX_DIR" ]]; then
    (cd "$SYFTBOX_DIR" && GOCACHE="$GO_CACHE_DIR" go run ./cmd/devstack stop --path "$SANDBOX_DIR") || true
  fi
  if (( RESET_FLAG )); then
    rm -rf "$SANDBOX_DIR"
  fi
  echo "Devstack stop complete."
}

bootstrap_biovault() {
  local email="$1"
  local client_dir="$SANDBOX_DIR/$email"
  local data_dir="$client_dir"
  local config_path="$client_dir/.syftbox/config.json"

  [[ -d "$client_dir" ]] || { echo "Client directory not found: $client_dir" >&2; exit 1; }
  require_file "$config_path"

  echo "Configuring BioVault for $email"
  HOME="$client_dir" \
  SYFTBOX_EMAIL="$email" \
  SYFTBOX_DATA_DIR="$data_dir" \
  SYFTBOX_CONFIG_PATH="$config_path" \
  "$BV_BIN" init --quiet "$email" >/dev/null
}

configure_syftbox_clients() {
  if (( ! CLIENT_MODE_EXPLICIT )); then
    return 0
  fi

  local mode
  mode="$(printf '%s' "${CLIENT_MODE:-go}" | tr '[:upper:]' '[:lower:]')"
  case "$mode" in
    ""|go)
      unset SBDEV_CLIENT_BIN SBDEV_CLIENT_MODE SBDEV_RUST_CLIENT_BIN
      return 0
      ;;
    embedded)
      # Don't spawn syftbox client daemons; we'll start `bv syftboxd` per client instead.
      EMBEDDED_MODE=1
      SKIP_CLIENT_DAEMONS=1
      unset SBDEV_CLIENT_BIN SBDEV_CLIENT_MODE SBDEV_RUST_CLIENT_BIN
      return 0
      ;;
    rust|mixed)
      ;;
    *)
      echo "Invalid --client-mode: $mode (expected go|rust|mixed|embedded)" >&2
      exit 1
      ;;
  esac

  local rust_bin
  rust_bin="$RUST_CLIENT_BIN"
  if [[ -z "$rust_bin" ]]; then
    rust_bin="$SYFTBOX_DIR/rust/target/release/syftbox-rs"
    if is_windows; then
      rust_bin="${rust_bin}.exe"
    fi
  fi

  rust_bin="$(abs_path "$rust_bin")"

  if (( ! SKIP_RUST_BUILD )); then
    require_bin cargo
    echo "Building Rust SyftBox client..."
    (cd "$SYFTBOX_DIR/rust" && cargo build --release)
  fi

  [[ -f "$rust_bin" ]] || { echo "Rust SyftBox client binary not found at $rust_bin" >&2; exit 1; }

  if [[ "$mode" == "rust" ]]; then
    export SBDEV_CLIENT_BIN="$rust_bin"
    unset SBDEV_CLIENT_MODE SBDEV_RUST_CLIENT_BIN
  else
    unset SBDEV_CLIENT_BIN
    export SBDEV_CLIENT_MODE="$mode"
    export SBDEV_RUST_CLIENT_BIN="$rust_bin"
  fi

  echo "Using SyftBox Rust client ($mode): $rust_bin"
}

start_syftboxd() {
  ensure_bv_binary

  for email in "${CLIENTS[@]}"; do
    local client_dir="$SANDBOX_DIR/$email"
    local data_dir="$client_dir"
    local config_path="$client_dir/.syftbox/config.json"

    [[ -d "$client_dir" ]] || { echo "Client directory not found: $client_dir" >&2; exit 1; }
    require_file "$config_path"

    echo "Starting syftboxd (embedded) for $email"
    HOME="$client_dir" \
    SYFTBOX_EMAIL="$email" \
    SYFTBOX_DATA_DIR="$data_dir" \
    SYFTBOX_CONFIG_PATH="$config_path" \
    BV_SYFTBOX_BACKEND=embedded \
    "$BV_BIN" syftboxd start >/dev/null
  done
}

start_stack() {
  require_bin go
  [[ -d "$SYFTBOX_DIR" ]] || { echo "Missing syftbox checkout at $SYFTBOX_DIR" >&2; exit 1; }

  mkdir -p "$SANDBOX_DIR"

  # Resolve client-mode first because it may force --skip-client-daemons (embedded mode).
  # This must happen before we build the sbdev argument list below.
  configure_syftbox_clients

  local args=(--path "$SANDBOX_DIR" --random-ports)
  (( RESET_FLAG )) && args+=(--reset)
  (( SKIP_SYNC_CHECK )) && args+=(--skip-sync-check)
  (( SKIP_CLIENT_DAEMONS )) && args+=(--skip-client-daemons)
  for email in "${CLIENTS[@]}"; do
    args+=(--client "$email")
  done

  # Export environment variables for the devstack command
  # These will be inherited by the Go process and passed to spawned clients
  export GOCACHE="$GO_CACHE_DIR"
  [[ -n "${SCENARIO_JAVA_HOME:-}" ]] && export JAVA_HOME="$SCENARIO_JAVA_HOME" JAVA_CMD="$SCENARIO_JAVA_HOME/bin/java"
  [[ -n "${SCENARIO_USER_PATH:-}" ]] && export PATH="$SCENARIO_USER_PATH"
  [[ -n "${NXF_DISABLE_JAVA_VERSION_CHECK:-}" ]] && export NXF_DISABLE_JAVA_VERSION_CHECK
  [[ -n "${NXF_IGNORE_JAVA_VERSION:-}" ]] && export NXF_IGNORE_JAVA_VERSION
  [[ -n "${NXF_OPTS:-}" ]] && export NXF_OPTS

  echo "Starting SyftBox devstack via syftbox/cmd/devstack..."
  (cd "$SYFTBOX_DIR" && go run ./cmd/devstack start "${args[@]}")

  if (( SKIP_KEYS )); then
    echo "Skipping key generation (--skip-keys)"
  else
    ensure_bv_binary
    for email in "${CLIENTS[@]}"; do
      bootstrap_biovault "$email"
    done
  fi

  if (( EMBEDDED_MODE )); then
    start_syftboxd
  fi

  echo ""
  echo "Devstack ready at $SANDBOX_DIR"
  for email in "${CLIENTS[@]}"; do
    printf '  - %s\n' "$email"
  done
}

print_state() {
  local state_path="$SANDBOX_DIR/relay/state.json"
  if [[ ! -f "$state_path" ]]; then
    echo "State file not found at $state_path"
    exit 1
  fi
  echo "Devstack state at $state_path:"
  python3 - <<'PY' "$state_path"
import json, sys, pathlib
path = pathlib.Path(sys.argv[1])
data = json.loads(path.read_text())
json.dump(data, sys.stdout, indent=2, sort_keys=True)
print()
PY
}

if [[ "$ACTION" == "stop" ]]; then
  stop_stack
  exit 0
fi

if [[ "$ACTION" == "status" ]]; then
  print_state
  exit 0
fi

start_stack
