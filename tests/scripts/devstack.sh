#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./devstack.sh [--clients email1,email2] [--sandbox DIR] [--reset] [--stop] [--status]

Starts or stops the SyftBox devstack (server + clients) without Docker using the
sbdev tool from the syftbox submodule. Defaults to the sandbox clients already
used in the scenario tests.

Options:
  --clients list   Comma-separated client emails (default: client1@sandbox.local,client2@sandbox.local)
  --sandbox DIR    Sandbox root path (default: ./sandbox)
  --reset          Remove any existing devstack state before starting (also removes sandbox on stop)
  --stop           Stop the devstack instead of starting it
  --status         Print the current devstack state (relay/state.json) and exit
  -h, --help       Show this message
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SYFTBOX_DIR="$ROOT_DIR/syftbox"
SANDBOX_DIR="${SANDBOX_DIR:-$ROOT_DIR/sandbox}"
GO_CACHE_DIR="$SYFTBOX_DIR/.gocache"
ACTION="start"
RESET_FLAG=0
RAW_CLIENTS=()

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
    --reset)
      RESET_FLAG=1
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

start_stack() {
  require_bin go
  ensure_bv_binary
  [[ -d "$SYFTBOX_DIR" ]] || { echo "Missing syftbox checkout at $SYFTBOX_DIR" >&2; exit 1; }

  mkdir -p "$SANDBOX_DIR"
  local args=(--path "$SANDBOX_DIR" --random-ports)
  (( RESET_FLAG )) && args+=(--reset)
  for email in "${CLIENTS[@]}"; do
    args+=(--client "$email")
  done

  # Pass through Java/Nextflow environment variables if available
  [[ -n "${SCENARIO_JAVA_HOME:-}" ]] && args+=(--env "JAVA_HOME=$SCENARIO_JAVA_HOME")
  [[ -n "${SCENARIO_JAVA_HOME:-}" ]] && args+=(--env "JAVA_CMD=$SCENARIO_JAVA_HOME/bin/java")
  [[ -n "${SCENARIO_USER_PATH:-}" ]] && args+=(--env "PATH=$SCENARIO_USER_PATH")
  [[ -n "${NXF_DISABLE_JAVA_VERSION_CHECK:-}" ]] && args+=(--env "NXF_DISABLE_JAVA_VERSION_CHECK=$NXF_DISABLE_JAVA_VERSION_CHECK")
  [[ -n "${NXF_IGNORE_JAVA_VERSION:-}" ]] && args+=(--env "NXF_IGNORE_JAVA_VERSION=$NXF_IGNORE_JAVA_VERSION")
  [[ -n "${NXF_OPTS:-}" ]] && args+=(--env "NXF_OPTS=$NXF_OPTS")

  echo "Starting SyftBox devstack via syftbox/cmd/devstack..."
  (cd "$SYFTBOX_DIR" && GOCACHE="$GO_CACHE_DIR" go run ./cmd/devstack start "${args[@]}")

  for email in "${CLIENTS[@]}"; do
    bootstrap_biovault "$email"
  done

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
