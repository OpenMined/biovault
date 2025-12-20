#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./datasite.sh [--count N] [--names email1,email2,...]

Options:
  -n, --count N          Generate N synthetic client emails (client{N}@sandbox.local).
  --names list           Comma-separated list of client emails to bootstrap.
  --reset                Stop all sandbox clients and delete the sandbox directory.
  -h, --help             Show this message.

Environment:
  SYFTBOX_SERVER_URL     Override the server URL (default: http://localhost:8080).
  SANDBOX_DIR            Override the sandbox directory (default: ./sandbox).
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SANDBOX_DIR="${SANDBOX_DIR:-$ROOT_DIR/sandbox}"
SYFTBOX_DIR="$ROOT_DIR/syftbox-sdk/syftbox"
CLI_DIR="$ROOT_DIR/cli"
SBENV_BIN="$ROOT_DIR/sbenv/sbenv"
SERVER_URL="${SYFTBOX_SERVER_URL:-http://localhost:8080}"
DEFAULT_DOMAIN="${SYFTBOX_TEST_DOMAIN:-sandbox.local}"

COUNT_REQUEST=0
declare -a RAW_NAME_SETS=()
RESET_FLAG=0
PROVISION_FLAG=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--count)
      [[ $# -lt 2 ]] && { echo "Missing value for $1" >&2; usage >&2; exit 1; }
      COUNT_REQUEST="$2"
      [[ ! "$COUNT_REQUEST" =~ ^[0-9]+$ ]] && { echo "Count must be numeric" >&2; exit 1; }
      shift
      PROVISION_FLAG=1
      ;;
    --names)
      [[ $# -lt 2 ]] && { echo "Missing value for --names" >&2; usage >&2; exit 1; }
      RAW_NAME_SETS+=("$2")
      shift
      PROVISION_FLAG=1
      ;;
    --reset)
      RESET_FLAG=1
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

trim() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s' "$value"
}

declare -a EMAILS=()

email_exists() {
  local needle="$1"
  local element
  for element in "${EMAILS[@]:-}"; do
    if [[ "$element" == "$needle" ]]; then
      return 0
    fi
  done
  return 1
}

add_email() {
  local raw norm
  raw="$(trim "$1")"
  [[ -z "$raw" ]] && return
  norm="$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]')"
  if ! email_exists "$norm"; then
    EMAILS+=("$norm")
  fi
}

stop_client() {
  local client_dir="$1"
  [[ -d "$client_dir" ]] || return
  if [[ -x "$SBENV_BIN" ]]; then
    (
      set +u
      cd "$client_dir"
      "$SBENV_BIN" stop >/dev/null 2>&1 || true
    )
  fi
  local pid_file="$client_dir/.syftbox/syftbox.pid"
  if [[ -f "$pid_file" ]]; then
    local pid
    pid="$(tr -d ' \n\r' < "$pid_file")"
    if [[ -n "$pid" ]] && ps -p "$pid" >/dev/null 2>&1; then
      kill "$pid" >/dev/null 2>&1 || true
    fi
  fi
}

reset_sandbox() {
  if [[ ! -d "$SANDBOX_DIR" ]]; then
    echo "Sandbox directory $SANDBOX_DIR not found (nothing to reset)."
    return
  fi
  echo "Stopping sandbox clients in $SANDBOX_DIR..."
  while IFS= read -r client_path; do
    [[ -z "$client_path" ]] && continue
    stop_client "$client_path"
  done < <(find "$SANDBOX_DIR" -mindepth 1 -maxdepth 1 -type d -print 2>/dev/null)
  echo "Removing $SANDBOX_DIR..."
  rm -rf "$SANDBOX_DIR"
  echo "Sandbox reset complete."
}

if (( PROVISION_FLAG == 0 && RESET_FLAG == 0 )); then
  PROVISION_FLAG=1
fi

if (( RESET_FLAG )); then
  reset_sandbox
  if (( PROVISION_FLAG == 0 )); then
    exit 0
  fi
fi

if (( PROVISION_FLAG )); then
  if ((${#RAW_NAME_SETS[@]})); then
    for block in "${RAW_NAME_SETS[@]}"; do
      IFS=',' read -r -a pieces <<< "$block"
      for piece in "${pieces[@]}"; do
        add_email "$piece"
      done
    done
  fi

  if (( COUNT_REQUEST > 0 )); then
    for idx in $(seq 1 "$COUNT_REQUEST"); do
      add_email "client${idx}@${DEFAULT_DOMAIN}"
    done
  fi

  if ((${#EMAILS[@]} == 0)); then
    for idx in 1 2; do
      add_email "client${idx}@${DEFAULT_DOMAIN}"
    done
  fi
fi

require_file() {
  [[ -f "$1" ]] || { echo "Missing required file: $1" >&2; exit 1; }
}

require_bin() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing required tool: $1" >&2; exit 1; }
}

if (( PROVISION_FLAG )); then
  require_bin go
  require_bin cargo
  require_bin curl
  require_file "$SBENV_BIN"
  [[ -d "$SYFTBOX_DIR" ]] || { echo "Missing syftbox checkout at $SYFTBOX_DIR" >&2; exit 1; }
  [[ -d "$CLI_DIR" ]] || { echo "Missing cli checkout at $CLI_DIR" >&2; exit 1; }

  mkdir -p "$SANDBOX_DIR"

  echo "Checking SyftBox server at $SERVER_URL..."
  curl -fsS --max-time 5 "$SERVER_URL" >/dev/null || {
    echo "Server is not reachable at $SERVER_URL" >&2
    exit 1
  }
else
  exit 0
fi

build_syftbox_binary() {
  local goos goarch target cgo ldflags version commit build_date
  goos="$(go env GOOS)"
  goarch="$(go env GOARCH)"
  target="$SYFTBOX_DIR/.out/syftbox_client_${goos}_${goarch}"
  if [[ -x "$target" ]]; then
    printf '%s\n' "$target"
    return
  fi

  mkdir -p "$SYFTBOX_DIR/.out"
  version="$(cd "$SYFTBOX_DIR" && git describe --tags --always 2>/dev/null || echo "dev")"
  commit="$(cd "$SYFTBOX_DIR" && git rev-parse --short HEAD 2>/dev/null || echo "unknown")"
  build_date="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  ldflags="-s -w"
  ldflags+=" -X github.com/openmined/syftbox/internal/version.Version=$version"
  ldflags+=" -X github.com/openmined/syftbox/internal/version.Revision=$commit"
  ldflags+=" -X github.com/openmined/syftbox/internal/version.BuildDate=$build_date"

  cgo=0
  [[ "$goos" == "darwin" ]] && cgo=1

  echo "Building SyftBox client binary for ${goos}/${goarch}..."
  (cd "$SYFTBOX_DIR" && GOOS="$goos" GOARCH="$goarch" CGO_ENABLED="$cgo" \
    go build -trimpath --tags "go_json nomsgpack" -ldflags "$ldflags" \
    -o "$target" ./cmd/client)
  printf '%s\n' "$target"
}

ensure_bv_binary() {
  local target="$CLI_DIR/target/release/bv"
  if [[ ! -x "$target" ]]; then
    echo "Compiling BioVault CLI (cargo build --release)..."
    (cd "$CLI_DIR" && cargo build --release >/dev/null)
  fi
  printf '%s\n' "$target"
}

SYFTBOX_BINARY_PATH="${SYFTBOX_BINARY_PATH:-$(build_syftbox_binary)}"
BV_BINARY_PATH="$(ensure_bv_binary)"

with_client_env() {
  local email="$1"; shift
  (
    set +u
    cd "$SANDBOX_DIR/$email"
    eval "$("$SBENV_BIN" activate --quiet)"
    export SYFTBOX_BINARY="${SYFTBOX_BINARY:-$SYFTBOX_BINARY_PATH}"
    "$@"
  )
}

bootstrap_client() {
  local email="$1"
  local client_dir="$SANDBOX_DIR/$email"
  rm -rf "$client_dir"
  mkdir -p "$client_dir"
  pushd "$client_dir" >/dev/null
  local args=(init --dev --server-url "$SERVER_URL" --email "$email")
  args+=(--binary "$SYFTBOX_BINARY_PATH")
  "$SBENV_BIN" "${args[@]}" >/dev/null
  popd >/dev/null
}

declare -a LAUNCHER_PIDS=()

start_client_daemon() {
  local email="$1"
  local client_dir="$SANDBOX_DIR/$email"
  local launcher_log="$client_dir/.syftbox/launcher.log"
  echo "Starting SyftBox daemon for $email"
  (
    set +u
    cd "$client_dir"
    eval "$("$SBENV_BIN" activate --quiet)"
    # Use SYFTBOX_BINARY_PATH (absolute path) to avoid circular symlink
    ln -sf "$SYFTBOX_BINARY_PATH" ./syftbox
    export SYFTBOX_BINARY="${SYFTBOX_BINARY:-$SYFTBOX_BINARY_PATH}"
    PATH="$PWD:$PATH"
    "$SBENV_BIN" start --skip-login-check >>"$launcher_log" 2>&1
  ) &
  LAUNCHER_PIDS+=("$!")
}

wait_for_client_ready() {
  local email="$1"
  local client_dir="$SANDBOX_DIR/$email"
  local pid_file="$client_dir/.syftbox/syftbox.pid"
  local log_file="$client_dir/.syftbox/daemon.log"
  local ready=0
  for attempt in $(seq 1 90); do
    if [[ -f "$pid_file" ]]; then
      local daemon_pid
      daemon_pid="$(tr -d ' \n\r' < "$pid_file")"
      if [[ -n "$daemon_pid" ]] && ps -p "$daemon_pid" >/dev/null 2>&1; then
        if [[ -f "$log_file" ]] && grep -q "full sync completed" "$log_file"; then
          ready=1
          break
        fi
      fi
    fi
    sleep 2
  done

  if (( ! ready )); then
    echo "Client $email failed to reach steady state." >&2
    [[ -f "$log_file" ]] && tail -n 40 "$log_file" >&2
    return 1
  fi
  echo "Client $email is syncing."
}

run_biovault_init() {
  local email="$1"
  echo "Configuring BioVault for $email"
  with_client_env "$email" "$BV_BINARY_PATH" init --quiet "$email" >/dev/null
  with_client_env "$email" "$BV_BINARY_PATH" config >/dev/null
}

write_shared_probe() {
  local owner="$1"
  local owner_dir="$SANDBOX_DIR/$owner/datasites/$owner/shared"
  ensure_shared_acl "$owner"
  mkdir -p "$owner_dir"
  local stamp
  stamp="$(date -u +%Y%m%dT%H%M%SZ)"
  local fname="sandbox-sync-${stamp}-${RANDOM}.txt"
  local fpath="$owner_dir/$fname"
  printf 'owner=%s\ncreated=%s\n' "$owner" "$stamp" >"$fpath"
  printf '%s\n' "$fpath"
}

ensure_shared_acl() {
  local owner="$1"
  local shared_dir="$SANDBOX_DIR/$owner/datasites/$owner/shared"
  mkdir -p "$shared_dir"
  local perm_file="$shared_dir/syft.pub.yaml"
  cat >"$perm_file" <<'EOF'
rules:
  - pattern: '**'
    access:
      admin: []
      read:
        - '*'
      write: []
EOF
}

wait_for_replica() {
  local path="$1"
  for attempt in $(seq 1 90); do
    [[ -f "$path" ]] && return 0
    sleep 1
  done
  return 1
}

echo "Provisioning sandbox clients in $SANDBOX_DIR:"
for email in "${EMAILS[@]}"; do
  echo "  - $email"
  bootstrap_client "$email"
  start_client_daemon "$email"
done

for email in "${EMAILS[@]}"; do
  wait_for_client_ready "$email"
  run_biovault_init "$email"
done

if ((${#EMAILS[@]} > 1)); then
  owner="${EMAILS[0]}"
  probe_file="$(write_shared_probe "$owner")"
  echo "Wrote sync probe at $probe_file"
  for email in "${EMAILS[@]:1}"; do
    replica="$SANDBOX_DIR/$email/datasites/$owner/shared/$(basename "$probe_file")"
    echo "Waiting for replication to $email..."
    if wait_for_replica "$replica"; then
      cmp -s "$probe_file" "$replica" || {
        echo "Content mismatch for $replica" >&2
        exit 1
      }
      echo "  ✓ $email received the probe."
    else
      echo "  ✗ $email did not receive $replica in time." >&2
      exit 1
    fi
  done
  echo "Shared file replicated across all clients."
else
  echo "Single client requested; skipping replication check."
fi

echo ""
echo "Sandbox ready. Active clients:"
for idx in "${!EMAILS[@]}"; do
  email="${EMAILS[$idx]}"
  pid="${LAUNCHER_PIDS[$idx]:-unknown}"
  printf '  %-30s launcher pid=%s\n' "$email" "$pid"
done
echo "Data root: $SANDBOX_DIR"
