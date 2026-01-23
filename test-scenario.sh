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
  cat <<'EOF' >&2
Usage: ./test-scenario.sh [options] <scenario.yaml>

Options:
  --client-mode MODE   SyftBox client: go|rust|mixed|embedded (default: rust)
  --sandbox DIR        Sandbox root (default: ./sandbox)
  --rust-client-bin P  Path to Rust client binary (optional)
  --skip-rust-build    Do not build Rust client (requires binary exists)
  -h, --help           Show this message

Examples:
  ./test-scenario.sh tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode go tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --sandbox sandbox-rs tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode embedded tests/scenarios/inbox-ping-pong.yaml
EOF
}

CLIENT_MODE="rust"
SANDBOX_DIR=""
RUST_CLIENT_BIN=""
SKIP_RUST_BUILD=0
SCENARIO=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --client-mode)
      [[ $# -lt 2 ]] && { echo "Missing value for --client-mode" >&2; usage; exit 1; }
      CLIENT_MODE="${2:-}"
      shift
      ;;
    --sandbox)
      [[ $# -lt 2 ]] && { echo "Missing value for --sandbox" >&2; usage; exit 1; }
      SANDBOX_DIR="${2:-}"
      shift
      ;;
    --rust-client-bin)
      [[ $# -lt 2 ]] && { echo "Missing value for --rust-client-bin" >&2; usage; exit 1; }
      RUST_CLIENT_BIN="${2:-}"
      shift
      ;;
    --skip-rust-build)
      SKIP_RUST_BUILD=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      SCENARIO="$1"
      ;;
  esac
  shift
done

if [[ -z "$SCENARIO" ]]; then
  usage
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building BioVault CLI (release)..."
(cd "$ROOT_DIR/cli" && cargo build --release)

echo "Running scenario: $SCENARIO"

# Preserve user JAVA_HOME/PATH for downstream shells (after sbenv activation)
export SCENARIO_JAVA_HOME="${JAVA_HOME:-}"
export SCENARIO_USER_PATH="$PATH"
if command -v nextflow >/dev/null 2>&1; then
  export BIOVAULT_BUNDLED_NEXTFLOW="$(command -v nextflow)"
fi

if [[ -n "$SANDBOX_DIR" ]]; then
  export SANDBOX_DIR="$SANDBOX_DIR"
fi

# Let tests/scripts/devstack.sh decide which syftbox client to run.
export BV_DEVSTACK_CLIENT_MODE="$CLIENT_MODE"
if [[ -n "$RUST_CLIENT_BIN" ]]; then
  export BV_DEVSTACK_RUST_CLIENT_BIN="$RUST_CLIENT_BIN"
fi
if (( SKIP_RUST_BUILD )); then
  export BV_DEVSTACK_SKIP_RUST_BUILD=1
fi

if [[ "${CLIENT_MODE}" == "embedded" ]]; then
  # Ensure BioVault-hosted SyftBox runs in embedded mode when started by devstack.sh.
  export BV_SYFTBOX_BACKEND=embedded
fi

if python3 -c 'import yaml' >/dev/null 2>&1; then
  python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
  exit 0
fi

python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
exit 0
