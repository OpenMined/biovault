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
  --docker             Force Docker mode for syqure runtime
  -h, --help           Show this message

Examples:
  ./test-scenario.sh tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode go tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --sandbox sandbox-rs tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --client-mode embedded tests/scenarios/inbox-ping-pong.yaml
  ./test-scenario.sh --docker tests/scenarios/syqure-distributed.yaml
EOF
}

CLIENT_MODE="rust"
SANDBOX_DIR=""
RUST_CLIENT_BIN=""
SKIP_RUST_BUILD=0
USE_DOCKER=0
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
    --docker)
      USE_DOCKER=1
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

# Syqure runtime setup (only when scenario needs it)
if [[ "$SCENARIO" == *"syqure"* ]]; then
  if [[ -n "${BV_SYQURE_DIR:-}" ]]; then
    SYQURE_DIR="$BV_SYQURE_DIR"
  elif [[ -d "$ROOT_DIR/../syqure" ]]; then
    SYQURE_DIR="$ROOT_DIR/../syqure"
  else
    SYQURE_DIR="$ROOT_DIR/syqure"
  fi
  SYQURE_BIN="$SYQURE_DIR/target/debug/syqure"

  if (( ! USE_DOCKER )); then
    # Preflight: require a bundle or bundled codon assets for native syqure.
    BUNDLE_OK=0
    if [[ -n "${SYQURE_BUNDLE_FILE:-}" && -f "${SYQURE_BUNDLE_FILE}" ]]; then
      BUNDLE_OK=1
    else
      if command -v rustc >/dev/null 2>&1; then
        HOST_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
        if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
          BUNDLE_OK=1
        fi
      fi
      if [[ -d "$SYQURE_DIR/bin/macos-arm64/codon" || -d "$SYQURE_DIR/bin/macos-x86_64/codon" || -d "$SYQURE_DIR/bin/linux-x86/codon" || -d "$SYQURE_DIR/bin/linux-arm64/codon" ]]; then
        BUNDLE_OK=1
      fi
      if [[ -d "$ROOT_DIR/../codon/install/lib/codon" ]]; then
        BUNDLE_OK=1
      fi
    fi

    if (( ! BUNDLE_OK )); then
      if [[ -d "$SYQURE_DIR/.git" && -f "$SYQURE_DIR/.gitmodules" ]]; then
        echo "Syqure bundle/assets missing; initializing submodules..."
        (cd "$SYQURE_DIR" && git submodule update --init --recursive)
      fi

      if [[ -x "$SYQURE_DIR/build_libs.sh" ]]; then
        echo "Building syqure bundle (build_libs.sh)..."
        (cd "$SYQURE_DIR" && ./build_libs.sh)
      fi

      # Re-check after attempted build.
      if [[ -n "${SYQURE_BUNDLE_FILE:-}" && -f "${SYQURE_BUNDLE_FILE}" ]]; then
        BUNDLE_OK=1
      else
        if command -v rustc >/dev/null 2>&1; then
          HOST_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
          if [[ -n "$HOST_TRIPLE" && -f "$SYQURE_DIR/bundles/${HOST_TRIPLE}.tar.zst" ]]; then
            BUNDLE_OK=1
          fi
        fi
      fi

      if (( ! BUNDLE_OK )); then
        echo "Syqure bundle/assets not found. Provide SYQURE_BUNDLE_FILE or build syqure bundles before running." >&2
        echo "Use --docker to force Docker mode if desired." >&2
        exit 1
      fi
    fi
  fi

  if (( USE_DOCKER )); then
    export BV_SYQURE_USE_DOCKER=1
    echo "Syqure mode: Docker"
  else
    # Native mode - build syqure if needed
    if [[ ! -x "$SYQURE_BIN" ]]; then
      if [[ -d "$SYQURE_DIR" ]]; then
        echo "Building syqure native binary..."
        (cd "$SYQURE_DIR" && cargo build) || {
          echo "Failed to build syqure. Use --docker flag for Docker mode." >&2
          exit 1
        }
      else
        echo "Syqure directory not found at $SYQURE_DIR. Use --docker flag for Docker mode." >&2
        exit 1
      fi
    fi
    if [[ -x "$SYQURE_BIN" ]]; then
      export SEQURE_NATIVE_BIN="$SYQURE_BIN"
      echo "Syqure mode: Native ($SYQURE_BIN)"
    else
      echo "Syqure binary not found at $SYQURE_BIN. Use --docker flag for Docker mode." >&2
      exit 1
    fi

    # Fast sanity check: ensure bundled Codon stdlib/plugins are present.
    if [[ "${BV_SYQURE_SKIP_PREFLIGHT:-0}" != "1" ]]; then
      echo "Syqure preflight check..."
      SYQURE_INFO_OUT="$("$SYQURE_BIN" info 2>/dev/null || true)"
      if [[ -z "$SYQURE_INFO_OUT" ]]; then
        echo "Syqure info failed; cannot verify bundle contents." >&2
        exit 1
      fi
      CODON_PATH="$(printf '%s\n' "$SYQURE_INFO_OUT" | sed -n 's/^  Codon path:[[:space:]]*//p' | head -n 1)"
      if [[ -z "$CODON_PATH" || ! -d "$CODON_PATH" ]]; then
        echo "Syqure bundle check failed: Codon path not found in info output." >&2
        exit 1
      fi
      if [[ ! -d "$CODON_PATH/stdlib" || -z "$(ls -A "$CODON_PATH/stdlib" 2>/dev/null)" ]]; then
        echo "Syqure bundle check failed: stdlib missing at $CODON_PATH/stdlib" >&2
        exit 1
      fi
      if [[ ! -d "$CODON_PATH/plugins/sequre" ]]; then
        echo "Syqure bundle check failed: sequre plugin missing at $CODON_PATH/plugins/sequre" >&2
        exit 1
      fi
      echo "Syqure preflight OK (stdlib + sequre plugin present)."
    fi
  fi
fi

if python3 -c 'import yaml' >/dev/null 2>&1; then
  python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
  exit 0
fi

python3 "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
exit 0
