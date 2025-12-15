#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROFILE="${CHAOS_BUILD_PROFILE:-release}"
if [[ "$PROFILE" == "release" ]]; then
  echo "Building BioVault CLI (release)..."
  (cd "$ROOT_DIR/cli" && cargo build --release)
  export BV_BIN="${BV_BIN:-$ROOT_DIR/cli/target/release/bv}"
else
  echo "Building BioVault CLI (dev)..."
  (cd "$ROOT_DIR/cli" && cargo build)
  export BV_BIN="${BV_BIN:-$ROOT_DIR/cli/target/debug/bv}"
fi

export SCENARIO_JAVA_HOME="${JAVA_HOME:-}"
export SCENARIO_USER_PATH="$PATH"

python3 "$ROOT_DIR/tests/scripts/run_chaos.py" "$@"
