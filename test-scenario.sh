#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: ./test-scenario.sh <scenario.yaml>" >&2
  echo "Example: ./test-scenario.sh tests/scenarios/inbox-ping-pong.yaml" >&2
  exit 1
fi

SCENARIO="$1"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building BioVault CLI (release)..."
(cd "$ROOT_DIR/cli" && cargo build --release)

echo "Running scenario: $SCENARIO"

# Preserve user JAVA_HOME/PATH for downstream shells (after sbenv activation)
export SCENARIO_JAVA_HOME="${JAVA_HOME:-}"
export SCENARIO_USER_PATH="$PATH"

uv run --with pyyaml "$ROOT_DIR/scripts/run_scenario.py" "$SCENARIO"
