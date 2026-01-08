#!/usr/bin/env bash
# lint.sh - Auto-fix + run quick tests (parallel, quiet on success)
# Usage: ./lint.sh [--check] [--test]
#   --check  Read-only mode for CI (no auto-fix)
#   --test   Also run tests (slower)
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR/cli"

RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m'

CHECK_MODE=0
RUN_TESTS=0
for arg in "$@"; do
  case "$arg" in
    --check) CHECK_MODE=1 ;;
    --test) RUN_TESTS=1 ;;
  esac
done

TMPDIR_LINT=$(mktemp -d)
trap "rm -rf $TMPDIR_LINT" EXIT

FAILED=0
PIDS=()
TASKS=()

run_task() {
  local name="$1"
  local outfile="$TMPDIR_LINT/$name.out"
  shift
  TASKS+=("$name")
  echo -e "${CYAN}→ $name${NC}"
  (
    if "$@" > "$outfile" 2>&1; then
      echo "0" > "$outfile.status"
    else
      echo "1" > "$outfile.status"
    fi
  ) &
  PIDS+=($!)
}

wait_all() {
  local i=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || true
    local name="${TASKS[$i]}"
    local outfile="$TMPDIR_LINT/$name.out"
    if [[ -f "$outfile.status" && "$(cat "$outfile.status")" != "0" ]]; then
      echo -e "${RED}✗ $name${NC}"
      cat "$outfile"
      echo ""
      FAILED=1
    fi
    i=$((i + 1))
  done
}

if [[ "$CHECK_MODE" -eq 1 ]]; then
  run_task "rust-fmt" cargo fmt --all -- --check
  run_task "rust-clippy" cargo clippy --all-targets --all-features --no-deps -- -D warnings
else
  run_task "rust-fmt" cargo fmt --all
  run_task "rust-clippy" cargo clippy --fix --allow-dirty --allow-staged --all-targets --all-features --no-deps -- -D warnings
fi

if [[ "$RUN_TESTS" -eq 1 ]]; then
  run_task "rust-test" cargo test
fi

wait_all

if [[ "$FAILED" -eq 0 ]]; then
  echo -e "${GREEN}✓ All checks passed${NC}"
else
  exit 1
fi
