#!/usr/bin/env bash
set -euo pipefail

# Simple test runner with fast/slow suites
# Usage:
#   ./test.sh            # fast (default)
#   ./test.sh --fast     # fast only
#   ./test.sh --slow     # slow only (features = slow-tests)
#   ./test.sh --all      # both fast then slow

MODE=${1:---fast}

cd cli

echo "==> cargo fmt"
cargo fmt

echo "==> cargo clippy"
cargo clippy --all-targets --all-features -q || true

run_fast() {
  echo "==> Running fast tests"
  cargo test
}

run_slow() {
  echo "==> Running slow tests (feature: slow-tests, ignored only)"
  BIOVAULT_NONINTERACTIVE=1 cargo test --features slow-tests -- --ignored
}

case "$MODE" in
  --fast)
    run_fast
    ;;
  --slow)
    run_slow
    ;;
  --all)
    run_fast
    run_slow
    ;;
  *)
    echo "Unknown option: $MODE" >&2
    echo "Usage: $0 [--fast|--slow|--all]" >&2
    exit 2
    ;;
esac
