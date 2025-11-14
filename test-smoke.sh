#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$SCRIPT_DIR"

SHOW_HELP=0
GENOTYPE_COUNT="${BIOVAULT_SMOKE_GENOTYPE_COUNT:-100}"
CLEAN_DATA=1
SKIP_BUILD=0

usage() {
  cat <<'USAGE'
Usage: ./test-smoke.sh [options]

Options:
  --genotype-count N  Generate N synthetic genotype files (default: env/BIOVAULT_SMOKE_GENOTYPE_COUNT or 100)
  --skip-build      Use existing bv binary (skip cargo build)
  --no-clean        Keep existing smoke client directory
  -h, --help        Show this help message

Environment overrides:
  BIOVAULT_SMOKE_CLIENT_EMAIL  Email for the smoke-test client (default: smoke-test@syftbox.net)
  BIOVAULT_SMOKE_GENOTYPE_COUNT Number of genotype files to synthesize (default: 100)
  TEST_CLIENTS_DIR             Root directory for client data (default: ./test-clients-local)
USAGE
}

while (($#)); do
  case "$1" in
    --genotype-count)
      shift
      GENOTYPE_COUNT="${1:-}"
      if [[ -z "$GENOTYPE_COUNT" || ! "$GENOTYPE_COUNT" =~ ^[0-9]+$ || "$GENOTYPE_COUNT" -lt 1 ]]; then
        echo "Error: --genotype-count requires a positive integer" >&2
        exit 1
      fi
      ;;
    --skip-build)
      SKIP_BUILD=1
      ;;
    --no-clean)
      CLEAN_DATA=0
      ;;
    -h|--help)
      SHOW_HELP=1
      ;;
    *)
      echo "Unknown option: $1" >&2
      SHOW_HELP=1
      ;;
  esac
  shift || true
done

if [[ "$SHOW_HELP" == "1" ]]; then
  usage
  if [[ "$#" -gt 0 ]]; then
    exit 1
  fi
  exit 0
fi

CLIENT_EMAIL="${BIOVAULT_SMOKE_CLIENT_EMAIL:-smoke-test@syftbox.net}"
TEST_CLIENTS_DIR="${TEST_CLIENTS_DIR:-$ROOT_DIR/test-clients-local}"

mkdir -p "$TEST_CLIENTS_DIR"
pushd "$TEST_CLIENTS_DIR" >/dev/null
TEST_CLIENTS_DIR_ABS="$(pwd)"
popd >/dev/null

CLIENT_PATH="$TEST_CLIENTS_DIR_ABS/$CLIENT_EMAIL"

echo "🧪 Smoke test client: $CLIENT_EMAIL"
echo "📁 Client data root: $TEST_CLIENTS_DIR_ABS"
echo "🧬 Genotype files to generate: $GENOTYPE_COUNT"

if [[ "$CLEAN_DATA" == "1" ]]; then
  echo "🧹 Cleaning previous smoke client directory..."
  rm -rf "$CLIENT_PATH"
fi

if [[ "$SKIP_BUILD" == "0" ]]; then
  echo "🔧 Building BioVault CLI (release)..."
  (cd "$ROOT_DIR/cli" && cargo build --release)
else
  echo "⏭  Skipping cargo build (requested)"
fi

echo "🚀 Running BRCA smoke test..."
(
  cd "$ROOT_DIR/cli"
  TEST_MODE="local" \
  BIOVAULT_SMOKE_CLIENT_EMAIL="$CLIENT_EMAIL" \
  TEST_CLIENTS_DIR="$TEST_CLIENTS_DIR_ABS" \
  BIOVAULT_SMOKE_GENOTYPE_COUNT="$GENOTYPE_COUNT" \
  cargo test --features e2e-tests test_single_client_brca_pipeline -- --ignored --nocapture
)

echo "✅ Smoke test completed successfully"
