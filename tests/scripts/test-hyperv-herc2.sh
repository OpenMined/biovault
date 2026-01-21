#!/usr/bin/env bash
# Test herc2 pipeline with Podman Hyper-V backend (VM copy mode)
# This tests the SSH file transfer approach that works around Hyper-V 9P mount issues.
#
# Prerequisites:
# - Podman installed (choco install podman-cli)
# - Podman machine initialized with Hyper-V:
#   $env:CONTAINERS_MACHINE_PROVIDER = "hyperv"
#   podman machine init
#   podman machine start
#
# Usage: ./win.ps1 tests/scripts/test-hyperv-herc2.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
CLI_DIR="$ROOT_DIR/cli"
BIOSCRIPT_DIR="$(cd "$ROOT_DIR/../bioscript" && pwd)"

echo "=== Hyper-V Podman Test ==="
echo "Testing SSH file copy mode for Hyper-V backend"
echo ""

# Check Podman is available
if ! command -v podman &> /dev/null; then
    echo "Error: podman not found. Install with: choco install podman-cli"
    exit 1
fi

# Check Podman machine is running
if ! podman machine list --format "{{.Running}}" | grep -q "true"; then
    echo "Error: Podman machine not running."
    echo "Start with: podman machine start"
    exit 1
fi

# Check if Hyper-V backend
PROVIDER=$(podman machine inspect --format "{{.VMType}}" 2>/dev/null || echo "unknown")
echo "Podman VM type: $PROVIDER"

# Force Hyper-V mode via environment variable
export CONTAINERS_MACHINE_PROVIDER="hyperv"
export BIOVAULT_CONTAINER_RUNTIME="podman"

echo "CONTAINERS_MACHINE_PROVIDER=$CONTAINERS_MACHINE_PROVIDER"
echo "BIOVAULT_CONTAINER_RUNTIME=$BIOVAULT_CONTAINER_RUNTIME"
echo ""

# Test data
GENOTYPE_DIR="$CLI_DIR/tests/data/genotype_files"
HERC2_PROJECT="$BIOSCRIPT_DIR/examples/herc2/herc2-classifier"

# Create temp directory for test
TEST_DIR="/tmp/hyperv-herc2-test-$$"
mkdir -p "$TEST_DIR"
echo "Test directory: $TEST_DIR"

# Set up minimal vault environment (needed by bv CLI)
export SYC_VAULT="$TEST_DIR/.syc"
export SYFTBOX_DATA_DIR="$TEST_DIR"
mkdir -p "$SYC_VAULT"

# Convert Git Bash path (/c/Users/...) to Windows path (C:/Users/...)
windows_path() {
    local p="$1"
    echo "$p" | sed 's|^/\([a-zA-Z]\)/|\1:/|'
}

GENOTYPE_DIR_WIN=$(windows_path "$GENOTYPE_DIR")

# Create minimal samplesheet (just 1 row for fast testing)
SAMPLESHEET="$TEST_DIR/samplesheet.csv"
echo "Using 1 participant for quick test..."
cat > "$SAMPLESHEET" <<EOF
participant_id,genotype_file
01,$GENOTYPE_DIR_WIN/01_23andme_grch36_standard_4col.txt
EOF

echo "Samplesheet contents:"
cat "$SAMPLESHEET"
echo ""

RESULTS_DIR="$TEST_DIR/results"
mkdir -p "$RESULTS_DIR"

# Build CLI if needed
BV_BIN="$CLI_DIR/target/release/bv"
if [[ ! -x "$BV_BIN" ]]; then
    echo "Building BioVault CLI..."
    (cd "$CLI_DIR" && cargo build --release)
fi

# Clean old nextflow dirs from project
if [[ -d "$HERC2_PROJECT/.nextflow" ]]; then
    echo "Cleaning old .nextflow directory..."
    rm -rf "$HERC2_PROJECT/.nextflow" 2>/dev/null || true
fi
if [[ -d "$HERC2_PROJECT/work" ]]; then
    echo "Cleaning old work directory..."
    rm -rf "$HERC2_PROJECT/work" 2>/dev/null || true
fi

echo ""
echo "=== Running herc2 pipeline with Hyper-V mode ==="
echo "Project:     $HERC2_PROJECT"
echo "Samplesheet: $SAMPLESHEET"
echo "Results:     $RESULTS_DIR"
echo ""

# Run bv - should detect Hyper-V and use VM copy mode
"$BV_BIN" run "$HERC2_PROJECT" \
    --results-dir "$RESULTS_DIR" \
    --set "inputs.participants=$SAMPLESHEET"

echo ""
echo "=== Test Complete ==="
echo "Results: $RESULTS_DIR"
ls -la "$RESULTS_DIR/" 2>/dev/null || echo "(results dir empty)"

# Cleanup option
echo ""
echo "To clean up: rm -rf $TEST_DIR"
