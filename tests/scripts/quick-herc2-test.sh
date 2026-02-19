#!/usr/bin/env bash
# Quick test for herc2 pipeline with Podman nested containers
# Runs just 1 participant to speed up testing
#
# Usage: ./tests/scripts/quick-herc2-test.sh [--all]
#   --all    Run all 13 participants (default: just 1)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
CLI_DIR="$ROOT_DIR/cli"
BIOSCRIPT_DIR="$(cd "$ROOT_DIR/../bioscript" && pwd)"

RUN_ALL=0
if [[ "${1:-}" == "--all" ]]; then
    RUN_ALL=1
fi

# Test data
GENOTYPE_DIR="$CLI_DIR/tests/data/genotype_files"
HERC2_PROJECT="$BIOSCRIPT_DIR/examples/herc2/herc2-classifier"

# Create temp directory for test
TEST_DIR="/tmp/quick-herc2-test-$$"
mkdir -p "$TEST_DIR"
echo "Test directory: $TEST_DIR"

# Set up minimal vault environment (needed by bv CLI)
export SBC_VAULT="$TEST_DIR/.sbc"
export SYFTBOX_DATA_DIR="$TEST_DIR"
mkdir -p "$SBC_VAULT"

# Convert Git Bash path (/c/Users/...) to Windows path (C:/Users/...)
# The BioVault CLI expects Windows paths and handles conversion to container format
windows_path() {
    local p="$1"
    # Convert /c/ to C:/
    echo "$p" | sed 's|^/\([a-zA-Z]\)/|\1:/|'
}

GENOTYPE_DIR_WIN=$(windows_path "$GENOTYPE_DIR")

# Create minimal samplesheet (just 1 row for fast testing)
# Use Windows-style paths - the BioVault CLI will extract these, mount them,
# and rewrite them to container-compatible format automatically
SAMPLESHEET="$TEST_DIR/samplesheet.csv"
if (( RUN_ALL )); then
    echo "Using all 13 participants..."
    cat > "$SAMPLESHEET" <<EOF
participant_id,genotype_file
01,$GENOTYPE_DIR_WIN/01_23andme_grch36_standard_4col.txt
02,$GENOTYPE_DIR_WIN/02_23andme_grch37_standard_4col.txt
03,$GENOTYPE_DIR_WIN/03_23andme_standard_4col.txt
04,$GENOTYPE_DIR_WIN/04_ancestrydna_grch37_v1_split_alleles_5col.txt
05,$GENOTYPE_DIR_WIN/05_ancestrydna_split_alleles_5col.txt
06,$GENOTYPE_DIR_WIN/06_decodeme_grch36_variation_6col.csv
07,$GENOTYPE_DIR_WIN/07_dynamicdna_grch38_extended_7col.txt
08,$GENOTYPE_DIR_WIN/08_familytreedna_grch37_standard_4col.csv
09,$GENOTYPE_DIR_WIN/09_genesforgood_standard_4col.txt
10,$GENOTYPE_DIR_WIN/10_livingdna_grch37_standard_4col.txt
11,$GENOTYPE_DIR_WIN/11_myheritage_grch37_standard_4col.csv
12,$GENOTYPE_DIR_WIN/12_unknown_no_header.txt
13,$GENOTYPE_DIR_WIN/13_unknown_standard_4col.txt
EOF
else
    echo "Using 1 participant for quick test..."
    cat > "$SAMPLESHEET" <<EOF
participant_id,genotype_file
01,$GENOTYPE_DIR_WIN/01_23andme_grch36_standard_4col.txt
EOF
fi

RESULTS_DIR="$TEST_DIR/results"
mkdir -p "$RESULTS_DIR"

# Build CLI if needed
BV_BIN="$CLI_DIR/target/release/bv"
if [[ ! -x "$BV_BIN" ]]; then
    echo "Building BioVault CLI..."
    (cd "$CLI_DIR" && cargo build --release)
fi

# Clean old nextflow dirs from project (permission issues when switching runtimes)
if [[ -d "$HERC2_PROJECT/.nextflow" ]]; then
    echo "Cleaning old .nextflow directory..."
    rm -rf "$HERC2_PROJECT/.nextflow" 2>/dev/null || true
fi
if [[ -d "$HERC2_PROJECT/work" ]]; then
    echo "Cleaning old work directory..."
    rm -rf "$HERC2_PROJECT/work" 2>/dev/null || true
fi

echo ""
echo "=== Running herc2 pipeline ==="
echo "Project:     $HERC2_PROJECT"
echo "Samplesheet: $SAMPLESHEET"
echo "Results:     $RESULTS_DIR"
echo ""

# Run bv directly (not through test harness)
# This uses Podman by default on Windows when Docker is not available
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
