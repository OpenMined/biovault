#!/usr/bin/env bash
set -euo pipefail

# Secure aggregation using Syqure MPC
# This runs on ALL parties (aggregator + clients) simultaneously

CODON_FILE="$(dirname "$0")/secure_aggregate.codon"

# Find syqure binary
SYQURE_BIN="${SEQURE_NATIVE_BIN:-}"
if [[ -z "$SYQURE_BIN" || ! -x "$SYQURE_BIN" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

    # Try common locations
    for candidate in \
        "${BIOVAULT_WORKSPACE:-}/../syqure/target/release/syqure" \
        "${BIOVAULT_WORKSPACE:-}/../syqure/target/debug/syqure" \
        "${SCRIPT_DIR}/../../../../../../../../../syqure/target/release/syqure" \
        "${SCRIPT_DIR}/../../../../../../../../../syqure/target/debug/syqure" \
        "/Users/madhavajay/dev/biovault-desktop/workspace2/syqure/target/release/syqure" \
        "/Users/madhavajay/dev/biovault-desktop/workspace2/syqure/target/debug/syqure"; do
        if [[ -x "$candidate" ]]; then
            SYQURE_BIN="$(cd "$(dirname "$candidate")" && pwd)/$(basename "$candidate")"
            break
        fi
    done
fi

if [[ -z "$SYQURE_BIN" || ! -x "$SYQURE_BIN" ]]; then
    echo "ERROR: syqure binary not found. Set SEQURE_NATIVE_BIN or build syqure." >&2
    exit 1
fi

# Get inputs from BioVault environment
PARTY_ID="${BV_INPUT_PARTY_ID:-}"
PARTY_EMAILS="${BV_INPUT_PARTY_EMAILS:-}"
COUNTS_PATH="${BV_INPUT_COUNTS:-}"
ARRAY_LENGTH_PATH="${BV_INPUT_ARRAY_LENGTH:-}"
OUTPUT_PATH="${BV_OUTPUT_AGGREGATED_COUNTS:-aggregated_counts.json}"

# Read array_length from file if it's a path (trim whitespace)
ARRAY_LENGTH="0"
echo "DEBUG: ARRAY_LENGTH_PATH='$ARRAY_LENGTH_PATH'"
if [[ -n "$ARRAY_LENGTH_PATH" ]]; then
    if [[ -f "$ARRAY_LENGTH_PATH" ]]; then
        ARRAY_LENGTH=$(cat "$ARRAY_LENGTH_PATH" | tr -d '[:space:]')
        echo "DEBUG: Read ARRAY_LENGTH='$ARRAY_LENGTH' from file"
    else
        echo "DEBUG: File not found at '$ARRAY_LENGTH_PATH'"
    fi
else
    echo "DEBUG: ARRAY_LENGTH_PATH is empty"
fi

# Get datasites root from environment
# SYFTBOX_DATA_DIR is the client's data folder (e.g., sandbox/client1@sandbox.local)
# Datasites are at {data_dir}/datasites/ (synced via SyftBox)
DATASITES_ROOT="${SYFTBOX_DATA_DIR:-}"
if [[ -n "$DATASITES_ROOT" ]]; then
    DATASITES_ROOT="${DATASITES_ROOT}/datasites"
fi

if [[ -z "$PARTY_ID" ]]; then
    echo "ERROR: BV_INPUT_PARTY_ID not set" >&2
    exit 1
fi

if [[ -z "$PARTY_EMAILS" ]]; then
    echo "ERROR: BV_INPUT_PARTY_EMAILS not set" >&2
    exit 1
fi

# Count parties
CP_COUNT=$(echo "$PARTY_EMAILS" | tr ',' '\n' | wc -l | tr -d ' ')

echo "=== Secure Aggregate Module ==="
echo "  Party ID: $PARTY_ID"
echo "  Party emails: $PARTY_EMAILS"
echo "  Counts path: $COUNTS_PATH"
echo "  Array length: $ARRAY_LENGTH"
echo "  Output path: $OUTPUT_PATH"
echo "  Datasites root: $DATASITES_ROOT"
echo "  Syqure binary: $SYQURE_BIN"
echo "==============================="

# Analyze codon file before running (validates syntax and shows cost estimate)
echo ""
ANALYZE_OUTPUT=$("$SYQURE_BIN" analyze "$CODON_FILE")
echo "$ANALYZE_OUTPUT"
echo ""

# Check if we can skip MHE setup (MPC only programs don't need it)
SKIP_MHE=""
if echo "$ANALYZE_OUTPUT" | grep -q "Can use --skip-mhe-setup:  Yes"; then
    SKIP_MHE="--skip-mhe-setup"
    echo "Note: Skipping MHE setup (MPC only)"
fi

# Set up Syqure environment for file transport
export SEQURE_TRANSPORT=file
export SEQURE_FILE_DIR="shared/flows/${BV_FLOW_NAME:-syqure-flow}/${BV_RUN_ID:-unknown}/_mpc"
export SEQURE_FILE_POLL_MS=50
export SEQURE_CP_COUNT="$CP_COUNT"
export SEQURE_PARTY_EMAILS="$PARTY_EMAILS"
export SEQURE_DATASITES_ROOT="$DATASITES_ROOT"
export SEQURE_FILE_KEEP=1
export SEQURE_FILE_DEBUG=1

# Pass data to codon program
export SEQURE_INPUT_COUNTS="$COUNTS_PATH"
export SEQURE_ARRAY_LENGTH="$ARRAY_LENGTH"
export SEQURE_OUTPUT_AGGREGATED="$OUTPUT_PATH"

# Run syqure with the codon file
# The party ID is passed as a command-line argument after --
exec "$SYQURE_BIN" $SKIP_MHE "$CODON_FILE" -- "$PARTY_ID"
