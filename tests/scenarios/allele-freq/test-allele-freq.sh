#!/usr/bin/env bash
set -euo pipefail

# Test script for allele-freq pipeline
# Generates synthetic genotype files and runs the pipeline

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR/allele-freq"
# Go up from tests/scenarios/allele-freq -> tests/scenarios -> tests -> biovault
BIOVAULT_DIR="${BIOVAULT_DIR:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"

# Configuration
FILE_COUNT="${FILE_COUNT:-100}"
SEED="${SEED:-42}"
OUTPUT_DIR="$SCRIPT_DIR/test-output"
DATA_DIR="$SCRIPT_DIR/test-data"
FORCE_REGEN="${FORCE_REGEN:-0}"
GENOSTATS_DB="${GENOSTATS_DB:-}"

# Per-group frequencies (empty = not enabled)
HERC2_FREQ=""
APOL1_FREQ=""
BRCA_FREQ=""
THAL_FREQ=""

info() { printf "\033[1;36m[allele-freq-test]\033[0m %s\n" "$1"; }
error() { printf "\033[1;31m[ERROR]\033[0m %s\n" "$1" >&2; }

# Build overlay JSON based on enabled groups with their frequencies
build_overlay_json() {
    local json='{'
    local first=1

    if [[ -n "$HERC2_FREQ" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"HERC2":{"alt_frequency":'$HERC2_FREQ',"variants":[
            {"rsid":"rs12913832","chromosome":"15","position":28120472,"reference":"A","alternates":["C","G"]}
        ]}'
        first=0
    fi

    if [[ -n "$APOL1_FREQ" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"APOL1":{"alt_frequency":'$APOL1_FREQ',"variants":[
            {"rsid":"rs60910145","chromosome":"22","position":36265988,"reference":"T","alternates":["C","G"]},
            {"rsid":"rs71785313","chromosome":"22","position":36266000,"reference":"I","alternates":["D"]},
            {"rsid":"rs73885319","chromosome":"22","position":36265860,"reference":"A","alternates":["G"]}
        ]}'
        first=0
    fi

    if [[ -n "$BRCA_FREQ" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"BRCA":{"alt_frequency":'$BRCA_FREQ',"variants":[
            {"rsid":"rs80357336","chromosome":"17","position":43045711,"reference":"G","alternates":["A","C","T"]},
            {"rsid":"rs80358650","chromosome":"13","position":32316463,"reference":"G","alternates":["A","C","T"]}
        ]}'
        first=0
    fi

    if [[ -n "$THAL_FREQ" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"THALASSEMIA":{"alt_frequency":'$THAL_FREQ',"variants":[
            {"rsid":"rs33985472","chromosome":"11","position":5225485,"reference":"T","alternates":["C"]},
            {"rsid":"rs33971634","chromosome":"11","position":5225660,"reference":"G","alternates":["A","C"]}
        ]}'
        first=0
    fi

    json+='}'
    echo "$json"
}

has_overlays() {
    [[ -n "$HERC2_FREQ" ]] || [[ -n "$APOL1_FREQ" ]] || [[ -n "$BRCA_FREQ" ]] || [[ -n "$THAL_FREQ" ]]
}

usage() {
    cat <<EOF
Usage: ./test-allele-freq.sh [OPTIONS]

Options:
    --count N        Number of synthetic files to generate (default: 100)
    --seed N         Random seed for synthetic data (default: 42)
    --force          Force regeneration of synthetic data
    --clean          Clean up test data and output before running
    --db PATH        Path to genostats.sqlite (auto-detected if not specified)
    -h, --help       Show this help

Overlay Variants (each with its own frequency 0.0-1.0):
    --herc2 FREQ     Eye color variants
    --apol1 FREQ     Kidney risk variants
    --brca FREQ      Cancer risk variants
    --thal FREQ      Thalassemia variants
    --all FREQ       All of the above at same FREQ

Variant Details:
    herc2:  rs12913832 (chr15:28120472)
    apol1:  rs60910145 (chr22:36265988)
            rs71785313 (chr22:36266000)
            rs73885319 (chr22:36265860)
    brca:   rs80357336 (chr17:43045711) - BRCA1
            rs80358650 (chr13:32316463) - BRCA2
    thal:   rs33985472 (chr11:5225485)
            rs33971634 (chr11:5225660)

Environment:
    BIOVAULT_DIR     Path to biovault repo (auto-detected)
    GENOSTATS_DB     Path to genostats.sqlite database

Examples:
    ./test-allele-freq.sh --count 10                      # No overlays
    ./test-allele-freq.sh --herc2 0.8                     # Eye color at 80%
    ./test-allele-freq.sh --apol1 0.5 --thal 0.5          # Two groups at 50% each
    ./test-allele-freq.sh --herc2 0.9 --brca 0.3          # Different frequencies
    ./test-allele-freq.sh --all 0.7 --count 20 --force    # All at 70%
EOF
}

CLEAN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --count) FILE_COUNT="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --force) FORCE_REGEN=1; shift ;;
        --clean) CLEAN=1; shift ;;
        --herc2) HERC2_FREQ="$2"; shift 2 ;;
        --apol1) APOL1_FREQ="$2"; shift 2 ;;
        --brca) BRCA_FREQ="$2"; shift 2 ;;
        --thal) THAL_FREQ="$2"; shift 2 ;;
        --all) HERC2_FREQ="$2"; APOL1_FREQ="$2"; BRCA_FREQ="$2"; THAL_FREQ="$2"; shift 2 ;;
        --db) GENOSTATS_DB="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# Clean if requested
if [[ "$CLEAN" == "1" ]]; then
    info "Cleaning test data and output..."
    rm -rf "$DATA_DIR" "$OUTPUT_DIR"
fi

mkdir -p "$DATA_DIR" "$OUTPUT_DIR"

# Check for bvs
if ! command -v bvs &>/dev/null; then
    info "bvs not found, attempting to install biosynth..."
    cargo install biosynth --locked || {
        error "Failed to install biosynth. Try: cargo install biosynth"
        exit 1
    }
fi

# Locate genostats.sqlite database if not specified
if [[ -z "$GENOSTATS_DB" ]]; then
    # Search in common locations
    SEARCH_PATHS=(
        "$SCRIPT_DIR/data/genostats.sqlite"
        "$BIOVAULT_DIR/data/genostats.sqlite"
        "$(dirname "$BIOVAULT_DIR")/data/genostats.sqlite"
        "$HOME/.biovault/genostats.sqlite"
        "./data/genostats.sqlite"
    )
    for path in "${SEARCH_PATHS[@]}"; do
        if [[ -f "$path" ]]; then
            GENOSTATS_DB="$path"
            break
        fi
    done
fi

if [[ -z "$GENOSTATS_DB" ]] || [[ ! -f "$GENOSTATS_DB" ]]; then
    error "genostats.sqlite not found!"
    error "Searched in:"
    for path in "${SEARCH_PATHS[@]}"; do
        error "  - $path"
    done
    error ""
    error "Please specify with --db or set GENOSTATS_DB environment variable"
    error "Or run 'bvs genostats -i <genotype_files>' to create one"
    exit 1
fi

info "Using database: $GENOSTATS_DB"

# Generate synthetic data if needed
EXISTING_COUNT=0
if [[ -d "$DATA_DIR" ]]; then
    EXISTING_COUNT=$(find "$DATA_DIR" -name "*.txt" 2>/dev/null | wc -l | tr -d ' ')
fi

if [[ "$FORCE_REGEN" == "1" ]] || [[ "$EXISTING_COUNT" -lt "$FILE_COUNT" ]]; then
    info "Generating $FILE_COUNT synthetic genotype files (seed=$SEED)..."

    rm -rf "$DATA_DIR"/*

    SYNTH_ARGS=(
        --sqlite "$GENOSTATS_DB"
        --output "$DATA_DIR/{id}.txt"
        --count "$FILE_COUNT"
        --threads 4
        --seed "$SEED"
    )

    if has_overlays; then
        OVERLAY_JSON=$(build_overlay_json)
        info "Including overlay variants:"
        [[ -n "$HERC2_FREQ" ]] && info "  herc2: $HERC2_FREQ"
        [[ -n "$APOL1_FREQ" ]] && info "  apol1: $APOL1_FREQ"
        [[ -n "$BRCA_FREQ" ]] && info "  brca: $BRCA_FREQ"
        [[ -n "$THAL_FREQ" ]] && info "  thal: $THAL_FREQ"
        SYNTH_ARGS+=(--variants-json "$OVERLAY_JSON")
    fi

    bvs synthetic "${SYNTH_ARGS[@]}"

    GENERATED=$(find "$DATA_DIR" -name "*.txt" | wc -l | tr -d ' ')
    info "Generated $GENERATED synthetic genotype files"
else
    info "Using existing synthetic data ($EXISTING_COUNT files)"
fi

# Create samplesheet CSV
SAMPLESHEET="$SCRIPT_DIR/test-samplesheet.csv"
info "Creating samplesheet: $SAMPLESHEET"

echo "participant_id,genotype_file" > "$SAMPLESHEET"
for f in "$DATA_DIR"/*.txt; do
    if [[ -f "$f" ]]; then
        basename_f=$(basename "$f" .txt)
        echo "$basename_f,$f" >> "$SAMPLESHEET"
    fi
done

SAMPLE_COUNT=$(tail -n +2 "$SAMPLESHEET" | wc -l | tr -d ' ')
info "Samplesheet has $SAMPLE_COUNT entries"

# Check for bv CLI
if ! command -v bv &>/dev/null; then
    info "Building biovault CLI..."
    (cd "$BIOVAULT_DIR/cli" && cargo build --release) || {
        error "Failed to build biovault CLI"
        exit 1
    }
    BV_CMD="$BIOVAULT_DIR/cli/target/release/bv"
else
    BV_CMD="bv"
fi

# Run the pipeline
info "Running allele-freq pipeline..."
info "  Pipeline: $PIPELINE_DIR"
info "  Files: $SAMPLE_COUNT"
info "  Output: $OUTPUT_DIR"

cd "$PIPELINE_DIR"

$BV_CMD run . \
    --participants "$SAMPLESHEET" \
    --results-dir "$OUTPUT_DIR" \
    2>&1 | tee "$OUTPUT_DIR/pipeline.log"

# Check results - bv run outputs to "results/" directory in the pipeline dir
RESULTS_DIR="$PIPELINE_DIR/results"
info "=== Results ==="
if [[ -f "$RESULTS_DIR/allele_freq.tsv" ]]; then
    LOCUS_COUNT=$(tail -n +2 "$RESULTS_DIR/allele_freq.tsv" | wc -l | tr -d ' ')
    info "Allele frequency table: $RESULTS_DIR/allele_freq.tsv"
    info "Total loci: $LOCUS_COUNT"

    info "Sample output (first 10 rows):"
    head -11 "$RESULTS_DIR/allele_freq.tsv" | column -t -s $'\t'

    # If overlay was used, show those specific rsids
    if has_overlays; then
        info ""
        info "Overlay variant frequencies:"

        # Build grep pattern based on enabled groups
        GREP_PATTERN=""
        [[ -n "$HERC2_FREQ" ]] && GREP_PATTERN="${GREP_PATTERN}|rs12913832"
        [[ -n "$APOL1_FREQ" ]] && GREP_PATTERN="${GREP_PATTERN}|rs60910145|rs71785313|rs73885319"
        [[ -n "$BRCA_FREQ" ]] && GREP_PATTERN="${GREP_PATTERN}|rs80357336|rs80358650"
        [[ -n "$THAL_FREQ" ]] && GREP_PATTERN="${GREP_PATTERN}|rs33985472|rs33971634"
        # Remove leading pipe
        GREP_PATTERN="${GREP_PATTERN#|}"

        grep -E "$GREP_PATTERN" "$RESULTS_DIR/allele_freq.tsv" | column -t -s $'\t' || true
    fi
else
    error "allele_freq.tsv not found!"
fi

if [[ -f "$RESULTS_DIR/vcf_conversion_results.tsv" ]]; then
    info ""
    info "VCF conversion stats: $RESULTS_DIR/vcf_conversion_results.tsv"
    head -6 "$RESULTS_DIR/vcf_conversion_results.tsv" | column -t -s $'\t'
fi

info ""
info "Done! Results in: $RESULTS_DIR"
