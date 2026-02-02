#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIOVAULT_DIR="${BIOVAULT_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

OUTPUT_DIR=""
SAMPLESHEET=""
FILE_COUNT="${FILE_COUNT:-100}"
SEED="${SEED:-42}"
FORCE_REGEN="${FORCE_REGEN:-0}"
GENOSTATS_DB="${GENOSTATS_DB:-}"

APOL1_FREQ=""
APOL1_HET_FREQ=""
APOL1_HOM_ALT_FREQ=""
THAL_FREQ=""
THAL_HET_FREQ=""
THAL_HOM_ALT_FREQ=""
NO_CALL_FREQ=""
NO_CALL_TOKEN=""

info() { printf "\033[1;36m[allele-freq-gen]\033[0m %s\n" "$1"; }
error() { printf "\033[1;31m[ERROR]\033[0m %s\n" "$1" >&2; }

usage() {
    cat <<EOF_USAGE
Usage: gen_allele_freq_data.sh --output-dir DIR [OPTIONS]

Options:
  --output-dir DIR   Output directory for generated data
  --samplesheet CSV  Samplesheet path (default: <output-dir>/samplesheet.csv)
  --count N          Number of synthetic files to generate (default: 100)
  --seed N           Random seed (default: 42)
  --force            Force regeneration
  --db PATH          Path to genostats.sqlite
  --apol1 FREQ       APOL1 overlay alt frequency (0.0-1.0)
  --apol1-het FREQ   APOL1 heterozygous frequency override (0.0-1.0)
  --apol1-hom-alt FREQ APOL1 homozygous alt frequency override (0.0-1.0)
  --thal FREQ        Thalassemia overlay alt frequency (0.0-1.0)
  --thal-het FREQ    Thalassemia heterozygous frequency override (0.0-1.0)
  --thal-hom-alt FREQ Thalassemia homozygous alt frequency override (0.0-1.0)
  --no-call-frequency FREQ   No-call frequency override (0.0-1.0)
  --no-call-token TOKEN      No-call token override (default: biosynth default)
EOF_USAGE
}

build_overlay_json() {
    local json='{'
    local first=1
    local apol1_mode="freq"
    local thal_mode="freq"

    if [[ -n "$APOL1_HET_FREQ" || -n "$APOL1_HOM_ALT_FREQ" ]]; then
        apol1_mode="genotype"
    fi
    if [[ -n "$THAL_HET_FREQ" || -n "$THAL_HOM_ALT_FREQ" ]]; then
        thal_mode="genotype"
    fi

    if [[ -n "$APOL1_FREQ" ]]; then
        apol1_mode="freq"
        [[ $first -eq 0 ]] && json+=','
        json+='"APOL1":{"alt_frequency":'$APOL1_FREQ',"variants":['
        json+='{"rsid":"rs60910145","chromosome":"22","position":36265988,"reference":"T","alternates":["C","G"]},'
        json+='{"rsid":"rs71785313","chromosome":"22","position":36266000,"reference":"I","alternates":["D"]},'
        json+='{"rsid":"rs73885319","chromosome":"22","position":36265860,"reference":"A","alternates":["G"]}'
        json+=']}'
        first=0
    fi

    if [[ -n "$THAL_FREQ" ]]; then
        thal_mode="freq"
        [[ $first -eq 0 ]] && json+=','
        json+='"THALASSEMIA":{"alt_frequency":'$THAL_FREQ',"variants":['
        json+='{"rsid":"rs33985472","chromosome":"11","position":5225485,"reference":"T","alternates":["C"]},'
        json+='{"rsid":"rs33971634","chromosome":"11","position":5225660,"reference":"G","alternates":["A","C"]}'
        json+=']}'
        first=0
    fi

    if [[ "$apol1_mode" == "genotype" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"APOL1":{"variants":['
        json+='{"rsid":"rs60910145","chromosome":"22","position":36265988,"reference":"T","alternates":["C","G"],'
        json+='"genotype_frequencies":{"het":'${APOL1_HET_FREQ:-0}',"hom_alt":'${APOL1_HOM_ALT_FREQ:-0}'}}'
        json+=','
        json+='{"rsid":"rs71785313","chromosome":"22","position":36266000,"reference":"I","alternates":["D"],'
        json+='"genotype_frequencies":{"het":'${APOL1_HET_FREQ:-0}',"hom_alt":'${APOL1_HOM_ALT_FREQ:-0}'}}'
        json+=','
        json+='{"rsid":"rs73885319","chromosome":"22","position":36265860,"reference":"A","alternates":["G"],'
        json+='"genotype_frequencies":{"het":'${APOL1_HET_FREQ:-0}',"hom_alt":'${APOL1_HOM_ALT_FREQ:-0}'}}'
        json+=']}'
        first=0
    fi

    if [[ "$thal_mode" == "genotype" ]]; then
        [[ $first -eq 0 ]] && json+=','
        json+='"THALASSEMIA":{"variants":['
        json+='{"rsid":"rs33985472","chromosome":"11","position":5225485,"reference":"T","alternates":["C"],'
        json+='"genotype_frequencies":{"het":'${THAL_HET_FREQ:-0}',"hom_alt":'${THAL_HOM_ALT_FREQ:-0}'}}'
        json+=','
        json+='{"rsid":"rs33971634","chromosome":"11","position":5225660,"reference":"G","alternates":["A","C"],'
        json+='"genotype_frequencies":{"het":'${THAL_HET_FREQ:-0}',"hom_alt":'${THAL_HOM_ALT_FREQ:-0}'}}'
        json+=']}'
        first=0
    fi

    json+='}'
    echo "$json"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --count) FILE_COUNT="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --force) FORCE_REGEN=1; shift ;;
        --db) GENOSTATS_DB="$2"; shift 2 ;;
        --apol1) APOL1_FREQ="$2"; shift 2 ;;
        --apol1-het) APOL1_HET_FREQ="$2"; shift 2 ;;
        --apol1-hom-alt) APOL1_HOM_ALT_FREQ="$2"; shift 2 ;;
        --thal) THAL_FREQ="$2"; shift 2 ;;
        --thal-het) THAL_HET_FREQ="$2"; shift 2 ;;
        --thal-hom-alt) THAL_HOM_ALT_FREQ="$2"; shift 2 ;;
        --no-call-frequency) NO_CALL_FREQ="$2"; shift 2 ;;
        --no-call-token) NO_CALL_TOKEN="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    error "--output-dir is required"
    usage
    exit 1
fi

if [[ -z "$SAMPLESHEET" ]]; then
    SAMPLESHEET="${OUTPUT_DIR}/samplesheet.csv"
fi

DATA_DIR="${OUTPUT_DIR}/genotypes"
mkdir -p "$DATA_DIR"

if ! command -v bvs &>/dev/null; then
    info "bvs not found, attempting to install biosynth..."
    cargo install biosynth --locked || {
        error "Failed to install biosynth. Try: cargo install biosynth"
        exit 1
    }
fi

if [[ -z "$GENOSTATS_DB" ]]; then
    SEARCH_PATHS=(
        "$SCRIPT_DIR/../scenarios/allele-freq/data/genostats.sqlite"
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

if [[ -n "$GENOSTATS_DB" && ! -f "$GENOSTATS_DB" ]]; then
    error "genostats.sqlite not found at: $GENOSTATS_DB"
    exit 1
fi

if [[ -n "$GENOSTATS_DB" ]]; then
    info "Using database: $GENOSTATS_DB"
else
    info "No genostats.sqlite found; bvs will auto-download to data/genostats.sqlite"
fi

EXISTING_COUNT=0
if [[ -d "$DATA_DIR" ]]; then
    EXISTING_COUNT=$(find "$DATA_DIR" -name "*.txt" 2>/dev/null | wc -l | tr -d ' ')
fi

if [[ "$FORCE_REGEN" == "1" ]] || [[ "$EXISTING_COUNT" -lt "$FILE_COUNT" ]]; then
    info "Generating $FILE_COUNT synthetic genotype files (seed=$SEED)..."
    rm -rf "$DATA_DIR"/*

    SYNTH_ARGS=(
        --output "$DATA_DIR/{id}.txt"
        --count "$FILE_COUNT"
        --threads 4
        --seed "$SEED"
    )

    if [[ -n "$GENOSTATS_DB" ]]; then
        SYNTH_ARGS+=(--sqlite "$GENOSTATS_DB")
    fi

    if [[ -n "$APOL1_FREQ" || -n "$THAL_FREQ" || -n "$APOL1_HET_FREQ" || -n "$APOL1_HOM_ALT_FREQ" || -n "$THAL_HET_FREQ" || -n "$THAL_HOM_ALT_FREQ" ]]; then
        OVERLAY_JSON=$(build_overlay_json)
        info "Including overlay variants"
        SYNTH_ARGS+=(--variants-json "$OVERLAY_JSON")
    fi
    if [[ -n "$NO_CALL_FREQ" ]]; then
        SYNTH_ARGS+=(--no-call-frequency "$NO_CALL_FREQ")
    fi
    if [[ -n "$NO_CALL_TOKEN" ]]; then
        SYNTH_ARGS+=(--no-call-token "$NO_CALL_TOKEN")
    fi

    bvs synthetic "${SYNTH_ARGS[@]}"

    GENERATED=$(find "$DATA_DIR" -name "*.txt" | wc -l | tr -d ' ')
    info "Generated $GENERATED synthetic genotype files"
else
    info "Using existing synthetic data ($EXISTING_COUNT files)"
fi

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
