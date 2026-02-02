#!/usr/bin/env bash
set -euo pipefail

BV_BIN="${BV_BIN:-bv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SRC_DIR="$(cd "${SCRIPT_DIR}/../../../allele-freq/syqure-allele-freq/assets/allele-freq" && pwd)"

if ! command -v "$BV_BIN" >/dev/null 2>&1; then
  echo "bv binary not found: ${BV_BIN}" >&2
  exit 1
fi

DATA_ROOT="${ALLELE_FREQ_DATA_DIR:-}"
if [[ -z "$DATA_ROOT" ]]; then
  if [[ -n "${BV_SYFTBOX_DATA_DIR:-}" ]]; then
    DATA_ROOT="${BV_SYFTBOX_DATA_DIR}/private/app_data/biovault/allele-freq-data"
  elif [[ -n "${BV_DATASITES_ROOT:-}" ]]; then
    DATA_ROOT="$(dirname "${BV_DATASITES_ROOT}")/private/app_data/biovault/allele-freq-data"
  fi
fi

if [[ -z "$DATA_ROOT" ]]; then
  echo "Unable to resolve allele-freq data root (set ALLELE_FREQ_DATA_DIR)." >&2
  exit 1
fi

SAMPLESHEET="${DATA_ROOT}/samplesheet.csv"
if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "Samplesheet not found: ${SAMPLESHEET}" >&2
  exit 1
fi

PIPELINE_ROOT_DIR="${BV_RESULTS_DIR:-${PWD}/results}"
PIPELINE_DIR="${PIPELINE_ROOT_DIR}/pipeline"
PIPELINE_RESULTS_DIR="${PIPELINE_DIR}/results"
PIPELINE_WORK_DIR="${PIPELINE_DIR}/work"
mkdir -p "$PIPELINE_RESULTS_DIR"
mkdir -p "$PIPELINE_WORK_DIR"

if [[ ! -f "$PIPELINE_DIR/project.yaml" ]]; then
  mkdir -p "$PIPELINE_DIR"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a "$PIPELINE_SRC_DIR"/ "$PIPELINE_DIR"/
  else
    cp -R "$PIPELINE_SRC_DIR"/. "$PIPELINE_DIR"/
  fi
fi

OUTPUT_NPZ="${BV_OUTPUT_ALLELE_FREQ_NPZ:-allele_freq.npz}"
OUTPUT_INDEX="${BV_OUTPUT_LOCUS_INDEX:-locus_index.json}"
OUTPUT_TSV="${BV_OUTPUT_ALLELE_FREQ_TSV:-allele_freq.tsv}"
OUTPUT_VCF_STATS="${BV_OUTPUT_VCF_RESULTS:-vcf_conversion_results.tsv}"
SAMPLESHEET_HASH_FILE="${PIPELINE_RESULTS_DIR}/samplesheet.sha256"

compute_samplesheet_hash() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$SAMPLESHEET" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$SAMPLESHEET" | awk '{print $1}'
  else
    SAMPLESHEET_PATH="$SAMPLESHEET" python3 - <<'PY'
import hashlib
import os
with open(os.environ["SAMPLESHEET_PATH"],"rb") as f:
    print(hashlib.sha256(f.read()).hexdigest())
PY
  fi
}

should_skip=0
if [[ "${ALLELE_FREQ_FORCE_RUN:-0}" != "1" ]]; then
  if [[ -f "${PIPELINE_RESULTS_DIR}/allele_freq.npz" && -f "${PIPELINE_RESULTS_DIR}/locus_index.json" ]]; then
    if [[ "${ALLELE_FREQ_SKIP_IF_DONE:-1}" == "1" ]]; then
      if [[ -f "$SAMPLESHEET_HASH_FILE" ]]; then
        current_hash="$(compute_samplesheet_hash)"
        stored_hash="$(cat "$SAMPLESHEET_HASH_FILE" 2>/dev/null || true)"
        if [[ -n "$stored_hash" && "$stored_hash" == "$current_hash" ]]; then
          should_skip=1
        fi
      fi
    fi
  fi
fi

echo "[gen-allele-freq] results_dir=${PIPELINE_RESULTS_DIR} skip=${should_skip}"

if [[ "$should_skip" != "1" ]]; then
  "$BV_BIN" run "$PIPELINE_DIR" \
    --participants "$SAMPLESHEET" \
    --results-dir "$PIPELINE_RESULTS_DIR"
  compute_samplesheet_hash > "$SAMPLESHEET_HASH_FILE"
else
  echo "[gen-allele-freq] skipping nextflow run (outputs present)"
fi

copy_out() {
  local name="$1"
  local dest="$2"
  local src="${PIPELINE_RESULTS_DIR}/${name}"
  if [[ -f "$src" ]]; then
    if [[ -f "$dest" && "$(realpath "$dest")" == "$(realpath "$src")" ]]; then
      return 0
    fi
    cp -f "$src" "$dest"
    return 0
  fi
  local alt="${PIPELINE_DIR}/results/${name}"
  if [[ -f "$alt" ]]; then
    cp -f "$alt" "$dest"
    return 0
  fi
  return 1
}

copy_out "allele_freq.npz" "$OUTPUT_NPZ"
copy_out "locus_index.json" "$OUTPUT_INDEX"
copy_out "allele_freq.tsv" "$OUTPUT_TSV"
copy_out "vcf_conversion_results.tsv" "$OUTPUT_VCF_STATS"
