#!/usr/bin/env bash
set -euo pipefail
echo "share-data: copying allele_freq TSV to output"
cp "$BV_INPUT_ALLELE_FREQ_TSV" "$BV_OUTPUT_ALLELE_FREQ_TSV"
echo "share-data: done ($(wc -l < "$BV_OUTPUT_ALLELE_FREQ_TSV") lines)"
