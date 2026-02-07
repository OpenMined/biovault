nextflow.enable.dsl=2


/*************************************************
 * WORKFLOW — wire the processes
 *************************************************/
workflow USER {
  take:
    context
    participants

  main:
    def assetsDir = System.getenv('BV_ASSETS_DIR') ?: "${projectDir}/assets"
    def assets_dir_ch = Channel.value(file(assetsDir))

    def per_participant = participants.map { record ->
      def refVersion = (record['ref_version'] ?: record['grch_version'] ?: 'GRCh38').toString()
      tuple(
        record['participant_id'],
        refVersion,
        record['reference_file'],
        record['reference_index'],
        record['aligned_file'],
        record['aligned_index']
      )
    }

    def version_ch = per_participant.map { pid, ref_version, ref, ref_index, aligned, aligned_index ->
      tuple(pid, ref_version)
    }

    def panel_ch = select_panel(assets_dir_ch, version_ch)

    def (regions_ch, panel_passthrough_ch) = build_regions(panel_ch)

    def call_inputs = per_participant.map { pid, ref_version, ref, ref_index, aligned, aligned_index ->
      tuple(pid, ref, ref_index, aligned, aligned_index)
    }

    // Detect chromosome naming convention (chrY vs Y) and adjust regions
    def detect_inputs = call_inputs.join(regions_ch)
    def adjusted_regions_ch = detect_chr_prefix(detect_inputs)

    def calls = call_sites(call_inputs.join(adjusted_regions_ch))

    def variants_tsv_ch = query_to_table(calls)

    def interpret_inputs = panel_passthrough_ch
      .join(variants_tsv_ch)
      .join(version_ch)
      .map { pid, panel_tsv, variants_tsv, ref_version ->
        tuple(pid, panel_tsv, variants_tsv, ref_version)
      }

    interpret_haplogroup(
      assets_dir_ch,
      interpret_inputs
    )
}


/*
 * Y-haplogroup mini-pipeline (DSL2)
 * 0) select_panel     — pick GRCh37/38 Y-marker panel from assets dir
 * 1) build_regions    — make comma-separated region list; keep panel passthrough
 * 2) call_sites       — bcftools mpileup+call over those regions (needs .fai/.crai)
 * 3) query_to_table   — compact TSV for Python post-processing
 * 4) interpret        — run Python script; copy pretty-named outputs to results_dir
 */

params.bcftools_container = 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
params.python_container   = 'python:3.11-slim'

/*************************************************
 * STEP 0 — select panel from assets directory
 *************************************************/
process select_panel {
  container params.bcftools_container
  errorStrategy 'terminate'

  tag { "select_panel" }

  input:
    path assets_dir
    tuple val(participant_id), val(ref_version)

  output:
    tuple val(participant_id), path("panel.tsv")

  shell:
  '''
  set -euo pipefail
  PANEL="!{assets_dir}/y_haplogroup_panel.!{ref_version}.tsv"
  echo "Using panel: $PANEL"
  if [ ! -s "$PANEL" ]; then
    echo "ERROR: Missing panel file: $PANEL" >&2
    exit 1
  fi
  cp "$PANEL" panel.tsv
  '''
}

/*************************************************
 * STEP 1 — build comma-separated regions list
 *************************************************/
process build_regions {
  container params.bcftools_container
  errorStrategy 'terminate'

  tag { "regions" }

  input:
    tuple val(participant_id), path(panel_tsv)

  output:
    tuple val(participant_id), path("regions.txt")
    tuple val(participant_id), path("panel.tsv")

  shell:
  '''
  set -euo pipefail
  # Column 2 is 'region' (e.g., chrY:20577481-20577481)
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $2}' !{panel_tsv} | paste -sd, - > regions.txt
  if [ ! -s regions.txt ]; then
    echo "ERROR: No regions parsed from panel.tsv" >&2
    exit 1
  fi
  cp !{panel_tsv} panel.tsv
  '''
}

/*************************************************
 * STEP 1.5 — detect chromosome naming convention
 * Checks if CRAM uses "chrY" or "Y" and adjusts regions
 *************************************************/
process detect_chr_prefix {
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
  errorStrategy 'terminate'

  input:
    tuple val(participant_id), path(ref_fa), path(ref_fa_fai), path(aln), path(aln_index), path(regions_txt)

  output:
    tuple val(participant_id), path("adjusted_regions.txt")

  shell:
  '''
  set -euo pipefail

  # Check if CRAM/BAM uses "chrY" or "Y" by inspecting header
  # Look for @SQ lines with SN:chrY or SN:Y
  HAS_CHR_PREFIX=false
  if samtools view -H "!{aln}" 2>/dev/null | head -n 100 | grep -q "SN:chrY"; then
    HAS_CHR_PREFIX=true
  fi

  # Read the original regions
  REGIONS=$(cat "!{regions_txt}")

  if [ "$HAS_CHR_PREFIX" = "true" ]; then
    # CRAM uses chrY, keep regions as-is
    echo "$REGIONS" > adjusted_regions.txt
    echo "Detected chrY naming convention, keeping original regions" >&2
  else
    # CRAM uses Y without prefix, strip "chr" from regions
    echo "$REGIONS" | sed 's/chr//g' > adjusted_regions.txt
    echo "Detected Y naming convention (no chr prefix), adjusted regions" >&2
  fi
  '''
}

/*************************************************
 * STEP 2 — bcftools mpileup + call on regions
 * Requires: FASTA + .fai; CRAM + .crai (or BAM + .bai)
 *************************************************/
process call_sites {
  container params.bcftools_container
  errorStrategy 'terminate'
  stageInMode 'symlink'

  input:
    tuple val(participant_id), path(ref_fa), path(ref_fa_fai), path(aln), path(aln_index), path(regions_txt)

  output:
    tuple val(participant_id), path("calls.vcf.gz"), path("calls.vcf.gz.csi")

  shell:
  '''
  set -euo pipefail

  # Ensure indexes are where htslib expects them
  # FASTA index
  if [ ! -e "!{ref_fa}.fai" ] && [ "!{ref_fa_fai}" != "!{ref_fa}.fai" ]; then
    cp -f "!{ref_fa_fai}" "!{ref_fa}.fai"
  fi

  # Alignment index (.crai for CRAM or .bai for BAM)
  if [ ! -e "!{aln}.crai" ] && [ ! -e "!{aln}.bai" ]; then
    case "!{aln_index}" in
      *.crai) if [ "!{aln_index}" != "!{aln}.crai" ]; then cp -f "!{aln_index}" "!{aln}.crai"; fi ;;
      *.bai)  if [ "!{aln_index}" != "!{aln}.bai" ]; then cp -f "!{aln_index}" "!{aln}.bai"; fi ;;
      *) echo "ERROR: Unknown alignment index extension: !{aln_index}" >&2; exit 1 ;;
    esac
  fi

  # Regions list (comma-separated)
  REGIONS=$(cat "!{regions_txt}")
  if [ -z "$REGIONS" ]; then
    echo "ERROR: Empty regions list" >&2
    exit 1
  fi

  # Call sites
  bcftools mpileup \
    -f "!{ref_fa}" \
    -r "$REGIONS" \
    -a AD,DP \
    -Ou "!{aln}" \
  | bcftools call -mv -A -Oz -o "calls.vcf.gz"

  bcftools index -f "calls.vcf.gz"
  '''
}

/*************************************************
 * STEP 3 — extract compact table for Python
 *************************************************/
process query_to_table {
  container params.bcftools_container
  errorStrategy 'terminate'

  input:
    tuple val(participant_id), path(vcf_gz), path(vcf_index)

  output:
    tuple val(participant_id), path("variants.tsv")

  shell:
  '''
  set -euo pipefail
  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t[%AD]\t%DP\n' \
    "!{vcf_gz}" \
    > "variants.tsv"
  '''
}

/*************************************************
 * STEP 4 — interpret_haplogroup in Python, write outputs
 *************************************************/
process interpret_haplogroup {
  container params.python_container
  errorStrategy 'terminate'
  publishDir params.results_dir, mode: 'copy'

  input:
    path assets_dir
    tuple val(participant_id), path(panel_tsv), path(variants_tsv), val(ref_version)

  output:
    path "report.txt"
    path "calls.tsv"

  shell:
  '''
  set -euo pipefail
  python3 --version >/dev/null

  python3 "!{assets_dir}/interpret_y_haplogroup.py" \
    --panel "!{panel_tsv}" \
    --variants "!{variants_tsv}" \
    --participant "!{participant_id}" \
    --ref-version "!{ref_version}" \
    --out-report "report.txt" \
    --out-variants "calls.tsv"

  # Pretty names (published by Nextflow)
  cp report.txt "!{participant_id}.y_haplogroup.report.txt"
  cp calls.tsv  "!{participant_id}.y_haplogroup.calls.tsv"
  '''
}
