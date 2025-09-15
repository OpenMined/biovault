nextflow.enable.dsl=2


/*************************************************
 * WORKFLOW — wire the processes
 *************************************************/
workflow USER {
  take:
    participant_id_ch
    ref_ch                // FASTA
    ref_index_ch          // FASTA .fai
    aligned_ch            // CRAM/BAM
    aligned_index_ch      // .crai/.bai
    ref_version           // 'GRCh38' or 'GRCh37'
    assets_dir_ch         // directory with y_haplogroup_panel.*.tsv + python script
    results_dir

  main:
    // 0) panel
    def panel_ch = select_panel(assets_dir_ch, ref_version)

    // 1) regions + passthrough panel (two distinct outputs)
    def (regions_ch, panel_passthrough_ch) = build_regions(panel_ch)

    // 2) calls (two outputs; use only the VCF)
    def (vcf_ch, vcf_index_ch) = call_sites(
      ref_ch,
      ref_index_ch,
      aligned_ch,
      aligned_index_ch,
      regions_ch
    )

    // 3) compact table (single output: variants.tsv)
    def variants_tsv_ch = query_to_table(vcf_ch)

    // 4) interpret
    interpret_haplogroup(
      participant_id_ch,
      ref_version,
      panel_passthrough_ch,
      variants_tsv_ch,
      assets_dir_ch,
      results_dir
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
    val  ref_version

  output:
    path "panel.tsv"

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
    path panel_tsv

  output:
    path "regions.txt"
    path "panel.tsv"

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
 * STEP 2 — bcftools mpileup + call on regions
 * Requires: FASTA + .fai; CRAM + .crai (or BAM + .bai)
 *************************************************/
process call_sites {
  container params.bcftools_container
  errorStrategy 'terminate'

  input:
    path ref_fa
    path ref_fa_fai
    path aln
    path aln_index
    path regions_txt

  output:
    path "calls.vcf.gz"
    path "calls.vcf.gz.csi"

  shell:
  '''
  set -euo pipefail

  # Ensure indexes are where htslib expects them
  # FASTA index
  if [ ! -e "!{ref_fa}.fai" ]; then
    ln -sf "!{ref_fa_fai}" "!{ref_fa}.fai"
  fi

  # Alignment index (.crai for CRAM or .bai for BAM)
  if [ ! -e "!{aln}.crai" ] && [ ! -e "!{aln}.bai" ]; then
    case "!{aln_index}" in
      *.crai) ln -sf "!{aln_index}" "!{aln}.crai" ;;
      *.bai)  ln -sf "!{aln_index}" "!{aln}.bai"  ;;
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
    path vcf_gz

  output:
    path "variants.tsv"

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

  input:
    val participant_id_ch
    val  ref_version
    path panel_tsv
    path variants_tsv
    path assets_dir
    val  results_dir

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
    --participant "!{participant_id_ch}" \
    --ref-version "!{ref_version}" \
    --out-report "report.txt" \
    --out-variants "calls.tsv"

  # Pretty names in results_dir
  mkdir -p "!{results_dir}"
  cp report.txt "!{results_dir}/!{participant_id_ch}.y_haplogroup.report.txt"
  cp calls.tsv  "!{results_dir}/!{participant_id_ch}.y_haplogroup.calls.tsv"

  printf "Y haplogroup summary written: %s\n" "!{results_dir}/!{participant_id_ch}.y_haplogroup.report.txt"
  '''
}
