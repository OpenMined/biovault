nextflow.enable.dsl=2
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowBroadcast

workflow USER {
    take:
      context
      participants
      clinvar_vcf
      clinvar_index

    main:
      def filterConfigPath = "${projectDir}/test/filter_config.json"
      if (!new File(filterConfigPath).exists()) {
        throw new IllegalArgumentException("Missing filter config at ${filterConfigPath}")
      }
      def filter_config_ch = Channel.value(file(filterConfigPath))

      def toChannel = { value ->
        if (value instanceof DataflowReadChannel || value instanceof DataflowBroadcast) {
          return value
        }
        return Channel.value(value)
      }

      def clinvar_vcf_ch = toChannel(clinvar_vcf).map { file(it) }
      def clinvar_index_ch = toChannel(clinvar_index).map { file(it) }
      def clinvar_ch = clinvar_vcf_ch
        .combine(clinvar_index_ch)
        .map { cv, ci -> tuple(cv, ci) }

      def per_participant = participants
        .splitCsv(header: true)
        .map { row ->
          def pid = row.participant_id
          def vcfValue = row.vcf_file
          def idxValue = row.vcf_index
          tuple(
            pid?.toString(),
            file(vcfValue.toString()),
            (idxValue && idxValue.toString().trim()) ? file(idxValue.toString()) : null
          )
        }

      // Ensure VCF is indexed
      def indexed_vcf = ensure_index(per_participant)

      // Intersect with ClinVar
      def intersected = intersect_clinvar(indexed_vcf, clinvar_ch)

      // Generate report
      def report = generate_report(filter_config_ch, intersected)

      // Output results
      report.msg
        .map { pid, msg -> "\n===== ClinVar Analysis for Participant: ${pid} =====\n${msg}\n====================================\n" }
        .view()
}

process ensure_index {
  container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
  stageInMode 'symlink'

  input:
  tuple val(participant_id), path(vcf), val(vcf_index)

  output:
  tuple val(participant_id), path("input.vcf.gz"), path("input.vcf.gz.tbi")

  shell:
  '''
  set -euo pipefail

  # Copy or compress input VCF
  if [[ "!{vcf}" == *.gz ]]; then
    cp "!{vcf}" input.vcf.gz
  else
    bgzip -c "!{vcf}" > input.vcf.gz
  fi

  # Index if not provided or invalid
  if [[ "!{vcf_index}" != "null" ]] && [[ -n "!{vcf_index}" ]] && [[ -f "!{vcf_index}" ]]; then
    cp "!{vcf_index}" input.vcf.gz.tbi 2>/dev/null || bcftools index -t input.vcf.gz
  else
    bcftools index -t input.vcf.gz
  fi
  '''
}

process intersect_clinvar {
  container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
  publishDir params.results_dir, mode: 'copy'
  stageInMode 'symlink'

  input:
  tuple val(participant_id), path(vcf), path(vcf_index)
  tuple path(clinvar_vcf), path(clinvar_index)

  output:
  tuple val(participant_id), path("clinvar_matched.vcf.gz"), path("clinvar_matched.vcf.gz.tbi"), path("clinvar_annotations.tsv")

  shell:
  '''
  set -euo pipefail

  # Use bcftools annotate to transfer ClinVar annotations to matching variants
  # This annotates our VCF with ClinVar INFO fields where positions match
  bcftools annotate \
    -a "!{clinvar_vcf}" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNVC,INFO/GENEINFO,INFO/RS,INFO/ALLELEID,INFO/MC,INFO/ORIGIN \
    -Oz -o clinvar_matched.vcf.gz \
    "!{vcf}"

  bcftools index -t clinvar_matched.vcf.gz

  # Extract annotations to TSV - only variants that have CLNSIG annotation (i.e., found in ClinVar)
  bcftools query \
    -i 'INFO/CLNSIG!="."' \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RS\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNVC\t%INFO/ALLELEID\n' \
    clinvar_matched.vcf.gz > clinvar_annotations.tsv || true

  # Create empty file if no matches
  if [ ! -s clinvar_annotations.tsv ]; then
    touch clinvar_annotations.tsv
  fi
  '''
}

process generate_report {
  container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
  publishDir params.results_dir, mode: 'copy'

  input:
  path filter_config
  tuple val(participant_id), path(matched_vcf), path(matched_vcf_index), path(annotations_tsv)

  output:
  path "clinvar_report.tsv"
  tuple val(participant_id), stdout, emit: msg

  shell:
  '''
  set -euo pipefail

  # Parse filter config for CLNSIG include patterns
  # Extract include patterns from JSON (simple grep/sed approach for shell-only)
  INCLUDE_PATTERNS=$(grep -A 20 '"include"' "!{filter_config}" | grep '"' | grep -v 'include\\|exclude\\|]' | sed 's/.*"\\([^"]*\\)".*/\\1/' | tr '\\n' '|' | sed 's/|$//')

  # Create header with URLs
  echo -e "CHROM\\tPOS\\tREF\\tALT\\tRS\\tGENEINFO\\tCLNSIG\\tCLNREVSTAT\\tCLNVC\\tdbSNP_URL\\tClinVar_URL" > clinvar_report.tsv

  # Filter and add URLs
  if [ -s "!{annotations_tsv}" ]; then
    awk -v patterns="$INCLUDE_PATTERNS" '
    BEGIN {
      FS="\\t"
      OFS="\\t"
      n = split(patterns, pat_arr, "|")
    }
    {
      clnsig = $7
      matched = 0
      for (i = 1; i <= n; i++) {
        if (tolower(clnsig) ~ tolower(pat_arr[i])) {
          matched = 1
          break
        }
      }
      if (matched) {
        rs = $5
        alleleid = $10
        dbsnp_url = (rs != "." && rs != "") ? "https://www.ncbi.nlm.nih.gov/snp/rs" rs : ""
        clinvar_url = (alleleid != "." && alleleid != "") ? "https://www.ncbi.nlm.nih.gov/clinvar/variation/" alleleid "/" : ""
        # Output: CHROM POS REF ALT RS GENEINFO CLNSIG CLNREVSTAT CLNVC dbSNP_URL ClinVar_URL
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, dbsnp_url, clinvar_url
      }
    }' "!{annotations_tsv}" >> clinvar_report.tsv
  fi

  # Count results
  TOTAL=$(tail -n +2 clinvar_report.tsv | wc -l | tr -d ' ')

  echo "Found ${TOTAL} clinically significant variants matching filter criteria"
  echo "Report saved to clinvar_report.tsv"

  # Show summary by CLNSIG
  if [ "$TOTAL" -gt 0 ]; then
    echo ""
    echo "Summary by clinical significance:"
    tail -n +2 clinvar_report.tsv | cut -f7 | sort | uniq -c | sort -rn
  fi
  '''
}
