nextflow.enable.dsl=2

workflow USER {
    take:
      context
      participants
    main:
      def assetsDir = System.getenv('BV_ASSETS_DIR') ?: "${projectDir}/assets"
      def assets_dir_ch = Channel.value(file(assetsDir))

      // (GRCh38_position, GRCh37_position)
      def rs12913832 = ['15:28120472-28120472', '15:28365618-28365618']

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

      def with_pos = per_participant.map { participant_id, ref_version, ref, ref_index, aligned, aligned_index ->
        def base_pos = pickPos(ref_version, rs12913832)
        tuple(participant_id, base_pos, ref, ref_index, aligned, aligned_index)
      }

      def region_ch = detect_region(with_pos).map { pid, region ->
        tuple(pid, region.toString().trim())
      }

      def call_inputs = region_ch.join(with_pos).map { pid, region, base_pos, ref, ref_index, aligned, aligned_index ->
        tuple(pid, region, ref, ref_index, aligned, aligned_index)
      }

      def calls = call_region(call_inputs)

      def interpreted = interpret_eyes(assets_dir_ch, calls)

      interpreted.msg
        .map { pid, msg -> "\n===== Eye Color Interpretation for Participant: ${pid} =====\n${msg}\n====================================\n" }
        .view()
}

// helper to choose correct coordinate
def pickPos(version, tuple) {
    def v = version?.toLowerCase()
    if(v in ['grch38','hg38']) return tuple[0]
    if(v in ['grch37','hg19']) return tuple[1]
    exit 1, "Unknown ref_version: ${version}"
}

process detect_region {
  // Lightweight header-only check
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

  input:
  tuple val(participant_id), val(base_pos), path(ref), path(ref_index), path(aligned), path(aligned_index)

  output:
  tuple val(participant_id), stdout, emit: region_out

  shell:
  '''
  set -euo pipefail

  BASE_POS="!{base_pos}"
  CHR="${BASE_POS%%:*}"

  # Extract first 50 header lines and collect SN: contig names
  # Example header lines: @SQ  SN:chr1  LN:248956422 ...
  if samtools view -H "!{aligned}" | head -n 50 \
    | awk -F'\\t' '$1=="@SQ"{for(i=1;i<=NF;i++){if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}}' \
    | awk -v want="chr${CHR}" 'BEGIN{ok=0} $0==want{ok=1} END{exit !ok}'
  then
    echo "chr${BASE_POS}"
    exit 0
  fi

  # Fallback: check FASTA index (tiny/fast)
  if awk -v want="chr${CHR}" 'BEGIN{ok=0} $1==want{ok=1} END{exit !ok}' "!{ref_index}"; then
    echo "chr${BASE_POS}"
    exit 0
  fi

  # Default: no prefix
  echo "${BASE_POS}"
  '''
}

process call_region {
    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
    publishDir params.results_dir, mode: 'copy'
    stageInMode 'symlink'

    input:
    tuple val(participant_id), val(base_pos), path(ref), path(ref_index), path(aligned), path(aligned_index)

    output:
    tuple val(participant_id), path("variants.vcf.gz"), path("variants.vcf.gz.csi"), path("snp.txt")

    shell:
    '''
    # Emit records even for homozygous-reference sites
    bcftools mpileup -Ou -f !{ref} -r !{base_pos} -a FORMAT/AD,FORMAT/DP !{aligned} \
      | bcftools call -m -A -Oz -o variants.vcf.gz

    # Use CSI index (robust for large references)
    bcftools index -f -c variants.vcf.gz

    # Query the variant and pass to python script
    bcftools query -r !{base_pos} -f "%REF\t%ALT\t[%GT]\t[%AD]\n" variants.vcf.gz > snp.txt
    '''
}

process interpret_eyes {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  input:
  path assets_dir
  tuple val(participant_id), path(vcf), path(vcf_index), path(snp)

  output:
  path "eye_color.txt"
  tuple val(participant_id), stdout, emit: msg
  shell:
  """
  cat !{snp} | python3 ./assets/eye_color.py > eye_color.txt
  cat eye_color.txt
  """
}
