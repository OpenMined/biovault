nextflow.enable.dsl=2

workflow USER {
    take:
      context
      participants
    main:
      def assetsDir = System.getenv('BV_ASSETS_DIR') ?: "${projectDir}/assets"
      def assets_dir_ch = Channel.value(file(assetsDir))

      // rs1805007 (MC1R R151C) coordinates:
      //   GRCh38: chr16:89919709
      //   GRCh37: chr16:89986117
      // (use 1-bp range like your eye example)
      def rs1805007 = ['16:89919709-89919709', '16:89986117-89986117']

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
        def base_pos = pickPos(ref_version, rs1805007)
        tuple(participant_id, base_pos, ref, ref_index, aligned, aligned_index)
      }

      def calls = call_region(with_pos)

      def out = interpret_redhair(assets_dir_ch, calls)

      out.msg
        .map { pid, msg -> "\n===== Red Hair (MC1R rs1805007) for Participant: ${pid} =====\n${msg}\n====================================\n" }
        .view()
}

// helper to choose correct coordinate
def pickPos(version, tuple) {
    def v = version?.toLowerCase()
    if(v in ['grch38','hg38']) return tuple[0]
    if(v in ['grch37','hg19']) return tuple[1]
    exit 1, "Unknown ref_version: ${version}"
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
  set -euo pipefail

  # Emit records even for homozygous-reference sites (-A)
  bcftools mpileup -Ou -f !{ref} -r !{base_pos} -a FORMAT/AD,FORMAT/DP !{aligned} \
    | bcftools call -m -A -Oz -o variants.vcf.gz

  # CSI index (robust on large refs)
  bcftools index -f -c variants.vcf.gz

  # Query one line: REF  ALT  GT  AD
  bcftools query -r !{base_pos} -f "%REF\t%ALT\t[%GT]\t[%AD]\n" variants.vcf.gz > snp.txt
  '''
}

process interpret_redhair {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'

  input:
  path assets_dir
  tuple val(participant_id), path(vcf), path(vcf_index), path(snp)

  output:
  path "red_hair.txt"
  tuple val(participant_id), stdout, emit: msg

  shell:
  """
  python3 !{assets_dir}/red_hair.py --participant '!{participant_id}' < !{snp} > red_hair.txt
  cat red_hair.txt
  """
}
