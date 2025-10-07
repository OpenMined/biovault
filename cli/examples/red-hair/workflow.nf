nextflow.enable.dsl=2

workflow USER {
    take:
      participant_id_ch
      ref_ch
      ref_index_ch
      aligned_ch
      aligned_index_ch
      ref_version
      assets_dir_ch
      results_dir
    main:
      // rs1805007 (MC1R R151C) coordinates:
      //   GRCh38: chr16:89919709
      //   GRCh37: chr16:89986117
      // (use 1-bp range like your eye example)
      def rs1805007 = ['16:89919709-89919709', '16:89986117-89986117']
      def base_pos = pickPos(ref_version, rs1805007)

      def call_region_ch = call_region(
        base_pos, ref_ch, ref_index_ch, aligned_ch, aligned_index_ch
      )

      def out = interpret_redhair(assets_dir_ch, call_region_ch, participant_id_ch, results_dir)

      out.msg
        .map { "\n===== Red Hair (MC1R rs1805007) for Participant: ${participant_id_ch.getVal()} =====\n${it}\n====================================\n" }
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

  input:
  val  base_pos
  path ref
  path ref_index
  path aligned
  path aligned_index

  output:
  tuple path("variants.vcf.gz"), path("variants.vcf.gz.csi"), path("snp.txt")

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
  tuple path(vcf), path(vcf_index), path(snp)
  val participant_id
  val results_dir

  output:
  path "red_hair.txt"
  stdout emit: msg

  shell:
  """
  python3 !{assets_dir}/red_hair.py --participant '!{participant_id}' < !{snp} > red_hair.txt
  cat red_hair.txt
  """
}
