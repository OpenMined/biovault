workflow USER {
    take:
      patient_id_ch
      ref_ch
      ref_index_ch
      aligned_ch
      aligned_index_ch
      ref_version
      assets_dir_ch
      results_dir
    main:
      // (GRCh38_position, GRCh37_position)
      def rs12913832 = ['chr15:28120472-28120472', 'chr15:28365618-28365618']
      def base_pos = pickPos(ref_version, rs12913832)

      def call_region_ch = call_region(
          base_pos, ref_ch, ref_index_ch, aligned_ch, aligned_index_ch
      )

      def out = interpret_eyes(assets_dir_ch, call_region_ch)

      out.msg
        .map { "\n===== Eye Color Interpretation for Patient: ${patient_id_ch.getVal()} =====\n${it}\n====================================\n" }
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
    val base_pos
    path ref
    path ref_index
    path aligned
    path aligned_index

    output:
    tuple path("variants.vcf.gz"), path("variants.vcf.gz.csi"), path("snp.txt")

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
  tuple path(vcf), path(vcf_index), path(snp)

  output:
  path "eye_color.txt"
  stdout emit: msg   // capture stdout as a channel
  shell:
  """
  cat !{snp} | python3 ./assets/eye_color.py > eye_color.txt
  cat eye_color.txt
  """
}
