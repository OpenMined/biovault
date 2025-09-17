workflow USER {
    take:
      participant_id_ch // example: MADHAVA
      snp_ch // example: NA07357.snp
      assets_dir_ch // example: /assets
      results_dir // example: /results/MADHAVA
    main:
      // your code here
      // Use .view to observe items flowing through channels
      participant_id_ch.view { "Participant ID: $it" }
      snp_ch.view            { "SNP file: $it" }
      assets_dir_ch.view     { "Assets Directory: $it" }
      // This is a plain value (not channel), so println works directly
      println "Results Directory: ${results_dir}"      

      def out = count_number_of_snps(assets_dir_ch, snp_ch)

      out.msg
        .map { "\n===== Number of SNPs: ${it} =====\n" }
        .view()
}

process count_number_of_snps {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  input:
  path assets_dir
  path snp_ch

  output:
  path "number_of_snps.txt"
  stdout emit: msg   // capture stdout as a channel
  shell:
  """
  cat !{snp_ch} | python3 ./assets/count_number_of_snps.py > number_of_snps.txt
  cat number_of_snps.txt
  """
}
