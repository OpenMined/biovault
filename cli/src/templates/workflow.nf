workflow USER {
    take:
      participant_id_ch // example: MADHAVA
      ref_ch // example: GRCh38_full_analysis_set_plus_decoy_hla.fa
      ref_index_ch // example: GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
      aligned_ch // example: NA06985.final.cram
      aligned_index_ch // example: NA06985.final.cram.crai
      ref_version // example:GRCh38
      assets_dir_ch // example: /assets
      results_dir // example: /results/MADHAVA
    main:
      // your code here
      // Use .view to observe items flowing through channels
      participant_id_ch.view { "Participant ID: $it" }
      ref_ch.view            { "Reference: $it" }
      ref_index_ch.view      { "Reference Index: $it" }
      aligned_ch.view        { "Aligned: $it" }
      aligned_index_ch.view  { "Aligned Index: $it" }
      assets_dir_ch.view     { "Assets Directory: $it" }
      // These are plain values (not channels), so println works directly
      println "Reference Version: ${ref_version}"
      println "Results Directory: ${results_dir}"      
}
