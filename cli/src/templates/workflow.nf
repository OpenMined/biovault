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
      // your code here
      println "Patient ID Channel: ${patient_id_ch}"
      println "Reference Channel: ${ref_ch}"
      println "Reference Index Channel: ${ref_index_ch}"
      println "Aligned Channel: ${aligned_ch}"
      println "Aligned Index Channel: ${aligned_index_ch}"
      println "Reference Version: ${ref_version}"
      println "Assets Directory Channel: ${assets_dir_ch}"
      println "Results Directory: ${results_dir}"
}