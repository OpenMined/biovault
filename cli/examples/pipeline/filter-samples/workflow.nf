nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext
        samplesheet   // File (CSV)
        data_dir      // Directory containing referenced files

    main:
        def assetsDir = file(context.params.assets_dir)
        def filterScript = file("${assetsDir}/filter_samples.py")
        def extensionValue = (context.params?.extension ?: '.txt').toString()
        def extension_ch = Channel.value(extensionValue)
        def script_ch = Channel.value(filterScript)
        def sheet_name = 'filtered_samplesheet.csv'
        def sheet_name_ch = Channel.value(sheet_name)
        filtered_sheet_ch = filter_samplesheet(
            samplesheet,
            data_dir,
            extension_ch,
            script_ch,
            sheet_name_ch
        )

    emit:
        filtered_sheet = filtered_sheet_ch
        filtered_data_dir = data_dir
}

process filter_samplesheet {
    tag { sheet.baseName }
    publishDir params.results_dir, mode: 'copy', pattern: 'filtered_samplesheet.csv', overwrite: true

    input:
        file sheet
        val data_dir_path
        val extension
        path filter_script
        val output_name

    output:
        path output_name

    script:
    """
    python3 ${filter_script} \
        --input "${sheet}" \
        --extension "${extension}" \
        --output "${output_name}" \
        --data-dir "${data_dir_path}"
    """
}
