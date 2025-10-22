nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext
        samplesheet   // File (CSV)

    main:
        def assetsDir = file(context.params.assets_dir)
        def filterScript = file("${assetsDir}/filter_samples.py")
        def extensionValue = (context.params?.extension ?: '.txt').toString()
        def extension_ch = Channel.value(extensionValue)
        def script_ch = Channel.value(filterScript)
        def output_name = 'filtered_samplesheet.csv'
        def output_ch = Channel.value(output_name)
        filtered_sheet_ch = filter_samplesheet(samplesheet, extension_ch, script_ch, output_ch)

    emit:
        filtered_sheet = filtered_sheet_ch
}

process filter_samplesheet {
    tag { sheet.baseName }
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        file sheet
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
        --output "${output_name}"
    """
}
