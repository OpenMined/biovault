nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext
        samplesheet   // File (CSV)

    main:
        def assetsDir = file(context.params.assets_dir)
        def countScript = file("${assetsDir}/count_lines.py")
        def script_ch = Channel.value(countScript)
        def output_name = 'line_counts.csv'
        def output_ch = Channel.value(output_name)
        counted_sheet_ch = annotate_line_counts(samplesheet, script_ch, output_ch)

    emit:
        counted_sheet = counted_sheet_ch
}

process annotate_line_counts {
    tag { sheet.baseName }
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        file sheet
        path count_script
        val output_name

    output:
        path output_name

    script:
    """
    python3 ${count_script} \
        --input "${sheet}" \
        --output "${output_name}"
    """
}
