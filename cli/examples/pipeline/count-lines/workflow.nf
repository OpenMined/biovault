nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext
        samplesheet   // File (CSV)
        data_dir      // Directory containing sample files

    main:
        def assetsDir = file(context.params.assets_dir)
        def countScript = file("${assetsDir}/count_lines.py")
        def script_ch = Channel.value(countScript)
        def output_name = 'line_counts.csv'
        def output_ch = Channel.value(output_name)
        counted_sheet_ch = annotate_line_counts(
            samplesheet,
            data_dir,
            script_ch,
            output_ch
        )

    emit:
        counted_sheet = counted_sheet_ch
        data_dir = data_dir
}

process annotate_line_counts {
    tag { sheet.baseName }
    publishDir params.results_dir, mode: 'copy', pattern: 'line_counts.csv', overwrite: true

    input:
        file sheet
        val data_dir_path
        path count_script
        val output_name

    output:
        path output_name

    script:
    """
    python3 ${count_script} \
        --input "${sheet}" \
        --assets-dir "${data_dir_path}" \
        --output "${output_name}"
    """
}
