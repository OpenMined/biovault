nextflow.enable.dsl=2

workflow USER {
    take:
        context       // BiovaultContext
        samplesheet   // File (CSV with line_count)

    main:
        def assetsDir = file(context.params.assets_dir)
        def sumScript = file("${assetsDir}/sum_counts.py")
        def script_ch = Channel.value(sumScript)
        def output_name = 'results.txt'
        def output_ch = Channel.value(output_name)
        result_summary_ch = sum_line_counts(samplesheet, script_ch, output_ch)

    emit:
        result_summary = result_summary_ch
}

process sum_line_counts {
    tag { sheet.baseName }
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        file sheet
        path sum_script
        val output_name

    output:
        path output_name

    script:
    """
    python3 ${sum_script} \
        --input "${sheet}" \
        --output "${output_name}"
    """
}
