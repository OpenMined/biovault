// BioVault allele-freq workflow v0.1.0
// Rust long-row pipeline (bvs emit-long + aggregate-long)

nextflow.enable.dsl=2

def EMIT_MAX_FORKS = params.emit_max_forks ?: (params.nextflow?.emit_max_forks ?: 0)
def AGG_THREADS = params.aggregate_threads ?: (params.nextflow?.aggregate_threads ?: 0)

workflow USER {
    take:
        context
        participants  // Channel emitting records with participant_id and input_file

    main:
        // Step 1: emit compact long rows (genotype or VCF)
        def participant_work_items = participants.map { record ->
            tuple(
                record.participant_id,
                file(record.genotype_file)
            )
        }
        def long_rows = emit_long(participant_work_items)

        // Step 2: aggregate long rows into allele frequencies (and matrix, if needed)
        def agg_out = aggregate_long(long_rows.collect())

    emit:
        allele_freq = agg_out.allele_freq  // allele_freq.tsv
        dosage_matrix = agg_out.dosage_matrix
}

process emit_long {
    container 'ghcr.io/openmined/biosynth:latest'
    containerOptions '--entrypoint ""'
    tag { participant_id }
    errorStrategy { params.nextflow?.error_strategy ?: 'ignore' }
    maxRetries { params.nextflow?.max_retries ?: 0 }
    maxForks EMIT_MAX_FORKS

    input:
        tuple val(participant_id), path(input_file)
    output:
        path "${participant_id}.bvlr"

    script:
    def inputFileName = input_file.getName()
    def isVcf = inputFileName.endsWith('.vcf') || inputFileName.endsWith('.vcf.gz')
    """
    #!/bin/bash
    set -euo pipefail

    if [ "${isVcf}" = "true" ]; then
        bvs emit-long \\
          --vcf "${inputFileName}" \\
          --output "${participant_id}.bvlr" \\
          --participant "${participant_id}" 2>&1
    else
        bvs emit-long \\
          --input "${inputFileName}" \\
          --output "${participant_id}.bvlr" \\
          --participant "${participant_id}" 2>&1
    fi
    """
}

process aggregate_long {
    container 'ghcr.io/openmined/biosynth:latest'
    containerOptions '--entrypoint ""'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'link'

    input:
        path bvlr_files

    output:
        path "allele_freq.tsv", emit: allele_freq
        path "dosage_matrix.tsv", emit: dosage_matrix

    script:
    def threads_arg = (AGG_THREADS > 0) ? "--threads ${AGG_THREADS}" : ""
    """
    set -euo pipefail

    # Build a stable input list for bvs to avoid ARG_MAX issues.
    # Resolve to absolute paths so bvs can read outside the work dir.
    printf "%s\\n" ${bvlr_files} | xargs -n1 readlink -f | sort > bvlr.list

    bvs aggregate-long \\
      --input-list "bvlr.list" \\
      --matrix-tsv "dosage_matrix.tsv" \\
      --allele-freq-tsv "allele_freq.tsv" \\
      ${threads_arg} 2>&1
    """
}
