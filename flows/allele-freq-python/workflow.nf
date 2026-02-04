// BioVault allele-freq workflow v0.1.0
// Converts genotypes to VCF, extracts allele counts, merges with frequency calculation

nextflow.enable.dsl=2

def BATCH_SIZE = 10

workflow USER {
    take:
        context
        participants  // Channel emitting records with participant_id and input_file

    main:
        def assetsDir = context.params.assets_dir
        if (!assetsDir) {
            throw new IllegalStateException("Missing assets directory in context params")
        }
        def assetsDirPath = file(assetsDir)
        def extract_script = file("${assetsDirPath}/extract_counts.py")
        def aggregate_batch_script = file("${assetsDirPath}/aggregate_batch.py")
        def merge_batches_script = file("${assetsDirPath}/merge_all_batches.py")
        def conversion_stats_script = file("${assetsDirPath}/aggregate_conversion_stats.py")
        def allele_freq_script = file("${assetsDirPath}/calc_allele_freq.py")
        def extract_script_ch = Channel.value(extract_script)
        def aggregate_batch_script_ch = Channel.value(aggregate_batch_script)
        def merge_batches_script_ch = Channel.value(merge_batches_script)
        def conversion_stats_script_ch = Channel.value(conversion_stats_script)
        def allele_freq_script_ch = Channel.value(allele_freq_script)
        // Step 1: Convert genotype files to VCF (or pass through existing VCF)
        def participant_work_items = participants.map { record ->
            tuple(
                record.participant_id,
                file(record.genotype_file)
            )
        }
        def vcf_files = genotype_to_vcf(participant_work_items)

        // Step 2: Extract counts from each VCF
        def per_file_counts = extract_counts(vcf_files.map { pid, vcf, stats -> tuple(pid, vcf) }, extract_script_ch)
        def snp_batched = per_file_counts.buffer(size: BATCH_SIZE, remainder: true)
            .map { batch ->
                def batch_id = UUID.randomUUID().toString().take(8)
                def batch_files = batch.collect { it[1] }
                tuple(batch_id, batch_files)
            }

        def snp_batch_counts = aggregate_batch(snp_batched, aggregate_batch_script_ch)

        merge_all_batches(
            snp_batch_counts.collect(),
            merge_batches_script_ch
        )
        def matrix_tsv_out = merge_all_batches.out.matrix_tsv
        def matrix_npz_out = merge_all_batches.out.matrix_npz
        def locus_index_out = merge_all_batches.out.locus_index
        def participants_out = merge_all_batches.out.participants
        def allele_freq_out = calc_allele_freq(matrix_tsv_out, allele_freq_script_ch)
        def conversion_stats = aggregate_conversion_stats(
            vcf_files.map { pid, vcf, stats -> stats }.collect(),
            conversion_stats_script_ch
        )

    emit:
        dosage_matrix_tsv = matrix_tsv_out  // dosage_matrix.tsv
        dosage_matrix_npz = matrix_npz_out  // dosage_matrix.npz
        locus_index = locus_index_out  // locus_index.txt
        participants = participants_out  // participants.txt
        allele_freq = allele_freq_out  // allele_freq.tsv
        vcf_results = conversion_stats
}

process genotype_to_vcf {
    container 'ghcr.io/openmined/biosynth:latest'
    containerOptions '--entrypoint ""'
    publishDir "${params.results_dir}/vcf", mode: 'copy', overwrite: true, pattern: '*.vcf.gz'
    publishDir "${params.results_dir}/vcf", mode: 'copy', overwrite: true, pattern: '*.vcf.log'
    tag { participant_id }
    errorStrategy { params.nextflow?.error_strategy ?: 'ignore' }
    maxRetries { params.nextflow?.max_retries ?: 0 }

    input:
        tuple val(participant_id), path(input_file)
    output:
        tuple val(participant_id), path("${participant_id}.vcf.gz"), path("${participant_id}_stats.tsv")

    script:
    def inputFileName = input_file.getName()
    def isVcf = inputFileName.endsWith('.vcf') || inputFileName.endsWith('.vcf.gz')
    """
    #!/bin/bash
    set -euo pipefail

    if [ "${isVcf}" = "true" ]; then
        if [[ "${inputFileName}" == *.vcf.gz ]]; then
            cp "${inputFileName}" "${participant_id}.vcf.gz"
        else
            gzip -c "${inputFileName}" > "${participant_id}.vcf.gz"
        fi
        echo -e "participant_id\\tfilename\\tvcf_rows\\tmissing_count\\tinferred_count\\tstatus" > "${participant_id}_stats.tsv"
        VCF_ROWS=\$(zcat "${participant_id}.vcf.gz" | grep -v "^#" | wc -l | tr -d ' ')
        echo -e "${participant_id}\\t${inputFileName}\\t\${VCF_ROWS}\\t0\\t0\\tpass_through" >> "${participant_id}_stats.tsv"
    else
        bvs genotype-to-vcf --input "${inputFileName}" --output "${participant_id}.vcf.gz" --gzip --missing-log "${participant_id}.vcf.log"
        if [ -f "${participant_id}.vcf.log" ]; then
            MISSING_COUNT=\$(grep -c "missing" "${participant_id}.vcf.log" || echo "0")
            INFERRED_COUNT=\$(grep -c "inferred" "${participant_id}.vcf.log" || echo "0")
        else
            MISSING_COUNT=0
            INFERRED_COUNT=0
        fi
        VCF_ROWS=\$(zcat "${participant_id}.vcf.gz" | grep -v "^#" | wc -l | tr -d ' ')
        echo -e "participant_id\\tfilename\\tvcf_rows\\tmissing_count\\tinferred_count\\tstatus" > "${participant_id}_stats.tsv"
        echo -e "${participant_id}\\t${inputFileName}\\t\${VCF_ROWS}\\t\${MISSING_COUNT}\\t\${INFERRED_COUNT}\\tconverted" >> "${participant_id}_stats.tsv"
    fi
    """
}

process extract_counts {
    container 'ghcr.io/openmined/bioscript:latest'
    tag { participant_id }
    errorStrategy { params.nextflow?.error_strategy ?: 'ignore' }
    maxRetries { params.nextflow?.max_retries ?: 0 }

    input:
        tuple val(participant_id), path(vcf_file)
        path extract_script

    output:
        tuple val(participant_id), path("snp_counts_${participant_id}.tsv")

    script:
    """
    python3 ${extract_script} \
      --vcf "${vcf_file}" \
      --participant "${participant_id}" \
      --output "snp_counts_${participant_id}.tsv"
    """
}

process aggregate_batch {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir "${params.results_dir}/batches", mode: 'copy', overwrite: true
    tag { batch_id }
    errorStrategy { params.nextflow?.error_strategy ?: 'ignore' }
    maxRetries { params.nextflow?.max_retries ?: 0 }

    input:
        tuple val(batch_id), path(count_files)
        path aggregate_batch_script

    output:
        path "batch_${batch_id}_counts.tsv"

    script:
    """
    python3 ${aggregate_batch_script} \
      --batch-id "${batch_id}" \
      --counts "${count_files}" \
      --output "batch_${batch_id}_counts.tsv"
    """
}

process merge_all_batches {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path batch_files
        path merge_batches_script

    output:
        path "dosage_matrix.tsv", emit: matrix_tsv
        path "dosage_matrix.npz", emit: matrix_npz
        path "locus_index.txt", emit: locus_index
        path "participants.txt", emit: participants

    script:
    """
    python3 ${merge_batches_script} \
      --batch-files "${batch_files}" \
      --matrix-tsv "dosage_matrix.tsv" \
      --npz "dosage_matrix.npz" \
      --participants "participants.txt" \
      --loci "locus_index.txt"
    """
}

process aggregate_conversion_stats {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path stats_files
        path conversion_stats_script

    output:
        path "vcf_conversion_results.tsv"

    script:
    """
    python3 ${conversion_stats_script} \
      --stats-files "${stats_files}" \
      --output "vcf_conversion_results.tsv"
    """
}

process calc_allele_freq {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path matrix
        path allele_freq_script

    output:
        path "allele_freq.tsv"

    script:
    """
    python3 ${allele_freq_script} \
      --matrix "${matrix}" \
      --output "allele_freq.tsv"
    """
}
