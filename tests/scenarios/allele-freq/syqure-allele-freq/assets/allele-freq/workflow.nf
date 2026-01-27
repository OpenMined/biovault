// BioVault allele-freq workflow v0.1.0
// Converts genotypes to VCF, extracts allele counts, merges with frequency calculation

nextflow.enable.dsl=2

def BATCH_SIZE = 10
def MODE = "dosage"  // "dosage" (0/1/2) or "carrier" (0/1)

workflow USER {
    take:
        context
        participants  // Channel emitting records with participant_id and input_file

    main:
        // Step 1: Convert genotype files to VCF (or pass through existing VCF)
        def participant_work_items = participants.map { record ->
            tuple(
                record.participant_id,
                file(record.genotype_file)
            )
        }
        def vcf_files = genotype_to_vcf(participant_work_items)

        // Step 2: Extract counts from each VCF
        def per_file_counts = extract_counts(vcf_files.map { pid, vcf, stats -> tuple(pid, vcf) })

        // Step 3: Batch the counts files and aggregate in parallel
        def batched = per_file_counts.buffer(size: BATCH_SIZE, remainder: true)
            .map { batch ->
                def batch_id = UUID.randomUUID().toString().take(8)
                tuple(batch_id, batch)
            }

        def batch_counts = aggregate_batch(batched)

        // Step 4: Final merge of all batches
        def final_result = merge_all_batches(batch_counts.collect())

        // Collect conversion stats
        def conversion_stats = aggregate_conversion_stats(
            vcf_files.map { pid, vcf, stats -> stats }.collect()
        )

    emit:
        allele_freq = final_result.map { it[0] }  // allele_freq.tsv
        allele_freq_npz = final_result.map { it[1] }  // allele_freq.npz
        locus_index = final_result.map { it[2] }  // locus_index.json
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

    output:
        path "counts_${participant_id}.tsv"

    script:
    def mode = MODE
    """
    #!/usr/bin/env python3
    import gzip

    mode = "${mode}"
    vcf_file = "${vcf_file}"
    participant_id = "${participant_id}"

    counts = {}  # locus -> [ac, an, rsid]

    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 10:
                continue

            chrom, pos, rsid, ref, alt = parts[0:5]
            sample = parts[9]

            locus = f"{chrom}:{pos}:{ref}:{alt}"
            gt = sample.split(':')[0]

            # Skip missing
            if gt in ('.', './.', '.|.'):
                continue

            # Parse genotype
            gt_clean = gt.replace('|', '/')
            alleles = gt_clean.split('/')
            if len(alleles) != 2 or '.' in alleles:
                continue

            a, b = int(alleles[0]), int(alleles[1])

            if mode == "dosage":
                ac = (1 if a > 0 else 0) + (1 if b > 0 else 0)
                an = 2
            else:  # carrier
                ac = 1 if (a > 0 or b > 0) else 0
                an = 1

            # Keep rsid if valid (not '.')
            rsid_val = rsid if rsid and rsid != '.' else ''

            if locus not in counts:
                counts[locus] = [0, 0, rsid_val]
            counts[locus][0] += ac
            counts[locus][1] += an
            # Keep first non-empty rsid we see
            if rsid_val and not counts[locus][2]:
                counts[locus][2] = rsid_val

    with open(f"counts_{participant_id}.tsv", 'w') as out:
        out.write("locus\\tac\\tan\\trsid\\n")
        for locus in sorted(counts.keys()):
            ac, an, rsid = counts[locus]
            out.write(f"{locus}\\t{ac}\\t{an}\\t{rsid}\\n")
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

    output:
        path "batch_${batch_id}_counts.tsv"

    script:
    """
    #!/usr/bin/env python3
    import os
    from collections import defaultdict

    batch_id = "${batch_id}"
    count_files = "${count_files}".split()

    counts = defaultdict(lambda: [0, 0, 0, ''])  # locus -> [ac, an, n_samples, rsid]

    for cf in count_files:
        if not os.path.exists(cf):
            continue
        with open(cf) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) < 4:
                    continue
                locus, ac, an, rsid = parts[0], int(parts[1]), int(parts[2]), parts[3] if len(parts) > 3 else ''
                counts[locus][0] += ac
                counts[locus][1] += an
                counts[locus][2] += 1
                # Keep first non-empty rsid
                if rsid and not counts[locus][3]:
                    counts[locus][3] = rsid

    with open(f"batch_{batch_id}_counts.tsv", 'w') as out:
        out.write("locus\\tac\\tan\\tn_samples\\trsid\\n")
        for locus in sorted(counts.keys()):
            ac, an, n, rsid = counts[locus]
            out.write(f"{locus}\\t{ac}\\t{an}\\t{n}\\t{rsid}\\n")
    """
}

process merge_all_batches {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path batch_files

    output:
        tuple path("allele_freq.tsv"), path("allele_freq.npz"), path("locus_index.json")

    script:
    """
    #!/usr/bin/env python3
    import os
    import json
    import numpy as np
    from collections import defaultdict

    batch_files = "${batch_files}".split()

    counts = defaultdict(lambda: [0, 0, 0, ''])  # locus -> [ac, an, n_samples, rsid]

    for bf in batch_files:
        if not os.path.exists(bf):
            continue
        with open(bf) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) < 5:
                    continue
                locus, ac, an, n, rsid = parts[0], int(parts[1]), int(parts[2]), int(parts[3]), parts[4] if len(parts) > 4 else ''
                counts[locus][0] += ac
                counts[locus][1] += an
                counts[locus][2] += n
                # Keep first non-empty rsid
                if rsid and not counts[locus][3]:
                    counts[locus][3] = rsid

    # Build canonical sorted locus index
    sorted_loci = sorted(counts.keys())
    n_loci = len(sorted_loci)

    # Build numpy arrays aligned to index
    ac_arr = np.zeros(n_loci, dtype=np.int64)
    an_arr = np.zeros(n_loci, dtype=np.int64)
    n_samples_arr = np.zeros(n_loci, dtype=np.int64)
    rsid_list = []

    for i, locus in enumerate(sorted_loci):
        ac, an, n, rsid = counts[locus]
        ac_arr[i] = ac
        an_arr[i] = an
        n_samples_arr[i] = n
        rsid_list.append(rsid)

    # Write TSV
    with open("allele_freq.tsv", 'w') as out:
        out.write("locus\\trsid\\tac\\tan\\taf\\tn_samples\\n")
        for i, locus in enumerate(sorted_loci):
            ac, an, n = int(ac_arr[i]), int(an_arr[i]), int(n_samples_arr[i])
            rsid = rsid_list[i]
            af = ac / an if an > 0 else 0.0
            out.write(f"{locus}\\t{rsid}\\t{ac}\\t{an}\\t{af:.6f}\\t{n}\\n")

    # Write numpy arrays (for secure aggregation)
    np.savez_compressed(
        "allele_freq.npz",
        ac=ac_arr,
        an=an_arr,
        n_samples=n_samples_arr
    )

    # Write locus index (canonical ordering for alignment)
    index_data = {
        "version": "1.0",
        "n_loci": n_loci,
        "loci": sorted_loci,
        "rsids": rsid_list
    }
    with open("locus_index.json", 'w') as f:
        json.dump(index_data, f)
    """
}

process aggregate_conversion_stats {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path stats_files

    output:
        path "vcf_conversion_results.tsv"

    script:
    """
    #!/usr/bin/env python3
    import os

    stats_files = "${stats_files}".split()

    with open("vcf_conversion_results.tsv", 'w') as out:
        out.write("participant_id\\tfilename\\tvcf_rows\\tmissing_count\\tinferred_count\\tstatus\\n")
        for sf in stats_files:
            if not os.path.exists(sf):
                continue
            with open(sf) as f:
                next(f)  # skip header
                for line in f:
                    out.write(line)
    """
}
