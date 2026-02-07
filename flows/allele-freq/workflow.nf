// BioVault allele-freq workflow v0.1.0
// Rust long-row pipeline (bvs emit-long + aggregate-long)

nextflow.enable.dsl=2

def envInt = { name, fallback ->
    def raw = System.getenv(name)
    if (raw && raw.isInteger()) {
        return raw.toInteger()
    }
    return fallback
}

def DEFAULT_EMIT_MAX_FORKS = envInt('BV_EMIT_MAX_FORKS', 10)
def DEFAULT_AGG_THREADS = envInt('BV_AGG_THREADS', 8)
def EMIT_MAX_FORKS = params.emit_max_forks ?: (params.nextflow?.emit_max_forks ?: (params.nextflow?.max_forks ?: DEFAULT_EMIT_MAX_FORKS))
def AGG_THREADS = params.aggregate_threads ?: (params.nextflow?.aggregate_threads ?: DEFAULT_AGG_THREADS)
def WIN_TMP_COPY = (System.getenv('BV_WIN_TMP_COPY') ?: '').toLowerCase() in ['1','true','yes','on']

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
    cpus { AGG_THREADS > 0 ? AGG_THREADS : 1 }

    input:
        path bvlr_files

    output:
        path "allele_freq.tsv", emit: allele_freq

    script:
    def threads_arg = (AGG_THREADS > 0) ? "--threads ${AGG_THREADS}" : ""
    def copy_block = '''
if [ "${BV_WIN_TMP_COPY:-}" = "1" ] || grep -q "^/c/Users/" bvlr.list; then
  # Copy inputs into container-local /tmp to avoid slow Windows mounts.
  tmpdir="/tmp/bv-agg-${NXF_TASK_ID:-0}"
  echo "[bv] Copying inputs into $tmpdir"
  mkdir -p "$tmpdir"
  xargs -I '{}' cp '{}' "$tmpdir/" < bvlr.list
  find "$tmpdir" -maxdepth 1 -name '*.bvlr' -print | sort > bvlr.list.local
  echo "[bv] Copy complete: $(wc -l < bvlr.list.local) file(s)"
else
  cp bvlr.list bvlr.list.local
fi
'''
    """
    set -euo pipefail

    # Build a stable input list for bvs to avoid ARG_MAX issues.
    # Resolve to absolute paths so bvs can read outside the work dir.
    printf "%s\\n" ${bvlr_files} | xargs -n1 readlink -f | sort > bvlr.list
    export RAYON_NUM_THREADS="${AGG_THREADS}"
    export OMP_NUM_THREADS="${AGG_THREADS}"
    export OPENBLAS_NUM_THREADS="${AGG_THREADS}"
    export MKL_NUM_THREADS="${AGG_THREADS}"
    export BVS_AGG_PROGRESS_EVERY="${'$'}{BVS_AGG_PROGRESS_EVERY:-250000}"

    ${copy_block}

    bvs aggregate-long \\
      --input-list "bvlr.list.local" \\
      --allele-freq-tsv "allele_freq.tsv" \\
      ${threads_arg} 2>&1
    """
}
