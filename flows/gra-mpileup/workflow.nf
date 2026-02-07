nextflow.enable.dsl=2

workflow USER {
    take:
      context
      participants
    main:
      def assetsDir = System.getenv('BV_ASSETS_DIR') ?: "${projectDir}/assets"
      def assets_dir_ch = Channel.value(file(assetsDir))

      // CYP11B1/CYP11B2 region on chr8 (GRCh38 and GRCh37 coordinates)
      // GRCh38: chr8:142869000-143030000
      // GRCh37: chr8:143950000-144110000
      def cyp11b_region = ['8:142869000-143030000', '8:143950000-144110000']

      def per_participant = participants.map { record ->
        def refVersion = (record['ref_version'] ?: record['grch_version'] ?: 'GRCh38').toString()
        tuple(
          record['participant_id'],
          refVersion,
          record['reference_file'],
          record['reference_index'],
          record['aligned_file'],
          record['aligned_index']
        )
      }

      def with_region = per_participant.map { participant_id, ref_version, ref, ref_index, aligned, aligned_index ->
        def base_region = pickRegion(ref_version, cyp11b_region)
        tuple(participant_id, base_region, ref, ref_index, aligned, aligned_index)
      }

      // Detect chr prefix and get adjusted region
      def region_ch = detect_region(with_region).map { pid, region ->
        tuple(pid, region.toString().trim())
      }

      // Join back with original data for depth and slice operations
      def analysis_inputs = region_ch.join(with_region).map { pid, region, base_region, ref, ref_index, aligned, aligned_index ->
        tuple(pid, region, ref, ref_index, aligned, aligned_index)
      }

      // Run depth analysis
      def depth_ch = compute_depth(analysis_inputs)

      // Slice region to BAM
      def slice_ch = slice_region(analysis_inputs)

      // Plot depth
      def plot_ch = plot_depth(assets_dir_ch, depth_ch)

      // Output results
      plot_ch.msg
        .map { pid, msg -> "\n===== GRA CYP11B1/B2 Analysis for Participant: ${pid} =====\n${msg}\n====================================\n" }
        .view()
}

// Helper to choose correct coordinate based on reference version
def pickRegion(version, tuple) {
    def v = version?.toLowerCase()
    if(v in ['grch38','hg38']) return tuple[0]
    if(v in ['grch37','hg19']) return tuple[1]
    exit 1, "Unknown ref_version: ${version}"
}

process detect_region {
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

  input:
  tuple val(participant_id), val(base_region), path(ref), path(ref_index), path(aligned), path(aligned_index)

  output:
  tuple val(participant_id), stdout, emit: region_out

  shell:
  '''
  set -euo pipefail

  BASE_REGION="!{base_region}"
  CHR="${BASE_REGION%%:*}"

  # Check CRAM/BAM header for chromosome naming convention
  # The aligned file's contigs determine what region format to use
  if samtools view -H "!{aligned}" 2>/dev/null | head -n 100 \
    | awk -F'\\t' '$1=="@SQ"{for(i=1;i<=NF;i++){if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}}' \
    | awk -v want="chr${CHR}" 'BEGIN{ok=0} $0==want{ok=1} END{exit !ok}'
  then
    echo "chr${BASE_REGION}"
    exit 0
  fi

  # Check if CRAM has the chromosome without prefix
  if samtools view -H "!{aligned}" 2>/dev/null | head -n 100 \
    | awk -F'\\t' '$1=="@SQ"{for(i=1;i<=NF;i++){if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}}' \
    | awk -v want="${CHR}" 'BEGIN{ok=0} $0==want{ok=1} END{exit !ok}'
  then
    echo "${BASE_REGION}"
    exit 0
  fi

  # Default: no prefix (base region already has no prefix)
  echo "${BASE_REGION}"
  '''
}

process compute_depth {
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
  publishDir params.results_dir, mode: 'copy'
  stageInMode 'symlink'

  input:
  tuple val(participant_id), val(region), path(ref), path(ref_index), path(aligned), path(aligned_index)

  output:
  tuple val(participant_id), val(region), path("depth.txt")

  shell:
  '''
  set -euo pipefail

  # Compute read depth for the CYP11B1/B2 region
  # Use --reference for CRAM files to ensure proper decoding
  samtools depth --reference "!{ref}" -r "!{region}" "!{aligned}" > depth.txt

  # Check if we got any data
  if [ ! -s depth.txt ]; then
    echo "WARNING: No depth data for region !{region}" >&2
    # Create empty file with header
    echo -e "chr\tpos\tdepth" > depth.txt
  fi
  '''
}

process slice_region {
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
  publishDir params.results_dir, mode: 'copy'
  stageInMode 'symlink'

  input:
  tuple val(participant_id), val(region), path(ref), path(ref_index), path(aligned), path(aligned_index)

  output:
  tuple val(participant_id), path("cyp11b_region.bam"), path("cyp11b_region.bam.bai")

  shell:
  '''
  set -euo pipefail

  # Slice the region from CRAM/BAM to a new BAM file
  samtools view -b -h -T "!{ref}" -o cyp11b_region.bam "!{aligned}" "!{region}"

  # Index the output BAM
  samtools index cyp11b_region.bam
  '''
}

process plot_depth {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'

  input:
  path assets_dir
  tuple val(participant_id), val(region), path(depth_txt)

  output:
  path "depth_plot.pdf"
  path "depth_plot.png"
  path "depth.txt"
  tuple val(participant_id), stdout, emit: msg

  shell:
  '''
  set -euo pipefail

  # Install dependencies to /tmp (container filesystem is read-only)
  pip install --no-cache-dir --target=/tmp/pylibs matplotlib numpy || {
    echo "pip install failed" >&2
    exit 1
  }
  export PYTHONPATH=/tmp/pylibs:${PYTHONPATH:-}
  export MPLCONFIGDIR=/tmp/matplotlib

  python3 !{assets_dir}/plot_depth.py \
    --input "!{depth_txt}" \
    --output-pdf "depth_plot.pdf" \
    --output-png "depth_plot.png" \
    --region "!{region}" \
    --participant "!{participant_id}"

  # Copy depth.txt to output only if not already named depth.txt
  if [ "!{depth_txt}" != "depth.txt" ]; then
    cp "!{depth_txt}" depth.txt
  fi
  echo "Generated coverage plot for !{region}"
  '''
}
