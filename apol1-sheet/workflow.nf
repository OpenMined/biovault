/*
 * Module: USER
 * Run the APOL1 classifier once per samplesheet row.
 * - Make assets_dir a value input so it’s available to every task.
 * - Map rows -> (pid, snp_file, age)
 * - No combine needed; process takes two inputs: value path + tuple.
 */

workflow USER {
  take:
    sample_sheet_ch   // kept for signature parity (unused)
    rows_ch           // emits maps: { participant_id, genotype_file_path, (optional) age }
    mapping_ch        // kept for signature parity (unused)
    assets_dir_ch     // MUST be a value channel (Channel.value(file(params.assets_dir)))
    results_dir       // pass-through value from the template

  main:
    /*
     * 1) Map rows -> (pid, snp_file, age)
     */
    def to_classify_ch = rows_ch.map { row ->
      def pid = row.participant_id as String
      def snp_file = file(row.genotype_file_path)

      if (!pid)               error "Row missing participant_id: ${row}"
      if (!snp_file.exists()) error "Genotype file not found for PID ${pid}: ${row.genotype_file_path}"

      def age =
        (row.age instanceof Number) ? (row.age as Number).intValue() :
        (row.age instanceof String && row.age.isInteger()) ? row.age.toInteger() :
        null

      tuple(pid, snp_file, age)
    }

    // DEBUG
    to_classify_ch.view { t -> "To classify: [${t[0]}, ${t[1]}, ${t[2]}]" }

    /*
     * 2) Run one task per row
     *    Pass assets_dir as a reusable value input + the per-row tuple
     */
    def apol1 = apol1_classifier(assets_dir_ch, to_classify_ch)

  emit:
    apol1_rsids = apol1.rsid_table
    apol1_genos = apol1.genotype_line
    apol1_evid  = apol1.evidence
}

process apol1_classifier {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  tag { pid }
  debug true

  input:
    // assets_dir must be a value channel so it’s reused for every task
    path assets_dir
    tuple val(pid), path(snp_file), val(age)

  output:
    path "apol1_${pid}.rsid.tsv",     emit: rsid_table
    path "apol1_${pid}.genotype.txt", emit: genotype_line
    path "apol1_${pid}.evidence.txt", emit: evidence

  script:
  """
  echo ">>> PROCESS START pid=${pid}"
  echo "assets_dir contents:"
  ls -la "${assets_dir}" || true
  echo "snp_file path: ${snp_file}"
  echo "age: ${age}"
  set -euo pipefail
  APOL1_OUT_PREFIX="apol1_${pid}" \\
    python3 "${assets_dir}/apol1_classifier.py" < "${snp_file}"
  echo ">>> PROCESS END pid=${pid}"
  """
}
