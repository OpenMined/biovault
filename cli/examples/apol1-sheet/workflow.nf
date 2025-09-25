/*
 * Module: USER (DSL2-safe; Python writer)
 * - Buffers rows once, fans out cleanly (no `into`)
 * - APOL1 classification per PID
 * - Left-join computed_genotype
 * - Adds computed_genotype to sheet
 * - Writes augmented CSV/TSV via Python (robust header mapping)
 */

import nextflow.Channel

nextflow.enable.dsl=2
params.results_dir = params.results_dir ?: 'results'

workflow USER {
  take:
    sample_sheet_ch   // value path to original sheet
    rows_ch           // queue of row maps (schema already applied)
    mapping_ch        // value { orig_headers, types, required, rename, defaults }
    assets_dir_ch     // value path to assets dir (has apol1_classifier.py)
    results_dir       // string path

  main:
    file(params.results_dir).mkdirs()

    //------------------------------------------------------------------
    // 0) Buffer rows ONCE, then replay
    //------------------------------------------------------------------
    def rows_all_v = rows_ch.collect()                   // VALUE: List<Map>
    rows_all_v.view { lst -> "Buffered rows: ${lst.size()}" }

    def rows_for_classify = rows_all_v.flatten()         // QUEUE: Map rows
    def rows_for_join     = rows_all_v.flatten()

    //------------------------------------------------------------------
    // 1) Build output header once
    //------------------------------------------------------------------
    def new_cols   = ['computed_genotype']
    def header_v   = mapping_ch.map { m -> (m.orig_headers + new_cols) as List }   // VALUE
    header_v.view { "Header (raw): $it" }

    // Pass rename map if provided
    def rename_v   = mapping_ch.map { m -> (m.containsKey('rename') && m.rename instanceof Map) ? m.rename : [:] }

    //------------------------------------------------------------------
    // 2) Prepare tuples to classify: (pid, snp_file)
    //------------------------------------------------------------------
    def to_classify_ch = rows_for_classify.map { row ->
      def snp_file = file(row.genotype_file_path)
      if (!snp_file.exists())
        error "Genotype file not found for PID ${row.participant_id}: ${row.genotype_file_path}"

      tuple(
        row.participant_id as String,
        snp_file
      )
    }
    to_classify_ch.view { "To classify: $it" }

    //------------------------------------------------------------------
    // 3) Run classifier
    //------------------------------------------------------------------
    def apol1 = apol1_classifier(assets_dir_ch, to_classify_ch)

    //------------------------------------------------------------------
    // 4) Extract (pid, computed_genotype) from apol1_<PID>.genotype.txt
    //------------------------------------------------------------------
    def apol1_calls = apol1.genotype_line.map { f ->
      def text = f.text?.trim() ?: ''
      def m = text =~ /^APOL1 genotype:\s*(\S+)/
      def computed = m ? m[0][1] : 'UNKNOWN'
      def pidm = f.name =~ /^apol1_(.+)\.genotype\.txt$/
      def pid = pidm ? pidm[0][1] : 'UNKNOWN_PID'
      tuple(pid, computed)
    }
    apol1_calls.view { "APOL1 call: $it" }

    //------------------------------------------------------------------
    // 5) Left-join: keep all input rows, attach computed_genotype
    //------------------------------------------------------------------
    def all_pids_ch = rows_for_join.map { r -> tuple(r.participant_id as String, r) }  // QUEUE

    def joined_rows_ch =
      all_pids_ch
        .map { pid, row -> tuple(pid, row, null) }
        .mix( apol1_calls.map { pid, computed -> tuple(pid, null, computed) } )
        .groupTuple(by: 0)
        .map { pid, rows, computed_list ->
          def row = rows.find { it != null }
          def computed = computed_list.find { it != null } ?: 'UNKNOWN'
          def out = new LinkedHashMap(row)
          out.computed_genotype = computed
          out
        }
    joined_rows_ch.view { "Joined row: $it" }

    //------------------------------------------------------------------
    // 6) Collect maps for Python writer
    //------------------------------------------------------------------
    def rows_maps_v = joined_rows_ch.collect()     // VALUE: List<Map>
    rows_maps_v.view { rows -> "Rows collected for writer: ${rows?.size() ?: 0}" }

    //------------------------------------------------------------------
    // 7) Detect base filename & delimiter from original sheet
    //------------------------------------------------------------------
    def sample_sheet_filename_v = sample_sheet_ch.map { f ->
      def n = f.name; def i = n.lastIndexOf('.'); i != -1 ? n[0..<i] : n
    } // VALUE

    def extension_v = sample_sheet_ch.map { f ->
      def firstLine = f.readLines().head()
      firstLine.contains('\t') ? 'tsv' : 'csv'
    } // VALUE

    def date_stamp   = java.time.LocalDate.now().format(java.time.format.DateTimeFormatter.ISO_DATE)
    def date_stamp_v = Channel.value(date_stamp)

    //------------------------------------------------------------------
    // 8) Write augmented sheet (Python)
    //------------------------------------------------------------------
    def csv_out    = write_sheet_py(header_v, rows_maps_v, rename_v, sample_sheet_filename_v, extension_v, date_stamp_v)

  emit:
    augmented_csv = csv_out.augmented_csv
    apol1_rsids   = apol1.rsid_table
    apol1_genos   = apol1.genotype_line
    apol1_evid    = apol1.evidence
}

process apol1_classifier {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  tag { pid }
  debug true

  input:
    path assets_dir
    tuple val(pid), path(snp_file)

  output:
    path "apol1_${pid}.rsid.tsv",     emit: rsid_table
    path "apol1_${pid}.genotype.txt", emit: genotype_line
    path "apol1_${pid}.evidence.txt", emit: evidence

  script:
  """
  echo "Processing PID: ${pid}, SNP file: \$(basename "${snp_file}")"
  set -euo pipefail
  APOL1_OUT_PREFIX="apol1_${pid}" \\
    python3 "${assets_dir}/apol1_classifier.py" < "${snp_file}"
  """
}

/*
 * Python-based sheet writer. No bash string gymnastics.
 * Inputs are VALUEs: header_list (List), rows_maps (List<Map>), rename (Map),
 *                    sample_sheet_filename (String), extension (String), date_stamp (String).
 */
process write_sheet_py {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  debug true

  input:
    val header_list
    val rows_maps
    val rename_map
    val sample_sheet_filename
    val extension
    val date_stamp

  output:
    path "${sample_sheet_filename}-${date_stamp}.${extension}", emit: augmented_csv

  script:
  // ---- SANITIZE to plain strings to avoid Groovy JSON recursion/stack overflows ----
  def header_clean = (header_list instanceof List) ? header_list.collect { it == null ? "" : it.toString().trim() } : []
  def rename_clean = (rename_map instanceof Map) ? rename_map.collectEntries { k,v ->
    [(k == null ? "" : k.toString().trim()) : (v == null ? "" : v.toString().trim())]
  } : [:]
  def rows_clean = (rows_maps instanceof List) ? rows_maps.collect { r ->
    (r instanceof Map) ? r.collectEntries { k,v ->
      [(k == null ? "" : k.toString().trim()) : (v == null ? "" : v.toString())]
    } : [:]
  } : []

  def header_json = groovy.json.JsonOutput.toJson(header_clean)
  def rows_json   = groovy.json.JsonOutput.toJson(rows_clean)
  def rename_json = groovy.json.JsonOutput.toJson(rename_clean)

  """
  set -euo pipefail

  # Write JSON payloads (plain strings only)
  cat > header.json <<EOF
  ${header_json}
EOF
  cat > rows.json <<EOF
  ${rows_json}
EOF
  cat > rename.json <<EOF
  ${rename_json}
EOF

  python3 - <<'PY'
import json, csv, sys

# Load inputs
with open('header.json','r') as fh:
    header = json.load(fh)
with open('rows.json','r') as fh:
    rows = json.load(fh)
with open('rename.json','r') as fh:
    rename = json.load(fh)

# Normalize header: strings, trimmed
header_norm = [ (h or '') for h in header ]
header_norm = [ str(h).strip() for h in header_norm ]

# Build a mapping from desired header -> actual row key (case-insensitive)
# Apply rename map FIRST if provided (logical header), THEN match case-insensitively in rows.
def build_keymap(sample_row):
    idx = { str(k).strip().lower(): k for k in sample_row.keys() }
    km = {}
    for h in header_norm:
        logical = rename.get(h, h) if isinstance(rename, dict) else h
        key = idx.get(str(logical).strip().lower())
        km[h] = key  # may be None
    return km

extension = '${extension}'
delimiter = '\\t' if extension.strip().lower() == 'tsv' else ','
out_path = f"{'${sample_sheet_filename}'}-{'${date_stamp}'}.{extension}"

with open(out_path, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=delimiter, quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header_norm)

    if not rows:
        print("DEBUG: No rows to write", file=sys.stderr)
    else:
        km = build_keymap(rows[0])

        # DEBUG
        print("DEBUG header_norm:", header_norm, file=sys.stderr)
        print("DEBUG rename_map:", rename, file=sys.stderr)
        print("DEBUG sample row keys:", list(rows[0].keys()), file=sys.stderr)
        print("DEBUG header->row key map:", km, file=sys.stderr)

        for r in rows:
            vals = []
            for h in header_norm:
                rk = km.get(h)
                v = r.get(rk, '') if rk is not None else ''
                vals.append('' if v is None else str(v))
            writer.writerow(vals)

        # DEBUG first 2 rendered rows
        preview = []
        for r in rows[:2]:
            preview.append([ '' if (r.get(km.get(h,''), '') is None) else str(r.get(km.get(h,''), ''))
                             for h in header_norm ])
        print("DEBUG first 2 rendered rows:", preview, file=sys.stderr)
PY
  """
}
