nextflow.enable.dsl=2

workflow USER {
  take:
    sample_sheet_ch // file channel: the original sample sheet
    rows_ch         // channel emitting Map per row (already validated/cast/renamed)
    mapping_ch      // value channel: [orig_headers, types, required, rename, defaults]
    assets_dir_ch   // (unused here)
    results_dir     // directory to write results (use --results_dir ...)

  main:
    // Derived columns to add
    def new_cols = ['bmi','age_group']

    // Header (original-after-rename + derived)
    def header_ch = mapping_ch.map { m -> (m.orig_headers + new_cols) as List }

    // Pair each row with mapping (so we can compute deriveds in header order)
    def paired_ch = rows_ch.combine(mapping_ch)

    // Emit LISTS of fields in header order (do NOT join here)
    def rows_list_ch = paired_ch.map { row, m ->
      // ----- derive new columns -----
      def weight = (row.weight instanceof BigDecimal) ? row.weight :
                   (row.weight != null ? new BigDecimal(row.weight.toString()) : null)

      def height_cm = (row.height instanceof BigDecimal) ? row.height :
                      (row.height != null ? new BigDecimal(row.height.toString()) : null)

      def height_m = (height_cm != null) ? height_cm.divide(new BigDecimal('100')) : null

      def bmi = (weight != null && height_m != null && height_m.compareTo(BigDecimal.ZERO) > 0)
                ? weight.divide(height_m.multiply(height_m), 4, java.math.RoundingMode.HALF_UP)
                : null

      def age = (row.age instanceof Number) ? row.age.intValue() :
                (row.age != null && row.age.toString().isInteger()
                  ? row.age.toString().toInteger()
                  : null)

      def age_group = (age == null) ? null :
                      (age < 18 ? 'child' :
                      (age < 30 ? '18-29' :
                      (age < 45 ? '30-44' :
                      (age < 60 ? '45-59' : '60+'))))

      // Build output map without mutating the original
      def out = new LinkedHashMap(row)
      out.bmi = bmi
      out.age_group = age_group

      // Return values in header order as a LIST
      def header = (m.orig_headers + new_cols) as List
      header.collect { h -> out.containsKey(h) ? out[h] : '' }
    }

    // Collect all rows (List<List<Object>> ideally)
    def body_rows_ch = rows_list_ch.collect()

    // Base filename (no extension)
    def sample_sheet_filename_ch = sample_sheet_ch.map { f ->
      def n = f.name
      def i = n.lastIndexOf('.')
      i != -1 ? n[0..<i] : n
    }

    // Sniff first line of the actual file to decide CSV vs TSV
    def extension_ch = sample_sheet_ch.map { f ->
      def firstLine = f.readLines().head()
      firstLine.contains('\t') ? 'tsv' : 'csv'
    }

    // Write out with the correct delimiter
    write_csv_or_tsv(header_ch, body_rows_ch, sample_sheet_filename_ch, extension_ch)
}

process write_csv_or_tsv {
  container 'alpine:latest'
  shell 'ash'
  publishDir params.results_dir, mode: 'copy'

  input:
    val header_list            // List<String>
    val body_rows              // EITHER List<List<Object>> OR a flat List<Object> (we'll normalize)
    val sample_sheet_filename  // String
    val extension              // 'csv' | 'tsv'

  output:
    path "${sample_sheet_filename}-updated.${extension}", emit: augmented_csv

  script:
  // Choose delimiter once
  def delim = (extension == 'tsv') ? '\t' : ','

  // Chunk a flat list of fields into rows of header-sized chunks
  def chunkRows = { rowsOrFlat, rowWidth ->
    if (!(rowsOrFlat instanceof List) || rowsOrFlat.isEmpty()) return []
    def first = rowsOrFlat[0]
    if (first instanceof List) return rowsOrFlat as List       // already List<List>
    // otherwise it's flat -> collate into rows
    return (rowsOrFlat as List).collate(rowWidth)
  }

  // Join + escape fields for CSV/TSV
  def toDelimited = { fields, d ->
    (fields as List).collect { v ->
      if (v == null) return ''
      def s = v.toString()
      (s.contains(d) || s.contains('"') || s.contains('\n'))
        ? "\"${s.replace('"','""')}\""
        : s
    }.join(d)
  }

  // Normalize body rows shape
  def normalized_rows = chunkRows(body_rows, header_list.size())

  def header = toDelimited(header_list as List, delim)
  def body   = (normalized_rows as List).collect { row -> toDelimited(row as List, delim) }.join('\n')

  """
cat > ${sample_sheet_filename}-updated.${extension} <<'EOF'
${header}
${body}
EOF
  """
}
