nextflow.enable.dsl=2

// ------------------------------------------------------------------
// Params
// ------------------------------------------------------------------

// samplesheet file (supports csv and tsv)
params.samplesheet     = params.samplesheet     ?: 'participants.csv'
// schema yaml file
params.schema_yaml     = params.schema_yaml     ?: 'schema.yaml'
// workflow file
params.work_flow_file  = params.work_flow_file  ?: null
// assets directory
params.assets_dir      = params.assets_dir      ?: null
// results directory
params.results_dir     = params.results_dir     ?: 'results'

// Channels for simple passthrough values/files
sample_sheet_ch = Channel.fromPath(params.samplesheet)
assets_dir_ch = Channel.value( file(params.assets_dir) )
work_flow_file  = Channel.fromPath(params.work_flow_file)

// ------------------------------------------------------------------
// Schema (YAML) loader
// schema.yaml keys supported:
//   required: [list of col names]
//   rename:   { csv_header_name: canonical_name, ... }
//   defaults: { col: default_value, ... }
//   types:    { col: string|int|float|bool|path, ... }
// ------------------------------------------------------------------

def yaml = new org.yaml.snakeyaml.Yaml()

/*
Example schema.yaml
required:
  - participant_id
  - genotype_file_path

rename:
  geno_path: genotype_file_path

defaults:
  age: 40

types:
  participant_id: string
  genotype_file_path: path
  weight: float
  height: float
  age: int
*/

def schema   = yaml.load(file(params.schema_yaml).text) ?: [:]
def sheetDir = file(params.samplesheet).toRealPath().normalize().parent

// ------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------

def castValue = { String v, String t ->
    if (v == null) return null
    switch ((t ?: 'string').toLowerCase()) {
        case 'int'   : return v.isInteger()    ? v.toInteger()   : (v ? v.toBigDecimal().intValue() : null)
        case 'float' : return v.isBigDecimal() ? v.toBigDecimal(): (v ? new BigDecimal(v)           : null)
        case 'bool'  : return ['1','true','yes','y'].contains(v.toLowerCase())
        case 'path'  :
            def raw = v.trim()
            def f = file(raw)
            if (f.isAbsolute() && f.exists()) {
                return f
            }
            def resolved = sheetDir.resolve(raw).normalize()
            return resolved
        default      : return v
    }
}

def applySchema = { Map row ->
    // 1) renames
    (schema.rename ?: [:]).each { from, to ->
        if (row.containsKey(from) && !row.containsKey(to)) {
            row[to] = row.remove(from)
        }
    }
    // 2) defaults
    (schema.defaults ?: [:]).each { k,v ->
        if (!row.containsKey(k) || row[k] == null || row[k].toString().trim() == '') {
            row[k] = v
        }
    }
    // 3) required
    (schema.required ?: []).each { k ->
        if (!row.containsKey(k) || row[k] == null || row[k].toString().trim() == '') {
            throw new IllegalArgumentException("Missing required field: ${k} in row: ${row}")
        }
    }
    // 4) types
    def types = (schema.types ?: [:])
    row.collectEntries { k,v -> [k, castValue(v?.toString(), types[k])] }
}

// ------------------------------------------------------------------
// Build channels
//  - mapping_ch: single value map with header/types/etc
//  - rows_src:   queue channel of row maps (schema-applied)
//  - rows_for_user: duplicated branch for the inner workflow
// ------------------------------------------------------------------

// Determine delimiter based on file extension
def delimiter = params.samplesheet.endsWith('.tsv') ? '\t' : ','

// Read header line safely
def header_line = file(params.samplesheet).newInputStream().withReader { it.readLine() }
if (!header_line) {
    throw new IllegalArgumentException("Samplesheet file appears empty or unreadable: ${params.samplesheet}")
}

def orig_headers_raw = header_line.split(delimiter).collect { it.trim() }
def rename_map       = (schema.rename ?: [:])
def renamed_headers  = orig_headers_raw.collect { h -> rename_map.getOrDefault(h, h) }

// mapping meta we pass inward (extend as needed)
mapping_ch = Channel.value([
    orig_headers: renamed_headers,   // original column order (after rename)
    types       : (schema.types ?: [:]),
    required    : (schema.required ?: []),
    rename      : rename_map,
    defaults    : (schema.defaults ?: [:])
])

// Build the rows source once (as you already have)
rows_src = Channel
  .fromPath(params.samplesheet)
  .splitCsv(header: true, sep: delimiter)
  .map { Map row ->
    def tidied = row.collectEntries { k,v -> [(k?.trim()), v] }
    applySchema(tidied)
  }

// ❌ remove this (it causes the DataflowBroadcast.into error):
// rows_src.into { rows_for_user /* , rows_for_stats, rows_for_logging, ... */ }

// ✅ just bind it (one consumer == no duplication needed)
def rows_for_user = rows_src

// ------------------------------------------------------------------
// Include and run inner workflow
// The inner workflow (workflow.nf) should define `workflow USER { ... }`
// and accept: sample_sheet_ch, rows_for_user, mapping_ch, assets_dir_ch, results_dir
// ------------------------------------------------------------------

include { USER } from "${params.work_flow_file}"

workflow {
    // ensure results dir exists (harmless if already present)
    file(params.results_dir).mkdirs()

    // hand off to inner workflow
    USER(
        sample_sheet_ch,
        rows_for_user,
        mapping_ch,
        assets_dir_ch,
        params.results_dir
    )
}
