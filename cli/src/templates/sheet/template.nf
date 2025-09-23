nextflow.enable.dsl=2

// samplesheet file (supports csv and tsv)
params.samplesheet    = params.samplesheet    ?: 'participants.csv'
// schema yaml file
params.schema_yaml  = params.schema_yaml  ?: 'schema.yaml'
// workflow file
params.work_flow_file = params.work_flow_file ?: null
// assets directory
params.assets_dir = params.assets_dir ?: null
// results directory
params.results_dir = params.results_dir ?: null

sample_sheet_ch    = Channel.fromPath(params.samplesheet)
assets_dir_ch    = Channel.fromPath(params.assets_dir)
work_flow_file   = Channel.fromPath(params.work_flow_file)

// nextflow run template.nf \
//     --samplesheet participants.csv \
//     --schema_yaml ./schema.yaml \
//     --work_flow_file ./workflow.nf \
//     --assets_dir ./assets \
//     --results_dir ./results \
//     -with-docker

/*
 * Load YAML schema
 * schema.yaml keys supported:
 *   required: [list of col names]
 *   rename:   { csv_header_name: canonical_name, ... }
 *   defaults: { col: default_value, ... }
 *   types:    { col: string|int|float|bool|path, ... }
 */
def yaml   = new org.yaml.snakeyaml.Yaml()

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

def schema = yaml.load(file(params.schema_yaml).text) ?: [:]

/*
 * Helpers
 */
def castValue = { String v, String t ->
    if (v == null) return null
    switch ((t ?: 'string').toLowerCase()) {
        case 'int'   : return v.isInteger()    ? v.toInteger()   : (v ? v.toBigDecimal().intValue() : null)
        case 'float' : return v.isBigDecimal() ? v.toBigDecimal(): (v ? new BigDecimal(v)           : null)
        case 'bool'  : return ['1','true','yes','y'].contains(v.toLowerCase())
        case 'path'  : return file(v) // NOTE: to stage on disk inside a process, pass again as a `path` input
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

/*
 * Build channels
 *
 * We read the header directly from the samplesheet file (first line), apply renames,
 * and expose it as a value channel. Rows are streamed separately.
 */

// Determine delimiter based on file extension
def delimiter = params.samplesheet.endsWith('.tsv') ? '\t' : ','

// read header line safely
def header_line = file(params.samplesheet).newInputStream().withReader { it.readLine() }
if (!header_line) {
    throw new IllegalArgumentException("Samplesheet file appears empty or unreadable: ${params.samplesheet}")
}

def orig_headers_raw = header_line.split(delimiter).collect { it.trim() }
def rename_map = (schema.rename ?: [:])
def renamed_headers = orig_headers_raw.collect { h -> rename_map.getOrDefault(h, h) }

// mapping meta we pass inward (extend as needed)
mapping_ch = Channel.value([
    orig_headers: renamed_headers,         // original column order (after rename)
    types       : (schema.types ?: [:]),
    required    : (schema.required ?: []),
    rename      : rename_map,
    defaults    : (schema.defaults ?: [:])
])

// rows stream (maps), tidy keys then apply schema
rows_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true, sep: delimiter)
    .map { Map row ->
        // tidy keys first
        def tidied = row.collectEntries { k,v -> [ (k?.trim()), v ] }
        applySchema(tidied)
    }

/*
 * Include and run your inner workflow
 * Implement whatever logic you want in workflow.nf, consuming:
 *   - rows_ch      (channel of row maps)
 *   - mapping_ch   (single map with headers/types/etc.)
 *   - params.results_dir
 */
include { USER } from "${params.work_flow_file}"

workflow {
    // ensure results dir exists (harmless if already present)
    file(params.results_dir).mkdirs()

    // hand off to inner workflow
    USER(
        sample_sheet_ch,
        rows_ch,
        mapping_ch,
        assets_dir_ch,
        params.results_dir
    )
}
