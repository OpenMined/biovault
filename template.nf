nextflow.enable.dsl=2

// patient id
params.patient_id = params.patient_id ?: null

// reference genome in fasta format
params.ref  = params.ref  ?: null
// reference genome index file
params.ref_index = params.ref_index ?: null
// reference genome version GRCh37 or GRCh38
params.ref_version = params.ref_version ?: null
// aligned reads in bam, cram etc
params.aligned = params.aligned ?: null
// aligned index file
params.aligned_index = params.aligned_index ?: null
// assets directory
params.assets_dir = params.assets_dir ?: null
// results directory
params.results_dir = params.results_dir ?: null
// workflow file
params.work_flow_file = params.work_flow_file ?: null

patient_id_ch    = Channel.value(params.patient_id)
ref_ch           = Channel.fromPath(params.ref)
ref_index_ch     = Channel.fromPath(params.ref_index)
aligned_ch       = Channel.fromPath(params.aligned)
aligned_index_ch = Channel.fromPath(params.aligned_index)
ref_version_ch   = Channel.value(params.ref_version)
assets_dir_ch    = Channel.fromPath(params.assets_dir)
work_flow_file   = Channel.fromPath(params.work_flow_file)

// EXAMPLE
// nextflow run template.nf \
//     --patient_id NA07357 \
//     --ref_version grch38 \
//     --ref ./data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
//     --ref_index ./data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
//     --aligned ./data/ERR3239283/NA07357.final.cram \
//     --aligned_index ./data/ERR3239283/NA07357.final.cram.crai \
//     --work_flow_file ./workflow.nf \
//     --assets_dir ./assets \
//     --results_dir ./results \
//     -with-docker

include { USER } from "${params.work_flow_file}"

workflow {
  USER(
    patient_id_ch,
    ref_ch,
    ref_index_ch,
    aligned_ch,
    aligned_index_ch,
    params.ref_version,
    assets_dir_ch,
    params.results_dir
  )
}
