nextflow.enable.dsl=2

// participant id
params.participant_id = params.participant_id ?: null

// snp file
params.snp = params.snp ?: null

// assets directory
params.assets_dir = params.assets_dir ?: null
// results directory
params.results_dir = params.results_dir ?: null
// workflow file
params.work_flow_file = params.work_flow_file ?: null

participant_id_ch    = Channel.value(params.participant_id)
snp_ch               = Channel.fromPath(params.snp)
assets_dir_ch        = Channel.fromPath(params.assets_dir)
work_flow_file       = Channel.fromPath(params.work_flow_file)

// EXAMPLE
// nextflow run template.nf \
//     --participant_id NA07357 \
//     --snp ./genome_23andMe_v4_Full.txt \
//     --work_flow_file ./workflow.nf \
//     --assets_dir ./assets \
//     --results_dir ./results \
//     -with-docker

include { USER } from "${params.work_flow_file}"

workflow {
  USER(
    participant_id_ch,
    snp_ch,
    assets_dir_ch,
    params.results_dir
  )
}