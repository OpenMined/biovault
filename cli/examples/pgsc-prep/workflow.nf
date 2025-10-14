workflow USER {
    take:
      participant_id_ch // example: 000000
      snp_ch            // example: 000000_carika.txt
      assets_dir_ch     // example: /assets
      results_dir       // example: /results/000000
    
    main:
      participant_id_ch.view { "Participant ID: $it" }
      snp_ch.view { "SNP file: $it" }
      assets_dir_ch.view { "Assets Directory: $it" }
      println "Results Directory: ${results_dir}"

      // Step 1: Convert SNP file to plink text format (.ped + .map)
      def plink_text = convert_to_plink_text(
        participant_id_ch,
        snp_ch,
        assets_dir_ch
      )

      // Step 2: Convert to plink binary format using plink
      def plink_binary = convert_to_binary(
        participant_id_ch,
        plink_text.ped,
        plink_text.map,
        results_dir
      )

      // Step 3: Generate instruction file
      generate_instructions(
        participant_id_ch,
        assets_dir_ch,
        results_dir
      )

      plink_binary.msg
        .map { "\n===== Plink Conversion Complete =====\n${it}\n====================================\n" }
        .view()
}

process convert_to_plink_text {
  container 'python:3.11-slim'
  
  input:
  val participant_id
  path snp_file
  path assets_dir

  output:
  path "*.ped", emit: ped
  path "*.map", emit: map

  script:
  def pid = participant_id.toString()
  """
  python3 ./assets/convert_to_plink.py \
    --input ${snp_file} \
    --participant-id '${pid}' \
    --output-prefix '${pid}'
  """
}

process convert_to_binary {
  container 'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1'
  publishDir params.results_dir, mode: 'copy', pattern: "*.{bed,bim,fam}"
  
  input:
  val participant_id
  path ped_file
  path map_file
  val results_dir

  output:
  path "*.bed"
  path "*.bim"
  path "*.fam"
  stdout emit: msg

  script:
  def ped_prefix = ped_file.baseName
  def pid = participant_id.toString()
  """
  set -e
  
  # Convert to binary format
  plink --file '${ped_prefix}' \
        --make-bed \
        --out '${pid}' \
        --allow-extra-chr \
        --output-chr MT
  
  echo "✓ Created plink binary files: ${pid}.bed, ${pid}.bim, ${pid}.fam"
  echo "✓ Total variants: \$(wc -l < ${pid}.bim)"
  echo ""
  echo "Next steps:"
  echo "1. Review the pgsc_calc_instructions.txt file in your results directory"
  echo "2. Run pgsc_calc with the generated plink files"
  """
}

process generate_instructions {
  container 'python:3.11-slim'
  publishDir params.results_dir, mode: 'copy'
  
  input:
  val participant_id
  path assets_dir
  val results_dir

  output:
  path "pgsc_calc_instructions.txt"

  script:
  def pid = participant_id.toString()
  """
  cp ./assets/README.md pgsc_calc_instructions.txt
  
  cat >> pgsc_calc_instructions.txt << 'EOF'

================================================================================
PARTICIPANT-SPECIFIC FILES
================================================================================

Participant ID: ${pid}
Generated files in: ${results_dir}

Files created:
  - ${pid}.bed (binary genotype file)
  - ${pid}.bim (variant information)
  - ${pid}.fam (sample information)

Example pgsc_calc command for this participant:
------------------------------------------------

nextflow run pgscatalog/pgsc_calc \\
  -profile docker \\
  --input ${results_dir}/${pid}.bed \\
  --pgs_id PGS000001 \\
  --outdir ${results_dir}/pgsc_calc_results

Replace PGS000001 with your desired PGS Catalog ID(s).

EOF
  """
}

