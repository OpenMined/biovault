nextflow.enable.dsl=2

// ---------- Parameters (override on CLI with --ref, --cram, etc.) ----------
params.ref  = params.ref  ?: null
params.cram = params.cram ?: null

// Set safe defaults without touching undefined params (prevents warnings)
if( !params.containsKey('region') ) params.region = "chr15:28110000-28130000"
if( !params.containsKey('pos')    ) params.pos    = "chr15:28120472-28120472"
if( !params.containsKey('prefix') ) params.prefix = "herc2_rs12913832"



// Ensure required inputs are provided
if( !params.ref || !params.cram ) {
    exit 1, "ERROR: Please provide both --ref <fasta> and --cram <cram> on the command line"
}

// ---------- Workflow ----------
workflow {
    // Stage FASTA + FAI
    Channel
        .fromPath(params.ref)
        .map { fasta -> tuple(fasta, file("${fasta}.fai")) }
        .set { ref_pair_ch }

    // Stage CRAM + CRAI
    Channel
        .fromPath(params.cram)
        .map { cram -> tuple(cram, file("${cram}.crai")) }
        .set { cram_pair_ch }

    // 1) Reference base at POS
    faidx(ref_pair_ch)

    // 2) Depth at POS
    depth(cram_pair_ch)

    // 3) Pileup base counts (A/C/G/T) at POS
    basecounts(cram_pair_ch, ref_pair_ch)

    // 4) Call variants in REGION (emit even 0/0 sites)
    call_region(cram_pair_ch, ref_pair_ch)
        .set { vcf_pair_ch }

    // 5) Interpret rs12913832 and print + write summary
    def out = interpret_eyes(vcf_pair_ch)

    out.msg
      .map { "\n===== Eye Color Interpretation =====\n${it}\n====================================\n" }
      .view()
}

// ---------- Processes ----------

process faidx {
    tag { "faidx:${params.pos}" }
    publishDir "results", mode: 'copy',
        saveAs: { fn -> "${params.prefix}.faidx.txt" }

    input:
    tuple path(ref), path(ref_fai)

    output:
    path "faidx.txt"
    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
    shell:
    '''
    samtools faidx !{ref} !{params.pos} > faidx.txt
    '''
}

process depth {
    tag { "depth:${params.pos}" }
    publishDir "results", mode: 'copy',
        saveAs: { fn -> "${params.prefix}.depth.txt" }

    input:
    tuple path(cram), path(crai)

    output:
    path "depth.txt"

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
    shell:
    '''
    # samtools depth doesn't use FASTA; just ensure CRAI is present
    samtools depth -r !{params.pos} -a -d 0 -Q 0 -G 0 !{cram} > depth.txt
    '''
}

process basecounts {
    tag { "basecounts:${params.pos}" }
    publishDir "results", mode: 'copy',
        saveAs: { fn -> "${params.prefix}.basecounts.txt" }

    input:
    tuple path(cram), path(crai)
    tuple path(ref),  path(ref_fai)

    output:
    path "basecounts.txt"

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'
    shell:
    '''
    set -euo pipefail
    samtools mpileup -f !{ref} -r !{params.pos} !{cram} \
    | awk '{
        b=$5; gsub("\\$", "", b); gsub("\\^.", "", b);
        nA=gsub(/A|a/, "&", b); nC=gsub(/C|c/, "&", b);
        nG=gsub(/G|g/, "&", b); nT=gsub(/T|t/, "&", b);
        print $1, $2, "ref="$3, "depth="$4, "A="nA, "C="nC, "G="nG, "T="nT
    }' > basecounts.txt
    '''
}

process call_region {
    tag { "bcftools:${params.region}" }
    publishDir "results", mode: 'copy',
        saveAs: { fn ->
            fn.endsWith('.vcf.gz.csi') ? "${params.prefix}.vcf.gz.csi" :
            fn.endsWith('.vcf.gz')     ? "${params.prefix}.vcf.gz"     : fn
        }

    input:
    tuple path(cram), path(crai)
    tuple path(ref),  path(ref_fai)

    output:
    tuple path("variants.vcf.gz"), path("variants.vcf.gz.csi")

    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
    shell:
    '''
    # Emit records even for homozygous-reference sites (-A)
    bcftools mpileup -Ou -f !{ref} -r !{params.region} -a FORMAT/AD,FORMAT/DP !{cram} \
      | bcftools call -m -A -Oz -o variants.vcf.gz

    # Use CSI index (robust for large references)
    bcftools index -f -c variants.vcf.gz
    '''
}

process interpret_eyes {
    tag { "interpret:${params.pos}" }
    publishDir "results", mode: 'copy',
        saveAs: { fn -> "${params.prefix}.eye_color.txt" }

    input:
    tuple path(vcf), path(vcf_index)

    output:
    path "eye_color.txt"
    stdout emit: msg   // capture stdout as a channel

    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
    shell:
    """
    set -euo pipefail

    # Query REF, ALT, GT, and AD (allele depths) at the exact site
    line=\$(bcftools query -r !{params.pos} -f "%REF\\t%ALT\\t[%GT]\\t[%AD]\\n" !{vcf} || true)

    if [ -z "\$line" ]; then
      text="The person at CRAM !{params.cram}
Has (unknown;unknown)."
    else
      ref=\$(echo "\$line" | cut -f1)
      alt=\$(echo "\$line" | cut -f2)
      gt=\$(echo  "\$line" | cut -f3)
      ad=\$(echo  "\$line" | cut -f4)

      # split genotype
      a1=\$(echo "\$gt" | sed 's@[|/]@\\t@' | cut -f1)
      a2=\$(echo "\$gt" | sed 's@[|/]@\\t@' | cut -f2)

      get_base () {
        idx="\$1"
        if [ "\$idx" = "0" ]; then
          echo "\$ref"
        elif [ "\$idx" = "." ] || [ -z "\$idx" ]; then
          echo "?"
        else
          echo "\$alt" | awk -v i="\$idx" 'BEGIN{FS=","} {print \$i}'
        fi
      }

      b1=\$(get_base "\$a1")
      b2=\$(get_base "\$a2")

      # Build labeled counts from AD (ref,alt1,alt2,...) and ALT list
      IFS=',' read -ra ad_arr  <<< "\$ad"
      IFS=',' read -ra alt_arr <<< "\$alt"

      count_str="Counts: \$ref: \${ad_arr[0]}"
      if [ \${#ad_arr[@]} -ge 2 ]; then
        for ((i=1; i<\${#ad_arr[@]}; i++)); do
          allele="\${alt_arr[\$((i-1))]:-ALT\$i}"
          count_str="\$count_str, \$allele: \${ad_arr[\$i]}"
        done
      fi

      text=\$(cat <<EOF
The person at CRAM !{params.cram}
Has (\$b1;\$b2).
\$count_str

Interpretation:
(A;A) yields brown eye color ~80% of the time.
(A;G) also tends toward brown.
(G;G) gives blue eye color ~99% of the time.
EOF
)
    fi

    # Write to file AND to stdout
    printf "%s\\n" "\$text" | tee eye_color.txt
    """
}
