#!/usr/bin/env python3
"""
Convert SNP genotype file to plink text format (.ped + .map)

Input format (tab-separated):
  rsid  chromosome  position  genotype  gs  baf  lrr
  rs9701055  1  630053  CC  0.6551  1  0.0887

Output:
  .map file: chromosome, rsid, genetic_distance (0), position
  .ped file: family_id, individual_id, paternal_id (0), maternal_id (0), 
             sex (0=unknown), phenotype (-9=missing), genotypes...
"""

import argparse
import sys
from pathlib import Path


def parse_genotype(genotype_str):
    """
    Convert genotype string (e.g., 'CC', 'AG', 'TT') to plink format.
    Plink uses space-separated alleles: 'CC' -> 'C C', 'AG' -> 'A G'
    Missing data: '--' or '00' -> '0 0'
    """
    if not genotype_str or genotype_str in ['--', '00', 'NA']:
        return '0 0'
    
    # Handle standard diploid genotypes
    if len(genotype_str) == 2:
        return f'{genotype_str[0]} {genotype_str[1]}'
    elif len(genotype_str) == 1:
        # Hemizygous (e.g., X or Y chromosome in males)
        return f'{genotype_str[0]} {genotype_str[0]}'
    else:
        # Invalid format, treat as missing
        return '0 0'


def convert_to_plink(input_file, participant_id, output_prefix):
    """Convert SNP file to plink text format."""
    
    map_entries = []
    genotypes = []
    
    print(f"Reading SNP file: {input_file}", file=sys.stderr)
    
    with open(input_file, 'r') as f:
        # Skip comment lines and find header
        for line in f:
            if line.startswith('#'):
                continue
            if line.strip() and not line.startswith('rsid'):
                # This is data, but no header found
                print("Warning: No header found, assuming standard format", file=sys.stderr)
                break
            if line.startswith('rsid'):
                # Found header, skip it
                break
        
        # Process data lines
        variant_count = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            
            rsid = parts[0]
            chrom = parts[1]
            position = parts[2]
            genotype = parts[3]
            
            # Handle chromosome naming
            # Plink uses 1-22, X, Y, XY, MT
            # Some files use 'chr1', convert to '1'
            chrom = chrom.replace('chr', '')
            
            # Map chromosome 23 to X, 24 to Y, 25 to XY, 26 to MT
            chrom_map = {
                '23': 'X',
                '24': 'Y', 
                '25': 'XY',
                '26': 'MT',
                'M': 'MT'
            }
            chrom = chrom_map.get(chrom, chrom)
            
            # Skip invalid chromosomes
            valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'XY', 'MT']
            if chrom not in valid_chroms:
                continue
            
            # Create .map entry: chrom, rsid, genetic_distance (0), position
            map_entries.append(f"{chrom}\t{rsid}\t0\t{position}")
            
            # Parse and store genotype
            plink_genotype = parse_genotype(genotype)
            genotypes.append(plink_genotype)
            
            variant_count += 1
            if variant_count % 100000 == 0:
                print(f"Processed {variant_count} variants...", file=sys.stderr)
    
    print(f"Total variants processed: {variant_count}", file=sys.stderr)
    
    # Write .map file
    map_file = f"{output_prefix}.map"
    print(f"Writing {map_file}...", file=sys.stderr)
    with open(map_file, 'w') as f:
        for entry in map_entries:
            f.write(entry + '\n')
    
    # Write .ped file
    # Format: family_id individual_id paternal_id maternal_id sex phenotype genotypes...
    ped_file = f"{output_prefix}.ped"
    print(f"Writing {ped_file}...", file=sys.stderr)
    with open(ped_file, 'w') as f:
        # Use participant_id as family and individual ID
        # 0 0 = unknown parents
        # 0 = unknown sex
        # -9 = missing phenotype
        ped_line = f"{participant_id}\t{participant_id}\t0\t0\t0\t-9"
        
        # Add all genotypes
        for gt in genotypes:
            ped_line += f"\t{gt}"
        
        f.write(ped_line + '\n')
    
    print(f"✓ Conversion complete!", file=sys.stderr)
    print(f"  .map file: {map_file} ({len(map_entries)} variants)", file=sys.stderr)
    print(f"  .ped file: {ped_file} (1 sample)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Convert SNP genotype file to plink text format'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input SNP file (tab-separated)'
    )
    parser.add_argument(
        '--participant-id',
        required=True,
        help='Participant/sample ID'
    )
    parser.add_argument(
        '--output-prefix',
        required=True,
        help='Output file prefix (will create .map and .ped files)'
    )
    
    args = parser.parse_args()
    
    try:
        convert_to_plink(args.input, args.participant_id, args.output_prefix)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

