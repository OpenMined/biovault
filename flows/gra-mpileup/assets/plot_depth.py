#!/usr/bin/env python3
"""
Plot read depth for CYP11B1/CYP11B2 region analysis (GRA detection).
Generates a coverage plot with highlighted duplication event region.
"""

import argparse
import sys
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Plot read depth for GRA analysis')
    parser.add_argument('--input', required=True, help='Input depth.txt file')
    parser.add_argument('--output', required=True, help='Output PNG file')
    parser.add_argument('--region', required=True, help='Genomic region (e.g., chr8:142869000-143030000)')
    parser.add_argument('--participant', default='Unknown', help='Participant ID')
    return parser.parse_args()

def parse_region(region_str):
    """Parse region string like 'chr8:142869000-143030000' into (chrom, start, end)"""
    chrom, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def load_depth(filename):
    """Load depth data from samtools depth output"""
    positions = []
    depths = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        pos = int(parts[1])
                        depth = int(parts[2])
                        positions.append(pos)
                        depths.append(depth)
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading depth file: {e}", file=sys.stderr)
        return np.array([]), np.array([])

    return np.array(positions), np.array(depths)

def main():
    args = parse_args()

    # Import matplotlib here to handle potential import errors gracefully
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"Error importing matplotlib: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse the region
    chrom, region_start, region_end = parse_region(args.region)

    # Load depth data
    positions, depths = load_depth(args.input)

    if len(positions) == 0:
        print("No depth data found, creating empty plot", file=sys.stderr)
        # Create empty plot with message
        fig, ax = plt.subplots(figsize=(14, 5))
        ax.text(0.5, 0.5, f'No coverage data for {args.region}',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        ax.set_xlabel(f'Genomic Position ({chrom})')
        ax.set_ylabel('Read Depth')
        ax.set_title(f'Coverage Plot for {chrom} - {args.participant}')
        plt.savefig(args.output, dpi=150, bbox_inches='tight')
        plt.close()
        return

    # Calculate statistics
    mean_depth = np.mean(depths)
    median_depth = np.median(depths)

    # Define the duplication event region (CYP11B1/B2 crossover)
    # This is approximately the middle portion where gene conversion/duplication occurs
    # For GRCh38: roughly 142880000-142920000
    dup_start = region_start + int((region_end - region_start) * 0.12)
    dup_end = region_start + int((region_end - region_start) * 0.75)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 5))

    # Highlight duplication event region (yellow background)
    ax.axvspan(dup_start, dup_end, alpha=0.3, color='yellow', label='Potential Duplication Region')

    # Plot coverage line
    ax.plot(positions, depths, color='blue', linewidth=0.5, label='Coverage')

    # Add reference lines
    ax.axhline(y=mean_depth, color='red', linestyle='--', linewidth=0.8, alpha=0.7, label=f'Mean: {mean_depth:.1f}x')

    # Calculate depth in duplication region vs flanking
    dup_mask = (positions >= dup_start) & (positions <= dup_end)
    flank_mask = ~dup_mask

    if np.any(dup_mask) and np.any(flank_mask):
        dup_mean = np.mean(depths[dup_mask])
        flank_mean = np.mean(depths[flank_mask])
        ratio = dup_mean / flank_mean if flank_mean > 0 else 0

        # Add annotation for ratio
        ax.text(0.02, 0.95, f'Dup/Flank ratio: {ratio:.2f}x',
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Labels and title
    ax.set_xlabel(f'Genomic Position ({chrom})', fontsize=11)
    ax.set_ylabel('Read Depth', fontsize=11)
    ax.set_title(f'Extended Coverage Plot for {chrom} with Duplication Event - {args.participant}', fontsize=12)

    # Legend
    ax.legend(loc='upper right', fontsize=9)

    # Format x-axis to show scientific notation
    ax.ticklabel_format(axis='x', style='scientific', scilimits=(0,0))

    # Set axis limits
    ax.set_xlim(positions.min() - 1000, positions.max() + 1000)
    ax.set_ylim(0, max(depths) * 1.1)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Save figure
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    plt.close()

    # Print summary
    print(f"Coverage analysis complete for {args.participant}", file=sys.stderr)
    print(f"  Region: {args.region}", file=sys.stderr)
    print(f"  Mean depth: {mean_depth:.1f}x", file=sys.stderr)
    print(f"  Median depth: {median_depth:.1f}x", file=sys.stderr)
    print(f"  Data points: {len(positions)}", file=sys.stderr)

if __name__ == '__main__':
    main()
