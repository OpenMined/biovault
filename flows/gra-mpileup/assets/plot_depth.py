#!/usr/bin/env python3
"""
Plot read depth for CYP11B1/CYP11B2 region analysis (GRA detection).
Generates a publication-quality coverage plot with highlighted duplication event region.
"""

import argparse
import sys
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Plot read depth for GRA analysis')
    parser.add_argument('--input', required=True, help='Input depth.txt file')
    parser.add_argument('--output-pdf', required=True, help='Output PDF file')
    parser.add_argument('--output-png', required=True, help='Output PNG file')
    parser.add_argument('--region', required=True, help='Genomic region (e.g., chr8:142869000-143030000)')
    parser.add_argument('--participant', default='Unknown', help='Participant ID')
    return parser.parse_args()

def parse_region(region_str):
    chrom, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def load_depth(filename):
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

def setup_publication_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'lines.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def format_position(x, pos):
    return f'{x/1e6:.2f}'

def main():
    args = parse_args()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, MultipleLocator
    except ImportError as e:
        print(f"Error importing matplotlib: {e}", file=sys.stderr)
        sys.exit(1)

    setup_publication_style()
    chrom, region_start, region_end = parse_region(args.region)
    positions, depths = load_depth(args.input)

    if len(positions) == 0:
        print("No depth data found, creating empty plot", file=sys.stderr)
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.text(0.5, 0.5, f'No coverage data for {args.region}',
                ha='center', va='center', fontsize=10, transform=ax.transAxes)
        ax.set_xlabel(f'Position on {chrom} (Mb)')
        ax.set_ylabel('Read Depth')
        plt.savefig(args.output_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1)
        plt.savefig(args.output_png, format='png', bbox_inches='tight', pad_inches=0.1)
        plt.close()
        return

    mean_depth = np.mean(depths)
    median_depth = np.median(depths)
    dup_start = region_start + int((region_end - region_start) * 0.12)
    dup_end = region_start + int((region_end - region_start) * 0.75)

    fig, ax = plt.subplots(figsize=(7, 3.5))

    ax.axvspan(dup_start, dup_end, alpha=0.15, color='#FFD700',
               label='Potential duplication region', zorder=1)

    window = max(1, len(positions) // 500)
    if window > 1:
        smoothed = np.convolve(depths, np.ones(window)/window, mode='valid')
        smooth_pos = positions[window//2:window//2+len(smoothed)]
        ax.fill_between(smooth_pos, 0, smoothed, alpha=0.3, color='#4A90D9', zorder=2)
        ax.plot(smooth_pos, smoothed, color='#2171B5', linewidth=0.8, zorder=3)
    else:
        ax.fill_between(positions, 0, depths, alpha=0.3, color='#4A90D9', zorder=2)
        ax.plot(positions, depths, color='#2171B5', linewidth=0.8, zorder=3)

    ax.axhline(y=mean_depth, color='#D62728', linestyle='--', linewidth=1.0,
               alpha=0.8, label=f'Mean depth: {mean_depth:.0f}×', zorder=4)

    dup_mask = (positions >= dup_start) & (positions <= dup_end)
    flank_mask = ~dup_mask
    if np.any(dup_mask) and np.any(flank_mask):
        dup_mean = np.mean(depths[dup_mask])
        flank_mean = np.mean(depths[flank_mask])
        ratio = dup_mean / flank_mean if flank_mean > 0 else 0
        ax.text(0.98, 0.95, f'Dup/Flank: {ratio:.2f}×',
                transform=ax.transAxes, fontsize=9, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor='#CCCCCC', linewidth=0.5))

    ax.set_xlabel(f'Position on {chrom} (Mb)')
    ax.set_ylabel('Read Depth (×)')
    ax.set_title(f'CYP11B1/CYP11B2 Coverage — {args.participant}', fontweight='medium')

    ax.xaxis.set_major_formatter(FuncFormatter(format_position))
    ax.set_xlim(positions.min(), positions.max())
    ax.set_ylim(0, max(depths) * 1.1)

    ax.yaxis.set_minor_locator(MultipleLocator(mean_depth / 4))
    ax.tick_params(which='minor', length=2)

    legend = ax.legend(loc='upper left', frameon=True, framealpha=0.95,
                       edgecolor='#CCCCCC', fancybox=False)
    legend.get_frame().set_linewidth(0.5)

    plt.tight_layout()
    plt.savefig(args.output_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(args.output_png, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close()

    print(f"Coverage analysis complete for {args.participant}", file=sys.stderr)
    print(f"  Region: {args.region}", file=sys.stderr)
    print(f"  Mean depth: {mean_depth:.1f}x", file=sys.stderr)
    print(f"  Median depth: {median_depth:.1f}x", file=sys.stderr)
    print(f"  Data points: {len(positions)}", file=sys.stderr)

if __name__ == '__main__':
    main()
