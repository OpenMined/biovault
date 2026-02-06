#!/usr/bin/env python3
"""
Plot read depth for CYP11B1/CYP11B2 region with gene structure annotation.
Shows exon positions and common GRA breakpoint regions.
"""

import argparse
import sys
import numpy as np

# GRCh38 coordinates for CYP11B1 and CYP11B2
# Source: NCBI Gene database
CYP11B1 = {
    'name': 'CYP11B1',
    'strand': '-',
    'start': 142876120,
    'end': 142879816,
    'exons': [
        (142879612, 142879816),  # Exon 1
        (142879066, 142879205),  # Exon 2
        (142878580, 142878716),  # Exon 3
        (142878169, 142878335),  # Exon 4
        (142877749, 142877879),  # Exon 5
        (142877407, 142877543),  # Exon 6
        (142877036, 142877208),  # Exon 7
        (142876613, 142876843),  # Exon 8
        (142876120, 142876406),  # Exon 9
    ]
}

CYP11B2 = {
    'name': 'CYP11B2',
    'strand': '-',
    'start': 142914143,
    'end': 142917843,
    'exons': [
        (142917639, 142917843),  # Exon 1
        (142917093, 142917232),  # Exon 2
        (142916607, 142916743),  # Exon 3
        (142916196, 142916362),  # Exon 4
        (142915776, 142915906),  # Exon 5
        (142915434, 142915570),  # Exon 6
        (142915063, 142915235),  # Exon 7
        (142914640, 142914870),  # Exon 8
        (142914143, 142914429),  # Exon 9
    ]
}

# Common breakpoint regions from GRAde paper (Table 1)
# Most common: E2-I2 (53%), then E3-I3, E4-I4, E5-I5
BREAKPOINT_LABELS = {
    'E2-I2': 'Intron 2 (most common ~53%)',
    'E3-I3': 'Exon 3 - Intron 3',
    'E4-I4': 'Exon 4 - Intron 4',
    'E5-I5': 'Exon 5 - Intron 5',
}

def parse_args():
    parser = argparse.ArgumentParser(description='Plot annotated read depth for GRA analysis')
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

def setup_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
        'font.size': 9,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 7,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'axes.linewidth': 0.8,
        'pdf.fonttype': 42,
    })

def format_position_mb(x, pos):
    return f'{x/1e6:.3f}'

def draw_gene(ax, gene, y_center, height, color, label_offset=0):
    """Draw gene structure with exons as boxes."""
    from matplotlib.patches import Rectangle, FancyBboxPatch

    # Draw gene body (thin line)
    ax.hlines(y_center, gene['start'], gene['end'], colors=color, linewidth=1.5, alpha=0.5)

    # Draw exons (thick boxes)
    for i, (exon_start, exon_end) in enumerate(gene['exons']):
        rect = Rectangle((exon_start, y_center - height/2),
                         exon_end - exon_start, height,
                         facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.8)
        ax.add_patch(rect)

    # Label gene
    ax.text(gene['start'] - 500, y_center + label_offset, gene['name'],
            fontsize=9, fontweight='bold', color=color, va='center', ha='right')

    # Arrow showing strand direction
    arrow_x = gene['end'] + 200 if gene['strand'] == '+' else gene['start'] - 200
    arrow_dx = 300 if gene['strand'] == '+' else -300
    ax.annotate('', xy=(arrow_x + arrow_dx, y_center), xytext=(arrow_x, y_center),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

def main():
    args = parse_args()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, MultipleLocator
        from matplotlib.patches import Rectangle
        import matplotlib.gridspec as gridspec
    except ImportError as e:
        print(f"Error importing matplotlib: {e}", file=sys.stderr)
        sys.exit(1)

    setup_style()
    chrom, region_start, region_end = parse_region(args.region)
    positions, depths = load_depth(args.input)

    if len(positions) == 0:
        print("No depth data found", file=sys.stderr)
        sys.exit(1)

    # Create figure with two panels
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 0.3, 2], hspace=0.05)

    # Top panel: Gene structure
    ax_genes = fig.add_subplot(gs[0])
    ax_spacer = fig.add_subplot(gs[1])
    ax_depth = fig.add_subplot(gs[2])

    # Draw genes
    draw_gene(ax_genes, CYP11B1, 0.7, 0.25, '#2171B5')  # Blue for B1
    draw_gene(ax_genes, CYP11B2, 0.3, 0.25, '#D62728')  # Red for B2

    ax_genes.set_xlim(region_start, region_end)
    ax_genes.set_ylim(0, 1)
    ax_genes.set_ylabel('')
    ax_genes.set_title(f'CYP11B1/CYP11B2 Region — {args.participant}', fontweight='bold', fontsize=12)
    ax_genes.set_xticks([])
    ax_genes.set_yticks([])
    ax_genes.spines['top'].set_visible(False)
    ax_genes.spines['right'].set_visible(False)
    ax_genes.spines['bottom'].set_visible(False)
    ax_genes.spines['left'].set_visible(False)

    # Add exon labels
    for i, (es, ee) in enumerate(CYP11B1['exons']):
        ax_genes.text((es+ee)/2, 0.92, f'E{i+1}', fontsize=6, ha='center', color='#2171B5')
    for i, (es, ee) in enumerate(CYP11B2['exons']):
        ax_genes.text((es+ee)/2, 0.08, f'E{i+1}', fontsize=6, ha='center', color='#D62728')

    # Spacer panel - show the intergenic region
    ax_spacer.set_xlim(region_start, region_end)
    ax_spacer.axvspan(CYP11B1['end'], CYP11B2['start'], alpha=0.1, color='gray')
    ax_spacer.text((CYP11B1['end'] + CYP11B2['start'])/2, 0.5,
                   f'Intergenic: {(CYP11B2["start"]-CYP11B1["end"])/1000:.1f} kb',
                   ha='center', va='center', fontsize=8, color='gray')
    ax_spacer.set_yticks([])
    ax_spacer.set_xticks([])
    for spine in ax_spacer.spines.values():
        spine.set_visible(False)

    # Bottom panel: Depth
    mean_depth = np.mean(depths)

    # Smooth the depth data
    window = max(1, len(positions) // 800)
    if window > 1:
        smoothed = np.convolve(depths, np.ones(window)/window, mode='valid')
        smooth_pos = positions[window//2:window//2+len(smoothed)]
    else:
        smoothed = depths
        smooth_pos = positions

    # Fill area
    ax_depth.fill_between(smooth_pos, 0, smoothed, alpha=0.3, color='#4A90D9')
    ax_depth.plot(smooth_pos, smoothed, color='#2171B5', linewidth=0.6)

    # Mean depth line
    ax_depth.axhline(y=mean_depth, color='#D62728', linestyle='--', linewidth=1.0,
                     alpha=0.7, label=f'Mean: {mean_depth:.0f}×')

    # Highlight gene regions
    ax_depth.axvspan(CYP11B1['start'], CYP11B1['end'], alpha=0.08, color='#2171B5', label='CYP11B1')
    ax_depth.axvspan(CYP11B2['start'], CYP11B2['end'], alpha=0.08, color='#D62728', label='CYP11B2')

    # Calculate coverage stats per gene
    b1_mask = (positions >= CYP11B1['start']) & (positions <= CYP11B1['end'])
    b2_mask = (positions >= CYP11B2['start']) & (positions <= CYP11B2['end'])
    intergenic_mask = (positions > CYP11B1['end']) & (positions < CYP11B2['start'])

    if np.any(b1_mask):
        b1_mean = np.mean(depths[b1_mask])
    else:
        b1_mean = 0
    if np.any(b2_mask):
        b2_mean = np.mean(depths[b2_mask])
    else:
        b2_mean = 0
    if np.any(intergenic_mask):
        inter_mean = np.mean(depths[intergenic_mask])
    else:
        inter_mean = mean_depth

    # Stats box
    stats_text = (f'CYP11B1: {b1_mean:.0f}×\n'
                  f'CYP11B2: {b2_mean:.0f}×\n'
                  f'Intergenic: {inter_mean:.0f}×\n'
                  f'B1/B2 ratio: {b1_mean/b2_mean:.2f}' if b2_mean > 0 else '')

    ax_depth.text(0.02, 0.98, stats_text, transform=ax_depth.transAxes,
                  fontsize=8, va='top', ha='left',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                           edgecolor='#CCCCCC', alpha=0.9))

    # Formatting
    ax_depth.set_xlabel(f'Position on {chrom} (Mb)', fontsize=10)
    ax_depth.set_ylabel('Read Depth (×)', fontsize=10)
    ax_depth.xaxis.set_major_formatter(FuncFormatter(format_position_mb))
    ax_depth.set_xlim(region_start, region_end)
    ax_depth.set_ylim(0, max(depths) * 1.1)
    ax_depth.spines['top'].set_visible(False)
    ax_depth.spines['right'].set_visible(False)

    ax_depth.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=7)

    plt.tight_layout()
    plt.savefig(args.output_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(args.output_png, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close()

    print(f"Generated annotated coverage plot", file=sys.stderr)
    print(f"  CYP11B1 mean depth: {b1_mean:.1f}x", file=sys.stderr)
    print(f"  CYP11B2 mean depth: {b2_mean:.1f}x", file=sys.stderr)

if __name__ == '__main__':
    main()
