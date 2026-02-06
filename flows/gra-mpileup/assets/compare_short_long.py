#!/usr/bin/env python3
"""
Compare short-read vs long-read coverage for CYP11B1/CYP11B2 region.
Useful for validating hybrid gene detection across platforms.

Usage:
  python compare_short_long.py \
    --short-depth short_depth.txt \
    --long-depth long_depth.txt \
    --output-pdf comparison.pdf \
    --output-png comparison.png \
    --participant "SAMPLE_ID"
"""

import argparse
import sys
import numpy as np

CYP11B1 = {'name': 'CYP11B1', 'start': 142876120, 'end': 142879816}
CYP11B2 = {'name': 'CYP11B2', 'start': 142914143, 'end': 142917843}

def parse_args():
    parser = argparse.ArgumentParser(description='Compare short vs long read coverage')
    parser.add_argument('--short-depth', required=True, help='Short-read depth.txt')
    parser.add_argument('--long-depth', required=True, help='Long-read depth.txt')
    parser.add_argument('--output-pdf', required=True, help='Output PDF')
    parser.add_argument('--output-png', required=True, help='Output PNG')
    parser.add_argument('--participant', default='Unknown', help='Participant ID')
    return parser.parse_args()

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
        print(f"Error reading {filename}: {e}", file=sys.stderr)
        return np.array([]), np.array([])
    return np.array(positions), np.array(depths)

def smooth(positions, depths, window=100):
    if len(depths) < window:
        return positions, depths
    smoothed = np.convolve(depths, np.ones(window)/window, mode='valid')
    smooth_pos = positions[window//2:window//2+len(smoothed)]
    return smooth_pos, smoothed

def calculate_stats(positions, depths, gene):
    mask = (positions >= gene['start']) & (positions <= gene['end'])
    if np.any(mask):
        return np.mean(depths[mask])
    return 0

def main():
    args = parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'pdf.fonttype': 42,
    })

    # Load data
    short_pos, short_depth = load_depth(args.short_depth)
    long_pos, long_depth = load_depth(args.long_depth)

    if len(short_pos) == 0 and len(long_pos) == 0:
        print("No data in either file", file=sys.stderr)
        sys.exit(1)

    # Smooth
    short_pos_s, short_depth_s = smooth(short_pos, short_depth)
    long_pos_s, long_depth_s = smooth(long_pos, long_depth)

    # Create figure
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 0.5], hspace=0.3)

    # Panel 1: Short reads
    ax_short = fig.add_subplot(gs[0, :])
    if len(short_pos_s) > 0:
        ax_short.fill_between(short_pos_s, 0, short_depth_s, alpha=0.3, color='#4A90D9')
        ax_short.plot(short_pos_s, short_depth_s, color='#2171B5', linewidth=0.6)
        mean_short = np.mean(short_depth)
        ax_short.axhline(y=mean_short, color='#D62728', linestyle='--', linewidth=0.8,
                        label=f'Mean: {mean_short:.0f}×')

    ax_short.axvspan(CYP11B1['start'], CYP11B1['end'], alpha=0.1, color='blue')
    ax_short.axvspan(CYP11B2['start'], CYP11B2['end'], alpha=0.1, color='red')
    ax_short.set_ylabel('Short-Read Depth (×)')
    ax_short.set_title(f'Short-Read Coverage — {args.participant}', fontweight='bold')
    ax_short.legend(loc='upper right', fontsize=8)

    # Panel 2: Long reads
    ax_long = fig.add_subplot(gs[1, :], sharex=ax_short)
    if len(long_pos_s) > 0:
        ax_long.fill_between(long_pos_s, 0, long_depth_s, alpha=0.3, color='#90D94A')
        ax_long.plot(long_pos_s, long_depth_s, color='#217B25', linewidth=0.6)
        mean_long = np.mean(long_depth)
        ax_long.axhline(y=mean_long, color='#D62728', linestyle='--', linewidth=0.8,
                       label=f'Mean: {mean_long:.0f}×')

    ax_long.axvspan(CYP11B1['start'], CYP11B1['end'], alpha=0.1, color='blue', label='CYP11B1')
    ax_long.axvspan(CYP11B2['start'], CYP11B2['end'], alpha=0.1, color='red', label='CYP11B2')
    ax_long.set_ylabel('Long-Read Depth (×)')
    ax_long.set_xlabel('Position on chr8 (Mb)')
    ax_long.set_title('Long-Read Coverage', fontweight='bold')
    ax_long.legend(loc='upper right', fontsize=8)

    # Format x-axis
    def format_mb(x, pos):
        return f'{x/1e6:.2f}'
    from matplotlib.ticker import FuncFormatter
    ax_long.xaxis.set_major_formatter(FuncFormatter(format_mb))

    # Panel 3: Comparison stats
    ax_stats = fig.add_subplot(gs[2, :])
    ax_stats.axis('off')

    # Calculate per-gene stats
    short_b1 = calculate_stats(short_pos, short_depth, CYP11B1)
    short_b2 = calculate_stats(short_pos, short_depth, CYP11B2)
    long_b1 = calculate_stats(long_pos, long_depth, CYP11B1)
    long_b2 = calculate_stats(long_pos, long_depth, CYP11B2)

    stats_text = f"""
┌─────────────────────┬────────────────┬────────────────┬─────────────┐
│ Metric              │ Short Reads    │ Long Reads     │ Difference  │
├─────────────────────┼────────────────┼────────────────┼─────────────┤
│ CYP11B1 mean depth  │ {short_b1:>10.1f}×   │ {long_b1:>10.1f}×   │ {long_b1-short_b1:>+9.1f}×  │
│ CYP11B2 mean depth  │ {short_b2:>10.1f}×   │ {long_b2:>10.1f}×   │ {long_b2-short_b2:>+9.1f}×  │
│ B1/B2 ratio         │ {short_b1/short_b2 if short_b2>0 else 0:>10.2f}    │ {long_b1/long_b2 if long_b2>0 else 0:>10.2f}    │             │
│ Total mean depth    │ {np.mean(short_depth) if len(short_depth)>0 else 0:>10.1f}×   │ {np.mean(long_depth) if len(long_depth)>0 else 0:>10.1f}×   │             │
└─────────────────────┴────────────────┴────────────────┴─────────────┘

INTERPRETATION:
• B1/B2 ratio ~1.0 suggests both genes present (normal or heterozygous hybrid)
• B1/B2 ratio >1.5 may indicate B1 duplication or hybrid gene extra copy
• Long reads provide better mapping in the high-similarity intergenic region
"""

    ax_stats.text(0.5, 0.5, stats_text, transform=ax_stats.transAxes,
                  fontsize=9, va='center', ha='center', family='monospace',
                  bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.9))

    plt.tight_layout()
    plt.savefig(args.output_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(args.output_png, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close()

    print("Generated comparison plot", file=sys.stderr)

if __name__ == '__main__':
    main()
