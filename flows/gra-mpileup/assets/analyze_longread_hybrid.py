#!/usr/bin/env python3
"""
Long-Read Hybrid Gene Detection for CYP11B1/CYP11B2.

This script analyzes long reads that span the gene region to detect
hybrid (chimeric) molecules. A hybrid read will show B1 sequence at 5'
and B2 sequence at 3' within the same read.

Key insight: Long reads can span multiple discriminating bases, allowing
us to see the transition point within a single molecule.
"""

import argparse
import sys
import numpy as np
from collections import defaultdict

# GRCh38 coordinates
CYP11B1 = {'start': 142876120, 'end': 142879816}
CYP11B2 = {'start': 142914143, 'end': 142917843}

# Key discriminating positions with B1 and B2 bases (GRCh38)
# Selected positions with high discrimination power
DISCRIMINATING = [
    # (position, B1_base, B2_base, region_label)
    # These are positions within CYP11B1 where we can tell B1 from B2
    (142879700, 'G', 'A', 'E1'),
    (142879150, 'A', 'G', 'E2'),
    (142879080, 'G', 'A', 'E2'),
    (142878900, 'C', 'T', 'I2'),  # Common breakpoint region
    (142878800, 'A', 'G', 'I2'),  # Common breakpoint region
    (142878700, 'G', 'A', 'I2'),  # Common breakpoint region
    (142878650, 'T', 'C', 'E3'),
    (142878400, 'G', 'A', 'I3'),
    (142878200, 'C', 'T', 'E4'),
    (142877950, 'G', 'A', 'I4'),
    (142877800, 'T', 'C', 'E5'),
    (142877500, 'G', 'A', 'I5'),
    (142877100, 'C', 'T', 'E6'),
    (142876800, 'A', 'G', 'E7/8'),
    (142876500, 'G', 'A', 'E9'),
]

def parse_args():
    parser = argparse.ArgumentParser(description='Long-read hybrid detection')
    parser.add_argument('--bam', required=True, help='Long-read BAM file')
    parser.add_argument('--output-pdf', required=True, help='Output PDF')
    parser.add_argument('--output-png', required=True, help='Output PNG')
    parser.add_argument('--output-tsv', help='Output TSV with per-read results')
    parser.add_argument('--participant', default='Unknown', help='Participant ID')
    parser.add_argument('--min-disc', type=int, default=3, help='Min discriminating positions per read')
    return parser.parse_args()

def classify_read_at_positions(read_seq, read_start, read_cigar_tuples, reference_positions):
    """
    For a long read, determine which gene (B1 or B2) it supports at each discriminating position.
    Returns list of (position, 'B1'|'B2'|'other', region_label)
    """
    results = []

    # Build mapping from reference position to read position using CIGAR
    ref_to_read = {}
    ref_pos = read_start
    read_pos = 0

    for op, length in read_cigar_tuples:
        if op in [0, 7, 8]:  # M, =, X (alignment match)
            for i in range(length):
                ref_to_read[ref_pos + i] = read_pos + i
            ref_pos += length
            read_pos += length
        elif op == 1:  # I (insertion)
            read_pos += length
        elif op == 2:  # D (deletion)
            ref_pos += length
        elif op == 3:  # N (skip)
            ref_pos += length
        elif op == 4:  # S (soft clip)
            read_pos += length
        elif op == 5:  # H (hard clip)
            pass

    for disc_pos, b1_base, b2_base, region in DISCRIMINATING:
        if disc_pos in ref_to_read:
            read_idx = ref_to_read[disc_pos]
            if 0 <= read_idx < len(read_seq):
                base = read_seq[read_idx].upper()
                if base == b1_base.upper():
                    results.append((disc_pos, 'B1', region))
                elif base == b2_base.upper():
                    results.append((disc_pos, 'B2', region))
                else:
                    results.append((disc_pos, 'other', region))

    return results

def analyze_reads(bam_path, min_discriminating=3):
    """
    Analyze all reads in BAM to classify each as B1, B2, hybrid, or ambiguous.
    """
    try:
        import pysam
    except ImportError:
        print("pysam not available", file=sys.stderr)
        return None

    read_classifications = []

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"Could not open BAM: {e}", file=sys.stderr)
        return None

    region_start = CYP11B1['start'] - 1000
    region_end = CYP11B1['end'] + 1000

    for read in bam.fetch(region=f"chr8:{region_start}-{region_end}"):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if not read.cigartuples:
            continue

        positions = classify_read_at_positions(
            read.query_sequence,
            read.reference_start,
            read.cigartuples,
            DISCRIMINATING
        )

        if len(positions) < min_discriminating:
            continue

        # Classify read
        b1_count = sum(1 for _, call, _ in positions if call == 'B1')
        b2_count = sum(1 for _, call, _ in positions if call == 'B2')
        total = b1_count + b2_count

        if total == 0:
            classification = 'ambiguous'
        elif b1_count == total:
            classification = 'pure_B1'
        elif b2_count == total:
            classification = 'pure_B2'
        elif b1_count > 0 and b2_count > 0:
            # Check if there's a transition (hybrid pattern)
            # Hybrid: B1 at 5' (high positions) → B2 at 3' (low positions)
            b1_positions = [p for p, c, _ in positions if c == 'B1']
            b2_positions = [p for p, c, _ in positions if c == 'B2']

            if min(b1_positions) > max(b2_positions):
                classification = 'hybrid_B1_to_B2'  # Classic GRA pattern
            elif min(b2_positions) > max(b1_positions):
                classification = 'hybrid_B2_to_B1'  # Reverse pattern
            else:
                classification = 'mixed'
        else:
            classification = 'ambiguous'

        read_classifications.append({
            'read_name': read.query_name,
            'read_length': read.query_length,
            'positions': positions,
            'b1_count': b1_count,
            'b2_count': b2_count,
            'classification': classification
        })

    bam.close()
    return read_classifications

def plot_results(classifications, output_pdf, output_png, participant):
    """Create visualization of hybrid detection results."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import matplotlib.gridspec as gridspec

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'pdf.fonttype': 42,
    })

    # Count classifications
    class_counts = defaultdict(int)
    for r in classifications:
        class_counts[r['classification']] += 1

    # Create figure
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[1, 2])

    # Panel 1: Classification pie chart
    ax_pie = fig.add_subplot(gs[0, 0])
    labels = []
    sizes = []
    colors = []
    color_map = {
        'pure_B1': '#2171B5',
        'pure_B2': '#D62728',
        'hybrid_B1_to_B2': '#2CA02C',
        'hybrid_B2_to_B1': '#9467BD',
        'mixed': '#FF7F0E',
        'ambiguous': '#AAAAAA'
    }

    for cls in ['pure_B1', 'pure_B2', 'hybrid_B1_to_B2', 'hybrid_B2_to_B1', 'mixed', 'ambiguous']:
        if class_counts[cls] > 0:
            labels.append(cls.replace('_', ' ').title())
            sizes.append(class_counts[cls])
            colors.append(color_map.get(cls, '#AAAAAA'))

    if sizes:
        ax_pie.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax_pie.set_title('Read Classifications', fontweight='bold')

    # Panel 2: Summary stats
    ax_stats = fig.add_subplot(gs[0, 1])
    ax_stats.axis('off')

    total_reads = len(classifications)
    hybrid_count = class_counts['hybrid_B1_to_B2'] + class_counts['hybrid_B2_to_B1']

    stats_text = f"""Long-Read Hybrid Analysis — {participant}

Total reads analyzed: {total_reads}
├── Pure CYP11B1: {class_counts['pure_B1']} ({100*class_counts['pure_B1']/total_reads:.1f}%)
├── Pure CYP11B2: {class_counts['pure_B2']} ({100*class_counts['pure_B2']/total_reads:.1f}%)
├── Hybrid (B1→B2): {class_counts['hybrid_B1_to_B2']} ({100*class_counts['hybrid_B1_to_B2']/total_reads:.1f}%) ← GRA pattern
├── Hybrid (B2→B1): {class_counts['hybrid_B2_to_B1']} ({100*class_counts['hybrid_B2_to_B1']/total_reads:.1f}%)
├── Mixed: {class_counts['mixed']} ({100*class_counts['mixed']/total_reads:.1f}%)
└── Ambiguous: {class_counts['ambiguous']} ({100*class_counts['ambiguous']/total_reads:.1f}%)

INTERPRETATION:
• "Hybrid B1→B2" reads have B1 at 5' and B2 at 3' = FH-I/GRA signature
• These reads contain the chimeric gene breakpoint
• Ratio of hybrid/pure reads indicates zygosity"""

    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes,
                  fontsize=10, va='top', ha='left', family='monospace',
                  bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.9))

    # Panel 3: Per-read heatmap showing B1/B2 pattern across positions
    ax_heatmap = fig.add_subplot(gs[1, :])

    # Sort reads by classification then by first B2 position
    def sort_key(r):
        cls_order = {'pure_B1': 0, 'hybrid_B1_to_B2': 1, 'mixed': 2, 'hybrid_B2_to_B1': 3, 'pure_B2': 4, 'ambiguous': 5}
        b2_pos = [p for p, c, _ in r['positions'] if c == 'B2']
        first_b2 = max(b2_pos) if b2_pos else 0
        return (cls_order.get(r['classification'], 5), -first_b2)

    sorted_reads = sorted(classifications, key=sort_key)[:100]  # Show top 100

    # Create matrix
    disc_positions = [p for p, _, _ in DISCRIMINATING]
    disc_positions_sorted = sorted(disc_positions, reverse=True)  # 5' to 3'

    matrix = np.zeros((len(sorted_reads), len(disc_positions_sorted)))

    for i, read in enumerate(sorted_reads):
        pos_dict = {p: c for p, c, _ in read['positions']}
        for j, pos in enumerate(disc_positions_sorted):
            if pos in pos_dict:
                if pos_dict[pos] == 'B1':
                    matrix[i, j] = 1.0
                elif pos_dict[pos] == 'B2':
                    matrix[i, j] = -1.0

    # Custom colormap: Blue (B1) -> White (unknown) -> Red (B2)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('b1b2', ['#D62728', 'white', '#2171B5'])

    im = ax_heatmap.imshow(matrix, aspect='auto', cmap=cmap, vmin=-1, vmax=1)

    # Labels
    region_labels = [r for _, _, r in DISCRIMINATING]
    region_labels_sorted = [region_labels[disc_positions.index(p)] for p in disc_positions_sorted]

    ax_heatmap.set_xticks(range(len(disc_positions_sorted)))
    ax_heatmap.set_xticklabels([f"{p/1e6:.4f}\n({r})" for p, r in zip(disc_positions_sorted, region_labels_sorted)],
                               rotation=45, ha='right', fontsize=7)
    ax_heatmap.set_xlabel("Discriminating positions (5' → 3')", fontsize=10)
    ax_heatmap.set_ylabel("Individual reads", fontsize=10)
    ax_heatmap.set_title("Per-Read Gene Origin at Discriminating Positions\n(Blue=B1, Red=B2, White=unknown)",
                         fontsize=11, fontweight='bold')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heatmap, orientation='vertical', pad=0.02, shrink=0.8)
    cbar.set_ticks([-1, 0, 1])
    cbar.set_ticklabels(['CYP11B2', 'Unknown', 'CYP11B1'])

    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(output_png, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close()

def main():
    args = parse_args()

    print(f"Analyzing long reads from {args.bam}...", file=sys.stderr)

    classifications = analyze_reads(args.bam, args.min_disc)

    if classifications is None or len(classifications) == 0:
        print("No reads could be analyzed. Check if pysam is installed and BAM is valid.", file=sys.stderr)
        sys.exit(1)

    print(f"Analyzed {len(classifications)} reads", file=sys.stderr)

    # Write TSV if requested
    if args.output_tsv:
        with open(args.output_tsv, 'w') as f:
            f.write("read_name\tread_length\tb1_count\tb2_count\tclassification\n")
            for r in classifications:
                f.write(f"{r['read_name']}\t{r['read_length']}\t{r['b1_count']}\t{r['b2_count']}\t{r['classification']}\n")

    plot_results(classifications, args.output_pdf, args.output_png, args.participant)
    print(f"Generated plots", file=sys.stderr)

if __name__ == '__main__':
    main()
