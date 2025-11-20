#!/usr/bin/env python3
"""
MSA Overview: Integrated Column-wise and Row-wise Analysis

This script creates a comprehensive overview of an MSA combining:
- Column-wise Shannon entropy (position-level variability)
- Row-wise sequence divergence (taxon-level divergence)

Usage:
    python plot_msa_overview.py <alignment_file.msa>

Output:
    - <gene_id>_msa_overview.png: Integrated visualization
    - <gene_id>_msa_overview.pdf: High-quality PDF version
"""

import numpy as np
from collections import Counter
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from typing import List, Tuple, Dict
import sys


def read_fasta(filepath: str) -> Dict[str, str]:
    """Read FASTA alignment file."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def trim_alignment(sequences: Dict[str, str], gap_threshold: float = 0.8) -> Tuple[Dict[str, str], List[int]]:
    """Trim alignment and return kept column indices."""
    if not sequences:
        return {}, []

    seq_list = list(sequences.values())
    seq_ids = list(sequences.keys())
    num_seqs = len(seq_list)
    alignment_length = len(seq_list[0])

    columns_to_keep = []
    for col_idx in range(alignment_length):
        gap_count = sum(1 for seq in seq_list if seq[col_idx] == '-')
        if gap_count / num_seqs <= gap_threshold:
            columns_to_keep.append(col_idx)

    trimmed_seqs = {}
    for seq_id, seq in zip(seq_ids, seq_list):
        trimmed_seqs[seq_id] = ''.join(seq[i] for i in columns_to_keep)

    return trimmed_seqs, columns_to_keep


def calculate_positional_entropy(sequences: Dict[str, str]) -> np.ndarray:
    """Calculate Shannon entropy for each position."""
    if not sequences:
        return np.array([])

    seq_list = list(sequences.values())
    alignment_length = len(seq_list[0])
    max_entropy = np.log2(20)
    entropies = []

    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in seq_list]
        aa_counts = Counter(aa for aa in column if aa != '-')

        if not aa_counts:
            entropies.append(0.0)
            continue

        total = sum(aa_counts.values())
        entropy = -sum((count/total) * np.log2(count/total) for count in aa_counts.values())
        entropies.append(entropy / max_entropy)

    return np.array(entropies)


def calculate_consensus(sequences: Dict[str, str]) -> str:
    """Calculate consensus sequence."""
    if not sequences:
        return ""

    seq_list = list(sequences.values())
    alignment_length = len(seq_list[0])
    consensus = []

    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in seq_list]
        aa_counts = Counter(aa for aa in column if aa != '-')
        consensus.append(aa_counts.most_common(1)[0][0] if aa_counts else '-')

    return ''.join(consensus)


def calculate_pairwise_distance(seq1: str, seq2: str) -> float:
    """Calculate proportion of differences between sequences."""
    valid_positions = 0
    differences = 0

    for aa1, aa2 in zip(seq1, seq2):
        if aa1 != '-' and aa2 != '-':
            valid_positions += 1
            if aa1 != aa2:
                differences += 1

    return differences / valid_positions if valid_positions > 0 else 0.0


def calculate_sequence_divergence(sequences: Dict[str, str]) -> pd.DataFrame:
    """Calculate divergence metrics for each sequence."""
    consensus = calculate_consensus(sequences)
    results = []

    for seq_id, seq in sequences.items():
        # Calculate mean pairwise distance
        distances = [calculate_pairwise_distance(seq, other_seq)
                    for other_id, other_seq in sequences.items() if other_id != seq_id]
        mean_pw_dist = np.mean(distances) if distances else 0.0

        # Distance to consensus
        dist_to_cons = calculate_pairwise_distance(seq, consensus)

        results.append({
            'sequence_id': seq_id,
            'mean_pairwise_distance': mean_pw_dist,
            'distance_to_consensus': dist_to_cons,
            'gap_fraction': seq.count('-') / len(seq)
        })

    return pd.DataFrame(results)


def plot_msa_overview(alignment_file: str):
    """
    Create comprehensive MSA overview with column and row-wise analyses.
    """
    gene_id = Path(alignment_file).stem

    print("=" * 80)
    print(f"MSA Overview Generator: {gene_id}")
    print("=" * 80)

    # Read and process alignment
    print("\nReading alignment...")
    sequences = read_fasta(alignment_file)
    print(f"  Sequences: {len(sequences)}")
    print(f"  Length: {len(list(sequences.values())[0])}")

    print("\nTrimming alignment (>80% gaps)...")
    trimmed_seqs, kept_cols = trim_alignment(sequences)
    print(f"  Trimmed length: {len(list(trimmed_seqs.values())[0])}")

    # Calculate metrics
    print("\nCalculating column-wise entropy...")
    entropy = calculate_positional_entropy(trimmed_seqs)

    print("Calculating row-wise divergence...")
    divergence_df = calculate_sequence_divergence(trimmed_seqs)
    divergence_df = divergence_df.sort_values('mean_pairwise_distance', ascending=False)

    # Create comprehensive figure
    print("\nCreating visualization...")
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[0.25, 0.5, 0.25],
                  width_ratios=[0.7, 0.3], hspace=0.3, wspace=0.3)

    # 1. Top: Positional entropy line plot
    ax_entropy = fig.add_subplot(gs[0, :])
    positions = np.arange(1, len(entropy) + 1)
    ax_entropy.plot(positions, entropy, color='steelblue', linewidth=1.5, alpha=0.7)
    ax_entropy.fill_between(positions, entropy, alpha=0.3, color='steelblue')
    ax_entropy.axhline(entropy.mean(), color='red', linestyle='--', linewidth=2,
                      label=f'Mean: {entropy.mean():.3f}', alpha=0.7)
    ax_entropy.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax_entropy.set_ylabel('Normalized\nShannon Entropy', fontsize=12, fontweight='bold')
    ax_entropy.set_title(f'Column-wise Analysis: Positional Entropy\n{gene_id}',
                        fontsize=14, fontweight='bold')
    ax_entropy.legend(loc='upper right', fontsize=10)
    ax_entropy.grid(alpha=0.3, linestyle='--')
    ax_entropy.set_ylim(0, 1)
    ax_entropy.axhspan(0, 0.2, alpha=0.08, color='green', zorder=0)
    ax_entropy.axhspan(0.2, 0.5, alpha=0.08, color='yellow', zorder=0)
    ax_entropy.axhspan(0.5, 1.0, alpha=0.08, color='red', zorder=0)

    # 2. Middle left: Heatmap of alignment (sample)
    ax_heatmap = fig.add_subplot(gs[1, 0])

    # Sample sequences for visualization (max 50)
    seq_ids = list(trimmed_seqs.keys())
    if len(seq_ids) > 50:
        # Sample: top 25 divergent + bottom 25 conserved
        sample_ids = (divergence_df.head(25)['sequence_id'].tolist() +
                     divergence_df.tail(25)['sequence_id'].tolist())
    else:
        sample_ids = divergence_df['sequence_id'].tolist()

    # Create numeric matrix for heatmap (amino acids as numbers)
    aa_to_num = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY-')}
    matrix = []
    for seq_id in sample_ids:
        seq = trimmed_seqs[seq_id]
        # Sample positions if too long
        if len(seq) > 500:
            step = len(seq) // 500
            seq_sample = seq[::step]
        else:
            seq_sample = seq
        matrix.append([aa_to_num.get(aa, 20) for aa in seq_sample])

    matrix = np.array(matrix)

    # Plot heatmap
    im = ax_heatmap.imshow(matrix, aspect='auto', cmap='tab20', interpolation='nearest')
    ax_heatmap.set_ylabel('Sequences (sorted by divergence)', fontsize=12, fontweight='bold')
    ax_heatmap.set_xlabel('Alignment Position (sampled)', fontsize=12, fontweight='bold')
    ax_heatmap.set_title('Alignment Heatmap\n(Top/Bottom divergent sequences)',
                        fontsize=13, fontweight='bold')
    ax_heatmap.set_yticks([])

    # 3. Middle right: Sequence divergence barplot
    ax_divbar = fig.add_subplot(gs[1, 1])
    top_n = min(30, len(divergence_df))
    top_div = divergence_df.head(top_n)

    colors = plt.cm.YlOrRd(np.linspace(0.3, 0.9, top_n))
    ax_divbar.barh(range(top_n), top_div['mean_pairwise_distance'], color=colors, edgecolor='black', linewidth=0.5)
    ax_divbar.set_yticks(range(top_n))
    ax_divbar.set_yticklabels([sid.split('_')[-1][:20] for sid in top_div['sequence_id']], fontsize=8)
    ax_divbar.set_xlabel('Mean Pairwise\nDistance', fontsize=11, fontweight='bold')
    ax_divbar.set_ylabel('')
    ax_divbar.set_title(f'Row-wise Analysis:\nTop {top_n} Divergent Sequences',
                       fontsize=13, fontweight='bold')
    ax_divbar.invert_yaxis()
    ax_divbar.grid(axis='x', alpha=0.3, linestyle='--')

    # Add value labels
    for i, (idx, row) in enumerate(top_div.iterrows()):
        ax_divbar.text(row['mean_pairwise_distance'] + 0.01, i,
                      f"{row['mean_pairwise_distance']:.3f}",
                      va='center', fontsize=7)

    # 4. Bottom left: Entropy distribution
    ax_entdist = fig.add_subplot(gs[2, 0])
    ax_entdist.hist(entropy, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax_entdist.axvline(entropy.mean(), color='red', linestyle='--', linewidth=2,
                      label=f'Mean: {entropy.mean():.3f}')
    ax_entdist.axvline(np.median(entropy), color='orange', linestyle='--', linewidth=2,
                      label=f'Median: {np.median(entropy):.3f}')
    ax_entdist.set_xlabel('Normalized Shannon Entropy', fontsize=11, fontweight='bold')
    ax_entdist.set_ylabel('Number of Positions', fontsize=11, fontweight='bold')
    ax_entdist.set_title('Distribution of Positional Entropy', fontsize=12, fontweight='bold')
    ax_entdist.legend(fontsize=9)
    ax_entdist.grid(axis='y', alpha=0.3)

    # 5. Bottom right: Divergence distribution
    ax_divdist = fig.add_subplot(gs[2, 1])
    ax_divdist.hist(divergence_df['mean_pairwise_distance'], bins=30,
                   color='coral', edgecolor='black', alpha=0.7)
    ax_divdist.axvline(divergence_df['mean_pairwise_distance'].mean(), color='darkred',
                      linestyle='--', linewidth=2,
                      label=f"Mean: {divergence_df['mean_pairwise_distance'].mean():.3f}")
    ax_divdist.axvline(divergence_df['mean_pairwise_distance'].median(), color='orange',
                      linestyle='--', linewidth=2,
                      label=f"Median: {divergence_df['mean_pairwise_distance'].median():.3f}")
    ax_divdist.set_xlabel('Mean Pairwise Distance', fontsize=11, fontweight='bold')
    ax_divdist.set_ylabel('Number of Sequences', fontsize=11, fontweight='bold')
    ax_divdist.set_title('Distribution of Sequence Divergence', fontsize=12, fontweight='bold')
    ax_divdist.legend(fontsize=9)
    ax_divdist.grid(axis='y', alpha=0.3)

    # Main title
    fig.suptitle(f'Complete MSA Overview: {gene_id}\n'
                f'Column-wise Entropy + Row-wise Divergence Analysis',
                fontsize=16, fontweight='bold', y=0.98)

    # Save outputs
    output_png = f"{gene_id}_msa_overview.png"
    output_pdf = f"{gene_id}_msa_overview.pdf"

    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')

    print(f"\nOutputs saved:")
    print(f"  {output_png}")
    print(f"  {output_pdf}")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("Summary Statistics")
    print("=" * 80)
    print(f"\nColumn-wise (Positional Entropy):")
    print(f"  Mean entropy: {entropy.mean():.4f}")
    print(f"  Median entropy: {np.median(entropy):.4f}")
    print(f"  Highly variable positions (>0.5): {(entropy > 0.5).sum()} / {len(entropy)}")

    print(f"\nRow-wise (Sequence Divergence):")
    print(f"  Mean divergence: {divergence_df['mean_pairwise_distance'].mean():.4f}")
    print(f"  Median divergence: {divergence_df['mean_pairwise_distance'].median():.4f}")
    print(f"  Highly divergent sequences (>0.3): {(divergence_df['mean_pairwise_distance'] > 0.3).sum()} / {len(divergence_df)}")

    print(f"\nMost divergent sequences:")
    for idx, row in divergence_df.head(5).iterrows():
        print(f"  {row['sequence_id']:60s} | Divergence: {row['mean_pairwise_distance']:.4f}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


def main():
    """Main execution function."""
    if len(sys.argv) != 2:
        print("Usage: python plot_msa_overview.py <alignment_file.msa>")
        sys.exit(1)

    alignment_file = sys.argv[1]

    if not Path(alignment_file).exists():
        print(f"Error: File '{alignment_file}' not found")
        sys.exit(1)

    plot_msa_overview(alignment_file)


if __name__ == '__main__':
    main()
