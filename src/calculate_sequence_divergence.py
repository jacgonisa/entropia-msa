#!/usr/bin/env python3
"""
Sequence Divergence Calculator for Multiple Sequence Alignments

This script calculates taxon-level divergence metrics for each sequence in an MSA,
complementing the column-wise Shannon entropy analysis with row-wise divergence analysis.

Usage:
    python calculate_sequence_divergence.py [--mode MODE]

Modes:
    - pairwise: Mean pairwise distance to all other sequences
    - consensus: Distance to consensus sequence
    - both: Calculate both metrics (default)

Output:
    - sequence_divergence_results.csv: Per-sequence divergence metrics
    - divergence_plot.png: Visualization of sequence divergence
"""

import numpy as np
from collections import Counter
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple, Dict
import glob
import argparse


def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Read a FASTA format alignment file.

    Args:
        filepath: Path to FASTA file

    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def trim_alignment(sequences: Dict[str, str], gap_threshold: float = 0.8) -> Dict[str, str]:
    """
    Trim alignment columns with gaps in more than gap_threshold of sequences.

    Args:
        sequences: Dictionary of sequence IDs to sequences
        gap_threshold: Fraction of sequences that must have gaps to remove column

    Returns:
        Dictionary of trimmed sequences
    """
    if not sequences:
        return {}

    seq_list = list(sequences.values())
    seq_ids = list(sequences.keys())
    num_seqs = len(seq_list)
    alignment_length = len(seq_list[0])

    # Identify columns to keep
    columns_to_keep = []
    for col_idx in range(alignment_length):
        gap_count = sum(1 for seq in seq_list if seq[col_idx] == '-')
        gap_fraction = gap_count / num_seqs
        if gap_fraction <= gap_threshold:
            columns_to_keep.append(col_idx)

    # Build trimmed sequences
    trimmed_seqs = {}
    for seq_id, seq in zip(seq_ids, seq_list):
        trimmed_seq = ''.join(seq[i] for i in columns_to_keep)
        trimmed_seqs[seq_id] = trimmed_seq

    return trimmed_seqs


def calculate_consensus(sequences: Dict[str, str]) -> str:
    """
    Calculate consensus sequence (most common amino acid at each position).

    Args:
        sequences: Dictionary of sequence IDs to sequences

    Returns:
        Consensus sequence string
    """
    if not sequences:
        return ""

    seq_list = list(sequences.values())
    alignment_length = len(seq_list[0])
    consensus = []

    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in seq_list]
        # Count amino acids (excluding gaps)
        aa_counts = Counter(aa for aa in column if aa != '-')
        if aa_counts:
            most_common = aa_counts.most_common(1)[0][0]
        else:
            most_common = '-'
        consensus.append(most_common)

    return ''.join(consensus)


def calculate_pairwise_distance(seq1: str, seq2: str) -> float:
    """
    Calculate pairwise distance between two sequences (proportion of differences).

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Proportion of differing positions (0-1)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have same length")

    # Only compare positions without gaps in either sequence
    valid_positions = 0
    differences = 0

    for aa1, aa2 in zip(seq1, seq2):
        if aa1 != '-' and aa2 != '-':
            valid_positions += 1
            if aa1 != aa2:
                differences += 1

    if valid_positions == 0:
        return 0.0

    return differences / valid_positions


def calculate_mean_pairwise_distance(seq_id: str, sequences: Dict[str, str]) -> float:
    """
    Calculate mean pairwise distance from one sequence to all others.

    Args:
        seq_id: ID of sequence to analyze
        sequences: Dictionary of all sequences

    Returns:
        Mean pairwise distance (0-1)
    """
    target_seq = sequences[seq_id]
    distances = []

    for other_id, other_seq in sequences.items():
        if other_id != seq_id:
            dist = calculate_pairwise_distance(target_seq, other_seq)
            distances.append(dist)

    return np.mean(distances) if distances else 0.0


def calculate_entropy_contribution(seq_id: str, sequences: Dict[str, str]) -> float:
    """
    Calculate how much a sequence contributes to overall positional entropy.

    Args:
        seq_id: ID of sequence to analyze
        sequences: Dictionary of all sequences

    Returns:
        Mean entropy contribution (higher = more divergent)
    """
    seq_list = list(sequences.values())
    target_seq = sequences[seq_id]
    alignment_length = len(target_seq)

    contributions = []

    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in seq_list]
        target_aa = target_seq[col_idx]

        if target_aa == '-':
            continue

        # Calculate how rare this amino acid is at this position
        aa_counts = Counter(aa for aa in column if aa != '-')
        if not aa_counts:
            continue

        total = sum(aa_counts.values())
        frequency = aa_counts.get(target_aa, 0) / total

        # Contribution is higher when amino acid is rarer
        # Using -log2(frequency) as information content
        if frequency > 0:
            contribution = -np.log2(frequency)
            contributions.append(contribution)

    return np.mean(contributions) if contributions else 0.0


def calculate_sequence_divergence(sequences: Dict[str, str], mode: str = 'both') -> pd.DataFrame:
    """
    Calculate divergence metrics for each sequence.

    Args:
        sequences: Dictionary of sequence IDs to sequences
        mode: 'pairwise', 'consensus', or 'both'

    Returns:
        DataFrame with divergence metrics per sequence
    """
    consensus = calculate_consensus(sequences)
    results = []

    for seq_id, seq in sequences.items():
        metrics = {
            'sequence_id': seq_id,
            'sequence_length': len(seq.replace('-', '')),
            'gap_count': seq.count('-'),
            'gap_fraction': seq.count('-') / len(seq)
        }

        if mode in ['pairwise', 'both']:
            metrics['mean_pairwise_distance'] = calculate_mean_pairwise_distance(seq_id, sequences)
            metrics['entropy_contribution'] = calculate_entropy_contribution(seq_id, sequences)

        if mode in ['consensus', 'both']:
            metrics['distance_to_consensus'] = calculate_pairwise_distance(seq, consensus)

        results.append(metrics)

    return pd.DataFrame(results)


def process_all_alignments(directory: str = '.', mode: str = 'both') -> Dict[str, pd.DataFrame]:
    """
    Process all .msa files and calculate sequence divergence.

    Args:
        directory: Directory containing .msa files
        mode: Divergence calculation mode

    Returns:
        Dictionary mapping gene IDs to divergence DataFrames
    """
    msa_files = sorted(glob.glob(f"{directory}/*.msa"))

    print(f"Found {len(msa_files)} alignment files")
    print("Processing alignments...\n")

    all_results = {}

    for msa_file in msa_files:
        filename = Path(msa_file).name
        gene_id = filename.replace('.msa', '')

        # Read and trim alignment
        sequences = read_fasta(msa_file)
        trimmed_seqs = trim_alignment(sequences, gap_threshold=0.8)

        # Calculate divergence
        divergence_df = calculate_sequence_divergence(trimmed_seqs, mode=mode)
        divergence_df['gene_id'] = gene_id

        all_results[gene_id] = divergence_df

        # Print summary
        mean_div = divergence_df['mean_pairwise_distance'].mean() if 'mean_pairwise_distance' in divergence_df else 0
        print(f"Processed: {gene_id:40s} | Seqs: {len(sequences):3d} | Mean divergence: {mean_div:.4f}")

    return all_results


def plot_divergence_summary(all_results: Dict[str, pd.DataFrame], output_file: str = 'divergence_plot.png'):
    """
    Create visualization of sequence divergence across genes.

    Args:
        all_results: Dictionary of gene IDs to divergence DataFrames
        output_file: Output filename
    """
    # Combine all results
    combined_df = pd.concat(all_results.values(), ignore_index=True)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Distribution of mean pairwise distances
    ax1 = axes[0, 0]
    if 'mean_pairwise_distance' in combined_df.columns:
        combined_df['mean_pairwise_distance'].hist(bins=50, ax=ax1, color='steelblue', edgecolor='black')
        ax1.set_xlabel('Mean Pairwise Distance', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Number of Sequences', fontsize=12, fontweight='bold')
        ax1.set_title('Distribution of Sequence Divergence\n(Mean Pairwise Distance)', fontsize=13, fontweight='bold')
        ax1.axvline(combined_df['mean_pairwise_distance'].median(), color='red', linestyle='--',
                   label=f"Median: {combined_df['mean_pairwise_distance'].median():.3f}", linewidth=2)
        ax1.legend()
        ax1.grid(alpha=0.3)

    # Plot 2: Distance to consensus
    ax2 = axes[0, 1]
    if 'distance_to_consensus' in combined_df.columns:
        combined_df['distance_to_consensus'].hist(bins=50, ax=ax2, color='coral', edgecolor='black')
        ax2.set_xlabel('Distance to Consensus', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of Sequences', fontsize=12, fontweight='bold')
        ax2.set_title('Distribution of Distance to Consensus', fontsize=13, fontweight='bold')
        ax2.axvline(combined_df['distance_to_consensus'].median(), color='darkred', linestyle='--',
                   label=f"Median: {combined_df['distance_to_consensus'].median():.3f}", linewidth=2)
        ax2.legend()
        ax2.grid(alpha=0.3)

    # Plot 3: Per-gene divergence summary
    ax3 = axes[1, 0]
    gene_summary = combined_df.groupby('gene_id')['mean_pairwise_distance'].mean().sort_values(ascending=False).head(20)
    gene_summary.plot(kind='barh', ax=ax3, color='teal')
    ax3.set_xlabel('Mean Pairwise Distance', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Gene ID', fontsize=12, fontweight='bold')
    ax3.set_title('Top 20 Genes by Mean Sequence Divergence', fontsize=13, fontweight='bold')
    ax3.grid(axis='x', alpha=0.3)

    # Plot 4: Entropy contribution
    ax4 = axes[1, 1]
    if 'entropy_contribution' in combined_df.columns:
        combined_df['entropy_contribution'].hist(bins=50, ax=ax4, color='mediumpurple', edgecolor='black')
        ax4.set_xlabel('Entropy Contribution', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Number of Sequences', fontsize=12, fontweight='bold')
        ax4.set_title('Distribution of Entropy Contribution\n(Higher = More Divergent)', fontsize=13, fontweight='bold')
        ax4.axvline(combined_df['entropy_contribution'].median(), color='darkviolet', linestyle='--',
                   label=f"Median: {combined_df['entropy_contribution'].median():.3f}", linewidth=2)
        ax4.legend()
        ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Calculate sequence-level divergence in MSAs')
    parser.add_argument('--mode', choices=['pairwise', 'consensus', 'both'], default='both',
                       help='Divergence calculation mode (default: both)')
    args = parser.parse_args()

    print("=" * 80)
    print("Sequence Divergence Calculator for MSAs")
    print("=" * 80)
    print(f"Mode: {args.mode}")
    print()

    # Process all alignments
    all_results = process_all_alignments('.', mode=args.mode)

    # Combine results
    combined_df = pd.concat(all_results.values(), ignore_index=True)

    # Save results
    output_csv = 'sequence_divergence_results.csv'
    combined_df.to_csv(output_csv, index=False)
    print(f"\n{'=' * 80}")
    print(f"Results saved to: {output_csv}")

    # Display summary
    print(f"\n{'=' * 80}")
    print("Summary Statistics")
    print("=" * 80)
    print(f"Total sequences analyzed: {len(combined_df)}")
    print(f"Total genes: {combined_df['gene_id'].nunique()}")

    if 'mean_pairwise_distance' in combined_df.columns:
        print(f"\nMean pairwise distance:")
        print(f"  Mean:   {combined_df['mean_pairwise_distance'].mean():.4f}")
        print(f"  Median: {combined_df['mean_pairwise_distance'].median():.4f}")
        print(f"  Max:    {combined_df['mean_pairwise_distance'].max():.4f}")

        # Identify most divergent sequences
        print(f"\n{'=' * 80}")
        print("Top 10 Most Divergent Sequences (by mean pairwise distance)")
        print("=" * 80)
        top_divergent = combined_df.nlargest(10, 'mean_pairwise_distance')
        for idx, row in top_divergent.iterrows():
            print(f"{row['sequence_id']:60s} | Gene: {row['gene_id']:20s} | Divergence: {row['mean_pairwise_distance']:.4f}")

    # Create visualization
    print(f"\n{'=' * 80}")
    print("Creating visualization...")
    print("=" * 80)
    plot_divergence_summary(all_results, 'divergence_plot.png')

    print(f"\n{'=' * 80}")
    print("Analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
