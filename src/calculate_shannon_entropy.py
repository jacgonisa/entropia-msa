#!/usr/bin/env python3
"""
Shannon Entropy Calculator for Multiple Sequence Alignments

This script calculates normalized Shannon entropy for amino acid alignments.
Alignments are trimmed to remove columns with >80% gaps before entropy calculation.

Usage:
    python calculate_shannon_entropy.py

Output:
    - shannon_entropy_results.csv: Detailed per-gene entropy statistics
    - entropy_plot.png: Visualization of genes ranked by mean normalized entropy
"""

import numpy as np
from collections import Counter
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple, Dict
import glob


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


def trim_alignment(sequences: Dict[str, str], gap_threshold: float = 0.8) -> List[str]:
    """
    Trim alignment columns with gaps in more than gap_threshold of sequences.

    Args:
        sequences: Dictionary of sequence IDs to sequences
        gap_threshold: Fraction of sequences that must have gaps to remove column (default: 0.8)

    Returns:
        List of trimmed sequences
    """
    if not sequences:
        return []

    # Convert to list for easier indexing
    seq_list = list(sequences.values())
    num_seqs = len(seq_list)
    alignment_length = len(seq_list[0])

    # Identify columns to keep
    columns_to_keep = []
    for col_idx in range(alignment_length):
        # Count gaps in this column
        gap_count = sum(1 for seq in seq_list if seq[col_idx] == '-')
        gap_fraction = gap_count / num_seqs

        # Keep column if gap fraction is below threshold
        if gap_fraction <= gap_threshold:
            columns_to_keep.append(col_idx)

    # Build trimmed sequences
    trimmed_seqs = []
    for seq in seq_list:
        trimmed_seq = ''.join(seq[i] for i in columns_to_keep)
        trimmed_seqs.append(trimmed_seq)

    return trimmed_seqs


def calculate_shannon_entropy(sequences: List[str], normalize: bool = True) -> Tuple[np.ndarray, float]:
    """
    Calculate Shannon entropy for each position in the alignment.

    For amino acids: H = -Σ(pi × log2(pi))
    Normalized by dividing by log2(20) = 4.321928... for 20 standard amino acids

    Args:
        sequences: List of aligned sequences (after trimming)
        normalize: Whether to normalize by maximum possible entropy (default: True)

    Returns:
        Tuple of (entropy array for each position, mean entropy)
    """
    if not sequences or len(sequences[0]) == 0:
        return np.array([]), 0.0

    alignment_length = len(sequences[0])
    num_seqs = len(sequences)

    # Maximum possible entropy for amino acids (20 standard amino acids)
    max_entropy = np.log2(20)

    entropies = []

    for col_idx in range(alignment_length):
        # Get all amino acids at this position
        column = [seq[col_idx] for seq in sequences]

        # Count amino acids (excluding gaps)
        aa_counts = Counter(aa for aa in column if aa != '-')

        if not aa_counts:
            # All gaps - entropy is 0
            entropies.append(0.0)
            continue

        # Calculate Shannon entropy
        total = sum(aa_counts.values())
        entropy = 0.0

        for count in aa_counts.values():
            if count > 0:
                p_i = count / total
                entropy -= p_i * np.log2(p_i)

        # Normalize if requested
        if normalize:
            entropy = entropy / max_entropy

        entropies.append(entropy)

    entropy_array = np.array(entropies)
    mean_entropy = np.mean(entropy_array) if len(entropy_array) > 0 else 0.0

    return entropy_array, mean_entropy


def process_all_alignments(directory: str = '.') -> pd.DataFrame:
    """
    Process all .msa files in the directory and calculate entropy statistics.

    Args:
        directory: Directory containing .msa files

    Returns:
        DataFrame with entropy statistics for each gene
    """
    results = []

    # Find all .msa files
    msa_files = sorted(glob.glob(f"{directory}/*.msa"))

    print(f"Found {len(msa_files)} alignment files")
    print("Processing alignments...\n")

    for msa_file in msa_files:
        filename = Path(msa_file).name
        # Extract gene ID and name from filename (format: OG0000103_H3.msa)
        gene_id = filename.replace('.msa', '')

        # Read alignment
        sequences = read_fasta(msa_file)
        num_sequences = len(sequences)
        original_length = len(list(sequences.values())[0]) if sequences else 0

        # Trim alignment
        trimmed_seqs = trim_alignment(sequences, gap_threshold=0.8)
        trimmed_length = len(trimmed_seqs[0]) if trimmed_seqs else 0

        # Calculate entropy
        entropy_array, mean_entropy = calculate_shannon_entropy(trimmed_seqs, normalize=True)

        # Calculate additional statistics
        median_entropy = np.median(entropy_array) if len(entropy_array) > 0 else 0.0
        std_entropy = np.std(entropy_array) if len(entropy_array) > 0 else 0.0
        max_entropy = np.max(entropy_array) if len(entropy_array) > 0 else 0.0
        min_entropy = np.min(entropy_array) if len(entropy_array) > 0 else 0.0

        results.append({
            'gene_id': gene_id,
            'num_sequences': num_sequences,
            'original_length': original_length,
            'trimmed_length': trimmed_length,
            'positions_removed': original_length - trimmed_length,
            'percent_removed': 100 * (original_length - trimmed_length) / original_length if original_length > 0 else 0,
            'mean_normalized_entropy': mean_entropy,
            'median_normalized_entropy': median_entropy,
            'std_normalized_entropy': std_entropy,
            'max_normalized_entropy': max_entropy,
            'min_normalized_entropy': min_entropy
        })

        print(f"Processed: {gene_id:40s} | Seqs: {num_sequences:3d} | "
              f"Length: {original_length:5d} -> {trimmed_length:5d} | "
              f"Mean entropy: {mean_entropy:.4f}")

    return pd.DataFrame(results)


def plot_entropy_results(df: pd.DataFrame, output_file: str = 'entropy_plot.png'):
    """
    Create visualization showing genes ranked by normalized Shannon entropy.

    Args:
        df: DataFrame with entropy results
        output_file: Output filename for plot
    """
    # Sort by mean normalized entropy
    df_sorted = df.sort_values('mean_normalized_entropy', ascending=False)

    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 1, figsize=(16, 12))

    # Plot 1: Bar plot of top genes by mean entropy
    top_n = min(50, len(df_sorted))
    top_genes = df_sorted.head(top_n)

    colors = plt.cm.viridis(np.linspace(0, 1, top_n))

    ax1 = axes[0]
    bars = ax1.barh(range(top_n), top_genes['mean_normalized_entropy'], color=colors)
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(top_genes['gene_id'], fontsize=8)
    ax1.set_xlabel('Mean Normalized Shannon Entropy', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Gene ID', fontsize=12, fontweight='bold')
    ax1.set_title(f'Top {top_n} Genes by Mean Normalized Shannon Entropy\n'
                  f'(Normalized by log₂(20) = 4.32 for 20 amino acids)',
                  fontsize=14, fontweight='bold', pad=20)
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 1.0)

    # Add value labels on bars
    for i, (idx, row) in enumerate(top_genes.iterrows()):
        ax1.text(row['mean_normalized_entropy'] + 0.01, i,
                f"{row['mean_normalized_entropy']:.3f}",
                va='center', fontsize=7)

    # Plot 2: Distribution of entropy across all genes
    ax2 = axes[1]

    # Histogram
    n, bins, patches = ax2.hist(df['mean_normalized_entropy'], bins=30,
                                 color='steelblue', alpha=0.7, edgecolor='black')

    # Color the bars by value
    for i, patch in enumerate(patches):
        patch.set_facecolor(plt.cm.viridis(bins[i]))

    ax2.set_xlabel('Mean Normalized Shannon Entropy', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution of Mean Normalized Shannon Entropy Across All Genes',
                  fontsize=14, fontweight='bold', pad=20)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    # Add vertical lines for statistics
    mean_val = df['mean_normalized_entropy'].mean()
    median_val = df['mean_normalized_entropy'].median()

    ax2.axvline(mean_val, color='red', linestyle='--', linewidth=2,
                label=f'Mean: {mean_val:.3f}')
    ax2.axvline(median_val, color='orange', linestyle='--', linewidth=2,
                label=f'Median: {median_val:.3f}')
    ax2.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")

    return fig


def main():
    """Main execution function."""
    print("=" * 80)
    print("Shannon Entropy Calculator for Amino Acid Alignments")
    print("=" * 80)
    print()

    # Process all alignments
    results_df = process_all_alignments('.')

    # Save results to CSV
    output_csv = 'shannon_entropy_results.csv'
    results_df.to_csv(output_csv, index=False)
    print(f"\n{'=' * 80}")
    print(f"Results saved to: {output_csv}")

    # Display summary statistics
    print(f"\n{'=' * 80}")
    print("Summary Statistics")
    print("=" * 80)
    print(f"Total genes analyzed: {len(results_df)}")
    print(f"\nMean normalized entropy statistics:")
    print(f"  Mean:   {results_df['mean_normalized_entropy'].mean():.4f}")
    print(f"  Median: {results_df['mean_normalized_entropy'].median():.4f}")
    print(f"  Std:    {results_df['mean_normalized_entropy'].std():.4f}")
    print(f"  Min:    {results_df['mean_normalized_entropy'].min():.4f}")
    print(f"  Max:    {results_df['mean_normalized_entropy'].max():.4f}")

    # Display top 10 genes with highest entropy
    print(f"\n{'=' * 80}")
    print("Top 10 Genes with Highest Mean Normalized Shannon Entropy")
    print("=" * 80)
    top_10 = results_df.nlargest(10, 'mean_normalized_entropy')
    for idx, row in top_10.iterrows():
        print(f"{row['gene_id']:40s} | Entropy: {row['mean_normalized_entropy']:.4f} | "
              f"Seqs: {row['num_sequences']:3.0f} | Length: {row['trimmed_length']:4.0f}")

    # Create visualization
    print(f"\n{'=' * 80}")
    print("Creating visualization...")
    print("=" * 80)
    plot_entropy_results(results_df, 'entropy_plot.png')

    print(f"\n{'=' * 80}")
    print("Analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
