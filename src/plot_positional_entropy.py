#!/usr/bin/env python3
"""
Positional Entropy Plotter for Multiple Sequence Alignments

This script generates a comprehensive multi-page PDF showing Shannon entropy
across alignment positions for each gene.

Usage:
    python plot_positional_entropy.py

Output:
    - positional_entropy_all_genes.pdf: Multi-page PDF with one plot per gene
"""

import numpy as np
from collections import Counter
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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


def trim_alignment(sequences: Dict[str, str], gap_threshold: float = 0.8) -> Tuple[List[str], List[int]]:
    """
    Trim alignment columns with gaps in more than gap_threshold of sequences.

    Args:
        sequences: Dictionary of sequence IDs to sequences
        gap_threshold: Fraction of sequences that must have gaps to remove column (default: 0.8)

    Returns:
        Tuple of (trimmed sequences, list of kept column indices)
    """
    if not sequences:
        return [], []

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

    return trimmed_seqs, columns_to_keep


def calculate_shannon_entropy(sequences: List[str], normalize: bool = True) -> np.ndarray:
    """
    Calculate Shannon entropy for each position in the alignment.

    For amino acids: H = -Σ(pi × log2(pi))
    Normalized by dividing by log2(20) = 4.321928... for 20 standard amino acids

    Args:
        sequences: List of aligned sequences (after trimming)
        normalize: Whether to normalize by maximum possible entropy (default: True)

    Returns:
        Array of entropy values for each position
    """
    if not sequences or len(sequences[0]) == 0:
        return np.array([])

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

    return np.array(entropies)


def plot_gene_entropy(gene_id: str, entropy_array: np.ndarray, num_sequences: int,
                      original_length: int, trimmed_length: int, ax: plt.Axes):
    """
    Create a plot of positional entropy for a single gene.

    Args:
        gene_id: Gene identifier
        entropy_array: Array of entropy values for each position
        num_sequences: Number of sequences in alignment
        original_length: Original alignment length before trimming
        trimmed_length: Alignment length after trimming
        ax: Matplotlib axes object to plot on
    """
    if len(entropy_array) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax.set_title(f'{gene_id}\nNo data available', fontsize=12, fontweight='bold')
        return

    positions = np.arange(1, len(entropy_array) + 1)

    # Calculate statistics
    mean_entropy = np.mean(entropy_array)
    median_entropy = np.median(entropy_array)
    max_entropy = np.max(entropy_array)
    std_entropy = np.std(entropy_array)

    # Create the plot
    ax.plot(positions, entropy_array, color='steelblue', linewidth=0.8, alpha=0.7)
    ax.fill_between(positions, entropy_array, alpha=0.3, color='steelblue')

    # Add horizontal lines for mean and median
    ax.axhline(mean_entropy, color='red', linestyle='--', linewidth=1.5,
               label=f'Mean: {mean_entropy:.3f}', alpha=0.8)
    ax.axhline(median_entropy, color='orange', linestyle='--', linewidth=1.5,
               label=f'Median: {median_entropy:.3f}', alpha=0.8)

    # Identify high entropy regions (>75th percentile)
    high_entropy_threshold = np.percentile(entropy_array, 75)
    high_entropy_positions = positions[entropy_array > high_entropy_threshold]
    high_entropy_values = entropy_array[entropy_array > high_entropy_threshold]

    # Highlight high entropy positions
    ax.scatter(high_entropy_positions, high_entropy_values, color='red',
               s=20, alpha=0.6, zorder=5, label=f'High variability (>75th percentile)')

    # Labels and formatting
    ax.set_xlabel('Alignment Position', fontsize=11, fontweight='bold')
    ax.set_ylabel('Normalized Shannon Entropy', fontsize=11, fontweight='bold')

    # Title with gene info
    title = f'{gene_id}\n'
    title += f'Sequences: {num_sequences} | Length: {original_length} → {trimmed_length} '
    title += f'({100*(original_length-trimmed_length)/original_length:.1f}% trimmed)'
    ax.set_title(title, fontsize=11, fontweight='bold', pad=15)

    # Add statistics box
    stats_text = f'Mean: {mean_entropy:.4f}\nMedian: {median_entropy:.4f}\n'
    stats_text += f'Max: {max_entropy:.4f}\nStd: {std_entropy:.4f}'
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0, len(entropy_array) + 1)

    # Add subtle background colors for entropy levels
    ax.axhspan(0, 0.2, alpha=0.05, color='green', zorder=0)
    ax.axhspan(0.2, 0.5, alpha=0.05, color='yellow', zorder=0)
    ax.axhspan(0.5, 1.0, alpha=0.05, color='red', zorder=0)


def create_entropy_pdf(output_file: str = 'positional_entropy_all_genes.pdf'):
    """
    Create a multi-page PDF with positional entropy plots for all genes.

    Args:
        output_file: Output PDF filename
    """
    # Find all .msa files
    msa_files = sorted(glob.glob("*.msa"))

    print("=" * 80)
    print("Positional Entropy Plotter")
    print("=" * 80)
    print(f"\nFound {len(msa_files)} alignment files")
    print(f"Generating plots for all genes...")
    print()

    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 100

    # Create PDF
    with PdfPages(output_file) as pdf:
        for idx, msa_file in enumerate(msa_files, 1):
            filename = Path(msa_file).name
            gene_id = filename.replace('.msa', '')

            print(f"[{idx}/{len(msa_files)}] Processing {gene_id}...")

            # Read and process alignment
            sequences = read_fasta(msa_file)
            num_sequences = len(sequences)
            original_length = len(list(sequences.values())[0]) if sequences else 0

            # Trim alignment
            trimmed_seqs, kept_indices = trim_alignment(sequences, gap_threshold=0.8)
            trimmed_length = len(trimmed_seqs[0]) if trimmed_seqs else 0

            # Calculate entropy
            entropy_array = calculate_shannon_entropy(trimmed_seqs, normalize=True)

            # Create figure
            fig, ax = plt.subplots(figsize=(14, 6))

            # Plot
            plot_gene_entropy(gene_id, entropy_array, num_sequences,
                            original_length, trimmed_length, ax)

            # Add page number
            fig.text(0.99, 0.01, f'Page {idx}/{len(msa_files)}',
                    ha='right', va='bottom', fontsize=8, style='italic', color='gray')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    print()
    print("=" * 80)
    print(f"PDF saved to: {output_file}")
    print(f"Total pages: {len(msa_files)}")
    print("=" * 80)


def main():
    """Main execution function."""
    create_entropy_pdf('positional_entropy_all_genes.pdf')


if __name__ == '__main__':
    main()
