#!/usr/bin/env python3
"""
Permutation Test for Shannon Entropy Significance

This script performs permutation tests to establish null distributions for Shannon entropy.
Useful for determining if observed entropy values are statistically significant.

Usage:
    python permutation_test_entropy.py

Output:
    - permutation_test_results.csv: Observed vs null entropy with p-values
    - permutation_distributions.png: Visualization of null distributions
"""

import numpy as np
from collections import Counter
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict
import glob


def read_fasta(filepath: str) -> Dict[str, str]:
    """Read a FASTA format alignment file."""
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


def trim_alignment(sequences: Dict[str, str], gap_threshold: float = 0.8) -> List[str]:
    """Trim alignment columns with gaps in more than gap_threshold of sequences."""
    if not sequences:
        return []

    seq_list = list(sequences.values())
    num_seqs = len(seq_list)
    alignment_length = len(seq_list[0])

    columns_to_keep = []
    for col_idx in range(alignment_length):
        gap_count = sum(1 for seq in seq_list if seq[col_idx] == '-')
        gap_fraction = gap_count / num_seqs

        if gap_fraction <= gap_threshold:
            columns_to_keep.append(col_idx)

    trimmed_seqs = []
    for seq in seq_list:
        trimmed_seq = ''.join(seq[i] for i in columns_to_keep)
        trimmed_seqs.append(trimmed_seq)

    return trimmed_seqs


def calculate_shannon_entropy(sequences: List[str], normalize: bool = True) -> np.ndarray:
    """Calculate Shannon entropy for each position in the alignment."""
    if not sequences or len(sequences[0]) == 0:
        return np.array([])

    alignment_length = len(sequences[0])
    max_entropy = np.log2(20)

    entropies = []

    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in sequences]
        aa_counts = Counter(aa for aa in column if aa != '-')

        if not aa_counts:
            entropies.append(0.0)
            continue

        total = sum(aa_counts.values())
        entropy = 0.0

        for count in aa_counts.values():
            if count > 0:
                p_i = count / total
                entropy -= p_i * np.log2(p_i)

        if normalize:
            entropy = entropy / max_entropy

        entropies.append(entropy)

    return np.array(entropies)


def permutation_test(sequences: List[str], n_permutations: int = 1000,
                     normalize: bool = True) -> Dict:
    """
    Perform permutation test for Shannon entropy.

    Strategy: For each permutation, shuffle amino acids within each column
    independently. This preserves the amino acid composition at each position
    but destroys any positional conservation patterns.

    Args:
        sequences: List of aligned sequences
        n_permutations: Number of permutations to perform
        normalize: Whether to normalize entropy values

    Returns:
        Dictionary containing observed and null statistics
    """
    if not sequences or len(sequences[0]) == 0:
        return None

    # Calculate observed entropy
    observed_entropy = calculate_shannon_entropy(sequences, normalize=normalize)
    observed_mean = np.mean(observed_entropy)
    observed_median = np.median(observed_entropy)
    observed_std = np.std(observed_entropy)

    alignment_length = len(sequences[0])
    max_entropy = np.log2(20) if normalize else 1.0

    # Perform permutations
    null_means = []
    null_medians = []
    null_position_entropies = []

    print(f"  Running {n_permutations} permutations...")

    for perm_idx in range(n_permutations):
        if (perm_idx + 1) % 100 == 0:
            print(f"    Permutation {perm_idx + 1}/{n_permutations}")

        perm_entropies = []

        for col_idx in range(alignment_length):
            # Get column
            column = [seq[col_idx] for seq in sequences]

            # Shuffle amino acids within column
            shuffled = list(np.random.permutation(column))

            # Calculate entropy
            aa_counts = Counter(aa for aa in shuffled if aa != '-')

            if not aa_counts:
                perm_entropies.append(0.0)
                continue

            total = sum(aa_counts.values())
            entropy = 0.0

            for count in aa_counts.values():
                if count > 0:
                    p_i = count / total
                    entropy -= p_i * np.log2(p_i)

            if normalize:
                entropy = entropy / max_entropy

            perm_entropies.append(entropy)

        # Store statistics for this permutation
        null_means.append(np.mean(perm_entropies))
        null_medians.append(np.median(perm_entropies))

        # Store first few permutations for visualization
        if len(null_position_entropies) < 10:
            null_position_entropies.append(perm_entropies)

    # Calculate p-value (one-tailed: observed > null)
    p_value = np.sum(np.array(null_means) >= observed_mean) / n_permutations

    # Empirical p-value adjustment (add 1 to numerator and denominator)
    p_value_adj = (np.sum(np.array(null_means) >= observed_mean) + 1) / (n_permutations + 1)

    return {
        'observed_mean': observed_mean,
        'observed_median': observed_median,
        'observed_std': observed_std,
        'observed_entropy': observed_entropy,
        'null_mean': np.mean(null_means),
        'null_median': np.median(null_medians),
        'null_std': np.std(null_means),
        'null_5th': np.percentile(null_means, 5),
        'null_95th': np.percentile(null_means, 95),
        'null_99th': np.percentile(null_means, 99),
        'null_distribution': null_means,
        'null_position_examples': null_position_entropies,
        'p_value': p_value,
        'p_value_adjusted': p_value_adj,
        'n_permutations': n_permutations
    }


def plot_permutation_results(results_dict: Dict[str, Dict], output_file: str):
    """
    Create visualization of permutation test results.

    Args:
        results_dict: Dictionary mapping gene IDs to permutation results
        output_file: Output filename for plot
    """
    n_genes = len(results_dict)

    # Create figure
    fig, axes = plt.subplots(n_genes, 2, figsize=(16, 5 * n_genes))

    if n_genes == 1:
        axes = axes.reshape(1, -1)

    for idx, (gene_id, results) in enumerate(results_dict.items()):
        ax1, ax2 = axes[idx]

        # Plot 1: Null distribution histogram
        ax1.hist(results['null_distribution'], bins=30, color='lightblue',
                edgecolor='black', alpha=0.7, label='Null distribution')

        # Add observed value
        ax1.axvline(results['observed_mean'], color='red', linewidth=3,
                   label=f"Observed: {results['observed_mean']:.4f}")

        # Add significance threshold
        ax1.axvline(results['null_95th'], color='orange', linewidth=2,
                   linestyle='--', label=f"95th percentile: {results['null_95th']:.4f}")

        ax1.set_xlabel('Mean Normalized Shannon Entropy', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax1.set_title(f'{gene_id}\nPermutation Test (n={results["n_permutations"]})',
                     fontsize=13, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(alpha=0.3)

        # Add p-value text
        p_val = results['p_value_adjusted']
        significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        ax1.text(0.98, 0.95, f'p = {p_val:.4f} {significance}',
                transform=ax1.transAxes, fontsize=12, fontweight='bold',
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        # Plot 2: Positional entropy comparison
        positions = np.arange(len(results['observed_entropy']))

        # Plot observed
        ax2.plot(positions, results['observed_entropy'], color='red',
                linewidth=2, alpha=0.8, label='Observed', zorder=5)

        # Plot a few null examples
        for i, null_entropy in enumerate(results['null_position_examples'][:5]):
            ax2.plot(positions, null_entropy, color='lightgray',
                    linewidth=0.5, alpha=0.5, zorder=1)

        # Plot mean of nulls
        if results['null_position_examples']:
            null_mean_pos = np.mean(results['null_position_examples'], axis=0)
            ax2.plot(positions, null_mean_pos, color='blue', linewidth=2,
                    linestyle='--', alpha=0.7, label='Null mean', zorder=3)

        ax2.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Normalized Shannon Entropy', fontsize=12, fontweight='bold')
        ax2.set_title(f'Positional Entropy: Observed vs Null',
                     fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(alpha=0.3)
        ax2.set_ylim(0, 1.0)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")


def main():
    """Main execution function."""
    print("=" * 80)
    print("Permutation Test for Shannon Entropy")
    print("=" * 80)
    print()

    # Find MSA files
    msa_files = sorted(glob.glob("*.msa"))

    if not msa_files:
        print("No .msa files found in current directory!")
        return

    print(f"Found {len(msa_files)} alignment files")
    print("Selecting first 5 for permutation testing (computational cost)...")
    print()

    results_dict = {}
    summary_data = []

    for msa_file in msa_files[:5]:
        gene_id = Path(msa_file).stem

        print(f"Processing: {gene_id}")

        # Read and trim alignment
        sequences = read_fasta(msa_file)
        trimmed_seqs = trim_alignment(sequences, gap_threshold=0.8)

        if not trimmed_seqs:
            print(f"  Skipping (no sequences after trimming)")
            continue

        print(f"  {len(sequences)} sequences, {len(trimmed_seqs[0])} positions after trimming")

        # Perform permutation test
        results = permutation_test(trimmed_seqs, n_permutations=1000, normalize=True)

        if results:
            results_dict[gene_id] = results

            # Store summary
            summary_data.append({
                'gene_id': gene_id,
                'num_sequences': len(sequences),
                'alignment_length': len(trimmed_seqs[0]),
                'observed_mean': results['observed_mean'],
                'observed_std': results['observed_std'],
                'null_mean': results['null_mean'],
                'null_std': results['null_std'],
                'null_95th': results['null_95th'],
                'null_99th': results['null_99th'],
                'p_value': results['p_value'],
                'p_value_adjusted': results['p_value_adjusted'],
                'significant_0.05': results['p_value_adjusted'] < 0.05,
                'significant_0.01': results['p_value_adjusted'] < 0.01,
            })

            print(f"  Observed mean: {results['observed_mean']:.4f}")
            print(f"  Null mean:     {results['null_mean']:.4f}")
            print(f"  p-value:       {results['p_value_adjusted']:.4f}")
            print()

    # Save summary results
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv('permutation_test_results.csv', index=False)
        print("=" * 80)
        print("Results saved to: permutation_test_results.csv")
        print("=" * 80)
        print()

        # Print summary
        print("Summary:")
        print(f"  Genes tested: {len(summary_df)}")
        print(f"  Significant at p < 0.05: {summary_df['significant_0.05'].sum()}")
        print(f"  Significant at p < 0.01: {summary_df['significant_0.01'].sum()}")
        print()

        # Create visualization
        print("Creating visualization...")
        plot_permutation_results(results_dict, 'permutation_distributions.png')

        print()
        print("=" * 80)
        print("Analysis complete!")
        print("=" * 80)


if __name__ == '__main__':
    main()
