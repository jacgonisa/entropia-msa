#!/usr/bin/env python3
"""
Shannon Entropy Violin Plots for Rhynchospora Kinetochore Proteins

Creates violin plots showing the distribution of entropy across alignment positions
for each protein, with jittered points and permutation-based null model.

Usage:
    python plot_entropy_violins.py

Output:
    - entropy_violins_all_proteins.pdf
    - entropy_violins_all_proteins.png
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


def calculate_null_entropy(sequences: List[str], n_permutations: int = 100) -> Dict:
    """
    Calculate null distribution by permuting amino acids within each column.

    This destroys positional conservation while preserving AA composition.
    """
    if not sequences or len(sequences[0]) == 0:
        return None

    alignment_length = len(sequences[0])
    max_entropy = np.log2(20)

    null_entropies_all = []

    for perm in range(n_permutations):
        perm_entropies = []

        for col_idx in range(alignment_length):
            # Get column and shuffle it
            column = [seq[col_idx] for seq in sequences]
            shuffled = np.random.permutation(column).tolist()

            # Calculate entropy of shuffled column
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

            perm_entropies.append(entropy / max_entropy)

        null_entropies_all.extend(perm_entropies)

    return {
        'null_values': null_entropies_all,
        'null_mean': np.mean(null_entropies_all),
        'null_median': np.median(null_entropies_all),
        'null_std': np.std(null_entropies_all)
    }


print("=" * 80)
print("Shannon Entropy Violin Plots for Rhynchospora Kinetochore Proteins")
print("=" * 80)
print()

# Load entropy data
print("Loading Shannon entropy results...")
entropy_df = pd.read_csv("shannon_entropy_results.csv")
entropy_df['orthogroup'] = entropy_df['gene_id'].str.split('_').str[0]
print(f"Loaded summary statistics for {len(entropy_df)} genes")
print()

# Load alignments and calculate positional entropy
print("Loading MSA files and calculating positional entropy...")
msa_files = sorted(glob.glob("../mafft_msa_core_proteins/*.msa"))

if not msa_files:
    print("ERROR: No MSA files found!")
    exit(1)

print(f"Found {len(msa_files)} alignment files")
print()

# Process all alignments
all_data = []
gene_order = []

for idx, msa_file in enumerate(msa_files, 1):
    gene_id = Path(msa_file).stem

    if idx % 20 == 0:
        print(f"  Processed {idx}/{len(msa_files)} alignments...")

    # Read and process alignment
    sequences = read_fasta(msa_file)
    trimmed_seqs = trim_alignment(sequences, gap_threshold=0.8)

    if not trimmed_seqs or len(trimmed_seqs[0]) == 0:
        continue

    # Calculate observed entropy
    entropy_array = calculate_shannon_entropy(trimmed_seqs, normalize=True)

    if len(entropy_array) == 0:
        continue

    # Store gene order
    gene_order.append(gene_id)

    # Add observed values to data
    for val in entropy_array:
        all_data.append({
            'gene_id': gene_id,
            'entropy': val,
            'type': 'observed'
        })

    # Calculate null distribution (only for every 10th gene to save time)
    if idx % 10 == 0:
        null_results = calculate_null_entropy(trimmed_seqs, n_permutations=50)
        if null_results:
            for val in null_results['null_values'][:len(entropy_array)]:  # Match sample size
                all_data.append({
                    'gene_id': gene_id,
                    'entropy': val,
                    'type': 'null'
                })

print(f"  Processed {len(msa_files)} alignments")
print(f"  Collected data for {len(gene_order)} genes")
print()

# Create DataFrame
df_plot = pd.DataFrame(all_data)

# Sort genes by mean observed entropy
gene_means = df_plot[df_plot['type'] == 'observed'].groupby('gene_id')['entropy'].mean().sort_values(ascending=False)
top_genes = gene_means.head(50).index.tolist()

# Filter to top 50 genes
df_plot = df_plot[df_plot['gene_id'].isin(top_genes)]

# Set gene order
df_plot['gene_id'] = pd.Categorical(df_plot['gene_id'], categories=top_genes, ordered=True)

print(f"Creating violin plot for top 50 genes by mean entropy...")
print()

# Create figure
fig, ax = plt.subplots(figsize=(20, 30))

# Plot violins
sns.violinplot(data=df_plot[df_plot['type'] == 'observed'],
               y='gene_id', x='entropy',
               color='steelblue', alpha=0.6, ax=ax, orient='h',
               inner=None, cut=0)

# Add jittered points for observed values
observed_data = df_plot[df_plot['type'] == 'observed']
for i, gene in enumerate(top_genes):
    gene_data = observed_data[observed_data['gene_id'] == gene]['entropy'].values

    # Add jitter to y-position
    y_positions = np.random.normal(i, 0.1, size=len(gene_data))

    ax.scatter(gene_data, y_positions, alpha=0.3, s=3, color='darkblue', zorder=3)

# Plot null distribution (if available)
null_data = df_plot[df_plot['type'] == 'null']
if len(null_data) > 0:
    # Get genes with null data
    genes_with_null = null_data['gene_id'].unique()

    for gene in genes_with_null:
        gene_null = null_data[null_data['gene_id'] == gene]['entropy'].values
        gene_idx = top_genes.index(gene)

        # Plot null distribution as light violin on the right side
        parts = ax.violinplot([gene_null], positions=[gene_idx],
                              vert=False, widths=0.7,
                              showmeans=False, showmedians=False, showextrema=False)

        for pc in parts['bodies']:
            pc.set_facecolor('orange')
            pc.set_alpha(0.3)
            pc.set_edgecolor('orange')
            pc.set_linewidth(1)

# Add mean line
mean_entropy = df_plot[df_plot['type'] == 'observed']['entropy'].mean()
ax.axvline(mean_entropy, color='red', linestyle='--', linewidth=2,
           label=f'Mean: {mean_entropy:.3f}', alpha=0.7, zorder=1)

# Formatting
ax.set_xlabel('Normalized Shannon Entropy', fontsize=16, fontweight='bold')
ax.set_ylabel('Gene ID', fontsize=16, fontweight='bold')
ax.set_title('Shannon Entropy Distribution Across Alignment Positions\nTop 50 Genes by Mean Entropy',
             fontsize=18, fontweight='bold', pad=20)

# Add background zones
ax.axvspan(0, 0.2, alpha=0.05, color='green', zorder=0)
ax.axvspan(0.2, 0.5, alpha=0.05, color='yellow', zorder=0)
ax.axvspan(0.5, 1.0, alpha=0.05, color='red', zorder=0)

ax.set_xlim(-0.05, 1.0)
ax.grid(axis='x', alpha=0.3, linestyle='--')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='steelblue', alpha=0.6, label='Observed distribution'),
    Patch(facecolor='orange', alpha=0.3, label='Null model (permuted)'),
    plt.Line2D([0], [0], color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_entropy:.3f}')
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=14, frameon=True)

plt.tight_layout()

# Save
print("Saving outputs...")
output_files = [
    'entropy_violins_all_proteins.png',
    'entropy_violins_all_proteins.pdf'
]

for output_file in output_files:
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")

print()
print("=" * 80)
print("Analysis complete!")
print("=" * 80)
print()
print("Summary:")
print(f"  Total genes analyzed: {len(gene_order)}")
print(f"  Genes plotted: {len(top_genes)}")
print(f"  Overall mean entropy: {mean_entropy:.4f}")
print(f"  Overall median entropy: {df_plot[df_plot['type'] == 'observed']['entropy'].median():.4f}")
print()
print("Top 5 genes by mean entropy:")
for i, gene in enumerate(gene_means.head(5).items(), 1):
    print(f"  {i}. {gene[0]:45s}  Mean: {gene[1]:.4f}")

plt.show()
