#!/usr/bin/env python3
"""
Integrated Heatmap with Shannon Entropy Boxplot and Permutation Testing

This script creates a heatmap of kinetochore protein copy numbers across species
with boxplots showing the distribution of Shannon entropy across alignment positions.
Includes permutation-based null model for significance testing.

Usage:
    python plot_heatmap_with_entropy_boxplot.py

Output:
    - kinetochore_heatmap_with_entropy_boxplot.png
    - kinetochore_heatmap_with_entropy_boxplot.pdf
    - kinetochore_heatmap_with_entropy_boxplot.svg
    - entropy_permutation_results.csv (null model statistics)
"""

import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colorbar import ColorbarBase
import numpy as np
from io import StringIO
from collections import Counter

# --- Configuration ---
manual_mapping = {
    "Rridleyi": ["Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr",
                 "Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr"],
    "Rradicans": ["Rhync_radicans.asm.bp.p_ctg.FINAL.chr"],
    "Rpubera": ["rhyPub2m.chr"],
    "Rriparia": ["Rhync_riparia_5519C.hap1.chr",
                 "Rhync_riparia_5519C.hap2.chr"],
    "Rcorymbosa": ["Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr",
                   "Rhync_corymbosa_6179D.asm.hic.FINAL.hap2.chr"],
}

print("=" * 80)
print("Integrated Heatmap with Shannon Entropy Boxplot")
print("=" * 80)
print()


def calculate_permuted_entropy(sequences, n_permutations=1000):
    """
    Calculate null distribution of entropy using permutations.

    For each permutation, shuffle amino acids within each column independently,
    then calculate entropy. This preserves amino acid frequencies but destroys
    positional conservation patterns.

    Args:
        sequences: List of aligned sequences
        n_permutations: Number of permutations (default: 1000)

    Returns:
        Dictionary with null distribution statistics
    """
    if not sequences or len(sequences[0]) == 0:
        return None

    alignment_length = len(sequences[0])
    max_entropy = np.log2(20)

    null_entropies = []

    for perm in range(n_permutations):
        perm_entropy = []

        for col_idx in range(alignment_length):
            # Get column
            column = [seq[col_idx] for seq in sequences]
            # Shuffle column
            shuffled = np.random.permutation(column)

            # Calculate entropy
            aa_counts = Counter(aa for aa in shuffled if aa != '-')
            if not aa_counts:
                perm_entropy.append(0.0)
                continue

            total = sum(aa_counts.values())
            entropy = 0.0
            for count in aa_counts.values():
                if count > 0:
                    p_i = count / total
                    entropy -= p_i * np.log2(p_i)

            perm_entropy.append(entropy / max_entropy)

        # Store mean entropy for this permutation
        null_entropies.append(np.mean(perm_entropy))

    return {
        'null_mean': np.mean(null_entropies),
        'null_median': np.median(null_entropies),
        'null_std': np.std(null_entropies),
        'null_95th': np.percentile(null_entropies, 95),
        'null_99th': np.percentile(null_entropies, 99),
        'null_distribution': null_entropies
    }


# --- Load Entropy Data ---
print("Loading Shannon entropy data...")
entropy_df = pd.read_csv("shannon_entropy_results.csv")
entropy_df['orthogroup'] = entropy_df['gene_id'].str.split('_').str[0]
entropy_df = entropy_df.set_index('orthogroup')
print(f"Loaded entropy data for {len(entropy_df)} genes")
print()

# Load raw alignments to get positional entropy distributions
print("Loading alignment files for positional entropy distributions...")

def read_fasta(filepath):
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

def trim_alignment(sequences, gap_threshold=0.8):
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

def calculate_shannon_entropy(sequences, normalize=True):
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

# Store positional entropy for each protein
import glob
positional_entropy = {}
permutation_results = {}

msa_files = sorted(glob.glob("../mafft_msa_core_proteins/*.msa"))
if not msa_files:
    print("Warning: No MSA files found in ../mafft_msa_core_proteins/")
    print("Proceeding with summary statistics only...")
else:
    print(f"Found {len(msa_files)} MSA files. Calculating positional entropy...")
    for msa_file in msa_files[:20]:  # Limit to first 20 for speed in example
        from pathlib import Path
        gene_id = Path(msa_file).stem
        og = gene_id.split('_')[0]

        sequences = read_fasta(msa_file)
        trimmed_seqs = trim_alignment(sequences, gap_threshold=0.8)
        entropy_array = calculate_shannon_entropy(trimmed_seqs, normalize=True)

        if len(entropy_array) > 0:
            positional_entropy[gene_id] = entropy_array

            # Calculate permutation null model for selected genes
            if len(permutation_results) < 5:  # Only for first 5 genes (computational cost)
                print(f"  Calculating null model for {gene_id}...")
                null_stats = calculate_permuted_entropy(trimmed_seqs, n_permutations=100)
                if null_stats:
                    permutation_results[gene_id] = null_stats

print(f"Loaded positional entropy for {len(positional_entropy)} genes")
print()

# Save permutation results
if permutation_results:
    perm_df = pd.DataFrame(permutation_results).T
    perm_df.to_csv('entropy_permutation_results.csv')
    print(f"Saved permutation results to entropy_permutation_results.csv")
    print()

# --- Data Loading ---
print("Loading protein copy number data...")
df = pd.read_csv("../hits_arabidopsis.kinetochore_label_summary_with_complex.tsv", sep="\t")

df["label"] = df["label"].replace("CenpA", "CENH3")
df.loc[df["orthogroup"] == "OG0000103", "label"] = "H3"
df['protein_og'] = df['label'] + ' / ' + df['orthogroup']

metadata_cols = ['label', 'orthogroup', 'n_ogs', 'split', 'gene_count', 'Total',
                 'is_arabidopsis_og', 'is_arabidopsis_best_og', 'Complex', 'protein_og']

count_cols = [c for c in df.columns if re.search(r'_chr|\.chr|chr', c)]
species_cols = [c for c in count_cols if c not in metadata_cols]

# --- Filtering ---
df.loc[df["label"].str.contains("NASP|CENH3|H3", case=False, na=False), "Complex"] = "Histones and NASP"
df = df[df["Complex"] != "Kinase"]
df = df[~df["label"].str.contains("Dad4", case=False, na=False)]
og_to_remove = ["OG0013292", "OG0015310", "OG0015704"]
df = df[~df["orthogroup"].str.contains("|".join(og_to_remove), case=False, na=False)]

print(f"Filtered to {len(df)} proteins")
print()

# --- Merge with Entropy Data ---
print("Merging entropy data with protein data...")
df = df.set_index("protein_og")
df = df.join(entropy_df[['mean_normalized_entropy', 'median_normalized_entropy',
                         'trimmed_length']], on='orthogroup', how='left')

matrix = df[species_cols].T
protein_complex = df["Complex"]
arabidopsis_og_status = df["is_arabidopsis_og"]
protein_entropy = df["mean_normalized_entropy"]

print(f"Merged entropy data: {protein_entropy.notna().sum()}/{len(protein_entropy)} proteins have entropy values")
print()

# --- Tree ---
print("Loading phylogenetic tree...")

def parse_newick_species(newick_file):
    with open(newick_file, 'r') as f:
        tree_str = f.read().strip()
    species = re.findall(r'([A-Z][a-z]+):', tree_str)
    seen = set()
    species_ordered = []
    for sp in species:
        if sp not in seen:
            seen.add(sp)
            species_ordered.append(sp)
    return species_ordered

tree_species = parse_newick_species("../../../Rhynchospora_GENESPACE_SpeciesTree_rooted_node_labels.txt")

def normalize_tree_name(name):
    return name[1:].lower() if name.startswith("R") else name.lower()

def normalize_matrix_name(name):
    name = name.lower()
    for p in ["rhynchospora_", "rhync_", "r_"]:
        if name.startswith(p):
            name = name[len(p):]
            break
    if name.startswith("r") and len(name) > 1 and name[1].isalpha():
        name = name[1:]
    m = re.match(r"([a-z]+)", name)
    return m.group(1) if m else name

def get_display_name(tree_name, matrix_name):
    epithet = normalize_tree_name(tree_name)
    hap_match = re.search(r"(hap[12])", matrix_name.lower())
    haplotype = f" ({hap_match.group(1)})" if hap_match else ""
    return f"R. {epithet}{haplotype}"

# --- Build mapping ---
final_mapping = {}
for sp in tree_species:
    norm_sp = normalize_tree_name(sp)
    candidates = manual_mapping.get(sp) or [s for s in matrix.index if normalize_matrix_name(s) == norm_sp]
    if candidates:
        for matrix_name in candidates:
            final_mapping[get_display_name(sp, matrix_name)] = [matrix_name]
    else:
        final_mapping[f"R. {norm_sp}"] = []

ordered_display_keys = []
ordered_matrix_indices = []

for sp_tree in tree_species:
    epithet = normalize_tree_name(sp_tree)
    matching_keys = sorted([k for k in final_mapping.keys() if k.startswith(f"R. {epithet}")])
    for key in matching_keys:
        if final_mapping[key] and final_mapping[key][0] in matrix.index:
            ordered_display_keys.append(key)
            ordered_matrix_indices.append(final_mapping[key][0])

matrix_plot = matrix.loc[ordered_matrix_indices]
species_labels = ordered_display_keys

print(f"Matrix: {len(matrix_plot)} species/haplotypes × {len(matrix_plot.columns)} proteins")
print()

# --- Colormap ---
c_discrete = ["#f8f8f8", "#cccccc", "#888888", "#444444", "#000000", "#e41a1c"]
max_count = int(matrix_plot.values.max())
bounds = [0, 1, 2, 3, 4, 11, max(12, max_count + 1)]
norm = BoundaryNorm(boundaries=bounds, ncolors=len(c_discrete))
cmap = ListedColormap(c_discrete)

unique_complexes = protein_complex.dropna().unique()
palette_complex = sns.color_palette("tab20", len(unique_complexes))
complex_to_color = dict(zip(unique_complexes, palette_complex))

protein_order_df = pd.DataFrame({
    "Complex": protein_complex.loc[matrix_plot.columns].fillna("ZZZ"),
    "prot_name": matrix_plot.columns
}).sort_values(by=["Complex", "prot_name"])

protein_order = protein_order_df["prot_name"].tolist()
matrix_plot = matrix_plot[protein_order]

col_colors = protein_complex[matrix_plot.columns].map(complex_to_color)
col_colors.name = 'Complex'

# --- Create Figure with Custom Layout ---
print("Creating integrated visualization with boxplots...")
print()

figsize_x, figsize_y = 50, 38

fig = plt.figure(figsize=(figsize_x, figsize_y))
gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    height_ratios=[0.08, 0.012, 0.908],
    width_ratios=[0.95, 0.05],
    hspace=0.02,
    wspace=0.02,
    bottom=0.25,
    top=0.975,
    left=0.05,
    right=0.75
)

# --- Create Entropy Boxplot Panel ---
ax_entropy = fig.add_subplot(gs[0, 0])

# Prepare data for boxplot
boxplot_data = []
boxplot_positions = []
boxplot_colors = []

for i, prot in enumerate(protein_order):
    gene_id_match = [k for k in positional_entropy.keys() if prot.split(' / ')[1] in k]
    if gene_id_match:
        entropy_vals = positional_entropy[gene_id_match[0]]
        boxplot_data.append(entropy_vals)
        boxplot_positions.append(i)
        boxplot_colors.append(complex_to_color.get(protein_complex.get(prot), 'gray'))
    else:
        # Use mean value as placeholder if positional data not available
        mean_val = protein_entropy.get(prot, 0)
        boxplot_data.append([mean_val])
        boxplot_positions.append(i)
        boxplot_colors.append(complex_to_color.get(protein_complex.get(prot), 'gray'))

# Create boxplot
bp = ax_entropy.boxplot(boxplot_data, positions=boxplot_positions,
                        widths=0.6, patch_artist=True,
                        showfliers=False,  # Hide outliers for cleaner look
                        boxprops=dict(linewidth=0.5),
                        whiskerprops=dict(linewidth=0.5),
                        capprops=dict(linewidth=0.5),
                        medianprops=dict(color='red', linewidth=1.5))

# Color boxes by complex
for patch, color in zip(bp['boxes'], boxplot_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Add overall mean line
all_entropy_vals = [val for sublist in boxplot_data for val in sublist]
mean_entropy = np.mean(all_entropy_vals)
ax_entropy.axhline(mean_entropy, color='darkred', linestyle='--', linewidth=2,
                   label=f'Overall Mean: {mean_entropy:.3f}', alpha=0.7)

# Add permutation null line if available
if permutation_results:
    null_means = [v['null_mean'] for v in permutation_results.values()]
    avg_null = np.mean(null_means)
    ax_entropy.axhline(avg_null, color='purple', linestyle=':', linewidth=2,
                       label=f'Null Model: {avg_null:.3f}', alpha=0.7)

ax_entropy.set_xlim(-0.5, len(protein_order) - 0.5)
ax_entropy.set_ylim(0, 0.5)
ax_entropy.set_ylabel('Normalized Shannon Entropy\n(Positional Distribution)',
                      fontsize=14, fontweight='bold')
ax_entropy.set_xticks([])
ax_entropy.legend(loc='upper right', fontsize=11, framealpha=0.9)
ax_entropy.grid(axis='y', alpha=0.3, linestyle='--')
ax_entropy.set_title('Shannon Entropy Distribution per Protein (Boxplots)',
                     fontsize=14, fontweight='bold', pad=10)

# Add background zones
ax_entropy.axhspan(0, 0.2, alpha=0.08, color='green', zorder=0)
ax_entropy.axhspan(0.2, 0.5, alpha=0.08, color='yellow', zorder=0)

# --- Calculate copy number variance ---
copy_number_cv = []
for prot in protein_order:
    copy_counts = matrix_plot[prot].values
    if copy_counts.sum() > 0:
        non_zero = copy_counts[copy_counts > 0]
        if len(non_zero) > 0 and non_zero.mean() > 0:
            cv = non_zero.std() / non_zero.mean()
        else:
            cv = 0
    else:
        cv = 0
    copy_number_cv.append(cv)

copy_number_cv = np.array(copy_number_cv)

# --- Create Heatmap ---
ax_heatmap = fig.add_subplot(gs[2, 0])

im = ax_heatmap.imshow(matrix_plot.values, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

ax_heatmap.set_xticks(np.arange(len(matrix_plot.columns)))
ax_heatmap.set_yticks(np.arange(len(matrix_plot.index)))
ax_heatmap.set_xticklabels(matrix_plot.columns, rotation=90, fontsize=14, ha='center', va='top')
ax_heatmap.set_yticklabels(species_labels, fontsize=12)
ax_heatmap.xaxis.tick_bottom()
ax_heatmap.tick_params(axis='x', which='major', pad=12, length=10)
plt.setp(ax_heatmap.xaxis.get_majorticklabels(), visible=True)

ax_heatmap.set_xticks([i + 0.5 for i in range(len(matrix_plot.columns))], minor=True)
ax_heatmap.set_yticks([i + 0.5 for i in range(len(matrix_plot.index))], minor=True)
ax_heatmap.grid(which='minor', color='lightgray', linestyle='-', linewidth=0.5)
ax_heatmap.tick_params(which='minor', length=0)

ax_heatmap.set_xlabel("Protein / Orthogroup", fontsize=15, fontweight='bold', labelpad=15)
ax_heatmap.set_ylabel("Species/Haplotype", fontsize=14, fontweight='bold')

xticklabels = ax_heatmap.get_xticklabels()
for i, label_obj in enumerate(xticklabels):
    prot_og = label_obj.get_text()
    complex_name = protein_complex.get(prot_og)
    color = complex_to_color.get(complex_name, "black")
    label_obj.set_color(color)
    label_obj.set_fontweight('bold')
    label_obj.set_visible(True)
    val = str(arabidopsis_og_status.get(prot_og)).strip().lower()
    if val in {"true", "1", "yes"}:
        ax_heatmap.plot(i, -0.5, marker="*", color="red", markersize=12,
                       markeredgecolor="red", clip_on=False, zorder=10)

# --- Add col_colors strip ---
ax_colors = fig.add_subplot(gs[1, 0])
col_colors_array = np.array([[complex_to_color.get(protein_complex.get(p), 'gray')
                              for p in protein_order]])
ax_colors.imshow(col_colors_array, aspect='auto', interpolation='nearest')
ax_colors.set_xticks([])
ax_colors.set_yticks([])
ax_colors.set_xlim(ax_heatmap.get_xlim())
ax_colors.set_ylabel('Complex', fontsize=10, fontweight='bold', rotation=0, ha='right', va='center')

# --- Custom colorbar ---
cbar_ax = fig.add_axes([0.02, 0.7, 0.012, 0.15])
cbar_labels = ['0', '1', '2', '3', '4-10', '>10']
cb = ColorbarBase(cbar_ax, cmap=cmap, norm=norm,
                  ticks=[0.5, 1.5, 2.5, 3.5, 7.5, 11],
                  orientation='vertical')
cb.set_ticklabels(cbar_labels)
cb.ax.tick_params(labelsize=11)
cb.set_label('Gene Copy Number', fontsize=13, fontweight='bold')

# --- Legend ---
sorted_complexes = sorted(complex_to_color.keys())
handles = [mpatches.Patch(color=complex_to_color[comp], label=comp) for comp in sorted_complexes]
labels = sorted_complexes

fig.legend(
    handles=handles,
    labels=labels,
    title="Protein Complex",
    title_fontsize=16,
    fontsize=13,
    loc="upper left",
    bbox_to_anchor=(0.77, 0.92),
    frameon=True,
    fancybox=True,
    shadow=True
)

# --- Main title ---
fig.suptitle("Kinetochore Protein Copy Number and Shannon Entropy Distribution\n★ = Arabidopsis OG",
             fontsize=18, fontweight='bold', y=0.985)

# --- Save outputs ---
print("Saving outputs...")
output_files = [
    "kinetochore_heatmap_with_entropy_boxplot.png",
    "kinetochore_heatmap_with_entropy_boxplot.pdf",
    "kinetochore_heatmap_with_entropy_boxplot.svg"
]

for output_file in output_files:
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")

print()
print("=" * 80)
print("Analysis complete!")
print("=" * 80)
print(f"\nOutputs saved to current directory")

if permutation_results:
    print(f"\nPermutation Test Results (null model):")
    for gene_id, stats in list(permutation_results.items())[:3]:
        observed = np.mean(positional_entropy[gene_id])
        print(f"\n  {gene_id}:")
        print(f"    Observed mean entropy: {observed:.4f}")
        print(f"    Null mean entropy:     {stats['null_mean']:.4f}")
        print(f"    Null 95th percentile:  {stats['null_95th']:.4f}")
        print(f"    Significant: {'Yes' if observed > stats['null_95th'] else 'No'} (p < 0.05)")

plt.show()
