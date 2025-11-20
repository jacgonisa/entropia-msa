#!/usr/bin/env python3
"""
Integrated Heatmap with Shannon Entropy Barplot

This script creates a heatmap of kinetochore protein copy numbers across species
with an accompanying barplot showing mean normalized Shannon entropy for each protein.

Usage:
    python plot_heatmap_with_entropy.py

Output:
    - kinetochore_heatmap_with_entropy.png: Integrated visualization
    - kinetochore_heatmap_with_entropy.pdf: High-quality PDF version
    - kinetochore_heatmap_with_entropy.svg: Vector graphics version
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
print("Integrated Heatmap with Shannon Entropy")
print("=" * 80)
print()

# --- Load Entropy Data ---
print("Loading Shannon entropy data...")
entropy_df = pd.read_csv("shannon_entropy_results.csv")
# Extract orthogroup from gene_id (format: OG0000103_H3)
entropy_df['orthogroup'] = entropy_df['gene_id'].str.split('_').str[0]
entropy_df = entropy_df.set_index('orthogroup')
print(f"Loaded entropy data for {len(entropy_df)} genes")
print()

# --- Data Loading ---
print("Loading protein copy number data...")
df = pd.read_csv("../hits_arabidopsis.kinetochore_label_summary_with_complex.tsv", sep="\t")

# 1. Apply user-requested label replacements FIRST to ensure consistency
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
# Merge entropy data
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
    """Parse species names from Newick format tree."""
    with open(newick_file, 'r') as f:
        tree_str = f.read().strip()

    # Extract species names (format: Rspecies:branch_length)
    import re
    species = re.findall(r'([A-Z][a-z]+):', tree_str)
    # Remove duplicates while preserving order
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

# --- Build ordered lists with matching indices ---
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

# --- Order proteins by Complex and name ---
protein_order_df = pd.DataFrame({
    "Complex": protein_complex.loc[matrix_plot.columns].fillna("ZZZ"),
    "prot_name": matrix_plot.columns
}).sort_values(by=["Complex", "prot_name"])

protein_order = protein_order_df["prot_name"].tolist()
matrix_plot = matrix_plot[protein_order]

col_colors = protein_complex[matrix_plot.columns].map(complex_to_color)
col_colors.name = 'Complex'

# --- Create Figure with Custom Layout ---
print("Creating integrated visualization...")
print()

figsize_x, figsize_y = 50, 38

# Create figure with GridSpec for custom layout
fig = plt.figure(figsize=(figsize_x, figsize_y))
gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    height_ratios=[0.08, 0.012, 0.908],  # Entropy bar, col_colors, heatmap (more space for heatmap)
    width_ratios=[0.95, 0.05],
    hspace=0.02,
    wspace=0.02,
    bottom=0.25,  # Extra space at bottom for labels
    top=0.975,
    left=0.05,
    right=0.75
)

# --- Create Entropy Barplot ---
ax_entropy = fig.add_subplot(gs[0, 0])

# Get entropy values in the same order as proteins
entropy_values = protein_entropy[protein_order].fillna(0)
x_positions = np.arange(len(entropy_values))

# Calculate copy number variance (coefficient of variation) across species for each protein
copy_number_cv = []
for prot in protein_order:
    copy_counts = matrix_plot[prot].values
    if copy_counts.sum() > 0:
        # Coefficient of variation: std / mean (only for non-zero values to avoid division issues)
        non_zero = copy_counts[copy_counts > 0]
        if len(non_zero) > 0 and non_zero.mean() > 0:
            cv = non_zero.std() / non_zero.mean()
        else:
            cv = 0
    else:
        cv = 0
    copy_number_cv.append(cv)

copy_number_cv = np.array(copy_number_cv)

# Color bars by complex
bar_colors = [complex_to_color.get(protein_complex.get(prot), 'gray') for prot in protein_order]

# Plot entropy bars
bars = ax_entropy.bar(x_positions, entropy_values, color=bar_colors,
                      edgecolor='black', linewidth=0.3, alpha=0.8, label='Shannon Entropy')

# Add horizontal line for mean entropy
mean_entropy = entropy_values[entropy_values > 0].mean()
ax_entropy.axhline(mean_entropy, color='red', linestyle='--', linewidth=2,
                   label=f'Mean Entropy: {mean_entropy:.3f}', alpha=0.7)

ax_entropy.set_xlim(-0.5, len(entropy_values) - 0.5)
ax_entropy.set_ylim(0, 0.5)
ax_entropy.set_ylabel('Mean Normalized\nShannon Entropy', fontsize=14, fontweight='bold', color='black')
ax_entropy.set_xticks([])
ax_entropy.tick_params(axis='y', labelcolor='black')
ax_entropy.grid(axis='y', alpha=0.3, linestyle='--')

# Add subtle background zones
ax_entropy.axhspan(0, 0.2, alpha=0.08, color='green', zorder=0)
ax_entropy.axhspan(0.2, 0.5, alpha=0.08, color='yellow', zorder=0)

# --- Add secondary y-axis for copy number variance ---
ax_cv = ax_entropy.twinx()

# Plot copy number CV as line with markers
line_cv = ax_cv.plot(x_positions, copy_number_cv, color='darkviolet', marker='o',
                     markersize=3, linewidth=1.5, alpha=0.7, label='Copy # Variance (CV)', zorder=5)

# Add horizontal line for mean CV
mean_cv = copy_number_cv[copy_number_cv > 0].mean()
ax_cv.axhline(mean_cv, color='darkviolet', linestyle=':', linewidth=2,
              label=f'Mean CV: {mean_cv:.2f}', alpha=0.5)

ax_cv.set_ylim(0, max(copy_number_cv) * 1.1 if max(copy_number_cv) > 0 else 1)
ax_cv.set_ylabel('Copy Number Variance\n(Coefficient of Variation)', fontsize=14,
                 fontweight='bold', color='darkviolet')
ax_cv.tick_params(axis='y', labelcolor='darkviolet')

# Combined legend
lines1, labels1 = ax_entropy.get_legend_handles_labels()
lines2, labels2 = ax_cv.get_legend_handles_labels()
ax_entropy.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11, framealpha=0.9)

ax_entropy.set_title('Shannon Entropy and Copy Number Variance per Protein',
                     fontsize=14, fontweight='bold', pad=10)

# --- Create Heatmap in our custom layout ---
ax_heatmap = fig.add_subplot(gs[2, 0])

# Plot heatmap
im = ax_heatmap.imshow(matrix_plot.values, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

# Set ticks and labels
ax_heatmap.set_xticks(np.arange(len(matrix_plot.columns)))
ax_heatmap.set_yticks(np.arange(len(matrix_plot.index)))
ax_heatmap.set_xticklabels(matrix_plot.columns, rotation=90, fontsize=14, ha='center', va='top')
ax_heatmap.set_yticklabels(species_labels, fontsize=12)
ax_heatmap.xaxis.tick_bottom()
ax_heatmap.tick_params(axis='x', which='major', pad=12, length=10)
plt.setp(ax_heatmap.xaxis.get_majorticklabels(), visible=True)

# Grid lines
ax_heatmap.set_xticks([i + 0.5 for i in range(len(matrix_plot.columns))], minor=True)
ax_heatmap.set_yticks([i + 0.5 for i in range(len(matrix_plot.index))], minor=True)
ax_heatmap.grid(which='minor', color='lightgray', linestyle='-', linewidth=0.5)
ax_heatmap.tick_params(which='minor', length=0)

# Force x-tick labels to be visible
ax_heatmap.set_xlabel("Protein / Orthogroup", fontsize=15, fontweight='bold', labelpad=15)
ax_heatmap.set_ylabel("Species/Haplotype", fontsize=14, fontweight='bold')

# --- X-tick labels: complex-color + star ---
# Re-get the xticklabels after setting them to ensure they're rendered
xticklabels = ax_heatmap.get_xticklabels()
for i, label_obj in enumerate(xticklabels):
    prot_og = label_obj.get_text()
    # Complex color
    complex_name = protein_complex.get(prot_og)
    color = complex_to_color.get(complex_name, "black")
    label_obj.set_color(color)
    label_obj.set_fontweight('bold')
    label_obj.set_visible(True)
    # Arabidopsis OG: small red star above the label
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
ax_colors.set_xlim(ax_heatmap.get_xlim())  # Match heatmap x-limits without sharex
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
fig.suptitle("Kinetochore Protein Copy Number and Shannon Entropy\n★ = Arabidopsis OG",
             fontsize=18, fontweight='bold', y=0.985)

# --- Save outputs ---
print("Saving outputs...")
output_files = [
    "kinetochore_heatmap_with_entropy.png",
    "kinetochore_heatmap_with_entropy.pdf",
    "kinetochore_heatmap_with_entropy.svg"
]

for output_file in output_files:
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")

print()
print("=" * 80)
print("Analysis complete!")
print("=" * 80)
print(f"\nOutputs saved to: {output_files[0].rsplit('/', 1)[0] if '/' in output_files[0] else '.'}")

print(f"\nTop 5 proteins by Shannon entropy:")
top_entropy = protein_entropy.sort_values(ascending=False).head(5)
for prot, ent in top_entropy.items():
    print(f"  {prot:45s}  Entropy: {ent:.4f}")

print(f"\nTop 5 proteins by copy number variance (CV):")
# Create CV dataframe for sorting
cv_data = pd.Series(copy_number_cv, index=protein_order)
top_cv = cv_data.sort_values(ascending=False).head(5)
for prot, cv in top_cv.items():
    print(f"  {prot:45s}  CV: {cv:.4f}")

print(f"\nProteins with high entropy AND high copy number variance:")
# Find proteins in top quartile for both metrics
entropy_threshold = protein_entropy.quantile(0.75)
cv_threshold = cv_data.quantile(0.75)
high_both = []
for prot in protein_order:
    ent = protein_entropy.get(prot, 0)
    cv = cv_data.get(prot, 0)
    if ent > entropy_threshold and cv > cv_threshold:
        high_both.append((prot, ent, cv))

if high_both:
    high_both.sort(key=lambda x: x[1] + x[2], reverse=True)  # Sort by sum
    for prot, ent, cv in high_both[:10]:
        print(f"  {prot:45s}  Entropy: {ent:.4f}  CV: {cv:.4f}")

plt.show()
