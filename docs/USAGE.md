# Entropia-MSA Usage Guide

## Table of Contents
- [Basic Usage](#basic-usage)
- [Script Descriptions](#script-descriptions)
- [Input Requirements](#input-requirements)
- [Output Files](#output-files)
- [Advanced Options](#advanced-options)

## Basic Usage

### Step 1: Calculate Shannon Entropy

Navigate to your directory containing `.msa` alignment files:

```bash
cd /path/to/your/alignments/
python /path/to/entropia-msa/src/calculate_shannon_entropy.py
```

This will:
1. Find all `.msa` files in the current directory
2. Trim alignments (remove columns with >80% gaps)
3. Calculate normalized Shannon entropy for each position
4. Generate summary statistics

**Outputs**:
- `shannon_entropy_results.csv`
- `entropy_plot.png`

### Step 2: Generate Positional Entropy Profiles

```bash
python /path/to/entropia-msa/src/plot_positional_entropy.py
```

This creates a multi-page PDF with one plot per gene showing entropy across all alignment positions.

**Output**:
- `positional_entropy_all_genes.pdf` (182 pages in example dataset)

### Step 3: Integrated Heatmap (Optional)

For projects with gene copy number data:

```bash
cd /path/to/analysis/directory/
# Ensure you have:
# - shannon_entropy_results.csv
# - copy_number_data.tsv
# - phylogenetic_tree.txt

python /path/to/entropia-msa/src/plot_heatmap_with_entropy.py
```

**Outputs**:
- `kinetochore_heatmap_with_entropy.png`
- `kinetochore_heatmap_with_entropy.pdf`
- `kinetochore_heatmap_with_entropy.svg`

## Script Descriptions

### 1. calculate_shannon_entropy.py

**Purpose**: Calculate Shannon entropy for all alignment files

**Key Functions**:
- `read_fasta()`: Parse FASTA alignment files
- `trim_alignment()`: Remove poorly aligned regions
- `calculate_shannon_entropy()`: Compute normalized entropy
- `process_all_alignments()`: Batch process all `.msa` files

**Parameters**:
- `gap_threshold`: Fraction of gaps to trigger column removal (default: 0.8)
- `normalize`: Whether to normalize by max entropy (default: True)

**Formula**:
```
H_normalized = -Σ(p_i × log₂(p_i)) / log₂(20)
```

### 2. plot_positional_entropy.py

**Purpose**: Generate detailed positional entropy profiles

**Features**:
- Line plot with filled area showing entropy values
- Mean and median lines
- High variability regions highlighted (>75th percentile)
- Statistics box with mean, median, max, std
- Background zones (green=conserved, yellow=moderate, red=variable)

**Output Format**: Multi-page PDF with one gene per page

### 3. plot_heatmap_with_entropy.py

**Purpose**: Integrate entropy with gene copy number heatmap

**Components**:
- **Top panel**: Entropy barplot + Copy number variance (CV) line
- **Middle strip**: Protein complex color-coding
- **Bottom panel**: Gene copy number heatmap

**Features**:
- Dual y-axes (entropy + coefficient of variation)
- Phylogeny-ordered species
- Complex-colored protein labels
- Red stars for Arabidopsis orthogroups

## Input Requirements

### MSA Files (.msa)

- **Format**: FASTA format with aligned sequences
- **Naming**: `OGXXXXXXX_ProteinName.msa` (e.g., `OG0000103_H3.msa`)
- **Content**: Pre-aligned amino acid sequences with gaps

Example:
```
>Species1_gene
MARTKQTARKSTGGKAPR---KQLATKAAR
>Species2_gene
MARTKQTARKSTGGKAPR---KQLATKAAR
>Species3_gene
MARTK-TARKST--KAPR---KQLATKAAR
```

### Copy Number Data (for integrated heatmap)

- **Format**: Tab-separated file
- **Required columns**:
  - `label`: Protein name
  - `orthogroup`: OG identifier
  - `Complex`: Protein complex assignment
  - Species columns with copy numbers

### Phylogenetic Tree (for integrated heatmap)

- **Format**: Newick format
- **Example**: `(((Species1:0.1,Species2:0.2):0.3,Species3:0.4):0.5);`

## Output Files

### shannon_entropy_results.csv

Columns:
- `gene_id`: Gene/OG identifier
- `num_sequences`: Number of sequences
- `original_length`: Length before trimming
- `trimmed_length`: Length after gap removal
- `positions_removed`: Number of columns removed
- `percent_removed`: Percentage of positions trimmed
- `mean_normalized_entropy`: Average entropy (0-1 scale)
- `median_normalized_entropy`: Median entropy
- `std_normalized_entropy`: Standard deviation
- `max_normalized_entropy`: Maximum entropy value
- `min_normalized_entropy`: Minimum entropy value

### entropy_plot.png

- Bar plot showing top 50 genes by mean entropy
- Histogram showing distribution across all genes
- Mean and median lines

### positional_entropy_all_genes.pdf

Multi-page PDF:
- One plot per gene
- X-axis: Alignment position
- Y-axis: Normalized Shannon entropy (0-1)
- Red dots: High variability positions (>75th percentile)
- Statistics box in upper right

### Integrated Heatmap

Three formats (PNG, PDF, SVG):
- **Top**: Entropy bars (left y-axis) + CV line (right y-axis)
- **Middle**: Complex color strip
- **Bottom**: Copy number heatmap (grayscale + red for high copy)
- **X-axis**: Protein names (rotated 90°, color-coded by complex)
- **Y-axis**: Species/haplotypes (phylogeny order)

## Advanced Options

### Custom Gap Threshold

Edit Line ~40 in `calculate_shannon_entropy.py`:

```python
# More stringent (remove columns with >70% gaps)
trimmed_seqs = trim_alignment(sequences, gap_threshold=0.7)

# More lenient (remove columns with >90% gaps)
trimmed_seqs = trim_alignment(sequences, gap_threshold=0.9)
```

### Change File Pattern

Edit Line ~120 in `calculate_shannon_entropy.py`:

```python
# Look for .fasta files instead
msa_files = sorted(glob.glob(f"{directory}/*.fasta"))

# Look in subdirectories
msa_files = sorted(glob.glob(f"{directory}/**/*.msa", recursive=True))
```

### Adjust Entropy Scale

Edit Line ~224 in `plot_heatmap_with_entropy.py`:

```python
# Change y-axis limits
ax_entropy.set_ylim(0, 0.3)  # Focus on 0-0.3 range
ax_entropy.set_ylim(0, 1.0)  # Full 0-1 range
```

### Filter by Complex

Edit Line ~70 in `plot_heatmap_with_entropy.py`:

```python
# Keep only specific complex
df = df[df["Complex"] == "CCAN"]

# Remove multiple complexes
df = df[~df["Complex"].isin(["Kinase", "NDC80"])]
```

### Custom Figure Size

Edit Line ~192 in `plot_heatmap_with_entropy.py`:

```python
# Wider figure
figsize_x, figsize_y = 60, 38

# Taller figure
figsize_x, figsize_y = 50, 45
```

## Tips and Best Practices

### 1. Alignment Quality
- Ensure alignments are high quality (use MAFFT, MUSCLE, or Clustal Omega)
- Pre-filter highly divergent sequences
- Check alignment manually for obviously misaligned regions

### 2. Gap Threshold Selection
- **0.7-0.8**: Standard, removes poorly aligned regions
- **0.9**: More lenient, keeps more positions
- **0.6**: More stringent, focuses on well-aligned core

### 3. Interpretation
- **Entropy < 0.1**: Highly conserved (e.g., active sites, structural motifs)
- **Entropy 0.1-0.3**: Moderately conserved
- **Entropy > 0.3**: Variable (e.g., rapidly evolving, lineage-specific)

### 4. Performance
- Processing time: ~1-2 seconds per alignment
- Memory: ~100 MB for 100 alignments
- For >1000 alignments: Process in batches

### 5. Visualization
- Use PDF for publication-quality figures
- Use SVG for further editing in Illustrator/Inkscape
- Use PNG for presentations and quick sharing

## Troubleshooting

### No output files generated
- Check that `.msa` files exist in current directory
- Verify file permissions
- Check Python environment has required packages

### Empty entropy values
- Alignments may have too many gaps
- Try lowering gap threshold
- Check alignment quality

### PDF too large
- Reduce DPI in plot_positional_entropy.py (Line ~180)
- Process fewer genes
- Use PNG instead of PDF

### Memory issues
- Process alignments in smaller batches
- Reduce figure DPI
- Close figures after saving: `plt.close()`

## Examples

See `data/example_alignments/` for sample data:
- Highly conserved: H3 histone
- Moderately variable: Nuf2
- Highly variable: Skp1

Run on examples:
```bash
cd data/example_alignments/
python ../../src/calculate_shannon_entropy.py
python ../../src/plot_positional_entropy.py
```
