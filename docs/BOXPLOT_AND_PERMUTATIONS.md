# Boxplot Visualization and Permutation Testing

## Overview

This document describes the advanced features for entropy analysis:
1. **Boxplot visualization** - Shows distribution of entropy across alignment positions
2. **Permutation testing** - Statistical significance testing with null models

## Why Boxplots Instead of Barplots?

### Barplot Limitations
- Shows only mean entropy (single summary statistic)
- Hides variability and distribution shape
- Cannot detect bimodal distributions

### Boxplot Advantages
- **Distribution shape**: See quartiles, spread, and skewness
- **Bimodal detection**: Identify conserved domains vs. disordered regions
- **Outlier identification**: Detect hypervariable positions
- **Inter-quartile range**: Understand consistency of conservation

### Biological Interpretation

#### Narrow Box (Low IQR)
- Uniformly conserved or variable across alignment
- Example: **Histone H3** - consistently low entropy (highly conserved)

#### Wide Box (High IQR)
- Mixed conservation patterns
- Example: **Multi-domain proteins** with conserved catalytic sites and variable linkers

#### Bimodal Distribution
- Two distinct regions in the protein
- **Low entropy mode**: Conserved domain (e.g., kinase domain)
- **High entropy mode**: Disordered region (e.g., regulatory tail)

## Permutation Testing

### Why Permutations?

**Question**: Is the observed entropy value significant, or just random noise?

**Problem**: For real proteins, entropy rarely approaches 1.0 because:
- Not all 20 amino acids occur at each position
- Biochemical constraints limit diversity
- Structural requirements favor certain residues

**Solution**: Permutation test establishes an empirical null distribution

### How It Works

```
For each permutation:
  1. For each alignment column:
     - Shuffle amino acids within the column
     - This preserves AA composition but destroys conservation patterns
  2. Calculate entropy for shuffled alignment
  3. Record mean entropy

After N permutations:
  - Compare observed entropy to null distribution
  - Calculate p-value: P(null ≥ observed)
```

### Interpretation

| p-value | Interpretation |
|---------|---------------|
| < 0.001 | Highly significant (***) - Strong conservation signal |
| < 0.01  | Very significant (**) - Clear conservation pattern |
| < 0.05  | Significant (*) - Moderate conservation |
| ≥ 0.05  | Not significant (ns) - Entropy consistent with random |

### Example Results

```
Gene: OG0000103_H3 (Histone H3)
  Observed mean entropy: 0.0095
  Null mean entropy:     0.4523
  95th percentile:       0.4689
  p-value:              < 0.001 (***)

Interpretation: H3 is HIGHLY conserved (far below null expectation)
```

```
Gene: OG0000122_Skp1 (Variable Skp1)
  Observed mean entropy: 0.3932
  Null mean entropy:     0.4201
  95th percentile:       0.4356
  p-value:               0.045 (*)

Interpretation: Skp1 shows moderate conservation (slightly below null)
```

## Usage

### 1. Boxplot Visualization (Integrated Heatmap)

```bash
cd /path/to/analysis/directory/
python /path/to/entropia-msa/src/plot_heatmap_with_entropy_boxplot.py
```

**Requirements**:
- `shannon_entropy_results.csv`
- MSA files in `../mafft_msa_core_proteins/`
- Copy number data (`.tsv`)
- Phylogenetic tree

**Outputs**:
- `kinetochore_heatmap_with_entropy_boxplot.{png,pdf,svg}`
- `entropy_permutation_results.csv` (null model stats for first 5 genes)

**Features**:
- Boxplots showing entropy distribution per protein
- Null model line (purple dashed)
- Overall mean line (red dashed)
- Color-coded by protein complex

### 2. Standalone Permutation Testing

```bash
cd /path/to/your/alignments/
python /path/to/entropia-msa/src/permutation_test_entropy.py
```

**Outputs**:
- `permutation_test_results.csv`: Summary statistics and p-values
- `permutation_distributions.png`: Visualization of null distributions

**Computation Time**:
- ~1-2 seconds per gene per 100 permutations
- Default: 1000 permutations
- For 5 genes: ~1-2 minutes

## Customization

### Adjust Number of Permutations

Edit `permutation_test_entropy.py` (Line ~250):

```python
# More permutations (more accurate p-values)
results = permutation_test(trimmed_seqs, n_permutations=10000)

# Fewer permutations (faster, less accurate)
results = permutation_test(trimmed_seqs, n_permutations=100)
```

### Process More Genes

Edit `plot_heatmap_with_entropy_boxplot.py` (Line ~144):

```python
# Process first 50 genes instead of 20
for msa_file in msa_files[:50]:
```

Edit `permutation_test_entropy.py` (Line ~260):

```python
# Test first 10 genes
for msa_file in msa_files[:10]:
```

### Change Significance Threshold

Default is one-tailed test (observed > null).

For two-tailed test, edit `permutation_test_entropy.py` (Line ~155):

```python
# Two-tailed p-value
p_value = np.sum(np.abs(null_means - null_mean) >= np.abs(observed_mean - null_mean)) / n_permutations
```

## Interpreting Boxplots

### Components

```
     ┌─────┐
  ───┤     │  ← 75th percentile (Q3)
     │  ━  │  ← Median (red line)
  ───┤     │  ← 25th percentile (Q1)
     └─────┘
      ╰─╯
   Whiskers (1.5 × IQR)
```

- **Box**: Inter-quartile range (IQR) - middle 50% of data
- **Median line**: Red line inside box
- **Whiskers**: Extend to 1.5 × IQR
- **Outliers**: Points beyond whiskers (hidden in default plot)

### Example Patterns

#### Pattern 1: Highly Conserved (e.g., H3)
```
Entropy
1.0 ┤
0.5 ┤
0.2 ┤
0.0 ┼ ▂  ← Tiny box near zero
```
**Interpretation**: Uniformly conserved across all positions

#### Pattern 2: Moderately Variable (e.g., Nuf2)
```
Entropy
1.0 ┤
0.5 ┤
0.3 ┤  ┌───┐
0.1 ┤──┤   │
0.0 ┼  └───┘
```
**Interpretation**: Most positions conserved, some variable

#### Pattern 3: Bimodal (e.g., Kinase with disordered tail)
```
Entropy
1.0 ┤
0.6 ┤    ╭───╮  ← High entropy mode
0.3 ┼────┴───┴
0.1 ┤  ╭─╮      ← Low entropy mode
0.0 ┼──┴─┴
```
**Interpretation**: Two distinct regions - conserved domain + variable region

## Advanced Analysis

### Identifying Domain Boundaries

Use positional entropy profiles (`plot_positional_entropy.py`) combined with boxplot statistics:

1. **High IQR proteins**: Likely have multiple domains
2. **Check positional plot**: Look for transitions from low→high entropy
3. **These transitions often mark domain boundaries**

### Detecting Intrinsically Disordered Regions (IDRs)

**IDRs show**:
- High entropy (>0.4)
- Low complexity (few amino acid types)
- Long stretches of variable positions

**Strategy**:
1. Identify proteins with high median entropy
2. Check positional plot for sustained high-entropy regions
3. Cross-reference with disorder predictors (IUPred, PONDR)

### Comparing Orthogroups

Use permutation results to ask:
- **Which orthogroups are under stronger purifying selection?**
  → Lower observed entropy, larger difference from null

- **Which show lineage-specific variation?**
  → High observed entropy, not significantly different from null

## Performance Optimization

### For Large Datasets (>100 genes)

1. **Batch processing**:
```python
# Process in chunks of 20
for chunk in range(0, len(msa_files), 20):
    batch = msa_files[chunk:chunk+20]
    # Process batch...
```

2. **Reduce permutations**:
```python
# Use 100 permutations for preliminary analysis
n_permutations = 100  # vs. default 1000
```

3. **Parallelize** (advanced):
```python
from multiprocessing import Pool

def process_gene(msa_file):
    # ... permutation test code ...
    return results

with Pool(4) as p:  # Use 4 cores
    results = p.map(process_gene, msa_files)
```

## Statistical Notes

### Multiple Testing Correction

When testing many genes, apply Bonferroni or FDR correction:

```python
from scipy.stats import false_discovery_control

# After running permutation tests
p_values = [result['p_value_adjusted'] for result in results_dict.values()]

# FDR correction (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
reject, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# Bonferroni correction (more conservative)
reject, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
```

### Power Analysis

Number of permutations needed for desired precision:

| Desired p-value precision | Min permutations |
|--------------------------|------------------|
| p = 0.05 ± 0.01 | 1,000 |
| p = 0.01 ± 0.001 | 10,000 |
| p = 0.001 ± 0.0001 | 100,000 |

## References

- **Shannon Entropy**: Shannon, C.E. (1948). Bell System Technical Journal
- **Permutation Tests**: Good, P. (2000). Permutation Tests: A Practical Guide
- **Multiple Testing**: Benjamini, Y. & Hochberg, Y. (1995). JRSS-B

## Troubleshooting

### Issue: All p-values = 1.0

**Cause**: Observed entropy higher than all null values
**Solution**: This is unusual - check that sequences are aligned correctly

### Issue: Permutation test very slow

**Cause**: Large alignment or many permutations
**Solutions**:
- Reduce n_permutations to 100
- Process fewer genes
- Use shorter alignments

### Issue: Memory error during permutations

**Cause**: Storing too many permutation results
**Solution**: Don't store all position-level permutations, only summary stats

```python
# Instead of storing full arrays
null_position_entropies.append(perm_entropies)  # Memory-intensive

# Just store summaries
null_means.append(np.mean(perm_entropies))  # Efficient
```
