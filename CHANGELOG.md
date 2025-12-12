# Changelog

All notable changes to Entropia-MSA will be documented in this file.

## [1.1.0] - 2025-11-20

### Added
- **Boxplot visualization** (`plot_heatmap_with_entropy_boxplot.py`)
  - Replaces barplot with boxplot showing entropy distribution across alignment positions
  - Reveals bimodal distributions (conserved domains vs disordered regions)
  - Shows quartiles, spread, and variability per protein
  - Includes permutation-based null model line for visual comparison

- **Permutation testing** (`permutation_test_entropy.py`)
  - Standalone tool for statistical significance testing
  - Establishes empirical null distributions via column shuffling
  - Calculates p-values for observed vs null entropy
  - Addresses biological reality that entropy rarely reaches 1.0
  - Identifies proteins under purifying selection

- **Comprehensive documentation** (`docs/BOXPLOT_AND_PERMUTATIONS.md`)
  - Explains why boxplots are superior to barplots
  - Details permutation testing methodology
  - Provides interpretation guidelines
  - Includes troubleshooting and performance tips

### Changed
- Original barplot version retained as `plot_heatmap_with_entropy.py`
- Both visualization options now available

### Improved
- Better statistical rigor with significance testing
- Enhanced ability to detect domain architecture
- More informative visualizations of entropy distribution

## [1.0.0] - 2025-11-20

### Initial Release

#### Core Features
- **Shannon entropy calculation** for amino acid MSAs
- **Gap trimming** (removes columns with >80% gaps)
- **Normalization** by log₂(20) = 4.32 for 20 amino acids
- **Positional entropy profiles** (multi-page PDF output)
- **Integrated heatmap** combining entropy with gene copy numbers

#### Tools Included
- `calculate_shannon_entropy.py` - Batch entropy calculation
- `plot_positional_entropy.py` - Detailed positional profiles
- `plot_heatmap_with_entropy.py` - Integrated visualization

#### Documentation
- Comprehensive README with installation and usage
- Detailed USAGE.md guide
- Example alignments included

#### Outputs
- CSV summary statistics
- High-quality PNG/PDF/SVG visualizations
- Multi-page positional entropy PDFs

#### Analysis Developed For
- Rhynchospora phylogenomics project
- Kinetochore protein evolution study
- 182 orthogroups across 32 species/haplotypes

---

## Upcoming Features (Planned)

### v1.2.0
- [ ] Sequence divergence analysis
- [ ] MSA overview plots
- [ ] Row-wise (sequence-wise) entropy metrics
- [ ] Parallelization for large datasets

### v1.3.0
- [ ] Interactive HTML visualizations
- [ ] Automated domain boundary detection
- [ ] Integration with disorder predictors
- [ ] Coevolution analysis

### v2.0.0
- [ ] Nucleotide sequence support (normalized by log₂(4))
- [ ] Codon-aware analysis
- [ ] dN/dS integration
- [ ] Web interface

---

## Version Numbering

- **Major (X.0.0)**: Breaking changes, new architecture
- **Minor (1.X.0)**: New features, backwards compatible
- **Patch (1.0.X)**: Bug fixes, documentation updates

---

## Contributors

- Javier González (@jacgonisa) - Lead developer

---

## Citation

If you use Entropia-MSA in your research, please cite:

```
González, J. (2025). Entropia-MSA: Shannon Entropy Analysis for Multiple Sequence Alignments.
GitHub: https://github.com/jacgonisa/entropia-msa
Version 1.1.0
```

---

## License

MIT License - See LICENSE file for details
