# plot_tmrca_bin_bars_v1.py

## Purpose
Creates overlapping bar plots from `tmrca_bins_frequencies.tsv`, comparing SNP or segment counts against gene/TE and TE subtype counts across TMRCA bins.

## What the script does
- Loads the bin frequency table produced by the TMRCA binning scripts.
- Builds paired bar charts for multiple feature comparisons (e.g., SNPs vs mRNA+TE, SNPs vs Copia+Gypsy, segments vs TE classes) using dual y-axes.
- Saves a series of PNG files with an `--out-prefix` for downstream visualization.

## How to run
```bash
python plot_tmrca_bin_bars_v1.py \
  --bins-file tmrca_bins_out_v3/tmrca_bins_frequencies.tsv \
  --out-prefix tmrca_bins_plots
```
Customize colors, font sizes, and figure sizes with the provided optional arguments.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
