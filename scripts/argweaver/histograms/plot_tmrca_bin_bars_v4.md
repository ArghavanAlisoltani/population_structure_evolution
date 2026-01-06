# plot_tmrca_bin_bars_v4.py

## Purpose
Creates dual-axis TMRCA bin plots with smoothed area overlays, optional y-axis breaks, and support for using upper bin edges as x-axis labels.

## What the script does
- Reads a `tmrca_bins_frequencies.tsv` table and selects the specified bin label column.
- Draws overlapping bar charts for a main series (e.g., SNPs or segments) and secondary series (genes, TEs, Copia, Gypsy), optionally adding smoothed area curves.
- Supports y-axis breaks, configurable legend placement, and smoothing windows for area plots.
- Writes PNG figures with filenames based on `--out-prefix` and plot suffixes.

## How to run
```bash
python plot_tmrca_bin_bars_v4.py \
  --input tmrca_bins_out_v4/tmrca_bins_frequencies.tsv \
  --out-prefix tmrca_bins_out_v4/tmrca_bins_barplots_v4_ \
  --col-bin bin_right \
  --col-snp-frac frac_snps \
  --col-seg-frac frac_segments \
  --col-gene-frac frac_genes \
  --col-te-frac frac_TEs_all \
  --col-copia-frac frac_TEs_Copia \
  --col-gypsy-frac frac_TEs_Gypsy \
  --make-area \
  --smooth-area-window 3 \
  --legend-loc right-out
```
Update column mappings and styling flags to match your binning output.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
