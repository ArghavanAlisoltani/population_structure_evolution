# TMRCA_bins_genes_TEs_SNPs_V2_fast_my_bin_selection.py

## Purpose
A vectorized, faster version of the TMRCA binning pipeline that assigns TMRCA values to genes, TEs, and SNPs, then bins frequencies with optional manual bin edges.

## What the script does
- Loads TMRCA segments plus mRNA, TE, and SNP annotations with customizable column names.
- Uses per-scaffold prefix sums to quickly compute average TMRCA for each gene and TE interval.
- Optionally reuses existing SNP TMRCA values or annotates SNPs via per-scaffold merges against segment starts/ends.
- Builds automatic log-scaled bins or uses provided comma-separated edges, then counts segments, SNPs, genes, all TEs, Copia TEs, and Gypsy TEs per bin.
- Writes per-feature TMRCA tables and `tmrca_bins_frequencies.tsv` into the chosen output directory.

## How to run
Provide the required inputs and optional binning choices:
```bash
python TMRCA_bins_genes_TEs_SNPs_V2_fast_my_bin_selection.py \
  --tmrca path/to/tmrca.tsv \
  --mrna path/to/mrna.tsv \
  --te path/to/te.tsv \
  --snps path/to/snps.tsv \
  --outdir tmrca_bins_out_v2_fast
```
Use `--bins` to set manual edges or `--nbins` to alter automatic bin counts.
