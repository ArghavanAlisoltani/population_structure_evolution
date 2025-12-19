# TMRCA_bins_genes_TEs_SNPs_V1.py

## Purpose
Annotates genes, TEs, and SNPs with average TMRCA values from segment data, then bins those values to summarize how often each feature type appears in each TMRCA range.

## What the script does
- Loads TMRCA segments plus mRNA, TE, and SNP annotations, with configurable column names.
- Assigns average TMRCA to each gene and TE by overlapping intervals, and optionally annotates SNPs by merging against segment starts.
- Builds either automatic log-scaled bins or user-specified bin edges.
- Counts how many segments, SNPs, genes, TEs, Copia elements, and Gypsy elements fall into each bin and writes per-feature tables and `tmrca_bins_frequencies.tsv` to the output directory.

## How to run
Supply the required inputs and optional binning parameters:
```bash
python TMRCA_bins_genes_TEs_SNPs_V1.py \
  --tmrca path/to/tmrca.tsv \
  --mrna path/to/mrna.tsv \
  --te path/to/te.tsv \
  --snps path/to/snps.tsv \
  --outdir tmrca_bins_out_v1
```
Add `--bins` for custom bin edges (comma-separated) or `--nbins` to change automatic bin counts.
