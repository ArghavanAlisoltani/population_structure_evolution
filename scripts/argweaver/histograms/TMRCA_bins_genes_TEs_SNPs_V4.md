# TMRCA_bins_genes_TEs_SNPs_V4.py

## Purpose
Computes gene-, TE-, and SNP-level TMRCA assignments and bins them, with options for manual bin edges and deduplicating overlapping TE intervals before summarizing.

## What the script does
- Loads TMRCA segments, SNP-level TMRCA, gene annotations, and TE annotations.
- Optionally collapses overlapping TEs by keeping the longest element per scaffold when `--dedup-overlapping-tes` is enabled.
- Assigns average TMRCA to genes and TEs, defines automatic or manual bins, and tallies segments, SNPs, genes, all TEs, Copia TEs, and Gypsy TEs per bin.
- Outputs a single `tmrca_bins_frequencies.tsv` with bin counts and fractions for each feature class.

## How to run
Run with the required inputs and optional flags:
```bash
python TMRCA_bins_genes_TEs_SNPs_V4.py \
  --tmrca path/to/tmrca.tsv \
  --snps path/to/snps_with_tmrca.tsv \
  --genes path/to/gene_annotations.tsv \
  --tes path/to/te_annotations.tsv \
  --outdir tmrca_bins_out_v4
```
Use `--manual-bins` for custom edges or `--dedup-overlapping-tes` to collapse overlapping TEs.
