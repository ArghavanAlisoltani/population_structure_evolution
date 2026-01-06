# TMRCA_bins_genes_TEs_SNPs_V3_automatic_bin.py

## Purpose
Automates TMRCA binning with the same vectorized overlap logic as V2 while keeping automatic or user-supplied bin edges for segment, SNP, gene, and TE summaries.

## What the script does
- Parses TMRCA segments, gene annotations, TE annotations, and SNP tables with flexible column names.
- Uses per-scaffold prefix sums to assign average TMRCA to each gene and TE interval and annotates SNPs via per-scaffold merges when needed.
- Builds automatic log-scaled bins (or accepts explicit edges) and counts how many segments, SNPs, genes, all TEs, Copia TEs, and Gypsy TEs fall into each bin.
- Writes the per-feature TMRCA outputs and consolidated `tmrca_bins_frequencies.tsv` to the chosen directory.

## How to run
Example invocation:
```bash
python TMRCA_bins_genes_TEs_SNPs_V3_automatic_bin.py \
  --tmrca path/to/tmrca.tsv \
  --mrna path/to/mrna.tsv \
  --te path/to/te.tsv \
  --snps path/to/snps.tsv \
  --outdir tmrca_bins_out_v3
```
Add `--bins` or `--nbins` to control binning behavior.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
