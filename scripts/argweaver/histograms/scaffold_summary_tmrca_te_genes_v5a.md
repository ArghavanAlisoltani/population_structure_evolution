# scaffold_summary_tmrca_te_genes_v5a.py

## Purpose
Summarizes TMRCA, TE, and gene features per scaffold, producing scaffold-level metrics and overlap statistics.

## What the script does
- Loads TMRCA segment, TE annotation, gene annotation, and optional SNP tables.
- Validates required columns, sanitizes column names, and computes interval unions, quantiles, and coverage metrics per scaffold.
- Reports scaffold-level summaries (e.g., segment counts, TE/gene coverage, quantiles) to TSV files in the output directory.

## How to run
```bash
python scaffold_summary_tmrca_te_genes_v5a.py \
  --tmrca path/to/tmrca.tsv \
  --tes path/to/te_annotations.tsv \
  --genes path/to/gene_annotations.tsv \
  --outdir scaffold_summary_out
```
Add optional inputs (e.g., `--snps`) or change separators as needed; the script creates the output directory if it does not exist.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
