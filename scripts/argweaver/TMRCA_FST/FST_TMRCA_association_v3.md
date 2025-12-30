# FST_TMRCA_association_v3.R

## Purpose
Merge per-site FST values with ARGweaver TMRCA segments, then generate correlation/regression plots and filtered output tables to explore FST–TMRCA relationships.

## Key inputs
- `fst_file`: per-site FST TSV (expects `Scaffold`/`POS` plus contrast columns).
- `tmrca_file`: TMRCA segment TSV (expects `CHROM`, `start_tmrca`, `end_tmrca`, `mean_tmrca`).
- Output paths: `out_merged_tsv`, `out_plot_dir`, and filter options near the middle of the script.

## How to run
1. Update `setwd()` and the `fst_file`/`tmrca_file` paths at the top of the script.
2. Run from R or via the command line:
   ```bash
   Rscript scripts/argweaver/TMRCA_FST/FST_TMRCA_association_v3.R
   ```

## Notes
- The merge uses interval overlap so each SNP is mapped to its containing TMRCA segment.
- Filtering affects the plots and filtered tables, while the merged TSV is preserved unchanged.
