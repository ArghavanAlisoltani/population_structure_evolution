# FST_staritified_meanTMRCA_v1.R

## Purpose
Analyze joint FST/TMRCA patterns by scanning for outliers and stratifying FST distributions across TMRCA bins.

## Key inputs
- `merged_file`: TSV produced by the merge step (expects `Scaffold`, `POS`, `mean_tmrca`, and FST contrast columns).
- Output directory: `out_dir`.
- Filter and plotting options in the `global_filter`, `opt1`, and `opt3` blocks.

## How to run
1. Update `setwd()` and `merged_file` to point to your merged FST+TMRCA table.
2. Execute:
   ```bash
   Rscript scripts/argweaver/TMRCA_FST/FST_staritified_meanTMRCA_v1.R
   ```

## Notes
- Option 1 generates quadrant/outlier plots; Option 3 bins TMRCA and compares FST distributions.
- The script expects contrast columns like `PROC3_vs_PROC4`; adjust the `contrasts` vector if needed.
