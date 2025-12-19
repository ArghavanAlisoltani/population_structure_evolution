# merge_FST_annotations.R

## Purpose
Merge per-site FST results with nearest annotation metadata, extract top signals for each comparison, and write a combined TSV for downstream visualization.【F:scripts/merge_FST_annotations.R†L1-L29】

## Key inputs
- Working directory set to the folder containing `fst_per_site_merged_rowmax_lt0.15.csv` and pairwise `PROC*_vs_*.weir.fst` files.【F:scripts/merge_FST_annotations.R†L4-L19】
- Annotation table path (`anno`) with scaffold and position columns to join against FST sites.【F:scripts/merge_FST_annotations.R†L15-L28】

## How to run
1. Install dependencies: `readr`, `readxl`, `dplyr`, `tibble`, `GenomicRanges`, and `data.table`.【F:scripts/merge_FST_annotations.R†L1-L13】
2. Update `setwd()`, input filenames, and annotation path near the top of the script.【F:scripts/merge_FST_annotations.R†L4-L28】
3. Execute:
   ```bash
   Rscript scripts/merge_FST_annotations.R
   ```
4. The script reads FST comparisons, filters high-scoring windows (e.g., ≥0.25), merges annotation by scaffold-position key, and writes `All_fst_merged_with_anno.tsv`.【F:scripts/merge_FST_annotations.R†L15-L28】

## Notes
- Edit the FST thresholds or columns used to define “top” windows to match your study; the current value is 0.25 for each comparison.【F:scripts/merge_FST_annotations.R†L20-L24】
- Ensure annotation columns are named consistently (`X.CHROM`, `POS`, etc.) or adjust the `paste`/`merge` logic accordingly.【F:scripts/merge_FST_annotations.R†L25-L28】
