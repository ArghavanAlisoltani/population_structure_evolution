# Visualize_FST_merged_v2.R

## Purpose
Plot three aligned FST comparison tracks over a selected scaffold window using a merged per-site CSV, while also loading overlapping gene annotations for the same region.【F:scripts/Visualize_FST_merged_v2.R†L1-L40】【F:scripts/Visualize_FST_merged_v2.R†L54-L79】

## Key inputs
- Working directory with `fst_per_site_merged_rowmax_lt0.15.csv` and annotation TSVs already merged (`All_fst_merged_with_anno.tsv`).【F:scripts/Visualize_FST_merged_v2.R†L6-L10】
- `in_csv`: Merged FST table; `comparisons`: three comparison column names to facet; region bounds (`scaffold`, `start_bp`, `end_bp`).【F:scripts/Visualize_FST_merged_v2.R†L11-L35】
- Annotation workbook path (`anno_xlsx`) to overlay gene intervals for the same scaffold/region.【F:scripts/Visualize_FST_merged_v2.R†L54-L69】

## How to run
1. Install `tidyverse`, `data.table`, and `readxl` packages.【F:scripts/Visualize_FST_merged_v2.R†L4-L8】【F:scripts/Visualize_FST_merged_v2.R†L55-L61】
2. Edit the input block to point to your merged CSV, comparison names, and target scaffold coordinates; adjust `setwd()` if needed.【F:scripts/Visualize_FST_merged_v2.R†L6-L23】
3. Execute:
   ```bash
   Rscript scripts/Visualize_FST_merged_v2.R
   ```
4. The script verifies the comparison columns, filters the region, reshapes to long format, and facets the three tracks with shared y-limits; gene intervals in the same window are prepared for overlay.【F:scripts/Visualize_FST_merged_v2.R†L25-L46】【F:scripts/Visualize_FST_merged_v2.R†L54-L79】

## Notes
- Negative or missing FST values are coerced to zero prior to plotting.【F:scripts/Visualize_FST_merged_v2.R†L36-L44】
- Overlapping feature counts are exported to `fst_415_overlapped.tsv` before plotting.【F:scripts/Visualize_FST_merged_v2.R†L7-L10】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
