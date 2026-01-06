# FST_visualization.R

## Purpose
Plot three aligned per-site FST tracks for a chosen scaffold region, using a merged CSV of FST comparisons and optional annotation overlap summaries.【F:scripts/FST_visualization.R†L1-L21】【F:scripts/FST_visualization.R†L48-L76】

## Key inputs
- Working directory pointing to the folder with `fst_per_site_merged_rowmax_lt0.15.csv` and annotation TSV files.【F:scripts/FST_visualization.R†L5-L9】
- `in_csv`: Merged per-site FST table with columns `CHROM`, `POS`, and the three comparison columns listed in `comparisons`.【F:scripts/FST_visualization.R†L10-L23】
- Region bounds (`scaffold`, `start_bp`, `end_bp`) to subset and plot.【F:scripts/FST_visualization.R†L10-L23】

## How to run
1. Install required R packages: `tidyverse`, `readr`, and `readxl` (loaded in the script).【F:scripts/FST_visualization.R†L4-L9】【F:scripts/FST_visualization.R†L48-L56】
2. Update the file paths and region parameters in the input block at the top of the script.【F:scripts/FST_visualization.R†L10-L23】
3. Execute:
   ```bash
   Rscript scripts/FST_visualization.R
   ```
4. The script filters the data to the requested scaffold window, reshapes to long format, and facets the three comparisons with consistent y-limits; annotation genes overlapping the region are also loaded for plotting.【F:scripts/FST_visualization.R†L35-L45】【F:scripts/FST_visualization.R†L52-L79】

## Notes
- The script guards against missing comparison columns and coerces negative/NA FST values to zero before plotting.【F:scripts/FST_visualization.R†L24-L45】
- Overlap counts are exported to `fst_415_overlapped.tsv` prior to plotting, based on the `overlaps_feature` flag in the merged annotation table.【F:scripts/FST_visualization.R†L6-L9】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
