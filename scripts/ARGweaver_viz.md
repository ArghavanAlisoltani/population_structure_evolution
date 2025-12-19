# ARGweaver_viz.R

## Purpose
Visualize a single ARGweaver TMRCA track for one scaffold region, optionally flagging young/old windows and labeling genes from an Excel annotation workbook. The script windows the TMRCA values, derives percentiles, and exports a PNG for the selected interval.

## Key inputs
- `tmrca_path`: Path to the TMRCA file with mean/mode/CI columns (ARGweaver output).【F:scripts/ARGweaver_viz.R†L12-L24】
- `ann_xlsx`: Excel workbook containing gene coordinates in the 1a/1b scaffold system.【F:scripts/ARGweaver_viz.R†L12-L24】
- Region parameters (`contig`, `x_start`, `x_end`) and window/statistic options (`win`, `tmrca_stat`, percentile cutoffs).【F:scripts/ARGweaver_viz.R†L16-L30】

## How to run
1. Install the required R packages (`data.table`, `dplyr`, `ggplot2`, `scales`, `readxl`, `ggrepel`).【F:scripts/ARGweaver_viz.R†L4-L10】
2. Edit the **USER SETTINGS** block near the top to point to your TMRCA file, annotation workbook, scaffold name, region coordinates, and window/percentile choices.【F:scripts/ARGweaver_viz.R†L12-L30】
3. From the repository root (or after updating any `setwd` you add), run:
   ```bash
   Rscript scripts/ARGweaver_viz.R
   ```
4. The script writes a region-specific PNG named from the scaffold, metric, and coordinates (e.g., `scaffold_10_tmrca_region_500000000_510000000_.png`).【F:scripts/ARGweaver_viz.R†L30-L32】

## Notes
- The helper `read_tm_one()` auto-detects 6- vs 7-column ARGweaver TMRCA files and adjusts zero-based starts when requested.【F:scripts/ARGweaver_viz.R†L33-L40】
- Percentile thresholds flag “young” and “old” windows based on the chosen statistic; gene labeling can be limited via `max_gene_rows`.【F:scripts/ARGweaver_viz.R†L22-L29】
