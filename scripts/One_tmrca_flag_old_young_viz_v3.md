# One_tmrca_flag_old_young_viz_v3.R

## Purpose
Visualize a single ARGweaver TMRCA file over a specified scaffold interval, flagging the youngest and oldest percentile windows and overlaying gene annotations from an Excel workbook.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L1-L33】

## Key inputs
- `tmrca_path`: Path to the TMRCA track; `ann_xlsx`: Annotation workbook with scaffold coordinates and gene IDs.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L12-L24】【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L50-L65】
- Region bounds (`contig`, `x_start`, `x_end`), window size (`win`), statistic (`tmrca_stat`), and percentile thresholds (`bottom_q`, `top_q`).【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L18-L33】

## How to run
1. Install the required packages: `data.table`, `dplyr`, `ggplot2`, `scales`, `readxl`, `ggrepel`.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L5-L11】
2. Update `setwd()` and the **USER SETTINGS** block with your file paths, scaffold, region, and window/percentile choices.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L12-L33】
3. Execute the script:
   ```bash
   Rscript scripts/One_tmrca_flag_old_young_viz_v3.R
   ```
4. The script reads the TMRCA track (handling 6- vs 7-column formats), adjusts zero-based coordinates, windows the region, and plots the flagged young/old intervals with optional gene labels; output PNG is named from the scaffold and coordinates.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L34-L48】【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L18-L33】

## Notes
- `read_annotation_1ab()` maps the annotation workbook to standardized columns and provides default gene IDs when missing.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L50-L65】
- Gene labels can be capped via `max_gene_rows` to avoid clutter on dense regions.【F:scripts/One_tmrca_flag_old_young_viz_v3.R†L27-L30】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
