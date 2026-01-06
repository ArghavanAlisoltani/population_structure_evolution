# Glob_tmrca_v2.R

## Purpose
Process multiple ARGweaver TMRCA files genome-wide: harmonize columns, convert zero-based starts, window the segments, flag youngest/oldest quantiles, and plot aggregated distributions with confidence intervals.【F:scripts/Glob_tmrca_v2.R†L1-L40】【F:scripts/Glob_tmrca_v2.R†L59-L75】

## Key inputs
- Working directory set to the folder containing the TMRCA track files (patterned by `glob`).【F:scripts/Glob_tmrca_v2.R†L5-L16】【F:scripts/Glob_tmrca_v2.R†L59-L63】
- `glob`: Filename pattern (e.g., `outargs_scaffold_10_*.tmrca.txt`) to gather per-scaffold TMRCA outputs.【F:scripts/Glob_tmrca_v2.R†L13-L20】【F:scripts/Glob_tmrca_v2.R†L59-L63】
- Window size, statistic choice, and percentile thresholds to classify young/old regions (`win`, `tmrca_stat`, `bottom_q`, `top_q`).【F:scripts/Glob_tmrca_v2.R†L14-L20】【F:scripts/Glob_tmrca_v2.R†L64-L70】

## How to run
1. Install `data.table`, `dplyr`, `ggplot2`, and `scales`.【F:scripts/Glob_tmrca_v2.R†L7-L12】
2. Update `setwd()` and the settings block (`glob`, window/stat parameters, output names) for your dataset.【F:scripts/Glob_tmrca_v2.R†L5-L20】
3. Execute:
   ```bash
   Rscript scripts/Glob_tmrca_v2.R
   ```
4. The script windows all matching files, summarizes per-window means and CI widths, flags percentile extremes, and writes both PNG and PDF plots named from the configured outputs.【F:scripts/Glob_tmrca_v2.R†L20-L22】【F:scripts/Glob_tmrca_v2.R†L73-L80】

## Notes
- `read_tm()` auto-detects 6- vs 7-column TMRCA inputs and converts zero-based starts when required.【F:scripts/Glob_tmrca_v2.R†L23-L38】
- `split_to_windows()` expands segment coverage across fixed windows, weighting by covered bases to derive mean TMRCA and CI widths.【F:scripts/Glob_tmrca_v2.R†L40-L57】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
