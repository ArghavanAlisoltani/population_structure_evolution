# concat_tmrca.R

## Purpose
Concatenate ARGweaver TMRCA outputs across scaffolds, export high-TMRCA subsets, and visualize TMRCA distributions by scaffold via violin plots; includes helpers for gene overlap mapping.【F:scripts/concat_tmrca.R†L1-L20】【F:scripts/concat_tmrca.R†L30-L66】【F:scripts/concat_tmrca.R†L70-L77】

## Key inputs
- Working directory pointing to the folder with `All_tmrca_60_run.txt` and any `*.tmrca.txt` files to merge.【F:scripts/concat_tmrca.R†L9-L33】
- Annotation workbook for downstream overlap mapping (see helper section).【F:scripts/concat_tmrca.R†L70-L77】

## How to run
1. Install required packages: `data.table`, `readxl`, `dplyr`, `stringr`, `GenomicRanges`, `ggplot2`, `scales`, and `tools`.【F:scripts/concat_tmrca.R†L1-L8】
2. Edit `setwd()` and input file names (e.g., `All_tmrca_60_run.txt`) to point to your TMRCA exports.【F:scripts/concat_tmrca.R†L9-L18】【F:scripts/concat_tmrca.R†L30-L33】
3. Execute:
   ```bash
   Rscript scripts/concat_tmrca.R
   ```
4. The script filters high-mean TMRCA segments, saves scaffold-specific subsets, builds a combined table from all `*.tmrca.txt` files, and produces log-scale violin plots saved as PNG/PDF.【F:scripts/concat_tmrca.R†L16-L19】【F:scripts/concat_tmrca.R†L30-L66】

## Notes
- Column names are harmonized whether the TMRCA files include a mode column or not; scaffolds are reordered by median TMRCA for clearer plotting.【F:scripts/concat_tmrca.R†L30-L50】
- Violin plots mark medians (black) and means (red) and rotate scaffold labels for readability.【F:scripts/concat_tmrca.R†L51-L66】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
