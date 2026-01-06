# FST_get_summary.R

## Purpose
Bin per-site vcftools FST values into predefined ranges for multiple pairwise comparisons and summarize counts/percentages across bins.【F:scripts/FST_get_summary.R†L8-L39】

## Key inputs
- Working directory set to where the per-site `*.weir.fst` files live (`setwd(...)`).【F:scripts/FST_get_summary.R†L6-L14】
- `files`: Vector of FST file names to include in the summary.【F:scripts/FST_get_summary.R†L8-L14】

## How to run
1. Ensure `data.table` and `ggplot2` are installed (loaded within the script).【F:scripts/FST_get_summary.R†L16-L20】
2. Update `setwd()` and the `files` vector to match your input paths/names.【F:scripts/FST_get_summary.R†L6-L14】
3. Execute from the repository root (or after changing directories):
   ```bash
   Rscript scripts/FST_get_summary.R
   ```
4. The script prints the binned summary table; customize downstream plotting as needed (a grouped bar plot is defined later in the script).

## Notes
- Negative FST values are coerced to zero before binning.【F:scripts/FST_get_summary.R†L28-L33】
- Bins are defined as `<0.05`, `0.05–0.1`, `0.1–0.15`, `0.15-0.25`, and `0.25-1` with percentages computed per comparison.【F:scripts/FST_get_summary.R†L32-L39】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
