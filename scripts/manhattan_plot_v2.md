# manhattan_plot_v2.R

## Purpose
Create Manhattan-style plots from tab-separated FST results by computing cumulative scaffold positions, cleaning negative/NaN values, and alternating point colors across scaffolds.【F:scripts/manhattan_plot_v2.R†L1-L38】

## Key inputs
- Input TSV with three columns: scaffold ID, position, and FST value (first two columns renamed to `CHROM` and `POS`).【F:scripts/manhattan_plot_v2.R†L8-L23】
- Plot title and file path supplied when calling `create_manhattan_plot(filepath, title)`.【F:scripts/manhattan_plot_v2.R†L5-L10】

## How to run
1. Install `tidyverse` and `patchwork` (loaded at the top).【F:scripts/manhattan_plot_v2.R†L1-L4】
2. From an interactive R session or via a small wrapper, source the script and call the helper function:
   ```r
   source("scripts/manhattan_plot_v2.R")
   plt <- create_manhattan_plot("path/to/fst.tsv", "My title")
   ggsave("fst_manhattan.png", plt, width = 10, height = 4, dpi = 300)
   ```
3. The function cleans negative FST values to zero, computes cumulative base-pair positions per scaffold, and returns a ggplot object ready for saving.【F:scripts/manhattan_plot_v2.R†L12-L39】

## Notes
- Scaffolds are ordered numerically based on digits in the ID (`CHROM_NUM`), enabling consistent alternating colors across chromosomes.【F:scripts/manhattan_plot_v2.R†L18-L33】【F:scripts/manhattan_plot_v2.R†L35-L39】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
