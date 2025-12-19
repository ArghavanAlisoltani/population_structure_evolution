# Summary_100_greedy_samples.R

## Purpose
Summarize per-site vcftools FST files by binning values across multiple comparisons, generate a stacked bar plot of bin percentages, and include utilities for reshaping pairwise FST outputs into a lower-triangle table.【F:scripts/Summary_100_greedy_samples.R†L8-L61】【F:scripts/Summary_100_greedy_samples.R†L64-L80】

## Key inputs
- Working directory containing the `*.weir.fst` files (`setwd(...)`).【F:scripts/Summary_100_greedy_samples.R†L6-L17】
- `files`: Vector of FST file names for the comparisons to summarize.【F:scripts/Summary_100_greedy_samples.R†L12-L17】

## How to run
1. Install required R packages: `data.table`, `ggplot2`, `dplyr`, `readr`, `stringr`, `purrr`, and `tidyr` (loaded within the script).【F:scripts/Summary_100_greedy_samples.R†L1-L5】【F:scripts/Summary_100_greedy_samples.R†L64-L73】
2. Adjust `setwd()` and the `files` list to match your environment and comparisons.【F:scripts/Summary_100_greedy_samples.R†L6-L17】
3. Execute:
   ```bash
   Rscript scripts/Summary_100_greedy_samples.R
   ```
4. The script bins FST values, writes a stacked bar plot (`fst_bin_summary_procs345_R.png`), and can reshape pairwise FST files into a lower-triangular matrix for reporting.【F:scripts/Summary_100_greedy_samples.R†L50-L61】【F:scripts/Summary_100_greedy_samples.R†L64-L80】

## Notes
- Negative FST values are set to zero before binning to avoid artifacts.【F:scripts/Summary_100_greedy_samples.R†L28-L33】
- Bin edges mirror those used in other summaries (`<0.05`, `0.05–0.1`, `0.1–0.15`, `0.15-0.25`, `0.25-1`).【F:scripts/Summary_100_greedy_samples.R†L34-L39】
