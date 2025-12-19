# Greedy_sample_selection.R

## Purpose
Select a target number of individuals that are maximally separated in PCA space using a greedy max–min heuristic followed by 1-swap polishing, and report the achieved minimum pairwise distance (radius).【F:scripts/Greedy_sample_selection.R†L1-L59】

## Key inputs
- User settings: `k_target`, `pc_rank`, `min_sep_q`, `relax`, and `seed` controlling selection size, PCA rank, and adaptive separation thresholds.【F:scripts/Greedy_sample_selection.R†L5-L15】【F:scripts/Greedy_sample_selection.R†L30-L45】
- Genotype matrix loaded from `../../snp_1490x25099_CharlesSEP2025.txt` (edit to your SNP table) and optional metadata CSV for labeling selections.【F:scripts/Greedy_sample_selection.R†L12-L33】【F:scripts/Greedy_sample_selection.R†L63-L76】

## How to run
1. Install required packages such as `mixOmics` (for PCA used later) and base R dependencies (`data.table` not needed; uses base functions).【F:scripts/Greedy_sample_selection.R†L63-L70】
2. Edit the working directory (`setwd()`), user settings, and input file paths near the top of the script to match your datasets.【F:scripts/Greedy_sample_selection.R†L5-L18】【F:scripts/Greedy_sample_selection.R†L12-L33】
3. Execute:
   ```bash
   Rscript scripts/Greedy_sample_selection.R
   ```
4. The script computes PC distances, runs greedy selection plus swap polishing, prints the selected ID count and minimum distance, and tags selected IDs in the merged metadata table.【F:scripts/Greedy_sample_selection.R†L30-L45】【F:scripts/Greedy_sample_selection.R†L46-L59】【F:scripts/Greedy_sample_selection.R†L63-L76】

## Notes
- `greedy_maxmin_adaptive()` seeds with the most isolated sample and relaxes the minimum separation if no candidates meet the threshold; `polish_1swap_maxmin()` iteratively improves the minimum distance via single swaps.【F:scripts/Greedy_sample_selection.R†L30-L59】
- Update the PCA rank or separation quantile to tune diversity; the script expects rownames on the genotype matrix to serve as IDs and echoes subsetting progress at each step.【F:scripts/Greedy_sample_selection.R†L12-L33】
