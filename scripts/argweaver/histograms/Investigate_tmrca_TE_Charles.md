# Investigate_tmrca_TE_Charles.R

## Purpose
Quickly inspects a specific TE on scaffold_1b and the surrounding TMRCA values, writing a subset of the data for review.

## What the script does
- Loads TE annotations and counts TE classes present in the dataset.
- Extracts two TE entries on `scaffold_1b` at the positions of interest and computes their lengths.
- Reads TMRCA segments, filters the region spanning 447,420,110–447,437,260 on scaffold_1b, and reports mean/median TMRCA.
- Loads per-TE TMRCA summaries, subsets the same region, and saves the filtered segment data to `tmrca_data_for_Charles/TMRCA_scaffold1b_447420110_447437260.txt`.

## How to run
1. Update the `setwd()` call and the file paths near the top of the script to match your local data locations.
2. Execute with R:
   ```bash
   Rscript Investigate_tmrca_TE_Charles.R
   ```
   The script expects the TE annotation TSV, TMRCA segment TSV, and TE-level TMRCA TSV to exist at the configured paths.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
