# Position_reassign_bb_to_b.R

## Purpose
Normalizes scaffold names and coordinates by shifting selected `*b` scaffolds forward by 1,000,000,000 bp and trimming the trailing `b` from their names.

## What the script does
- Reads `all.tmrca.txt` and labels columns as chromosome, TMRCA start/end, mean, and confidence intervals.
- Forces start/end coordinates to numeric values for safety.
- Adds 1,000,000,000 to `start_tmrca` and `end_tmrca` for rows whose `CHROM` matches a predefined list of scaffolds ending in `b` or `bb`.
- Drops the final `b` in the scaffold name for those rows, recomputes segment length, and writes the corrected table to `all_tmrca_corrected_position.tsv`.

## How to run
1. Update the working directory and input filename at the top of the script if needed.
2. Execute with R:
   ```bash
   Rscript Position_reassign_bb_to_b.R
   ```
   The script expects `all.tmrca.txt` in the working directory and writes the corrected TSV alongside it.
