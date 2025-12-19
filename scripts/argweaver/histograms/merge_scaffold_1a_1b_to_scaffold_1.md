# merge_scaffold_1a_1b_to_scaffold_1.py

## Purpose
Merges separate scaffold_1a and scaffold_1b coordinates back into a single scaffold_1 by applying an offset to 1b coordinates and renaming both to a unified scaffold name.

## What the script does
- Reads a TSV with chromosome, start, end, and optional segment length columns (auto-detected).
- Adds the provided offset to start/end positions for rows matching the 1b scaffold, renames both 1a and 1b to the output scaffold name, and recomputes segment length when present.
- Prints a brief range summary for the merged scaffold and writes the adjusted table to the specified output file.

## How to run
Supply the input file, output file, and offset where scaffold_1b begins on the merged scaffold:
```bash
python merge_scaffold_1a_1b_to_scaffold_1.py \
  --in path/to/input.tsv \
  --out merged_scaffold1.tsv \
  --offset 100000000 \
  --scaf1a scaffold_1a --scaf1b scaffold_1b --scaf-out scaffold_1
```
Use `--offset-is-1based` if the offset provided refers to a 1-based coordinate.
