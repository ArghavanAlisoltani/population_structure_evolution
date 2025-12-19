# plot_tmrca_genomewide_ci_v2.py

## Purpose
Plots genome-wide TMRCA segment means and confidence intervals, producing both raw and windowed views with options for bar or line rendering.

## What the script does
- Reads TMRCA segment tables and optionally scaffold length information to scale genomic coordinates.
- Sorts scaffolds, computes windowed statistics (mean/CI) per scaffold, and draws segment rectangles without false links between scaffolds.
- Exports RAW and WINDOWED plots (lines/CI and optional bars) with axis formatting helpers for gigabase positions.

## How to run
```bash
python plot_tmrca_genomewide_ci_v2.py \
  --tmrca path/to/tmrca_segments.tsv \
  --out-prefix genome_tmrca \
  --window-bp 1000000
```
Provide scaffold lengths with `--scaffold-lengths` to scale axes; adjust plotting style via flags such as `--hide-bars` or `--plot-windowed-only`.
