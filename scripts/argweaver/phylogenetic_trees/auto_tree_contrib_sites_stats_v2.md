# auto_tree_contrib_sites_stats_v2.R

## Purpose
Auto-select the ARGweaver window containing a query position, identify the contributing tree, and export site/genotype/trait summary tables for one or more replicates.

## Key inputs
- `--smc_dir`: directory containing `.smc` files.
- `--sites_dir`: directory containing matching `.sites` files.
- `--scaffold` and `--position`: coordinates used to select the window.
- `--pheno`: phenotype/metadata table (with ID column matching haplotype labels).
- Optional: `--replicates`, `--trait`, `--outbase`, group column names (`--id_col`, `--proc_col`, `--site_col`, `--mum_col`).

## How to run
```bash
Rscript scripts/argweaver/phylogenetic_trees/auto_tree_contrib_sites_stats_v2.R \
  --smc_dir=smc_files \
  --sites_dir=sites_for_tree \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --trait=C13 \
  --replicates=0,10,20 \
  --outbase=tree_contrib_stats
```

## Notes
- Outputs are written under `<outbase>/<scaffold>_<position>/rep*/` and include selected tree info, contributing site positions, and genotype summaries.
- Requires `data.table` and expects ARGweaver file naming in the `outargs_<scaffold>_<start>_<end>.<rep>.*` format.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
