# auto_tree_informative_sites_stats_v1.R

## Purpose
Auto-select an ARGweaver tree interval for a query position and compute summary statistics for informative sites, including genotype and trait summaries.

## Key inputs
- `--smc_dir`: directory containing `.smc` files.
- `--sites_dir`: directory containing `.sites` files.
- `--scaffold` and `--position`: coordinates used to select the window.
- `--pheno`: phenotype/metadata table with a matching ID column.
- Optional: `--replicates`, `--trait`, `--outbase`, and column names (`--id_col`, `--proc_col`, `--site_col`, `--mum_col`).

## How to run
```bash
Rscript scripts/argweaver/phylogenetic_trees/auto_tree_informative_sites_stats_v1.R \
  --smc_dir=smc_files \
  --sites_dir=sites_for_tree \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --trait=C13 \
  --replicates=0,10,20 \
  --outbase=tree_stats_out
```

## Notes
- Produces per-position and pooled genotype/allele summaries under `<outbase>/<scaffold>_<position>/rep*/`.
- Requires `data.table` and ARGweaver file naming of `outargs_<scaffold>_<start>_<end>.<rep>.*`.
