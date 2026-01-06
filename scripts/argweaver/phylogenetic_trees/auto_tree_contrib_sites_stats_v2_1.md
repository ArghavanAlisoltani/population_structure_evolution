# auto_tree_contrib_sites_stats_v2_1.R

## Purpose
A variant of the contributing-site summary script that auto-selects an ARGweaver window and reports the sites/genotypes supporting the selected tree interval for one or more replicates.

## Key inputs
- `--smc_dir`: directory containing `.smc` files.
- `--sites_dir`: directory containing `.sites` files.
- `--scaffold` and `--position`: query coordinates for selecting the window.
- `--pheno`: phenotype/metadata table.
- Optional: `--replicates`, `--trait`, `--outbase`, and metadata column names (`--id_col`, `--proc_col`, `--site_col`, `--mum_col`).

## How to run
```bash
Rscript scripts/argweaver/phylogenetic_trees/auto_tree_contrib_sites_stats_v2_1.R \
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
- Outputs land in `<outbase>/<scaffold>_<position>/rep*/` with tree, site, and genotype summaries.
- Designed for ARGweaver outputs named `outargs_<scaffold>_<start>_<end>.<rep>.*`.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
