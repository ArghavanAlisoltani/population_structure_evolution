# auto_tree_contrib_sites_stats_per_position_v3.R guide

## Purpose
Locate the ARGweaver tree interval that spans a single genomic position, then summarize which sites contribute to that tree and produce per-position allele/genotype counts (overall and by metadata groups).

## What the script does
- Parses CLI flags in `--flag=value` format.
- Scans `outargs_<scaffold>_<start>_<end>.<rep>.smc` files to find the window containing the requested position.
- Loads the matching `.smc` and `.sites` files for each requested replicate.
- Reads TREE intervals from the `.smc` file and assigns each site to a tree index.
- Identifies the tree interval that overlaps the query position and writes diagnostics.
- Builds per-position genotype tables for sites that map to the selected tree.
- Writes allele/genotype counts overall and per group (provenance/site/mum) when metadata columns exist.

## Requirements
- R with `data.table` installed.
- ARGweaver `.smc` and `.sites` files following the naming convention:
  - `outargs_<scaffold>_<start>_<end>.<rep>.smc`
  - `outargs_<scaffold>_<start>_<end>.<rep>.sites`

## Key inputs (CLI)
- `--smc_dir`: directory containing `.smc` files.
- `--sites_dir`: directory containing `.sites` files.
- `--scaffold`: scaffold/contig ID to query.
- `--position`: 1-based position to locate within a tree interval.
- `--pheno`: phenotype/metadata table (tab-delimited).
- `--replicates`: comma-separated list of replicate IDs (default: `0`).
- `--outbase`: output base directory (default: `tree_contrib_stats_per_position`).
- `--recursive`: search recursively in input dirs (`true`/`false`, default: `true`).

### Metadata column overrides
- `--id_col` (default: `codg`)
- `--proc_col` (default: `proc`)
- `--site_col` (default: `site`)
- `--mum_col` (default: `mum`)
- `--trait` (default: `C13`, optional; included in long output if present)

### Other flags
- `--write_long`: write per-sample long table (`true`/`false`, default: `true`).

## Outputs
Within `--outbase/<scaffold>_<position>/rep<rep>/`:
- `trees_all_with_site_counts.tsv`: tree intervals with site counts.
- `selected_tree_info.tsv`: interval metadata for the selected tree.
- `selected_tree.nwk`: Newick string for the selected tree.
- `selected_tree_site_positions.tsv`: positions assigned to the selected tree.
- `rep_summary_diagnostics.tsv`: per-replicate diagnostics.
- `genotypes_long_per_position.tsv`: per-sample genotype table (when `--write_long=true`).
- `allele_overall_by_position.tsv`: allele counts per position.
- `genotype_overall_by_position.tsv`: genotype counts per position.
- `het_hom_overall_by_position.tsv`: heterozygous vs homozygous counts.
- `allele_by_provenance_by_position.tsv` + `genotype_by_provenance_by_position.tsv` (if `proc` available).
- `allele_by_site_by_position.tsv` + `genotype_by_site_by_position.tsv` (if `site` available).
- `allele_by_mum_by_position.tsv` + `genotype_by_mum_by_position.tsv` (if `mum` available).

## Run examples
### Example: single replicate (default)
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/auto_tree_contrib_sites_stats_per_position_v3.R \
  --smc_dir=smc_files \
  --sites_dir=sites_files \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --trait=C13 \
  --outbase=tree_contrib_stats
```

### Example: multiple replicates and custom metadata columns
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/auto_tree_contrib_sites_stats_per_position_v3.R \
  --smc_dir=smc_files \
  --sites_dir=sites_files \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --replicates=0,10,20,30 \
  --id_col=sample_id \
  --proc_col=provenance \
  --site_col=collection_site \
  --mum_col=maternal_id \
  --outbase=tree_contrib_stats
```

### Example: disable recursive search and long output
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/auto_tree_contrib_sites_stats_per_position_v3.R \
  --smc_dir=smc_files \
  --sites_dir=sites_files \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --recursive=false \
  --write_long=false
```
