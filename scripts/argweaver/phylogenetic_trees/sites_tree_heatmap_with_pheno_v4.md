# sites_tree_heatmap_with_pheno_v4.R run guide

Builds a neighbor-joining tree from an ARGweaver `.sites` file, joins phenotype data, and plots allele/phenotype heatmaps alongside the tree.

## Usage
```bash
Rscript sites_tree_heatmap_with_pheno_v4.R --sites <FILE.sites> --pheno <PHENO.txt> --positions <POS_LIST> [options]
```

## Options
- `--sites` (required): ARGweaver `.sites` file containing sequences and positions.
- `--pheno` (required): Phenotype table with `codg`, `proc`, `site`, and `mum` columns.
- `--positions` (required): Comma-separated genomic positions to extract alleles for.
- `--trait`: Phenotype column to visualize (default: `C13`).
- `--outprefix`: Output prefix for tree/plot files (default: `sites_tree_pheno`).
- `--hide_labels`: Hide tree tip labels (`true`/`false`; default: `true`).
- `--tree_scale`: Tree scaling mode `distance` or `cladogram` (default: `distance`).

The script exits with a usage message unless the sites file, phenotype file, and allele positions are provided.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
