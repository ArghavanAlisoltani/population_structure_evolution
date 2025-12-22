# plot_flanks_middle_sidebars_xbreak.R run guide

Plots a focal ARGweaver tree interval with flanks, phenotype/allele sidebars, and optional x-axis break handling.

## Usage
```bash
Rscript plot_flanks_middle_sidebars_xbreak.R --smc <FILE.smc> --sites <FILE.sites> --pheno <PHENO.txt> --position <POS> [options]
```

## Options
- `--smc` (required): ARGweaver `.smc` file to scan.
- `--sites` (required): Matching `.sites` file for allele data.
- `--pheno` (required): Phenotype table for proc/site/trait sidebars.
- `--position` (required): Focal base position.
- `--flank`: Number of flanking tree intervals on each side (default: `1`).
- `--trait`: Phenotype column to plot (default: `C13`).
- `--allele_positions`: Comma-separated allele positions (defaults to `--position`).
- `--outdir`: Output directory (default: `trees_out`).
- `--outprefix`: Filename prefix (default: `pos_flanks_mid_sidebars`).
- `--tree_scale`: Tree scaling mode `distance` or `cladogram` (default: `distance`).
- `--ladderize`: Ladderize trees before plotting (default: `true`).
- `--tip_labels_mid`: Show tip labels on the middle tree (default: `false`).
- `--tip_size_mid`: Tip label size on the middle tree (default: `1.1`).
- `--xbreak`: X-axis break mode `none`, `auto`, or `manual` (default: `none`).
- `--xbreak_at`: Break point to use when `--xbreak=manual`.
- `--xbreak_shrink`: Fractional shrink applied at the break (default: `0.03`).
- `--xbreak_trigger`: Distance ratio that triggers an automatic break (default: `6`).
- `--xbreak_q`: Quantile used for auto break detection (default: `0.90`).
- `--width`: Figure width in inches (default: `18`).
- `--height`: Figure height in inches (default: `7`).
- `--dpi`: Output DPI (default: `300`).
- `--format`: Output format `pdf`, `png`, or `both` (default: `png`).
- `--w_left`: Relative width for the left flank tree (default: `4`).
- `--w_mid`: Relative width for the middle block (default: `10`).
- `--w_right`: Relative width for the right flank tree (default: `4`).
- `--w_mid_tree`: Width for the middle tree inside the middle block (default: `5`).
- `--w_proc`: Width for the proc sidebar (default: `0.35`).
- `--w_site`: Width for the site sidebar (default: `0.35`).
- `--w_trait`: Width for the trait sidebar (default: `0.45`).
- `--w_allele`: Width for the allele sidebar (default: `0.9`).

The script stops with a usage message unless the SMC, sites, phenotype files, and focal position are provided.
