# plot_flanks_with_middle_sidebars.R run guide

Plots a focal ARGweaver tree interval with flanking trees and middle-panel phenotype/allele sidebars.

## Usage
```bash
Rscript plot_flanks_with_middle_sidebars.R --smc <FILE.smc> --sites <FILE.sites> --pheno <PHENO.txt> --position <POS> [options]
```

## Options
- `--smc` (required): ARGweaver `.smc` file to scan.
- `--sites` (required): Matching `.sites` file for allele data.
- `--pheno` (required): Phenotype table for sidebar annotations.
- `--position` (required): Focal base position.
- `--flank`: Number of flanking tree intervals on each side (default: `1`).
- `--trait`: Phenotype column to plot (default: `C13`).
- `--allele_positions`: Comma-separated allele positions (defaults to `--position`).
- `--outdir`: Output directory (default: `trees_out`).
- `--outprefix`: Filename prefix (default: `pos_flanks_mid_sidebars`).
- `--tree_scale`: Tree scaling mode `distance` or `cladogram` (default: `distance`).
- `--ladderize`: Ladderize trees before plotting (default: `true`).
- `--tip_labels_mid`: Show tip labels on the middle tree (default: `false`).
- `--tip_size_mid`: Tip label size on the middle tree (default: `1.3`).
- `--width`: Figure width in inches (default: `18`).
- `--height`: Figure height in inches (default: `8`).
- `--dpi`: Output DPI (default: `300`).
- `--format`: Output format `pdf`, `png`, or `both` (default: `pdf`).
- `--w_left`: Relative width for the left flank tree (default: `4`).
- `--w_mid`: Relative width for the middle block (default: `10`).
- `--w_right`: Relative width for the right flank tree (default: `4`).
- `--w_mid_tree`: Width for the middle tree inside the middle block (default: `5`).
- `--w_proc`: Width for the proc sidebar (default: `0.35`).
- `--w_site`: Width for the site sidebar (default: `0.35`).
- `--w_trait`: Width for the trait sidebar (default: `0.45`).
- `--w_allele`: Width for the allele sidebar (default: `0.9`).
- `--x_unit`: Genomic unit for x-axis labels `bp` or `Mb` (default: `bp`).
- `--bp_labels`: Show base-pair labels on the middle x-axis (default: `true`).

The script stops with a usage message unless the SMC, sites, phenotype files, and focal position are provided.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
