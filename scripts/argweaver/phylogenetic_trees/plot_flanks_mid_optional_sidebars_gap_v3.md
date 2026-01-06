# plot_flanks_mid_optional_sidebars_gap_v3.R run guide

Plots an ARGweaver SMC interval that contains a focal position, plus optional flanking trees and phenotype/allele sidebars.

## Usage
```bash
Rscript plot_flanks_mid_optional_sidebars_gap_v3.R --smc <FILE.smc> --position <POS> [options]
```

## Options
- `--smc` (required): ARGweaver `.smc` file to scan.
- `--sites`: Matching `.sites` file for allele sidebars.
- `--pheno`: Phenotype table for sidebar annotations.
- `--position` (required): Focal base position.
- `--flank`: Number of flanking tree intervals on each side (default: `1`).
- `--trait`: Phenotype column to plot (default: `C13`).
- `--allele_positions`: Comma-separated allele positions (default: uses `--position`).
- `--outdir`: Output directory (default: `trees_out`).
- `--outprefix`: Filename prefix (default: `pos_flanks_mid_sidebars`).
- `--tree_scale`: Tree scaling mode `distance` or `cladogram` (default: `distance`).
- `--ladderize`: Ladderize trees before plotting (default: `true`).
- `--xlab_mode`: X-axis labeling mode `distance` or `time_before_present` (default: `time_before_present`).
- `--tip_labels_hit`: Show tip labels on the hit tree (default: `true`).
- `--tip_size_hit`: Tip label size for the hit tree (default: `0.9`).
- `--tip_labels_flanks`: Show tip labels on flank trees (default: `false`).
- `--tip_size_flanks`: Tip label size for flank trees (default: `0.6`).
- `--panel_gap`: Gap (points) between panels (default: `25`).
- `--show_sidebars`: Master toggle for all sidebars (default: `true`).
- `--show_proc`: Include processing location sidebar (default: `true`).
- `--show_site`: Include site sidebar (default: `true`).
- `--show_trait`: Include trait sidebar (default: `true`).
- `--show_allele`: Include allele sidebar (default: `true`).
- `--width`: Figure width in inches (default: `22`).
- `--height`: Figure height in inches (default: `7`).
- `--dpi`: Output DPI (default: `300`).
- `--format`: Output format `png`, `pdf`, or `both` (default: `png`).
- `--w_other`: Relative width for non-hit flanks (default: `4`).
- `--w_midblock`: Width for the middle block containing the hit tree (default: `10`).
- `--w_hit_tree`: Width for the hit tree portion of the middle block (default: `5`).
- `--w_proc`: Width for the proc sidebar (default: `0.35`).
- `--w_site`: Width for the site sidebar (default: `0.35`).
- `--w_trait`: Width for the trait sidebar (default: `0.45`).
- `--w_allele`: Width for the allele sidebar (default: `0.9`).
- `--mid_spacer`: Optional spacer width inside the middle block (default: `0`).

Sidebars require `--pheno` and/or `--sites` when the corresponding toggles are enabled; the script aborts if the required files are missing.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
