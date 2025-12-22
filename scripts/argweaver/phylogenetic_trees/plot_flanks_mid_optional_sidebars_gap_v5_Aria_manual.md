# plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R run guide

Manually chooses a center tree and flanking intervals from an ARGweaver `.smc` file, with optional phenotype and allele sidebars.

## Usage
```bash
Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R --smc <FILE.smc> --position <POS> --positions <POS_LIST> [options]
```

## Options
- `--smc` (required): ARGweaver `.smc` file containing tree intervals.
- `--sites`: Matching `.sites` file for allele sidebars (optional but needed when `--show_allele=true`).
- `--pheno`: Phenotype table for sidebar annotations (optional but needed when proc/site/trait sidebars are shown).
- `--position` (required): Center position used to identify the middle tree.
- `--positions`: Comma-separated allele positions for the allele sidebar (defaults to the center position when omitted).
- `--trait`: Phenotype column to plot (default: `C13`).
- `--n_flanks`: Number of flanking trees to include on each side (default: `1`).
- `--outprefix`: Output filename prefix (default: `flanks_mid_tree`).
- `--out_png`: Write PNG output (default: `true`).
- `--out_pdf`: Write PDF output (default: `false`).
- `--show_mid_sidebars`: Show phenotype/allele sidebars next to the middle tree (default: `true`).
- `--show_proc`: Include processing location sidebar (default: `true`).
- `--show_site`: Include site sidebar (default: `true`).
- `--show_trait`: Include trait sidebar (default: `true`).
- `--show_allele`: Include allele sidebar (default: `true`).
- `--show_tip_labels_mid`: Show tip labels on the middle tree (default: `true`).
- `--tip_label_size`: Tip label font size (default: `2.0`).
- `--font_title`: Plot title font size (default: `10`).
- `--font_axis_text`: Axis text font size (default: `8`).
- `--font_axis_title`: Axis title font size (default: `9`).
- `--font_legend_title`: Legend title font size (default: `9`).
- `--font_legend_text`: Legend text font size (default: `8`).
- `--panel_gap`: Gap between tree panels (default: `0.5`).
- `--panel_margin_pt`: Margin (points) around each panel (default: `6`).
- `--sidebar_gap`: Gap between proc/site/trait sidebars (default: `0.05`).
- `--gap_before_allele`: Gap before the allele sidebar (default: `0.05`).
- `--w_tree`: Relative width for the tree panel (default: `4.5`).
- `--w_proc`: Relative width for the proc sidebar (default: `0.6`).
- `--w_site`: Relative width for the site sidebar (default: `0.6`).
- `--w_trait`: Relative width for the trait sidebar (default: `0.8`).
- `--w_allele`: Relative width for the allele sidebar (default: `1.2`).
- `--w_gap_outer`: Override gap width between panels (defaults to `--panel_gap`).
- `--w_gap_sidebar`: Override gap width between sidebars (defaults to `--sidebar_gap`).
- `--w_gap_allele`: Override gap width before the allele panel (defaults to `--gap_before_allele`).

The script prints a usage message and exits when the SMC file, center position, or allele position list is missing.
