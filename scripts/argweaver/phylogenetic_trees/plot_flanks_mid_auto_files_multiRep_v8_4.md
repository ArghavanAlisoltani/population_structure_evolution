# plot_flanks_mid_auto_files_multiRep_v8_4.R run guide

This script scans an ARGweaver output directory for SMC/SITES windows that overlap a requested scaffold position, extracts the hit window (plus optional flanks), and generates composite tree + sidebar plots for the requested replicate IDs.

## Example
```bash
Rscript plot_flanks_mid_auto_files_multiRep_v8_4.R \
  --smc_dir=smc_files \
  --sites_dir=sites_for_tree \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --trait=C13 \
  --n_up=0 --n_down=0 \
  --replicates=0,10,20,30,40,50 \
  --outbase=plots_v8 \
  --show_mid_sidebars=true \
  --show_tip_labels_mid=true \
  --sidebars_on_right=true \
  --allele_near_tree=true \
  --w_tree_mid=1.5 \
  --panel_gap=0.25 \
  --gap_before_allele=0 \
  --w_proc=0.45 --w_site=0.45 --w_trait=0.45 --w_allele=0.45 \
  --font_title=11 --font_axis_text=9 --font_axis_title=10 --tip_label_size=2.2 \
  --pdf=false --png=true --width=22 --height=11
```

## Usage
```bash
Rscript plot_flanks_mid_auto_files_multiRep_v8_4.R --smc_dir <DIR> --scaffold <SCAFFOLD_ID> --position <POS> [options]
```

## Options
- `--smc_dir` (required): Directory containing `outargs_*.<rep>.smc` files to scan.
- `--sites_dir` (optional): Directory containing matching `.sites` files.
- `--scaffold` / `--scaffold_id`: Scaffold/contig ID to query.
- `--position`: Focal base position (integer) used to select the hit window.
- `--pheno`: Phenotype table used for sidebar annotations (default: empty string).
- `--trait`: Phenotype column to display (default: `C13`).
- `--n_up`: Number of upstream flank windows (default: `--n_flanks` value).
- `--n_down`: Number of downstream flank windows (default: `--n_flanks` value).
- `--n_flanks`: Symmetric flank count if `--n_up/--n_down` are omitted (default: `2`).
- `--replicates`: Comma/space-separated list of replicate IDs to plot (default: `0`).
- `--outbase`: Base output directory for per-position folders (default: `plots`).
- `--show_mid_sidebars` / `--show_sidebars`: Toggle phenotype/allele sidebars for the center panel (default: `true`).
- `--show_tip_labels_mid`: Show tip labels on the middle tree (default: `true`).
- `--show_proc`: Include processing location sidebar (default: `true`).
- `--show_site`: Include site sidebar (default: `true`).
- `--show_trait`: Include phenotype trait sidebar (default: `true`).
- `--show_allele`: Include allele sidebar (default: `true`).
- `--sidebars_on_right`: Place sidebars to the right of the tree block (default: `true`).
- `--allele_near_tree`: Position allele sidebar adjacent to the tree (default: `true`).
- `--panel_gap`: Gap between panels (default: `0.25`).
- `--gap_before_allele`: Extra gap before allele sidebar (default: `0.01`).
- `--w_proc`: Relative width for the proc sidebar (default: `0.10`).
- `--w_site`: Relative width for the site sidebar (default: `0.10`).
- `--w_trait`: Relative width for the trait sidebar (default: `0.16`).
- `--w_allele`: Relative width for the allele sidebar (default: `0.16`).
- `--w_tree_mid`: Relative width for the middle tree panel (default: `1.0`).
- `--w_mid_panel`: Weight multiplier for the middle panel in the patchwork layout (default: `1.0`).
- `--font_title` / `--title_size`: Plot title font size (default: `11`).
- `--font_axis_text` / `--axis_text_size`: Axis text font size (default: `9`).
- `--font_axis_title` / `--axis_title_size`: Axis title font size (default: `10`).
- `--tip_label_size` / `--tiplab_size`: Tip label font size (default: `2.2`).
- `--linewidth_tree`: Tree line width (default: `0.35`).
- `--axis_title`: X-axis label text (default: `Generations before present (ARGweaver units)`).
- `--axis_label_mid_only`: Show x-axis labels only on the middle panel (default: `true`).
- `--panel_margin_pt`: Margin (points) around each panel (default: `6`).
- `--pdf`: Write PDF output (default: `true`).
- `--png`: Write PNG output (default: `true`).
- `--width`: Figure width in inches (default: `22`).
- `--height`: Figure height in inches (default: `11`).
- `--dpi`: Rasterization DPI (default: `300`).
- `--positions`: Comma-separated allele positions; falls back to `--position` when omitted.

The script errors if `--smc_dir`, `--scaffold`, or `--position` are missing.
