ARGweaver phylogenetic tree plotting helpers
===========================================
This folder contains three R scripts that turn ARGweaver .smc/.sites output and phenotype tables into publication-ready phylogenetic plots. Each script can be run directly with `Rscript` and supports command-line arguments to customize the region, flanking windows, and annotation panels.

Scripts and how to run them
--------------------------
Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R \
  --smc=smc_files/outargs_scaffold_4_900000001_1050000000.0.smc \
  --sites=sites_for_tree/outargs_scaffold_4_900000001_1050000000.0.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=983057685 --trait=C13 \
  --n_flanks=0 --outprefix=plots/scaffold_4_pos983057685 \
  --out_png=true \
  --out_pdf=false \
  --show_mid_sidebars=true \
  --show_proc=true \
  --show_site=true \
  --show_trait=true \
  --show_allele=true \
  --show_tip_labels_mid=true \
  --tip_label_size=1.8 \
  --font_title=15 \
  --font_axis_text=10 \
  --font_axis_title=12 \
  --font_legend_title=10 \
  --font_legend_text=12 \
  --panel_gap=0.4 \
  --panel_margin_pt=6 \
  --sidebar_gap=0 \
  --gap_before_allele=0 \
  --w_tree=4.5 \
  --w_proc=0.7 \
  --w_site=0.7 \
  --w_trait=0.7 \
  --w_allele=0.7 \
  --w_gap_outer=0 \
  --w_gap_sidebar=0 \
  --w_gap_allele=0



Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R \
  --smc=outargs_scaffold_8_1050000001_1200000000.10a.smc \
  --sites=outargs_scaffold_8_1050000001_1200000000.10a.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=45802607 --trait=C13 \
  --n_flanks=0 --outprefix=plots/scaffold_8_1095802607_8b_pos45802607.10a \
  --out_png=true \
  --out_pdf=false \
  --show_mid_sidebars=true \
  --show_proc=true \
  --show_site=true \
  --show_trait=true \
  --show_allele=true \
  --show_tip_labels_mid=true \
  --tip_label_size=1.8 \
  --font_title=15 \
  --font_axis_text=10 \
  --font_axis_title=12 \
  --font_legend_title=10 \
  --font_legend_text=12 \
  --panel_gap=0.4 \
  --panel_margin_pt=6 \
  --sidebar_gap=0 \
  --gap_before_allele=0 \
  --w_tree=4.5 \
  --w_proc=0.7 \
  --w_site=0.7 \
  --w_trait=0.7 \
  --w_allele=0.7 \
  --w_gap_outer=0 \
  --w_gap_sidebar=0 \
  --w_gap_allele=0


Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R \
  --smc=outargs_scaffold_1b_1200000001_1350000000.0b.smc \
  --sites=outargs_scaffold_1b_1200000001_1350000000.0b.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=211706312 --trait=lDECLb \
  --n_flanks=0 --outprefix=plots/scaffold_1b_pos211706312_lDECLb \
  --out_png=true \
  --out_pdf=false \
  --show_mid_sidebars=true \
  --show_proc=true \
  --show_site=true \
  --show_trait=true \
  --show_allele=false \
  --show_tip_labels_mid=true \
  --tip_label_size=1.8 \
  --font_title=15 \
  --font_axis_text=10 \
  --font_axis_title=12 \
  --font_legend_title=10 \
  --font_legend_text=12 \
  --panel_gap=0.4 \
  --panel_margin_pt=6 \
  --sidebar_gap=0 \
  --gap_before_allele=0 \
  --w_tree=4.5 \
  --w_proc=0.7 \
  --w_site=0.7 \
  --w_trait=0.7 \
  --w_allele=0.7 \
  --w_gap_outer=0 \
  --w_gap_sidebar=0 \
  --w_gap_allele=0

Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R \
  --smc=outargs_scaffold_1b_1200000001_1350000000.0b.smc \
  --sites=outargs_scaffold_1b_1200000001_1350000000.0b.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=211245637 --trait=lDECLb \
  --n_flanks=0 --outprefix=plots/scaffold_1b_pos211245637_lDECLb \
  --out_png=true \
  --out_pdf=false \
  --show_mid_sidebars=true \
  --show_proc=true \
  --show_site=true \
  --show_trait=true \
  --show_allele=true \
  --show_tip_labels_mid=true \
  --tip_label_size=1.8 \
  --font_title=15 \
  --font_axis_text=10 \
  --font_axis_title=12 \
  --font_legend_title=10 \
  --font_legend_text=12 \
  --panel_gap=0.4 \
  --panel_margin_pt=6 \
  --sidebar_gap=0 \
  --gap_before_allele=0 \
  --w_tree=4.5 \
  --w_proc=0.7 \
  --w_site=0.7 \
  --w_trait=0.7 \
  --w_allele=0.7 \
  --w_gap_outer=0 \
  --w_gap_sidebar=0 \
  --w_gap_allele=0

Rscript plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R \
  --smc=smc_files/outargs_scaffold_31_000000001_013421136.0.smc \
  --sites=sites_for_tree/outargs_scaffold_31_000000001_013421136.0.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=4124403 --trait=lDECLb \
  --n_flanks=0 --outprefix=plots/scaffold_31_pos4124403_lDECLb \
  --out_png=true \
  --out_pdf=false \
  --show_mid_sidebars=true \
  --show_proc=true \
  --show_site=true \
  --show_trait=true \
  --show_allele=true \
  --show_tip_labels_mid=true \
  --tip_label_size=1.8 \
  --font_title=15 \
  --font_axis_text=10 \
  --font_axis_title=12 \
  --font_legend_title=10 \
  --font_legend_text=12 \
  --panel_gap=0.4 \
  --panel_margin_pt=6 \
  --sidebar_gap=0 \
  --gap_before_allele=0 \
  --w_tree=4.5 \
  --w_proc=0.7 \
  --w_site=0.7 \
  --w_trait=0.7 \
  --w_allele=0.7 \
  --w_gap_outer=0 \
  --w_gap_sidebar=0 \
  --w_gap_allele=0

#####################
# v8 automatic file finder
######################
Dec_16_2025 % Rscript plot_flanks_mid_auto_files_multiRep_v8_4.R \
  --smc_dir=smc_files \
  --sites_dir=sites_for_tree \
  --scaffold=scaffold_4 \
  --position=983057685 \
  --pheno=PHENO_Charles_6_2025.txt \
  --trait=C13 \
  --n_up=0 --n_down=0 \
  --replicates=0,10,20,30,40,50 \
  --outbase=plots_v8_4 \
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



- **plot_flanks_middle_sidebars_xbreak.R** – Plots the tree covering a focal SNP with configurable upstream/downstream flanks, optional manual crossover lines, phenotype barplots, and allele sidebars. Outputs trees plus per-figure metadata into `--outdir` with prefix `--outprefix`.
  - Example:
    ```
    Rscript plot_flanks_middle_sidebars_xbreak.R \
      --smc=outargs_scaffold_4_900000001_1050000000.0.smc \
      --sites=sites_for_tree/outargs_scaffold_4_900000001_1050000000.50.sites \
      --pheno=PHENO_Charles_6_2025.txt \
      --position=983057685 \
      --flank=2 \
      --trait=C13 \
      --allele_positions=983057685 \
      --xbreak=manual --xbreak_at=1800 --xbreak_shrink=0.002 \
      --outdir=scaffold4_plots_xbreak \
      --outprefix=scaffold4_pos983057685_1800
    ```

- **plot_flanks_with_middle_sidebars.R** – Similar to the script above but without crossover lines. Generates a tree for the focal position, adds allele sidebars, and can control units/labels in the x-axis.
  - Example:
    ```
    Rscript plot_flanks_with_middle_sidebars.R \
      --smc=smc_files/outargs_scaffold_6_000000001_150000000.30.smc \
      --sites=sites_for_tree/outargs_scaffold_6_000000001_150000000.30.sites \
      --pheno=PHENO_Charles_6_2025.txt \
      --position=72516188 \
      --flank=2 \
      --trait=lcamphenedwb \
      --allele_positions=72516188 \
      --outdir=scaffold6_plots_sidebar \
      --outprefix=scaffold6_flanks2_30_pos983057685_midSidebars \
      --format=png \
      --x_unit=bp \
      --bp_labels=true \
      --tip_labels_mid=false
    ```

Rscript plot_flanks_mid_optional_sidebars_gap_v6.R \
  --smc=smc_files/outargs_scaffold_4_900000001_1050000000.0.smc \
  --sites=sites_for_tree/outargs_scaffold_4_900000001_1050000000.0.sites \
  --pheno=PHENO_Charles_6_2025.txt \
  --position=983057685 \
  --n_up=1 --n_down=1 \
  --trait=C13 \
  --outprefix=v6_scaffold4_pos983057685_0 \
  --show_sidebars=true \
  --show_tiplab_mid=true \
  --panel_gap=0.25 \
  --w_gap_before_allele=0 \
  --w_proc=0.4 --w_site=0.4 --w_trait=0.4 --w_allele=0.4 \
  --title_size=11 --axis_text_size=9 --axis_title_size=10 --tiplab_size=2.2 \
  --pdf=true --png=true --width=22 --height=12


- **sites_tree_heatmap_with_pheno_v4.R** – Builds a Neighbor-Joining tree from a `.sites` file, overlays allele calls for selected positions as a heatmap, and joins phenotype values. Writes both the tree (`.nwk`) and the plotted figure with the chosen prefix.
  - Example:
    ```
    Rscript sites_tree_heatmap_with_pheno_v4.R \
      --sites=sites_for_tree/outargs_scaffold_4_900000001_1050000000.50.sites \
      --pheno=PHENO_Charles_6_2025.txt \
      --positions=983057685,983057700 \
      --trait=C13 \
      --outprefix=scaffold4_pos983057685_heatmap \
      --tree_scale=distance \
      --hide_labels=true
    ```

How the phylogenetic trees are built
------------------------------------
1. **Convert .sites into per-sample sequences (alignment)** – A `.sites` file is site-major (one row per genomic position with an allele string across samples). The script transposes this into sample-major sequences and stores them as an `ape::DNAbin` alignment.
2. **Compute a JC69 pairwise distance matrix** – `dist.dna(aln, model="JC69", pairwise.deletion=TRUE)` estimates distances between every pair of sequences, ignoring sites with missing/ambiguous data in either sample for that pair.
3. **Build a Neighbor-Joining (NJ) tree** – `nj()`/`bionj()` constructs a distance-based tree that approximately minimizes total branch length.

Common caveats
--------------
- Few sites can make the NJ topology unstable because many distances tie.
- Non-ACGT characters (e.g., 0/1) may be treated as ambiguous; swap in a Hamming-distance version if your data use binary alleles.
- Nearly identical sequences produce very short branches that can visually collapse.
