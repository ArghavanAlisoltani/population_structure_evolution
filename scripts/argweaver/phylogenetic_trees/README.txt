ARGweaver phylogenetic tree plotting helpers
===========================================
This folder contains three R scripts that turn ARGweaver .smc/.sites output and phenotype tables into publication-ready phylogenetic plots. Each script can be run directly with `Rscript` and supports command-line arguments to customize the region, flanking windows, and annotation panels.

Scripts and how to run them
---------------------------
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
