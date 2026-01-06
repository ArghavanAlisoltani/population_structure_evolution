# Script content guide

This guide lists the analysis scripts in this repository, grouped by folder, with a brief description and the primary input file(s) they expect. Inputs are described at a high level (file type or dataset), since many scripts use hard-coded paths or CLI arguments within the script itself.

## `scripts/` (top-level)
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/ARGweaver_viz.R` | Visualizes TMRCA tracks for a selected scaffold region, highlighting young/old windows and gene overlays. | ARGweaver TMRCA tables and gene annotation tracks. |
| `scripts/Glob_tmrca_v2.R` | Aggregates TMRCA tracks genome-wide, windows them, and plots distributions. | Multiple ARGweaver TMRCA track files. |
| `scripts/One_tmrca_flag_old_young_viz_v3.R` | Windows a single TMRCA file, labels youngest/oldest windows, and plots summaries. | One ARGweaver TMRCA track file. |
| `scripts/Greedy_sample_selection.R` | Selects maximally distinct individuals in PC space using a greedy strategy. | PCA/PC-space coordinates table. |
| `scripts/Summary_100_greedy_samples.R` | Summarizes greedy sample selection results across runs. | Output tables from `Greedy_sample_selection.R`. |
| `scripts/subsampling_Greedy_100.R` | Runs repeated greedy subsampling to evaluate robustness. | PCA coordinates and sampling configuration. |
| `scripts/concat_tmrca.R` | Concatenates TMRCA tracks, derives subsets, and plots violin summaries. | Multiple TMRCA track files (per scaffold). |
| `scripts/manhattan_plot_v2.R` | Builds Manhattan-style plots from FST summaries. | FST window summary table with scaffold coordinates. |
| `scripts/merge_FST_annotations.R` | Merges FST windows with external annotation tables. | FST windows table plus annotation tables (genes/TE/GWAS). |
| `scripts/Visualize_FST_merged_v2.R` | Visualizes merged FST + annotation tables. | Output of `merge_FST_annotations.R`. |
| `scripts/FST_get_summary.R` | Computes descriptive statistics for processed FST windows. | FST window summary tables. |
| `scripts/FST_visualization.R` | Plots FST summaries and significance patterns. | FST window summary tables. |
| `scripts/never_used_prep_genotype_pipeline_v1.sh` | Prototype (unused) pipeline to prepare genotype inputs. | Raw genotype/VCF inputs (paths set in script). |

## `scripts/EDTA`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/EDTA/insertion_time_v7.py` | Summarizes LTR insertion times and TE family comparisons. | EDTA/LTR retriever `.list` outputs. |
| `scripts/EDTA/insertion_time_v8.py` | Newer version of insertion-time summaries. | EDTA/LTR retriever `.list` outputs. |
| `scripts/EDTA/insertion_time_viz_v4.py` | Generates insertion-time plots/figures. | Output tables from insertion-time scripts. |
| `scripts/EDTA/ltr_identity.sh` | Computes LTR identity metrics. | LTR retriever outputs (FASTA/GFF or list files). |
| `scripts/EDTA/concat_passlists.sh` | Concatenates EDTA pass lists for combined processing. | EDTA passlist files. |

## `scripts/BetaScan`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/BetaScan/pre_processing_scripts.sh` | Pre-filters VCFs and converts to PLINK/PHYLIP/FASTA alignments for downstream selection scans. | VCFs, plus PLINK/seqmagick outputs. |

## `scripts/FST`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/FST/FST_calc_v1.sh` | Computes pairwise FST from VCFs (older wrapper). | VCFs per scaffold/region. |
| `scripts/FST/FST_calc_v2.sh` | Computes pairwise FST from VCFs (updated wrapper). | VCFs per scaffold/region. |
| `scripts/FST/fst_windows_v4.py` | Windows FST, benchmarks thresholds, and runs enrichment models. | Per-site FST outputs and window config files. |
| `scripts/FST/merged_fst_proc3-5.R` | Merges FST results across processing batches. | Batch-level FST window tables. |
| `scripts/FST/summarize_sig_windows_by_chrom.py` | Summarizes significant FST windows by chromosome. | FST window table with significance flags. |
| `scripts/FST/te_vs_fst_windows.py` | Compares TE density vs. FST windows (baseline). | FST windows + TE annotation/density tables. |
| `scripts/FST/te_vs_fst_windows.v7.py` | TE vs. FST comparison (v7). | FST windows + TE annotation/density tables. |
| `scripts/FST/te_vs_fst_windows.v8.py` | TE vs. FST comparison (v8). | FST windows + TE annotation/density tables. |
| `scripts/FST/te_vs_fst_windows.v9.py` | TE vs. FST comparison (v9). | FST windows + TE annotation/density tables. |
| `scripts/FST/te_enrichment_and_correlations_v1.py` | Evaluates TE enrichment and correlations with FST signals. | FST windows + TE annotation tables. |
| `scripts/FST/te_qvalue_cor_and_manhattan_v2.py` | Computes TE q-values vs. FST and plots Manhattan-style summaries. | FST windows + TE annotation tables. |
| `scripts/FST/fst_te_ranked_correlations_sigFST_v3.py` | Ranks TE–FST correlations for significant windows. | Significant FST windows + TE annotation data. |
| `scripts/FST/fst_te_enrichment_two_groups.py` | TE enrichment tests split by two sample groups. | FST windows + group labels + TE annotations. |
| `scripts/FST/fst_te_enrichment_two_groups_v1.py` | Earlier two-group TE enrichment variant. | FST windows + group labels + TE annotations. |
| `scripts/FST/tag_gwas_in_sigFST_windows_v1.py` | Tags significant FST windows with GWAS hits (v1). | Significant FST windows + GWAS hit table. |
| `scripts/FST/tag_gwas_in_sigFST_windows_v2.py` | Tags significant FST windows with GWAS hits (v2). | Significant FST windows + GWAS hit table. |
| `scripts/FST/tag_gwas_in_sigFST_windows_v3.py` | Tags significant FST windows with GWAS hits (v3). | Significant FST windows + GWAS hit table. |
| `scripts/FST/tag_gwas_in_sigFST_windows_v4.py` | Tags significant FST windows with GWAS hits (v4). | Significant FST windows + GWAS hit table. |
| `scripts/FST/add_gypsy_copia_enrichment_v1.py` | Adds gypsy/copia TE enrichment columns to summaries. | FST windows + TE class annotation table. |
| `scripts/FST/add_gypsy_copia_enrichment_v2.py` | Updated gypsy/copia enrichment columns. | FST windows + TE class annotation table. |
| `scripts/FST/add_mtag_gwas_to_windows_and_summaries_v2.py` | Integrates MTAG GWAS hits into FST windows/summaries. | FST windows + MTAG GWAS hit table. |
| `scripts/FST/gene_vs_fst_window_v1.py` | Compares gene features vs. FST window signals. | FST windows + gene annotation table. |
| `scripts/FST/annotate_mrna_with_fst_sig_and_go_enrichment.py` | Annotates mRNAs overlapping significant windows and runs GO enrichment. | FST windows + mRNA annotations + GO term mapping. |
| `scripts/FST/te_class_enrichment_panel_v1.py` | Panel plot of TE class enrichment across windows. | TE class enrichment summary tables. |
| `scripts/FST/pos_relabel_TEanno_1ab.sh` | Re-labels TE annotation scaffolds to 1a/1b naming. | TE annotation GFF/GTF. |

### `scripts/FST/manhattan`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/FST/manhattan/manhattan_windows_poi_single_v1.py` | Builds a single-panel Poisson Manhattan plot. | FST window table with Poisson results. |
| `scripts/FST/manhattan/manhattan_windows_poi_threepanel_v2.py` | Builds three-panel Poisson Manhattan plots (v2). | FST window table with Poisson results. |
| `scripts/FST/manhattan/manhattan_windows_poi_threepanel_v3.py` | Builds three-panel Poisson Manhattan plots (v3). | FST window table with Poisson results. |

## `scripts/GWAS`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots.py` | Joins GWAS hits with TMRCA summaries and plots. | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v2.py` | GWAS/TMRCA join + plots (v2). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v3.py` | GWAS/TMRCA join + plots (v3). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v4.py` | GWAS/TMRCA join + plots (v4). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v5.py` | GWAS/TMRCA join + plots (v5). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v5_1.py` | GWAS/TMRCA join + plots (v5.1). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v6.py` | GWAS/TMRCA join + plots (v6). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v7.py` | GWAS/TMRCA join + plots (v7). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v9.py` | GWAS/TMRCA join + plots (v9). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v10.py` | GWAS/TMRCA join + plots (v10). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v12.py` | GWAS/TMRCA join + plots (v12). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v13.py` | GWAS/TMRCA join + plots (v13). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/GWAS_TMRCA_join_and_plots_v14.py` | GWAS/TMRCA join + plots (v14). | GWAS hit table + TMRCA windows. |
| `scripts/GWAS/combine_gwas_hits_pval_maf.py` | Combines GWAS hits, p-values, and MAF filters. | GWAS summary statistics table(s). |

## `scripts/LP_assembly_annotations`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/LP_assembly_annotations/De_duplicate_genes_tissues_v1.R` | Deduplicates gene calls across tissues. | Gene annotation tables per tissue. |
| `scripts/LP_assembly_annotations/modify_manipulate_gene_TE_files.R` | General manipulation of gene/TE annotation files. | Gene annotation and TE annotation tables. |
| `scripts/LP_assembly_annotations/pos_relabel_TEanno_1ab.sh` | Relabels TE coordinates to 1a/1b scaffolds. | TE annotation GFF/GTF. |
| `scripts/LP_assembly_annotations/te_landscape_report.py` | Generates TE landscape summaries from annotations. | TE annotation outputs (GFF/GTF + family tables). |

### `scripts/LP_assembly_annotations/TE`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/LP_assembly_annotations/TE/gff3_te_identity_to_age_v1.py` | Converts TE identity metrics in GFF3 to age estimates. | TE GFF3 with identity fields. |
| `scripts/LP_assembly_annotations/TE/Extract_TE_family_ID_GFF3.sh` | Extracts TE family IDs from GFF3 records. | TE GFF3 annotation file. |
| `scripts/LP_assembly_annotations/TE/extract_parant.sh` | Extracts parent feature references from TE GFF3. | TE GFF3 annotation file. |

## `scripts/python_circus`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/python_circus/circos_multilayer_python_1ab_FIX.py` | Builds multi-layer Circos plots for scaffold 1a/1b data. | Circos data tracks (FST/TE/TMRCA) configured in script. |
| `scripts/python_circus/circos_multilayer_python_1ab_FIX_v2.py` | Updated Circos multi-layer plotting (v2). | Circos data tracks configured in script. |
| `scripts/python_circus/circos_multilayer_python_1ab_FIX_v3.py` | Updated Circos multi-layer plotting (v3). | Circos data tracks configured in script. |

## `scripts/singer`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/singer/run_vcf_prep.sh` | Prepares VCFs for Singer pipeline. | Raw VCFs (paths set in script). |
| `scripts/singer/Singer_single_v5.sh` | Runs Singer analysis on a single sample. | Prepared VCFs + Singer configuration. |

## `scripts/snp_viz`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/snp_viz/SNP_VIZ_v1.R` | Visualizes SNP metrics across scaffolds. | SNP summary tables per scaffold. |

## `scripts/argweaver`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/ARGweaver.sh` | Wrapper to launch ARGweaver runs. | ARGweaver config + VCF/FASTA inputs. |
| `scripts/argweaver/rename_from_map.sh` | Renames scaffolds/IDs using mapping files. | ID mapping table + target files. |
| `scripts/argweaver/run_sed_loop_replace.sh` | Batch renaming using `sed` across files. | List of files + mapping patterns. |
| `scripts/argweaver/split_vcf_with_breakpoints.sh` | Splits VCFs at breakpoint coordinates. | VCF + breakpoint list. |
| `scripts/argweaver/vcf_shifted_position_breakpoints_prep.sh` | Shifts VCF positions based on breakpoints. | VCF + breakpoint list. |
| `scripts/argweaver/region_sanity_check.sh` | Validates region definitions for ARGweaver. | Region definition table. |
| `scripts/argweaver/batch_script_gen.sh` | Generates batch scripts for scaffold runs. | Template parameters in script. |
| `scripts/argweaver/find_breakpoints.R` | Detects scaffold breakpoints. | Coordinate or alignment tables. |
| `scripts/argweaver/Win_150mb_generator_position_shifted.R` | Generates shifted 150 Mb windows. | Scaffold length table. |
| `scripts/argweaver/make_balanced_windows_v1.R` | Creates balanced genomic windows. | Scaffold length table or genome index. |
| `scripts/argweaver/allele_viz.R` | Visualizes allele distributions across windows. | Allele frequency tables from ARGweaver outputs. |
| `scripts/argweaver/R_summary_violin_plot.R` | Violin summaries of TMRCA per scaffold. | TMRCA summary tables. |
| `scripts/argweaver/window_tmrca_v1.py` | Windows ARGweaver TMRCA tracks. | ARGweaver TMRCA track files. |
| `scripts/argweaver/win_compare.py` | Compares windowing strategies for TMRCA tracks. | Windowed TMRCA outputs. |
| `scripts/argweaver/tmrca_windowing_enrichment_v1.py` | Enrichment analyses on windowed TMRCA tracks. | Windowed TMRCA + annotation tables. |
| `scripts/argweaver/length_matched_fisher_tmrca.py` | Length-matched Fisher tests on TMRCA segments. | TMRCA segment tables + annotations. |
| `scripts/argweaver/Length_Matched_Fisher_mRNA_V5.py` | Length-matched Fisher tests for mRNA overlaps. | TMRCA windows + mRNA annotations. |
| `scripts/argweaver/Length-matched_Fisher_TE_v1.py` | Length-matched Fisher tests for TE overlaps. | TMRCA windows + TE annotations. |
| `scripts/argweaver/tmrca_Length_Matched_Fisher_TE_V1.py` | Length-matched Fisher tests for TE overlaps (variant). | TMRCA windows + TE annotations. |
| `scripts/argweaver/staritified_fisher_segment.length_v1.py` | Stratified Fisher tests by segment length. | TMRCA segment tables + annotations. |
| `scripts/argweaver/segment_length–stratified_Fisher_tests_v2.py` | Stratified Fisher tests by segment length (v2). | TMRCA segment tables + annotations. |
| `scripts/argweaver/Annotate_Short_TMRCA_Clusters_V1.py` | Labels short-TMRCA clusters. | TMRCA segment tables. |
| `scripts/argweaver/Annotate_GWAS_with_TMRCA_V2.py` | Annotates GWAS variants with overlapping TMRCA segments. | GWAS hits + TMRCA segment tables. |
| `scripts/argweaver/GWAS_coverage_by_TMRCA_thresholds_V1.py` | GWAS coverage across TMRCA thresholds (v1). | GWAS hits + TMRCA windows. |
| `scripts/argweaver/GWAS_coverage_by_TMRCA_thresholds_V2_fix.py` | GWAS coverage across TMRCA thresholds (v2). | GWAS hits + TMRCA windows. |
| `scripts/argweaver/analyze_tmrca_genes_enrichment.py` | Enrichment of gene sets in TMRCA-defined regions. | TMRCA windows + gene annotation sets. |
| `scripts/argweaver/Check_gene_overlapes.R` | Checks gene overlaps with TMRCA windows. | TMRCA windows + gene annotations. |

### `scripts/argweaver/SNP_wise`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/SNP_wise/SNP_wise.py` | Annotates SNPs with TMRCA/TE/mRNA overlaps. | SNP table + TMRCA windows + TE/mRNA annotations. |
| `scripts/argweaver/SNP_wise/Annotate_SNPs_with_TMRCA_TE_mRNA_V1.py` | SNP-wise TMRCA/TE/mRNA annotation (v1). | SNP table + TMRCA windows + TE/mRNA annotations. |
| `scripts/argweaver/SNP_wise/5_Annotate_GWAS_with_TMRCA_V2.py` | Annotates GWAS SNPs with TMRCA windows. | GWAS SNP table + TMRCA windows. |
| `scripts/argweaver/SNP_wise/GWAS_coverage_by_TMRCA_thresholds_V2_fix.py` | GWAS coverage across TMRCA thresholds (SNP-wise). | GWAS hits + TMRCA windows. |
| `scripts/argweaver/SNP_wise/fisher_snps_tmrca.R` | Fisher tests for SNP overlap with TMRCA windows. | SNP table + TMRCA windows. |

### `scripts/argweaver/Summary_TE_mRNA_TMRCA`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/scaffold_summary_tmrca_te_genes_v2.py` | Summarizes scaffold TMRCA with TE/gene overlaps (v2). | TMRCA windows + TE/gene annotations. |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/scaffold_summary_tmrca_te_genes_v3.py` | Summarizes scaffold TMRCA with TE/gene overlaps (v3). | TMRCA windows + TE/gene annotations. |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/scaffold_summary_tmrca_te_genes_v3a.py` | Summarizes scaffold TMRCA with TE/gene overlaps (v3a). | TMRCA windows + TE/gene annotations. |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/scaffold_summary_tmrca_te_genes_v4.py` | Summarizes scaffold TMRCA with TE/gene overlaps (v4). | TMRCA windows + TE/gene annotations. |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/scaffold_summary_tmrca_te_genes_v5.py` | Summarizes scaffold TMRCA with TE/gene overlaps (v5). | TMRCA windows + TE/gene annotations. |
| `scripts/argweaver/Summary_TE_mRNA_TMRCA/test.sh` | Example shell wrapper for summary scripts. | Script parameters set inside wrapper. |

### `scripts/argweaver/barplots_bins`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/barplots_bins/TMRCA_bins_genes_TEs_SNPs_V2_fast.py` | Builds binned barplots of TMRCA vs. genes/TEs/SNPs. | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/barplots_bins/TMRCA_bins_genes_TEs_SNPs_V3.py` | Binned barplots (v3). | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/barplots_bins/plot_tmrca_bin_bars_v1.py` | Plots TMRCA bin bar charts. | Binned TMRCA summary table. |

### `scripts/argweaver/histograms`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/histograms/Investigate_tmrca_TE_Charles.R` | Exploratory histograms for TMRCA/TE data. | TMRCA windows + TE annotation tables. |
| `scripts/argweaver/histograms/Position_reassign_bb_to_b.R` | Reassigns positions between scaffolds (bb to b). | Coordinate reassignment table. |
| `scripts/argweaver/histograms/TMRCA_bins_genes_TEs_SNPs_V1.py` | Histogram/binning of TMRCA vs. genes/TEs/SNPs (v1). | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/histograms/TMRCA_bins_genes_TEs_SNPs_V2_fast_my_bin_selection.py` | Binned histograms with custom bin selection. | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/histograms/TMRCA_bins_genes_TEs_SNPs_V3_automatic_bin.py` | Binned histograms with automatic binning. | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/histograms/TMRCA_bins_genes_TEs_SNPs_V4.py` | Binned histograms (v4). | TMRCA windows + gene/TE/SNP tables. |
| `scripts/argweaver/histograms/TMRCA_mRNA_TE_Histograms_v1.R` | Histogram plots for mRNA and TE overlaps. | TMRCA windows + mRNA/TE annotations. |
| `scripts/argweaver/histograms/merge_scaffold_1a_1b_to_scaffold_1.py` | Merges scaffold 1a/1b data for histogram analysis. | Scaffold-level TMRCA/annotation tables. |
| `scripts/argweaver/histograms/plot_tmrca_bin_bars_v1.py` | Plots bar charts from binned TMRCA tables (v1). | Binned TMRCA summary table. |
| `scripts/argweaver/histograms/plot_tmrca_bin_bars_v3.py` | Plots bar charts from binned TMRCA tables (v3). | Binned TMRCA summary table. |
| `scripts/argweaver/histograms/plot_tmrca_bin_bars_v4.py` | Plots bar charts from binned TMRCA tables (v4). | Binned TMRCA summary table. |
| `scripts/argweaver/histograms/plot_tmrca_genomewide_ci_v2.py` | Plots genome-wide TMRCA confidence intervals. | Genome-wide TMRCA summary table. |
| `scripts/argweaver/histograms/scaffold_summary_tmrca_te_genes_v5a.py` | Scaffold summary with TE/gene overlaps (v5a). | TMRCA windows + TE/gene annotations. |

### `scripts/argweaver/manipulate_tmrca`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/manipulate_tmrca/1_concat.sh` | Concatenates per-scaffold TMRCA outputs. | Per-scaffold TMRCA files. |
| `scripts/argweaver/manipulate_tmrca/2_check_concat.R` | QC checks on concatenated TMRCA output. | Concatenated TMRCA table. |
| `scripts/argweaver/manipulate_tmrca/Position_reassign_bb_to_b.R` | Reassigns scaffold positions (bb to b). | Coordinate reassignment table. |
| `scripts/argweaver/manipulate_tmrca/merge_scaffold_1a_1b_to_scaffold_1.py` | Merges scaffold 1a/1b TMRCA outputs. | TMRCA outputs for 1a and 1b. |

### `scripts/argweaver/other_scaffolds`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/other_scaffolds/run.sh` | Example ARGweaver run script for other scaffolds. | ARGweaver inputs configured in script. |

### `scripts/argweaver/phylogenetic_trees`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/phylogenetic_trees/auto_tree_contrib_sites_stats_v2.R` | Summarizes contributing sites for tree building. | Site-level tree contribution tables. |
| `scripts/argweaver/phylogenetic_trees/auto_tree_contrib_sites_stats_v2_1.R` | CLI workflow to summarize contributing sites from ARGweaver `.smc`/`.sites` files. | ARGweaver `.smc` and `.sites` files. |
| `scripts/argweaver/phylogenetic_trees/auto_tree_informative_sites_stats_v1.R` | Summarizes informative sites for trees. | Site-level informative-sites tables. |
| `scripts/argweaver/phylogenetic_trees/alignment/GPT_alignment_viz.R` | Extracts selected `.sites` positions, outputs FASTA, and plots base alignment tiles. | ARGweaver `.sites` files with selected positions. |
| `scripts/argweaver/phylogenetic_trees/plot_flanks_mid_auto_files_multiRep_v8_4.R` | Plots trees across flanks/middle with multiple replicates. | Tree output files per region. |
| `scripts/argweaver/phylogenetic_trees/plot_flanks_mid_optional_sidebars_gap_v3.R` | Plots trees with optional sidebars (v3). | Tree output files per region. |
| `scripts/argweaver/phylogenetic_trees/plot_flanks_mid_optional_sidebars_gap_v5_Aria_manual.R` | Plots trees with optional sidebars (v5). | Tree output files per region. |
| `scripts/argweaver/phylogenetic_trees/plot_flanks_middle_sidebars_xbreak.R` | Plots trees with middle sidebars and x-breaks. | Tree output files per region. |
| `scripts/argweaver/phylogenetic_trees/plot_flanks_with_middle_sidebars.R` | Plots trees with middle sidebars. | Tree output files per region. |
| `scripts/argweaver/phylogenetic_trees/sites_tree_heatmap_with_pheno_v4.R` | Heatmap of tree sites with phenotypes. | Site-level tree data + phenotype table. |

### `scripts/argweaver/phylogenetic_trees/tree_sites`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/phylogenetic_trees/tree_sites/build_tree_from_sites_dir.py` | Builds trees from per-site alignment directories. | Directory of per-site alignment files. |
| `scripts/argweaver/phylogenetic_trees/tree_sites/compute_pdist_from_sites_dir.py` | Computes pairwise distances from per-site alignments. | Directory of per-site alignment files. |
| `scripts/argweaver/phylogenetic_trees/tree_sites/concat_sites_rep0_scaffold_order.R` | Concatenates replicate-0 `.sites` files into scaffold-ordered FASTA/PHYLIP alignments. | ARGweaver `.sites` files. |
| `scripts/argweaver/phylogenetic_trees/tree_sites/fasttree.sh` | Runs FastTreeMP on concatenated FASTA alignments. | FASTA alignment from concatenated `.sites`. |
| `scripts/argweaver/phylogenetic_trees/tree_sites/iqtree.sh` | Runs IQ-TREE model selection and support on concatenated PHYLIP alignments. | PHYLIP alignment from concatenated `.sites`. |

### `scripts/argweaver/py_viz_slide`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/py_viz_slide/tmrca_pipeline_v1.py` | TMRCA windowing + visualization pipeline (v1). | ARGweaver TMRCA outputs + annotations. |
| `scripts/argweaver/py_viz_slide/tmrca_pipeline_v2.py` | TMRCA windowing + visualization pipeline (v2). | ARGweaver TMRCA outputs + annotations. |
| `scripts/argweaver/py_viz_slide/tmrca_pipeline_v3.py` | TMRCA windowing + visualization pipeline (v3). | ARGweaver TMRCA outputs + annotations. |
| `scripts/argweaver/py_viz_slide/tmrca_pipeline_v4.py` | TMRCA windowing + visualization pipeline (v4). | ARGweaver TMRCA outputs + annotations. |
| `scripts/argweaver/py_viz_slide/tmrca_pipeline_v5.py` | TMRCA windowing + visualization pipeline (v5). | ARGweaver TMRCA outputs + annotations. |

### `scripts/argweaver/visualization`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/visualization/make_tmrca_arg_deck.py` | Builds a slide deck/figure set for TMRCA results. | TMRCA summary tables + plot outputs. |
| `scripts/argweaver/visualization/plot_tmrca_genomewide_with_ci.R` | Plots genome-wide TMRCA with confidence intervals. | Genome-wide TMRCA summary table. |
| `scripts/argweaver/visualization/plot_tmrca_genomewide_with_ci_v2.R` | Plots genome-wide TMRCA with confidence intervals (v2). | Genome-wide TMRCA summary table. |

### `scripts/argweaver/TMRCA_FST`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/argweaver/TMRCA_FST/FST_TMRCA_association_v3.R` | Merges per-site FST with TMRCA segments and produces regression/summary plots. | Per-site FST table + TMRCA segment table. |
| `scripts/argweaver/TMRCA_FST/FST_staritified_meanTMRCA_v1.R` | Stratifies TMRCA classes and compares high-FST contrasts with scatter/bar summaries. | Merged FST + TMRCA table. |
| `scripts/argweaver/TMRCA_FST/combined_FST_procs_tmrca_v1.R` | Runs joint outlier/quadrant scans and stratified TMRCA–FST analyses. | Merged FST + TMRCA table. |

## `scripts/phylo_tree`
| Script | What it does | Input file(s) |
| --- | --- | --- |
| `scripts/phylo_tree/vcf2aln.sh` | Converts VCFs to PLINK, PHYLIP, and FASTA alignments for phylogenetic inference. | VCF files plus PLINK/seqmagick outputs. |


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
