# Scripts directory overview

This repository contains many analysis utilities for population structure and evolutionary genomics.  The notes below summarize what each script or helper in `scripts/` and its subfolders is meant to do so you can quickly locate the right tool.

## Top-level R and shell scripts
- **ARGweaver_viz.R** – Visualize TMRCA tracks for a selected scaffold region, flag young/old windows, and overlay gene annotations.
- **Glob_tmrca_v2.R** – Read multiple TMRCA tracks, window them genome-wide, flag extremes, and plot aggregated distributions.
- **One_tmrca_flag_old_young_viz_v3.R** – Window a single TMRCA file, label youngest/oldest windows, and produce summary plots.
- **Greedy_sample_selection.R** – Select maximally distinct individuals in PC space using a greedy max–min strategy with adjustable separation thresholds.
- **Summary_100_greedy_samples.R** – Summarize greedy sample selections, including per-run overlap and diversity metrics.
- **subsampling_Greedy_100.R** – Run repeated greedy subsampling to evaluate robustness of the 100-sample selection.
- **concat_tmrca.R** – Concatenate TMRCA tracks, derive high-TMRCA subsets, and generate violin plots per scaffold.
- **manhattan_plot_v2.R** – Create Manhattan-style plots from FST tables with cumulative scaffold coordinates.
- **merge_FST_annotations.R** – Combine FST window summaries with external annotations for downstream visualization.
- **Visualize_FST_merged_v2.R** – Plot merged FST/annotation tables, highlighting significant windows or features.
- **FST_get_summary.R** – Compute descriptive summaries from processed FST windows.
- **never_used_prep_genotype_pipeline_v1.sh** – Unused prototype shell workflow for preparing genotype inputs.

## FST analyses (`scripts/FST`)
These scripts focus on FST calculation, windowing, enrichment, and annotation.
- **FST_calc_v1.sh / FST_calc_v2.sh** – Shell wrappers to compute pairwise FST from VCFs across scaffolds.
- **fst_windows_v4.py** – Benchmark multiple window sizes/overlaps, compute counts above FST thresholds, run enrichment models (binomial/Poisson/etc.), and generate PowerPoint reports.
- **merged_fst_proc3-5.R** – Merge FST results across processing batches 3–5 for consolidated analyses.
- **summarize_sig_windows_by_chrom.py** – Aggregate significant FST windows per chromosome/scaffold.
- **te_vs_fst_windows.py / v7.py / v8.py / v9.py** – Compare TE densities against FST windows across successive versions of the analysis.
- **te_enrichment_and_correlations_v1.py** – Evaluate TE enrichment and correlations with FST signals.
- **te_qvalue_cor_and_manhattan_v2.py** – Compute TE q-values versus FST windows and plot Manhattan-style summaries.
- **fst_te_ranked_correlations_sigFST_v3.py** – Rank TE–FST correlations focusing on significant windows.
- **fst_te_enrichment_two_groups.py / v1.py** – TE enrichment tests split by two sample groupings.
- **tag_gwas_in_sigFST_windows_v1-4.py** – Tag significant FST windows with GWAS hits (several iterations).
- **add_gypsy_copia_enrichment_v1.py / v2.py** – Add gypsy/copia TE enrichment columns to FST summaries.
- **add_mtag_gwas_to_windows_and_summaries_v2.py** – Integrate MTAG GWAS hits into FST window tables and summaries.
- **gene_vs_fst_window_v1.py** – Compare gene features against FST window signals.
- **annotate_mrna_with_fst_sig_and_go_enrichment.py** – Mark mRNAs overlapping significant FST windows and run GO enrichment on significant vs. non-significant sets.
- **te_class_enrichment_panel_v1.py** – Panel plot summarizing TE class enrichment across windows.
- **pos_relabel_TEanno_1ab.sh** – Reassign scaffold labels in TE annotations to match 1a/1b naming.

## ARGweaver analyses (`scripts/argweaver`)
ARGweaver utilities for TMRCA processing, enrichment, and visualization.
- **ARGweaver.sh** – Wrapper to launch ARGweaver runs with project-specific settings.
- **rename_from_map.sh / run_sed_loop_replace.sh** – Helpers to rename scaffolds/IDs using mapping files.
- **split_vcf_with_breakpoints.sh / vcf_shifted_position_breakpoints_prep.sh** – Split or shift VCFs at predefined breakpoints for region-specific analyses.
- **region_sanity_check.sh** – Validate region definitions before running ARGweaver.
- **win_compare.py** – Compare different windowing strategies on TMRCA tracks.
- **window_tmrca_v1.py** – Window TMRCA outputs and export summaries.
- **tmrca_windowing_enrichment_v1.py** – Perform enrichment analyses on windowed TMRCA tracks.
- **length_matched_fisher_tmrca.py** – Length-matched Fisher tests on TMRCA segments.
- **Length_Matched_Fisher_mRNA_V5.py** – Apply length-matched Fisher testing to mRNA overlaps.
- **Length-matched_Fisher_TE_v1.py / tmrca_Length_Matched_Fisher_TE_V1.py** – Length-matched Fisher tests focused on TE overlaps with TMRCA windows.
- **segment_length–stratified_Fisher_tests_v2.py / staritified_fisher_segment.length_v1.py** – Segment-length–stratified Fisher testing variants.
- **Annotate_Short_TMRCA_Clusters_V1.py** – Label short TMRCA clusters for downstream reporting.
- **GWAS_coverage_by_TMRCA_thresholds_V1.py / V2_fix.py** – Quantify GWAS hit coverage across TMRCA threshold bins.
- **Annotate_GWAS_with_TMRCA_V2.py** – Annotate GWAS variants with overlapping TMRCA segments.
- **analyze_tmrca_genes_enrichment.py** – Assess enrichment of gene sets within TMRCA-defined regions.
- **allele_viz.R** – Visualize allele frequency distributions across ARGweaver windows.
- **R_summary_violin_plot.R** – Violin-plot summaries of TMRCA by scaffold.
- **make_balanced_windows_v1.R** – Generate balanced genomic windows for fair enrichment testing.
- **find_breakpoints.R / breakpoints.txt / changed_region.txt** – Identify or document scaffold breakpoints for splitting analyses.
- **Win_150mb_generator_position_shifted.R / scaffold_windows_150Mb.tsv** – Generate shifted 150 Mb windows for scaffolds.
- **batch_script_gen.sh** – Template generator for batch processing scripts.
- **ARGweaver/other_scaffolds/run.sh** – Example run script for additional scaffolds.

### ARGweaver subfolders
- **manipulate_tmrca/** – Tools for concatenating and QC-ing TMRCA outputs (`1_concat.sh`, `2_check_concat.R`, `Position_reassign_bb_to_b.R`).
- **SNP_wise/** – SNP-level TMRCA annotation and enrichment (`Annotate_SNPs_with_TMRCA_TE_mRNA_V1.py`, `SNP_wise.py`, `fisher_snps_tmrca.R`).
- **Summary_TE_mRNA_TMRCA/** – Scaffold-level summaries combining TMRCA with TE and gene annotations (several `scaffold_summary_tmrca_te_genes_v*.py` versions plus `test.sh`).
- **py_viz_slide/** – TMRCA windowing + visualization pipelines with successive versions (`tmrca_pipeline_v1.py`–`v5.py` and `readme.txt`).
- **other_scaffolds/** – Supplemental run scripts for non-primary scaffolds.

## EDTA insertion-time analyses (`scripts/EDTA`)
- **insertion_time_v7.py / insertion_time_v8.py** – Summaries and visualizations of LTR retriever `.list` files, including rescaled insertion times and TE family comparisons.
- **insertion_time_viz_v4.py** – Plot-focused insertion time visualization pipeline.
- **ltr_identity.sh** – Shell helper to derive LTR identity metrics.
- **concat_passlists.sh** – Concatenate EDTA pass lists for combined processing.
- **list.txt** – Example list of EDTA output files.

## LP assembly annotation helpers (`scripts/LP_assembly_annotations`)
- **De_duplicate_genes_tissues_v1.R** – Collapse redundant gene calls across tissues.
- **modify_manipulate_gene_TE_files.R** – General wrangling of gene/TE annotation tables.
- **pos_relabel_TEanno_1ab.sh** – Relabel TE annotation coordinates to the 1a/1b scaffold convention.
- **TE/** – TE-specific helpers:
  - **gff3_te_identity_to_age_v1.py** – Convert TE identity metrics in GFF3 to estimated ages.
  - **Extract_TE_family_ID_GFF3.sh** – Pull TE family IDs from GFF3 annotations.
  - **extract_parant.sh** – Extract parent feature references from TE GFF3 files.
  - **te_landscape_report.py** – Build TE landscape summaries from annotation outputs.

## Singer workflows (`scripts/singer`)
- **run_vcf_prep.sh** – Prepare VCFs for Singer analysis pipeline.
- **Singer_single_v5.sh** – Single-sample Singer execution wrapper.

## SNP visualization (`scripts/snp_viz`)
- **SNP_VIZ_v1.R** – Visualization of SNP metrics across scaffolds.

## Argweaver helper folders
- **argweaver/** – (see above) primary ARGweaver analysis suite.
- **EDTA/** – (see above) TE insertion-time analyses.
- **FST/** – (see above) FST enrichment and annotation scripts.
- **LP_assembly_annotations/** – (see above) assembly/annotation utilities.
- **singer/** – (see above) Singer pipeline helpers.
- **snp_viz/** – (see above) SNP visualization.



## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
