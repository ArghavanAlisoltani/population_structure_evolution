# vcf_hap_box_plots_v1.R guide

## Purpose
Build haplotype groups from selected SNP positions in a VCF, merge with phenotype data, and generate boxplots/jitter plots of trait values by haplotype group.

## What the script does
- Ensures the VCF is bgzip + tabix indexed (or creates a temporary indexed copy).
- Extracts genotypes at the requested positions using `bcftools query`.
- Converts per-sample genotypes into haplotype strings across the selected positions.
- Collapses identical haplotype strings into haplotype groups and counts unique individuals per group.
- Merges haplotype membership with phenotype traits.
- Generates a boxplot + jitter plot per haplotype group and saves plot data.

## Requirements
- R packages: `data.table`, `ggplot2`, `stringr`, `optparse`.
- Command-line tools in PATH: `bcftools` (required), `bgzip` and `tabix` (if VCF is not already `.vcf.gz` with `.tbi`).

## Key inputs (CLI)
- `--vcf`: input VCF (`.vcf` or `.vcf.gz`).
- `--pheno`: phenotype table (TSV/CSV).
- `--pheno_sep`: delimiter for the phenotype file (default: tab).
- `--id_col`: ID column in the phenotype file (default: `codg`).
- `--scaffold`: scaffold/contig name (e.g., `scaffold_4`).
- `--positions`: comma-separated list of positions (1-based).
- `--traits`: comma-separated traits to plot (first trait is used for y-axis).
- `--outdir`: output directory base.

## Outputs
Within a scaffold/position-specific directory (e.g., `vcf_hap_box_out/scaffold_4_983057685_983057724`):
- `vcf_selected_positions_raw.tsv`: extracted GT table.
- `hap_code_to_sequence.tsv`: haplotype code mapping.
- `haplotype_counts_unique_individuals.tsv`: unique individual counts per haplotype.
- `sample_to_haplotype_groups.tsv`: sample-level haplotype membership.
- `plot_data_<trait>_hapgroups.tsv`: merged plot data.
- `box_jitter_<trait>_by_HAPLOTYPES_vcf.png` (and/or `.pdf`).

## Run examples
### Basic example
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/vcf_hap_box_plots_v1.R \
  --vcf Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
  --pheno phenotypes.tsv \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057694,983057707,983057714,983057718,983057724 \
  --traits C13
```

### Example with custom output and multiple traits
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/vcf_hap_box_plots_v1.R \
  --vcf Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
  --pheno phenotypes.tsv \
  --pheno_sep , \
  --id_col SampleID \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057694 \
  --traits C13,Height,Diameter \
  --outdir hap_boxplots \
  --png TRUE --pdf TRUE
```


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
