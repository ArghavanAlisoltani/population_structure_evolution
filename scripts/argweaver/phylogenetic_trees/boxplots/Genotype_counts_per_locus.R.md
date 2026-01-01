# Genotype_counts_per_locus.R guide

## Purpose
Summarize per-variant genotype counts (AA/AB/BB/missing) from a VCF and write a TSV report with counts per locus.

## What the script does
- Reads a VCF (optionally gzipped) with `vcfR::read.vcfR`.
- Optionally filters to a single scaffold/contig (`only_chrom`).
- Extracts the GT matrix for all samples.
- Recodes genotypes into `AA`, `AB`, `BB`, or missing (`NA`), treating phased/unphased formats equivalently.
- Computes counts per variant: `n_AA`, `n_AB` (heterozygotes), `n_BB`, `n_missing`, `n_called`.
- Writes a TSV with variant metadata (CHROM/POS/ID/REF/ALT) and counts.

## Key inputs
- `vcf_file`: path to a `.vcf` or `.vcf.gz`.
- `out_tsv`: output TSV path.
- `only_chrom`: set to a specific scaffold ID to filter, or `NULL` for all.

## Outputs
- `genotype_counts_per_locus.tsv` (or `out_tsv` value) with columns:
  - `CHROM`, `POS`, `ID`, `REF`, `ALT`
  - `n_AA`, `n_AB`, `n_BB`, `n_missing`, `n_called`
  - `n_total_samples`

## Run examples
### Example: run as-is after editing input paths
```r
# In R or RStudio
source("scripts/argweaver/phylogenetic_trees/boxplots/Genotype_counts_per_locus.R")
```

### Example: run from command line (after editing inputs)
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/Genotype_counts_per_locus.R
```

### Example: minimal edits for a specific scaffold
```r
# inside the script
vcf_file <- "Imputed_whole_panel_Esteban_Soms_shared.vcf.gz"
only_chrom <- "scaffold_4"
out_tsv <- "scaffold4_genotype_counts.tsv"
```
