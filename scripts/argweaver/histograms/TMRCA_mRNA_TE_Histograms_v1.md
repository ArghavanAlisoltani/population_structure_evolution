# TMRCA_mRNA_TE_Histograms_v1.R

## Purpose
Builds histogram-style summaries of SNP TMRCA bins and overlays mRNA/TE overlap proportions, producing TSV and PNG outputs.

## What the script does
- Reads SNP annotations that already include `mean_tmrca_at_pos`, `mrna_id`, and TE overlap indicators.
- Bins SNPs into predefined TMRCA intervals, counting SNPs per bin and the number/proportion overlapping mRNA and TE features.
- Saves the bin summary to TSV and generates three plots: SNP counts per bin, SNP counts with mRNA proportions, and SNP counts with TE proportions.

## How to run
1. Adjust the `setwd()` and input/output filenames at the top of the script to match your data.
2. Execute with R:
   ```bash
   Rscript TMRCA_mRNA_TE_Histograms_v1.R
   ```
   The script writes one TSV and three PNG figures in the working directory.
