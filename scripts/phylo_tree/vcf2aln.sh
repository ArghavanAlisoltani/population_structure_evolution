#Pre-filter VCF (biallelic SNPs only)
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz \
  -o filtered_whole_Esteban_imputed.vcf.gz \
  Imputed_whole_panel_Esteban_Soms_shared.vcf


#Index it:
