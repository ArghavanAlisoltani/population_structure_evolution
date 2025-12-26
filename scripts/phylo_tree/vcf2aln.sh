#Pre-filter VCF (biallelic SNPs only)
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz \
  -o filtered.vcf.gz \
  Imputed_whole_panel_Esteban_Soms_shared.vcf

