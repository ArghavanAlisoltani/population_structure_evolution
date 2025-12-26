#Pre-filter VCF (biallelic SNPs only)
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz \
  -o filtered_whole_Esteban_imputed.vcf.gz \
  Imputed_whole_panel_Esteban_Soms_shared.vcf


#Index it:
bcftools index filtered_whole_Esteban_imputed.vcf.gz 

#Convert VCF → PLINK (keeps sample order stable)
plink2 \                      
  --vcf filtered_whole_Esteban_imputed.vcf.gz \
  --allow-extra-chr \
  --snps-only just-acgt \
  --max-alleles 2 \
  --make-bed \
  --out panel

#Convert VCF → PLINK (keeps sample order stable)
plink2 \
  --bfile panel  --allow-extra-chr \
  --maf 0.05 \
  --geno 0.1 \
  --indep-pairwise 50 5 0.2 \
  --out ld
