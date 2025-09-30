module load bcftools/1.8
bcftools view -Oz -o split_poly_s100_scaffolds.vcf.gz \
  split_poly_s100_scaffolds.vcf

# (optional but smart) ensure position-sorted before indexing
bcftools sort -Oz -o split_poly_s100_scaffolds.vcf.gz \
 split_poly_s100_scaffolds.vcf.vcf

# 1) CSI indexes (not TBI)
module load bcftools/1.8
bcftools index -c -f split_poly_s100_scaffolds.vcf.gz
#or
#bcftools index -c split_poly_s100_scaffolds.vcf.gz
