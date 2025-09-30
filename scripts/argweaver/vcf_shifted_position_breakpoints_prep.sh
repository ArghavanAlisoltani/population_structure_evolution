module load bcftools/1.8
bcftools view -Oz -o renamed_1ab_samples_1489_merged_sorted_tworules.vcf.gz \
  renamed_1ab_samples_1489_merged_sorted_tworules.vcf

# (optional but smart) ensure position-sorted before indexing
bcftools sort -Oz -o renamed_1ab_samples_1489_merged_sorted_tworules.vcf.gz \
  renamed_1ab_samples_1489_merged_sorted_tworules.vcf

# 1) CSI indexes (not TBI)
module load bcftools/1.8
bcftools index -c -f renamed_1ab_samples_1489_merged_sorted_tworules.vcf.gz
#or
bcftools index -c renamed_1ab_samples_1489_merged_sorted_tworules.vcf.gz
