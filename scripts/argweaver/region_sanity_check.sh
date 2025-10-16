bcftools view "/scratch/arghavan/LP/subsampling/poly_s100_All_1a1b_renamed.vcf.gz" "scaffold_2:300000000-350000000" | wc -l

bcftools view "/scratch/arghavan/LP/subsampling/split_poly_s100_scaffolds.sorted.vcf.gz" "scaffold_2b:1-30000000" | wc -l

bcftools view "/scratch/arghavan/LP/subsampling/split_poly_s100_scaffolds.sorted.vcf.gz" "scaffold_2b:150000000-180000000" | wc -l

bcftools view "/scratch/arghavan/LP/subsampling/split_poly_s100_scaffolds.sorted.vcf.gz" 

# check 1 tmrca value per 6 million 
bcftools view "/scratch/arghavan/LP/subsampling/poly_s100_All_1a1b_renamed.vcf.gz" "scaffold_5:600000000-606310940" | wc -l
#or

