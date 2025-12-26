#1-Pre-filter VCF (biallelic SNPs only)
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz \
  -o filtered_whole_Esteban_imputed.vcf.gz \
  Imputed_whole_panel_Esteban_Soms_shared.vcf


#2- Index it:
bcftools index filtered_whole_Esteban_imputed.vcf.gz 

#3- Convert VCF → PLINK (keeps sample order stable)
plink2 \                      
  --vcf filtered_whole_Esteban_imputed.vcf.gz \
  --allow-extra-chr \
  --snps-only just-acgt \
  --max-alleles 2 \
  --make-bed \
  --out panel

# 4-Optional but STRONGLY recommended filtering two steps
# This avoids pathological trees and speeds everything up.plink2 \
  --bfile panel  --allow-extra-chr \
  --maf 0.05 \
  --geno 0.1 \
  --indep-pairwise 50 5 0.2 \
  --out ld


  plink2 \
  --bfile panel  --allow-extra-chr \
  --extract ld.prune.in \
  --make-bed \
  --out panel_pruned

# 5- PLINK2 → PHYLIP alignment for iqtree2
  plink2 \
  --bfile panel_pruned --allow-extra-chr \
  --export phylip \
  --out alignment

# 6- PHYLIP → FASTA alignment for fasttree
#pip install seqmagick

seqmagick convert alignment.phy alignment.fasta

