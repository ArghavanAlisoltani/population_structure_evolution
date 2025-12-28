#
bcftools view -m2 -M2 -v snps -r scaffold_4 Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
| bcftools +fill-tags -Oz -o scaffold_4.biallelic.tags.vcf.gz -- -t AN,AC,AF
bcftools index -f scaffold_4.biallelic.tags.vcf.gz
