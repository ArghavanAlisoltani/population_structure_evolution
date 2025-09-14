# 1) ARGweaver MCMC on the region
arg-sample \
  --vcf "path/to/s100_All_1a1b_renamed.filtered.vcf.gz" \
  --region "scaffold_1a:750000001-900000000" \
  -N 10000 --randseed 1753216431 -r 1.6e-8 -m 1.8e-8 \
  --ntimes 20 --maxtime 2e5 -c 20 -n 50 \
  -o "/outargs_scaffold_1a_750000001_900000000"

# 2) Posterior-mean TMRCA 
arg-extract-tmrca "/outargs_scaffold_1a_750000001_900000000.%d.smc.gz" \
> "/scaffold_1a_750000001_900000000/outargs_scaffold_1a_750000001_900000000.tmrca.txt"
