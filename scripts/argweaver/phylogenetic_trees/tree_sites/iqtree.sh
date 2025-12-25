iqtree2 \
  -s ALL_rep0_scaffoldOrder.phy \
  -m MFP+ASC \
  -B 1000 \
  --alrt 1000 \
  --keep-ident \
  --safe \
  -nt 32 \
  -pre ALL_taxa_rep0_iqtree_tree2
