# IQ-TREE 2 phylogenetic inference with model selection and support values.
# Option notes:
#   -s: input alignment file (PHYLIP format).
#   -m MFP+ASC: use ModelFinder Plus to pick the best substitution model and
#     apply ascertainment bias correction for alignments that exclude invariant sites.
#   -B 1000: run 1000 ultrafast bootstrap replicates for branch support.
#   --alrt 1000: run 1000 SH-aLRT replicates for additional branch support.
#   --keep-ident: retain identical sequences in the alignment (do not deduplicate).
#   --safe: enable safer but slower search settings to avoid numerical issues.
#   -nt 32: use 32 CPU threads.
#   -pre: output file prefix for all IQ-TREE result files.
iqtree2 \
  -s ALL_rep0_scaffoldOrder.phy \
  -m MFP+ASC \
  -B 1000 \
  --alrt 1000 \
  --keep-ident \
  --safe \
  -nt 32 \
  -pre ALL_taxa_rep0_iqtree_tree2
