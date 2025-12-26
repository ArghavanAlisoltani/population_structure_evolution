# FastTreeMP phylogenetic inference for nucleotide alignments.
# Option notes:
#   -gtr: use the GTR substitution model.
#   -nt: interpret input as nucleotides (DNA) rather than amino acids.
#   ALL_rep0_scaffoldOrder.fa: input alignment in FASTA format.
#   > ALL_rep0_fasttree.nwk: write the Newick tree output to this file.
FastTreeMP -gtr -nt ALL_rep0_scaffoldOrder.fa > ALL_rep0_fasttree.nwk
