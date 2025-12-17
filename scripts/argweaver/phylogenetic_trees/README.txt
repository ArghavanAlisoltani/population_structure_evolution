How the phylogenetic tree is constructed in this script
Step 1: Convert .sites into per-sample sequences (alignment)

Your .sites file is “site-major”:

Each row = one genomic position

The allele string contains one character per sample

The script transposes that into “sample-major” sequences:

For each sample, it collects its allele at every site in the file

That becomes a DNA sequence string like "ACGTACG..."

This is stored as an ape::DNAbin alignment object (aln).

Step 2: Compute a pairwise distance matrix (JC69)

dist.dna(aln, model="JC69") computes distances between every pair of sequences.

For two samples, it first estimates the fraction of differing sites (p)

Under Jukes–Cantor (JC69), it converts that to an evolutionary distance:

𝑑=−3/4ln(

Key detail: pairwise.deletion=TRUE means for each pair, sites with missing/ambiguous data in either sample are ignored.

Output d is a symmetric distance object: sample×sample.

Step 3: Build a Neighbor-Joining (NJ) tree from the distances

nj(d) builds a tree that approximately minimizes total branch length given those distances.

Conceptually:

It starts with all samples as separate tips (“star”)

Iteratively finds the “best” pair of taxa/clusters to join based on the NJ criterion

Creates an internal node, assigns branch lengths

Repeats until only two clusters remain

So this tree is a distance-based tree, not a likelihood/Bayesian coalescent tree.

Important caveats (why a tree can look “weird”)

Even with correct code, trees can look odd if:

The file is small (few sites): distances have many ties → NJ topology unstable

Alleles are not A/C/G/T (e.g., 0/1): as.DNAbin + dist.dna may treat them as ambiguous/missing

Branch lengths near zero (very similar sequences): tree collapses visually (comb-like)

If your .sites uses 0/1, tell me and I’ll give you a version that constructs the tree using Hamming distance on characters (treating 0/1 as valid states) instead of dist.dna(JC69).
