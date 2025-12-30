# iqtree.sh

## Purpose
Run IQ-TREE 2 with model selection and support metrics on a nucleotide alignment in PHYLIP format.

## Key inputs
- `ALL_rep0_scaffoldOrder.phy`: input PHYLIP alignment.
- `iqtree2`: must be installed and in `PATH`.

## How to run
```bash
bash scripts/argweaver/phylogenetic_trees/tree_sites/iqtree.sh
```

## Notes
- Outputs are written with the `ALL_taxa_rep0_iqtree_tree2` prefix in the current directory.
- Adjust thread count (`-nt`) and model settings as needed for your environment.
