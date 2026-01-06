# fasttree.sh

## Purpose
Run FastTreeMP on a nucleotide alignment to infer a phylogenetic tree in Newick format.

## Key inputs
- `ALL_rep0_scaffoldOrder.fa`: input FASTA alignment.
- `FastTreeMP`: must be installed and in `PATH`.

## How to run
```bash
bash scripts/argweaver/phylogenetic_trees/tree_sites/fasttree.sh
```

## Notes
- The script outputs `ALL_rep0_fasttree.nwk` in the current working directory.
- Update the input/output filenames or options if you use different alignment names.


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
