# concat_sites_rep0_scaffold_order.R

## Purpose
Concatenate ARGweaver replicate-0 `.sites` files into a single alignment, ordered by scaffold and window coordinates, then export FASTA/PHYLIP outputs for phylogenetic inference.

## Key inputs
- `--sites_dir`: directory containing `*.0.sites` files.
- `--outprefix`: prefix for output alignment files (FASTA/PHYLIP).
- Optional flags: `--recursive`, `--none_last`, `--fasta`, `--phylip`.

## How to run
```bash
Rscript scripts/argweaver/phylogenetic_trees/tree_sites/concat_sites_rep0_scaffold_order.R \
  --sites_dir=sites_for_tree \
  --outprefix=ALL_rep0_concat \
  --recursive=true --none_last=true --fasta=true --phylip=true
```

## Notes
- The script enforces consistent haplotype order across all `.sites` files.
- Outputs include `<outprefix>.fa`, `<outprefix>.phy`, and a TSV listing files in concatenation order.
