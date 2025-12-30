# vcf2aln.sh

## Purpose
Convert an input VCF to filtered, LD-pruned alignments suitable for phylogenetic inference with IQ-TREE or FastTree.

## Key inputs
- `Imputed_whole_panel_Esteban_Soms_shared.vcf`: source VCF.
- `bcftools`, `plink2`, and `seqmagick` installed and available in `PATH`.

## How to run
```bash
bash scripts/phylo_tree/vcf2aln.sh
```

## Notes
- The script creates intermediate PLINK datasets (`panel`, `panel_pruned`) and exports `alignment.phy` and `alignment.fasta`.
- Update filenames or filtering thresholds (`--maf`, `--geno`, `--indep-pairwise`) to match your dataset.
