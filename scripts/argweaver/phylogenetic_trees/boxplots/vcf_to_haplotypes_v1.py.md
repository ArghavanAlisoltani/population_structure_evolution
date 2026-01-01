# vcf_to_haplotypes_v1.py guide

## Purpose
Read a VCF and print two haplotype sequences per sample (FASTA format) for a specified scaffold and position list.

## What the script does
- Parses a VCF line-by-line without external dependencies.
- Collects sample names from the VCF header.
- For each matching scaffold/position, extracts the `GT` field and converts allele indices to bases.
- Builds two haplotype strings per sample and prints them in FASTA format.

## Key inputs (inside the script)
The script is configured by editing the `__main__` block:
- `input_vcf`: VCF path.
- `scaffold`: scaffold/contig string to match.
- `positions`: list of integer positions (1-based).

## Outputs
- FASTA printed to stdout with two records per sample (`<sample>_1` and `<sample>_2`).

## Run examples
### Example: edit inputs and run
```python
# In scripts/argweaver/phylogenetic_trees/boxplots/vcf_to_haplotypes_v1.py
input_vcf = "Imputed_whole_panel_Esteban_Soms_shared.vcf"
scaffold = "scaffold_4"
positions = [983057685, 983057688, 983057694]
```

```bash
python3 scripts/argweaver/phylogenetic_trees/boxplots/vcf_to_haplotypes_v1.py > haplotypes.fasta
```

### Example: redirect to a file
```bash
python3 scripts/argweaver/phylogenetic_trees/boxplots/vcf_to_haplotypes_v1.py \
  > scaffold4_haplotypes.fasta
```
