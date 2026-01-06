# vcf_haplotype_groups_v1.py guide

## Purpose
Extract phased haplotype strings from selected SNP positions in a VCF, group identical haplotypes, and write summary tables plus FASTA outputs.

## What the script does
- Parses sites from either `--sites_file`, `--scaffold` + `--positions`, or defaults.
- Uses `bcftools query` to pull genotypes for the requested sites.
- Converts GT fields into two per-sample haplotype strings across the site list.
- Writes per-sample haplotypes in FASTA format.
- Groups identical haplotypes, assigns haplotype IDs, and writes summary tables.
- Reports missing sites and outputs paths to generated files.

## Requirements
- Python 3.
- `bcftools` available in PATH.
- VCF should be bgzip + tabix indexed (`.vcf.gz` with `.tbi`/`.csi`).

## Key inputs (CLI)
- `--vcf`: input VCF (`.vcf.gz`).
- `--outdir`: output directory (default: `hap_out`).
- `--scaffold`: scaffold/contig (required with `--positions`).
- `--positions`: comma-separated list of 1-based positions.
- `--sites_file`: file with `CHROM POS` per line (overrides `--positions`).
- `--missing_char`: placeholder for missing/invalid alleles (default `N`).
- `--drop_any_missing`: drop any haplotype containing missing chars from grouping.

## Outputs
In `--outdir`:
- `sites.regions.txt`: region list for bcftools query.
- `loci_used.tsv`: requested loci with status.
- `haplotypes_per_sample.fasta`: per-sample haplotypes.
- `haplotype_groups.tsv`: haplotype group summary with counts and frequencies.
- `haplotypes_unique.fasta`: unique haplotype sequences.
- `sample_haplotype_calls.tsv`: per-sample haplotype assignments.

## Run examples
### Example with scaffold + positions
```bash
python3 scripts/argweaver/phylogenetic_trees/boxplots/vcf_haplotype_groups_v1.py \
  --vcf Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057694,983057707,983057714,983057718,983057724 \
  --outdir hap_scaffold4
```

### Example using a sites file
```bash
cat > sites.txt <<'SITES'
scaffold_4 983057685
scaffold_4 983057688
scaffold_4 983057694
SITES

python3 scripts/argweaver/phylogenetic_trees/boxplots/vcf_haplotype_groups_v1.py \
  --vcf Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
  --sites_file sites.txt \
  --outdir hap_sites
```

### Example dropping haplotypes with missing data
```bash
python3 scripts/argweaver/phylogenetic_trees/boxplots/vcf_haplotype_groups_v1.py \
  --vcf Imputed_whole_panel_Esteban_Soms_shared.vcf.gz \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057694 \
  --drop_any_missing
```


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
