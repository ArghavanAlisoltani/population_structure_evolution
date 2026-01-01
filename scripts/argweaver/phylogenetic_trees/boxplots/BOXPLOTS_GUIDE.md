# Boxplots pipeline guide (ArgWeaver phylogenetic trees)

This folder contains small utilities for summarizing genotype data at selected loci and visualizing phenotype distributions by inferred haplotype groups.

## Contents

- `vcf_to_haplotypes.py`: Minimal Python extractor that converts selected VCF positions into two haplotype FASTA sequences per sample.
- `Genotype_counts_per_locus.R`: R script to count genotype categories (AA/AB/BB/missing) per variant and write a TSV.
- `vcf_hap_box_plots_v1.R`: Rscript CLI that builds haplotype groups from selected VCF positions, merges phenotypes, and produces box/jitter plots.

---

## 1) `vcf_to_haplotypes.py`

### What it does
- Reads a VCF text file and extracts genotypes for a **single scaffold** at **user-provided positions**.
- Builds two haplotype strings per sample (phased or unphased GTs are supported).
- Prints FASTA output to stdout (two records per sample: `sample_1`, `sample_2`).

### Key logic
- Parses the `#CHROM` header to collect sample IDs.
- For each matching variant, extracts the `GT` field and maps `0` to `REF` and `1` to `ALT` (missing `.` becomes `N`).
- Uses `|` for phased and `/` for unphased genotypes.

### Configuration (edit in file)
At the bottom of the script:
```python
input_vcf = "vcf.header.txt"
scaffold = "scaffold_17"
positions = [21197, 21238, 21240, 21245, 21297, 21300, 21302]
```

### Example usage
```bash
python scripts/argweaver/phylogenetic_trees/boxplots/vcf_to_haplotypes.py > haplotypes.fasta
```

---

## 2) `Genotype_counts_per_locus.R`

### What it does
- Loads a VCF (optionally a single chromosome/scaffold).
- Extracts genotype calls and recodes them into `AA`, `AB`, `BB`, or missing.
- Writes a per-variant TSV with genotype counts and metadata.

### Key logic
- Uses `vcfR::read.vcfR()` to load VCF data.
- Normalizes phased/unphased GTs and recodes to `AA/AB/BB`.
- Counts called/missing and heterozygotes per variant.

### Configuration (edit in file)
```r
vcf_file <- "Imputed_whole_panel_Esteban_Soms_shared.vcf.gz"
out_tsv  <- "genotype_counts_per_locus.tsv"
only_chrom <- NULL  # or "scaffold_4"
```

### Example usage
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/Genotype_counts_per_locus.R
```

Outputs:
- `genotype_counts_per_locus.tsv` (per-variant summary)

---

## 3) `vcf_hap_box_plots_v1.R`

### What it does
- Ensures the VCF is bgzip/tabix-indexed (creates a local copy if needed).
- Extracts selected positions using `bcftools query`.
- Builds haplotype sequences per sample and groups individuals by haplotype.
- Merges a phenotype table and generates box/jitter plots by haplotype group.

### Key outputs
- `vcf_selected_positions_raw.tsv`: VCF rows for the selected positions.
- `hap_code_to_sequence.tsv`: Mapping of Hxx codes to haplotype sequences.
- `haplotype_counts_unique_individuals.tsv`: Unique-individual counts per haplotype.
- `sample_to_haplotype_groups.tsv`: Sample/haplotype membership.
- `plot_data_<trait>_hapgroups.tsv`: Plot-ready table.
- `box_jitter_<trait>_by_HAPLOTYPES_vcf.(png|pdf)`

### Dependencies
- R packages: `data.table`, `ggplot2`, `stringr`, `optparse`
- CLI tools: `bcftools`, `bgzip`, `tabix` (from htslib)

### Example usage
```bash
Rscript scripts/argweaver/phylogenetic_trees/boxplots/vcf_hap_box_plots_v1.R \
  --vcf data/input.vcf.gz \
  --pheno data/phenotypes.tsv \
  --pheno_sep '\t' \
  --id_col codg \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057701 \
  --traits C13,C14 \
  --outdir results/vcf_hap_boxplots \
  --png TRUE \
  --pdf FALSE
```

### Notes on phenotype matching
- The script creates `sample_id` by **extracting digits** from the phenotype ID column if present.
- If no digits are present, it uses the raw ID string.
- Ensure the phenotype IDs match the VCF sample IDs (or at least share matching numeric identifiers).

---

## Suggested workflow
1. (Optional) Run `Genotype_counts_per_locus.R` to check genotype balance across loci.
2. Use `vcf_hap_box_plots_v1.R` to summarize haplotypes and create phenotype boxplots.
3. If you only need FASTA haplotypes for a small set of positions, run `vcf_to_haplotypes.py`.

