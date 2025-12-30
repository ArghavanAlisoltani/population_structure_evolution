# pre_processing_scripts.sh

## Purpose
Generate a biallelic SNP-only VCF for a target scaffold, add basic site tags (AN/AC/AF), and index the result for downstream BetaScan-style analyses.

## Key inputs
- `Imputed_whole_panel_Esteban_Soms_shared.vcf.gz`: the source VCF.
- `scaffold_4`: scaffold or region to extract.
- `bcftools`: required for filtering and tagging.

## How to run
```bash
bash scripts/BetaScan/pre_processing_scripts.sh
```

## Notes
- The script writes `scaffold_4.biallelic.tags.vcf.gz` and creates a `.csi` index with `bcftools index`.
- Update the scaffold name and input VCF if you want to process a different region or dataset.
