# never_used_prep_genotype_pipeline_v1.sh

## Purpose
Draft SLURM pipeline (not used in production) for preparing genotype data: filter VCFs, convert to PLINK2 formats, apply missingness and per-provenance MAF filters, thin by LD, and output pruned datasets and log files.【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L1-L33】【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L40-L120】

## Key inputs
- SLURM resources and output paths defined in the header (`#SBATCH` directives and `OUTDIR`).【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L1-L33】
- Input VCF (`VCF_IN`), metadata table (`META`) with IID/PROC columns, and thresholds for missingness, MAF, LD pruning, and KING relatedness.【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L13-L33】【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L64-L120】

## How to run
1. Ensure dependencies are available on the cluster: `bcftools`, `tabix`, `plink2`, and `king` (used later in the script).【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L40-L120】
2. Update the user input block for your VCF, metadata, output directory, thread count, and thresholds (`MAX_MISSING_SNP`, `MAX_MISSING_IND`, `PER_PROC_MAF`, LD parameters).【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L13-L33】【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L64-L120】
3. Submit to SLURM:
   ```bash
   sbatch scripts/never_used_prep_genotype_pipeline_v1.sh
   ```
4. The script logs progress to `${OUTDIR}/logs/run.log`, performs sequential filtering (biallelic-only, missingness, per-proc MAF intersection), LD pruning, and prepares pruned PGEN outputs; adjust or remove KING-relatedness steps if unused.【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L40-L120】

## Notes
- File names and thresholds are placeholders from the draft; review each step before adopting in production.
- The script builds intermediate SNP lists per provenance to enforce within-group MAF, then intersects them to keep SNPs common to all groups.【F:scripts/never_used_prep_genotype_pipeline_v1.sh†L64-L99】


## Additional notes
- These steps assume paths and filenames can be adjusted to match your environment.
- When re-running, consider versioning outputs (e.g., suffixes) to avoid overwriting prior results.
- Record software versions and key parameters alongside outputs for reproducibility.
