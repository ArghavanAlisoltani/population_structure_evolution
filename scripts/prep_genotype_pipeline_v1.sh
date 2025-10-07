#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=8
#SBATCH --mail-user=arghavana85@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --output="/scratch/arghavan/LP/subsampling/out_prep.txt"
#SBATCH --error="/scratch/arghavan/LP/singer/scaffold_1a_single/err.txt"
#SBATCH --job-name="prep"

# ------------------ USER INPUTS ------------------
VCF_IN="/scratch/arghavan/LP/subsampling/poly_s100_All_1a1b_renamed.vcf.gz"        # your VCF (bgzipped + tabix)
META="samples_metadata.tsv"            # 2 cols: IID<TAB>PROC  (proc/provenance)
OUTDIR="prep_out"
THREADS=8

# Thresholds (edit to match the paper; defaults are sane placeholders)
MAX_MISSING_SNP=0.10     # per-SNP missingness (e.g., 0.10 = keep SNPs with callrate ≥0.90)
MAX_MISSING_IND=0.10     # per-individual missingness
PER_PROC_MAF=0.05        # MAF threshold to enforce within each provenance

# LD thinning (indep-pairwise window_kb step_kb r2)
LD_WIN_KB=50
LD_STEP_KB=5
LD_R2=0.2

# KING relatedness cutoff (0.125 ~ 2nd-degree)
KING_CUTOFF=0.125
# -------------------------------------------------

mkdir -p "${OUTDIR}"/{logs,lists,tmp}

echo "[$(date)] Index VCF if needed..." | tee -a "${OUTDIR}/logs/run.log"
if [ ! -f "${VCF_IN}.tbi" ]; then
  tabix -p vcf "${VCF_IN}"
fi

# 1) Keep strictly biallelic SNPs (drop non-SNPs, multiallelics)
echo "[$(date)] Keep biallelic SNPs only..." | tee -a "${OUTDIR}/logs/run.log"
bcftools view -m2 -M2 -v snps -Oz -o "${OUTDIR}/step1.biallelic.snps.vcf.gz" "${VCF_IN}"
tabix -p vcf "${OUTDIR}/step1.biallelic.snps.vcf.gz"

# 2) Convert to PLINK2 PGEN
echo "[$(date)] Convert to PLINK2..." | tee -a "${OUTDIR}/logs/run.log"
plink2 \
  --vcf "${OUTDIR}/step1.biallelic.snps.vcf.gz" 'dosage=DS' 'vcf-half-call=m' \
  --double-id \
  --threads ${THREADS} \
  --make-pgen \
  --out "${OUTDIR}/step2.biallelic"

# 3) Missingness filtering (per-SNP and per-individual)
echo "[$(date)] Apply missingness filters..." | tee -a "${OUTDIR}/logs/run.log"
plink2 \
  --pfile "${OUTDIR}/step2.biallelic" \
  --geno ${MAX_MISSING_SNP} \
  --mind ${MAX_MISSING_IND} \
  --threads ${THREADS} \
  --make-pgen \
  --out "${OUTDIR}/step3.missingness"

# 4) Per-provenance MAF filter
#    Rule: a SNP is kept only if it has MAF ≥ PER_PROC_MAF in *each* provenance
#    (Change the intersection logic below if you want “kept if MAF≥X in at least one/two provenances”.)
echo "[$(date)] Build per-provenance SNP lists at MAF≥${PER_PROC_MAF}..." | tee -a "${OUTDIR}/logs/run.log"

# Expect META file: IID<TAB>PROC (header allowed)
# Create one keep list per PROC
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2}' "${META}" > "${OUTDIR}/tmp/iid_proc.tsv"
cut -f2 "${OUTDIR}/tmp/iid_proc.tsv" | sort -u > "${OUTDIR}/lists/proc.list"

> "${OUTDIR}/lists/maf_pass.snplist"  # will be overwritten after intersection

while read PROC; do
  awk -v p="${PROC}" '$2==p{print $1,$1}' OFS="\t" "${OUTDIR}/tmp/iid_proc.tsv" > "${OUTDIR}/lists/keep_${PROC}.txt"
  plink2 \
    --pfile "${OUTDIR}/step3.missingness" \
    --keep "${OUTDIR}/lists/keep_${PROC}.txt" \
    --maf ${PER_PROC_MAF} \
    --threads ${THREADS} \
    --write-snplist \
    --out "${OUTDIR}/lists/maf_${PROC}"
done < "${OUTDIR}/lists/proc.list"

# Intersect all per-PROC snp lists (require SNP to pass in every PROC)
echo "[$(date)] Intersect per-PROC MAF-passing SNP lists..." | tee -a "${OUTDIR}/logs/run.log"
cp "$(head -n1 ${OUTDIR}/lists/proc.list | sed "s|^|${OUTDIR}/lists/maf_|;s|$|.snplist|")" "${OUTDIR}/lists/maf_pass.snplist"

tail -n +2 "${OUTDIR}/lists/proc.list" | while read PROC; do
  comm -12 \
    <(sort "${OUTDIR}/lists/maf_pass.snplist") \
    <(sort "${OUTDIR}/lists/maf_${PROC}.snplist") \
    > "${OUTDIR}/lists/maf_pass.tmp"
  mv "${OUTDIR}/lists/maf_pass.tmp" "${OUTDIR}/lists/maf_pass.snplist"
done

# Apply the joint per-PROC MAF filter
plink2 \
  --pfile "${OUTDIR}/step3.missingness" \
  --extract "${OUTDIR}/lists/maf_pass.snplist" \
  --threads ${THREADS} \
  --make-pgen \
  --out "${OUTDIR}/step4.perproc_maf"

# 5) LD thinning (indep-pairwise)
echo "[$(date)] LD thinning (indep-pairwise ${LD_WIN_KB} ${LD_STEP_KB} ${LD_R2})..." | tee -a "${OUTDIR}/logs/run.log"
plink2 \
  --pfile "${OUTDIR}/step4.perproc_maf" \
  --indep-pairwise ${LD_WIN_KB} ${LD_STEP_KB} ${LD_R2} \
  --threads ${THREADS} \
  --out "${OUTDIR}/lists/ld"
# Pruned dataset (for structure/PCs/etc.)
plink2 \
  --pfile "${OUTDIR}/step4.perproc_maf" \
  --extract "${OUTDIR}/lists/ld.prune.in" \
  --threads ${THREADS} \
  --make-pgen \
  --out "${OUTDIR}/step5.ldpruned"

# 6) Relatedness pruning (KING; drop one of each pair with kinship > cutoff)
echo "[$(date)] KING relatedness pruning (cutoff ${KING_CUTOFF})..." | tee -a "${OUTDIR}/logs/run.log"
plink2 \
  --pfile "${OUTDIR}/step5.ldpruned" \
  --king-cutoff ${KING_CUTOFF} \
  --threads ${THREADS} \
  --make-pgen \
  --out "${OUTDIR}/final.ldpruned.king${KING_CUTOFF}"

echo "[$(date)] Done."
echo "Final dataset: ${OUTDIR}/final.ldpruned.king${KING_CUTOFF}.pgen/.pvar/.psam"
