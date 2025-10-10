#!/bin/bash
#SBATCH -p bigmem
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks-per-node=8
#SBATCH --mail-user=arghavana85@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --output="/scratch/arghavan/LP/FSTcalc/proc/fst_100_greedy/out.txt"
#SBATCH --error="/scratch/arghavan/LP/FSTcalc/proc/fst_100_greedy/err.txt"
#SBATCH --job-name="FSTallLP"

module load vcftools/0.1.16
module load bcftools

# === INPUTS ===
VCF="/scratch/arghavan/LP/vcf_v1/Aria_processed/S10b_1ab/renamed_1ab_samples_1489_merged_sorted_tworules.vcf.gz"
POPDIR="lists_100_greedy_perproc/"
OUTDIR="fst_out"
WINDOW_SIZE=50000
WINDOW_STEP=10000

# filtering thresholds (applied AFTER subsetting each population)
CALLRATE=0.90            # keep variants with >=90% non-missing calls (per-site)
MAF_POP=0.01             # per-population MAF threshold (removes monomorphic)
# (biallelic restriction is enforced explicitly)

mkdir -p "$OUTDIR"/{vcf_by_pop,fst_per_site,fst_windowed,merged_for_pairs,tmp}

echo "[*] CSI-indexing the main VCF"
#bcftools index -f -c "$VCF"

# --- 1) Find population files (one sample ID per line) ---
mapfile -t POP_FILES < <(find "$POPDIR" -maxdepth 1 -type f -name "proc*_samples.txt" | sort)
if (( ${#POP_FILES[@]} == 0 )); then
  echo "No population files found in $POPDIR (pattern: proc*_samples.txt)"; exit 1
fi
echo "[*] Found ${#POP_FILES[@]} population files:"
printf '  - %s\n' "${POP_FILES[@]}"

# --- 2) Create per-population VCFs (subset by sample list) ---
for POP in "${POP_FILES[@]}"; do
  base=$(basename "$POP" .txt)
  rawvcf="$OUTDIR/vcf_by_pop/${base}.vcf.gz"
  echo "[*] Subsetting VCF for $base"
  bcftools view -S "$POP" -Oz -o "$rawvcf" "$VCF"
  bcftools index -f -c "$rawvcf"

  # --- 2b) FILTER the per-pop VCF ---
  #  - keep biallelic SNPs only
  #  - per-site callrate >= 0.90
  #  - per-pop MAF >= 0.01 (also drops monomorphic: MAF=0)
  echo "[*] Filtering (biallelic SNPs, callrate >=${CALLRATE}, MAF >=${MAF_POP}) for $base"
  biavcf="$OUTDIR/tmp/${base}.biallelic.snps.vcf.gz"
  bcftools view -m2 -M2 -v snps -Oz -o "$biavcf" "$rawvcf"
  bcftools index -f -c "$biavcf"

  # vcftools needs an uncompressed recode step; then re-compress
  vcftools --gzvcf "$biavcf" \
           --max-missing "$CALLRATE" \
           --maf "$MAF_POP" \
           --recode --recode-INFO-all \
           --out "$OUTDIR/tmp/${base}.filtered"

  bgzip -f "$OUTDIR/tmp/${base}.filtered.recode.vcf"
  tabix -f -p vcf "$OUTDIR/tmp/${base}.filtered.recode.vcf.gz"

  # replace the original per-pop VCF with the filtered one (same directory)
  mv -f "$OUTDIR/tmp/${base}.filtered.recode.vcf.gz"  "$OUTDIR/vcf_by_pop/${base}.filtered.vcf.gz"
  mv -f "$OUTDIR/tmp/${base}.filtered.recode.vcf.gz.tbi" "$OUTDIR/vcf_by_pop/${base}.filtered.vcf.gz.tbi"
done

# --- 3) Build all unique pairwise combinations ---
pairs=()
for ((i=0; i<${#POP_FILES[@]}; i++)); do
  for ((j=i+1; j<${#POP_FILES[@]}; j++)); do
    pairs+=("${POP_FILES[$i]}:::${POP_FILES[$j]}")
  done
done

# --- 4) For each pair, compute per-site FST and windowed FST ---
# NOTE: As before, this runs FST from the full-cohort VCF using --weir-fst-pop lists.
# If you prefer to compute FST ONLY on sites that pass per-pop filters in BOTH groups,
# you can switch to using the pair-specific merged, filtered VCFs in step 5.
for P in "${pairs[@]}"; do
  A="${P%%:::*}"; B="${P##*:::}"
  Aname=$(basename "$A" .txt)
  Bname=$(basename "$B" .txt)
  tag="${Aname}_vs_${Bname}"

  echo "[*] FST (per-site) for $tag"
  vcftools --gzvcf "$VCF" \
           --weir-fst-pop "$A" \
           --weir-fst-pop "$B" \
           --out "$OUTDIR/fst_per_site/${tag}"

  echo "[*] FST (windowed) for $tag  (win=${WINDOW_SIZE}, step=${WINDOW_STEP})"
  vcftools --gzvcf "$VCF" \
           --weir-fst-pop "$A" \
           --weir-fst-pop "$B" \
           --fst-window-size "$WINDOW_SIZE" \
           --fst-window-step "$WINDOW_STEP" \
           --out "$OUTDIR/fst_windowed/${tag}"
done

# --- 5) (Optional) Run FST on merged, *filtered* per-pop VCFs instead of the full-cohort VCF ---
# Set DO_MERGED=true if you want FST computed only on intersecting, filtered loci.
DO_MERGED=false
if $DO_MERGED; then
  for P in "${pairs[@]}"; do
    A="${P%%:::*}"; B="${P##*:::}"
    Aname=$(basename "$A" .txt)
    Bname=$(basename "$B" .txt)
    tag="${Aname}_vs_${Bname}"

    Avcf="$OUTDIR/vcf_by_pop/${Aname}.filtered.vcf.gz"
    Bvcf="$OUTDIR/vcf_by_pop/${Bname}.filtered.vcf.gz"
    Mvcf="$OUTDIR/merged_for_pairs/${tag}.filtered.merged.vcf.gz"

    echo "[*] Merging FILTERED sub-VCFs for $tag"
    bcftools merge -Oz -o "$Mvcf" "$Avcf" "$Bvcf"
    bcftools index -f -c "$Mvcf"

    echo "[*] FST (per-site) on merged FILTERED VCF for $tag"
    vcftools --gzvcf "$Mvcf" \
             --weir-fst-pop "$A" \
             --weir-fst-pop "$B" \
             --out "$OUTDIR/fst_per_site/${tag}_MERGED_FILTERED"

    echo "[*] FST (windowed) on merged FILTERED VCF for $tag"
    vcftools --gzvcf "$Mvcf" \
             --weir-fst-pop "$A" \
             --weir-fst-pop "$B" \
             --fst-window-size "$WINDOW_SIZE" \
             --fst-window-step "$WINDOW_STEP" \
             --out "$OUTDIR/fst_windowed/${tag}_MERGED_FILTERED"
  done
fi

echo "[✓] All done. Per-pop filtered VCFs: $OUTDIR/vcf_by_pop/*.filtered.vcf.gz"
echo "[✓] Per-site FST: $OUTDIR/fst_per_site/*.weir.fst ;  Windowed: $OUTDIR/fst_windowed/*.windowed.weir.fst"
