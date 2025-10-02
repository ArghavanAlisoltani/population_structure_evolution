#!/bin/bash
#SBATCH -p bigmem              # check below for Different Queues
#SBATCH -t 10:00:00             # Walltime/duration of the job
#SBATCH -N 1                   # Number of Nodes
#SBATCH --mem-per-cpu=2G       # Memory per node in GB needed for a job. Also s$
#SBATCH --ntasks-per-node=8    # Number of Cores (Processors)
#SBATCH --mail-user=arghavana85@gmail.com  # Designate e$
#SBATCH --mail-type=FAIL     # Events options are job BEGIN, END, NONE, FAIL, R$
#SBATCH --output="/scratch/arghavan/LP/FSTcalc/proc/out.txt"  # Path for output must alread$
#SBATCH --error="/scratch/arghavan/LP/FSTcalc/proc/err.txt"   # Path for errors must alread$
#SBATCH --job-name="FST_all_LP"       # Name of job

module load vcftools/0.1.16

# === INPUTS ===
VCF="/scratch/arghavan/LP/vcf_v1/Aria_processed/S10_taxa_filter/renamed_samples_1489_merged_sorted_tworules.vcf"                      # your single input VCF/BCF with ALL samples
POPDIR="all_population/"                       # directory containing your population lists (*.txt)
OUTDIR="fst_out_all"                         # all outputs go here
WINDOW_SIZE=50000                        # 50 kb windows (adjust as needed)
WINDOW_STEP=10000                        # 10 kb step (adjust as needed)

mkdir -p "$OUTDIR"/{vcf_by_pop,fst_per_site,fst_windowed,merged_for_pairs}

# --- 0) Index main VCF with CSI (handles long chromosomes/scaffolds) ---
# (works for .vcf.gz or .bcf input)
echo "[*] CSI-indexing the main VCF"
bcftools index -f -c "$VCF"

# --- 1) Find population files (one sample ID per line) ---
#    Uses every *.txt in POPDIR. Your attached lists should match this.
mapfile -t POP_FILES < <(find "$POPDIR" -maxdepth 1 -type f -name "PROC*_one_per_mum__subset100.txt" | sort)
if (( ${#POP_FILES[@]} == 0 )); then
  echo "No population files found in $POPDIR (pattern: PROC*_one_per_mum__subset100.txt)"; exit 1
fi

echo "[*] Found ${#POP_FILES[@]} population files:"
printf '  - %s\n' "${POP_FILES[@]}"

# --- 2) Create per-population VCFs (subset by sample list) ---
# NOTE: vcftools doesn't require this (you can pass --weir-fst-pop files directly),
# but you asked to generate pop-specific VCFs first; we do both (subset + FST from full VCF).
for POP in "${POP_FILES[@]}"; do
  base=$(basename "$POP" .txt)
  outvcf="$OUTDIR/vcf_by_pop/${base}.vcf.gz"
  echo "[*] Subsetting VCF for $base"
  bcftools view -S "$POP" -Oz -o "$outvcf" "$VCF"
  bcftools index -f -c "$outvcf"
done

# --- 3) Build all unique pairwise combinations ---
pairs=()
for ((i=0; i<${#POP_FILES[@]}; i++)); do
  for ((j=i+1; j<${#POP_FILES[@]}; j++)); do
    pairs+=("${POP_FILES[$i]}:::${POP_FILES[$j]}")
  done
done

# --- 4) For each pair, compute per-site FST and windowed FST ---
# Preferred: compute directly from the full cohort VCF using the two --weir-fst-pop files.
# (This preserves identical variant sets and avoids merge artifacts.)
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

# --- 5) (Optional) If you *also* want to run FST on a merged pair-specific VCF: ---
# This is unnecessary, but provided to match "new VCFs before pairwise" literally.
DO_MERGED=false
if $DO_MERGED; then
  for P in "${pairs[@]}"; do
    A="${P%%:::*}"; B="${P##*:::}"
    Aname=$(basename "$A" .txt)
    Bname=$(basename "$B" .txt)
    tag="${Aname}_vs_${Bname}"
    Avcf="$OUTDIR/vcf_by_pop/${Aname}.vcf.gz"
    Bvcf="$OUTDIR/vcf_by_pop/${Bname}.vcf.gz"
    Mvcf="$OUTDIR/merged_for_pairs/${tag}.vcf.gz"

    echo "[*] Merging sub-VCFs for $tag"
    bcftools merge -Oz -o "$Mvcf" "$Avcf" "$Bvcf"
    bcftools index -f -c "$Mvcf"

    echo "[*] FST (per-site) on merged VCF for $tag"
    vcftools --gzvcf "$Mvcf" \
             --weir-fst-pop "$A" \
             --weir-fst-pop "$B" \
             --out "$OUTDIR/fst_per_site/${tag}_MERGED"

    echo "[*] FST (windowed) on merged VCF for $tag"
    vcftools --gzvcf "$Mvcf" \
             --weir-fst-pop "$A" \
             --weir-fst-pop "$B" \
             --fst-window-size "$WINDOW_SIZE" \
             --fst-window-step "$WINDOW_STEP" \
             --out "$OUTDIR/fst_windowed/${tag}_MERGED"
  done
fi

echo "[✓] All done. Per-site FST: $OUTDIR/fst_per_site/*.weir.fst ;  Windowed: $OUTDIR/fst_windowed/*.windowed.weir.fst"
