#!/usr/bin/env bash
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --mem-per-cpu=256G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=arghavana85@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --output="/scratch/arghavan/LP/singer/scaffold_1a_single/out.txt"
#SBATCH --error="/scratch/arghavan/LP/singer/scaffold_1a_single/err.txt"
#SBATCH --job-name="singer_scaffold_1a_single"

# ---- Env ----
export PATH="$PATH:/scratch/arghavan/mybin/singer-0.1.8-beta-linux-x86_64/"
module load bcftools
module load miniconda3/4.7.12.1
source activate /scratch/arghavan/mybin/python3.9

export LC_ALL=C

############################################
# 1) Inputs & run parameters (EDIT HERE)
############################################
# Accepts .vcf.gz or .vcf; we’ll always create a plain .vcf for SINGER.
SRC_VCF="/scratch/arghavan/LP/subsampling/"

CHR="scaffold_1a"
WIN=1000000          # window size (bp); bump to 2e6–5e6 if too sparse
MIN_SNPS=50          # skip windows with fewer than this many SNPs

# SINGER settings
MUT="1.6e-8"
NE="2.0e5"
RATIO="1"
NSAMP="50"           # be conservative per window; raise later if stable
THIN="20"
POLAR="0.5"

OUTDIR="/scratch/arghavan/LP/singer/scaffold_1a_single/scaffold_1a_windows"
mkdir -p "${OUTDIR}"/{logs,vcf,bed,windows,trees}

LOG_MAIN="${OUTDIR}/run.log"
exec > >(tee -a "$LOG_MAIN") 2>&1

############################################
# 2) Prepare VCFs
############################################
# Ensure index for gz source (used by bcftools slicing/stats)
if [[ "$SRC_VCF" =~ \.vcf\.gz$ ]]; then
  bcftools index -f "$SRC_VCF"
  SRC_GZ="$SRC_VCF"
elif [[ "$SRC_VCF" =~ \.vcf$ ]]; then
  # If plain .vcf was provided, bgzip + index it so we can query fast
  echo "[`date`] bgzip + index source VCF for queries..."
  SRC_GZ="${OUTDIR}/$(basename "${SRC_VCF}").gz"
  bgzip -c "$SRC_VCF" > "$SRC_GZ"
  bcftools index -f "$SRC_GZ"
else
  echo "ERROR: SRC_VCF must be .vcf.gz or .vcf (got: $SRC_VCF)"; exit 1
fi

# Always create a plain .vcf that SINGER will read
WORK_VCF="${OUTDIR}/vcf/$(basename "${SRC_GZ%.gz}")"   # *.vcf
if [[ ! -s "$WORK_VCF" ]]; then
  echo "[`date`] Writing plain .vcf for SINGER: $WORK_VCF"
  # (Optional) restrict to CHR here to make WORK_VCF smaller:
  bcftools view -r "$CHR" -Ov -o "$WORK_VCF" "$SRC_GZ"
else
  echo "[`date`] Using existing plain .vcf for SINGER: $WORK_VCF"
fi

############################################
# 3) Preflight checks
############################################
# Sample count
NSAMPLES=$(bcftools query -l "$SRC_GZ" | wc -l | awk '{print $1}')
echo "[`date`] Samples in VCF: $NSAMPLES"

# Rough phasing check (look for '|' in GT); SINGER prefers phased data
PHASED_PROP=$(bcftools query -r "$CHR" -f '[%GT\n]' "$SRC_GZ" 2>/dev/null \
  | head -n 200000 | awk 'BEGIN{t=0;p=0} {t++; if(index($0,"|")>0)p++} END{if(t==0)print 0; else print p/t}')
echo "[`date`] Rough phased GT proportion (first ~200k genotypes on $CHR): $PHASED_PROP"

# Contig length: header length if present, else max position from gz
CONTIG_LEN=$(bcftools view -h "$SRC_GZ" \
  | awk -v c="$CHR" 'match($0,/^##contig=<ID=([^,>]+)(,length=([0-9]+))?/,a){if(a[1]==c && a[3]!=""){print a[3]; exit}}')
if [[ -z "${CONTIG_LEN:-}" ]]; then
  echo "[`date`] Header lacks length for $CHR; inferring from variants..."
  CONTIG_LEN=$(bcftools view -r "$CHR" -H "$SRC_GZ" | awk 'max<$2{max=$2} END{print max+0}')
  if [[ -z "$CONTIG_LEN" || "$CONTIG_LEN" -eq 0 ]]; then
    echo "ERROR: Could not infer length for $CHR (no records?)."; exit 1
  fi
fi
echo "[`date`] $CHR length = $CONTIG_LEN bp"

############################################
# 4) Build windows
############################################
BED="${OUTDIR}/bed/${CHR}_${WIN}.bed"
awk -v c="$CHR" -v L="$CONTIG_LEN" -v W="$WIN" 'BEGIN{OFS="\t"}{for(s=1; s<=L; s+=W){e=s+W-1; if(e>L)e=L; print c,s,e}}' > "$BED"
NWIN=$(wc -l < "$BED" | awk '{print $1}')
echo "[`date`] Windows to process: $NWIN (size ~"$(printf "%'d" $WIN)" bp, min SNPs $MIN_SNPS)"

############################################
# 5) Run SINGER window-by-window (sequential)
############################################
WIN_I=0
while read -r C S E; do
  ((WIN_I++))
  PREFIX="${OUTDIR}/windows/${C}.${S}-${E}"
  LOG_W="${OUTDIR}/logs/${C}.${S}-${E}.log"

  # Count SNPs in this window; skip if too few
  NSNPS=$(bcftools view -r "${C}:${S}-${E}" -H "$SRC_GZ" | wc -l | awk '{print $1}')
  if (( NSNPS < MIN_SNPS )); then
    echo "[`date`] [$WIN_I/$NWIN] Skipping ${C}:${S}-${E} — only $NSNPS SNPs (< $MIN_SNPS)."
    continue
  fi

  echo "[`date`] [$WIN_I/$NWIN] SINGER ${C}:${S}-${E} with $NSNPS SNPs..."
  # Run SINGER (the singer_master shell calls the compiled 'singer' with auto-debug)
  set +e
  singer_master \
    -m "$MUT" \
    -vcf "$WORK_VCF" \
    -output "$PREFIX" \
    -start "$S" \
    -end "$E" \
    -Ne "$NE" \
    -ratio "$RATIO" \
    -n "$NSAMP" \
    -thin "$THIN" \
    -polar "$POLAR" \
    > "$LOG_W" 2>&1
  SINGER_RC=$?
  set -e

  # Verify that at least the first sample’s files exist
  NODES0="${PREFIX}_nodes_0.txt"
  BRANCH0="${PREFIX}_branches_0.txt"
  MUT0="${PREFIX}_mutations_0.txt"

  if (( SINGER_RC != 0 )) || [[ ! -s "$NODES0" || ! -s "$BRANCH0" || ! -s "$MUT0" ]]; then
    echo "[`date`] [$WIN_I/$NWIN] SINGER failed for ${C}:${S}-${E} (rc=$SINGER_RC or missing outputs). See $LOG_W. Skipping conversion."
    continue
  fi

  # Convert all posterior samples for this window to a single .trees
  echo "[`date`] [$WIN_I/$NWIN] Converting to tskit..."
  set +e
  convert_to_tskit \
    -input "$PREFIX" \
    -output "${OUTDIR}/trees/${C}.${S}-${E}" \
    -start 0 -end $((NSAMP-1)) -step 1 \
    >> "$LOG_W" 2>&1
  CONV_RC=$?
  set -e

  if (( CONV_RC != 0 )) || [[ ! -s "${OUTDIR}/trees/${C}.${S}-${E}.trees" ]]; then
    echo "[`date`] [$WIN_I/$NWIN] Conversion failed for ${C}:${S}-${E} (rc=$CONV_RC). See $LOG_W."
  else
    echo "[`date`] [$WIN_I/$NWIN] DONE: ${OUTDIR}/trees/${C}.${S}-${E}.trees"
  fi

done < "$BED"

echo "[`date`] All windows processed. Trees (per-window) in: ${OUTDIR}/trees"
echo "[`date`] Tip: once you have many per-window .trees, we can stitch adjacent windows with tskit in Python (optional)."
