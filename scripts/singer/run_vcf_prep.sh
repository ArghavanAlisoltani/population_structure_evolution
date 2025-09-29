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
#SBATCH --job-name="singer_prep"

module load bcftools/1.8

# ------------ EDIT THESE ------------
SRC_VCF="/scratch/arghavan/LP/subsampling/s100_All_1a1b_renamed.filtered.vcf.gz"  # .vcf or .vcf.gz
OUTDIR="/scratch/arghavan/LP/singer/prep"

# ------------------------------------

mkdir -p "${OUTDIR}"

echo "[`date`] Input VCF: $SRC_VCF"
# Normalize input -> bgzipped + indexed for bcftools queries
if [[ "$SRC_VCF" =~ \.vcf\.gz$ ]]; then
  SRC_GZ="$SRC_VCF"
  bcftools index -f "$SRC_GZ"
elif [[ "$SRC_VCF" =~ \.vcf$ ]]; then
  echo "[`date`] bgzip + index input .vcf -> .vcf.gz"
  SRC_GZ="${OUTDIR}/source.vcf.gz"
  bgzip -c "$SRC_VCF" > "$SRC_GZ"
  bcftools index -f "$SRC_GZ"
else
  echo "ERROR: SRC_VCF must end with .vcf or .vcf.gz" >&2
  exit 1
fi

# 1) Grab the original header (text)
HDR_ORIG="${OUTDIR}/header.original.txt"
bcftools view -h "$SRC_GZ" > "$HDR_ORIG"

# 2) Build a "base" header WITHOUT any ##contig lines and WITHOUT scaffold_1
#    (We’ll reinsert fresh ##contig lines with computed lengths.)
HDR_BASE="${OUTDIR}/header.base.txt"
grep -v '^##contig=' "$HDR_ORIG" > "$HDR_BASE"

# 3) Compute max position (length proxy) for each contig EXCEPT scaffold_1
#    Using a single pass over variant lines (fast and memory-light).
#    Output: contig\tmaxpos  (only contigs that actually have records)
LEN_TSV="${OUTDIR}/contig_lengths.tsv"
bcftools view -H "$SRC_GZ" \
| awk '$1!="scaffold_1"{ if($2>a[$1]) a[$1]=$2 } END{ for(c in a){ print c"\t"a[c] } }' \
| sort -k1,1V \
> "$LEN_TSV"

if [[ ! -s "$LEN_TSV" ]]; then
  echo "ERROR: Did not find any variant records (after excluding scaffold_1)." >&2
  exit 1
fi

# 4) Turn the TSV into fresh ##contig lines with length=
HDR_CONTIG="${OUTDIR}/header.contigs.txt"
awk -v OFS="" '{print "##contig=<ID=",$1,",length=",$2,">"}' "$LEN_TSV" > "$HDR_CONTIG"

# 5) Stitch new header = base header + fresh contig lines
HDR_NEW="${OUTDIR}/header.rebuild.txt"
cat "$HDR_BASE" "$HDR_CONTIG" > "$HDR_NEW"

# 6) Reheader the bgzipped VCF; write a clean bgzipped copy + index
REHEADERED_GZ="${OUTDIR}/reheadered.vcf.gz"
bcftools reheader -h "$HDR_NEW" -o "$REHEADERED_GZ" "$SRC_GZ"
bcftools index -f "$REHEADERED_GZ"

# 7) Also write a plain .vcf for SINGER (SINGER prefers uncompressed .vcf)
SINGER_VCF="${OUTDIR}/singer_input.vcf"
bcftools view -Ov -o "$SINGER_VCF" "$REHEADERED_GZ"

# 8) Manifest
{
  echo "SOURCE_GZ=$SRC_GZ"
  echo "REHEADERED_GZ=$REHEADERED_GZ"
  echo "SINGER_VCF=$SINGER_VCF"
  echo "LEN_TSV=$LEN_TSV"
} > "${OUTDIR}/manifest.txt"

echo "[`date`] Done."
echo "  - Reheadered (bgz): $REHEADERED_GZ (indexed)"
echo "  - Plain VCF for SINGER: $SINGER_VCF"
echo "  - Contig lengths (from data): $LEN_TSV"
