# Inputs / outputs
IN="readable.TEanno.cds.tsv"
OUT="readable.TEanno.cds.scaf1ab.tsv"

# Parameters
OFFSET=1427634029

awk -v OFS="\t" \
    -v CHR="scaffold_1" -v A="scaffold_1a" -v B="scaffold_1b" \
    -v OFFSET="$OFFSET" '
BEGIN { shift = OFFSET - 1 }

NF==0 { next }                        # skip empty lines

$1 == CHR {
  start = $4 + 0
  end   = $5 + 0

  if (start >= OFFSET) {
    # move to scaffold_1b and shift coordinates
    $1 = B
    $4 = start - shift
    $5 = end   - shift
  } else {
    # keep on scaffold_1a, coordinates unchanged
    $1 = A
  }
  print; next
}

# already-split or other scaffolds: pass through unchanged
{ print }
' "$IN" > "$OUT"
