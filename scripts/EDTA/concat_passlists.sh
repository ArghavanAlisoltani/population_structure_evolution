#!/usr/bin/env bash
# concat_passlists.sh
# Usage: bash concat_passlists.sh list.txt all_insertion_times.concat.tsv

set -euo pipefail

LIST_FILE="${1:-list.txt}"
OUT="${2:-all_insertion_times.concat.tsv}"

# start fresh
: > "$OUT"

first=1
while IFS= read -r f || [[ -n "${f:-}" ]]; do
  # skip blanks and comment lines
  [[ -z "${f// }" || "${f:0:1}" == "#" ]] && continue
  [[ ! -f "$f" ]] && { echo "WARN: missing file: $f" >&2; continue; }

  # pick a reader (supports .gz too)
  if [[ "$f" == *.gz ]]; then
    reader="zcat"
  else
    reader="cat"
  fi

  if (( first )); then
    # write header once; normalize "5_TSD 3_TSD" and add SourceFile column
    $reader "$f" | awk 'BEGIN{FS=OFS="\t"}
      NR==1 { gsub(/5_TSD 3_TSD/,"5_TSD\t3_TSD"); print $0, "SourceFile"; next }
      { print $0, src }' src="$f" >> "$OUT"
    first=0
  else
    # append rows from next files, skipping their header
    $reader "$f" | awk 'BEGIN{FS=OFS="\t"} NR>1 { print $0, src }' src="$f" >> "$OUT"
  fi
done < "$LIST_FILE"

echo "Done: wrote $(wc -l < "$OUT") lines to $OUT"

