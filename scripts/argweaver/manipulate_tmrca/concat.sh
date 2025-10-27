#!/usr/bin/env bash
set -euo pipefail

LIST="list.txt"
OUT="all.tmrca.txt"

# empty/initialize output
: > "$OUT"

first_file=1  # set to 0 if you do NOT want to keep only the first header

while IFS= read -r f; do
  # skip blanks and comments
  [[ -z "${f// }" || "$f" =~ ^[[:space:]]*# ]] && continue

  if [[ -f "$f" ]]; then
    if (( first_file )); then
      cat -- "$f" >> "$OUT"      # keep header from the very first file (if present)
      first_file=0
    else
      # append but drop the first line (useful if each file repeats the same header)
      awk 'NR>1' "$f" >> "$OUT"
      # If there is NO header to drop, replace the line above with: cat -- "$f" >> "$OUT"
    fi
  else
    echo "WARNING: missing file: $f" >&2
  fi
done < "$LIST"

echo "Wrote: $OUT"
