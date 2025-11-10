#!/usr/bin/env bash
# rename_from_map.sh
# Usage:
#   bash rename_from_map.sh path/to/file_name.txt         # real run
#   bash rename_from_map.sh path/to/file_name.txt --dry   # dry-run (no changes)

set -euo pipefail

MAPFILE="${1:-file_name.txt}"
DRY="${2:-}"

if [[ ! -f "$MAPFILE" ]]; then
  echo "Mapping file not found: $MAPFILE" >&2
  exit 1
fi

# Command to perform (mv) or print (dry-run)
do_mv() {
  if [[ "$DRY" == "--dry" || "$DRY" == "--dry-run" ]]; then
    printf 'DRY mv %q -> %q\n' "$1" "$2"
  else
    mv -v -- "$1" "$2"
  fi
}

# Expect TAB-delimited file: OLD<TAB>NEW
# Lines starting with # or blank lines are ignored.
while IFS=$'\t' read -r old new || [[ -n "${old:-}" ]]; do
  # skip empty / commented lines
  [[ -z "${old:-}" || "${old:0:1}" == "#" ]] && continue

  if [[ ! -e "$old" ]]; then
    printf 'MISS: source not found: %s\n' "$old" >&2
    continue
  fi

  if [[ -e "$new" ]]; then
    printf 'SKIP: target exists: %s -> %s\n' "$old" "$new" >&2
    continue
  fi

  # create target directory if needed
  mkdir -p -- "$(dirname -- "$new")"

  do_mv "$old" "$new"
done < "$MAPFILE"
