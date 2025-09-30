#!/usr/bin/env bash
set -euo pipefail

# --------- Config (edit as needed) ---------
IN="${1:-poly_s100_All_1a1b_renamed.vcf}"     # input VCF
OUT="${2:-split_poly_s100_scaffolds.vcf}"     # output VCF
BP="${3:-breakpoints.txt}"                    # 2 cols: scaffold <TAB> breakpoint
QC="${4:-breakpoints.qc.tsv}"                 # optional; if present, used for contig lengths

# Notes:
#   - If QC is not present/readable, lengths are omitted in ##contig lines (still valid).
#   - Use the Option A R script to generate breakpoints.qc.tsv & breakpoints.txt first.

if [[ ! -r "$BP" ]]; then
  echo "ERROR: Breakpoints file not found/readable: $BP" >&2
  exit 1
fi

awk -v OFS="\t" -v BP="$BP" -v QC="$QC" '
function ends_with_alpha(s,    c) {
  c = substr(s, length(s), 1)
  return (c ~ /[A-Za-z]/)
}
function makeA(s) { return s "a" }
function makeB(s) { return s "b" }

function load_breakpoints(bpfile,    line,sc,bp) {
  while ((getline line < bpfile) > 0) {
    if (line ~ /^[[:space:]]*$/) continue
    split(line, t, /\t/)
    sc = t[1]; bp = t[2] + 0
    if (sc != "" && bp > 0) {
      brk[sc] = bp
      # Precompute derived names
      Aname[sc] = makeA(sc)
      Bname[sc] = makeB(sc)
      need_contig[Aname[sc]] = 1
      need_contig[Bname[sc]] = 1
      contig_done[sc] = 0
    }
  }
  close(bpfile)
}

# Try to load lengths from QC file (len_a / len_b)
# QC header expected to have columns: scaffold, len_a, len_b (among others)
function load_lengths(qcfile,   hdr,idx,i,line,cols,sca,la,lb) {
  if (qcfile == "" || system("[ -r \"" qcfile "\" ]") != 0) {
    return
  }
  # read header
  if ((getline hdr < qcfile) <= 0) {
    close(qcfile); return
  }
  n = split(hdr, cols, /\t/)
  for (i=1;i<=n;i++) {
    if (cols[i] == "scaffold") idx["scaffold"]=i
    else if (cols[i] == "len_a") idx["len_a"]=i
    else if (cols[i] == "len_b") idx["len_b"]=i
  }
  if (!("scaffold" in idx) || !("len_a" in idx) || !("len_b" in idx)) {
    # header not as expected—skip lengths
    close(qcfile); return
  }
  while ((getline line < qcfile) > 0) {
    split(line, t, /\t/)
    sca = t[idx["scaffold"]]
    la  = t[idx["len_a"]] + 0
    lb  = t[idx["len_b"]] + 0
    if (sca in brk) {
      lenA[sca] = la
      lenB[sca] = lb
      have_len[sca] = 1
    }
  }
  close(qcfile)
}

BEGIN {
  load_breakpoints(BP)
  load_lengths(QC)
  injected = 0
}

# We suppress original contig lines for scaffolds we split, and print two new contigs instead.
# If lengths are known from QC, emit length=; otherwise omit it.
/^##contig=<ID=/ {
  # Check if this contig is one we split
  match($0, /ID=([^,>]+)/, m)
  if (m[1] in brk) {
    sc = m[1]
    if (!contig_done[sc]) {
      a = Aname[sc]; b = Bname[sc]
      if (have_len[sc]) {
        printf("##contig=<ID=%s,length=%d>\n", a, lenA[sc])
        printf("##contig=<ID=%s,length=%d>\n", b, lenB[sc])
      } else {
        printf("##contig=<ID=%s>\n", a)
        printf("##contig=<ID=%s>\n", b)
      }
      contig_done[sc] = 1
    }
    next
  }
  # Non-split contig → pass through
  print; next
}

# Right before #CHROM, inject any missing new contigs (some VCFs lack per-contig headers).
/^#CHROM/ {
  if (!injected) {
    for (sc in brk) {
      if (!contig_done[sc]) {
        a = Aname[sc]; b = Bname[sc]
        if (have_len[sc]) {
          printf("##contig=<ID=%s,length=%d>\n", a, lenA[sc])
          printf("##contig=<ID=%s,length=%d>\n", b, lenB[sc])
        } else {
          printf("##contig=<ID=%s>\n", a)
          printf("##contig=<ID=%s>\n", b)
        }
        contig_done[sc] = 1
      }
    }
    injected = 1
  }
  print; next
}

# Other headers
/^#/ { print; next }

# Body lines: rewrite CHROM and POS based on breakpoints
{
  chrom = $1
  pos   = $2 + 0

  if (chrom in brk) {
    bp = brk[chrom]
    if (pos < bp) {
      $1 = Aname[chrom]            # keep POS unchanged
    } else {
      $1 = Bname[chrom]
      $2 = (pos - bp) + 1          # shift so breakpoint maps to POS=1
      if ($2 < 1) $2 = 1           # safety
    }
  }
  print
}
' "$IN" > "$OUT"

echo "Done → $OUT"
