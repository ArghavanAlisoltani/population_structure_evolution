#!/usr/bin/env python3
'''
Run example
python merge_scaffold_1a_1b_to_scaffold_1.py \
  --in all_tmrca_split.tsv \
  --out all_tmrca_as_scaffold1.tsv \
  --offset 1427634029 \
  --scaf1a scaffold_1a --scaf1b scaffold_1b --scaf-out scaffold_1

'''


import argparse
import sys
import pandas as pd

CAND_CHR = ["CHROM", "chrom", "scaffold", "seqid"]
CAND_START = ["start_tmrca", "start", "Start", "seg_start"]
CAND_END   = ["end_tmrca", "end", "End", "seg_end"]
CAND_LEN   = ["seg_length", "length", "len"]

def pick_col(cols, candidates, required=True, label="column"):
    for c in candidates:
        if c in cols: return c
    if required:
        sys.exit(f"ERROR: Could not find {label}. Looked for: {candidates}")
    return None

def main():
    ap = argparse.ArgumentParser(
        description="Merge scaffold_1a and scaffold_1b back to scaffold_1 using a coordinate offset.")
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV")
    ap.add_argument("--out", dest="out", required=True, help="Output TSV")
    ap.add_argument("--offset", type=int, required=True,
                    help="Coordinate of the FIRST base of 1b on original scaffold_1")
    ap.add_argument("--scaf1a", default="scaffold_1a", help="Name of the 1a scaffold in input (default: scaffold_1a)")
    ap.add_argument("--scaf1b", default="scaffold_1b", help="Name of the 1b scaffold in input (default: scaffold_1b)")
    ap.add_argument("--scaf-out", default="scaffold_1", help="Output scaffold name (default: scaffold_1)")
    ap.add_argument("--offset-is-1based", action="store_true",
                    help="If set, use (OFFSET-1) when adding to 1b rows (for 0-based input).")
    ap.add_argument("--sep", default="\t", help="Input delimiter (default: tab)")
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep=args.sep, dtype=str)
    # Auto-detect columns
    cols = list(df.columns)
    chr_col   = pick_col(cols, CAND_CHR,  True, "chromosome column")
    start_col = pick_col(cols, CAND_START,True, "start column")
    end_col   = pick_col(cols, CAND_END,  True, "end column")
    len_col   = pick_col(cols, CAND_LEN,  False, "length column")

    # Cast to numeric for math
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col]   = pd.to_numeric(df[end_col],   errors="coerce")

    # Work on a copy
    out = df.copy()

    # Compute the delta to add for 1b
    delta = args.offset - 1 if args.offset_is_1based else args.offset

    # Shift 1b rows
    mask_1b = out[chr_col] == args.scaf1b
    out.loc[mask_1b, start_col] = out.loc[mask_1b, start_col] + delta
    out.loc[mask_1b, end_col]   = out.loc[mask_1b, end_col]   + delta
    out.loc[mask_1b, chr_col]   = args.scaf_out

    # Rename 1a rows to scaffold_1 (no shift)
    mask_1a = out[chr_col] == args.scaf1a
    out.loc[mask_1a, chr_col] = args.scaf_out

    # Recompute seg_length if present
    if len_col and len_col in out.columns:
        out[len_col] = (out[end_col] - out[start_col]).astype("Int64")

    # Optional sanity summary to stderr
    def rng(s): 
        return (int(s.min()) if s.notna().any() else None, int(s.max()) if s.notna().any() else None)
    m = out[out[chr_col] == args.scaf_out]
    if not m.empty:
        r_start = rng(m[start_col])
        r_end   = rng(m[end_col])
        print(f"[INFO] {args.scaf_out} after merge: start range {r_start}, end range {r_end}", file=sys.stderr)

    # Save
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] Wrote: {args.out}")

if __name__ == "__main__":
    main()
