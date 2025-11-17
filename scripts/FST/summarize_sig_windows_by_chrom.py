'''
Run Example
python summarize_sig_windows_by_chrom.py \
  -i fst_window_counts_tests_win1000000_ov200000.txt \
  -o sig_window_counts_by_CHROM.tsv \
  --q 0.05
'''

#!/usr/bin/env python3
import argparse, sys
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True, help="Input TSV with window stats")
    ap.add_argument("-o","--output", required=True, help="Output TSV summary")
    ap.add_argument("--q", type=float, default=0.05, help="q-value threshold (default 0.05)")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    if "CHROM" not in df.columns:
        sys.exit("Input must contain CHROM column.")

    # Which q_poi_* columns are present?
    cand = ["q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"]
    comps = [c for c in cand if c in df.columns]
    if not comps:
        sys.exit("No q_poi_* columns found (expected one or more of q_poi_3vs4, q_poi_3vs5, q_poi_4vs5).")

    # Coerce q columns numeric; non-numeric -> NaN -> treated as non-significant
    for c in comps:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Boolean significance per comparison
    for c in comps:
        df[f"sig_{c}"] = df[c] <= args.q

    sig_cols = [f"sig_{c}" for c in comps]
    # How many comparisons significant at each row (0..len(comps))
    df["sig_count"] = df[sig_cols].sum(axis=1)

    # Summaries per CHROM
    gb = df.groupby("CHROM", dropna=False)

    out = gb.apply(lambda g: pd.Series({
        "n_windows": len(g),
        **{f"n_sig_{c.replace('q_poi_','')}": int(g[f"sig_{c}"].sum()) for c in comps},
        "n_sig_any": int((g["sig_count"] >= 1).sum()),
        "n_sig_atleast2": int((g["sig_count"] >= 2).sum()) if len(comps) >= 2 else 0,
        "n_sig_all": int((g["sig_count"] == len(comps)).sum()) if len(comps) >= 2 else 0
    })).reset_index()

    # Optional overall totals row
    tot = pd.DataFrame([{
        "CHROM": "__TOTAL__",
        "n_windows": int(len(df)),
        **{f"n_sig_{c.replace('q_poi_','')}": int(df[f"sig_{c}"].sum()) for c in comps},
        "n_sig_any": int((df["sig_count"] >= 1).sum()),
        "n_sig_atleast2": int((df["sig_count"] >= 2).sum()) if len(comps) >= 2 else 0,
        "n_sig_all": int((df["sig_count"] == len(comps)).sum()) if len(comps) >= 2 else 0
    }])

    out = pd.concat([out, tot], ignore_index=True)
    out.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()

