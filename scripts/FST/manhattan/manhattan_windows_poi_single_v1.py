#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Manhattan plot for window-level Poisson q-values.
- X = cumulative genomic coordinate (bp) across numerically ordered scaffolds
- Y = -log10(q)
- Significant (q < --q) = black dots
- Non-significant       = gray dots
- No labels/annotations

Input must contain at least: CHROM, WIN_START, WIN_END, and one Poisson q column
like q_poi_3vs4 (set with --qcol). CSV or TSV accepted (auto-detected).
"""

import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

def scaffold_sort_key(chrom: str):
    """
    Sort like scaffold_1, scaffold_1a, scaffold_1b, scaffold_2, ...
    Falls back to name if no digits.
    """
    s = str(chrom)
    m = re.search(r'(\d+)([A-Za-z]*)', s)
    if m:
        num = int(m.group(1))
        suf = m.group(2) or ""
        return (num, suf, s)
    return (10**12, "", s)

def autodf(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")

def build_cumpos(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "WIN_MID" not in df.columns:
        df["WIN_MID"] = ((pd.to_numeric(df["WIN_START"], errors="coerce")
                          + pd.to_numeric(df["WIN_END"], errors="coerce")) // 2)
    df["WIN_MID"] = pd.to_numeric(df["WIN_MID"], errors="coerce")

    chroms = sorted(df["CHROM"].dropna().unique().tolist(), key=scaffold_sort_key)
    offsets = {}
    cursor = 0
    for ch in chroms:
        sub = df[df["CHROM"] == ch]
        # use max end as scaffold length proxy
        max_end = pd.to_numeric(sub["WIN_END"], errors="coerce").max()
        max_end = 0 if pd.isna(max_end) else int(max_end)
        offsets[ch] = cursor
        cursor += max_end + 1  # +1 spacer

    df["offset"] = df["CHROM"].map(offsets).fillna(0).astype(np.int64)
    df["cumpos"] = df["offset"] + df["WIN_MID"].fillna(0).astype(np.int64)
    return df

def plot_one(df: pd.DataFrame, qcol: str, qthr: float, title: str, outpng: str, outpdf: str | None = None):
    if qcol not in df.columns:
        raise SystemExit(f"ERROR: q column '{qcol}' not found. Columns: {list(df.columns)[:10]} ...")

    q = pd.to_numeric(df[qcol], errors="coerce")
    # treat NA as non-sig (q=1)
    q = q.fillna(1.0).clip(lower=1e-300, upper=1.0)
    ml10 = -np.log10(q)

    sig = q < qthr
    nonsig = ~sig

    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.scatter(df.loc[nonsig, "cumpos"], ml10.loc[nonsig], s=6, c="#bdbdbd", linewidths=0, alpha=0.7)
    ax.scatter(df.loc[sig,    "cumpos"], ml10.loc[sig],    s=8, c="black",  linewidths=0, alpha=0.8)

    ax.set_ylabel(r"$-\log_{10}(q)$")
    ax.set_xlabel("Genomic coordinate (bp)")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)

    # Make x-axis show Mb for readability, still "genome size" (not scaffold labels)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/1e6:.0f} Mb"))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=12, integer=False))

    fig.tight_layout()
    fig.savefig(outpng, dpi=220)
    if outpdf:
        fig.savefig(outpdf)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Window table (CSV/TSV) with CHROM, WIN_START, WIN_END, q columns")
    ap.add_argument("-o", "--out",   required=True, help="Output image prefix (will save .png and .pdf)")
    ap.add_argument("--qcol",        default="q_poi_3vs4", help="Poisson q-value column to plot (e.g., q_poi_3vs4)")
    ap.add_argument("--q",           type=float, default=0.05, help="Significance threshold (default 0.05)")
    args = ap.parse_args()

    df = autodf(args.input)
    need = ["CHROM", "WIN_START", "WIN_END"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")

    df = df[["CHROM","WIN_START","WIN_END", args.qcol]].copy()
    df["WIN_START"] = pd.to_numeric(df["WIN_START"], errors="coerce")
    df["WIN_END"]   = pd.to_numeric(df["WIN_END"], errors="coerce")
    df["CHROM"]     = df["CHROM"].astype(str)

    df = df.dropna(subset=["WIN_START","WIN_END"])
    df = build_cumpos(df)

    title = f"Windows Manhattan (Poisson q): {args.qcol}  —  q<{args.q}"
    outpng = f"{args.out}.png"
    outpdf = f"{args.out}.pdf"
    plot_one(df, args.qcol, args.q, title, outpng, outpdf)
    print(f"Saved: {outpng}\nSaved: {outpdf}")

if __name__ == "__main__":
    main()

