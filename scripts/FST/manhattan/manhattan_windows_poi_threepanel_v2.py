#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Three-panel Manhattan plot for Poisson q-values per window.

- Panels: q_poi_3vs4, q_poi_3vs5, q_poi_4vs5 (or any you pass via --qcols)
- X = cumulative genomic coordinate across numerically ordered scaffolds
- X-axis labels in Gb (not scaffold names)
- Two colors alternating by scaffold: black / blue
- Significance (q < --q) via alpha: 1.0 (sig), 0.5 (non-sig)
- No labels/annotations for peaks

Input must contain: CHROM, WIN_START, WIN_END, and the q-value columns.
CSV or TSV supported (auto-detected).

Example:
  python manhattan_windows_poi_threepanel_v2.py \
    -i fst_window_counts_tests_win1000000_ov200000.csv \
    -o manhattan_poi_3pan_win1mb \
    --qcols q_poi_3vs4,q_poi_3vs5,q_poi_4vs5 \
    --q 0.05
"""

import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def autodf(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")

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

def build_cumpos(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str], dict]:
    df = df.copy()
    # midpoints
    if "WIN_MID" not in df.columns:
        df["WIN_MID"] = ((pd.to_numeric(df["WIN_START"], errors="coerce")
                          + pd.to_numeric(df["WIN_END"], errors="coerce")) // 2)
    df["WIN_MID"] = pd.to_numeric(df["WIN_MID"], errors="coerce")
    df["WIN_END"] = pd.to_numeric(df["WIN_END"], errors="coerce")

    # scaffold order & offsets
    chroms = sorted(df["CHROM"].dropna().unique().tolist(), key=scaffold_sort_key)
    offsets = {}
    cursor = 0
    for ch in chroms:
        sub = df[df["CHROM"] == ch]
        max_end = sub["WIN_END"].max()
        max_end = 0 if pd.isna(max_end) else int(max_end)
        offsets[ch] = cursor
        cursor += max_end + 1  # +1 spacer

    df["offset"] = df["CHROM"].map(offsets).fillna(0).astype(np.int64)
    df["cumpos"] = df["offset"] + df["WIN_MID"].fillna(0).astype(np.int64)
    return df, chroms, offsets

def plot_three_panels(df: pd.DataFrame, qcols: list[str], qthr: float,
                      outprefix: str, dot_size: float = 7.0):
    # Build cumulative coordinate once
    df = df.copy()
    df, chroms, offsets = build_cumpos(df)

    # color map (alternate by scaffold: black / blue)
    alt_colors = ["black", "blue"]
    chrom_color = {ch: alt_colors[i % 2] for i, ch in enumerate(chroms)}

    # Prepare figure
    n_pan = len(qcols)
    fig, axes = plt.subplots(1, n_pan, figsize=(5.0 * n_pan, 4.6), sharey=True)

    if n_pan == 1:
        axes = [axes]

    for ax, qcol in zip(axes, qcols):
        if qcol not in df.columns:
            raise SystemExit(f"ERROR: q column '{qcol}' not found in input.")

        # Coerce q, compute -log10(q)
        q = pd.to_numeric(df[qcol], errors="coerce").fillna(1.0).clip(lower=1e-300, upper=1.0)
        ml10 = -np.log10(q)
        sig_mask = (q < qthr)

        # Plot per scaffold in alternating colors, alpha by significance
        for ch in chroms:
            sub_idx = (df["CHROM"] == ch)
            if not sub_idx.any():
                continue
            col = chrom_color[ch]

            # non-significant (alpha=0.5)
            idx_ns = sub_idx & (~sig_mask)
            if idx_ns.any():
                ax.scatter(df.loc[idx_ns, "cumpos"], ml10.loc[idx_ns],
                           s=dot_size, c=col, alpha=0.5, linewidths=0)

            # significant (alpha=1.0)
            idx_s = sub_idx & sig_mask
            if idx_s.any():
                ax.scatter(df.loc[idx_s, "cumpos"], ml10.loc[idx_s],
                           s=dot_size, c=col, alpha=1.0, linewidths=0)

        # axes cosmetics
        ax.set_title(qcol, fontsize=11)
        ax.grid(axis="y", alpha=0.25)
        # X in Gb
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/1e9:.1f} Gb"))

    axes[0].set_ylabel(r"$-\log_{10}(q)$", fontsize=11)
    for ax in axes:
        ax.set_xlabel("Genomic coordinate (Gb)", fontsize=10)

    fig.suptitle(f"Window Manhattan (Poisson q), q<{qthr}", y=0.98, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_png = f"{outprefix}.png"
    out_pdf = f"{outprefix}.pdf"
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved: {out_png}\nSaved: {out_pdf}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True,
                    help="Window table (CSV/TSV) with CHROM, WIN_START, WIN_END, q columns")
    ap.add_argument("-o", "--out", required=True, help="Output image prefix (without extension)")
    ap.add_argument("--qcols", default="q_poi_3vs4,q_poi_3vs5,q_poi_4vs5",
                    help="Comma-separated q columns to plot (default: q_poi_3vs4,q_poi_3vs5,q_poi_4vs5)")
    ap.add_argument("--q", type=float, default=0.05, help="Significance threshold (default 0.05)")
    ap.add_argument("--size", type=float, default=7.0, help="Dot size (default 7)")
    args = ap.parse_args()

    df = autodf(args.input)
    need = ["CHROM", "WIN_START", "WIN_END"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")

    qcols = [c.strip() for c in args.qcols.split(",") if c.strip()]
    df = df[["CHROM", "WIN_START", "WIN_END"] + qcols].copy()
    df["CHROM"] = df["CHROM"].astype(str)

    # ensure numeric window coords
    df["WIN_START"] = pd.to_numeric(df["WIN_START"], errors="coerce")
    df["WIN_END"]   = pd.to_numeric(df["WIN_END"], errors="coerce")
    df = df.dropna(subset=["WIN_START", "WIN_END"])

    plot_three_panels(df, qcols, args.q, args.out, dot_size=args.size)

if __name__ == "__main__":
    main()

