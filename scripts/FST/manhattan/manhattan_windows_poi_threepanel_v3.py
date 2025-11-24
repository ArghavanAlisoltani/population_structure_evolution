#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Three-panel Manhattan plot for Poisson q-values per window.

- Panels: default q_poi_3vs4, q_poi_3vs5, q_poi_4vs5 (customizable via --qcols)
- X = cumulative genomic coordinate (bp) across numerically ordered scaffolds
- X-axis labeled in Gb (not scaffold names)
- Scaffolds colored alternately (black / blue)
- Significance by alpha only: --alpha_sig (default 1.0), --alpha_nonsig (default 0.5)
- Orientation: --orientation horizontal|vertical (default: horizontal)
- No red dots, no annotations

Input must contain: CHROM, WIN_START, WIN_END, and the q columns.
CSV or TSV accepted (auto-detected).
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
    Order like: scaffold_1, scaffold_1a, scaffold_1b, scaffold_2, ...
    Falls back to name if no digits.
    """
    s = str(chrom)
    m = re.search(r'(\d+)([A-Za-z]*)', s)
    if m:
        return (int(m.group(1)), m.group(2) or "", s)
    return (10**12, "", s)

def build_cumpos(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    df = df.copy()
    if "WIN_MID" not in df.columns:
        df["WIN_MID"] = ((pd.to_numeric(df["WIN_START"], errors="coerce")
                          + pd.to_numeric(df["WIN_END"],   errors="coerce")) // 2)
    df["WIN_MID"] = pd.to_numeric(df["WIN_MID"], errors="coerce")
    df["WIN_END"] = pd.to_numeric(df["WIN_END"], errors="coerce")

    chroms = sorted(df["CHROM"].dropna().unique().tolist(), key=scaffold_sort_key)

    offsets = {}
    cursor = 0
    for ch in chroms:
        sub = df[df["CHROM"] == ch]
        max_end = sub["WIN_END"].max()
        max_end = 0 if pd.isna(max_end) else int(max_end)
        offsets[ch] = cursor
        cursor += max_end + 1  # spacer

    df["offset"] = df["CHROM"].map(offsets).fillna(0).astype(np.int64)
    df["cumpos"] = df["offset"] + df["WIN_MID"].fillna(0).astype(np.int64)
    return df, chroms

def plot_panels(df: pd.DataFrame, qcols: list[str], qthr: float,
                outprefix: str, size: float,
                orientation: str, alpha_sig: float, alpha_nonsig: float):
    df = df.copy()
    df, chroms = build_cumpos(df)

    # alternate scaffold colors
    alt_colors = ["black", "blue"]
    chrom_color = {ch: alt_colors[i % 2] for i, ch in enumerate(chroms)}

    n_pan = len(qcols)
    if orientation.lower() == "vertical":
        fig, axes = plt.subplots(n_pan, 1, figsize=(8, 3.8 * n_pan), sharex=True, sharey=True)
    else:
        fig, axes = plt.subplots(1, n_pan, figsize=(5.2 * n_pan, 4.6), sharey=True)
    if n_pan == 1:
        axes = [axes]

    for ax, qcol in zip(axes, qcols):
        if qcol not in df.columns:
            raise SystemExit(f"ERROR: q column '{qcol}' not found in input.")

        q = pd.to_numeric(df[qcol], errors="coerce").fillna(1.0).clip(lower=1e-300, upper=1.0)
        ml10 = -np.log10(q)
        sig = (q < qthr)

        # per scaffold: same color, alpha encodes significance
        for ch in chroms:
            mask = (df["CHROM"] == ch)
            if not mask.any():
                continue
            col = chrom_color[ch]

            # non-significant
            idx_ns = mask & (~sig)
            if idx_ns.any():
                ax.scatter(df.loc[idx_ns, "cumpos"], ml10.loc[idx_ns],
                           s=size, c=col, alpha=alpha_nonsig, linewidths=0)

            # significant
            idx_s = mask & sig
            if idx_s.any():
                ax.scatter(df.loc[idx_s, "cumpos"], ml10.loc[idx_s],
                           s=size, c=col, alpha=alpha_sig, linewidths=0)

        ax.set_title(qcol, fontsize=11)
        ax.grid(axis="y", alpha=0.25)
        # X in Gb
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/1e9:.1f} Gb"))

    axes[0].set_ylabel(r"$-\log_{10}(q)$", fontsize=11)
    if orientation.lower() == "vertical":
        axes[-1].set_xlabel("Genomic coordinate (Gb)", fontsize=10)
    else:
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
    ap.add_argument("-o", "--out", required=True, help="Output image prefix (no extension)")
    ap.add_argument("--qcols", default="q_poi_3vs4,q_poi_3vs5,q_poi_4vs5",
                    help="Comma-separated q columns to plot (default 3 panels)")
    ap.add_argument("--q", type=float, default=0.05, help="Significance threshold (default 0.05)")
    ap.add_argument("--size", type=float, default=7.0, help="Dot size (default 7)")
    ap.add_argument("--orientation", choices=["horizontal","vertical"], default="horizontal",
                    help="Panel arrangement (default: horizontal)")
    ap.add_argument("--alpha_sig", type=float, default=1.0, help="Alpha for significant points (default 1.0)")
    ap.add_argument("--alpha_nonsig", type=float, default=0.5, help="Alpha for non-significant points (default 0.5)")
    args = ap.parse_args()

    df = autodf(args.input)
    need = ["CHROM", "WIN_START", "WIN_END"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")

    qcols = [c.strip() for c in args.qcols.split(",") if c.strip()]
    df = df[["CHROM","WIN_START","WIN_END"] + qcols].copy()
    df["CHROM"] = df["CHROM"].astype(str)
    df["WIN_START"] = pd.to_numeric(df["WIN_START"], errors="coerce")
    df["WIN_END"]   = pd.to_numeric(df["WIN_END"],   errors="coerce")
    df = df.dropna(subset=["WIN_START","WIN_END"])

    plot_panels(df, qcols, args.q, args.out, args.size, args.orientation,
                args.alpha_sig, args.alpha_nonsig)

if __name__ == "__main__":
    main()

