#!/usr/bin/env python3
"""
Manhattan plot (LG-ordered) for MTAG/GWAS-style summary stats.

Requirements:
  pip install pandas numpy matplotlib

Example:
  python3 manhattan_lg_mtag.py \
    --infile sumstat.txt \
    --out manhattan_LG.png \
    --threshold 6
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="Input table (TSV by default).")
    ap.add_argument("--out", required=True, help="Output figure path (png/pdf/jpg).")
    ap.add_argument("--sep", default="\t", help="Field separator (default: tab). Use ',' for CSV.")
    ap.add_argument("--p_col", default="P_value", help="P-value column name.")
    ap.add_argument("--trait_col", default="TRAIT", help="Trait column name.")
    ap.add_argument("--lg_col", default="LG_group", help="LG group column name (e.g., LG1..LG12).")
    ap.add_argument("--scaf_col", default="Scaffold", help="Scaffold identifier column name.")
    ap.add_argument("--scaf_order_col", default="Scaffold_order_per_LG", help="Scaffold order within LG column.")
    ap.add_argument("--pos_col", default="POS", help="Position column name (within scaffold).")
    ap.add_argument("--len_col", default="Length", help="Scaffold length column name.")
    ap.add_argument("--threshold", type=float, default=6.0, help="Significance threshold in -log10(p).")
    ap.add_argument("--title", default="", help="Optional plot title.")
    ap.add_argument("--dpi", type=int, default=300, help="Output DPI (default 300).")
    ap.add_argument("--fig_w", type=float, default=18.0, help="Figure width in inches.")
    ap.add_argument("--fig_h", type=float, default=5.5, help="Figure height in inches.")
    ap.add_argument("--bg_size", type=float, default=3.0, help="Point size for background (non-significant) dots.")
    ap.add_argument("--sig_size", type=float, default=18.0, help="Point size for significant dots.")
    ap.add_argument("--max_points_bg", type=int, default=0,
                    help="Optional: downsample background points to this many (0 = no downsampling).")
    return ap.parse_args()


def lg_to_num(x: str) -> int:
    if pd.isna(x):
        return 10**9
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else 10**9


def distinct_colors(n: int):
    # tab20 gives 20 distinct-ish colors; if more, fall back to HSV sampling
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        return [cmap(i) for i in range(n)]
    cmap = plt.get_cmap("hsv")
    return [cmap(i / n) for i in range(n)]


def main():
    a = parse_args()

    usecols = [a.p_col, a.trait_col, a.lg_col, a.scaf_col, a.scaf_order_col, a.pos_col, a.len_col]
    # Read and keep only needed columns (if some are missing, handle below)
    df = pd.read_csv(a.infile, sep=a.sep, dtype=str, low_memory=False)
    missing = [c for c in usecols if c not in df.columns]
    if missing:
        # allow missing Length or Scaffold_order_per_LG; other columns are required
        required = {a.p_col, a.trait_col, a.lg_col, a.scaf_col, a.pos_col}
        miss_req = [c for c in missing if c in required]
        if miss_req:
            raise SystemExit(f"ERROR: missing required columns: {miss_req}\nFound columns: {list(df.columns)}")

    # Convert numeric columns
    df[a.p_col] = pd.to_numeric(df[a.p_col], errors="coerce")
    df[a.pos_col] = pd.to_numeric(df[a.pos_col], errors="coerce")
    if a.scaf_order_col in df.columns:
        df[a.scaf_order_col] = pd.to_numeric(df[a.scaf_order_col], errors="coerce")
    else:
        df[a.scaf_order_col] = np.nan

    if a.len_col in df.columns:
        df[a.len_col] = pd.to_numeric(df[a.len_col], errors="coerce")
    else:
        df[a.len_col] = np.nan

    # Drop rows missing essentials
    df = df.dropna(subset=[a.p_col, a.pos_col, a.lg_col, a.scaf_col])

    # -log10(p) with floor for zeros
    p_floor = 1e-300
    p = df[a.p_col].to_numpy()
    p = np.where(p <= 0, p_floor, p)
    df["mlogp"] = -np.log10(p)

    # LG numeric ordering
    df["lg_num"] = df[a.lg_col].map(lg_to_num).astype(int)

    # If scaffold length missing, estimate as max POS per scaffold within LG
    if df[a.len_col].isna().all():
        est_len = (
            df.groupby(["lg_num", a.scaf_col], as_index=False)[a.pos_col]
              .max()
              .rename(columns={a.pos_col: a.len_col})
        )
        df = df.merge(est_len, on=["lg_num", a.scaf_col], how="left", suffixes=("", "_est"))

    # Prepare unique scaffold table to compute offsets within each LG
    scaf_cols = ["lg_num", a.lg_col, a.scaf_col, a.scaf_order_col, a.len_col]
    scaf_tbl = df[scaf_cols].drop_duplicates()

    # If scaffold order missing, build a stable order by scaffold id
    if scaf_tbl[a.scaf_order_col].isna().all():
        # sort scaffold IDs alphanumerically within each LG
        scaf_tbl = scaf_tbl.sort_values(["lg_num", a.scaf_col]).copy()
        scaf_tbl[a.scaf_order_col] = (
            scaf_tbl.groupby("lg_num").cumcount().astype(int) + 1
        )

    # Make sure lengths exist; if still missing, estimate again
    if scaf_tbl[a.len_col].isna().any():
        # try using max POS as length (per LG, scaffold)
        est_len2 = (
            df.groupby(["lg_num", a.scaf_col], as_index=False)[a.pos_col]
              .max()
              .rename(columns={a.pos_col: a.len_col})
        )
        scaf_tbl = scaf_tbl.drop(columns=[a.len_col]).merge(est_len2, on=["lg_num", a.scaf_col], how="left")

    scaf_tbl[a.len_col] = scaf_tbl[a.len_col].fillna(0).astype(float)

    # Sort scaffolds within LG by scaffold_order_per_LG
    scaf_tbl = scaf_tbl.sort_values(["lg_num", a.scaf_order_col, a.scaf_col]).copy()

    # Offset within LG (cumulative sum of scaffold lengths)
    scaf_tbl["offset_in_lg"] = (
        scaf_tbl.groupby("lg_num")[a.len_col].cumsum() - scaf_tbl[a.len_col]
    )

    # Total length per LG and global offset per LG
    lg_len = scaf_tbl.groupby(["lg_num", a.lg_col], as_index=False)[a.len_col].sum().rename(columns={a.len_col: "lg_len"})
    lg_len = lg_len.sort_values("lg_num").copy()
    lg_len["lg_global_offset"] = lg_len["lg_len"].cumsum() - lg_len["lg_len"]

    # Merge offsets back
    df = df.merge(
        scaf_tbl[["lg_num", a.scaf_col, a.scaf_order_col, "offset_in_lg"]],
        on=["lg_num", a.scaf_col, a.scaf_order_col],
        how="left"
    ).merge(
        lg_len[["lg_num", "lg_global_offset", a.lg_col, "lg_len"]],
        on=["lg_num", a.lg_col],
        how="left"
    )

    # Final x coordinate: within-scaffold pos + scaffold offset + LG offset
    df["x"] = df["lg_global_offset"] + df["offset_in_lg"] + df[a.pos_col].astype(float)

    # Sort (for any optional downsampling or consistent rendering)
    df = df.sort_values(["lg_num", a.scaf_order_col, a.pos_col])

    # Split significant vs background
    sig = df["mlogp"] >= a.threshold
    df_sig = df.loc[sig].copy()
    df_bg = df.loc[~sig].copy()

    # Optional downsampling of background for speed/size
    if a.max_points_bg and a.max_points_bg > 0 and len(df_bg) > a.max_points_bg:
        df_bg = df_bg.sample(n=a.max_points_bg, random_state=1).sort_values("x")

    # Alternate greys by LG parity (like grey75/grey88)
    grey75 = (0.75, 0.75, 0.75, 1.0)
    grey88 = (0.88, 0.88, 0.88, 1.0)

    # Build figure
    fig, ax = plt.subplots(figsize=(a.fig_w, a.fig_h))

    # Background points (no legend)
    bg_even = df_bg["lg_num"] % 2 == 0
    ax.scatter(df_bg.loc[~bg_even, "x"], df_bg.loc[~bg_even, "mlogp"],
               s=a.bg_size, c=[grey88], linewidths=0, alpha=1.0, rasterized=True)
    ax.scatter(df_bg.loc[bg_even, "x"], df_bg.loc[bg_even, "mlogp"],
               s=a.bg_size, c=[grey75], linewidths=0, alpha=1.0, rasterized=True)

    # Significant points colored by trait (legend = traits only)
    if len(df_sig) > 0:
        traits = sorted(df_sig[a.trait_col].fillna("NA").unique().tolist())
        colors = distinct_colors(len(traits))
        trait2color = dict(zip(traits, colors))

        for tr in traits:
            d = df_sig[df_sig[a.trait_col].fillna("NA") == tr]
            ax.scatter(d["x"], d["mlogp"],
                       s=a.sig_size, c=[trait2color[tr]],
                       linewidths=0, alpha=0.95, label=tr, rasterized=True)

    # Threshold line
    ax.axhline(a.threshold, linestyle="--", linewidth=1.2, alpha=0.9)

    # X ticks at centers of each LG region, labels "LG1..LG12"
    lg_len["lg_start"] = lg_len["lg_global_offset"]
    lg_len["lg_end"] = lg_len["lg_global_offset"] + lg_len["lg_len"]
    lg_len["lg_mid"] = (lg_len["lg_start"] + lg_len["lg_end"]) / 2.0

    ax.set_xticks(lg_len["lg_mid"].to_numpy())
    ax.set_xticklabels(lg_len[a.lg_col].astype(str).to_list(), rotation=0)

    # Vertical boundaries between LGs
    for xline in lg_len["lg_end"].to_numpy()[:-1]:
        ax.axvline(xline, linewidth=0.6, alpha=0.25)

    # Labels / style
    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("Linkage Group")
    if a.title:
        ax.set_title(a.title)

    # Y limits with some headroom
    ymax = max(df["mlogp"].max(), a.threshold) if len(df) else a.threshold
    ax.set_ylim(0, ymax + 0.5)

    ax.grid(axis="y", linewidth=0.6, alpha=0.25)
    ax.set_axisbelow(True)

    # Legend (traits only) on the right
    if len(df_sig) > 0:
        ax.legend(title=a.trait_col, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, borderaxespad=0.0)

    fig.tight_layout()
    outpath = Path(a.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=a.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()

