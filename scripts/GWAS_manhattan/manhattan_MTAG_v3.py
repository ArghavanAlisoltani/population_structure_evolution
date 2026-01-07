#!/usr/bin/env python3
"""
Manhattan plot (LG-ordered) with:
- Background (mlogp < threshold): alternating grey by LG parity
- Significant (mlogp >= threshold): colored by TRAIT
- Point SHAPES by Symbol column:
    Single trait -> circle
    MTAG         -> triangle

Legend:
  Organized by Symbol first, then Trait ("Symbol | Trait") for significant points.

NEW:
  Control point sizes per SYMBOL separately (and optionally per significance):
    --size_single_bg, --size_mtag_bg
    --size_single_sig, --size_mtag_sig

Example:
  python3 manhattan_lg_symbol_sizes.py \
    --infile sumstat.tsv \
    --out manhattan.png \
    --threshold 6 \
    --p_scale raw \
    --symbol_col Simbol \
    --size_single_bg 8 \
    --size_mtag_bg 10 \
    --size_single_sig 60 \
    --size_mtag_sig 110
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sep", default="\t")

    ap.add_argument("--p_col", default="P_value")
    ap.add_argument("--p_scale", choices=["raw", "mlog10"], default="raw")

    ap.add_argument("--trait_col", default="TRAIT")
    ap.add_argument("--symbol_col", default="Simbol")

    ap.add_argument("--lg_col", default="LG_group")
    ap.add_argument("--scaf_col", default="Scaffold")
    ap.add_argument("--scaf_order_col", default="Scaffold_order_per_LG")
    ap.add_argument("--pos_col", default="POS")
    ap.add_argument("--len_col", default="Length")

    ap.add_argument("--threshold", type=float, default=6.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--fig_w", type=float, default=18.0)
    ap.add_argument("--fig_h", type=float, default=5.5)

    # Per-symbol sizes (background + significant)
    ap.add_argument("--size_single_bg", type=float, default=6.0, help="Circle size for background (Single trait).")
    ap.add_argument("--size_mtag_bg", type=float, default=7.0, help="Triangle size for background (MTAG).")
    ap.add_argument("--size_single_sig", type=float, default=55.0, help="Circle size for significant (Single trait).")
    ap.add_argument("--size_mtag_sig", type=float, default=85.0, help="Triangle size for significant (MTAG).")

    ap.add_argument("--max_points_bg", type=int, default=0,
                    help="Downsample background points to this many (0 = no downsampling).")
    ap.add_argument("--legend_max_entries", type=int, default=120,
                    help="Max legend entries (Symbol|Trait combos). 0 = no limit.")
    return ap.parse_args()


def lg_to_num(x) -> int:
    if pd.isna(x):
        return 10**9
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else 10**9


def distinct_colors(n: int):
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        return [cmap(i) for i in range(n)]
    cmap = plt.get_cmap("hsv")
    return [cmap(i / n) for i in range(n)]


def compute_mlogp(pvals: np.ndarray, scale: str) -> np.ndarray:
    if scale == "mlog10":
        return pvals.astype(float)
    p_floor = 1e-300
    p = pvals.astype(float)
    p = np.where(np.isnan(p), np.nan, p)
    p = np.where(p <= 0, p_floor, p)
    return -np.log10(p)


def normalize_symbol(s) -> str:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return "NA"
    t = str(s).strip().lower()
    if t in {"single", "singletrait", "single trait", "single-trait"}:
        return "Single trait"
    if t == "mtag":
        return "MTAG"
    return str(s).strip()


def main():
    a = parse_args()

    df = pd.read_csv(a.infile, sep=a.sep, dtype=str, low_memory=False)

    required = [a.p_col, a.trait_col, a.symbol_col, a.lg_col, a.scaf_col, a.pos_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}\nFound columns: {list(df.columns)}")

    # numeric conversions
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

    df = df.dropna(subset=[a.p_col, a.pos_col, a.lg_col, a.scaf_col, a.symbol_col])

    df["mlogp"] = compute_mlogp(df[a.p_col].to_numpy(), a.p_scale)
    df = df.dropna(subset=["mlogp"])

    df["lg_num"] = df[a.lg_col].map(lg_to_num).astype(int)
    df["symbol_norm"] = df[a.symbol_col].map(normalize_symbol)

    # marker + size maps
    marker_map = {"Single trait": "o", "MTAG": "^", "NA": "o"}
    size_bg_map = {"Single trait": a.size_single_bg, "MTAG": a.size_mtag_bg, "NA": a.size_single_bg}
    size_sig_map = {"Single trait": a.size_single_sig, "MTAG": a.size_mtag_sig, "NA": a.size_single_sig}
    symbol_order = ["Single trait", "MTAG"]

    # estimate Length if missing
    if df[a.len_col].isna().all():
        est_len = (
            df.groupby(["lg_num", a.scaf_col], as_index=False)[a.pos_col]
              .max()
              .rename(columns={a.pos_col: a.len_col})
        )
        df = df.merge(est_len, on=["lg_num", a.scaf_col], how="left")

    scaf_cols = ["lg_num", a.lg_col, a.scaf_col, a.scaf_order_col, a.len_col]
    scaf_tbl = df[scaf_cols].drop_duplicates()

    if scaf_tbl[a.scaf_order_col].isna().all():
        scaf_tbl = scaf_tbl.sort_values(["lg_num", a.scaf_col]).copy()
        scaf_tbl[a.scaf_order_col] = scaf_tbl.groupby("lg_num").cumcount().astype(int) + 1

    if scaf_tbl[a.len_col].isna().any():
        est_len2 = (
            df.groupby(["lg_num", a.scaf_col], as_index=False)[a.pos_col]
              .max()
              .rename(columns={a.pos_col: a.len_col})
        )
        scaf_tbl = scaf_tbl.drop(columns=[a.len_col]).merge(est_len2, on=["lg_num", a.scaf_col], how="left")

    scaf_tbl[a.len_col] = scaf_tbl[a.len_col].fillna(0).astype(float)
    scaf_tbl = scaf_tbl.sort_values(["lg_num", a.scaf_order_col, a.scaf_col]).copy()

    scaf_tbl["offset_in_lg"] = scaf_tbl.groupby("lg_num")[a.len_col].cumsum() - scaf_tbl[a.len_col]

    lg_len = (scaf_tbl.groupby(["lg_num", a.lg_col], as_index=False)[a.len_col]
              .sum()
              .rename(columns={a.len_col: "lg_len"})
              .sort_values("lg_num")
              .copy())
    lg_len["lg_global_offset"] = lg_len["lg_len"].cumsum() - lg_len["lg_len"]

    df = (
        df.merge(
            scaf_tbl[["lg_num", a.scaf_col, a.scaf_order_col, "offset_in_lg"]],
            on=["lg_num", a.scaf_col, a.scaf_order_col],
            how="left"
        )
        .merge(
            lg_len[["lg_num", "lg_global_offset", a.lg_col, "lg_len"]],
            on=["lg_num", a.lg_col],
            how="left"
        )
    )

    df["x"] = df["lg_global_offset"] + df["offset_in_lg"] + df[a.pos_col].astype(float)
    df = df.sort_values(["lg_num", a.scaf_order_col, a.pos_col])

    sig_mask = df["mlogp"] >= a.threshold
    df_sig = df.loc[sig_mask].copy()
    df_bg = df.loc[~sig_mask].copy()

    if a.max_points_bg and a.max_points_bg > 0 and len(df_bg) > a.max_points_bg:
        df_bg = df_bg.sample(n=a.max_points_bg, random_state=1).sort_values("x")

    # trait colors (from significant points)
    if len(df_sig) > 0:
        df_sig[a.trait_col] = df_sig[a.trait_col].fillna("NA").astype(str)
        traits = sorted(df_sig[a.trait_col].unique().tolist())
        colors = distinct_colors(len(traits))
        trait2color = dict(zip(traits, colors))
    else:
        traits, trait2color = [], {}

    fig, ax = plt.subplots(figsize=(a.fig_w, a.fig_h))

    grey75 = (0.75, 0.75, 0.75, 1.0)
    grey88 = (0.88, 0.88, 0.88, 1.0)

    # background with per-symbol sizes + markers
    for sym in sorted(df_bg["symbol_norm"].unique().tolist()):
        m = marker_map.get(sym, "o")
        s_bg = size_bg_map.get(sym, a.size_single_bg)
        d = df_bg[df_bg["symbol_norm"] == sym]
        even = (d["lg_num"] % 2 == 0)
        ax.scatter(d.loc[~even, "x"], d.loc[~even, "mlogp"],
                   s=s_bg, c=[grey88], marker=m, linewidths=0, alpha=1.0, rasterized=True)
        ax.scatter(d.loc[even, "x"], d.loc[even, "mlogp"],
                   s=s_bg, c=[grey75], marker=m, linewidths=0, alpha=1.0, rasterized=True)

    # significant points: plot by (symbol, trait) in desired order
    if len(df_sig) > 0:
        df_sig["symbol_norm"] = df_sig["symbol_norm"].astype(str)
        present_symbols = [s for s in symbol_order if s in set(df_sig["symbol_norm"])]
        others = sorted(set(df_sig["symbol_norm"]) - set(present_symbols))
        present_symbols += others

        for sym in present_symbols:
            msym = marker_map.get(sym, "o")
            ssig = size_sig_map.get(sym, a.size_single_sig)
            d_sym = df_sig[df_sig["symbol_norm"] == sym]
            for tr in sorted(d_sym[a.trait_col].unique().tolist()):
                d = d_sym[d_sym[a.trait_col] == tr]
                ax.scatter(d["x"], d["mlogp"],
                           s=ssig, c=[trait2color.get(tr, (0, 0, 0, 1))],
                           marker=msym, linewidths=0, alpha=0.95, rasterized=True)

    ax.axhline(a.threshold, linestyle="--", linewidth=1.2, alpha=0.9)

    # LG ticks
    lg_len["lg_start"] = lg_len["lg_global_offset"]
    lg_len["lg_end"] = lg_len["lg_global_offset"] + lg_len["lg_len"]
    lg_len["lg_mid"] = (lg_len["lg_start"] + lg_len["lg_end"]) / 2.0
    ax.set_xticks(lg_len["lg_mid"].to_numpy())
    ax.set_xticklabels(lg_len[a.lg_col].astype(str).to_list(), rotation=0)

    # LG separators
    for xline in lg_len["lg_end"].to_numpy()[:-1]:
        ax.axvline(xline, linewidth=0.6, alpha=0.25)

    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("Linkage Group")
    if a.title:
        ax.set_title(a.title)

    ymax = max(df["mlogp"].max(), a.threshold) if len(df) else a.threshold
    ax.set_ylim(0, ymax + 0.5)
    ax.grid(axis="y", linewidth=0.6, alpha=0.25)
    ax.set_axisbelow(True)

    # Legend: Symbol first, then Trait (significant combos), with marker + color + size hint
    if len(df_sig) > 0:
        present_symbols = [s for s in symbol_order if s in set(df_sig["symbol_norm"])]
        others = sorted(set(df_sig["symbol_norm"]) - set(present_symbols))
        present_symbols += others

        handles, labels = [], []
        for sym in present_symbols:
            msym = marker_map.get(sym, "o")
            # make legend marker reflect the significant size (scaled to markersize units)
            # markersize is in points; scatter 's' is area. This is a reasonable visual mapping.
            legend_ms = max(5, np.sqrt(size_sig_map.get(sym, a.size_single_sig)))
            d_sym = df_sig[df_sig["symbol_norm"] == sym]
            for tr in sorted(d_sym[a.trait_col].unique().tolist()):
                col = trait2color.get(tr, (0, 0, 0, 1))
                handles.append(
                    Line2D([0], [0], marker=msym, linestyle="None",
                           markerfacecolor=col, markeredgecolor=col,
                           markersize=legend_ms)
                )
                labels.append(f"{sym} | {tr}")

        if a.legend_max_entries and a.legend_max_entries > 0:
            handles = handles[:a.legend_max_entries]
            labels = labels[:a.legend_max_entries]

        ax.legend(handles, labels,
                  title=f"{a.symbol_col} | {a.trait_col}",
                  loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, borderaxespad=0.0)

    fig.tight_layout()
    outpath = Path(a.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=a.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()

