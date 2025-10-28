#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Correlations/regressions: FST q_poi_* vs TE/TE-class q_poi_* (23 classes)
+ Manhattan plots for each TE/TE-class q_poi column.

This version is robust to:
 - mixed-type columns (reads as strings, then coerces needed cols)
 - constant vectors (skips regression & correlation gracefully)

Usage:
python te_qvalue_cor_and_manhattan_v2.py \
  --in fst_windows_with_TEcounts__augmented.tsv \
  --outdir te_qval_cor_manhattan_out \
  --fst_qcols q_poi_3vs4 q_poi_3vs5 q_poi_4vs5
"""

import os
import re
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress

# ---- Exact TE count column names (23), as in your table ----
TE_COUNT_COLS = [
    "TE_count",
    "TE_B2_SINE_retrotransposon",
    "TE_CACTA_TIR_transposon",
    "TE_CR1_LINE_retrotransposon",
    "TE_Copia_LTR_retrotransposon",
    "TE_Gypsy_LTR_retrotransposon",
    "TE_ID_SINE_retrotransposon",
    "TE_Jockey_LINE_retrotransposon",
    "TE_L1_LINE_retrotransposon",
    "TE_L2_LINE_retrotransposon",
    "TE_LINE_element",
    "TE_LTR_retrotransposon",
    "TE_Mutator_TIR_transposon",
    "TE_PIF_Harbinger_TIR_transposon",
    "TE_Penelope_retrotransposon",
    "TE_RTE_LINE_retrotransposon",
    "TE_Tad1_LINE_retrotransposon",
    "TE_Tc1_Mariner_TIR_transposon",
    "TE_hAT_TIR_transposon",
    "TE_helitron",
    "TE_long_terminal_repeat",
    "TE_repeat_fragment",
    "TE_repeat_region",
    "TE_target_site_duplication",
]

# FST q-value columns
FST_Q_COLS_DEFAULT = ["q_poi_3vs4", "q_poi_3vs5", "q_poi_4vs5"]

# ---- Utils ----
def bh_qvalues(p):
    """Benjamini–Hochberg FDR for a 1D array-like (ignores NaN)."""
    p = np.asarray(p, dtype=float)
    n = np.sum(~np.isnan(p))
    if n == 0:
        return np.full_like(p, np.nan, dtype=float)
    # rank NaNs to the end
    p_nonan = np.where(np.isnan(p), np.inf, p)
    order = np.argsort(p_nonan)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(p) + 1, dtype=float)
    q = p * (n / ranks)
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q_adj = np.empty_like(q)
    q_adj[order] = np.minimum(q_sorted, 1.0)
    q_adj[np.isnan(p)] = np.nan
    return q_adj

def parse_scaffold_key(chrom):
    """
    Turn CHROM like 'scaffold_1a', 'scaffold_1b', 'scaffold_2', ... into a sort key.
    Orders by numeric part first, then optional letter. Falls back to lexical.
    """
    if not isinstance(chrom, str):
        chrom = str(chrom)
    m = re.search(r'(\d+)\s*([A-Za-z]?)$', chrom)
    if m:
        num = int(m.group(1))
        letter = m.group(2).lower() if m.group(2) else ''
        letter_rank = ord(letter) - ord('a') + 1 if letter else 0
        return (num, letter_rank, chrom)
    m2 = re.search(r'scaffold[_-]?(\d+)\s*([A-Za-z]?)', chrom, flags=re.I)
    if m2:
        num = int(m2.group(1))
        letter = m2.group(2).lower() if m2.group(2) else ''
        letter_rank = ord(letter) - ord('a') + 1 if letter else 0
        return (num, letter_rank, chrom)
    return (10**9, 10**9, chrom)

def build_cumulative_pos(df, chrom_col="CHROM", pos_col="WIN_MID", end_col="WIN_END"):
    """
    Compute cumulative genomic position for plotting across scaffolds in numeric order.
    offset(scaff) = sum of max(end) of all previous scaffolds in order.
    Returns a Series of cumulative positions and a DataFrame of scaffold ticks.
    """
    scafs = (
        df[[chrom_col, end_col]]
        .groupby(chrom_col, as_index=False)
        .max()
        .assign(_key=lambda d: d[chrom_col].map(parse_scaffold_key))
        .sort_values("_key")
        .drop(columns="_key")
    )
    offsets = {}
    cum = 0.0
    ticks = []
    for _, row in scafs.iterrows():
        chrom = row[chrom_col]
        length = float(row[end_col])
        offsets[chrom] = cum
        ticks.append((chrom, cum + length/2.0))
        cum += length
    cumpos = df[pos_col].astype(float).values + df[chrom_col].map(offsets).astype(float).values
    ticks_df = pd.DataFrame(ticks, columns=[chrom_col, "tick_mid"])
    return pd.Series(cumpos, index=df.index), ticks_df

def make_reg_plot(x, y, r, p, slope, intercept, out_png, xlabel, ylabel, title, draw_line=True):
    plt.figure(figsize=(5.2, 4.2))
    plt.scatter(x, y, s=6, alpha=0.6)
    if draw_line and np.isfinite(slope) and np.isfinite(intercept):
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ys = slope * xs + intercept
        plt.plot(xs, ys, linewidth=1.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    txt = f"r = {r:.3f}" if np.isfinite(r) else "r = NA"
    txt += f"\np = {p:.2e}" if np.isfinite(p) else "\np = NA"
    plt.gca().text(0.02, 0.98, txt, transform=plt.gca().transAxes,
                   ha="left", va="top",
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.5", alpha=0.85))
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

def manhattan_plot(df, chrom_col, cumpos_col, y_col, out_png,
                   alpha_sig=1.0, alpha_nonsig=0.25, qcut=0.05,
                   colors=("#888888", "#1f77b4")):
    """
    Per-TE Manhattan: y = -log10(q), x = cumulative position, colored by alternating scaffold.
    Alpha encodes significance (q<qcut → alpha_sig, else alpha_nonsig).
    """
    d = df.copy()
    if y_col not in d.columns:
        raise ValueError(f"missing column: {y_col}")
    q = pd.to_numeric(d[y_col], errors="coerce")
    y = -np.log10(q.clip(lower=1e-300))
    d["_y"] = y
    d["_is_sig"] = (q < qcut).astype(int)

    scafs = d[[chrom_col, cumpos_col]].groupby(chrom_col, as_index=False).agg({"cumpos": "median"})
    scafs = scafs.sort_values("cumpos").reset_index(drop=True)
    scafs["parity"] = scafs.index % 2
    parity_map = dict(zip(scafs[chrom_col], scafs["parity"]))
    d["_parity"] = d[chrom_col].map(parity_map).fillna(0).astype(int)

    plt.figure(figsize=(11, 3.6))
    for par in (0, 1):
        sub = d[d["_parity"] == par]
        sub_sig = sub[sub["_is_sig"] == 1]
        sub_non = sub[sub["_is_sig"] == 0]
        if len(sub_non):
            plt.scatter(sub_non[cumpos_col], sub_non["_y"], s=8, alpha=alpha_nonsig,
                        c=colors[par], edgecolors="none")
        if len(sub_sig):
            plt.scatter(sub_sig[cumpos_col], sub_sig["_y"], s=10, alpha=alpha_sig,
                        c=colors[par], edgecolors="none")

    plt.ylabel(r"$-\log_{10}(q)$")
    plt.xlabel("Genomic position (concatenated scaffolds)")
    plt.title(y_col)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close()

def parse_args():
    """Parse CLI arguments for correlation and Manhattan plotting utilities."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True,
                    help="fst_windows_with_TEcounts__augmented.tsv")
    ap.add_argument("--outdir", default="te_qval_cor_manhattan_out",
                    help="Output directory")
    ap.add_argument("--fst_qcols", nargs="+", default=FST_Q_COLS_DEFAULT,
                    help="FST q-value columns to compare (default: q_poi_3vs4 q_poi_3vs5 q_poi_4vs5)")
    ap.add_argument("--qcut", type=float, default=0.05,
                    help="q-value threshold for significance in Manhattan")
    ap.add_argument("--only_first_n_te", type=int, default=None,
                    help="Optional: limit number of TE classes when testing/plotting")
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    INFILE = args.infile
    OUTDIR = args.outdir
    FST_Q_COLS = args.fst_qcols
    Q_CUT = args.qcut
    N_TE = args.only_first_n_te

    os.makedirs(OUTDIR, exist_ok=True)
    fig_corr_root = os.path.join(OUTDIR, "figures", "corr")
    fig_manh_dir = os.path.join(OUTDIR, "figures", "manhattan")
    os.makedirs(fig_corr_root, exist_ok=True)
    os.makedirs(fig_manh_dir, exist_ok=True)

    # --- Robust read: everything as string, then we'll coerce what we need ---
    df = pd.read_csv(INFILE, sep="\t", dtype=str, low_memory=False)

    # Ensure essential columns exist
    need = ["CHROM", "WIN_START", "WIN_END", "WIN_MID", "WIN_LEN"] + FST_Q_COLS
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    # Create numeric versions we need
    for col in ["WIN_START", "WIN_END", "WIN_MID", "WIN_LEN"] + FST_Q_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # TE q_poi columns
    TE_Q_COLS = [f"q_poi_{c}" for c in TE_COUNT_COLS]
    if N_TE is not None:
        TE_Q_COLS = TE_Q_COLS[:N_TE]
    missing_teq = [c for c in TE_Q_COLS if c not in df.columns]
    if missing_teq:
        raise SystemExit(f"Missing TE q-p columns in input: {missing_teq}")
    for col in TE_Q_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Build cumulative x for manhattan (use numeric WIN_MID/END we just coerced)
    df["cumpos"], ticks_df = build_cumulative_pos(df, chrom_col="CHROM", pos_col="WIN_MID", end_col="WIN_END")

    # -----------------------------
    # 1) Correlations & regressions: FST q_poi_* vs TE q_poi_* (3 x 23)
    # -----------------------------
    rows = []
    for fst_q in FST_Q_COLS:
        y_all = df[fst_q]
        q_dir = os.path.join(fig_corr_root, fst_q)
        os.makedirs(q_dir, exist_ok=True)

        for te_q in TE_Q_COLS:
            x_all = df[te_q]
            ok = x_all.notna() & y_all.notna()
            n = int(ok.sum())
            r = p = slope = intercept = p_slope = np.nan
            note = ""

            if n >= 3:
                x = x_all[ok].values.astype(float)
                y = y_all[ok].values.astype(float)
                # handle constant arrays
                x_var = float(np.nanvar(x))
                y_var = float(np.nanvar(y))
                if x_var == 0.0 or y_var == 0.0:
                    note = "constant_x" if x_var == 0.0 else "constant_y"
                    draw_line = False
                else:
                    try:
                        r, p = pearsonr(x, y)
                    except Exception:
                        r = np.nan; p = np.nan
                    try:
                        lr = linregress(x, y)
                        slope, intercept, p_slope = lr.slope, lr.intercept, lr.pvalue
                        draw_line = True
                    except Exception:
                        slope = intercept = p_slope = np.nan
                        draw_line = False

                # Plot (show points; draw line only if estimated)
                out_png = os.path.join(q_dir, f"{fst_q}__vs__{te_q}.png")
                make_reg_plot(
                    x, y, r, p, slope, intercept, out_png,
                    xlabel=f"{te_q}",
                    ylabel=f"{fst_q}",
                    title=f"{fst_q} vs {te_q}",
                    draw_line=draw_line
                )
            else:
                note = "n<3"

            rows.append({
                "fst_q": fst_q,
                "te_q": te_q,
                "n": n,
                "pearson_r": r,
                "pearson_p": p,
                "slope": slope,
                "intercept": intercept,
                "p_slope": p_slope,
                "note": note,
                "q_r_pair": np.nan,
                "q_r_global": np.nan
            })

    corr_df = pd.DataFrame(rows)

    # BH per FST q column (within its 23 tests)
    for fst_q in FST_Q_COLS:
        mask = corr_df["fst_q"] == fst_q
        corr_df.loc[mask, "q_r_pair"] = bh_qvalues(corr_df.loc[mask, "pearson_p"].values)

    # BH across all 3×23 tests
    corr_df["q_r_global"] = bh_qvalues(corr_df["pearson_p"].values)

    out_corr = os.path.join(OUTDIR, "te_fst_q_correlations.tsv")
    corr_df.to_csv(out_corr, sep="\t", index=False)

    # -----------------------------
    # 2) Manhattan plots per TE q_poi column
    # -----------------------------
    for te_q in TE_Q_COLS:
        out_png = os.path.join(fig_manh_dir, f"manhattan_{te_q}.png")
        dplot = df[["CHROM", "cumpos", te_q]].copy()
        manhattan_plot(
            dplot.rename(columns={te_q: te_q}),
            chrom_col="CHROM",
            cumpos_col="cumpos",
            y_col=te_q,
            out_png=out_png,
            alpha_sig=1.0,
            alpha_nonsig=0.25,
            qcut=Q_CUT,
            colors=("#808080", "#1f77b4")
        )

    print("[✓] Done.")
    print("Wrote:")
    print("  ", out_corr)
    print("  ", os.path.join(OUTDIR, "figures"))

