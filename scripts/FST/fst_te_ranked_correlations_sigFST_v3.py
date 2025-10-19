#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter to Poisson-significant FST windows per comparison,
then compute ranked (Spearman) correlations:
  A) FST q  vs TE/TE-class Poisson q
  B) FST counts (C0.25_*) vs TE/TE-class counts

Robust to constant inputs and degenerate cases. Plots omit the fit line when invalid.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

COMPS = ["3vs4", "3vs5", "4vs5"]

def find_te_q_columns(df, prefix="q_poi_TE_"):
    return [c for c in df.columns if c.startswith(prefix)]

def find_te_count_columns(df):
    cols = []
    for c in df.columns:
        if c == "TE_count":
            cols.append(c)
        elif c.startswith("TE_") and not c.startswith("TEdens_") and c != "TE_density":
            cols.append(c)
    # de-duplicate while preserving order
    seen, out = set(), []
    for c in cols:
        if c not in seen:
            out.append(c); seen.add(c)
    return out

def to_num(arr_like):
    return pd.to_numeric(arr_like, errors="coerce").to_numpy()

def valid_xy(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return m, n, "too_few_points"
    ux = np.unique(x[m])
    uy = np.unique(y[m])
    if ux.size < 2:
        return m, n, "constant_x"
    if uy.size < 2:
        return m, n, "constant_y"
    return m, n, None

def spearman_safe(x, y):
    m, n, why = valid_xy(x, y)
    if why is not None:
        return np.nan, np.nan, n, why
    rho, p = spearmanr(x[m], y[m])
    return float(rho), float(p), n, None

def ols_line_closed_form(x, y):
    """Return slope, intercept using cov/var on finite values; None if invalid."""
    m, n, why = valid_xy(x, y)
    if why is not None:
        return None
    xx, yy = x[m], y[m]
    vx = np.var(xx)
    if vx == 0 or not np.isfinite(vx):
        return None
    cov = np.cov(xx, yy, ddof=0)[0, 1]
    b1 = cov / vx
    b0 = float(np.mean(yy) - b1 * np.mean(xx))
    if not (np.isfinite(b1) and np.isfinite(b0)):
        return None
    return b1, b0, float(np.min(xx)), float(np.max(xx))

def scatter_with_optional_line(x, y, xlab, ylab, title, annot, out_png, xlim=None, ylim=None, line_params=None):
    plt.figure(figsize=(7.5, 6.0), dpi=150)
    ax = plt.gca()
    ax.scatter(x, y, alpha=0.35, s=16)  # no linewidths/edgecolors to avoid alias issues
    if line_params is not None:
        b1, b0, xmin, xmax = line_params
        xs = np.array([xmin, xmax], dtype=float)
        ys = b1 * xs + b0
        ax.plot(xs, ys, linewidth=1.5)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab); ax.set_title(title)
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    ax.text(0.03, 0.97, annot, transform=ax.transAxes,
            va="top", ha="left",
            bbox=dict(fc="white", ec="0.6", alpha=0.9))
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def compute_q_vs_q(sig_df, comp, te_q_cols, outdir, fst_q_col):
    rows = []
    plotdir = Path(outdir, "qq_plots"); plotdir.mkdir(parents=True, exist_ok=True)

    fq = to_num(sig_df[fst_q_col])
    for tecol in te_q_cols:
        tq = to_num(sig_df[tecol])
        rho, p, n, note = spearman_safe(fq, tq)
        rows.append(dict(fst_q=fst_q_col, te_q=tecol, n=n, spearman_rho=rho, p_value=p, note=note))
        # plot
        annot = f"ρ = {np.nan if np.isnan(rho) else round(rho,3)}\n" \
                f"p = {np.nan if np.isnan(p) else f'{p:.2e}'}\n" \
                f"n = {n}" + (f"\n{note}" if note else "")
        line_params = None if note else ols_line_closed_form(tq, fq)
        scatter_with_optional_line(
            x=tq, y=fq,
            xlab=tecol, ylab=fst_q_col,
            title=f"{comp}: {fst_q_col} vs {tecol}",
            annot=annot, out_png=plotdir / f"{comp}__{fst_q_col}__vs__{tecol}.png",
            xlim=(0,1), ylim=(0,1),
            line_params=line_params
        )

    pd.DataFrame(rows).sort_values(["note","p_value"], na_position="last") \
        .to_csv(Path(outdir, f"qq_correlations_{comp}.tsv"), sep="\t", index=False)

def compute_counts(sig_df, comp, te_count_cols, outdir, fst_count_col):
    rows = []
    plotdir = Path(outdir, "count_plots"); plotdir.mkdir(parents=True, exist_ok=True)

    fy = to_num(sig_df[fst_count_col])
    for tecol in te_count_cols:
        tx = to_num(sig_df[tecol])
        rho, p, n, note = spearman_safe(fy, tx)
        rows.append(dict(fst_count=fst_count_col, te_count=tecol, n=n, spearman_rho=rho, p_value=p, note=note))
        annot = f"ρ = {np.nan if np.isnan(rho) else round(rho,3)}\n" \
                f"p = {np.nan if np.isnan(p) else f'{p:.2e}'}\n" \
                f"n = {n}" + (f"\n{note}" if note else "")
        # Axis limits from finite data
        m = np.isfinite(tx) & np.isfinite(fy)
        xlim = (float(np.min(tx[m])) if m.any() else None, float(np.max(tx[m])) if m.any() else None)
        ylim = (float(np.min(fy[m])) if m.any() else None, float(np.max(fy[m])) if m.any() else None)
        line_params = None if note else ols_line_closed_form(tx, fy)
        scatter_with_optional_line(
            x=tx, y=fy,
            xlab=tecol, ylab=fst_count_col,
            title=f"{comp}: {fst_count_col} vs {tecol}",
            annot=annot, out_png=plotdir / f"{comp}__{fst_count_col}__vs__{tecol}.png",
            xlim=None if None in xlim else xlim,
            ylim=None if None in ylim else ylim,
            line_params=line_params
        )
    pd.DataFrame(rows).sort_values(["note","p_value"], na_position="last") \
        .to_csv(Path(outdir, f"count_correlations_{comp}.tsv"), sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser(description="Rank correlations between FST and TE metrics using only FST-poisson-significant windows.")
    ap.add_argument("-i", "--input", required=True, help="Input TSV (fst_windows_with_TEcounts__augmented.tsv)")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory")
    ap.add_argument("--fst_q_threshold", type=float, default=0.05,
                    help="Threshold on FST Poisson q to define significant windows (default: 0.05)")
    ap.add_argument("--fst_count_choice", choices=["C0", "C0.25"], default="C0.25",
                    help="Which FST count column family to use for the count correlations (default: C0.25)")
    ap.add_argument("--min_rows", type=int, default=10,
                    help="Minimum number of significant windows required to compute correlations (default: 10)")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    te_q_cols = find_te_q_columns(df, prefix="q_poi_TE_")
    te_count_cols = find_te_count_columns(df)

    # sanity
    for comp in COMPS:
        if f"q_poi_{comp}" not in df.columns:
            raise SystemExit(f"Missing column q_poi_{comp}")
        if f"{args.fst_count_choice}_{comp}" not in df.columns:
            raise SystemExit(f"Missing column {args.fst_count_choice}_{comp}")

    # per-comparison
    for comp in COMPS:
        subdir = Path(outdir, comp); subdir.mkdir(parents=True, exist_ok=True)
        fst_q_col = f"q_poi_{comp}"
        fst_count_col = f"{args.fst_count_choice}_{comp}"

        sig = df.loc[pd.to_numeric(df[fst_q_col], errors="coerce") <= args.fst_q_threshold].copy()
        sig.to_csv(Path(subdir, f"sig_windows_{comp}.tsv"), sep="\t", index=False)

        if len(sig) < args.min_rows:
            # emit empty tables for completeness
            pd.DataFrame(columns=["fst_q","te_q","n","spearman_rho","p_value","note"]).to_csv(
                Path(subdir, f"qq_correlations_{comp}.tsv"), sep="\t", index=False)
            pd.DataFrame(columns=["fst_count","te_count","n","spearman_rho","p_value","note"]).to_csv(
                Path(subdir, f"count_correlations_{comp}.tsv"), sep="\t", index=False)
            print(f"[WARN] {comp}: only {len(sig)} significant rows; skipping correlations.")
            continue

        compute_q_vs_q(sig, comp, te_q_cols, subdir, fst_q_col)
        compute_counts(sig, comp, te_count_cols, subdir, fst_count_col)

    with open(Path(outdir, "README.txt"), "w") as fh:
        fh.write(
f"""Comparisons: {', '.join(COMPS)}
FST-q threshold: {args.fst_q_threshold}
FST counts used: {args.fst_count_choice}_*

Per comparison directory:
  sig_windows_<comp>.tsv
  qq_correlations_<comp>.tsv  (Spearman + notes: constant_x/constant_y/too_few_points)
  qq_plots/*.png              (fit line omitted if invalid)
  count_correlations_<comp>.tsv
  count_plots/*.png

Notes:
- TE q columns matched by '^q_poi_TE_'.
- TE count columns include TE_count and all 'TE_*' except density columns.
- All numeric inputs are coerced with errors='coerce'.
"""
        )

if __name__ == "__main__":
    main()

