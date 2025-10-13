#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-TE per-window enrichment (Poisson & NegBin) + correlations vs FST q-values.

Run guide

python te_enrichment_and_correlations.py \
  --in fst_windows_with_TEcounts.tsv \
  --outdir te_teclass_enrichment_out \
  --q_fst_colset q_poi_3vs4 q_poi_3vs5 q_poi_4vs5


Inputs (TSV; tab-separated) with exact columns (subset shown):
  CHROM  WIN_START  WIN_END  WIN_MID  WIN_LEN
  q_poi_3vs4  q_poi_3vs5  q_poi_4vs5
  TE_count  TE_B2_SINE_retrotransposon  TE_CACTA_TIR_transposon  ...
  TE_target_site_duplication
(23 TE/TE-class count columns listed by the user)

Outputs:
  OUTDIR/
    fst_windows_with_TEcounts__augmented.tsv     # original + p/q for each TE (Poisson+NB)
    te_enrichment_summary.tsv                    # one row per TE with lambda, alpha, dispersion diag
    te_qpoi_correlations.tsv                     # 69 rows (3 q_poi cols × 23 TE cols) with r/p/q and regression stats
    figures/                                     # regression PNGs organized per q_poi column
"""

import os
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import poisson, nbinom, pearsonr, linregress

# -------------------------
# Edit TE columns here only if needed (must match your file EXACTLY)
# -------------------------
TE_COLS = [
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

# Default FST q-value columns (can override via CLI)
Q_POI_COLS_DEFAULT = ["q_poi_3vs4", "q_poi_3vs5", "q_poi_4vs5"]

# -------------------------

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True, help="Input TSV (v9 table)")
    ap.add_argument("--outdir", default="te_teclass_enrichment_out", help="Output directory")
    ap.add_argument("--q_fst_colset", nargs="+", default=Q_POI_COLS_DEFAULT,
                    help="Which q_poi columns to use for correlations (default: 3vs4,3vs5,4vs5)")
    ap.add_argument("--min_alpha", type=float, default=1e-8, help="Minimum NB alpha to avoid degenerate NB")
    ap.add_argument("--plot_subset", type=int, default=None,
                    help="Optional: only plot first N TE columns (for quick runs)")
    return ap.parse_args()

def bh_qvalues(p):
    """Benjamini-Hochberg FDR for a 1D array-like (ignores NaN)."""
    p = np.asarray(p, dtype=float)
    n = np.sum(~np.isnan(p))
    if n == 0: return np.full_like(p, np.nan, dtype=float)
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(p) + 1, dtype=float)
    q = p * (n / ranks)
    # enforce monotonicity
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q_adj = np.empty_like(q)
    q_adj[order] = np.minimum(q_sorted, 1.0)
    q_adj[np.isnan(p)] = np.nan
    return q_adj

def safe_nb_params(mu_vec, y_vec, min_alpha=1e-8):
    """
    Estimate a global NB alpha using method-of-moments:
      Var(Y) ≈ mean((Y - mu)^2)  and  Var(Y) ≈ mean(mu) + alpha * mean(mu^2)
      => alpha ≈ max( (Var - mean(mu)) / mean(mu^2), min_alpha )
    Robust to windows with zero length (mu=0) and NaNs.
    """
    mu = np.asarray(mu_vec, dtype=float)
    y  = np.asarray(y_vec,  dtype=float)
    ok = np.isfinite(mu) & np.isfinite(y)
    if not np.any(ok):
        return max(min_alpha, 1e-8)
    mu_ok = mu[ok]; y_ok = y[ok]
    var_obs = np.mean((y_ok - mu_ok) ** 2)
    m1 = np.mean(mu_ok)
    m2 = np.mean(mu_ok ** 2)
    if m2 <= 0:
        return max(min_alpha, 1e-8)
    alpha = (var_obs - m1) / max(m2, 1e-12)
    if not np.isfinite(alpha):
        alpha = min_alpha
    return float(max(alpha, min_alpha))

def per_te_enrichment(df, te_col, min_alpha=1e-8):
    """
    For a given TE count column:
      - compute global lambda = sum(counts)/sum(length)
      - per-window Poisson p = P[Y >= y_i | mu_i=lambda*L_i] (upper tail)
      - estimate NB alpha by method-of-moments from mu_i and y_i
      - per-window NB p = P[Y >= y_i | NB(mu_i, alpha)]
      - BH q-values across windows (separately for Poisson and NB)
    Returns dict of new columns (Series) and a small summary dict.
    """
    y = pd.to_numeric(df[te_col], errors="coerce").fillna(0).astype(float).values
    L = pd.to_numeric(df["WIN_LEN"], errors="coerce").fillna((df["WIN_END"] - df["WIN_START"] + 1)).clip(lower=1).astype(float).values

    # Global rate
    total_counts = float(np.nansum(y))
    total_bases  = float(np.nansum(L))
    lam = total_counts / total_bases if total_bases > 0 else np.nan
    mu = lam * L if np.isfinite(lam) else np.full_like(L, np.nan)

    # Poisson p-values (upper-tail): p = 1 - CDF(y-1; mu)
    p_poi = np.full_like(y, np.nan, dtype=float)
    ok = np.isfinite(mu)
    if np.any(ok):
        p_poi[ok] = 1.0 - poisson.cdf(np.maximum(0, y[ok]) - 1, mu[ok])

    # NB parameters: r (size) = 1/alpha; p = r/(r+mu)
    alpha = safe_nb_params(mu, y, min_alpha=min_alpha)
    r = 1.0 / alpha if alpha > 0 else np.inf

    p_nb = np.full_like(y, np.nan, dtype=float)
    if np.isfinite(r) and np.any(ok):
        # SciPy nbinom: #successes r, probability p; mean = r*(1-p)/p
        # With mean=mu and alpha=1/r: p = r/(r+mu)
        p_param = r / (r + mu[ok])
        # guard: p must be in (0,1)
        p_param = np.clip(p_param, 1e-12, 1 - 1e-12)
        p_nb_vals = 1.0 - nbinom.cdf(np.maximum(0, y[ok]) - 1, r, p_param)
        p_nb[ok] = p_nb_vals

    # BH q-values across windows (per TE column)
    q_poi = bh_qvalues(p_poi)
    q_nb  = bh_qvalues(p_nb)

    # Wrangle column suffix-friendly name
    suf = te_col

    outcols = {
        f"p_poi_{suf}": p_poi,
        f"q_poi_{suf}": q_poi,
        f"p_nb_{suf}":  p_nb,
        f"q_nb_{suf}":  q_nb,
    }
    summary = {
        "te_col": te_col,
        "lambda_per_bp": lam,
        "alpha_nb": alpha,
        "total_counts": total_counts,
        "total_bases": total_bases
    }
    return outcols, summary

def make_reg_plot(x, y, r, p, slope, intercept, out_png, xlabel, ylabel, title):
    plt.figure(figsize=(5.2, 4.2))
    # scatter
    plt.scatter(x, y, s=6, alpha=0.5)
    # regression line
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ys = slope * xs + intercept
    plt.plot(xs, ys, linewidth=1.5)
    # labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # annotation: Pearson r and p-value
    txt = f"r = {r:.3f}\np = {p:.2e}"
    plt.gca().text(0.02, 0.98, txt, transform=plt.gca().transAxes,
                   ha="left", va="top", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.5", alpha=0.8))
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def main():
    args = parse_args()
    infile = args.infile
    outdir = args.outdir
    q_cols = args.q_fst_colset
    min_alpha = args.min_alpha
    plot_subset = args.plot_subset

    os.makedirs(outdir, exist_ok=True)
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    df = pd.read_csv(infile, sep="\t")
    # Basic presence checks
    need = ["WIN_START","WIN_END","WIN_LEN"] + q_cols
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit(f"Missing required columns: {miss}")

    # -----------------------------
    # 1) Per-TE per-window enrichment (Poisson & NB)
    # -----------------------------
    all_summ = []
    newcols = {}

    te_cols_iter = TE_COLS[:plot_subset] if plot_subset else TE_COLS
    for te in te_cols_iter:
        if te not in df.columns:
            raise SystemExit(f"TE column missing in input: {te}")
        outcols, summ = per_te_enrichment(df, te, min_alpha=min_alpha)
        newcols.update(outcols)
        all_summ.append(summ)

    # Append to original dataframe (preserve original columns first)
    for k, v in newcols.items():
        df[k] = v

    # Write augmented table
    out_aug = os.path.join(outdir, "fst_windows_with_TEcounts__augmented.tsv")
    df.to_csv(out_aug, sep="\t", index=False)

    # Per-TE summary table (lambda, alpha, totals)
    summ_df = pd.DataFrame(all_summ, columns=["te_col","lambda_per_bp","alpha_nb","total_counts","total_bases"])
    summ_df.to_csv(os.path.join(outdir, "te_enrichment_summary.tsv"), sep="\t", index=False)

    # -----------------------------
    # 2) Correlations & simple regressions: q_poi_* vs TE counts (all 3×23 = 69)
    # -----------------------------
    rows = []
    # For BH FDR we’ll compute per-pair and global later
    for qcol in q_cols:
        if qcol not in df.columns:
            raise SystemExit(f"Requested q_poi column missing: {qcol}")
        y = pd.to_numeric(df[qcol], errors="coerce")

        pair_figdir = os.path.join(figdir, qcol)
        os.makedirs(pair_figdir, exist_ok=True)

        for te in TE_COLS:
            x = pd.to_numeric(df[te], errors="coerce")
            # drop rows where either is NaN
            ok = x.notna() & y.notna()
            n = int(ok.sum())
            if n < 3:
                r = np.nan; p = np.nan; slope = np.nan; intercept = np.nan; p_slope = np.nan
            else:
                r, p = pearsonr(x[ok], y[ok])
                lr = linregress(x[ok].values, y[ok].values)
                slope, intercept, p_slope = lr.slope, lr.intercept, lr.pvalue

                # Plot
                out_png = os.path.join(pair_figdir, f"{qcol}__vs__{te}.png")
                make_reg_plot(x[ok].values, y[ok].values, r, p, slope, intercept,
                              out_png,
                              xlabel=f"{te} (count per window)",
                              ylabel=f"{qcol} (FST Poisson q)",
                              title=f"{qcol} vs {te}")

            rows.append({
                "qcol": qcol, "te_col": te, "n": n,
                "pearson_r": r, "pearson_p": p,
                "slope": slope, "intercept": intercept, "p_slope": p_slope,
                # placeholders for q-values (filled below)
                "q_r_pair": np.nan, "q_r_global": np.nan
            })

    corr_df = pd.DataFrame(rows)

    # BH per pair (within each qcol across its 23 tests)
    for qcol in q_cols:
        mask = corr_df["qcol"] == qcol
        corr_df.loc[mask, "q_r_pair"] = bh_qvalues(corr_df.loc[mask, "pearson_p"].values)

    # BH global across all 69
    corr_df["q_r_global"] = bh_qvalues(corr_df["pearson_p"].values)

    out_corr = os.path.join(outdir, "te_qpoi_correlations.tsv")
    corr_df.to_csv(out_corr, sep="\t", index=False)

    print("[✓] Done")
    print("Wrote:")
    print("  ", out_aug)
    print("  ", os.path.join(outdir, "te_enrichment_summary.tsv"))
    print("  ", out_corr)
    print("  ", figdir)

if __name__ == "__main__":
    main()

