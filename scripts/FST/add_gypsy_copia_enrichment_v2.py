#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gypsy vs Copia enrichment tied to significant FST windows
(Exact column names from your v9 table)

Inputs (TSV; tab-separated) must contain at least:
  CHROM, WIN_START, WIN_END, WIN_MID, WIN_LEN,
  q_poi_3vs4, q_poi_3vs5, q_poi_4vs5,
  TE_Gypsy_LTR_retrotransposon, TE_Copia_LTR_retrotransposon,
  TEdens_Gypsy_LTR_retrotransposon, TEdens_Copia_LTR_retrotransposon

Outputs:
  - fst_windows_with_TEcounts__gypsy_copia.tsv  (adds Gypsy/Copia counts/dens + is_sig flags)
  - gypsy_copia_fisher_presence.tsv             (2x2 Fisher tests per pair and TE)
  - gypsy_copia_glm_interaction.tsv             (Poisson & NB interaction results)
  - figures/*.png                                (optional quick plots)
"""

import os
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
import statsmodels.api as sm
from statsmodels.genmod.families import Poisson, NegativeBinomial
from statsmodels.genmod.families.links import log as LogLink

PAIRS = ["3vs4","3vs5","4vs5"]

def parse_args():
    """Parse command-line options for the enrichment workflow."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True,
                    help="Path to v9 output TSV (fst_windows_with_TEcounts.tsv)")
    ap.add_argument("--outdir", default="gypsy_copia_addon",
                    help="Output directory")
    ap.add_argument("--qcut", type=float, default=0.05,
                    help="q-value cutoff for FST Poisson significance")
    ap.add_argument("--make_plots", action="store_true",
                    help="Emit small summary plots per pair")
    return ap.parse_args()

def check_required_columns(df):
    """Ensure the v9 table contains the columns needed for Gypsy/Copia analyses."""
    required = [
        "CHROM","WIN_START","WIN_END","WIN_MID","WIN_LEN",
        "q_poi_3vs4","q_poi_3vs5","q_poi_4vs5",
        "TE_Gypsy_LTR_retrotransposon","TE_Copia_LTR_retrotransposon",
        "TEdens_Gypsy_LTR_retrotransposon","TEdens_Copia_LTR_retrotransposon"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        # We allow the TEdens_* to be missing (we'll fallback to count/WIN_LEN)
        fallback_ok = all(c in ["TEdens_Gypsy_LTR_retrotransposon","TEdens_Copia_LTR_retrotransposon"] for c in missing)
        if not fallback_ok:
            raise ValueError(f"Missing required columns: {missing}")

def add_gypsy_copia_columns(df):
    """Create Gypsy/Copia count & density columns with fallbacks when densities are absent."""
    df = df.copy()
    # counts
    df["Gypsy_count"] = pd.to_numeric(df["TE_Gypsy_LTR_retrotransposon"], errors="coerce").fillna(0).astype(float)
    df["Copia_count"] = pd.to_numeric(df["TE_Copia_LTR_retrotransposon"],  errors="coerce").fillna(0).astype(float)
    # densities
    win_len = pd.to_numeric(df["WIN_LEN"], errors="coerce").fillna((df["WIN_END"]-df["WIN_START"]+1)).clip(lower=1).astype(float)
    if "TEdens_Gypsy_LTR_retrotransposon" in df.columns:
        gy_den = pd.to_numeric(df["TEdens_Gypsy_LTR_retrotransposon"], errors="coerce")
    else:
        gy_den = np.nan
    if "TEdens_Copia_LTR_retrotransposon" in df.columns:
        co_den = pd.to_numeric(df["TEdens_Copia_LTR_retrotransposon"], errors="coerce")
    else:
        co_den = np.nan

    df["Gypsy_density"] = np.where(np.isfinite(gy_den), gy_den, df["Gypsy_count"]/win_len)
    df["Copia_density"] = np.where(np.isfinite(co_den), co_den, df["Copia_count"]/win_len)
    return df

def get_sig_mask(df, pair, qcut):
    """Return a boolean mask where Poisson q-values for ``pair`` are below ``qcut``."""
    col = f"q_poi_{pair}"
    q = pd.to_numeric(df[col], errors="coerce")
    return (q < qcut)

def fisher_presence(sig_mask, has_feat):
    """Fisher's exact test for presence (any>0) between sig vs non-sig windows."""
    a = int((sig_mask & has_feat).sum())
    b = int((sig_mask & ~has_feat).sum())
    c = int((~sig_mask & has_feat).sum())
    d = int((~sig_mask & ~has_feat).sum())
    if (a+b == 0) or (c+d == 0):
        return {"a":a,"b":b,"c":c,"d":d,"odds":np.nan,"p":np.nan}
    odds, p = fisher_exact([[a,b],[c,d]], alternative="two-sided")
    return {"a":a,"b":b,"c":c,"d":d,"odds":odds,"p":p}

def glm_interaction_counts(df, pair, family="poisson", qcut=0.05):
    """Fit Gypsy-vs-Copia GLMs with a significance interaction term."""
    d = df.copy()
    # ensure numeric
    d["Gypsy_count"] = pd.to_numeric(d["Gypsy_count"], errors="coerce").fillna(0).astype(float)
    d["Copia_count"] = pd.to_numeric(d["Copia_count"], errors="coerce").fillna(0).astype(float)
    length = pd.to_numeric(d["WIN_LEN"], errors="coerce").fillna((d["WIN_END"]-d["WIN_START"]+1)).clip(lower=1).astype(float)
    sig = get_sig_mask(d, pair, qcut).astype(int)

    long = pd.DataFrame({
        "count": np.concatenate([d["Gypsy_count"].values, d["Copia_count"].values]).astype(float),
        "te_is_gypsy": np.concatenate([np.ones(len(d), dtype=int), np.zeros(len(d), dtype=int)]),
        "is_sig": np.concatenate([sig.values, sig.values]).astype(int),
        "length": np.concatenate([length.values, length.values]).astype(float)
    })
    # model matrix
    X = pd.DataFrame({
        "intercept": 1.0,
        "te_is_gypsy": long["te_is_gypsy"].astype(int),  # Copia=0 (reference), Gypsy=1
        "is_sig": long["is_sig"].astype(int),
        "int_te_sig": (long["te_is_gypsy"] * long["is_sig"]).astype(int)
    })
    off = np.log(long["length"].clip(lower=1.0))
    y = long["count"].astype(float)

    fam = Poisson(LogLink()) if family=="poisson" else NegativeBinomial(LogLink())
    fit = sm.GLM(y, X, family=fam, offset=off).fit()

    # dispersion (Pearson chi^2 / df)
    mu = fit.fittedvalues
    resid = y - mu
    if family=="poisson":
        var = mu
    else:
        var = fit.family.variance(mu)  # NB variance from family object
    ok = var > 0
    dfree = max(1, ok.sum() - X.shape[1])
    dispersion = float(np.sum((resid[ok]**2)/var[ok]) / dfree)

    # interaction stats
    b = fit.params["int_te_sig"]; se = fit.bse["int_te_sig"]
    if se > 0 and np.isfinite(b):
        z = b/se
        p = 2*(1 - 0.5*(1+math.erf(abs(z)/math.sqrt(2))))
    else:
        z = np.nan; p = np.nan
    rr = float(np.exp(b)) if np.isfinite(b) else np.nan

    # main effects (optional context)
    b_gy = fit.params["te_is_gypsy"]; se_gy = fit.bse["te_is_gypsy"]
    p_gy = 2*(1 - 0.5*(1+math.erf(abs(b_gy/se_gy)/math.sqrt(2)))) if se_gy>0 else np.nan
    b_sig = fit.params["is_sig"]; se_sig = fit.bse["is_sig"]
    p_sig = 2*(1 - 0.5*(1+math.erf(abs(b_sig/se_sig)/math.sqrt(2)))) if se_sig>0 else np.nan

    return {
        "dispersion": dispersion,
        "beta_gypsy_vs_copia": float(b_gy), "p_gypsy_vs_copia": float(p_gy),
        "beta_is_sig": float(b_sig), "p_is_sig": float(p_sig),
        "beta_int": float(b), "se_int": float(se), "z_int": float(z),
        "p_int": float(p), "RR_int": float(rr)
    }

def save_quick_plot(df, pair, outdir, qcut):
    """Small scatter of counts by FST significance for sanity."""
    sig = get_sig_mask(df, pair, qcut).values
    x_jitter = (sig.astype(int) + np.random.uniform(-0.07,0.07,len(sig)))
    fig, ax = plt.subplots(1,1, figsize=(7,4))
    ax.scatter(x_jitter, df["Gypsy_count"], s=6, alpha=0.4, label="Gypsy")
    ax.scatter(x_jitter, df["Copia_count"], s=6, alpha=0.4, label="Copia")
    ax.set_xticks([0,1]); ax.set_xticklabels(["Not-sig","Sig"])
    ax.set_ylabel("TE count per window")
    ax.set_title(f"Gypsy vs Copia counts by FST significance — {pair}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"gypsy_copia_counts_{pair}.png"), dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    args   = parse_args()
    INFILE = args.infile
    OUTDIR = args.outdir
    Q_CUT  = args.qcut
    MAKE_P = args.make_plots

    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, "figures"), exist_ok=True)

    df = pd.read_csv(INFILE, sep="\t")
    check_required_columns(df)
    df = add_gypsy_copia_columns(df)

    # Add FST significance flags (Poisson q)
    for pair in PAIRS:
        df[f"is_sig_{pair}"] = get_sig_mask(df, pair, Q_CUT).astype(int)

    # Save augmented table
    out_aug = os.path.join(OUTDIR, "fst_windows_with_TEcounts__gypsy_copia.tsv")
    df.to_csv(out_aug, sep="\t", index=False)

    # 1) Fisher presence tests (Gypsy, Copia) per pair
    fisher_rows = []
    for pair in PAIRS:
        sig = df[f"is_sig_{pair}"] == 1
        has_gy = (pd.to_numeric(df["Gypsy_count"], errors="coerce").fillna(0) > 0)
        has_co = (pd.to_numeric(df["Copia_count"], errors="coerce").fillna(0) > 0)

        g = fisher_presence(sig, has_gy); g.update({"pair":pair,"te":"Gypsy"})
        c = fisher_presence(sig, has_co); c.update({"pair":pair,"te":"Copia"})
        fisher_rows.extend([g,c])

    fisher_df = pd.DataFrame(fisher_rows, columns=["pair","te","a","b","c","d","odds","p"])
    # Re-compute columns ordering if needed (we set after since we updated dicts)
    if fisher_df.isnull().all(axis=None):
        fisher_df = pd.DataFrame(fisher_rows)
    # Add BH per pair (2 tests per pair)
    fisher_df["q_bh_perpair"] = (fisher_df.groupby("pair")["p"]
                                 .transform(lambda p: (p * p.notna().sum() / p.rank(method="min")).clip(upper=1.0)))
    fisher_out = os.path.join(OUTDIR, "gypsy_copia_fisher_presence.tsv")
    fisher_df.to_csv(fisher_out, sep="\t", index=False)

    # 2) Poisson & NB GLM with interaction per pair
    glm_rows = []
    for pair in PAIRS:
        res_p = glm_interaction_counts(df, pair, family="poisson", qcut=Q_CUT)
        res_n = glm_interaction_counts(df, pair, family="negbin",  qcut=Q_CUT)
        res_p.update({"pair":pair,"family":"poisson"})
        res_n.update({"pair":pair,"family":"negbin"})
        glm_rows.extend([res_p, res_n])

    glm_df = pd.DataFrame(glm_rows)[[
        "pair","family","dispersion",
        "beta_gypsy_vs_copia","p_gypsy_vs_copia",
        "beta_is_sig","p_is_sig",
        "beta_int","se_int","z_int","p_int","RR_int"
    ]].copy()

    # BH-FDR for interaction p across the 3 pairs within each family
    for fam in ["poisson","negbin"]:
        mask = glm_df["family"]==fam
        pseries = glm_df.loc[mask, "p_int"]
        q = (pseries * pseries.notna().sum() / pseries.rank(method="min")).clip(upper=1.0)
        glm_df.loc[mask, "q_int_bh_across_pairs"] = q

    glm_out = os.path.join(OUTDIR, "gypsy_copia_glm_interaction.tsv")
    glm_df.to_csv(glm_out, sep="\t", index=False)

    # Optional quick plots
    if MAKE_P:
        for pair in PAIRS:
            save_quick_plot(df, pair, os.path.join(OUTDIR,"figures"), Q_CUT)

    print(f"[✓] Wrote:\n  {out_aug}\n  {fisher_out}\n  {glm_out}")

