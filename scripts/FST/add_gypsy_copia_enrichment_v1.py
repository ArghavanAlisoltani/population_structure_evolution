#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gypsy vs Copia enrichment tied to significant FST windows
---------------------------------------------------------

This script augments the v9 output table (fst_windows_with_TEcounts.tsv) by:
  * deriving Gypsy and Copia counts/densities per window
  * comparing enrichment in significant vs non-significant FST windows
      - Fisher presence test (Gypsy separately; Copia separately)
      - Poisson GLM (counts ~ te_type + is_sig + te_type:is_sig) with log(window_length) offset
      - Negative Binomial GLM for overdispersion robustness
  * writing results to tidy TSVs and producing a small summary plot per pair (optional)

Inputs
------
A TSV like your v9 output, expected to contain at least:
  CHROM, WIN_START, WIN_END, WIN_MID,
  q_poi_3vs4, q_poi_3vs5, q_poi_4vs5,
  ... and TE per-class columns such as TEcnt_* (from v9), OR at least the class-wise counts
      (e.g., TEcnt_Gypsy, TEcnt_Copia) — if not present, we infer from TEcnt_* columns by name.

Usage
-----
python add_gypsy_copia_enrichment.py \
  --in fst_windows_with_TEcounts.tsv \
  --outdir gypsy_copia_addon \
  --qcut 0.05

"""

import os
import re
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True,
                    help="Path to fst_windows_with_TEcounts.tsv (from v9)")
    ap.add_argument("--outdir", default="gypsy_copia_addon",
                    help="Output directory for TSVs/figures")
    ap.add_argument("--qcut", type=float, default=0.05,
                    help="q-value cutoff for FST Poisson significance")
    ap.add_argument("--make_plots", action="store_true",
                    help="Emit small summary plots per pair")
    return ap.parse_args()

def find_class_columns(df, prefix="TEcnt_"):
    """Return list of TE count columns with the given prefix."""
    return [c for c in df.columns if c.startswith(prefix)]

def derive_gypsy_copia(df):
    """
    Create Gypsy/Copia count & density columns.

    Strategy:
      1) If TEcnt_Gypsy / TEcnt_Copia already exist, use them.
      2) Else, sum any TEcnt_* columns whose name includes 'gypsy' or 'copia' (case-insensitive).
      3) Densities are bp-normalized if TEbp_* exist; otherwise count/length.

    Returns df with added:
      'Gypsy_count', 'Copia_count', 'Gypsy_density', 'Copia_density'
    """
    df = df.copy()

    # Window length
    if not {"WIN_START","WIN_END"}.issubset(df.columns):
        raise ValueError("Input must contain WIN_START and WIN_END")
    win_len = (df["WIN_END"] - df["WIN_START"] + 1).clip(lower=1)

    # Prefer explicit columns if present
    if "TEcnt_Gypsy" in df.columns: gypsy_cnt = df["TEcnt_Gypsy"].fillna(0).astype(float)
    else:
        # Sum any TEcnt_* whose suffix contains 'gypsy'
        cnt_cols = find_class_columns(df, "TEcnt_")
        gypsy_like = [c for c in cnt_cols if re.search(r"gypsy", c, flags=re.I)]
        if not gypsy_like:
            raise ValueError("Could not find any TE columns that look like Gypsy (e.g., 'TEcnt_Gypsy').")
        gypsy_cnt = df[gypsy_like].fillna(0).sum(axis=1).astype(float)

    if "TEcnt_Copia" in df.columns: copia_cnt = df["TEcnt_Copia"].fillna(0).astype(float)
    else:
        cnt_cols = find_class_columns(df, "TEcnt_")
        copia_like = [c for c in cnt_cols if re.search(r"copia", c, flags=re.I)]
        if not copia_like:
            raise ValueError("Could not find any TE columns that look like Copia (e.g., 'TEcnt_Copia').")
        copia_cnt = df[copia_like].fillna(0).sum(axis=1).astype(float)

    # Densities
    # If bp-overlap columns exist, use them; else counts per base
    if "TEbp_Gypsy" in df.columns:
        gypsy_den = (df["TEbp_Gypsy"].fillna(0).astype(float) / win_len)
    else:
        # fall back to count density
        gypsy_den = gypsy_cnt / win_len

    if "TEbp_Copia" in df.columns:
        copia_den = (df["TEbp_Copia"].fillna(0).astype(float) / win_len)
    else:
        copia_den = copia_cnt / win_len

    df["Gypsy_count"]   = gypsy_cnt
    df["Copia_count"]   = copia_cnt
    df["Gypsy_density"] = gypsy_den
    df["Copia_density"] = copia_den
    return df

def get_sig_mask(df, pair, qcut):
    col = f"q_poi_{pair}"
    if col not in df.columns:
        raise ValueError(f"Missing FST q-value column: {col}")
    q = pd.to_numeric(df[col], errors="coerce")
    return (q < qcut)

def fisher_presence(sig_mask, has_feat):
    """
    Fisher's exact test for presence (any>0) between sig vs non-sig windows.
    Returns dict with a,b,c,d, odds, p.
      a = sig & present
      b = sig & absent
      c = non-sig & present
      d = non-sig & absent
    """
    a = int((sig_mask & has_feat).sum())
    b = int((sig_mask & ~has_feat).sum())
    c = int((~sig_mask & has_feat).sum())
    d = int((~sig_mask & ~has_feat).sum())
    if (a+b == 0) or (c+d == 0):
        return {"a":a,"b":b,"c":c,"d":d,"odds":np.nan,"p":np.nan}
    odds, p = fisher_exact([[a,b],[c,d]], alternative="two-sided")
    return {"a":a,"b":b,"c":c,"d":d,"odds":odds,"p":p}

def glm_interaction_counts(df, pair, family="poisson"):
    """
    Long-format GLM:
      count ~ te_type (Gypsy/Copia) + is_sig + te_type:is_sig
      offset = log(window_length)

    Returns dict with:
      disp (Pearson χ²/df), beta_int, se_int, p_int (Wald), RR_int (exp(beta_int)),
      plus overall coef table for reference.
    """
    # build long frame
    length = (df["WIN_END"] - df["WIN_START"] + 1).clip(lower=1).astype(float)
    sig = get_sig_mask(df, pair, qcut=1.0)  # we'll pass true qcut elsewhere; here just getting column exists
    sig = (pd.to_numeric(df[f"q_poi_{pair}"], errors="coerce") < Q_CUT).astype(int)

    long = pd.DataFrame({
        "count": np.concatenate([df["Gypsy_count"].values, df["Copia_count"].values]).astype(float),
        "te_type": np.array(["Gypsy"]*len(df) + ["Copia"]*len(df)),
        "is_sig": np.concatenate([sig.values, sig.values]).astype(int),
        "length": np.concatenate([length.values, length.values]).astype(float)
    })

    X = pd.get_dummies(long[["te_type","is_sig"]], drop_first=True)
    # columns: te_type_Gypsy? (Copia dropped) → we want Copia as reference; so re-encode:
    # Make baseline te_type=Copia, so add explicit coding:
    X = pd.DataFrame({
        "intercept": 1.0,
        "te_is_gypsy": (long["te_type"]=="Gypsy").astype(int),
        "is_sig": long["is_sig"].astype(int),
        "int_te_sig": ((long["te_type"]=="Gypsy") & (long["is_sig"]==1)).astype(int)
    })

    off = np.log(long["length"].clip(lower=1.0))
    y = long["count"].astype(float)

    fam = Poisson(LogLink()) if family=="poisson" else NegativeBinomial(LogLink())
    model = sm.GLM(y, X, family=fam, offset=off)
    fit   = model.fit()

    # Pearson dispersion
    mu = fit.fittedvalues
    resid = y - mu
    if family=="poisson":
        var = mu
    else:
        # NB variance from statsmodels is mu + alpha*mu^2 ; alpha = 1/theta in some paramizations.
        # statsmodels stores 'scale' differently; use fitted var from family:
        var = fit.family.variance(mu)
    ok = var > 0
    dfree = max(1, ok.sum() - X.shape[1])
    disp = float(np.sum((resid[ok]**2)/var[ok]) / dfree)

    # interaction is coefficient on int_te_sig (Gypsy x Sig)
    b = fit.params["int_te_sig"]; se = fit.bse["int_te_sig"]
    z = np.nan
    p = np.nan
    if se > 0 and np.isfinite(b) and np.isfinite(se):
        z = b/se
        # two-sided normal approx
        p = 2*(1 - 0.5*(1+math.erf(abs(z)/math.sqrt(2))))
    rr = float(np.exp(b)) if np.isfinite(b) else np.nan

    # also return main effect contrasts for context
    out = {
        "dispersion": disp,
        "beta_int": float(b), "se_int": float(se), "z_int": float(z), "p_int": float(p), "RR_int": rr,
        "beta_gypsy_vs_copia": float(fit.params["te_is_gypsy"]),
        "p_gypsy_vs_copia": float(2*(1 - 0.5*(1+math.erf(abs(fit.params["te_is_gypsy"]/fit.bse["te_is_gypsy"])/math.sqrt(2))))) if fit.bse["te_is_gypsy"]>0 else np.nan,
        "beta_is_sig": float(fit.params["is_sig"]),
        "p_is_sig": float(2*(1 - 0.5*(1+math.erf(abs(fit.params["is_sig"]/fit.bse["is_sig"])/math.sqrt(2))))) if fit.bse["is_sig"]>0 else np.nan
    }
    return out

def save_plots(df, pair, outdir, qcut):
    # Simple dot plot of Gypsy/Copia counts per window (optional)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    sig = get_sig_mask(df, pair, qcut)
    ax.scatter(np.where(sig, 1, 0)+np.random.uniform(-0.05,0.05,len(df)), df["Gypsy_count"], s=6, alpha=0.4, label="Gypsy")
    ax.scatter(np.where(sig, 1, 0)+np.random.uniform(-0.05,0.05,len(df)), df["Copia_count"], s=6, alpha=0.4, label="Copia")
    ax.set_xticks([0,1]); ax.set_xticklabels(["Not-sig","Sig"])
    ax.set_ylabel("Count per window")
    ax.set_title(f"Gypsy vs Copia counts by FST significance — {pair}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"gypsy_copia_counts_{pair}.png"), dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    args = parse_args()
    INFILE  = args.infile
    OUTDIR  = args.outdir
    Q_CUT   = args.qcut
    MAKE_P  = args.make_plots

    os.makedirs(OUTDIR, exist_ok=True)

    df = pd.read_csv(INFILE, sep="\t")
    # Derive Gypsy/Copia columns (counts + densities)
    df = derive_gypsy_copia(df)

    # Add/overwrite is_sig flags per pair (based on q_poi_<pair>)
    for pair in PAIRS:
        col = f"q_poi_{pair}"
        if col not in df.columns:
            raise SystemExit(f"Missing FST q-value column: {col}")
        df[f"is_sig_{pair}"] = (pd.to_numeric(df[col], errors="coerce") < Q_CUT).astype(int)

    # Save augmented table with Gypsy/Copia columns
    out_aug = os.path.join(OUTDIR, "fst_windows_with_TEcounts__gypsy_copia.tsv")
    df.to_csv(out_aug, sep="\t", index=False)

    # 1) Presence Fisher tests (Gypsy and Copia, each: sig vs non)
    fisher_rows = []
    for pair in PAIRS:
        sig = df[f"is_sig_{pair}"] == 1

        # presence flags
        has_gypsy = (df["Gypsy_count"] > 0)
        has_copia = (df["Copia_count"] > 0)

        g = fisher_presence(sig, has_gypsy)
        c = fisher_presence(sig, has_copia)
        g.update({"pair":pair, "te":"Gypsy"})
        c.update({"pair":pair, "te":"Copia"})
        fisher_rows.append(g); fisher_rows.append(c)

    fisher_df = pd.DataFrame(fisher_rows)[["pair","te","a","b","c","d","odds","p"]].copy()
    # BH-FDR per pair across the two tests (Gypsy/Copia)
    fisher_df["q_bh_perpair"] = (fisher_df
                                 .groupby("pair")["p"]
                                 .transform(lambda p: (p*len(p)/p.rank(method="min")).clip(upper=1.0)))
    fisher_out = os.path.join(OUTDIR, "gypsy_copia_fisher_presence.tsv")
    fisher_df.to_csv(fisher_out, sep="\t", index=False)

    # 2) Rate comparison with interaction: Poisson and NB
    glm_rows = []
    for pair in PAIRS:
        res_p = glm_interaction_counts(df, pair, family="poisson")
        res_n = glm_interaction_counts(df, pair, family="negbin")
        res_p.update({"pair":pair, "family":"poisson"})
        res_n.update({"pair":pair, "family":"negbin"})
        glm_rows.extend([res_p, res_n])

    glm_df = pd.DataFrame(glm_rows)[[
        "pair","family","dispersion",
        "beta_gypsy_vs_copia","p_gypsy_vs_copia",
        "beta_is_sig","p_is_sig",
        "beta_int","se_int","z_int","p_int","RR_int"
    ]].copy()

    # BH-FDR across the 3 pairs per family for interaction p
    for fam in ["poisson","negbin"]:
        mask = (glm_df["family"]==fam)
        p = glm_df.loc[mask, "p_int"]
        q = (p * p.notna().sum() / p.rank(method="min")).clip(upper=1.0)
        glm_df.loc[mask, "q_int_bh_across_pairs"] = q

    glm_out = os.path.join(OUTDIR, "gypsy_copia_glm_interaction.tsv")
    glm_df.to_csv(glm_out, sep="\t", index=False)

    # Optional quick plots
    if MAKE_P:
        for pair in PAIRS:
            save_plots(df, pair, OUTDIR, Q_CUT)

    print(f"[✓] Wrote:\n  {out_aug}\n  {fisher_out}\n  {glm_out}")
