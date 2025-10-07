#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Windows jointly enriched for high FST and TE (two major TE groups)

Inputs
------
1) FST windows table (exact columns you listed):
   CHROM,WIN_START,WIN_END,C0_3vs4,C0.25_3vs4,C0_3vs5,C0.25_3vs5,C0_4vs5,C0.25_4vs5,
   R_3vs4,R_3vs5,R_4vs5,p_glob_3vs4,lambda_3vs4,disp_bin_3vs4,disp_poi_3vs4,rho_bb_3vs4,
   theta_nb_3vs4,p_bin_3vs4,q_bin_3vs4,p_poi_3vs4,q_poi_3vs4,p_bb_3vs4,q_bb_3vs4,
   p_nb_3vs4,q_nb_3vs4, ... (same for 3vs5, 4vs5), WIN_MID

2) TE TSV with headers:
   seqid, source, sequence_ontology, start, end, score, strand

What it outputs
---------------
- fst_windows_with_TEcounts_twoGroups.tsv  (per window, all TE counts, densities, and 2-group tallies)
- stats_twoGroups.tsv                      (Wilcoxon / GLM / Fisher summaries per pair & feature)
- joint_enrichment_twoGroups.tsv           (windows passing FST q<cut and TE-high criteria)
- figures/manhattan_q_poi_{pair}.png       (−log10 q from FST Poisson, 3-row panels)
- figures/violin_density_{pair}.png        (TE density by sig vs non-sig, two groups)
- figures/violin_counts_{pair}.png         (TE counts by sig vs non-sig, two groups)
"""

import os, math, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from statsmodels.api import GLM
from statsmodels.genmod.families import Poisson, NegativeBinomial
from statsmodels.genmod.families.links import log as LogLink
from scipy.stats import mannwhitneyu, fisher_exact

# ---------------- user knobs ----------------
FST_CSV   = "fst_window_counts_tests_win10000000_ov2000000.csv"   # change if needed
TE_TSV    = os.path.expanduser("~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv")
OUT_DIR   = "two_groups_out"
Q_CUT     = 0.05                         # significance for FST windows (Poisson q_<pair>)
ALPHA     = 0.05                         # for TE tests display
TOP_TE_Q  = 0.10                         # top 10% TE density to call "TE-high" (fallback criterion)

# Define the two major TE groups by keywords found in `sequence_ontology`
# Edit these lists to match your file's labels:
GROUPS = {
    "LTR_retro": ["LTR", "Gypsy", "Copia", "ERV", "retrotransposon"],
    "DNA_TIR_Helitron": ["TIR", "CACTA", "hAT", "Mutator", "PIF", "MULE", "Helitron", "Helitron_NA"]
}
# -------------- end knobs -------------------

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "figures"), exist_ok=True)

# -------- helpers --------
def load_fst(fname):
    df = pd.read_csv(fname)
    # required columns check (abbrev)
    must = ["CHROM","WIN_START","WIN_END","q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"]
    missing = [c for c in must if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required FST columns: {missing}")
    return df

def load_te(fname):
    te = pd.read_csv(fname, sep="\t", header=0)
    te = te.rename(columns=str.strip)
    # enforce expected header names
    need = ["seqid","sequence_ontology","start","end"]
    miss = [c for c in need if c not in te.columns]
    if miss:
        raise ValueError(f"TE file missing columns: {miss}")
    te["start"] = te["start"].astype(int)
    te["end"]   = te["end"].astype(int)
    return te

def scaffold_order(chroms):
    # order like scaffold_1a, 1b, 2, 3, ...
    def keyer(s):
        # extract integer and trailing letter if present
        base = s
        if s.startswith("scaffold_"): base = s[len("scaffold_"):]
        # split “1a” to (1, 'a'); “10” to (10, '')
        num = ""
        suf = ""
        for ch in base:
            if ch.isdigit(): num += ch
            else: suf += ch
        n = int(num) if num else 0
        return (n, suf)
    return sorted(list(set(chroms)), key=keyer)

def add_genome_x(df, chrom_col="CHROM", pos_col="WIN_MID"):
    # build cumulative offsets for Manhattan x
    order = scaffold_order(df[chrom_col].unique())
    offsets = {}
    acc = 0
    for c in order:
        length = df.loc[df[chrom_col]==c, pos_col].max()
        offsets[c] = acc
        acc += length
    x = df[pos_col] + df[chrom_col].map(offsets).astype(int)
    df = df.assign(x_coord = x/1e9)  # in Gb
    return df, order

def label_two_groups(so):
    s = str(so).lower()
    for g, keys in GROUPS.items():
        for k in keys:
            if k.lower() in s:
                return g
    return "other"

def interval_overlap_count(te_df, windows_df):
    """
    Count TE per window (all + per sequence_ontology + two major groups).
    Assumes per-scaffold processing for speed.
    """
    out = []
    # precompute ontology and group
    te_df = te_df.copy()
    te_df["te_len"] = te_df["end"] - te_df["start"] + 1
    te_df["group2"] = te_df["sequence_ontology"].apply(label_two_groups)

    for chrom in scaffold_order(windows_df["CHROM"].unique()):
        w = windows_df[windows_df["CHROM"]==chrom].copy()
        t = te_df[te_df["seqid"]==chrom].copy()
        if t.empty:
            w["TE_total"] = 0
            w["TE_bp"]    = 0
            w["TE_density"] = 0.0
            # per 2-group zeros
            for g in GROUPS.keys():
                w[f"TEcnt_{g}"] = 0
                w[f"TEbp_{g}"]  = 0
                w[f"TEdens_{g}"]= 0.0
            out.append(w); continue

        # For faster overlap, sort and sweep
        t = t.sort_values("start").reset_index(drop=True)
        t_starts = t["start"].values
        t_ends   = t["end"].values
        t_len    = t["te_len"].values
        t_so     = t["sequence_ontology"].values
        t_grp    = t["group2"].values

        # Build small index by start bins (1 Mb bins) to cut comparisons
        bin_size = 1_000_000
        bins = defaultdict(list)
        for i,(a,b) in enumerate(zip(t_starts, t_ends)):
            b0 = a // bin_size
            b1 = b // bin_size
            for bidx in range(b0, b1+1):
                bins[bidx].append(i)

        # compute overlaps
        per_so_keys = sorted(t["sequence_ontology"].dropna().astype(str).unique())
        for col in per_so_keys:
            w[f"TEcnt_{col}"] = 0
        w["TE_total"]=0; w["TE_bp"]=0

        for g in GROUPS.keys():
            w[f"TEcnt_{g}"]=0; w[f"TEbp_{g}"]=0

        W = w.to_dict("records")
        for rec in W:
            a = int(rec["WIN_START"]); b = int(rec["WIN_END"])
            # candidate TE indices from start-bin slice
            cand = set()
            for bidx in range(a // bin_size, b // bin_size + 1):
                cand.update(bins.get(bidx, []))

            te_bp = 0; te_cnt = 0
            grp_bp = {g:0 for g in GROUPS.keys()}
            grp_ct = {g:0 for g in GROUPS.keys()}
            so_ct  = defaultdict(int)

            for i in cand:
                aa = t_starts[i]; bb = t_ends[i]
                if bb < a or aa > b: 
                    continue
                ov = min(b,bb) - max(a,aa) + 1
                if ov <= 0: 
                    continue
                te_cnt += 1
                te_bp  += ov
                g = t_grp[i]
                if g in grp_bp:
                    grp_bp[g] += ov
                    grp_ct[g] += 1
                so_ct[t_so[i]] += 1

            rec["TE_total"] = te_cnt
            rec["TE_bp"]    = te_bp
            win_len = (b - a + 1)
            rec["TE_density"] = te_bp / win_len if win_len>0 else 0.0
            for g in GROUPS.keys():
                rec[f"TEcnt_{g}"] = grp_ct[g]
                rec[f"TEbp_{g}"]  = grp_bp[g]
                rec[f"TEdens_{g}"]= grp_bp[g]/win_len if win_len>0 else 0.0
            # per-class counts
            for so in so_ct:
                rec[f"TEcnt_{so}"] = so_ct[so]

        w = pd.DataFrame(W)
        out.append(w)

    res = pd.concat(out, ignore_index=True)
    return res

def is_sig_fst(row, pair):
    return (row[f"q_poi_{pair}"] < Q_CUT) if f"q_poi_{pair}" in row else False

def glm_counts(y, group_sig, offset_len, family="poisson"):
    # model: count ~ is_sig + offset(log(length))
    X = pd.DataFrame({"intercept":1.0, "is_sig": group_sig.astype(int)})
    off = np.log(np.clip(offset_len.astype(float), 1.0, None))
    fam = Poisson(LogLink()) if family=="poisson" else NegativeBinomial(LogLink())
    m = GLM(y.astype(float), X, family=fam, offset=off).fit()
    beta = m.params["is_sig"]; se = m.bse["is_sig"]
    rr = math.exp(beta)
    # Wald p
    z  = beta/se if se>0 else np.nan
    p  = 2*(1 - 0.5*(1+math.erf(abs(z)/math.sqrt(2)))) if np.isfinite(z) else np.nan
    return {"rr": rr, "beta": beta, "se": se, "p": p}

# -------- main workflow --------
fst = load_fst(FST_CSV)
te  = load_te(TE_TSV)

# TE overlaps per window
tew = interval_overlap_count(te, fst[["CHROM","WIN_START","WIN_END","WIN_MID"]])

# Merge back to fst table
aug = fst.merge(tew, on=["CHROM","WIN_START","WIN_END","WIN_MID"], how="left")
for pair in ["3vs4","3vs5","4vs5"]:
    aug[f"is_sig_{pair}"] = aug.apply(is_sig_fst, axis=1, pair=pair)

# Save the augmented table
out_main = os.path.join(OUT_DIR, "fst_windows_with_TEcounts_twoGroups.tsv")
aug.to_csv(out_main, sep="\t", index=False)

# ------------- stats and figures per pair -------------
stats_rows = []
pairs = ["3vs4","3vs5","4vs5"]
features_density = ["TE_density"] + [f"TEdens_{g}" for g in GROUPS.keys()]
features_counts  = ["TE_total"]  + [f"TEcnt_{g}"  for g in GROUPS.keys()]

for pair in pairs:
    is_sig = aug[f"is_sig_{pair}"].values
    length = (aug["WIN_END"] - aug["WIN_START"] + 1).clip(lower=1)

    # Wilcoxon (density) + effect size (rank-biserial)
    for feat in features_density:
        x = aug.loc[ is_sig, feat].astype(float)
        y = aug.loc[~is_sig, feat].astype(float)
        if len(x)>0 and len(y)>0:
            u, p = mannwhitneyu(x, y, alternative="two-sided")
            # rank-biserial = 1 - 2*U/(n1*n2)
            rbes = 1 - 2*u/(len(x)*len(y))
        else:
            p = np.nan; u = np.nan; rbes = np.nan
        stats_rows.append({"pair":pair,"type":"wilcoxon_density","feature":feat,"p":p,"U":u,
                           "n_sig":len(x),"n_non":len(y),"med_sig":np.median(x) if len(x) else np.nan,
                           "med_non":np.median(y) if len(y) else np.nan,"rbes":rbes})

    # GLM Poisson & NegBin on counts with offset
    for feat in features_counts:
        ycount = aug[feat].astype(float)
        res_p  = glm_counts(ycount, aug[f"is_sig_{pair}"], length, family="poisson")
        res_nb = glm_counts(ycount, aug[f"is_sig_{pair}"], length, family="negbin")
        stats_rows.append({"pair":pair,"type":"glm_poisson_counts","feature":feat, **res_p})
        stats_rows.append({"pair":pair,"type":"glm_negbin_counts","feature":feat, **res_nb})

    # Fisher (any TE present)
    for feat in features_counts:
        any_te = (aug[feat].fillna(0) > 0).astype(int)
        a = int(((is_sig==1)&(any_te==1)).sum())
        b = int(((is_sig==1)&(any_te==0)).sum())
        c = int(((is_sig==0)&(any_te==1)).sum())
        d = int(((is_sig==0)&(any_te==0)).sum())
        if (a+b>0) and (c+d>0):
            od, p = fisher_exact([[a,b],[c,d]], alternative="two-sided")
        else:
            od, p = (np.nan, np.nan)
        stats_rows.append({"pair":pair,"type":"fisher_presence","feature":feat,
                           "odds":od,"p":p,"table":json.dumps({"sig":[a,b],"nonsig":[c,d]})})

# BH-FDR within each (pair,type,feature)
stats = pd.DataFrame(stats_rows)
def bh(df, pcol="p"):
    p = df[pcol].astype(float)
    m = p.notna().sum()
    ranks = p.rank(method="min")
    q = p * m / ranks
    df["q_bh"] = q.clip(upper=1.0)
    return df
stats = stats.groupby(["pair","type","feature"], group_keys=False).apply(bh).reset_index(drop=True)
stats_out = os.path.join(OUT_DIR, "stats_twoGroups.tsv")
stats.to_csv(stats_out, sep="\t", index=False)

# -------- joint enrichment: FST q<cut AND TE-high ----------
joint_rows = []
for pair in pairs:
    sig_mask = aug[f"is_sig_{pair}"]
    # choose a TE metric to rank: total density and each group density
    for feat in ["TE_density"] + [f"TEdens_{g}" for g in GROUPS.keys()]:
        thr = aug[feat].quantile(1 - TOP_TE_Q)
        tehigh = aug[feat] >= thr
        joint = aug[sig_mask & tehigh].copy()
        joint["pair"] = pair; joint["te_feature"] = feat; joint["te_thr"] = thr
        joint_rows.append(joint)
joint = pd.concat(joint_rows, ignore_index=True)
joint_out = os.path.join(OUT_DIR, "joint_enrichment_twoGroups.tsv")
joint_cols = ["pair","te_feature","CHROM","WIN_START","WIN_END","WIN_MID","TE_total","TE_bp","TE_density",
              "TEcnt_LTR_retro","TEcnt_DNA_TIR_Helitron","TEdens_LTR_retro","TEdens_DNA_TIR_Helitron",
              "q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"]
joint[joint_cols].to_csv(joint_out, sep="\t", index=False)

# ------------- plots: Manhattan (FST q) -------------
def manhattan_q(df, pair, out_png, title):
    sub = df.copy()
    sub["q"] = sub[f"q_poi_{pair}"].astype(float)
    sub["mlogq"] = -np.log10(sub["q"].clip(lower=1e-300))
    sub["sig"] = sub["q"] < Q_CUT
    sub["WIN_MID"] = ((sub["WIN_START"]+sub["WIN_END"])//2).astype(int)
    sub, order = add_genome_x(sub, "CHROM", "WIN_MID")
    # alternating colors by scaffold
    scaffold_to_color = {c: ("#9ecae1" if i%2==0 else "#bdbdbd") for i,c in enumerate(order)}
    colr = sub["CHROM"].map(scaffold_to_color)
    alpha = np.where(sub["sig"], 1.0, 0.3)

    plt.figure(figsize=(12, 7))
    ax = plt.gca()
    ax.scatter(sub["x_coord"], sub["mlogq"], s=8, c=colr, alpha=alpha, linewidths=0)
    ax.axhline(-math.log10(Q_CUT), color="k", lw=1)
    ax.set_xlabel("Genome position (Gb, scaffolds concatenated)")
    ax.set_ylabel("−log10(q)")
    ax.set_title(title)
    # annotate top 5
    top = sub.nsmallest(5, f"q_poi_{pair}")
    for _,r in top.iterrows():
        ax.scatter(r["x_coord"], r["mlogq"], s=30, c="red", zorder=5)
        ax.text(r["x_coord"], r["mlogq"]+0.2,
                f"{r['CHROM']}_{r['WIN_START']}-{r['WIN_END']}",
                color="red", fontsize=8, rotation=45, ha="left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300); plt.close()

for pair in pairs:
    manhattan_q(aug, pair,
                os.path.join(OUT_DIR,"figures",f"manhattan_q_poi_{pair}.png"),
                f"FST Poisson q — {pair}")

# -------- violin/box: density & counts (two groups) ----------
def violin_by_sig(df, pair, feats, kind="density"):
    fig, axes = plt.subplots(len(feats), 1, figsize=(10, 3.2*len(feats)), sharex=True)
    if len(feats)==1: axes=[axes]
    for ax, feat in zip(axes, feats):
        data = [df.loc[df[f"is_sig_{pair}"], feat].values,
                df.loc[~df[f"is_sig_{pair}"], feat].values]
        ax.violinplot(data, showmedians=True)
        ax.set_title(f"{pair} — {feat}")
        ax.set_xticks([1,2]); ax.set_xticklabels(["Sig", "Not"])
        ax.set_ylabel(feat)
    axes[-1].set_xlabel("Window class")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR,"figures",f"violin_{kind}_{pair}.png"), dpi=300)
    plt.close()

for pair in pairs:
    violin_by_sig(aug, pair, features_density, "density")
    violin_by_sig(aug, pair, features_counts,  "counts")

print(f"Done. Tables in {OUT_DIR}")

