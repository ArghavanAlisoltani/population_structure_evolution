#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TMRCA near-TE enrichment with LENGTH-MATCHED Fisher tests (Top 5% vs background).

This is a TE-focused clone of your original gene-based script, preserving the
same structure and flow while swapping in EDTA TE annotations. It performs:

1) Reads a tab-delimited TMRCA file with EXACT columns:
   CHROM, start_tmrca, end_tmrca, mean_tmrca, lower_CI, upper_CI,
   overlap_genes, overlap_gene_n, nearest_genes

   • Computes segment length = (end_tmrca - start_tmrca) + 1
   • Splits rows into TOP 5% vs BACKGROUND by mean_tmrca (configurable)

2) Reads an EDTA TE table (TSV) with columns:
   seqid, source, sequence_ontology, start, end, score, strand

   • Builds two flavors of TE sets:
       a) Any TE (all rows)
       b) Per-class TE sets inferred from `sequence_ontology`.
          e.g., "Copia_LTR_retrotransposon" → class "Copia";
                bare terms like "LTR_retrotransposon" are kept as-is.

3) Computes edge-to-edge *distance to nearest TE* (0 if any overlap), for:
     • Any TE, and
     • Each TE class.

4) Fisher exact tests (two-sided) for TOP5 vs BACKGROUND inside/outside a
   set of TE-distance thresholds (e.g., 0, 1k, 2k, 5k, 10k, 50k).

5) NEW: LENGTH-MATCHED Fisher tests (identical structure to your script)
   • Compute exact length frequencies, keep only lengths with freq ≥ MIN_LEN_FREQ.
   • For EACH kept length L, run the SAME Fisher set using ONLY rows with seg_len == L
     (so Top and Background are matched on length).
   • Return tidy tables with one row per (length × distance × TE class).

Outputs (in OUTDIR)
-------------------
- fisher_overall_TE.tsv                     (overall across all lengths, any + per class)
- fisher_by_length_TE.tsv                   (one row per length×distance×class)
- length_frequencies.tsv                    (all exact lengths and counts)
- summary.txt                               (cutoff and basic counts)

Notes
-----
- Coordinates are treated as 0-based half-open internally when computing distances.
- Dependencies: pandas, numpy, scipy, statsmodels.
- You can run this next to your original pipeline without changing it.

Author: you + ChatGPT
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =========================
# USER SETTINGS
# =========================
IN_TMRCA_TSV = "annotated_tmrca_4_GPT_13columns.tsv"   # required columns above
IN_EDTA_TSV  = "col8_readable.TEanno.cds.scaf1ab.txt"  # EDTA table

# Define Top 5% on mean_tmrca (set DROP_ZERO_FOR_CUTOFF=True to ignore zeros)
TOP_PERCENTILE         = 95.0
DROP_ZERO_FOR_CUTOFF   = True

# Distance thresholds (in bp) to use in Fisher tests (0 means “overlaps a TE”)
DIST_THRESH = [0, 1_000, 2_000, 5_000, 10_000, 50_000]

# Minimum exact length frequency required to include a length in the
# length-matched analysis (helps stability when counts are tiny)
MIN_LEN_FREQ = 10

# Output directory
OUTDIR = Path("tmrca_TE_fisher_out")
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================
# Loading
# =========================

def load_tmrca_table(path: str | Path) -> pd.DataFrame:
    need = ["CHROM","start_tmrca","end_tmrca","mean_tmrca","lower_CI","upper_CI",
            "overlap_genes","overlap_gene_n","nearest_genes"]
    df = pd.read_csv(path, sep="	", dtype={"CHROM":str})
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")
    # enforce numeric on coords and tmrca
    for c in ["start_tmrca","end_tmrca","mean_tmrca","lower_CI","upper_CI","overlap_gene_n"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["CHROM","start_tmrca","end_tmrca","mean_tmrca"]).copy()
    df["start_tmrca"] = df["start_tmrca"].astype(np.int64)
    df["end_tmrca"]   = df["end_tmrca"].astype(np.int64)
    df["seg_len"]     = (df["end_tmrca"] - df["start_tmrca"] + 1).clip(lower=0)
    return df


def load_edta(path: str | Path) -> pd.DataFrame:
    need = ["seqid","source","sequence_ontology","start","end","score","strand"]
    ed = pd.read_csv(path, sep="	")
    miss = [c for c in need if c not in ed.columns]
    if miss:
        raise ValueError(f"EDTA missing required columns: {miss}")
    ed = ed.rename(columns={"seqid":"CHROM"}).copy()
    for c in ["start","end"]:
        ed[c] = pd.to_numeric(ed[c], errors="coerce")
    ed = ed.dropna(subset=["CHROM","start","end"]).copy()
    ed["start"] = ed["start"].astype(np.int64)
    ed["end"]   = ed["end"].astype(np.int64)
    ed = ed[ed["end"] > ed["start"]].copy()
    # derive TE class
    ed["te_class"] = ed["sequence_ontology"].apply(infer_te_class)
    return ed


def infer_te_class(so: str) -> str:
    if isinstance(so, str):
        if so.endswith("_LTR_retrotransposon"):
            return so.split("_LTR_retrotransposon")[0] or "LTR_retrotransposon"
        if so.endswith("_TIR_transposon"):
            return so.split("_TIR_transposon")[0] or "TIR_transposon"
        if "_non-LTR_retrotransposon" in so:
            return so.split("_non-LTR_retrotransposon")[0] or "non-LTR_retrotransposon"
        return so
    return "unknown"

# =========================
# Distance to nearest TE (0 if overlap)
# =========================

def _prepare_interval_index(ed: pd.DataFrame) -> Dict[str, np.ndarray]:
    idx: Dict[str, np.ndarray] = {}
    for chrom, d in ed.groupby("CHROM", sort=False):
        arr = d[["start","end"]].to_numpy(np.int64)
        arr = arr[np.argsort(arr[:,0])]
        idx[chrom] = arr
    return idx


def compute_nearest_dist(df: pd.DataFrame, idx: Dict[str, np.ndarray],
                          start_col="start_tmrca", end_col="end_tmrca") -> pd.Series:
    out = np.full(len(df), np.nan, dtype=float)
    # For speed, iterate per chrom; use a moving pointer heuristic
    pos = {c:0 for c in idx.keys()}
    for chrom, d in df.groupby("CHROM", sort=False):
        arr = idx.get(chrom)
        if arr is None or len(arr) == 0:
            continue
        j = 0
        for row_idx, row in d.sort_values(start_col).iterrows():
            s = int(row[start_col]); e = int(row[end_col])
            # advance j while intervals end before s
            while j < len(arr) and arr[j,1] < s:
                j += 1
            # check nearby intervals j-1, j, j+1, j+2 for min distance
            best = np.inf
            for k in (j-2, j-1, j, j+1, j+2):
                if 0 <= k < len(arr):
                    a,b = int(arr[k,0]), int(arr[k,1])
                    if b >= s and a <= e:  # overlap
                        best = 0.0
                        break
                    gap = min(abs(s - b), abs(a - e))
                    if gap < best:
                        best = gap
            out[df.index.get_loc(row_idx)] = float(best)
    return pd.Series(out, index=df.index)

# =========================
# Fisher tests
# =========================

def fisher_table(sub: pd.DataFrame, flag_col: str,
                 dist_col: str, thresholds: List[int | float], te_class: str) -> pd.DataFrame:
    rows = []
    D = sub[np.isfinite(sub[dist_col].values)].copy()
    for thr in thresholds:
        top_in  = int((D[flag_col] & (D[dist_col] <= thr)).sum())
        top_out = int((D[flag_col] & (D[dist_col] >  thr)).sum())
        bkg_in  = int((~D[flag_col] & (D[dist_col] <= thr)).sum())
        bkg_out = int((~D[flag_col] & (D[dist_col] >  thr)).sum())
        # Fisher two-sided
        _, p = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="two-sided")
        top_rate = (top_in  / max(1, (top_in + top_out))) if (top_in+top_out)>0 else np.nan
        bkg_rate = (bkg_in  / max(1, (bkg_in + bkg_out))) if (bkg_in+bkg_out)>0 else np.nan
        rr = (top_rate / bkg_rate) if (np.isfinite(top_rate) and np.isfinite(bkg_rate) and bkg_rate>0) else np.nan
        rows.append({
            "te_class": te_class,
            "threshold_bp": thr,
            "TOP_in": top_in, "TOP_out": top_out,
            "BKG_in": bkg_in, "BKG_out": bkg_out,
            "TOP_rate": top_rate, "BKG_rate": bkg_rate,
            "rate_ratio": rr, "fisher_p": p
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_q"] = multipletests(out["fisher_p"].values, method="fdr_bh")[1]
    return out

# =========================
# Main
# =========================

def main():
    df = load_tmrca_table(IN_TMRCA_TSV)
    ed = load_edta(IN_EDTA_TSV)

    # Define Top5% vs Background on mean_tmrca
    vals = df["mean_tmrca"].replace([np.inf, -np.inf], np.nan)
    if DROP_ZERO_FOR_CUTOFF:
        vals = vals.mask(vals <= 0, np.nan)
    cutoff = float(np.nanpercentile(vals, TOP_PERCENTILE))
    df["TOP5"] = df["mean_tmrca"] >= cutoff

    # Distance to Any TE
    any_idx = _prepare_interval_index(ed)
    df["nearest_TE_any"] = compute_nearest_dist(df, any_idx)

    # Per-class distances
    by_class = {}
    for cl, dcl in ed.groupby("te_class", sort=False):
        idx = _prepare_interval_index(dcl)
        col = f"nearest_TE_{cl}"
        df[col] = compute_nearest_dist(df, idx)
        by_class[cl] = col

    # Save exact length frequencies
    len_freq = (df["seg_len"].value_counts().rename_axis("seg_len").reset_index(name="count")
                .sort_values("seg_len"))
    len_freq.to_csv(OUTDIR/"length_frequencies.tsv", sep="	", index=False)

    # OVERALL (across all lengths)
    all_rows = []
    # Any TE
    all_rows.append(fisher_table(df, "TOP5", "nearest_TE_any", DIST_THRESH, te_class="Any_TE"))
    # Per-class
    for cl, col in by_class.items():
        all_rows.append(fisher_table(df, "TOP5", col, DIST_THRESH, te_class=cl))
    overall = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    if not overall.empty:
        overall.to_csv(OUTDIR/"fisher_overall_TE.tsv", sep="	", index=False)

    # LENGTH-MATCHED (exact length equality)
    keep_lengths = set(len_freq.loc[len_freq["count"] >= MIN_LEN_FREQ, "seg_len"].astype(int).tolist())
    lm_rows = []
    for L in sorted(keep_lengths):
        sub = df[df["seg_len"] == L]
        lm_rows.append(fisher_table(sub, "TOP5", "nearest_TE_any", DIST_THRESH, te_class="Any_TE").assign(seg_len=L))
        for cl, col in by_class.items():
            lm_rows.append(fisher_table(sub, "TOP5", col, DIST_THRESH, te_class=cl).assign(seg_len=L))
    lm = pd.concat(lm_rows, ignore_index=True) if lm_rows else pd.DataFrame()
    if not lm.empty:
        # Order columns for readability
        col_order = ["te_class","seg_len","threshold_bp","TOP_in","TOP_out","BKG_in","BKG_out",
                     "TOP_rate","BKG_rate","rate_ratio","fisher_p","fdr_q"]
        lm = lm[col_order]
        lm.to_csv(OUTDIR/"fisher_by_length_TE.tsv", sep="	", index=False)

    # Summary
    with open(OUTDIR/"summary.txt", "w") as f:
        f.write("TMRCA near-TE enrichment with length-matched Fisher tests
")
        f.write(f"Top percentile = {TOP_PERCENTILE:.1f}%  (cutoff on mean_tmrca = {cutoff:.6g})
")
        f.write(f"Total rows: {df.shape[0]:,}
")
        f.write(f"Kept lengths (freq ≥ {MIN_LEN_FREQ}): {len(keep_lengths)}
")
        if keep_lengths:
            f.write(f"Min kept length: {min(keep_lengths)}  Max kept length: {max(keep_lengths)}
")
        f.write(f"Overall Fisher table: fisher_overall_TE.tsv
")
        f.write(f"Length-matched Fisher table: fisher_by_length_TE.tsv
")

    print("Done.")
    print("Outputs written under:", OUTDIR.resolve())


if __name__ == "__main__":
    main()
