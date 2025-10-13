#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Length-matched Fisher exact tests for TMRCA segments vs. Transposable Elements (TEs)

YOU ASKED:
- Use the ORIGINAL TMRCA segments (no windowing).
- Split segments into TOP 5% vs BACKGROUND by mean_tmrca.
- For EACH segment length with frequency ≥ MIN_LEN_FREQ, run Fisher tests
  comparing TOP vs BACKGROUND on proximity to TEs:
    (a) “Any TE” at multiple distance thresholds.
    (b) Per TE CLASS at multiple distance thresholds.

INPUTS
------
1) TMRCA segments (tab-delimited) with EXACT columns:
   CHROM  start_tmrca  end_tmrca  mean_tmrca  lower_CI  upper_CI
   overlap_genes  overlap_gene_n  nearest_genes
   (Other columns are ignored.)

2) TE annotation file: `col8_readable.TEanno.cds.scaf1ab.tsv`
   with headers exactly:
   seqid  source  sequence_ontology  start  end  score  strand

   - We treat `seqid` as chromosome.
   - We treat `sequence_ontology` as the TE CLASS label (e.g., Gypsy_LTR_retrotransposon).
   - Coordinates are assumed 1-based inclusive (as in your example).

OUTPUTS (in OUTDIR)
-------------------
- fisher_anyTE_overall.tsv
- fisher_anyTE_by_length.tsv
- fisher_by_class_overall.tsv
- fisher_by_class_by_length.tsv
- length_frequencies.tsv
- summary.txt

NOTES
-----
• “Distance” = 0 if a segment overlaps ≥1 TE (edge-to-edge overlap).
  Else it is the minimum gap (in bp) to the nearest TE (or nearest TE of a class).
• Fisher tables have one row per distance threshold (e.g., 0, 1k, 2k, 5k, 10k, 50k).
• Length-matched tables subset to a single exact segment length at a time,
  but only for lengths with frequency ≥ MIN_LEN_FREQ (to keep tests stable).
• FDR (BH) is computed within each table.

Author: you + ChatGPT
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from bisect import bisect_left
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =========================
# USER SETTINGS
# =========================
TMRCA_TSV   = "annotated_tmrca_4_GPT_13columns.tsv"        # original segments (not windows)
TE_TSV      = "col8_readable.TEanno.cds.scaf1ab.tsv"       # your TE file (headers given)

TOP_PERCENTILE       = 95.0        # Top 5% by mean_tmrca
DROP_ZERO_FOR_CUTOFF = True        # ignore mean_tmrca<=0 when computing the cutoff

DIST_THRESH          = [0, 1_000, 2_000, 5_000, 10_000, 50_000]
MIN_LEN_FREQ         = 300         # only run length-matched Fisher for lengths with ≥ this many segments

OUTDIR = Path("tmrca_TE_length_matched_out")
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================
# IO HELPERS
# =========================
def load_tmrca(path: str | Path) -> pd.DataFrame:
    need = ["CHROM","start_tmrca","end_tmrca","mean_tmrca",
            "lower_CI","upper_CI","overlap_genes","overlap_gene_n","nearest_genes"]
    df = pd.read_csv(path, sep="\t", dtype={"CHROM":str})
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"TMRCA file missing columns: {miss}")
    for c in ["start_tmrca","end_tmrca","mean_tmrca"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["CHROM","start_tmrca","end_tmrca","mean_tmrca"]).copy()
    df["start_tmrca"] = df["start_tmrca"].astype(np.int64)
    df["end_tmrca"]   = df["end_tmrca"].astype(np.int64)
    df["seg_len"]     = (df["end_tmrca"] - df["start_tmrca"] + 1).clip(lower=0)
    return df

def load_te(path: str | Path) -> pd.DataFrame:
    # seqid  source  sequence_ontology  start  end  score  strand
    te = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    need = ["seqid","source","sequence_ontology","start","end","score","strand"]
    miss = [c for c in need if c not in te.columns]
    if miss:
        raise ValueError(f"TE file missing columns: {miss}")
    te = te[["seqid","sequence_ontology","start","end"]].copy()
    te.columns = ["CHROM","TE_CLASS","START","END"]
    te["CHROM"] = te["CHROM"].astype(str)
    te["START"] = pd.to_numeric(te["START"], errors="coerce").astype("Int64")
    te["END"]   = pd.to_numeric(te["END"],   errors="coerce").astype("Int64")
    te = te.dropna(subset=["CHROM","START","END"]).copy()
    te["START"] = te["START"].astype(np.int64)
    te["END"]   = te["END"].astype(np.int64)
    te = te[te["END"] >= te["START"]].copy()
    return te.sort_values(["CHROM","START","END"]).reset_index(drop=True)

# =========================
# DISTANCE KERNELS
# =========================
def _nearest_gap_to_sorted_features(seg_s:int, seg_e:int,
                                    f_starts:np.ndarray, f_ends:np.ndarray) -> int:
    """
    Given a segment [seg_s, seg_e] and arrays of feature starts/ends (both sorted by start),
    return 0 if any overlap, else the minimal edge gap in bp.
    """
    # position to insert seg_s among starts
    i = bisect_left(f_starts, seg_s)

    # check a small neighborhood for overlap (i-1, i, i+1)
    best = np.inf
    for k in (i-2, i-1, i, i+1, i+2):
        if 0 <= k < f_starts.size:
            fs, fe = f_starts[k], f_ends[k]
            if fe >= seg_s and fs <= seg_e:
                return 0
            # compute gap
            gap = min(abs(seg_s - fe), abs(fs - seg_e))
            if gap < best: best = gap
    return int(best if np.isfinite(best) else np.inf)

def compute_nearest_distances_any_and_by_class(tm: pd.DataFrame,
                                               te: pd.DataFrame) -> tuple[pd.Series, dict[str, pd.Series]]:
    """
    For each TMRCA segment, compute:
      - distance to nearest TE of ANY class
      - distance to nearest TE of EACH class (returned as dict of Series)

    Returns:
      dist_any  : pd.Series (length = len(tm))
      dist_by_c : dict {class_name: pd.Series}
    """
    dist_any = np.full(len(tm), np.nan, dtype=float)
    dist_by_class: dict[str, np.ndarray] = {}

    # map index for assignment
    idx_map = {i: pos for pos, i in enumerate(tm.index)}

    for chrom, dfc in tm.groupby("CHROM", sort=False):
        te_c = te[te["CHROM"] == chrom]
        if te_c.empty:
            continue

        # ANY TE arrays
        arr_starts_any = te_c["START"].to_numpy(np.int64)
        arr_ends_any   = te_c["END"].to_numpy(np.int64)

        # class-specific arrays (pre-sorted by START already)
        class_groups = {cls: sub[["START","END"]].to_numpy(np.int64)
                        for cls, sub in te_c.groupby("TE_CLASS", sort=False)}

        for i, row in dfc[["start_tmrca","end_tmrca"]].iterrows():
            s, e = int(row["start_tmrca"]), int(row["end_tmrca"])
            # ANY
            d_any = _nearest_gap_to_sorted_features(s, e, arr_starts_any, arr_ends_any)
            dist_any[idx_map[i]] = d_any

            # BY CLASS
            for cls, arr in class_groups.items():
                if cls not in dist_by_class:
                    dist_by_class[cls] = np.full(len(tm), np.nan, dtype=float)
                d_cls = _nearest_gap_to_sorted_features(s, e, arr[:,0], arr[:,1])
                dist_by_class[cls][idx_map[i]] = d_cls

    # wrap into pandas objects
    dist_any_s = pd.Series(dist_any, index=tm.index, name="nearest_TE_any_bp")
    dist_by_class_s = {cls: pd.Series(vals, index=tm.index, name=f"nearest_TE_{cls}_bp")
                       for cls, vals in dist_by_class.items()}
    return dist_any_s, dist_by_class_s

# =========================
# FISHER HELPERS
# =========================
def fisher_table(mask_top: pd.Series, dist: pd.Series,
                 thresholds: list[int|float]) -> pd.DataFrame:
    """
    Build multi-threshold Fisher table (TOP vs BACKGROUND × within vs outside threshold).
    """
    D = dist.copy()
    ok = np.isfinite(D.values)
    mask_top = mask_top & ok
    mask_bkg = (~mask_top) & ok

    rows=[]
    for thr in thresholds:
        top_in  = int((mask_top & (D <= thr)).sum())
        top_out = int((mask_top & (D >  thr)).sum())
        bkg_in  = int((mask_bkg & (D <= thr)).sum())
        bkg_out = int((mask_bkg & (D >  thr)).sum())

        _, p = fisher_exact([[top_in, top_out], [bkg_in, bkg_out]], alternative="two-sided")
        top_rate = top_in / max(1, top_in + top_out)
        bkg_rate = bkg_in / max(1, bkg_in + bkg_out)
        rr = (top_rate / bkg_rate) if bkg_rate > 0 else np.inf

        rows.append(dict(distance_bp=thr,
                         TOP_in=top_in, TOP_out=top_out,
                         BKG_in=bkg_in, BKG_out=bkg_out,
                         TOP_rate=top_rate, BKG_rate=bkg_rate,
                         rate_ratio=rr, fisher_p=p))
    out = pd.DataFrame(rows)
    out["fdr_q"] = multipletests(out["fisher_p"].values, method="fdr_bh")[1]
    return out

# =========================
# MAIN
# =========================
def main():
    tm = load_tmrca(TMRCA_TSV)
    te = load_te(TE_TSV)

    # --- define TOP5 vs BACKGROUND on mean_tmrca
    vals = tm["mean_tmrca"].replace([np.inf, -np.inf], np.nan)
    if DROP_ZERO_FOR_CUTOFF:
        vals = vals.mask(vals <= 0, np.nan)
    cutoff = float(np.nanpercentile(vals, TOP_PERCENTILE))
    tm["TOP5"] = tm["mean_tmrca"] >= cutoff

    # --- compute nearest distances (ANY TE + per-class)
    dist_any, dist_by_class = compute_nearest_distances_any_and_by_class(tm, te)
    tm["nearest_TE_any_bp"] = dist_any

    # --- save length frequencies (for transparency)
    len_freq = tm["seg_len"].value_counts().sort_index()
    len_freq.rename_axis("seg_length_bp").rename("count").to_csv(
        OUTDIR / "length_frequencies.tsv", sep="\t"
    )

    # ======================
    # ANY-TE Fisher (overall)
    # ======================
    any_overall = fisher_table(tm["TOP5"], tm["nearest_TE_any_bp"], DIST_THRESH)
    any_overall.to_csv(OUTDIR / "fisher_anyTE_overall.tsv", sep="\t", index=False)

    # ======================
    # ANY-TE Fisher (by exact segment length, freq ≥ MIN_LEN_FREQ)
    # ======================
    keep_lengths = sorted(len_freq[len_freq >= MIN_LEN_FREQ].index.astype(int).tolist())
    rows=[]
    for L in keep_lengths:
        idx = (tm["seg_len"] == L)
        if idx.sum() == 0:
            continue
        tbl = fisher_table(tm.loc[idx, "TOP5"], tm.loc[idx, "nearest_TE_any_bp"], DIST_THRESH)
        tbl.insert(0, "seg_length_bp", int(L))
        tbl.insert(1, "N_rows", int(idx.sum()))
        rows.append(tbl)
    any_bylen = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["seg_length_bp","N_rows","distance_bp","TOP_in","TOP_out",
                 "BKG_in","BKG_out","TOP_rate","BKG_rate","rate_ratio","fisher_p","fdr_q"])
    any_bylen.to_csv(OUTDIR / "fisher_anyTE_by_length.tsv", sep="\t", index=False)

    # ======================
    # BY-CLASS Fisher (overall)
    # ======================
    # Build one table per TE CLASS: use class-specific nearest distances
    class_rows = []
    for cls, dist_s in dist_by_class.items():
        tbl = fisher_table(tm["TOP5"], dist_s, DIST_THRESH)
        tbl.insert(0, "TE_class", cls)
        class_rows.append(tbl)
    byclass_overall = pd.concat(class_rows, ignore_index=True) if class_rows else pd.DataFrame(
        columns=["TE_class","distance_bp","TOP_in","TOP_out","BKG_in","BKG_out",
                 "TOP_rate","BKG_rate","rate_ratio","fisher_p","fdr_q"])
    # global FDR across classes×thresholds
    if not byclass_overall.empty:
        byclass_overall["fdr_q_global"] = multipletests(byclass_overall["fisher_p"].values, method="fdr_bh")[1]
    byclass_overall.to_csv(OUTDIR / "fisher_by_class_overall.tsv", sep="\t", index=False)

    # ======================
    # BY-CLASS Fisher (by exact segment length, freq ≥ MIN_LEN_FREQ)
    # ======================
    rows=[]
    for L in keep_lengths:
        idx_len = (tm["seg_len"] == L)
        if idx_len.sum() == 0:
            continue
        for cls, dist_s in dist_by_class.items():
            tbl = fisher_table(tm.loc[idx_len, "TOP5"], dist_s.loc[idx_len], DIST_THRESH)
            tbl.insert(0, "seg_length_bp", int(L))
            tbl.insert(1, "N_rows", int(idx_len.sum()))
            tbl.insert(2, "TE_class", cls)
            rows.append(tbl)
    byclass_bylen = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["seg_length_bp","N_rows","TE_class","distance_bp","TOP_in","TOP_out",
                 "BKG_in","BKG_out","TOP_rate","BKG_rate","rate_ratio","fisher_p","fdr_q"])
    if not byclass_bylen.empty:
        byclass_bylen["fdr_q_global"] = multipletests(byclass_bylen["fisher_p"].values, method="fdr_bh")[1]
    byclass_bylen.to_csv(OUTDIR / "fisher_by_class_by_length.tsv", sep="\t", index=False)

    # summary
    with open(OUTDIR / "summary.txt","w") as fh:
        fh.write(f"Top percentile: {TOP_PERCENTILE:.1f}% (cutoff on mean_tmrca = {cutoff:.6g})\n")
        fh.write(f"Total segments: {tm.shape[0]:,}\n")
        fh.write(f"Kept lengths for length-matched tests (freq ≥ {MIN_LEN_FREQ}): {len(keep_lengths)}\n")
        if keep_lengths:
            fh.write(f"Smallest kept length: {min(keep_lengths)} bp; Largest kept length: {max(keep_lengths)} bp\n")
        fh.write(f"Any-TE Fisher: fisher_anyTE_overall.tsv, fisher_anyTE_by_length.tsv\n")
        fh.write(f"Per-class Fisher: fisher_by_class_overall.tsv, fisher_by_class_by_length.tsv\n")

    print("Done. Outputs in:", OUTDIR.resolve())

if __name__ == "__main__":
    main()
