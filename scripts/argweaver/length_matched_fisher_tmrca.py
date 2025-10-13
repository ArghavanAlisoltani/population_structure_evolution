#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TMRCA near-gene enrichment with LENGTH-MATCHED Fisher tests.

What this script does
---------------------
1) Reads a tab-delimited TMRCA file with EXACT columns:
   CHROM, start_tmrca, end_tmrca, mean_tmrca, lower_CI, upper_CI,
   overlap_genes, overlap_gene_n, nearest_genes

   • We compute segment length = (end_tmrca - start_tmrca) + 1
   • We split rows into TOP 5% vs BACKGROUND by mean_tmrca (configurable)
   • We need a distance to the nearest gene:
        - If a numeric column 'nearest_gene_dist' already exists, we use it.
        - Otherwise we compute it from an annotation workbook (see settings).

2) “Overall” Fisher exact tests (TOP vs BACKGROUND × within distance vs outside)
   for a set of gene-distance thresholds (e.g. 0, 1k, 2k, 5k, 10k, 50k).

3) NEW: LENGTH-MATCHED Fisher tests
   • Compute exact length frequencies, keep only lengths with freq ≥ MIN_LEN_FREQ.
   • For EACH kept length L, run the SAME Fisher set as in (2) using ONLY rows
     with seg_len == L (so Top and Background are matched on length).
   • Return one tidy table with one row per (length, distance).

Outputs (in OUTDIR)
-------------------
- fisher_overall.tsv                        (overall across all lengths)
- fisher_by_length.tsv                      (one row per length×distance)
- length_frequencies.tsv                    (all exact lengths and counts)
- summary.txt                               (cutoff and basic counts)

You can run this script on your existing file without touching your previous,
already-working pipeline. It only depends on pandas/scipy/statsmodels.

Author: you + ChatGPT
"""

from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =========================
# USER SETTINGS
# =========================
IN_TMRCA_TSV = "annotated_tmrca_4_GPT_13columns.tsv"  # your file with exact columns above

# If your TMRCA file does NOT already contain a numeric distance column named
# 'nearest_gene_dist', set ANNO_XLSX to an Excel workbook with columns:
#   contig_1ab, start1ab_pos, end1ab_pos  (+ a gene ID column is ignored here)
ANNO_XLSX   = "Aria_curated_annotation_1ab.xlsx"   # or None if 'nearest_gene_dist' exists
ANNO_SHEET  = None                                  # None = first sheet

# Define Top 5% on mean_tmrca (set DROP_ZERO_FOR_CUTOFF=True to ignore zeros)
TOP_PERCENTILE         = 95.0
DROP_ZERO_FOR_CUTOFF   = True

# Distance thresholds (in bp) to use in Fisher tests (0 means “overlaps a gene”)
DIST_THRESH = [0, 1_000, 2_000, 5_000, 10_000, 50_000]

# Minimum exact length frequency required to include a length in the
# length-matched Fisher panel
MIN_LEN_FREQ = 300

# Output directory
OUTDIR = Path("tmrca_len_matched_out")
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================
# Helpers
# =========================
def load_tmrca_table(path: str | Path) -> pd.DataFrame:
    need = ["CHROM","start_tmrca","end_tmrca","mean_tmrca",
            "lower_CI","upper_CI","overlap_genes","overlap_gene_n","nearest_genes"]
    df = pd.read_csv(path, sep="\t", dtype={"CHROM":str})
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

def load_annotation_distance(xlsx: str | Path, sheet=None) -> pd.DataFrame:
    """
    Returns a tidy gene table with columns:
      CHROM, GENE_START, GENE_END
    """
    if sheet is None:
        sheet = pd.ExcelFile(xlsx).sheet_names[0]
    ann = pd.read_excel(xlsx, sheet_name=sheet)
    req = ["contig_1ab","start1ab_pos","end1ab_pos"]
    miss = [c for c in req if c not in ann.columns]
    if miss:
        raise ValueError(f"Annotation missing columns: {miss}")
    g = ann[["contig_1ab","start1ab_pos","end1ab_pos"]].copy()
    g.columns = ["CHROM","GENE_START","GENE_END"]
    g["CHROM"]      = g["CHROM"].astype(str)
    g["GENE_START"] = pd.to_numeric(g["GENE_START"], errors="coerce").astype("Int64")
    g["GENE_END"]   = pd.to_numeric(g["GENE_END"], errors="coerce").astype("Int64")
    g = g.dropna(subset=["CHROM","GENE_START","GENE_END"]).copy()
    g["GENE_START"] = g["GENE_START"].astype(np.int64)
    g["GENE_END"]   = g["GENE_END"].astype(np.int64)
    return g.sort_values(["CHROM","GENE_START","GENE_END"]).reset_index(drop=True)

def compute_nearest_gene_distance(df: pd.DataFrame, genes: pd.DataFrame) -> pd.Series:
    """
    Compute edge-to-edge distance from each TMRCA segment to the nearest gene on the same scaffold.
    Distance = 0 if any overlap.
    """
    out = np.full(len(df), np.nan, dtype=float)
    by_chrom = df.groupby("CHROM", sort=False)
    for chrom, d in by_chrom:
        g = genes[genes["CHROM"] == chrom]
        if g.empty:
            continue
        # sort for a light sweep
        d = d.sort_values("start_tmrca")
        g = g.sort_values("GENE_START")
        gi = 0
        garr = g[["GENE_START","GENE_END"]].to_numpy(np.int64)
        for idx, (s,e) in d[["start_tmrca","end_tmrca"]].iterrows():
            # advance gene pointer to the first gene with start beyond (s - typical gene size cushion)
            while gi+1 < len(garr) and garr[gi+1,1] < s:
                gi += 1
            # check neighbors gi and gi+1
            candidates = []
            for k in (gi-1, gi, gi+1, gi+2):
                if 0 <= k < len(garr):
                    gs, ge = garr[k]
                    # overlap distance = 0, else min gap
                    if ge >= s and gs <= e:
                        candidates.append(0)
                    else:
                        gap = min(abs(s - ge), abs(gs - e))
                        candidates.append(gap)
            out[df.index.get_loc(idx)] = float(min(candidates)) if candidates else np.nan
    return pd.Series(out, index=df.index, name="nearest_gene_dist")

def fisher_table(sub: pd.DataFrame, flag_col: str,
                 dist_col: str, thresholds: list[int | float]) -> pd.DataFrame:
    """
    Count and Fisher test at multiple distance thresholds.
    Returns tidy table with counts/rates and p + BH–FDR q.
    """
    rows = []
    # use only rows with a finite distance
    D = sub[np.isfinite(sub[dist_col].values)].copy()
    for thr in thresholds:
        top_in  = int((D[flag_col] & (D[dist_col] <= thr)).sum())
        top_out = int((D[flag_col] & (D[dist_col] >  thr)).sum())
        bkg_in  = int((~D[flag_col] & (D[dist_col] <= thr)).sum())
        bkg_out = int((~D[flag_col] & (D[dist_col] >  thr)).sum())
        # Fisher two-sided
        _, p = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="two-sided")
        top_rate = top_in / max(1, (top_in + top_out))
        bkg_rate = bkg_in / max(1, (bkg_in + bkg_out))
        rr = (top_rate / bkg_rate) if bkg_rate > 0 else np.inf
        rows.append({
            "distance_bp": thr,
            "TOP_in": top_in, "TOP_out": top_out,
            "BKG_in": bkg_in, "BKG_out": bkg_out,
            "TOP_rate": top_rate, "BKG_rate": bkg_rate,
            "rate_ratio": rr, "fisher_p": p
        })
    out = pd.DataFrame(rows)
    out["fdr_q"] = multipletests(out["fisher_p"].values, method="fdr_bh")[1]
    return out

# =========================
# Main
# =========================
def main():
    df = load_tmrca_table(IN_TMRCA_TSV)

    # If a numeric 'nearest_gene_dist' already exists, use it; else compute
    if "nearest_gene_dist" in df.columns and pd.api.types.is_numeric_dtype(df["nearest_gene_dist"]):
        dist_col = "nearest_gene_dist"
    else:
        if ANNO_XLSX is None:
            raise ValueError(
                "No numeric 'nearest_gene_dist' column found in the input and ANNO_XLSX is None.\n"
                "Either add a numeric 'nearest_gene_dist' column or provide the annotation workbook "
                "so distances can be computed."
            )
        genes = load_annotation_distance(ANNO_XLSX, sheet=ANNO_SHEET)
        df["nearest_gene_dist"] = compute_nearest_gene_distance(df, genes)
        dist_col = "nearest_gene_dist"

    # Define Top5% vs Background on mean_tmrca
    vals = df["mean_tmrca"].replace([np.inf, -np.inf], np.nan)
    if DROP_ZERO_FOR_CUTOFF:
        vals = vals.mask(vals <= 0, np.nan)
    cutoff = float(np.nanpercentile(vals, TOP_PERCENTILE))
    df["TOP5"] = df["mean_tmrca"] >= cutoff

    # Save length frequencies (for transparency)
    len_freq = df["seg_len"].value_counts().sort_index()
    len_freq.rename_axis("seg_length_bp").rename("count").to_csv(
        OUTDIR / "length_frequencies.tsv", sep="\t"
    )

    # 1) OVERALL Fisher (all lengths pooled)
    overall = fisher_table(df, "TOP5", dist_col, DIST_THRESH)
    overall.to_csv(OUTDIR / "fisher_overall.tsv", sep="\t", index=False)

    # 2) LENGTH-MATCHED Fisher
    #    keep only exact lengths with >= MIN_LEN_FREQ
    keep_lengths = sorted(len_freq[len_freq >= MIN_LEN_FREQ].index.astype(int).tolist())
    all_rows = []
    all_p = []
    for L in keep_lengths:
        subL = df[df["seg_len"] == L].copy()
        if subL.empty:
            continue
        tbl = fisher_table(subL, "TOP5", dist_col, DIST_THRESH)
        tbl.insert(0, "seg_length_bp", int(L))
        tbl.insert(1, "N_rows", subL.shape[0])
        all_rows.append(tbl)
        all_p.extend(tbl["fisher_p"].tolist())

    if all_rows:
        bylen = pd.concat(all_rows, ignore_index=True)
        # Global BH–FDR across ALL length×distance tests
        bylen["fdr_q_global"] = multipletests(bylen["fisher_p"].values, method="fdr_bh")[1]
        bylen.to_csv(OUTDIR / "fisher_by_length.tsv", sep="\t", index=False)
    else:
        bylen = pd.DataFrame(columns=["seg_length_bp","N_rows","distance_bp","TOP_in","TOP_out",
                                      "BKG_in","BKG_out","TOP_rate","BKG_rate","rate_ratio",
                                      "fisher_p","fdr_q","fdr_q_global"])
        bylen.to_csv(OUTDIR / "fisher_by_length.tsv", sep="\t", index=False)

    # summary note
    with open(OUTDIR / "summary.txt", "w") as f:
        f.write(f"Top percentile = {TOP_PERCENTILE:.1f}%  (cutoff on mean_tmrca = {cutoff:.6g})\n")
        f.write(f"Total rows: {df.shape[0]:,}\n")
        f.write(f"Kept lengths (freq ≥ {MIN_LEN_FREQ}): {len(keep_lengths)}\n")
        if keep_lengths:
            f.write(f"Min kept length: {min(keep_lengths)}  Max kept length: {max(keep_lengths)}\n")
        f.write(f"Overall Fisher table: fisher_overall.tsv\n")
        f.write(f"Length-matched Fisher table: fisher_by_length.tsv\n")

    print("Done.")
    print("Outputs written under:", OUTDIR.resolve())

if __name__ == "__main__":
    main()
