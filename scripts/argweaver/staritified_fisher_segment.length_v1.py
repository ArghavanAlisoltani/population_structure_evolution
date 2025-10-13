#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Windowed TMRCA pipeline + NEW segment-length–stratified Fisher tests.

WHAT’S NEW vs. your last script
-------------------------------
1) Keeps ALL the windowed analysis you already had (weighted means per window,
   TOP5% on window means, near-gene Fisher on windows, Manhattan/line plots).

2) Adds a SECOND analysis that works at the **segment level** (raw ARGweaver
   intervals), exactly as you requested:
   • Compute segment length = end_tmrca - start_tmrca + 1
   • Define Top 5% **on segments** by mean_tmrca (you can choose to drop zeros)
   • Determine whether each segment overlaps ≥1 gene (uses `overlap_gene_n` if
     present; otherwise computes overlap from the same annotation workbook)
   • Tabulate the frequency of every exact segment length and keep only those
     lengths with count ≥ MIN_LEN_FREQ (default 300)
   • For EACH such length L, run a Fisher exact test (Top vs Background ×
     Overlap vs Not) and output a per-length results table with BH–FDR q values.

INPUT (tab-delimited) — EXACT column names expected (you already have these):
  CHROM  start_tmrca  end_tmrca  mean_tmrca  lower_CI  upper_CI
  overlap_genes  overlap_gene_n  nearest_genes   <-- used if present (overlap flag)

Annotation workbook (optional but recommended for distances/overlap if needed):
  Excel with columns: contig_1ab, start1ab_pos, end1ab_pos, and a gene ID col (default "NR.ID")

Outputs (written to OUTDIR):
  - tmrca_windows_win{WIN}_ov{OVERLAP}.tsv            (window table + flags)
  - fisher_near_gene_windows.tsv                      (window Fisher, various radii)
  - binom_enrichment.tsv                              (per-scaffold binomial on overlap=0)
  - seglen_hist.png, seglen_boxplot.png               (diagnostics on raw segments)
  - seglen_fisher_by_exact_length.tsv                 (NEW: per-length Fisher on segments)
  - manhattan_linear.png / manhattan_log.png          (windows)
  - line_with_ci_linear.png / line_with_ci_log.png    (windows)
  - README_analysis.txt

Tune everything in USER SETTINGS below.
"""

from pathlib import Path
import os, re, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import fisher_exact, binom
from statsmodels.stats.multitest import multipletests

# ------------------------
# USER SETTINGS
# ------------------------
# pwd all_tmrca/all_tmrca
IN_TMRCA_TSV = "annotated_tmrca_4_GPT_13columns.tsv"  # segments (columns listed above)
ANNO_XLSX    = None     # set to None to skip annotation
ANNO_SHEET   = None                                    # None => first sheet
ANNO_ID_COL  = "NR.ID"

# Windowing (uniform tiling; OVERLAP_BP = WIN_BP - STEP_BP)
WIN_BP       = 1_000_000
OVERLAP_BP   = 100_000
MIN_COVER_BP = 10_000           # drop windows with tiny covered bp

# Define Top 5% on WINDOW means:
TOP_PCT_WINDOWS = 95.0          # percentile of W_CENTER (after optional zero drop)
REMOVE_ZERO_FOR_WINDOW_CUTOFF = True

# Distances (bp) for near-gene Fisher on windows (0 = overlap)
DIST_THRESH = [0, 1_000, 2_000, 5_000, 10_000, 50_000]

# NEW: Segment-level Top 5% and per-length Fisher
TOP_PCT_SEGMENTS = 95.0         # percentile on mean_tmrca across raw segments
REMOVE_ZERO_FOR_SEG_CUTOFF = True
MIN_LEN_FREQ     = 300          # only test exact lengths with at least this many segments

# Plots
FIG_FONT   = "DejaVu Sans"
FIG_W, FIG_H = 14, 4
SCATTER_SIZE = 10
LABEL_SIZE   = 9
TITLE_SIZE   = 16
AXIS_LABEL   = 13
TICK_SIZE    = 10
ALPHA_BELOW_Q3 = 0.35
POINT_ALPHA    = 0.9
TOP_N_ANNOTATE = 10

OUTDIR = Path("tmrca_windowed_out")
OUTDIR.mkdir(parents=True, exist_ok=True)
rcParams["font.family"] = FIG_FONT

# ------------------------
# Helper functions
# ------------------------
def scaffold_sort_key(chrom: str):
    s = str(chrom)
    m = re.search(r"(\d+)([A-Za-z]?)$", s)
    if not m: return (10**9, 9, s)
    num = int(m.group(1)); let = m.group(2).lower()
    return (num, (ord(let) - ord('a')) if let else 9, s)

def windows_for_scaffold(xmin: int, xmax: int, win: int, overlap: int) -> np.ndarray:
    step = win - overlap
    if step <= 0:
        raise ValueError("OVERLAP_BP must be smaller than WIN_BP (step > 0).")
    starts = np.arange(xmin, xmax+1, step, dtype=np.int64)
    ends   = np.minimum(starts + win - 1, xmax)
    return np.vstack([starts, ends]).T

def weighted_window_stats(seg_s, seg_e, seg_mean, seg_lo, seg_hi, w_start, w_end):
    bp_cov = 0; num_m = num_lo = num_hi = 0.0; nseg = 0
    for s,e,m,lo,hi in zip(seg_s, seg_e, seg_mean, seg_lo, seg_hi):
        ov = max(0, min(e, w_end) - max(s, w_start) + 1)
        if ov > 0:
            bp_cov += ov
            num_m  += m  * ov
            num_lo += lo * ov
            num_hi += hi * ov
            nseg   += 1
    if bp_cov == 0:
        return 0, np.nan, np.nan, np.nan, 0
    return bp_cov, num_m/bp_cov, num_lo/bp_cov, num_hi/bp_cov, nseg

def read_annotation(xlsx, sheet=None, id_col="NR.ID"):
    if xlsx is None: return None
    if sheet is None: sheet = pd.ExcelFile(xlsx).sheet_names[0]
    ann = pd.read_excel(xlsx, sheet_name=sheet)
    req = ["contig_1ab","start1ab_pos","end1ab_pos"]
    miss = [c for c in req if c not in ann.columns]
    if miss: raise ValueError(f"Annotation missing columns: {miss}")
    if id_col not in ann.columns:
        ann[id_col] = [f"GENE_{i+1}" for i in range(len(ann))]
    out = ann[["contig_1ab","start1ab_pos","end1ab_pos", id_col]].copy()
    out.columns = ["CHROM","GENE_START","GENE_END","GENE_ID"]
    out["CHROM"]  = out["CHROM"].astype(str)
    out["GENE_START"] = pd.to_numeric(out["GENE_START"], errors="coerce").astype("Int64")
    out["GENE_END"]   = pd.to_numeric(out["GENE_END"],   errors="coerce").astype("Int64")
    out = out.dropna(subset=["CHROM","GENE_START","GENE_END"])
    out["GENE_START"] = out["GENE_START"].astype(np.int64)
    out["GENE_END"]   = out["GENE_END"].astype(np.int64)
    return out.sort_values(["CHROM","GENE_START","GENE_END"]).reset_index(drop=True)

def gene_overlap_and_distance(genes_df: pd.DataFrame, chrom: str, w_start: int, w_end: int):
    if genes_df is None: return (np.nan, "", np.nan)
    sub = genes_df[genes_df["CHROM"] == chrom]
    if sub.empty: return (0, "", np.nan)
    # overlap
    ov = sub[(sub["GENE_END"] >= w_start) & (sub["GENE_START"] <= w_end)]
    n_ov = int(ov.shape[0]); ids = ";".join(map(str, ov["GENE_ID"].tolist())) if n_ov>0 else ""
    if n_ov > 0: return (n_ov, ids, 0)
    # edge-to-edge distance
    left  = (w_start - sub["GENE_END"]).clip(lower=0)
    right = (sub["GENE_START"] - w_end).clip(lower=0)
    dist = int(pd.concat([left, right], axis=1).min(axis=1).min())
    return (0, "", dist)

def bh_fdr(p):
    p = np.asarray(p, float)
    mask = np.isfinite(p)
    q = np.full_like(p, np.nan, dtype=float)
    if mask.sum():
        q[mask] = multipletests(p[mask], method="fdr_bh")[1]
    return q

# ------------------------
# Load segments
# ------------------------
need = ["CHROM","start_tmrca","end_tmrca","mean_tmrca","lower_CI","upper_CI"]
seg = pd.read_csv(IN_TMRCA_TSV, sep="\t", dtype={"CHROM":str})
missing = [c for c in need if c not in seg.columns]
if missing:
    raise ValueError(f"Input missing columns: {missing}")

for c in need[1:]:
    seg[c] = pd.to_numeric(seg[c], errors="coerce")
seg = seg.dropna(subset=["CHROM","start_tmrca","end_tmrca","mean_tmrca"]).copy()
seg["start_tmrca"] = seg["start_tmrca"].astype(np.int64)
seg["end_tmrca"]   = seg["end_tmrca"].astype(np.int64)
seg["seg_len"]     = (seg["end_tmrca"] - seg["start_tmrca"] + 1).clip(lower=0)

# Prefer provided overlap flag if available
if "overlap_gene_n" in seg.columns:
    seg["SEG_OVERLAP_GENE"] = (pd.to_numeric(seg["overlap_gene_n"], errors="coerce").fillna(0) > 0)
else:
    seg["SEG_OVERLAP_GENE"] = np.nan  # will be filled if ANNO_XLSX set

# ------------------------
# Read annotation once
# ------------------------
genes = read_annotation(ANNO_XLSX, sheet=ANNO_SHEET, id_col=ANNO_ID_COL) if ANNO_XLSX else None

# If segment overlap flag missing, compute it quickly from genes
if genes is not None and seg["SEG_OVERLAP_GENE"].isna().any():
    msk = seg["SEG_OVERLAP_GENE"].isna()
    seg_tmp = seg.loc[msk, ["CHROM","start_tmrca","end_tmrca"]].copy()
    seg_tmp["SEG_OVERLAP_GENE"] = False
    for chrom, gsub in genes.groupby("CHROM", sort=False):
        ssub = seg_tmp[seg_tmp["CHROM"]==chrom]
        if ssub.empty: continue
        # simple interval sweep (naive; OK for one-off)
        idxs = ssub.index.values
        for i in idxs:
            a,b = int(ssub.at[i,"start_tmrca"]), int(ssub.at[i,"end_tmrca"])
            ov = gsub[(gsub["GENE_END"] >= a) & (gsub["GENE_START"] <= b)]
            if not ov.empty:
                seg_tmp.at[i,"SEG_OVERLAP_GENE"] = True
    seg.loc[msk, "SEG_OVERLAP_GENE"] = seg_tmp["SEG_OVERLAP_GENE"].values

# ------------------------
# Segment-length diagnostics (plots)
# ------------------------
def seglen_hist(df: pd.DataFrame, out_png: Path):
    plt.figure(figsize=(8,3.2))
    x = df["seg_len"].clip(lower=1).astype(float).values
    plt.hist(x, bins=60, log=False, color="#8da0cb")
    plt.xscale("log")
    plt.xlabel("Segment length (bp, log scale)")
    plt.ylabel("Count")
    plt.title("TMRCA segment lengths (log-scale histogram)")
    plt.tight_layout(); plt.savefig(out_png, dpi=200); plt.close()

seglen_hist(seg, OUTDIR / "seglen_hist.png")

# ------------------------
# Build windows & compute weighted means (WINDOW analysis)
# ------------------------
scaffolds = sorted(seg["CHROM"].unique(), key=scaffold_sort_key)
rows = []
for ch in scaffolds:
    d = seg[seg["CHROM"]==ch].sort_values(["start_tmrca","end_tmrca"])
    if d.empty: continue
    xmin, xmax = int(d["start_tmrca"].min()), int(d["end_tmrca"].max())
    win_arr = windows_for_scaffold(xmin, xmax, WIN_BP, OVERLAP_BP)
    ss = d["start_tmrca"].to_numpy(np.int64)
    ee = d["end_tmrca"].to_numpy(np.int64)
    mm = d["mean_tmrca"].to_numpy(float)
    ll = d["lower_CI"].to_numpy(float)
    uu = d["upper_CI"].to_numpy(float)

    for ws, we in win_arr:
        bp_cov, wmean, wlo, whi, nseg = weighted_window_stats(ss,ee,mm,ll,uu, ws, we)
        n_g, ids, ndist = gene_overlap_and_distance(genes, ch, int(ws), int(we))
        if bp_cov >= MIN_COVER_BP:
            rows.append({
                "CHROM": ch, "WIN_START": int(ws), "WIN_END": int(we),
                "WIN_MID": int((ws+we)//2),
                "bp_covered": int(bp_cov), "n_segments": int(nseg),
                "W_CENTER": wmean, "W_LO": wlo, "W_HI": whi,
                "NGENE_OVERLAP": n_g, "overlap_gene_ids": ids,
                "NEAREST_GENE_DIST": ndist
            })

win = pd.DataFrame(rows)
if win.empty:
    raise RuntimeError("No windows passed coverage threshold; adjust MIN_COVER_BP / window settings.")

# Top 5% cutoff on window means
vals_win = win["W_CENTER"].replace([np.inf,-np.inf], np.nan)
if REMOVE_ZERO_FOR_WINDOW_CUTOFF:
    vals_win = vals_win.mask(vals_win <= 0, np.nan)
cutoff_win = float(np.nanpercentile(vals_win, TOP_PCT_WINDOWS))
win["TOP5_WIN"] = win["W_CENTER"] >= cutoff_win
win["CUTOFF_95pct_WINDOW"] = cutoff_win

# ------------------------
# Near-gene Fisher on WINDOWS (TOP vs BKG) for DIST_THRESH
# ------------------------
fisher_rows = []
subw = win[pd.notna(win["NEAREST_GENE_DIST"])].copy()
for D in DIST_THRESH:
    top_in  = int((subw["TOP5_WIN"] & (subw["NEAREST_GENE_DIST"] <= D)).sum())
    top_out = int((subw["TOP5_WIN"] & (subw["NEAREST_GENE_DIST"] >  D)).sum())
    bkg_in  = int(((~subw["TOP5_WIN"]) & (subw["NEAREST_GENE_DIST"] <= D)).sum())
    bkg_out = int(((~subw["TOP5_WIN"]) & (subw["NEAREST_GENE_DIST"] >  D)).sum())
    _, pval = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="two-sided")
    top_rate = top_in / max(1, (top_in + top_out))
    bkg_rate = bkg_in / max(1, (bkg_in + bkg_out))
    rr = (top_rate / bkg_rate) if bkg_rate>0 else math.inf
    fisher_rows.append(dict(distance_bp=D, TOP_in=top_in, TOP_out=top_out,
                            BKG_in=bkg_in, BKG_out=bkg_out,
                            TOP_rate=top_rate, BKG_rate=bkg_rate,
                            rate_ratio=rr, fisher_p=pval))
fisher_win = pd.DataFrame(fisher_rows)
fisher_win["fdr_q"] = bh_fdr(fisher_win["fisher_p"].to_numpy())
fisher_win.to_csv(OUTDIR / "fisher_near_gene_windows.tsv", sep="\t", index=False)

# ------------------------
# Per-scaffold binomial (windows; overlap=0)
# ------------------------
def per_scaffold_binomial(win_df: pd.DataFrame, flag_col: str, near_col: str) -> pd.DataFrame:
    rows = []
    for chrom, d in win_df.groupby("CHROM", sort=False):
        top = d[d[flag_col]]; bkg = d[~d[flag_col]]
        if d.empty or bkg.empty: continue
        p0 = float((bkg[near_col]).mean())  # background rate on that scaffold
        K  = int((top[near_col]).sum())
        n  = int(top.shape[0])
        if n == 0: continue
        p = 1.0 - binom.cdf(K-1, n, p0) if p0 > 0 else (1.0 if K==0 else 0.0)
        rows.append({"CHROM": chrom, "p_binom": p, "n_top": n, "k_top": K, "p0_bkg": p0})
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["CHROM","p_binom","q_binom","n_top","k_top","p0_bkg"])
    out["q_binom"] = bh_fdr(out["p_binom"].to_numpy())
    return out.sort_values("q_binom")

win["OVERLAP_BOOL"] = (win["NGENE_OVERLAP"] > 0)
binom_df = per_scaffold_binomial(win, flag_col="TOP5_WIN", near_col="OVERLAP_BOOL")
binom_df.to_csv(OUTDIR / "binom_enrichment.tsv", sep="\t", index=False)

# ------------------------
# SAVE window table
# ------------------------
win_out = OUTDIR / f"tmrca_windows_win{WIN_BP}_ov{OVERLAP_BP}.tsv"
win.to_csv(win_out, sep="\t", index=False)

# ------------------------
# NEW: Segment-length–stratified Fisher (SEGMENT analysis)
# ------------------------
# Define Top 5% on segments by mean_tmrca (optionally drop zeros before percentile)
vals_seg = seg["mean_tmrca"].replace([np.inf,-np.inf], np.nan)
if REMOVE_ZERO_FOR_SEG_CUTOFF:
    vals_seg = vals_seg.mask(vals_seg <= 0, np.nan)
cutoff_seg = float(np.nanpercentile(vals_seg, TOP_PCT_SEGMENTS))
seg["TOP5_SEG"] = seg["mean_tmrca"] >= cutoff_seg

# Ensure an overlap flag for segments exists (boolean)
if "SEG_OVERLAP_GENE" not in seg.columns:
    seg["SEG_OVERLAP_GENE"] = False  # should not happen (we set it above)

# Count exact length frequencies and keep only frequent ones
len_counts = seg["seg_len"].value_counts()
keep_lengths = set(len_counts[len_counts >= MIN_LEN_FREQ].index.tolist())

rows_len = []
for L in sorted(keep_lengths):
    subL = seg[seg["seg_len"] == L]
    if subL.empty: continue
    top_in  = int(((subL["TOP5_SEG"]) &  (subL["SEG_OVERLAP_GENE"])).sum())
    top_out = int(((subL["TOP5_SEG"]) & ~(subL["SEG_OVERLAP_GENE"])).sum())
    bkg_in  = int((~(subL["TOP5_SEG"]) &  (subL["SEG_OVERLAP_GENE"])).sum())
    bkg_out = int((~(subL["TOP5_SEG"]) & ~(subL["SEG_OVERLAP_GENE"])).sum())

    # Guard tiny tables (shouldn’t happen with MIN_LEN_FREQ, but be safe)
    if (top_in + top_out == 0) or (bkg_in + bkg_out == 0):
        p_two = np.nan; p_one = np.nan
        rr = np.nan
    else:
        # two-sided Fisher
        _, p_two = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="two-sided")
        # one-sided (enrichment of overlaps in TOP)
        _, p_one = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="greater")
        # overlap rate ratio
        rt = top_in / max(1, (top_in + top_out))
        rb = bkg_in / max(1, (bkg_in + bkg_out))
        rr = (rt / rb) if rb > 0 else math.inf

    rows_len.append({
        "seg_length_bp": int(L),
        "N_segments": int(subL.shape[0]),
        "TOP_in": top_in, "TOP_out": top_out,
        "BKG_in": bkg_in, "BKG_out": bkg_out,
        "top_overlap_rate": (top_in / max(1, (top_in+top_out))),
        "bkg_overlap_rate": (bkg_in / max(1, (bkg_in+bkg_out))),
        "rate_ratio": rr,
        "fisher_two_sided_p": p_two,
        "fisher_one_sided_p": p_one
    })

seglen_fisher = pd.DataFrame(rows_len)
if not seglen_fisher.empty:
    seglen_fisher["q_two_sided"] = bh_fdr(seglen_fisher["fisher_two_sided_p"].to_numpy())
    seglen_fisher["q_one_sided"] = bh_fdr(seglen_fisher["fisher_one_sided_p"].to_numpy())
seglen_fisher = seglen_fisher.sort_values("seg_length_bp")
seglen_fisher.to_csv(OUTDIR / "seglen_fisher_by_exact_length.tsv", sep="\t", index=False)

# OPTIONAL: a quick boxplot to show segment lengths vs TOP/BKG (log y)
def seglen_top_box(df: pd.DataFrame, cutoff: float, out_png: Path):
    d = df.copy()
    d["group"] = np.where(d["mean_tmrca"] >= cutoff, "TOP5%", "BACKGROUND")
    plt.figure(figsize=(5,6))
    d.boxplot(column="seg_len", by="group")
    plt.yscale("log")
    plt.title("")
    plt.suptitle("")
    plt.ylabel("Segment length (bp, log)")
    plt.tight_layout(); plt.savefig(out_png, dpi=200); plt.close()

seglen_top_box(seg, cutoff_seg, OUTDIR / "seglen_boxplot.png")

# ------------------------
# Plotting (WINDOWS): Manhattan + Line w/ CI
# ------------------------
def build_concat_positions(win_df: pd.DataFrame):
    chroms = sorted(win_df["CHROM"].unique(), key=scaffold_sort_key)
    offsets, xticks, xtlabs = {}, [], []
    cursor = 0
    for ch in chroms:
        sub = win_df[win_df["CHROM"]==ch]
        if sub.empty: continue
        max_mid = int(sub["WIN_MID"].max())
        offsets[ch] = cursor
        xticks.append(cursor + max_mid//2)
        xtlabs.append(ch)
        cursor += max_mid + 1
    return offsets, list(zip(xtlabs, xticks))

offsets, chrom_ticks = build_concat_positions(win)
q75 = win["W_CENTER"].quantile(0.75)
alphas = np.where(win["W_CENTER"] < q75, ALPHA_BELOW_Q3, POINT_ALPHA)

def concat_x(df):
    return df["WIN_MID"] + df["CHROM"].map(offsets).fillna(0).astype(np.int64)

def manhattan_scatter(df, ycol, log10=False, out_png="scatter.png", title=""):
    cmap = {ch: ("#4C78A8" if (i % 2 == 0) else "#A0A0A0")
            for i, ch in enumerate(sorted(df["CHROM"].unique(), key=scaffold_sort_key))}
    colors = df["CHROM"].map(cmap).values
    X = concat_x(df); Y = df[ycol].astype(float)
    topN = df.nlargest(TOP_N_ANNOTATE, ycol).copy()
    topN_x = concat_x(topN); topN_y = topN[ycol].astype(float)

    fig, ax = plt.subplots(figsize=(18,5))
    ax.scatter(X, Y, s=SCATTER_SIZE, c=colors, alpha=alphas, linewidths=0)
    ax.scatter(topN_x, topN_y, s=SCATTER_SIZE*2.2, c="red", alpha=1.0, linewidths=0, zorder=3)
    for _, r in topN.iterrows():
        label = f"{r['CHROM']}:{int(r['WIN_START'])}-{int(r['WIN_END'])}"
        ax.text(int(r['WIN_MID'])+offsets[r['CHROM']], float(r[ycol])*1.02,
                label, fontsize=LABEL_SIZE-2, color="red", rotation=35, ha="left", va="bottom")
    ax.set_xlabel("Genome (scaffolds concatenated)", fontsize=AXIS_LABEL)
    ax.set_ylabel(("W_CENTER (generations, log10)" if log10 else "W_CENTER (generations)"),
                  fontsize=AXIS_LABEL)
    ax.set_xticks([t for _,t in chrom_ticks])
    ax.set_xticklabels([c for c,_ in chrom_ticks], rotation=90, fontsize=TICK_SIZE)
    if log10: ax.set_yscale("log")
    ax.set_title(title, fontsize=TITLE_SIZE)
    fig.tight_layout(); fig.savefig(out_png, dpi=220); plt.close(fig)

manhattan_scatter(
    win, "W_CENTER", log10=False,
    out_png=OUTDIR / f"manhattan_linear.png",
    title=f"Windowed TMRCA Manhattan (win={WIN_BP:,}, overlap={OVERLAP_BP:,})"
)
manhattan_scatter(
    win, "W_CENTER", log10=True,
    out_png=OUTDIR / f"manhattan_log.png",
    title=f"Windowed TMRCA Manhattan (log10) (win={WIN_BP:,}, overlap={OVERLAP_BP:,})"
)

def line_with_ci(df, out_png, title, use_log=False):
    d = df.sort_values(["CHROM","WIN_START"]).copy()
    d["x"] = concat_x(d)
    y  = d["W_CENTER"].astype(float).values
    lo = d["W_LO"].astype(float).values
    hi = d["W_HI"].astype(float).values
    if use_log:
        y  = np.log10(np.clip(y, 1e-9, None))
        lo = np.log10(np.clip(lo, 1e-9, None))
        hi = np.log10(np.clip(hi, 1e-9, None))
        ylabel = "W_CENTER (log10 generations)"
    else:
        ylabel = "W_CENTER (generations)"
    fig, ax = plt.subplots(figsize=(20,5))
    ax.fill_between(d["x"], lo, hi, alpha=0.15, step="pre", color="#999999")
    ax.plot(d["x"], y, lw=0.9, color="black")
    ax.set_xlabel("Genome (scaffolds concatenated)", fontsize=AXIS_LABEL)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL)
    ax.set_xticks([t for _,t in chrom_ticks])
    ax.set_xticklabels([c for c,_ in chrom_ticks], rotation=90, fontsize=TICK_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE)
    fig.tight_layout(); fig.savefig(out_png, dpi=220); plt.close(fig)

line_with_ci(win, OUTDIR / "line_with_ci_linear.png",
             f"Windowed mean TMRCA with CI (win={WIN_BP:,}, overlap={OVERLAP_BP:,})", use_log=False)
line_with_ci(win, OUTDIR / "line_with_ci_log.png",
             f"Windowed mean TMRCA with CI (log10) (win={WIN_BP:,}, overlap={OVERLAP_BP:,})", use_log=True)

# ------------------------
# README
# ------------------------
with open(OUTDIR / "README_analysis.txt", "w") as fh:
    fh.write(
        "All enrichment and Manhattan/line plots use coverage-weighted WINDOW means\n"
        f"(WIN_BP={WIN_BP:,}, OVERLAP_BP={OVERLAP_BP:,}, MIN_COVER_BP={MIN_COVER_BP:,}).\n"
        f"Window Top5 cutoff (95th percentile) = {cutoff_win:.6g}\n"
        f"Segment Top5 cutoff (95th percentile) = {cutoff_seg:.6g}\n"
        f"Segment-length Fisher ran for exact lengths with freq ≥ {MIN_LEN_FREQ}.\n"
        "Window Fisher tables report enrichment of near-gene windows at distances D={DIST_THRESH}.\n"
        "Per-scaffold binomial compares Top vs Background overlap rates within each scaffold.\n"
    )

print("\nDone.")
print("Outputs:", OUTDIR.resolve())
print("  - Windows TSV:", win_out.name)
print("  - Window Fisher:", "fisher_near_gene_windows.tsv")
print("  - Per-scaffold binomial:", "binom_enrichment.tsv")
print("  - Segment-length Fisher:", "seglen_fisher_by_exact_length.tsv")
print("  - Plots: manhattan_*.png, line_with_ci_*.png, seglen_*.png")
