#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Windowed TMRCA pipeline (segment-length aware)

INPUT (tab-delimited) — EXACT column names expected:
  CHROM   start_tmrca  end_tmrca  mean_tmrca  lower_CI  upper_CI
  overlap_genes  overlap_gene_n  nearest_genes     <-- not used here, kept for compatibility

PLUS (recommended) an annotation workbook to compute overlaps/distances to genes:
  Excel sheet with columns: contig_1ab, start1ab_pos, end1ab_pos, and a gene ID column (default "NR.ID")

What this script does (all metrics are computed on FIXED, OVERLAPPING WINDOWS — not on raw segments):
  1) Slide a window of size WIN_BP with step STEP_BP along each scaffold.
  2) For each window, compute LENGTH-WEIGHTED means from all TMRCA segments that overlap the window:
        - w_mean_tmrca  (weighted by bp overlap)
        - w_lower_CI    (weighted)
        - w_upper_CI    (weighted)
        - bp_covered (total bp overlapped by segments inside the window)
        - n_segments  (# overlapping segments)
  3) Compute gene annotation for windows (if ANNO_XLSX provided):
        - n_overlap_genes, overlap_gene_ids
        - nearest_gene_dist (0 if any overlap, otherwise min distance)
  4) Compute the 95th percentile of w_mean_tmrca across windows and flag the TOP5% windows.
     Save this FLAG in the TSV so the split is traceable.
  5) Fisher-style enrichment on distance-to-gene for TOP vs BACKGROUND windows, for thresholds:
        DIST_THRESH = [0, 1k, 2k, 5k, 10k, 50k]  (0 counts overlaps)
  6) Visualizations (scaffolds concatenated on x):
        - Manhattan-style scatter of w_mean_tmrca (linear and log10 scales)
        - Line plot across genome with CI band (uses w_lower_CI/w_upper_CI)
        - Top-N windows highlighted in red and annotated
  7) Write all derived windows to a TSV.

Customize everything in the USER SETTINGS block.
"""

# =========================
# USER SETTINGS
# =========================
IN_TMRCA_TSV = "annotated_tmrca_4_GPT_13columns.tsv"   # your segments file (exact cols listed above)
ANNO_XLSX    = "Aria_curated_annotation_both_1ab_1.xlsx" # set to None if you don't want gene distances
ANNO_SHEET   = None                                     # None -> first sheet
ANNO_ID_COL  = "NR.ID"                                  # column to use as gene ID

# Windowing
WIN_BP   = 1_000_000     # window size in bp
STEP_BP  = 100_000       # step between windows (overlap = WIN_BP - STEP_BP) — must be >0

# Top group definition on WINDOW means
TOP_PCT  = 0.05          # Top 5% of windows by w_mean_tmrca

# Distance thresholds (in bp) for Fisher comparisons (TOP vs BACKGROUND, on windows)
DIST_THRESH = [0, 1_000, 2_000, 5_000, 10_000, 50_000]

# Threshold to define "old" bp for optional comparisons (not used in Fisher here)
THR_HI = 50_000

# Plot controls
FIG_FONT = "DejaVu Sans"
SCATTER_SIZE = 10
LABEL_SIZE   = 9
TITLE_SIZE   = 16
AXIS_LABEL   = 13
TICK_SIZE    = 10
ALPHA_LOW    = 0.35      # alpha for windows below Q3 (optional dimming)
POINT_ALPHA  = 0.9
TOP_N_ANNOTATE = 10      # top-N windows to color red and annotate in scatter

# Output folder
OUTDIR = "tmrca_windowed_out"

# =========================
# IMPORTS
# =========================
import os, re, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from matplotlib import rcParams

# Fonts
rcParams["font.family"] = FIG_FONT

# =========================
# HELPERS
# =========================
def scaffold_sort_key(chrom: str):
    """Order scaffolds: 1a, 1b, 2, 3, ..., 10, 11, ...  (numbers first, letter suffix a<b<none)"""
    s = str(chrom)
    m = re.search(r"(\d+)([A-Za-z]?)$", s)
    if not m:
        return (10**9, 9, s)
    num = int(m.group(1))
    let = m.group(2).lower()
    let_rank = (ord(let) - ord('a')) if let else 9
    return (num, let_rank, s)

def windows_for_scaffold(xmin: int, xmax: int, win: int, step: int) -> np.ndarray:
    if step <= 0:
        raise ValueError("STEP_BP must be > 0 (overlap must be smaller than window size).")
    starts = np.arange(xmin, xmax+1, step, dtype=np.int64)
    ends   = np.minimum(starts + win - 1, xmax)
    return np.vstack([starts, ends]).T

def weighted_window_stats(seg_s, seg_e, seg_mean, seg_lo, seg_hi, w_start, w_end):
    """Compute length-weighted mean over the overlap between [w_start,w_end] and each segment.
       Returns (bp_covered, w_mean, w_lo, w_hi, n_segments).
    """
    bp_cov = 0
    num_m  = 0.0
    num_lo = 0.0
    num_hi = 0.0
    nseg   = 0
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
    return (bp_cov,
            num_m / bp_cov,
            num_lo / bp_cov,
            num_hi / bp_cov,
            nseg)

def read_annotation(xlsx, sheet=None, id_col="NR.ID"):
    if xlsx is None:
        return None
    if sheet is None:
        sheet = pd.ExcelFile(xlsx).sheet_names[0]
    ann = pd.read_excel(xlsx, sheet_name=sheet)
    req = ["contig_1ab","start1ab_pos","end1ab_pos"]
    miss = [c for c in req if c not in ann.columns]
    if miss:
        raise ValueError(f"Annotation missing columns: {miss}")
    if id_col not in ann.columns:
        ann[id_col] = [f"GENE_{i+1}" for i in range(len(ann))]
    out = ann[["contig_1ab","start1ab_pos","end1ab_pos", id_col]].copy()
    out.columns = ["chrom","gstart","gend","gene_id"]
    out["chrom"]  = out["chrom"].astype(str)
    out["gstart"] = pd.to_numeric(out["gstart"], errors="coerce").astype("Int64")
    out["gend"]   = pd.to_numeric(out["gend"],   errors="coerce").astype("Int64")
    out = out.dropna(subset=["chrom","gstart","gend"])
    out["gstart"] = out["gstart"].astype(np.int64)
    out["gend"]   = out["gend"].astype(np.int64)
    return out.sort_values(["chrom","gstart","gend"]).reset_index(drop=True)

def gene_overlap_and_distance(genes_df: pd.DataFrame, chrom: str, w_start: int, w_end: int):
    """Return (n_overlap, overlap_ids_str, nearest_dist_bp). Distance=0 if overlap present, else min edge-to-edge distance."""
    if genes_df is None:
        return (np.nan, "", np.nan)
    sub = genes_df[genes_df["chrom"] == chrom]
    if sub.empty:
        return (0, "", np.nan)
    # overlaps
    ov = sub[(sub["gend"] >= w_start) & (sub["gstart"] <= w_end)]
    n_ov = int(ov.shape[0])
    ids  = ";".join(map(str, ov["gene_id"].tolist())) if n_ov>0 else ""
    if n_ov > 0:
        return (n_ov, ids, 0)
    # distance (min of left/right gap)
    left_gap  = (w_start - sub["gend"]).clip(lower=0)
    right_gap = (sub["gstart"] - w_end).clip(lower=0)
    dist = int(pd.concat([left_gap, right_gap], axis=1).min(axis=1).min())
    return (0, "", dist)

# =========================
# LOAD TMRCA SEGMENTS
# =========================
os.makedirs(OUTDIR, exist_ok=True)
seg = pd.read_csv(IN_TMRCA_TSV, sep="\t", dtype={"CHROM":str})
need = ["CHROM","start_tmrca","end_tmrca","mean_tmrca","lower_CI","upper_CI"]
miss = [c for c in need if c not in seg.columns]
if miss:
    raise ValueError(f"Input {IN_TMRCA_TSV} missing columns: {miss}")

for c in need[1:]:
    seg[c] = pd.to_numeric(seg[c], errors="coerce")
seg = seg.dropna(subset=["CHROM","start_tmrca","end_tmrca","mean_tmrca"]).copy()
seg["start_tmrca"] = seg["start_tmrca"].astype(np.int64)
seg["end_tmrca"]   = seg["end_tmrca"].astype(np.int64)

# =========================
# BUILD WINDOWED METRICS
# =========================
scaffolds = sorted(seg["CHROM"].unique().tolist(), key=scaffold_sort_key)

genes = read_annotation(ANNO_XLSX, sheet=ANNO_SHEET, id_col=ANNO_ID_COL)

rows = []
for ch in scaffolds:
    d = seg[seg["CHROM"]==ch].sort_values(["start_tmrca","end_tmrca"])
    if d.empty: 
        continue
    xmin, xmax = int(d["start_tmrca"].min()), int(d["end_tmrca"].max())
    win_arr = windows_for_scaffold(xmin, xmax, WIN_BP, STEP_BP)
    ss = d["start_tmrca"].to_numpy(np.int64)
    ee = d["end_tmrca"].to_numpy(np.int64)
    mm = d["mean_tmrca"].to_numpy(float)
    ll = d["lower_CI"].to_numpy(float)
    uu = d["upper_CI"].to_numpy(float)

    for ws, we in win_arr:
        bp_cov, wmean, wlo, whi, nseg = weighted_window_stats(ss,ee,mm,ll,uu, ws, we)
        n_g, ids, ndist = gene_overlap_and_distance(genes, ch, int(ws), int(we))
        rows.append({
            "CHROM": ch,
            "WIN_START": int(ws),
            "WIN_END":   int(we),
            "WIN_MID":   int((ws+we)//2),
            "bp_covered": int(bp_cov),
            "n_segments": int(nseg),
            "w_mean_tmrca": wmean,
            "w_lower_CI": wlo,
            "w_upper_CI": whi,
            "n_overlap_genes": n_g,
            "overlap_gene_ids": ids,
            "nearest_gene_dist": ndist
        })

win = pd.DataFrame(rows)
if win.empty:
    raise RuntimeError("No windows were created. Check WIN_BP/STEP_BP and input coordinates.")

# flag top 5% by window mean (non-NaN)
q95 = win["w_mean_tmrca"].quantile(1.0 - TOP_PCT, interpolation="linear")
win["FLAG"] = np.where(win["w_mean_tmrca"] >= q95, "TOP5%", "BACKGROUND")

# save window table
out_tsv = os.path.join(OUTDIR, f"tmrca_windows_win{WIN_BP}_step{STEP_BP}.tsv")
win.to_csv(out_tsv, sep="\t", index=False)

# =========================
# FISHER ENRICHMENT (TOP vs BKG) by distance thresholds on WINDOWS
# =========================
fisher_rows = []
sub = win[pd.notna(win["nearest_gene_dist"])].copy()
for D in DIST_THRESH:
    top_in  = int((sub["FLAG"].eq("TOP5%") & (sub["nearest_gene_dist"] <= D)).sum())
    top_out = int((sub["FLAG"].eq("TOP5%") & (sub["nearest_gene_dist"] >  D)).sum())
    bkg_in  = int((sub["FLAG"].eq("BACKGROUND") & (sub["nearest_gene_dist"] <= D)).sum())
    bkg_out = int((sub["FLAG"].eq("BACKGROUND") & (sub["nearest_gene_dist"] >  D)).sum())
    _, pval = fisher_exact([[top_in, top_out],[bkg_in, bkg_out]], alternative="two-sided")
    top_rate = top_in / max(1, (top_in + top_out))
    bkg_rate = bkg_in / max(1, (bkg_in + bkg_out))
    rr = (top_rate / bkg_rate) if bkg_rate>0 else math.inf
    fisher_rows.append(dict(distance_bp=D, TOP_in=top_in, TOP_out=top_out,
                            BKG_in=bkg_in, BKG_out=bkg_out,
                            TOP_rate=top_rate, BKG_rate=bkg_rate,
                            rate_ratio=rr, fisher_p=pval))
fisher_df = pd.DataFrame(fisher_rows)
fisher_df.to_csv(os.path.join(OUTDIR, "fisher_near_gene_windows.tsv"), sep="\t", index=False)

# =========================
# PLOTTING HELPERS
# =========================
def build_concat_positions(win_df: pd.DataFrame):
    chroms = sorted(win_df["CHROM"].unique(), key=scaffold_sort_key)
    offsets, xticks, xtlabs = {}, [], []
    cursor = 0
    for i, ch in enumerate(chroms):
        sub = win_df[win_df["CHROM"]==ch]
        if sub.empty: continue
        max_mid = int(sub["WIN_MID"].max())
        offsets[ch] = cursor
        xticks.append(cursor + max_mid//2)
        xtlabs.append(ch)
        cursor += max_mid + 1
    return offsets, list(zip(xtlabs, xticks))

offsets, chrom_ticks = build_concat_positions(win)

def concat_x(df):
    return df["WIN_MID"] + df["CHROM"].map(offsets).fillna(0).astype(np.int64)

# Optional dimming: alpha=ALPHA_LOW for windows below the 75th percentile
q75 = win["w_mean_tmrca"].quantile(0.75)
alphas = np.where(win["w_mean_tmrca"] < q75, ALPHA_LOW, POINT_ALPHA)

# =========================
# MANHATTAN-style SCATTER (linear and log10)
# =========================
def manhattan_scatter(df, ycol, log10=False, out_png="scatter.png", title=""):
    cmap = {ch: ("#4C78A8" if (i % 2 == 0) else "#A0A0A0")
            for i, ch in enumerate(sorted(df["CHROM"].unique(), key=scaffold_sort_key))}
    colors = df["CHROM"].map(cmap).values
    X = concat_x(df)
    Y = df[ycol].astype(float)

    # Top-N by ycol
    topN = df.nlargest(TOP_N_ANNOTATE, ycol).copy()
    topN_x = concat_x(topN); topN_y = topN[ycol].astype(float)

    fig, ax = plt.subplots(figsize=(18,5))
    ax.scatter(X, Y, s=SCATTER_SIZE, c=colors, alpha=alphas, linewidths=0)
    ax.scatter(topN_x, topN_y, s=SCATTER_SIZE*2.2, c="red", alpha=1.0, linewidths=0, zorder=3)
    for _, r in topN.iterrows():
        label = f"{r['CHROM']}:{int(r['WIN_START'])}-{int(r['WIN_END'])}"
        ax.text(int(r['WIN_MID'])+offsets[r['CHROM']], float(r[ycol])+0.02*(np.nanmax(Y) if not log10 else 1),
                label, fontsize=LABEL_SIZE-2, color="red", rotation=35, ha="left", va="bottom")

    ax.set_xlabel("Genome (scaffolds concatenated)", fontsize=AXIS_LABEL)
    ax.set_ylabel(("w_mean_tmrca (generations, log10)" if log10 else "w_mean_tmrca (generations)"),
                  fontsize=AXIS_LABEL)
    ax.set_xticks([t for _,t in chrom_ticks])
    ax.set_xticklabels([c for c,_ in chrom_ticks], rotation=90, fontsize=TICK_SIZE)
    if log10:
        ax.set_yscale("log")
    ax.set_title(title, fontsize=TITLE_SIZE)
    fig.tight_layout(); fig.savefig(out_png, dpi=220); plt.close(fig)

manhattan_scatter(
    win, "w_mean_tmrca", log10=False,
    out_png=os.path.join(OUTDIR, f"manhattan_linear_win{WIN_BP}_step{STEP_BP}.png"),
    title=f"Windowed TMRCA Manhattan plot (linear)  win={WIN_BP:,}, step={STEP_BP:,}"
)

manhattan_scatter(
    win, "w_mean_tmrca", log10=True,
    out_png=os.path.join(OUTDIR, f"manhattan_log_win{WIN_BP}_step{STEP_BP}.png"),
    title=f"Windowed TMRCA Manhattan plot (log10)  win={WIN_BP:,}, step={STEP_BP:,}"
)

# =========================
# LINE PLOT WITH CI BAND (across concatenated genome)
# =========================
def line_with_ci(df, out_png, title):
    d = df.sort_values(["CHROM","WIN_START"]).copy()
    d["x"] = concat_x(d)
    fig, ax = plt.subplots(figsize=(20,5))
    ax.fill_between(d["x"], d["w_lower_CI"], d["w_upper_CI"], alpha=0.15, step="pre", color="#999999")
    ax.plot(d["x"], d["w_mean_tmrca"], lw=0.9, color="black")
    ax.set_xlabel("Genome (scaffolds concatenated)", fontsize=AXIS_LABEL)
    ax.set_ylabel("w_mean_tmrca (generations)", fontsize=AXIS_LABEL)
    ax.set_xticks([t for _,t in chrom_ticks])
    ax.set_xticklabels([c for c,_ in chrom_ticks], rotation=90, fontsize=TICK_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE)
    fig.tight_layout(); fig.savefig(out_png, dpi=220); plt.close(fig)

line_with_ci(
    win,
    out_png=os.path.join(OUTDIR, f"line_with_ci_win{WIN_BP}_step{STEP_BP}.png"),
    title=f"Windowed mean TMRCA with CI band  (win={WIN_BP:,}, step={STEP_BP:,})"
)

# =========================
# NOTES (how windows change Fisher)
# =========================
with open(os.path.join(OUTDIR, "README_analysis.txt"), "w") as fh:
    fh.write(
        "This analysis averages mean TMRCA on fixed overlapping windows using bp-overlap weights.\n"
        "All downstream statistics (TOP5% split, Fisher near-gene enrichment, plots) use these window means.\n"
        "Distance thresholds [0,1k,2k,5k,10k,50k] are now computed from window intervals to the nearest gene;\n"
        "distance=0 means the window overlaps ≥1 gene (n_overlap_genes>0). This is meaningful on windows because\n"
        "windows have uniform length and even coverage, reducing segment-length biases present in raw intervals.\n"
        f"95th percentile cutoff on window means (TOP {int(TOP_PCT*100)}%) = {q95:.3f}\n"
        f"Windows written to: {os.path.basename(out_tsv)}\n"
    )

print("\nDone.")
print("Outputs in:", os.path.abspath(OUTDIR))
print(" - Windows TSV:", os.path.basename(out_tsv))
print(" - Fisher table:", "fisher_near_gene_windows.tsv")
print(" - Plots: manhattan_linear_*.png, manhattan_log_*.png, line_with_ci_*.png")

