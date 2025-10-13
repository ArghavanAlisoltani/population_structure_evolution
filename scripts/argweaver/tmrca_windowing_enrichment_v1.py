#!/usr/bin/env python3
"""
Windowed TMRCA + gene enrichment with binomial tests, figures, and robust handling
of short segments. Designed for ARGweaver outputs on GBS data.

Input TMRCA (tab-delimited):
  - 6 cols:  CHROM, start, end, tmrca_mean, tmrca_lo, tmrca_hi
  - or 7 cols: CHROM, start, end, tmrca_mean, tmrca_mode, tmrca_lo, tmrca_hi
    (mode is ignored unless you switch to it explicitly)

Gene annotation (CSV/TSV):
  Required cols: CHROM, GENE_START, GENE_END, GENE_ID   (1-based, inclusive)

Outputs:
  - tmrca_windows_win{WIN}_ov{OVERLAP}.tsv      (full window table with flags)
  - manhattan_linear.png / manhattan_log.png    (scatter by scaffold)
  - line_with_ci.png                             (coverage-weighted mean + ribbon)
  - seglen_hist.png / seglen_boxplot.png         (segment length diagnostics)
  - binom_enrichment.tsv                         (per-scaffold binomial tests)
"""

from pathlib import Path
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

# ------------------------
# User settings
# ------------------------
TMRCA_PATH   = "All_tmrca_60_run.txt"          # your big tmrca file (segments across scaffolds)
GENES_PATH   = "genes_1ab.tsv"                 # TSV/CSV with CHROM,GENE_START,GENE_END,GENE_ID
OUTDIR       = Path("tmrca_windowed_out")
OUTDIR.mkdir(exist_ok=True, parents=True)

# ARGweaver coordinate assumption:
ARG_ZERO_BASED = True     # set True if starts are 0-based half-open; False if already 1-based inclusive

# Windowing
WIN_BP       = 1_000_000  # window size (bp)
OVERLAP_BP   = 100_000    # exact overlap between successive windows (bp)
MIN_COVER_BP = 10_000     # min covered bp inside a window to keep it (avoid windows with almost no segments)

# Which TMRCA center to use for window averaging
CENTER = "mean"           # "mean" or "mode"; CI is always from lo/hi

# Percentile cut (after excluding zeros if desired)
REMOVE_ZERO_FOR_CUTOFF = True
TOP_PCT   = 95.0          # 95th percentile = top 5%

# Distances for "near a gene" (including 0 => strict overlap)
DIST_THRESH = [0, 1_000, 2_000, 5_000, 10_000, 50_000]

# Plot knobs
FIG_W = 14
FIG_H = 4
FONTSIZE = 12
ALPHA_BELOW_Q3 = 0.35     # make dots below Q3 less visible
DOT_SIZE = 8

# ------------------------
# Helpers
# ------------------------
def read_tmrca(path: str, zero_based=True) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+|\t", engine="python", header=None)
    if df.shape[1] == 6:
        df.columns = ["CHROM","start","end","tmrca_mean","tmrca_lo","tmrca_hi"]
        df["tmrca_mode"] = np.nan
    elif df.shape[1] == 7:
        df.columns = ["CHROM","start","end","tmrca_mean","tmrca_mode","tmrca_lo","tmrca_hi"]
    else:
        raise ValueError(f"Unexpected columns ({df.shape[1]}).")
    # coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    if zero_based:
        df["start"] = df["start"].astype("int64") + 1
    df["end"]   = df["end"].astype("int64")
    # numeric TMRCA
    for c in ["tmrca_mean","tmrca_mode","tmrca_lo","tmrca_hi"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").replace([np.inf,-np.inf], np.nan)
    # length and CI width
    df["seg_len"] = (df["end"] - df["start"] + 1).clip(lower=0)
    df["ci_width"] = (df[["tmrca_lo","tmrca_hi"]].max(axis=1) - df[["tmrca_lo","tmrca_hi"]].min(axis=1))
    # which center
    if CENTER == "mode" and df["tmrca_mode"].notna().any():
        df["center"] = df["tmrca_mode"].fillna(df["tmrca_mean"])
    else:
        df["center"] = df["tmrca_mean"]
    # clean negatives/zeros safely
    for c in ["center","tmrca_mean","tmrca_lo","tmrca_hi"]:
        df[c] = df[c].astype(float)
        df[c] = df[c].mask(df[c] < 0, 0.0)
    df = df.dropna(subset=["CHROM","start","end","center"])
    df = df[df["end"] >= df["start"]].copy()
    return df.reset_index(drop=True)

def read_genes(path: str) -> pd.DataFrame:
    sep = "," if path.lower().endswith(".csv") else "\t"
    g = pd.read_csv(path, sep=sep)
    need = ["CHROM","GENE_START","GENE_END","GENE_ID"]
    miss = [c for c in need if c not in g.columns]
    if miss: raise ValueError(f"Missing gene columns: {miss}")
    g = g.copy()
    g["GENE_START"] = pd.to_numeric(g["GENE_START"], errors="coerce").astype("Int64")
    g["GENE_END"]   = pd.to_numeric(g["GENE_END"],   errors="coerce").astype("Int64")
    g = g.dropna(subset=["CHROM","GENE_START","GENE_END"])
    g["GENE_START"] = g["GENE_START"].astype("int64")
    g["GENE_END"]   = g["GENE_END"].astype("int64")
    g = g[g["GENE_END"] >= g["GENE_START"]].copy()
    return g.reset_index(drop=True)

def scaffold_sort_key(chrom: str):
    # put 1a, 1b first, then numeric
    if chrom.endswith("_1a"): return (-1, chrom)
    if chrom.endswith("_1b"): return (0, chrom)
    m = re.search(r"(\d+)", chrom)
    return (int(m.group(1)) if m else 10**12, chrom)

def build_windows(df: pd.DataFrame, win_bp: int, overlap_bp: int) -> pd.DataFrame:
    """Return a DataFrame of windows per scaffold: CHROM, WSTART, WEND, WMID"""
    assert overlap_bp < win_bp, "OVERLAP_BP must be smaller than WIN_BP"
    rows = []
    for chrom, d in df.groupby("CHROM", sort=False):
        if d.empty: continue
        start_min = int(d["start"].min())
        end_max   = int(d["end"].max())
        step = win_bp - overlap_bp
        starts = np.arange(start_min, end_max + 1, step, dtype=np.int64)
        ends   = np.minimum(starts + win_bp - 1, end_max)
        mids   = (starts + ends) // 2
        rows.append(pd.DataFrame({"CHROM": chrom, "WSTART": starts, "WEND": ends, "WMID": mids}))
    if not rows: return pd.DataFrame(columns=["CHROM","WSTART","WEND","WMID"])
    out = pd.concat(rows, ignore_index=True)
    # order scaffolds
    chroms = sorted(out["CHROM"].unique(), key=scaffold_sort_key)
    out["CHROM"] = pd.Categorical(out["CHROM"], categories=chroms, ordered=True)
    out = out.sort_values(["CHROM","WSTART"]).reset_index(drop=True)
    return out

def window_weighted_means(segments: pd.DataFrame, windows: pd.DataFrame, min_cover_bp=0) -> pd.DataFrame:
    """Coverage-weighted mean for center, lo, hi within each window."""
    out = []
    for chrom, w in windows.groupby("CHROM", sort=False):
        seg = segments[segments["CHROM"] == chrom]
        if seg.empty: 
            continue
        # two-pointer sweep
        seg = seg.sort_values("start").reset_index(drop=True)
        s_idx = 0
        for _, row in w.iterrows():
            ws, we = int(row.WSTART), int(row.WEND)
            cover = 0
            w_num_center = 0.0
            w_num_loww   = 0.0
            w_num_highw  = 0.0
            # advance seg pointer
            while s_idx < len(seg) and seg.loc[s_idx, "end"] < ws:
                s_idx += 1
            j = s_idx
            while j < len(seg) and seg.loc[j, "start"] <= we:
                ss = int(seg.loc[j, "start"]); ee = int(seg.loc[j, "end"])
                ovl = max(0, min(we, ee) - max(ws, ss) + 1)
                if ovl > 0:
                    c  = float(seg.loc[j, "center"])
                    lo = float(seg.loc[j, "tmrca_lo"])
                    hi = float(seg.loc[j, "tmrca_hi"])
                    w_num_center += ovl * c
                    w_num_loww   += ovl * lo
                    w_num_highw  += ovl * hi
                    cover += ovl
                j += 1
            if cover >= min_cover_bp:
                out.append({
                    "CHROM": chrom,
                    "WSTART": ws, "WEND": we, "WMID": int((ws+we)//2),
                    "COVER_BP": cover,
                    "W_CENTER": w_num_center/cover if cover>0 else np.nan,
                    "W_LO":     w_num_loww/cover   if cover>0 else np.nan,
                    "W_HI":     w_num_highw/cover  if cover>0 else np.nan
                })
    if not out: 
        return pd.DataFrame(columns=["CHROM","WSTART","WEND","WMID","COVER_BP","W_CENTER","W_LO","W_HI"])
    return pd.DataFrame(out)

def add_near_gene_counts(win: pd.DataFrame, genes: pd.DataFrame, dists: list[int]) -> pd.DataFrame:
    """
    For each window and each distance threshold D, mark is_near_D (>=1 gene)
    and count genes (Ngenes_D) where dist <= D (dist=0 means overlap).
    """
    win = win.copy()
    for D in dists:
        win[f"is_near_{D}"] = False
        win[f"Ngenes_{D}"]  = 0

    for chrom, w in win.groupby("CHROM", sort=False):
        g = genes[genes["CHROM"] == chrom]
        if w.empty: continue
        if g.empty:
            continue
        # sort for fast sweep
        g = g.sort_values("GENE_START").reset_index(drop=True)
        # naive but readable scan per window (ok if window count ~ tens of thousands)
        for idx, row in w.iterrows():
            ws, we = int(row.WSTART), int(row.WEND)
            gm = g[((g["GENE_END"] >= ws - max(dists)) & (g["GENE_START"] <= we + max(dists)))]
            if gm.empty:
                continue
            # distances (0 if overlap)
            # distance = min(|gene_start - window_end|, |window_start - gene_end|) when disjoint; else 0
            d_left  = np.maximum(0, ws - gm["GENE_END"])
            d_right = np.maximum(0, gm["GENE_START"] - we)
            dist = np.minimum(d_left, d_right).to_numpy()
            for D in dists:
                mask = dist <= D
                c = int(np.count_nonzero(mask))
                if c > 0:
                    win.at[idx, f"is_near_{D}"] = True
                    win.at[idx, f"Ngenes_{D}"]  = c
    return win

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, float)
    ok = np.isfinite(pvals)
    q = np.full_like(pvals, np.nan, dtype=float)
    if ok.sum():
        q[ok] = multipletests(pvals[ok], method="fdr_bh")[1]
    return q

def per_scaffold_binomial(win: pd.DataFrame, flag_col: str, near_col: str) -> pd.DataFrame:
    """
    For each scaffold, test whether the fraction of Top windows that are near genes
    is higher than the fraction among Background windows on the same scaffold.
    One-sided binomial with null p0 = background rate on that scaffold.
    """
    rows = []
    for chrom, d in win.groupby("CHROM", sort=False):
        top = d[d[flag_col]]
        bkg = d[~d[flag_col]]
        if d.empty or bkg.empty:
            continue
        # background rate on this scaffold
        p0 = float(bkg[near_col].mean())
        K  = int(top[near_col].sum())   # successes in Top
        n  = int(top.shape[0])          # trials in Top
        if n == 0:
            continue
        # P(X >= K | n, p0)
        p = 1.0 - binom.cdf(K-1, n, p0) if p0 > 0 else (1.0 if K == 0 else 0.0)
        rows.append({"CHROM": chrom, "p_binom": p, "n_top": n, "k_top": K, "p0_bkg": p0})
    if not rows:
        return pd.DataFrame(columns=["CHROM","p_binom","q_binom","n_top","k_top","p0_bkg"])
    out = pd.DataFrame(rows).reset_index(drop=True)
    out["q_binom"] = bh_fdr(out["p_binom"].to_numpy())
    return out.sort_values("q_binom")

# ------------------------
# Plots
# ------------------------
def seglen_hist(df: pd.DataFrame, out_png: Path):
    plt.figure(figsize=(8,3.2))
    x = df["seg_len"].clip(lower=1).astype(float).values
    plt.hist(x, bins=60, log=False)
    plt.xscale("log")
    plt.xlabel("Segment length (bp, log scale)")
    plt.ylabel("Count")
    plt.title("TMRCA segment lengths (log-scale histogram)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200); plt.close()

def seglen_top_box(df: pd.DataFrame, cutoff: float, out_png: Path):
    d = df.copy()
    d["group"] = np.where(d["center"] >= cutoff, "TOP5%", "BACKGROUND")
    plt.figure(figsize=(5,6))
    d.boxplot(column="seg_len", by="group")
    plt.yscale("log")
    plt.title("")
    plt.suptitle("")
    plt.ylabel("Segment length (bp, log)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200); plt.close()

def make_offsets(w: pd.DataFrame) -> tuple[pd.DataFrame, dict, list, list]:
    chroms = list(w["CHROM"].cat.categories) if hasattr(w["CHROM"], "cat") else sorted(w["CHROM"].unique(), key=scaffold_sort_key)
    offsets = {}
    cursor = 0
    xticks, xtlabs = [], []
    augmented = []
    for ch in chroms:
        sub = w[w["CHROM"] == ch]
        if sub.empty: continue
        offsets[ch] = cursor
        max_pos = int(sub["WMID"].max())
        xticks.append(cursor + max_pos//2)
        xtlabs.append(ch)
        tmp = sub.copy()
        tmp["XPOS"] = tmp["WMID"] + cursor
        augmented.append(tmp)
        cursor += max_pos + 1
    if not augmented:
        return w.copy(), offsets, [], []
    aug = pd.concat(augmented, ignore_index=True)
    return aug, offsets, xticks, xtlabs

def manhattan_scatter(win: pd.DataFrame, out_png: Path, use_log=False, annotate_topN=10):
    w = win.copy()
    # hide zeros
    w = w[w["W_CENTER"] > 0].copy()
    # alpha rule: below Q3 => lower alpha
    q3 = np.nanquantile(w["W_CENTER"], 0.75) if not use_log else np.nanquantile(np.log10(w["W_CENTER"].clip(lower=1e-9)), 0.75)
    w, offsets, xticks, xtlabs = make_offsets(w)
    # color by scaffold alternating
    chroms = list(w["CHROM"].cat.categories) if hasattr(w["CHROM"], "cat") else sorted(w["CHROM"].unique(), key=scaffold_sort_key)
    colA, colB = "#1f77b4", "#7f7f7f"
    colors = {ch: (colA if i%2==0 else colB) for i, ch in enumerate(chroms)}
    # top10 by value
    order = w.sort_values("W_CENTER", ascending=False).head(annotate_topN)
    plt.figure(figsize=(FIG_W, FIG_H))
    for ch, d in w.groupby("CHROM", sort=False):
        y = d["W_CENTER"].to_numpy(float)
        if use_log: y = np.log10(np.clip(y, 1e-9, None))
        a = np.where((d["W_CENTER"] if not use_log else np.log10(np.clip(d["W_CENTER"],1e-9,None))) < q3, ALPHA_BELOW_Q3, 0.9)
        plt.scatter(d["XPOS"], y, s=DOT_SIZE, c=colors[ch], alpha=a, edgecolors="none")
    # annotate top
    yo = order["W_CENTER"].to_numpy(float)
    if use_log: yo = np.log10(np.clip(yo, 1e-9, None))
    plt.scatter(order["XPOS"], yo, s=40, c="red", zorder=5)
    for _, r in order.iterrows():
        yy = math.log10(r["W_CENTER"]) if use_log else r["W_CENTER"]
        plt.text(r["XPOS"], yy, f"{r['CHROM']}:{int(r['WSTART'])}-{int(r['WEND'])}",
                 color="red", fontsize=8, rotation=45, ha="left", va="bottom")
    plt.xticks(xticks, xtlabs, rotation=90, fontsize=9)
    plt.ylabel("Windowed TMRCA mean" + (" (log10)" if use_log else " (generations)"))
    plt.xlabel("Genomic position (scaffolds concatenated)")
    plt.title(f"Windowed TMRCA Manhattan plot (win={WIN_BP:,} bp, overlap={OVERLAP_BP:,} bp)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220); plt.close()

def line_with_ci(win: pd.DataFrame, out_png: Path, use_log=False):
    w = win.copy()
    w, offsets, xticks, xtlabs = make_offsets(w)
    # sort by XPOS for line continuity
    w = w.sort_values("XPOS")
    y  = w["W_CENTER"].to_numpy(float)
    lo = w["W_LO"].to_numpy(float); hi = w["W_HI"].to_numpy(float)
    if use_log:
        y  = np.log10(np.clip(y, 1e-9, None))
        lo = np.log10(np.clip(lo, 1e-9, None))
        hi = np.log10(np.clip(hi, 1e-9, None))
    plt.figure(figsize=(FIG_W, FIG_H))
    plt.plot(w["XPOS"], y, lw=0.8, c="black")
    plt.fill_between(w["XPOS"], lo, hi, color="grey", alpha=0.2, step="mid")
    # faint scaffold bands
    for i, (ch, d) in enumerate(w.groupby("CHROM", sort=False)):
        if i % 2 == 1:
            plt.axvspan(d["XPOS"].min(), d["XPOS"].max(), color="#000000", alpha=0.03)
    plt.xticks(xticks, xtlabs, rotation=90, fontsize=9)
    plt.ylabel("Windowed TMRCA mean" + (" (log10)" if use_log else " (generations)"))
    plt.xlabel("Genomic position (scaffolds concatenated)")
    plt.title("Windowed TMRCA line plot with CI ribbon")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220); plt.close()

# ------------------------
# Main
# ------------------------
def main():
    # 1) Read
    seg = read_tmrca(TMRCA_PATH, zero_based=ARG_ZERO_BASED)
    genes = read_genes(GENES_PATH)

    # 2) Segment-length diagnostics
    seglen_hist(seg, OUTDIR / "seglen_hist.png")

    # 3) Windows
    wins = build_windows(seg, WIN_BP, OVERLAP_BP)
    if wins.empty:
        print("No windows produced."); return

    # 4) Coverage-weighted window means
    wagg = window_weighted_means(seg, wins, min_cover_bp=MIN_COVER_BP)
    if wagg.empty:
        print("No windows passed coverage threshold."); return

    # 5) Top 5% cutoff on windowed mean
    vals = wagg["W_CENTER"].replace([np.inf,-np.inf], np.nan)
    if REMOVE_ZERO_FOR_CUTOFF:
        vals = vals.mask(vals <= 0, np.nan)
    cutoff = float(np.nanpercentile(vals, TOP_PCT))
    wagg["TOP5"] = wagg["W_CENTER"] >= cutoff
    wagg["CUTOFF_95pct"] = cutoff

    # 6) Top10 windows for plots
    wagg = wagg.sort_values(["CHROM","WSTART"]).reset_index(drop=True)

    # 7) Add gene proximity/overlap counts per threshold
    wnear = add_near_gene_counts(wagg, genes, DIST_THRESH)

    # 8) Binomial enrichment per scaffold (using “near 0” == overlap)
    enrich = per_scaffold_binomial(wnear, flag_col="TOP5", near_col="is_near_0")
    enrich.to_csv(OUTDIR / "binom_enrichment.tsv", sep="\t", index=False)

    # 9) Save the window table
    out_tsv = OUTDIR / f"tmrca_windows_win{WIN_BP}_ov{OVERLAP_BP}.tsv"
    wnear.to_csv(out_tsv, sep="\t", index=False)

    # 10) Segment-length boxplot using Top5 cutoff derived from WINDOWED values
    seglen_top_box(seg, cutoff=cutoff, out_png=OUTDIR / "seglen_boxplot.png")

    # 11) Manhattan (linear + log)
    manhattan_scatter(wnear, OUTDIR / "manhattan_linear.png", use_log=False, annotate_topN=10)
    manhattan_scatter(wnear, OUTDIR / "manhattan_log.png",    use_log=True,  annotate_topN=10)

    # 12) Line plots with CI ribbons (linear + log)
    line_with_ci(wnear, OUTDIR / "line_with_ci_linear.png", use_log=False)
    line_with_ci(wnear, OUTDIR / "line_with_ci_log.png",    use_log=True)

    print(f"Done.\nWindows table: {out_tsv}\nEnrichment: {OUTDIR/'binom_enrichment.tsv'}")
    print(f"Figures in {OUTDIR}")

if __name__ == "__main__":
    main()
