# tmrca_pipeline_v3.py
# ------------------------------------------------------------
# End-to-end TMRCA pipeline for ARGweaver 6-col tracks:
#   • Read: chrom start end tmrca_mean tmrca_lo tmrca_hi  (whitespace/TSV)
#   • Global & per-scaffold summaries (all vs zeros-removed)
#   • Coverage-weighted windowing with overlap
#   • TWO Manhattan plots (linear & log10). Zeros dropped.
#       - Alternating scaffold colors
#       - Points with y_mean < THRESH_FADE are transparent (alpha=ALPHA_FAINT)
#       - Top-K windows highlighted in red
#   • NEW: Line plot with confidence band (ribbon) from the raw TMRCA file
#       - Choose scaffold and optional region
#       - Linear or log10 y-scale
#   • PowerPoint deck with all figures and summary tables
# ------------------------------------------------------------

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pptx import Presentation
from pptx.util import Inches

# ============ USER SETTINGS ============
INPUT_FILE = "All_tmrca_60_run.txt"   # 6 columns, whitespace/TSV
OUTDIR     = "tmrca_reports"

# Coordinates for display (doesn't affect windowing)
ZERO_BASED_START = True

# Windowing for Manhattan
WIN_BP  = 1_000_000      # window size (bp)
STEP_BP = 100_000        # step (bp). Overlap = WIN_BP - STEP_BP

# Manhattan plot controls
TOP_K             = 10        # highlight the top-K windows by mean TMRCA
THRESH_FADE       = 500.0     # points with y_mean < this get a low alpha
ALPHA_FAINT       = 0.12      # alpha for faint points (< threshold)
ALPHA_STRONG      = 0.90      # alpha for points >= threshold
POINT_SIZE        = 12        # non-top points
POINT_SIZE_TOP    = 46        # red top-K points
ALT_SCAFF_COLORS  = ("#4C78A8", "#A0A0A0")  # alternating scaffold colors (blue/gray)

# NEW: confidence-band line plot (raw TMRCA track)
CI_SCAFFOLD   = "scaffold_2"   # which scaffold to draw
CI_X_START    = None           # None = auto (min). Otherwise an int bp coordinate.
CI_X_END      = None           # None = auto (max)
CI_Y_SCALE    = "log"          # "log" or "linear"
CI_MIN_POSVAL = 1e-9           # to avoid 0 on log axis

np.random.seed(13)
os.makedirs(OUTDIR, exist_ok=True)

# ============ HELPERS ============
def parse_scaffold_order(chrom: str):
    """
    Sort key so labels go: 1a, 1b, 2, 3, ..., 27, 31, 44, etc.
    Works with names like 'scaffold_1a', 'scaffold_2', '1b', '20', etc.
    """
    m = re.search(r'(\d+)([A-Za-z]?)$', chrom)
    if not m:  # unknown -> send to end
        return (10**9, 10**9, chrom)
    num = int(m.group(1))
    letter = m.group(2).lower()
    # order a(0), b(1), c(2), ..., ''(9 so number-without-letter after letters)
    letter_rank = (ord(letter) - ord('a')) if letter else 9
    return (num, letter_rank, chrom)

def window_reduce(df_scaf: pd.DataFrame, win_bp: int, step_bp: int) -> pd.DataFrame:
    """Coverage-weighted mean of tmrca_mean per sliding window in one scaffold."""
    start_min = int(max(0, df_scaf["start"].min()))
    end_max   = int(df_scaf["end"].max())
    starts = np.array([start_min], dtype=int) if end_max - start_min < win_bp \
             else np.arange(start_min, end_max - win_bp + 1, step_bp, dtype=int)

    d = df_scaf.sort_values(["start","end"]).reset_index(drop=True)
    s_arr = d["start"].to_numpy()
    e_arr = d["end"].to_numpy()
    y_arr = d["tmrca_mean"].to_numpy()

    rows, j0 = [], 0
    for wstart in starts:
        wend = wstart + win_bp
        while j0 < len(d) and e_arr[j0] <= wstart:
            j0 += 1
        j, cov, num = j0, 0, 0.0
        while j < len(d) and s_arr[j] < wend:
            ovl = max(0, min(e_arr[j], wend) - max(s_arr[j], wstart))
            if ovl > 0:
                cov += ovl
                num += ovl * y_arr[j]
            j += 1
        mean_y = (num / cov) if cov > 0 else np.nan
        rows.append((wstart, wend, mean_y, cov))
    return pd.DataFrame(rows, columns=["w_start","w_end","y_mean","bp_covered"])

def plot_manhattan(win_all, scaffolds_ordered, chrom_ticks, yscale, out_png, title_suffix):
    """Scatter without bands; fade points below THRESH_FADE; highlight top-K in red."""
    y = win_all["y_mean"].to_numpy()
    alphas = np.where(y >= THRESH_FADE, ALPHA_STRONG, ALPHA_FAINT)

    # Colors alternate by scaffold
    chrom2col = {sc: ALT_SCAFF_COLORS[i % 2] for i, sc in enumerate(scaffolds_ordered)}
    colors = np.array([chrom2col[c] for c in win_all["chrom"]], dtype=object)

    # Top-K by y_mean (positive only)
    topK = win_all[win_all["y_mean"] > 0].nlargest(TOP_K, "y_mean")
    mask_top = win_all.index.isin(topK.index)

    fig, ax = plt.subplots(figsize=(18, 4.5))
    # scatter: non-top
    ax.scatter(win_all.loc[~mask_top, "mid_global"],
               win_all.loc[~mask_top, "y_mean"],
               s=POINT_SIZE, c=colors[~mask_top], alpha=alphas[~mask_top], linewidths=0)
    # scatter: top-K
    ax.scatter(win_all.loc[mask_top, "mid_global"],
               win_all.loc[mask_top, "y_mean"],
               s=POINT_SIZE_TOP, c="red", alpha=1.0, linewidths=0, zorder=3)

    if yscale == "log":
        ax.set_yscale("log")
        ylab = "Windowed TMRCA mean (generations, log10)"
    else:
        ylab = "Windowed TMRCA mean (generations)"

    ax.set_ylabel(ylab)
    ax.set_xlabel("Genomic position (scaffolds concatenated)")
    ax.set_xticks([t for _, t in chrom_ticks])
    ax.set_xticklabels([c for c, _ in chrom_ticks], rotation=90, fontsize=9)
    ax.set_title(f"Windowed TMRCA Manhattan plot (win={WIN_BP:,} bp, step={STEP_BP:,} bp)\n"
                 f"Top {TOP_K} windows in red {title_suffix}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)

def plot_ci_track(raw_df: pd.DataFrame, scaffold: str,
                  x_start: int | None, x_end: int | None,
                  yscale: str, out_png: str):
    """
    Draw a step line (tmrca_mean) with confidence ribbon (lo-hi) for one scaffold/region.
    """
    d = raw_df.loc[raw_df["chrom"] == scaffold, ["start","end","tmrca_mean","tmrca_lo","tmrca_hi"]].copy()
    if d.empty:
        raise ValueError(f"No rows found for scaffold '{scaffold}'")

    # Region clamp
    xmin = d["start"].min() if x_start is None else int(x_start)
    xmax = d["end"].max()   if x_end   is None else int(x_end)
    d = d[(d["end"] >= xmin) & (d["start"] <= xmax)].copy()
    if d.empty:
        raise ValueError(f"No rows in requested interval {scaffold}:{xmin}-{xmax}")

    # Clamp segment bounds to [xmin, xmax]
    d["xs"] = np.maximum(d["start"], xmin)
    d["xe"] = np.minimum(d["end"],   xmax)
    # Avoid zeros on log scale
    eps = CI_MIN_POSVAL
    d["lo"] = np.maximum(d["tmrca_lo"].astype(float), eps)
    d["hi"] = np.maximum(d["tmrca_hi"].astype(float), eps)
    d["ym"] = np.maximum(d["tmrca_mean"].astype(float), eps)

    # Mb for x
    d["xs_mb"] = d["xs"] / 1e6
    d["xe_mb"] = d["xe"] / 1e6

    fig, ax = plt.subplots(figsize=(9, 4))
    # Ribbon (lo-hi)
    for _, r in d.iterrows():
        ax.fill_between([r["xs_mb"], r["xe_mb"]], [r["lo"], r["lo"]], [r["hi"], r["hi"]],
                        color="lightgrey", alpha=0.6, linewidth=0)
    # Step line for mean
    for _, r in d.iterrows():
        ax.plot([r["xs_mb"], r["xe_mb"]], [r["ym"], r["ym"]], color="black", linewidth=1.0)

    if yscale == "log":
        ax.set_yscale("log")
        ylab = "TMRCA (mean, log10 generations)"
    else:
        ylab = "TMRCA (mean, generations)"

    ax.set_xlabel(f"{scaffold} position (Mb)")
    ax.set_ylabel(ylab)
    ax.set_title(f"{scaffold}:{xmin:,}-{xmax:,}   (ARGweaver TMRCA with 95% CI)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)

# ============ 1) READ DATA ============
cols = ["chrom","start","end","tmrca_mean","tmrca_lo","tmrca_hi"]
tm = pd.read_csv(INPUT_FILE, sep=r"\s+|\t", engine="python", header=None, names=cols)
for c in ["start","end","tmrca_mean","tmrca_lo","tmrca_hi"]:
    tm[c] = pd.to_numeric(tm[c], errors="coerce")
tm = tm.dropna(subset=["chrom","start","end","tmrca_mean"]).reset_index(drop=True)
tm["start_disp"] = tm["start"] + 1 if ZERO_BASED_START else tm["start"]
tm["end_disp"]   = tm["end"]

# ============ 2) SUMMARIES ============
def summarize_series(x: pd.Series, label: str) -> pd.DataFrame:
    x = pd.to_numeric(x, errors="coerce")
    return pd.DataFrame([{
        "n": int(x.size),
        "n_nonzero": int((x>0).sum()),
        "min": float(x.min()),
        "q1": float(x.quantile(0.25)),
        "median": float(x.median()),
        "mean": float(x.mean()),
        "q3": float(x.quantile(0.75)),
        "max": float(x.max()),
        "sd": float(x.std())
    }], index=[label])

summ_all = summarize_series(tm["tmrca_mean"], "all_values")
summ_no0 = summarize_series(tm.loc[tm["tmrca_mean"]>0, "tmrca_mean"], "nonzero_only")
summ_table = pd.concat([summ_all, summ_no0])
summ_table.to_csv(os.path.join(OUTDIR, "summary_global_tmrca.csv"))

per_scaf = tm.groupby("chrom")["tmrca_mean"].agg(
    n="size",
    n_nonzero=lambda s: (s>0).sum(),
    min="min",
    q1=lambda s: s.quantile(0.25),
    median="median",
    mean="mean",
    q3=lambda s: s.quantile(0.75),
    max="max",
    sd="std"
).reset_index()
per_scaf.to_csv(os.path.join(OUTDIR, "summary_by_scaffold.csv"), index=False)

# ============ 3) WINDOWING ============
scaffolds_ordered = sorted(tm["chrom"].unique().tolist(), key=parse_scaffold_order)

win_tables = []
offset = 0
chrom_ticks = []
xpos_min = offset

for chrom in scaffolds_ordered:
    d = tm.loc[tm["chrom"] == chrom, ["start","end","tmrca_mean"]].copy()
    if d.empty:
        continue
    w = window_reduce(d, WIN_BP, STEP_BP)
    w["chrom"] = chrom
    w["mid"] = (w["w_start"] + w["w_end"]) / 2.0
    w["mid_global"] = w["mid"] + offset
    win_tables.append(w)

    xpos_max = offset + d["end"].max()
    chrom_ticks.append((chrom, (xpos_min + xpos_max) / 2.0))
    offset = xpos_max + 1
    xpos_min = offset

win_all = pd.concat(win_tables, ignore_index=True)

# keep only informative windows: coverage>0, finite mean, and mean>0 (drop zeros)
win_all = win_all[(win_all["bp_covered"] > 0) & (win_all["y_mean"].notna()) & (win_all["y_mean"] > 0)].copy()
if win_all.empty:
    raise RuntimeError("No windows with positive mean and coverage. Adjust WIN_BP/STEP_BP or check input.")

win_csv = os.path.join(OUTDIR, f"tmrca_windows_{WIN_BP}_{STEP_BP}.csv")
win_all.to_csv(win_csv, index=False)

# ============ 4) PLOTS ============
# Manhattan (linear)
out_lin = os.path.join(OUTDIR, f"manhattan_tmrca_LINEAR_win{WIN_BP}_step{STEP_BP}.png")
plot_manhattan(win_all, scaffolds_ordered, chrom_ticks, yscale="linear", out_png=out_lin,
               title_suffix="(linear scale)")
# Manhattan (log)
out_log = os.path.join(OUTDIR, f"manhattan_tmrca_LOG_win{WIN_BP}_step{STEP_BP}.png")
plot_manhattan(win_all, scaffolds_ordered, chrom_ticks, yscale="log", out_png=out_log,
               title_suffix="(log scale)")

# NEW: CI line plot on raw TMRCA for chosen scaffold/region
ci_png = os.path.join(OUTDIR, f"tmrca_ci_track_{CI_SCAFFOLD}.png")
plot_ci_track(tm, CI_SCAFFOLD, CI_X_START, CI_X_END, CI_Y_SCALE, ci_png)

# Optional: genome-wide histogram (zeros removed)
vals = tm["tmrca_mean"].replace(0, np.nan).dropna()
fig, ax = plt.subplots(figsize=(7,4))
ax.hist(np.log10(vals), bins=60, color="#747474", edgecolor="none")
ax.set_xlabel("log10(TMRCA mean, generations)")
ax.set_ylabel("Count")
ax.set_title("Genome-wide distribution of local TMRCA (zeros removed)")
fig.tight_layout()
hist_path = os.path.join(OUTDIR, "tmrca_hist_log10.png")
fig.savefig(hist_path, dpi=220)
plt.close(fig)

# ============ 5) TOP-K TABLE ============
topK = win_all.nlargest(TOP_K, "y_mean").copy()
topK["w_start_Mb"] = topK["w_start"]/1e6
topK["w_end_Mb"]   = topK["w_end"]/1e6
topK_csv = os.path.join(OUTDIR, f"top_{TOP_K}_windows.csv")
topK[["chrom","w_start","w_end","w_start_Mb","w_end_Mb","y_mean","bp_covered"]].to_csv(topK_csv, index=False)

# ============ 6) SLIDES ============
prs = Presentation()

# Title
s0 = prs.slides.add_slide(prs.slide_layouts[0])
s0.shapes.title.text = "ARGweaver TMRCA summary"
s0.placeholders[1].text = (
    f"Input file: {os.path.basename(INPUT_FILE)}\n"
    f"Windows: {WIN_BP:,} bp  |  Step: {STEP_BP:,} bp\n"
    f"Top {TOP_K} windows highlighted in red\n"
    f"Fade rule: alpha={ALPHA_FAINT} if mean < {THRESH_FADE} generations"
)

# Global stats
s1 = prs.slides.add_slide(prs.slide_layouts[5])
s1.shapes.title.text = "Global TMRCA (generations): all vs zeros-removed"
rows, cols = summ_table.shape[0]+1, summ_table.shape[1]+1
tbl = s1.shapes.add_table(rows, cols, Inches(0.5), Inches(1.2), Inches(9), Inches(1.4)).table
tbl.cell(0,0).text = ""
for j, c in enumerate(summ_table.columns, start=1):
    tbl.cell(0,j).text = c
for i, (idx, row) in enumerate(summ_table.iterrows(), start=1):
    tbl.cell(i,0).text = str(idx)
    for j, c in enumerate(summ_table.columns, start=1):
        v = row[c]
        tbl.cell(i,j).text = f"{v:,.3g}" if isinstance(v,(int,float)) and not pd.isna(v) else str(v)

# Per-scaffold stats table (first 20 rows for readability)
s2 = prs.slides.add_slide(prs.slide_layouts[5])
s2.shapes.title.text = "Per-scaffold summary (generations)"
sub = per_scaf.sort_values("chrom", key=lambda s: s.map(parse_scaffold_order)).head(20).copy()
rows, cols = len(sub)+1, len(sub.columns)
tbl2 = s2.shapes.add_table(rows, cols, Inches(0.5), Inches(1.2), Inches(9), Inches(3.5)).table
for j,c in enumerate(sub.columns):
    tbl2.cell(0,j).text = c
for i, (_, r) in enumerate(sub.iterrows(), start=1):
    for j, c in enumerate(sub.columns):
        val = r[c]
        tbl2.cell(i,j).text = f"{val:,.3g}" if isinstance(val,(int,float)) and not pd.isna(val) else str(val)

# Manhattan LINEAR
s3 = prs.slides.add_slide(prs.slide_layouts[5])
s3.shapes.title.text = "Windowed TMRCA Manhattan (linear scale)"
s3.shapes.add_picture(out_lin, Inches(0.5), Inches(1.2), width=Inches(9))

# Manhattan LOG
s4 = prs.slides.add_slide(prs.slide_layouts[5])
s4.shapes.title.text = "Windowed TMRCA Manhattan (log scale)"
s4.shapes.add_picture(out_log, Inches(0.5), Inches(1.2), width=Inches(9))

# NEW: CI ribbon plot slide
s5 = prs.slides.add_slide(prs.slide_layouts[5])
s5.shapes.title.text = f"TMRCA with 95% CI — {CI_SCAFFOLD}"
s5.shapes.add_picture(ci_png, Inches(0.5), Inches(1.2), width=Inches(9))

# Histogram slide
s6 = prs.slides.add_slide(prs.slide_layouts[5])
s6.shapes.title.text = "Genome-wide distribution of local TMRCA (zeros removed)"
s6.shapes.add_picture(hist_path, Inches(0.5), Inches(1.2), width=Inches(9))

pptx_path = os.path.join(OUTDIR, "tmrca_summary_slides_v3.pptx")
prs.save(pptx_path)

print("Wrote:")
print("  -", win_csv)
print("  -", topK_csv)
print("  -", out_lin)
print("  -", out_log)
print("  -", ci_png)
print("  -", hist_path)
print("  -", pptx_path)
