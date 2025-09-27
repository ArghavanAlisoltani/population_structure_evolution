# tmrca_pipeline_manhattan.py
# ------------------------------------------------------------
# Read ARGweaver 6-col TMRCA track, window with overlap, and draw
# Manhattan-style plots with:
#   • scaffolds ordered 1a, 1b, 2, 3, … (numeric then letter)
#   • TWO figures: log10 Y and linear Y
#   • hide zeros; all points shown but those < Q3 are fainter (alpha=0.5)
#   • top-K windows highlighted in red
# Produces tables as before.
# ------------------------------------------------------------

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    from pptx import Presentation
    from pptx.util import Inches
except Exception:
    Presentation = None

# --------------------
# USER SETTINGS
# --------------------
INPUT_FILE = "All_tmrca_60_run.txt"   # tab-delimited: chrom start end mean lo hi
OUTDIR     = "tmrca_reports"

ZERO_BASED_START = True                # for display columns (not used in plots)
WIN_BP  = 1_000_000                    # window size (bp)
STEP_BP = 100_000                      # step / overlap
TOP_K   = 10                           # highlight top-K windows

# --------------------
# helpers
# --------------------
def parse_scaffold_order(chrom: str):
    """
    Return a sortable tuple: (number, letter_rank, original)
    Accepts things like 'scaffold_1a', 'scaffold_1b', 'scaffold_2', 'scaffold_27', '1a', '2', etc.
    Letters are ranked: a < b < c < '' (no letter goes AFTER letters of same number).
    """
    m = re.search(r'(\d+)([A-Za-z]?)$', chrom)
    if not m:
        return (float('inf'), float('inf'), chrom)
    num = int(m.group(1))
    letter = m.group(2).lower()
    # order a(0), b(1), c(2), …, ''(9) so 1a, 1b, then 1 (if exists)
    letter_rank = ord(letter) - ord('a') if letter else 9
    return (num, letter_rank, chrom)

def window_reduce(df_scaf: pd.DataFrame, win_bp: int, step_bp: int) -> pd.DataFrame:
    """Coverage-weighted mean of tmrca_mean per sliding window within a scaffold."""
    start_min = int(max(0, df_scaf["start"].min()))
    end_max   = int(df_scaf["end"].max())
    if end_max - start_min < win_bp:
        starts = np.array([start_min], dtype=int)
    else:
        starts = np.arange(start_min, end_max - win_bp + 1, step_bp, dtype=int)

    d = df_scaf.sort_values(["start","end"])
    s_arr = d["start"].to_numpy()
    e_arr = d["end"].to_numpy()
    y_arr = d["tmrca_mean"].to_numpy()

    rows, j0 = [], 0
    for wstart in starts:
        wend = wstart + win_bp
        while j0 < len(d) and e_arr[j0] <= wstart:
            j0 += 1
        j = j0
        cov = 0
        num = 0.0
        while j < len(d) and s_arr[j] < wend:
            ovl = max(0, min(e_arr[j], wend) - max(s_arr[j], wstart))
            if ovl > 0:
                cov += ovl
                num += ovl * y_arr[j]
            j += 1
        y = num/cov if cov > 0 else np.nan
        rows.append((wstart, wend, y, cov))
    return pd.DataFrame(rows, columns=["w_start","w_end","y_mean","bp_covered"])

def make_manhattan(win_all: pd.DataFrame, scaffolds_ordered, chrom_ticks,
                   out_png, yscale="log", title_suffix=""):
    """Draw Manhattan. Fainter alpha for points < per-scaffold Q3; >=Q3 alpha=0.9."""
    # Compute per-scaffold Q3 on window means (ignore zeros)
    q3 = (
        win_all[win_all["y_mean"] > 0]
        .groupby("chrom")["y_mean"]
        .quantile(0.75)
        .to_dict()
    )
    # alpha rule: < Q3 => 0.5, >= Q3 => 0.9
    alphas = np.where(
        win_all["y_mean"].to_numpy() >= win_all["chrom"].map(q3).to_numpy(),
        0.9, 0.5
    )

    # base colors alternating per scaffold
    base_colors = ["#4C78A8", "#A0A0A0"]  # blue, gray
    chrom_to_color = {chrom: base_colors[i % 2] for i, chrom in enumerate(scaffolds_ordered)}
    colors = np.array([chrom_to_color[c] for c in win_all["chrom"]], dtype=object)

    # Top-K in red
    topK = win_all[win_all["y_mean"] > 0].nlargest(TOP_K, "y_mean")
    mask_top = win_all.index.isin(topK.index)

    # Plot
    fig, ax = plt.subplots(figsize=(14, 4))
    # background bands
    for i, chrom in enumerate(scaffolds_ordered):
        scaf_mask = win_all["chrom"] == chrom
        if not scaf_mask.any():
            continue
        x_min = win_all.loc[scaf_mask, "mid_global"].min() - WIN_BP/2
        x_max = win_all.loc[scaf_mask, "mid_global"].max() + WIN_BP/2
        if i % 2 == 1:
            ax.axvspan(x_min, x_max, color="#000000", alpha=0.05, lw=0, zorder=0)

    # non-top points
    base_size = 14
    ax.scatter(
        win_all.loc[~mask_top, "mid_global"],
        win_all.loc[~mask_top, "y_mean"],
        s=base_size,
        c=colors[~mask_top],
        alpha=alphas[~mask_top],
        linewidths=0
    )
    # top-K
    ax.scatter(
        win_all.loc[mask_top, "mid_global"],
        win_all.loc[mask_top, "y_mean"],
        s=40, c="red", alpha=1.0, zorder=3, linewidths=0
    )

    if yscale == "log":
        ax.set_yscale("log")
        ylab = "Windowed TMRCA mean (generations, log10)"
    else:
        ylab = "Windowed TMRCA mean (generations)"

    ax.set_ylabel(ylab)
    ax.set_xlabel("Genomic position (scaffolds concatenated)")
    ax.set_xticks([t for _, t in chrom_ticks])
    ax.set_xticklabels([c for c, _ in chrom_ticks], rotation=90, fontsize=8)
    ax.set_title(
        f"Windowed TMRCA Manhattan plot (win={WIN_BP:,} bp, step={STEP_BP:,} bp)\n"
        f"Top {TOP_K} windows in red{title_suffix}"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)

# --------------------
# 1) Read data
# --------------------
os.makedirs(OUTDIR, exist_ok=True)
cols = ["chrom","start","end","tmrca_mean","tmrca_lo","tmrca_hi"]
tm = pd.read_csv(INPUT_FILE, sep=r"\s+|\t", engine="python", header=None, names=cols)
for c in ["start","end","tmrca_mean","tmrca_lo","tmrca_hi"]:
    tm[c] = pd.to_numeric(tm[c], errors="coerce")
tm = tm.dropna(subset=["chrom","start","end","tmrca_mean"]).reset_index(drop=True)
tm["start_disp"] = tm["start"] + 1 if ZERO_BASED_START else tm["start"]
tm["end_disp"]   = tm["end"]

# --------------------
# 2) Windowing per scaffold (coverage-weighted mean)
# --------------------
# Order scaffolds 1a, 1b, 2, 3, ...
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

# keep only informative windows: coverage >0, mean not NaN, and mean > 0
win_all = win_all[(win_all["bp_covered"] > 0) & (win_all["y_mean"].notna()) & (win_all["y_mean"] > 0)].copy()
if win_all.empty:
    raise RuntimeError("No windows with positive mean and coverage. Adjust WIN_BP/STEP_BP or check input.")

win_all.to_csv(os.path.join(OUTDIR, f"tmrca_windows_{WIN_BP}_{STEP_BP}.csv"), index=False)

# --------------------
# 3) Summaries (unchanged from earlier; saved for completeness)
# --------------------
def summarize_series(x: pd.Series, label: str) -> pd.DataFrame:
    desc = {
        "n": int(x.size),
        "n_nonzero": int((x>0).sum()),
        "min": float(x.min()),
        "q1": float(x.quantile(0.25)),
        "median": float(x.median()),
        "mean": float(x.mean()),
        "q3": float(x.quantile(0.75)),
        "max": float(x.max()),
        "sd": float(x.std())
    }
    return pd.DataFrame([desc], index=[label])

summ_all = summarize_series(tm["tmrca_mean"], "all_values")
summ_no0 = summarize_series(tm.loc[tm["tmrca_mean"]>0, "tmrca_mean"], "nonzero_only")
pd.concat([summ_all, summ_no0]).to_csv(os.path.join(OUTDIR, "summary_global_tmrca.csv"))

# --------------------
# 4) Manhattan plots (two versions)
# --------------------
# LOG scale
make_manhattan(
    win_all, scaffolds_ordered, chrom_ticks,
    out_png=os.path.join(OUTDIR, f"manhattan_tmrca_log_win{WIN_BP}_step{STEP_BP}.png"),
    yscale="log",
    title_suffix=" (log scale)"
)
# LINEAR scale
make_manhattan(
    win_all, scaffolds_ordered, chrom_ticks,
    out_png=os.path.join(OUTDIR, f"manhattan_tmrca_linear_win{WIN_BP}_step{STEP_BP}.png"),
    yscale="linear",
    title_suffix=" (linear scale)"
)

print("Done. Plots written to:", OUTDIR)

