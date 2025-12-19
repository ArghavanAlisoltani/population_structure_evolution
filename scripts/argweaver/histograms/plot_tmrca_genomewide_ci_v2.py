#!/usr/bin/env python3
"""
V2 — Genome-wide TMRCA plotting with clean segment drawing (no false links)
and separate RAW vs WINDOWED outputs (each with line/CI and bar options).
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

# --------------------- helpers ---------------------

def pick_col(cols, candidates, required=True, label="column"):
    lowmap = {c.lower(): c for c in cols}
    for name in candidates:
        if name.lower() in lowmap:
            return lowmap[name.lower()]
    if required:
        raise SystemExit(f"ERROR: Could not find {label}. Looked for: {candidates}")
    return None

def parse_scaffold_key(name: str):
    import re
    m = re.search(r'(\d+)', name)
    if m:
        num = int(m.group(1))
        suf = name[m.end():]
        return (num, suf)
    return (10**9, name)

def gb_format(x, pos):
    return f"{x/1e9:.2f}"

def compute_window_stats(sub_df, start_col, end_col, val_col, win_bp, sc_len):
    n = max(1, int(np.ceil(sc_len / win_bp)))
    starts = np.arange(n, dtype=np.int64) * win_bp
    ends = np.minimum(starts + win_bp, sc_len)
    mids = (starts + ends) / 2.0

    num = np.zeros(n, dtype=np.float64)
    den = np.zeros(n, dtype=np.float64)

    s = sub_df[start_col].to_numpy(np.int64)
    e = sub_df[end_col].to_numpy(np.int64)
    v = sub_df[val_col].to_numpy(np.float64)

    for ss, ee, vv in zip(s, e, v):
        if ee <= ss:
            continue
        w0 = ss // win_bp
        w1 = (ee - 1) // win_bp
        for wi in range(int(w0), int(w1) + 1):
            left = max(ss, wi * win_bp)
            right = min(ee, (wi + 1) * win_bp)
            ol = right - left
            if ol > 0:
                num[wi] += vv * ol
                den[wi] += ol

    with np.errstate(divide="ignore", invalid="ignore"):
        wmean = np.where(den > 0, num / den, np.nan)

    return mids, wmean, starts, ends

def add_scaffold_axis(ax, scaf_df, x_mode):
    if x_mode == "gb":
        ax.xaxis.set_major_formatter(FuncFormatter(gb_format))
        ax.set_xlabel("Genome position (Gb)")
    else:
        ax.set_xticks(scaf_df["mid_cum"].to_numpy())
        ax.set_xticklabels(scaf_df["scaffold"].astype(str).tolist(), rotation=90)
        ax.set_xlabel("Scaffolds (ordered: 1a, 1b, 2, …)")

def draw_boundaries(ax, scaf_df):
    for b in scaf_df["cum_offset"].to_numpy():
        if b > 0:
            ax.axvline(b, color="#E6E6E6", linewidth=0.5, zorder=0)

# --------------------- drawing primitives ---------------------

def plot_raw_line(ax, df, c_scaf, c_start, c_end, c_mean,
                  line_color="black", line_width=0.6,
                  max_segments=None, min_seg_len_bp=0):
    segs = []
    counter = 0
    for scaf, sub in df.sort_values([c_scaf, c_start]).groupby(c_scaf, observed=True):
        if min_seg_len_bp > 0:
            sub = sub[(sub[c_end] - sub[c_start]) >= min_seg_len_bp]
        for _, r in sub.iterrows():
            segs.append([(r["start_cum"], r[c_mean]), (r["end_cum"], r[c_mean])])
            counter += 1
            if max_segments and counter >= max_segments:
                break
        if max_segments and counter >= max_segments:
            break
    if segs:
        lc = LineCollection(segs, colors=line_color, linewidths=line_width)
        ax.add_collection(lc)
        ax.set_xlim(df["start_cum"].min(), df["end_cum"].max())

def plot_raw_bars(ax, df, c_start, c_end, c_mean,
                  face_color="black", alpha=0.8,
                  min_seg_len_bp=0, max_bars=None):
    patches = []
    count = 0
    for _, r in df.sort_values(c_start).iterrows():
        w = r["end_cum"] - r["start_cum"]
        if w <= 0 or w < min_seg_len_bp:
            continue
        rect = Rectangle((r["start_cum"], 0), width=w, height=r[c_mean])
        patches.append(rect)
        count += 1
        if max_bars and count >= max_bars:
            break
    if patches:
        pc = PatchCollection(patches, facecolor=face_color, alpha=alpha, edgecolor="none")
        ax.add_collection(pc)

def build_window_series(df, scaf_df, c_scaf, c_start, c_end, c_mean, c_lo, c_hi,
                        win_bp, ci_agg="mean"):
    xs, ys, los, his = [], [], [], []
    for _, row in scaf_df.iterrows():
        scaf = row["scaffold"]
        sc_len = int(row["scaffold_len"])
        offset = float(row["cum_offset"])
        sub = df[df[c_scaf] == scaf].sort_values(c_start)
        if sub.empty:
            xs += [np.array([np.nan])]; ys += [np.array([np.nan])]
            los += [np.array([np.nan])]; his += [np.array([np.nan])]
            continue

        mids, wmean, _, _ = compute_window_stats(sub, c_start, c_end, c_mean, win_bp, sc_len)

        if ci_agg == "mean":
            _, wlo, _, _ = compute_window_stats(sub, c_start, c_end, c_lo, win_bp, sc_len)
            _, whi, _, _ = compute_window_stats(sub, c_start, c_end, c_hi, win_bp, sc_len)
        else:
            # min/max across contributing segments per window
            n = len(mids)
            wlo = np.full(n, np.nan, dtype=float)
            whi = np.full(n, np.nan, dtype=float)
            s = sub[c_start].to_numpy(np.int64)
            e = sub[c_end].to_numpy(np.int64)
            lo = sub[c_lo].to_numpy(np.float64)
            hi = sub[c_hi].to_numpy(np.float64)
            for ss, ee, l, h in zip(s, e, lo, hi):
                if ee <= ss:
                    continue
                w0 = int(ss // win_bp)
                w1 = int((e - 1) // win_bp)
                rng = slice(w0, w1 + 1)
                # initialize if NaN, else update
                if np.any(np.isnan(wlo[rng])):
                    wlo[rng] = np.where(np.isnan(wlo[rng]), l, np.minimum(wlo[rng], l))
                else:
                    wlo[rng] = np.minimum(wlo[rng], l)
                if np.any(np.isnan(whi[rng])):
                    whi[rng] = np.where(np.isnan(whi[rng]), h, np.maximum(whi[rng], h))
                else:
                    whi[rng] = np.maximum(whi[rng], h)

        xs.append(offset + mids); ys.append(wmean); los.append(wlo); his.append(whi)
        xs.append(np.array([np.nan])); ys.append(np.array([np.nan]))
        los.append(np.array([np.nan])); his.append(np.array([np.nan]))

    x = np.concatenate(xs); y = np.concatenate(ys)
    lo = np.concatenate(los); hi = np.concatenate(his)
    return x, y, lo, hi

# --------------------- main ---------------------

def main():
    ap = argparse.ArgumentParser(description="Genome-wide TMRCA plotter (V2, fixed)")
    ap.add_argument("--tmrca", required=True, help="Input TSV with segments & CIs")
    ap.add_argument("--out-prefix", default="tmrca_v2", help="Output prefix")

    # Axis + theme
    ap.add_argument("--x-mode", choices=["gb", "scaffold"], default="gb")
    ap.add_argument("--ylog", action="store_true")
    ap.add_argument("--fig-width", type=float, default=14.0)
    ap.add_argument("--fig-height", type=float, default=5.0)
    ap.add_argument("--dpi", type=int, default=300)

    # RAW options
    ap.add_argument("--make-raw-line", action="store_true")
    ap.add_argument("--make-raw-bar", action="store_true")
    ap.add_argument("--line-color", default="black")
    ap.add_argument("--line-width", type=float, default=0.6)
    ap.add_argument("--min-seg-len-bp", type=int, default=0)
    ap.add_argument("--max-raw-segs", type=int, default=None)
    ap.add_argument("--max-raw-bars", type=int, default=None)
    ap.add_argument("--bar-alpha", type=float, default=0.85)
    ap.add_argument("--bar-color", default="black")

    # WINDOW options
    ap.add_argument("--window-bp", type=int, default=0)
    ap.add_argument("--make-window-line", action="store_true")
    ap.add_argument("--make-window-bar", action="store_true")
    ap.add_argument("--ci-agg", choices=["mean", "range"], default="mean")
    ap.add_argument("--ci-alpha", type=float, default=0.35)
    ap.add_argument("--ci-color", default="lightgray")

    args = ap.parse_args()

    # ---------- read & standardize ----------
    df = pd.read_csv(args.tmrca, sep="\t", dtype={})
    cols = df.columns.tolist()

    c_scaf = pick_col(cols, ["CHROM","scaffold","seqid"], label="scaffold")
    c_start = pick_col(cols, ["start_tmrca","start","seg_start"], label="start")
    c_end   = pick_col(cols, ["end_tmrca","end","seg_end"], label="end")
    c_mean  = pick_col(cols, ["mean_tmrca","tmrca_mean"], label="mean_tmrca")
    c_lo    = pick_col(cols, ["Low_CI","low_ci","lower_ci"], label="Low_CI")
    c_hi    = pick_col(cols, ["UP_CI","up_ci","upper_ci","hi_ci"], label="UP_CI")

    for c in [c_start, c_end, c_mean, c_lo, c_hi]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=[c_scaf, c_start, c_end, c_mean])

    # scaffold ordering & cumulative coords
    scaf_len = df.groupby(c_scaf)[c_end].max().rename("scaffold_len").reset_index()
    scaf_len["order_key"] = scaf_len[c_scaf].map(parse_scaffold_key)
    scaf_len = scaf_len.sort_values("order_key").drop(columns="order_key")
    scaf_order = scaf_len[c_scaf].tolist()

    scaf_len["cum_offset"] = scaf_len["scaffold_len"].shift(fill_value=0).cumsum().shift(fill_value=0)
    scaf_len["mid_cum"] = scaf_len["cum_offset"] + scaf_len["scaffold_len"]/2.0
    scaf_len = scaf_len.rename(columns={c_scaf: "scaffold"})

    df[c_scaf] = pd.Categorical(df[c_scaf], categories=scaf_order, ordered=True)
    df = df.merge(scaf_len[["scaffold","cum_offset"]], left_on=c_scaf, right_on="scaffold", how="left")
    df["start_cum"] = df["cum_offset"] + df[c_start]
    df["end_cum"]   = df["cum_offset"] + df[c_end]

    # ---------- RAW: line ----------
    if args.make_raw_line:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        plot_raw_line(ax, df, c_scaf, c_start, c_end, c_mean,
                      line_color=args.line_color, line_width=args.line_width,
                      max_segments=args.max_raw_segs, min_seg_len_bp=args.min_seg_len_bp)
        draw_boundaries(ax, scaf_len)
        add_scaffold_axis(ax, scaf_len, args.x_mode)
        ax.set_ylabel("TMRCA (generations)")
        if args.ylog:
            ax.set_yscale("log")
        ax.set_title("Genome-wide TMRCA — raw segments (line)")
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}_raw_line.png", dpi=args.dpi)
        fig.savefig(f"{args.out_prefix}_raw_line.pdf")

    # ---------- RAW: bars ----------
    if args.make_raw_bar:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        plot_raw_bars(ax, df, c_start, c_end, c_mean,
                      face_color=args.bar_color, alpha=args.bar_alpha,
                      min_seg_len_bp=args.min_seg_len_bp, max_bars=args.max_raw_bars)
        draw_boundaries(ax, scaf_len)
        add_scaffold_axis(ax, scaf_len, args.x_mode)
        ax.set_ylabel("TMRCA (generations)")
        if args.ylog:
            ax.set_yscale("log")
        ax.set_title("Genome-wide TMRCA — raw segments (bars)")
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}_raw_bar.png", dpi=args.dpi)
        fig.savefig(f"{args.out_prefix}_raw_bar.pdf")

    # ---------- WINDOWED series ----------
    need_window = args.make_window_line or args.make_window_bar
    if need_window:
        if args.window_bp <= 0:
            raise SystemExit("ERROR: --window-bp must be > 0 for windowed plots.")
        x, y, lo, hi = build_window_series(
            df, scaf_len, c_scaf, c_start, c_end, c_mean, c_lo, c_hi,
            win_bp=args.window_bp, ci_agg=args.ci_agg
        )

    # ---------- WINDOW: line + CI ----------
    if args.make_window_line:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        ax.fill_between(x, lo, hi, where=~np.isnan(x), color=args.ci_color, alpha=args.ci_alpha)
        ax.plot(x, y, color=args.line_color, linewidth=args.line_width)
        draw_boundaries(ax, scaf_len)
        add_scaffold_axis(ax, scaf_len, args.x_mode)
        ax.set_ylabel("TMRCA (generations)")
        if args.ylog:
            ax.set_yscale("log")
        ax.set_title(f"Genome-wide TMRCA — windowed ({args.window_bp:,} bp) line + CI")
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}_win{args.window_bp}_line.png", dpi=args.dpi)
        fig.savefig(f"{args.out_prefix}_win{args.window_bp}_line.pdf")

    # ---------- WINDOW: bars ----------
    if args.make_window_bar:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        patches = []
        for _, row in scaf_len.iterrows():
            scaf = row["scaffold"]; sc_len = int(row["scaffold_len"]); offset = float(row["cum_offset"])
            sub = df[df[c_scaf] == scaf]
            if sub.empty:
                continue
            mids, wmean, wstarts, wends = compute_window_stats(sub, c_start, c_end, c_mean, args.window_bp, sc_len)
            for xm, ym, ws, we in zip(offset + mids, wmean, wstarts, wends):
                if np.isnan(ym):
                    continue
                rect = Rectangle((offset + ws, 0), width=(we - ws), height=ym)
                patches.append(rect)
        if patches:
            pc = PatchCollection(patches, facecolor=args.bar_color, alpha=args.bar_alpha, edgecolor="none")
            ax.add_collection(pc)

        draw_boundaries(ax, scaf_len)
        add_scaffold_axis(ax, scaf_len, args.x_mode)
        ax.set_ylabel("TMRCA (generations)")
        if args.ylog:
            ax.set_yscale("log")
        ax.set_title(f"Genome-wide TMRCA — windowed ({args.window_bp:,} bp) bars")
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}_win{args.window_bp}_bar.png", dpi=args.dpi)
        fig.savefig(f"{args.out_prefix}_win{args.window_bp}_bar.pdf")

    print("[OK] Done.")

if __name__ == "__main__":
    main()

