#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm


'''
Run examples
cd /Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025

python GWAS_TMRCA_join_and_plots_v7.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025/To_merge_with_TMRCA_MTAG_SIgnificant_sumstat.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_MTAG_trait_tmrca_results_v7" \
  --allpoints-color "#000000" \
 --regline --show-lm --show-r --show-p \
  --stats-outside --stats-panel --stats-position top --stats-panel-height 1.6 \
  --min-n-trait 5  \
  --top-percentiles 5,10,25

cd /Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files

python GWAS_TMRCA_join_and_plots_v7.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files/combined_pval_maf_hits/combined_hits_full_merge_with_tmrca.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_single_trait_tmrca_results_v7" \
  --allpoints-color "#000000" \
 --regline --show-lm --show-r --show-p \
  --stats-outside --stats-panel --stats-position top --stats-panel-height 1.6 \
  --min-n-trait 5 \
  --top-percentiles 5,10,25
  
'''
# Optional SciPy for Pearson p-values
try:
    from scipy import stats as sstats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

# ============================ utilities ============================

def ensure_dir(p): Path(p).mkdir(parents=True, exist_ok=True)

def pick_col(cols, candidates, required=True, label="column"):
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    if required:
        raise SystemExit(f"ERROR: Could not find {label}. Looked for: {candidates}")
    return None

def merge_asof_per_chrom(gwas, tmrca, scaf_col, pos_col,
                         t_scaf_col, t_start_col, t_end_col):
    """
    Interval-aware join per scaffold.
    If SNP position is NOT inside [start,end], set ALL TMRCA fields to NaN.
    """
    out = []
    for scaf, gsub in gwas.groupby(scaf_col, sort=False):
        tsub = tmrca[tmrca[t_scaf_col] == scaf]
        gsub = gsub.copy()
        if tsub.empty:
            out.append(gsub)
            continue

        gsub[pos_col] = pd.to_numeric(gsub[pos_col], errors="coerce").astype("float64")
        tsub = tsub.copy()
        tsub[t_start_col] = pd.to_numeric(tsub[t_start_col], errors="coerce").astype("float64")
        tsub[t_end_col]   = pd.to_numeric(tsub[t_end_col], errors="coerce").astype("float64")

        gsub.sort_values(pos_col, inplace=True)
        tsub.sort_values(t_start_col, inplace=True)

        merged = pd.merge_asof(
            gsub, tsub,
            left_on=pos_col, right_on=t_start_col,
            direction="backward", allow_exact_matches=True
        )

        inside = merged[pos_col].between(merged[t_start_col], merged[t_end_col], inclusive="both")
        tmrca_cols = [c for c in tsub.columns if c != t_scaf_col]
        merged.loc[~inside, tmrca_cols] = np.nan
        out.append(merged)

    return pd.concat(out, ignore_index=True) if out else gwas.copy()

def trait_summary(df, trait_col, tmrca_col):
    """NaN-aware trait summaries; also report n_total and n_with_tmrca."""
    def one(g):
        x = g[tmrca_col]
        has = x.notna().any()
        return pd.Series({
            "n_total": len(g),
            "n_with_tmrca": x.notna().sum(),
            "median": np.nanmedian(x) if has else np.nan,
            "q1": np.nanpercentile(x, 25) if has else np.nan,
            "q3": np.nanpercentile(x, 75) if has else np.nan,
            "q4": np.nanpercentile(x, 100) if has else np.nan,
            "min": np.nanmin(x) if has else np.nan,
            "max": np.nanmax(x) if has else np.nan,
        })
    return df.groupby(trait_col, dropna=False, observed=True).apply(one).reset_index()

def ols_line(x, y):
    msk = np.isfinite(x) & np.isfinite(y)
    if msk.sum() < 2:
        return None, None
    m, b = np.polyfit(x[msk], y[msk], 1)
    return m, b

def pearson_r_p(x, y):
    msk = np.isfinite(x) & np.isfinite(y)
    if msk.sum() < 2:
        return None, None
    if _HAVE_SCIPY:
        r, p = sstats.pearsonr(x[msk], y[msk])
        return float(r), float(p)
    r = np.corrcoef(x[msk], y[msk])[0, 1]
    return float(r), None

def format_stats_text(x, y, show_lm, show_r, show_p):
    lines = []
    m, b = ols_line(x, y)
    r, p = pearson_r_p(x, y)
    if show_lm and (m is not None):
        lines.append(f"y = {m:.3g} x + {b:.3g}")
    if show_r and (r is not None):
        lines.append(f"r = {r:.3f}")
    if show_p:
        lines.append(f"p = {p:.2e}" if p is not None else "p = (requires SciPy)")
    return "\n".join(lines)

def add_reg_line(ax, x, y, color="black", lw=1.2, alpha=0.8):
    m, b = ols_line(x, y)
    if m is None:
        return
    xs = np.linspace(np.nanmin(x[np.isfinite(x)]), np.nanmax(x[np.isfinite(x)]), 200)
    ax.plot(xs, m*xs + b, color=color, linewidth=lw, alpha=alpha)

def make_axes(figsize=(7,5), stats_strip=False, strip_pos="top", strip_height_in=1.0):
    """
    Create fig/ax, optionally with a separate 'strip' axes for stats text that
    does NOT shrink the plotting area. Strip size is in inches.
    """
    if not stats_strip:
        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax, None

    w_in, h_in = figsize
    if strip_pos in ("top", "bottom"):
        fig_h = h_in + strip_height_in
        fig = plt.figure(figsize=(w_in, fig_h))
        if strip_pos == "top":
            gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[strip_height_in, h_in], hspace=0.02)
            ax_strip = fig.add_subplot(gs[0, 0]); ax_main = fig.add_subplot(gs[1, 0])
        else:
            gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[h_in, strip_height_in], hspace=0.02)
            ax_main = fig.add_subplot(gs[0, 0]); ax_strip = fig.add_subplot(gs[1, 0])
    else:
        fig_w = w_in + strip_height_in
        fig = plt.figure(figsize=(fig_w, h_in))
        gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[w_in, strip_height_in], wspace=0.02)
        ax_main = fig.add_subplot(gs[0, 0]); ax_strip = fig.add_subplot(gs[0, 1])

    ax_strip.axis("off")
    return fig, ax_main, ax_strip

# ============================ plotting helpers ============================

def scatter_basic(x, y, ax, title, xlabel, ylabel,
                  color="black", alpha=0.6, s=10,
                  regline=False, show_lm=False, show_r=False, show_p=False,
                  stats_outside=False, stats_position="top", stats_margin=0.18,
                  stats_panel=False, stats_panel_ax=None):
    ax.scatter(x, y, s=s, c=color, alpha=alpha, edgecolors="none")
    if regline:
        add_reg_line(ax, x, y, color=color)

    if show_lm or show_r or show_p:
        txt = format_stats_text(x, y, show_lm, show_r, show_p)
        if stats_outside and stats_panel and (stats_panel_ax is not None):
            stats_panel_ax.text(0.5, 0.5, txt, ha="center", va="center",
                                fontsize=10, color=color,
                                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=color, alpha=0.95))
        elif stats_outside:
            fig = ax.figure
            margin = float(np.clip(stats_margin, 0.05, 0.35))
            if stats_position == "top":
                plt.subplots_adjust(top=1.0 - margin)
                fig.text(0.5, 1.0 - (margin/2.0), txt, ha="center", va="center",
                         fontsize=9, color=color,
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, alpha=0.9))
            elif stats_position == "bottom":
                plt.subplots_adjust(bottom=margin)
                fig.text(0.5, margin/2.0, txt, ha="center", va="center",
                         fontsize=9, color=color,
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, alpha=0.9))
            else:
                plt.subplots_adjust(right=1.0 - margin)
                fig.text(1.0 - (margin/2.0), 0.5, txt, ha="center", va="center",
                         fontsize=9, color=color,
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, alpha=0.9))
        else:
            ax.text(0.02, 0.98, txt, transform=ax.transAxes,
                    ha="left", va="top", fontsize=9, color=color,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.75))

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.15)

def scatter_by_trait(df, xcol, ycol, trait_col, ax, title,
                     regline=False, show_lm=False, show_r=False, show_p=False,
                     per_trait_stats=False,
                     stats_outside=False, stats_position="top", stats_margin=0.18,
                     stats_panel=False, stats_panel_ax=None):
    traits = df[trait_col].astype(str).fillna("NA").unique().tolist()
    cmap = cm.get_cmap("tab20", max(10, len(traits)))

    # Overall stats across filtered subset
    if (show_lm or show_r or show_p) and not per_trait_stats and not df.empty:
        txt = format_stats_text(df[xcol].to_numpy(), df[ycol].to_numpy(), show_lm, show_r, show_p)
        if stats_outside and stats_panel and (stats_panel_ax is not None):
            stats_panel_ax.text(0.5, 0.5, txt, ha="center", va="center",
                                fontsize=10, color="black",
                                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="black", alpha=0.95))
        elif stats_outside:
            fig = ax.figure
            margin = float(np.clip(stats_margin, 0.05, 0.35))
            if stats_position == "top":
                plt.subplots_adjust(top=1.0 - margin)
                fig.text(0.5, 1.0 - (margin/2.0), txt, ha="center", va="center",
                         fontsize=9, color="black",
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9))
            elif stats_position == "bottom":
                plt.subplots_adjust(bottom=margin)
                fig.text(0.5, margin/2.0, txt, ha="center", va="center",
                         fontsize=9, color="black",
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9))
            else:
                plt.subplots_adjust(right=1.0 - margin)
                fig.text(1.0 - (margin/2.0), 0.5, txt, ha="center", va="center",
                         fontsize=9, color="black",
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", alpha=0.9))
        else:
            ax.text(0.02, 0.98, txt, transform=ax.transAxes, ha="left", va="top",
                    fontsize=9, color="black",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.75))

    for i, t in enumerate(traits):
        sub = df[df[trait_col].astype(str).fillna("NA") == t]
        col = cmap(i)
        ax.scatter(sub[xcol], sub[ycol], s=10, alpha=0.7, color=col, label=t, edgecolors="none")
        if regline:
            add_reg_line(ax, sub[xcol].to_numpy(), sub[ycol].to_numpy(), color=col)

        if per_trait_stats and (show_lm or show_r or show_p):
            xv = sub[xcol].to_numpy(); yv = sub[ycol].to_numpy()
            m, b = ols_line(xv, yv); r, p = pearson_r_p(xv, yv)
            bits = []
            if show_lm and (m is not None): bits.append(f"m={m:.2g}")
            if show_r and (r is not None):  bits.append(f"r={r:.2f}")
            if show_p: bits.append(f"p={p:.1e}" if p is not None else "p=?")
            if bits: ax.scatter([], [], color=col, label=f"{t} ({', '.join(bits)})")

    ax.set_title(title)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.legend(markerscale=1.5, fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.15)

def filter_by_min_points(df, trait_col, xcol, ycol, min_n):
    if min_n <= 0 or df.empty:
        return df, pd.Series(dtype=int)
    dvalid = df[[trait_col, xcol, ycol]].dropna()
    counts = dvalid.groupby(trait_col, observed=True).size().sort_values(ascending=False)
    keep = counts[counts >= min_n].index
    return df[df[trait_col].isin(keep)], counts

# ============================ new: percent-in-top-X ============================

def compute_top_thresholds(tmrca_means, percents):
    """
    percents: iterable of P (e.g., 5, 10) meaning top P percent.
    Returns dict {P: cutoff_value} where cutoff = quantile at 1 - P/100.
    """
    tm = pd.to_numeric(tmrca_means, errors="coerce").dropna().to_numpy()
    if tm.size == 0:
        return {p: np.nan for p in percents}
    return {p: float(np.quantile(tm, 1 - p/100.0)) for p in percents}

def summarize_snps_in_top(annotated, trait_col, mean_col, cutoffs):
    """
    cutoffs: dict {P: cutoff_value}
    Returns:
      overall_df, by_trait_df, annotated_with_flags
    """
    rows_overall = []
    rows_trait = []

    n_total = len(annotated)
    n_with  = annotated[mean_col].notna().sum()

    for P, cutoff in cutoffs.items():
        flag = annotated[mean_col] >= cutoff
        flag = flag.fillna(False)

        nin_top = int(flag.sum())
        pct_all = (nin_top / n_total * 100.0) if n_total > 0 else np.nan
        pct_with = (nin_top / n_with * 100.0) if n_with > 0 else np.nan

        rows_overall.append({
            "top_percent": P,
            "cutoff_mean_tmrca": cutoff,
            "n_snps_total": n_total,
            "n_snps_with_tmrca": int(n_with),
            "n_snps_in_top": nin_top,
            "pct_of_all_snps": pct_all,
            "pct_of_snps_with_tmrca": pct_with,
        })

        for tr, g in annotated.groupby(trait_col, dropna=False, observed=True):
            n_t_total = len(g)
            n_t_with  = g[mean_col].notna().sum()
            f2 = g[mean_col].ge(cutoff).fillna(False)
            n_t_top = int(f2.sum())
            rows_trait.append({
                "trait": tr,
                "top_percent": P,
                "cutoff_mean_tmrca": cutoff,
                "n_snps_total": n_t_total,
                "n_snps_with_tmrca": int(n_t_with),
                "n_snps_in_top": n_t_top,
                "pct_of_all_snps": (n_t_top / n_t_total * 100.0) if n_t_total > 0 else np.nan,
                "pct_of_snps_with_tmrca": (n_t_top / n_t_with * 100.0) if n_t_with > 0 else np.nan,
            })

        colname = f"in_top_{int(P)}pct"
        annotated[colname] = flag

    overall_df = pd.DataFrame(rows_overall).sort_values("top_percent")
    by_trait_df = pd.DataFrame(rows_trait).sort_values(["top_percent","trait"])
    return overall_df, by_trait_df, annotated

# ============================ main ============================

def main():
    ap = argparse.ArgumentParser(
        description="Annotate GWAS with TMRCA & summarize/plot (v7: TMRCA on x-axis; top-percentile coverage)."
    )
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--tmrca", required=True)
    ap.add_argument("--outdir", default="gwas_tmrca_out")

    # Columns
    ap.add_argument("--trait-col", default="Significant_Trait")
    ap.add_argument("--scaf-col", default="Scaffold")
    ap.add_argument("--pos-col",  default="position", type=str)
    ap.add_argument("--beta-col", default="BETA")
    ap.add_argument("--maf-col",  default="MAF")

    # Figure sizing
    ap.add_argument("--fig-w", type=float, default=7.0)
    ap.add_argument("--fig-h", type=float, default=5.0)

    # Regression & stats toggles
    ap.add_argument("--regline", action="store_true")
    ap.add_argument("--show-lm", action="store_true")
    ap.add_argument("--show-r", action="store_true")
    ap.add_argument("--show-p", action="store_true")
    ap.add_argument("--per-trait-stats", action="store_true")

    # Outside stats + strip panel
    ap.add_argument("--stats-outside", action="store_true")
    ap.add_argument("--stats-position", default="top", choices=["top","bottom","right"])
    ap.add_argument("--stats-margin", type=float, default=0.18)
    ap.add_argument("--stats-panel", action="store_true")
    ap.add_argument("--stats-panel-height", type=float, default=1.0)

    # Appearance
    ap.add_argument("--allpoints-color", default="#000000")
    ap.add_argument("--min-n-trait", type=int, default=0)
    ap.add_argument("--dpi", type=int, default=300)

    # NEW: thresholds
    ap.add_argument("--top-percentiles", default="5",
                    help="Comma-separated list of percentages (e.g., '5,10,25').")

    args = ap.parse_args()
    ensure_dir(args.outdir)

    # ----------- load
    gwas = pd.read_csv(args.gwas, sep="\t", dtype={})
    tmrca = pd.read_csv(args.tmrca, sep="\t", dtype={})

    gcols = gwas.columns.tolist()
    tcols = tmrca.columns.tolist()

    c_trait = pick_col(gcols, [args.trait_col], label="trait")
    c_scaf  = pick_col(gcols, [args.scaf_col], label="scaffold (GWAS)")
    c_pos   = pick_col(gcols, [args.pos_col], label="position (GWAS)")
    c_beta  = pick_col(gcols, [args.beta_col], label="BETA")
    c_maf   = pick_col(gcols, [args.maf_col], label="MAF")

    t_scaf  = pick_col(tcols, ["CHROM","scaffold","seqid"], label="scaffold (TMRCA)")
    t_start = pick_col(tcols, ["start_tmrca","start","seg_start"], label="start_tmrca")
    t_end   = pick_col(tcols, ["end_tmrca","end","seg_end"], label="end_tmrca")
    t_mean  = pick_col(tcols, ["mean_tmrca","tmrca_mean"], label="mean_tmrca")
    t_lo    = pick_col(tcols, ["Low_CI","low_ci","lower_ci"], label="Low_CI")
    t_hi    = pick_col(tcols, ["UP_CI","up_ci","upper_ci"], label="UP_CI")
    t_len   = pick_col(tcols, ["seg_length","length"], required=False, label="seg_length")

    for c in [c_pos, c_beta, c_maf]: gwas[c] = pd.to_numeric(gwas[c], errors="coerce")
    for c in [t_start, t_end, t_mean, t_lo, t_hi]: tmrca[c] = pd.to_numeric(tmrca[c], errors="coerce")
    if t_len: tmrca[t_len] = pd.to_numeric(tmrca[t_len], errors="coerce")

    gwas = gwas.dropna(subset=[c_scaf, c_pos])

    keep_cols = [t_scaf, t_start, t_end, t_mean, t_lo, t_hi] + ([t_len] if t_len else [])
    annotated = merge_asof_per_chrom(
        gwas, tmrca[keep_cols].copy(),
        scaf_col=c_scaf, pos_col=c_pos,
        t_scaf_col=t_scaf, t_start_col=t_start, t_end_col=t_end
    )

    rename_map = {t_scaf:"tmrca_scaffold", t_start:"start_tmrca", t_end:"end_tmrca",
                  t_mean:"mean_tmrca", t_lo:"Low_CI", t_hi:"UP_CI"}
    if t_len: rename_map[t_len] = "seg_length"
    annotated = annotated.rename(columns=rename_map)

    # ===================== NEW: percent of SNPs in top-X% =====================
    try:
        top_list = [float(x.strip()) for x in str(args.top_percentiles).split(",") if x.strip() != ""]
    except Exception:
        raise SystemExit("ERROR: --top-percentiles must be numbers like '5' or '5,10,25'.")
    if any((p <= 0 or p >= 100) for p in top_list):
        raise SystemExit("ERROR: --top-percentiles values must be >0 and <100.")

    cutoffs = compute_top_thresholds(tmrca_means=tmrca[t_mean], percents=top_list)
    overall_top, by_trait_top, annotated = summarize_snps_in_top(
        annotated, trait_col=c_trait, mean_col="mean_tmrca", cutoffs=cutoffs
    )

    overall_path = Path(args.outdir, "snps_in_top_tmrca_percentiles_overall.tsv")
    bytrait_path = Path(args.outdir, "snps_in_top_tmrca_percentiles_by_trait.tsv")
    overall_top.to_csv(overall_path, sep="\t", index=False)
    by_trait_top.to_csv(bytrait_path, sep="\t", index=False)

    # ===================== main outputs =====================
    out_annot = Path(args.outdir) / "gwas_with_tmrca.tsv"
    annotated.to_csv(out_annot, sep="\t", index=False)

    trait_summ = trait_summary(annotated, trait_col=c_trait, tmrca_col="mean_tmrca")
    out_summ = Path(args.outdir) / "tmrca_trait_summary.tsv"
    trait_summ.to_csv(out_summ, sep="\t", index=False)

    def savefig(fig, name, dpi=300):
        fig.tight_layout()
        fig.savefig(Path(args.outdir, f"{name}.png"), dpi=dpi)
        fig.savefig(Path(args.outdir, f"{name}.pdf"))

    base_size = (args.fig_w, args.fig_h)

    # -------- BETA vs TMRCA (TMRCA on x axis) – all points
    dfp = annotated.dropna(subset=["mean_tmrca", c_beta]).copy()
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_basic(dfp["mean_tmrca"].to_numpy(), dfp[c_beta].to_numpy(), ax,
                  title="mean TMRCA vs BETA",
                  xlabel="mean TMRCA", ylabel="BETA",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position,
                  stats_margin=args.stats_margin,
                  stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_BETA_all", dpi=args.dpi)

    # -------- |BETA| vs TMRCA (TMRCA on x axis) – all points
    dfp["absBETA"] = dfp[c_beta].abs()
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_basic(dfp["mean_tmrca"].to_numpy(), dfp["absBETA"].to_numpy(), ax,
                  title="mean TMRCA vs |BETA|",
                  xlabel="mean TMRCA", ylabel="|BETA|",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position,
                  stats_margin=args.stats_margin,
                  stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_absBETA_all", dpi=args.dpi)

    # -------- TMRCA vs MAF – all points (TMRCA on x axis for consistency)
    dfm = annotated.dropna(subset=["mean_tmrca", c_maf]).copy()
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_basic(dfm["mean_tmrca"].to_numpy(), dfm[c_maf].to_numpy(), ax,
                  title="mean TMRCA vs MAF",
                  xlabel="mean TMRCA", ylabel="MAF",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position,
                  stats_margin=args.stats_margin,
                  stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_MAF_all", dpi=args.dpi)

    # -------- BY-TRAIT PLOTS (TMRCA on x axis)

    # mean TMRCA vs BETA by trait
    dfc = annotated.dropna(subset=["mean_tmrca", c_beta]).copy()
    dfc_bt, counts_bt = filter_by_min_points(dfc, c_trait, "mean_tmrca", c_beta, args.min_n_trait)
    if args.min_n_trait > 0:
        counts_bt.to_csv(Path(args.outdir, "counts_TMRCA_vs_BETA_byTrait.tsv"),
                         sep="\t", header=["n"])
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_by_trait(dfc_bt, xcol="mean_tmrca", ycol=c_beta,
                     trait_col=c_trait, ax=ax,
                     title=f"mean TMRCA vs BETA (by trait; min_n={args.min_n_trait})",
                     regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position,
                     stats_margin=args.stats_margin,
                     stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_BETA_byTrait", dpi=args.dpi)

    # mean TMRCA vs |BETA| by trait
    dfc_abs = dfc.copy()
    dfc_abs["absBETA"] = dfc_abs[c_beta].abs()
    dfc_abs_f, counts_abt = filter_by_min_points(dfc_abs, c_trait, "mean_tmrca", "absBETA", args.min_n_trait)
    if args.min_n_trait > 0:
        counts_abt.to_csv(Path(args.outdir, "counts_TMRCA_vs_absBETA_byTrait.tsv"),
                          sep="\t", header=["n"])
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_by_trait(dfc_abs_f, xcol="mean_tmrca", ycol="absBETA",
                     trait_col=c_trait, ax=ax,
                     title=f"mean TMRCA vs |BETA| (by trait; min_n={args.min_n_trait})",
                     regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position,
                     stats_margin=args.stats_margin,
                     stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_absBETA_byTrait", dpi=args.dpi)

    # mean TMRCA vs MAF by trait
    dfc_maf = annotated.dropna(subset=["mean_tmrca", c_maf]).copy()
    dfc_maf_f, counts_maf = filter_by_min_points(dfc_maf, c_trait, "mean_tmrca", c_maf, args.min_n_trait)
    if args.min_n_trait > 0:
        counts_maf.to_csv(Path(args.outdir, "counts_TMRCA_vs_MAF_byTrait.tsv"),
                          sep="\t", header=["n"])
    fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                  strip_pos=args.stats_position,
                                  strip_height_in=args.stats_panel_height)
    scatter_by_trait(dfc_maf_f, xcol="mean_tmrca", ycol=c_maf,
                     trait_col=c_trait, ax=ax,
                     title=f"mean TMRCA vs MAF (by trait; min_n={args.min_n_trait})",
                     regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position,
                     stats_margin=args.stats_margin,
                     stats_panel=args.stats_panel, stats_panel_ax=ax_strip)
    savefig(fig, "TMRCA_vs_MAF_byTrait", dpi=args.dpi)

    print("[OK] Wrote:")
    print(f"  {out_annot}")
    print(f"  {overall_path}")
    print(f"  {bytrait_path}")
    print(f"  {out_summ}")
    print(f"  figures & counts under {args.outdir}/")

if __name__ == "__main__":
    main()
