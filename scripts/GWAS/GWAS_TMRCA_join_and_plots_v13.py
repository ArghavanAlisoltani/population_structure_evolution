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

python GWAS_TMRCA_join_and_plots_v13.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025/To_merge_with_TMRCA_MTAG_SIgnificant_sumstat.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_MTAG_trait_tmrca_results_v13" \
  --regline --show-lm --show-r --show-p \
  --stats-outside --stats-panel --stats-position top --stats-panel-height 1.6 \
  --min-n-trait 5 \
  --exclude-traits WWD \
  --top-percentiles 5,10,25 \
  --highlight-traits C13 \
  --trait-alpha-default 0.2 \
  --trait-alpha-highlight 0.9 \
  --trait-colors "HT:red,LDECL:orange,C13:green" \
  --stats-traits C13 \         
  --stats-trait-xq 50 \
  --stats-trait-yq 98 \
  --tmrca-regfilter 0 \
  --tmrca-log10

###two trait highlighted in MTAG
python GWAS_TMRCA_join_and_plots_v10.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025/To_merge_with_TMRCA_MTAG_SIgnificant_sumstat.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_MTAG_trait_tmrca_results_v10" \
  --regline --show-lm --show-r --show-p \
  --stats-outside --stats-panel --stats-position top --stats-panel-height 1.6 \
  --min-n-trait 5 \
  --top-percentiles 5,10,25 \
  --highlight-traits C13 \
  --trait-alpha-default 1 \
  --trait-alpha-highlight 1 \
  --trait-colors "HT:red,LDECL:orange,C13:green" \
  --stats-traits "C13, LDECL" \
  --exclude-traits WWD \
  --stats-trait-xq 60 \
  --stats-trait-yq 85

### single trait

cd /Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025

python GWAS_TMRCA_join_and_plots_v13.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files/combined_pval_maf_hits/combined_hits_full_merge_with_tmrca.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_single_trait_tmrca_results_v13_filtered_regression" \
  --regline --show-lm --show-r --show-p \
  --stats-outside --stats-panel --stats-position top --stats-panel-height 1.6 \
  --min-n-trait 5 \
  --top-percentiles 5,10,25 \
  --tmrca-regfilter 0 \    
  --highlight-traits HT30 \  
  --trait-alpha-default 0.4 \
  --trait-alpha-highlight 1 \                       
  --trait-colors "HT30:red,LDECL:orange,C13:green" \
  --stats-traits HT30 \
  --stats-trait-xq 89 \
  --stats-trait-yq 29 \
  --tmrca-log10

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

        inside = merged[pos_col].between(merged[t_start_col],
                                         merged[t_end_col],
                                         inclusive="both")
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
    r = np.corrcoef(x[msk])[0, 1]
    return float(r), None

def format_stats_text(x, y, show_lm, show_r, show_p):
    lines = []
    if x.size < 2:
        return ""
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
    if x.size < 2:
        return
    m, b = ols_line(x, y)
    if m is None:
        return
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
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
            gs = fig.add_gridspec(nrows=2, ncols=1,
                                  height_ratios=[strip_height_in, h_in],
                                  hspace=0.02)
            ax_strip = fig.add_subplot(gs[0, 0])
            ax_main = fig.add_subplot(gs[1, 0])
        else:
            gs = fig.add_gridspec(nrows=2, ncols=1,
                                  height_ratios=[h_in, strip_height_in],
                                  hspace=0.02)
            ax_main = fig.add_subplot(gs[0, 0])
            ax_strip = fig.add_subplot(gs[1, 0])
    else:
        fig_w = w_in + strip_height_in
        fig = plt.figure(figsize=(fig_w, h_in))
        gs = fig.add_gridspec(nrows=1, ncols=2,
                              width_ratios=[w_in, strip_height_in],
                              wspace=0.02)
        ax_main = fig.add_subplot(gs[0, 0])
        ax_strip = fig.add_subplot(gs[0, 1])

    ax_strip.axis("off")
    return fig, ax_main, ax_strip

def parse_trait_list(s):
    return [x.strip() for x in str(s).split(",") if x.strip() != ""]

def parse_trait_colors(s):
    """
    'HT:red,C13:green,WWD:blue' -> {'HT':'red', 'C13':'green', 'WWD':'blue'}
    """
    out = {}
    s = str(s).strip()
    if not s:
        return out
    for part in s.split(","):
        if ":" not in part:
            continue
        k, v = part.split(":", 1)
        k = k.strip()
        v = v.strip()
        if k and v:
            out[k] = v
    return out

# ============================ top-percentile summaries ============================

def compute_top_thresholds(tmrca_means, percents):
    tm = pd.to_numeric(tmrca_means, errors="coerce").dropna().to_numpy()
    if tm.size == 0:
        return {p: np.nan for p in percents}
    return {p: float(np.quantile(tm, 1 - p/100.0)) for p in percents}

def summarize_snps_in_top(annotated, trait_col, mean_col, cutoffs):
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
            "cutoff_tmrca": cutoff,
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
                "cutoff_tmrca": cutoff,
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

# ============================ low-TMRCA (below global mean) ============================

def summarize_snps_below_mean(annotated, trait_col, mean_col, global_mean):
    col_flag = "below_global_mean_tmrca"

    tm = annotated[mean_col]
    flag = tm < global_mean
    flag = flag & tm.notna()
    annotated[col_flag] = flag

    n_total = len(annotated)
    n_with  = tm.notna().sum()
    n_below = int(flag.sum())

    rows_overall = [{
        "threshold_type": "global_mean",
        "threshold_tmrca": global_mean,
        "n_snps_total": n_total,
        "n_snps_with_tmrca": int(n_with),
        "n_snps_below": n_below,
        "pct_of_all_snps": (n_below / n_total * 100.0) if n_total > 0 else np.nan,
        "pct_of_snps_with_tmrca": (n_below / n_with * 100.0) if n_with > 0 else np.nan,
    }]

    rows_trait = []
    for tr, g in annotated.groupby(trait_col, dropna=False, observed=True):
        tm_t = g[mean_col]
        flag_t = g[col_flag]
        n_t_total = len(g)
        n_t_with  = tm_t.notna().sum()
        n_t_below = int(flag_t.sum())
        rows_trait.append({
            "trait": tr,
            "threshold_type": "global_mean",
            "threshold_tmrca": global_mean,
            "n_snps_total": n_t_total,
            "n_snps_with_tmrca": int(n_t_with),
            "n_snps_below": n_t_below,
            "pct_of_all_snps": (n_t_below / n_t_total * 100.0) if n_t_total > 0 else np.nan,
            "pct_of_snps_with_tmrca": (n_t_below / n_t_with * 100.0) if n_t_with > 0 else np.nan,
        })

    overall_df = pd.DataFrame(rows_overall)
    by_trait_df = pd.DataFrame(rows_trait).sort_values(["trait"])
    return overall_df, by_trait_df, annotated

# ============================ plotting helpers ============================

def scatter_basic(x, y, ax, title, xlabel, ylabel,
                  color="black", alpha=0.6, s=10,
                  regline=False, show_lm=False, show_r=False, show_p=False,
                  reg_filter=None,
                  stats_outside=False, stats_position="top", stats_margin=0.18,
                  stats_panel=False, stats_panel_ax=None):
    m_plot = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[m_plot], y[m_plot], s=s, c=color, alpha=alpha, edgecolors="none")

    mask = m_plot.copy()
    if reg_filter is not None:
        mask &= (x > reg_filter)
    x_reg = x[mask]
    y_reg = y[mask]

    if regline and x_reg.size >= 2:
        add_reg_line(ax, x_reg, y_reg, color=color, alpha=min(1.0, alpha+0.2))

    if show_lm or show_r or show_p:
        txt = format_stats_text(x_reg, y_reg, show_lm, show_r, show_p)
        if txt:
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
                     reg_filter=None,
                     per_trait_stats=False,
                     stats_traits=None,
                     highlight_traits=None,
                     alpha_default=0.25, alpha_highlight=0.9,
                     trait_colors=None,
                     stats_outside=False, stats_position="top", stats_margin=0.18,
                     stats_panel=False, stats_panel_ax=None,
                     xlabel_override=None,
                     stats_trait_xq=80.0, stats_trait_yq=95.0):
    if highlight_traits is None:
        highlight_traits = []
    if stats_traits is None:
        stats_traits = []
    if trait_colors is None:
        trait_colors = {}

    traits = df[trait_col].astype(str).fillna("NA").unique().tolist()
    cmap = cm.get_cmap("tab20", max(10, len(traits)))

    stats_trait_xq = float(np.clip(stats_trait_xq, 0.0, 100.0))
    stats_trait_yq = float(np.clip(stats_trait_yq, 0.0, 100.0))

    # overall stats
    if (show_lm or show_r or show_p) and not per_trait_stats and not df.empty:
        x_all = df[xcol].to_numpy()
        y_all = df[ycol].to_numpy()
        m_all = np.isfinite(x_all) & np.isfinite(y_all)
        if reg_filter is not None:
            m_all &= (x_all > reg_filter)
        x_reg_all = x_all[m_all]
        y_reg_all = y_all[m_all]
        txt = format_stats_text(x_reg_all, y_reg_all, show_lm, show_r, show_p)
        if txt:
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
        if sub.empty:
            continue

        x = sub[xcol].to_numpy()
        y = sub[ycol].to_numpy()

        col = trait_colors.get(str(t), cmap(i))
        alpha = alpha_highlight if str(t) in highlight_traits else alpha_default

        m_plot = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[m_plot], y[m_plot], s=10, alpha=alpha, color=col,
                   label=str(t), edgecolors="none")

        mask = m_plot.copy()
        if reg_filter is not None:
            mask &= (x > reg_filter)
        x_reg = x[mask]
        y_reg = y[mask]

        if regline and x_reg.size >= 2:
            add_reg_line(ax, x_reg, y_reg, color=col, alpha=min(1.0, alpha + 0.2))

        if per_trait_stats and (show_lm or show_r or show_p):
            txt_short = []
            if x_reg.size >= 2:
                m, b = ols_line(x_reg, y_reg); r, p = pearson_r_p(x_reg, y_reg)
                if show_lm and (m is not None): txt_short.append(f"m={m:.2g}")
                if show_r and (r is not None):  txt_short.append(f"r={r:.2f}")
                if show_p: txt_short.append(f"p={p:.1e}" if p is not None else "p=?")
            label = str(t)
            if txt_short:
                label = f"{t} ({', '.join(txt_short)})"
            ax.scatter([], [], color=col, alpha=alpha, label=label)

        if (str(t) in stats_traits) and (show_lm or show_r or show_p) and x_reg.size >= 2:
            txt = format_stats_text(x_reg, y_reg, show_lm, show_r, show_p)
            if txt:
                try:
                    x_pos = np.nanpercentile(x, stats_trait_xq)
                    y_pos = np.nanpercentile(y, stats_trait_yq)
                except Exception:
                    x_pos = np.nanmean(x)
                    y_pos = np.nanmax(y)
                ax.text(x_pos, y_pos, txt,
                        fontsize=8, color=col,
                        ha="left", va="bottom",
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, alpha=0.9))

    ax.set_title(title)
    ax.set_xlabel(xlabel_override if xlabel_override is not None else xcol)
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

# ============================ plotting wrapper ============================

def make_all_plots(annotated, c_trait, c_beta, c_maf, xcol, x_label, args,
                   label_suffix="", fname_suffix="",
                   reg_filter=None,
                   highlight_traits=None, trait_alpha_default=0.25,
                   trait_alpha_highlight=0.9,
                   trait_colors=None, stats_traits=None,
                   exclude_traits=None):
    def savefig(fig, name, dpi=300):
        final = f"{name}{fname_suffix}" if fname_suffix else name
        fig.tight_layout()
        fig.savefig(Path(args.outdir, f"{final}.png"), dpi=dpi)
        fig.savefig(Path(args.outdir, f"{final}.pdf"))

    base_size = (args.fig_w, args.fig_h)

    # ---- TMRCA vs BETA (all points)
    dfp = annotated.dropna(subset=[xcol, c_beta]).copy()
    if not dfp.empty:
        fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                      strip_pos=args.stats_position,
                                      strip_height_in=args.stats_panel_height)
        scatter_basic(
            dfp[xcol].to_numpy(), dfp[c_beta].to_numpy(), ax,
            title=f"{x_label} vs BETA{label_suffix}",
            xlabel=x_label, ylabel="BETA",
            color=args.allpoints_color, alpha=0.6, regline=args.regline,
            show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
            reg_filter=reg_filter,
            stats_outside=args.stats_outside, stats_position=args.stats_position,
            stats_margin=args.stats_margin,
            stats_panel=args.stats_panel, stats_panel_ax=ax_strip
        )
        savefig(fig, "TMRCA_vs_BETA_all", dpi=args.dpi)

    # ---- TMRCA vs |BETA|
    if not dfp.empty:
        dfp["absBETA"] = dfp[c_beta].abs()
        fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                      strip_pos=args.stats_position,
                                      strip_height_in=args.stats_panel_height)
        scatter_basic(
            dfp[xcol].to_numpy(), dfp["absBETA"].to_numpy(), ax,
            title=f"{x_label} vs |BETA|{label_suffix}",
            xlabel=x_label, ylabel="|BETA|",
            color=args.allpoints_color, alpha=0.6, regline=args.regline,
            show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
            reg_filter=reg_filter,
            stats_outside=args.stats_outside, stats_position=args.stats_position,
            stats_margin=args.stats_margin,
            stats_panel=args.stats_panel, stats_panel_ax=ax_strip
        )
        savefig(fig, "TMRCA_vs_absBETA_all", dpi=args.dpi)

    # ---- TMRCA vs MAF
    dfm = annotated.dropna(subset=[xcol, c_maf]).copy()
    if not dfm.empty:
        fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                      strip_pos=args.stats_position,
                                      strip_height_in=args.stats_panel_height)
        scatter_basic(
            dfm[xcol].to_numpy(), dfm[c_maf].to_numpy(), ax,
            title=f"{x_label} vs MAF{label_suffix}",
            xlabel=x_label, ylabel="MAF",
            color=args.allpoints_color, alpha=0.6, regline=args.regline,
            show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
            reg_filter=reg_filter,
            stats_outside=args.stats_outside, stats_position=args.stats_position,
            stats_margin=args.stats_margin,
            stats_panel=args.stats_panel, stats_panel_ax=ax_strip
        )
        savefig(fig, "TMRCA_vs_MAF_all", dpi=args.dpi)

    exclude_traits = set(exclude_traits or [])

    # ---- BY TRAIT: TMRCA vs BETA
    dfc = annotated.dropna(subset=[xcol, c_beta]).copy()
    if not dfc.empty:
        dfc_bt, counts_bt = filter_by_min_points(dfc, c_trait,
                                                 xcol, c_beta,
                                                 args.min_n_trait)
        if exclude_traits:
            dfc_bt = dfc_bt[~dfc_bt[c_trait].astype(str).isin(exclude_traits)]
        if args.min_n_trait > 0 and fname_suffix == "":
            counts_bt.to_csv(Path(args.outdir,
                                  "counts_TMRCA_vs_BETA_byTrait.tsv"),
                             sep="\t", header=["n"])
        if not dfc_bt.empty:
            fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                          strip_pos=args.stats_position,
                                          strip_height_in=args.stats_panel_height)
            scatter_by_trait(
                dfc_bt, xcol=xcol, ycol=c_beta, trait_col=c_trait, ax=ax,
                title=f"{x_label} vs BETA (by trait; min_n={args.min_n_trait}){label_suffix}",
                regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                reg_filter=reg_filter,
                per_trait_stats=args.per_trait_stats,
                stats_traits=stats_traits,
                highlight_traits=highlight_traits,
                alpha_default=trait_alpha_default,
                alpha_highlight=trait_alpha_highlight,
                trait_colors=trait_colors,
                stats_outside=args.stats_outside, stats_position=args.stats_position,
                stats_margin=args.stats_margin,
                stats_panel=args.stats_panel, stats_panel_ax=ax_strip,
                xlabel_override=x_label,
                stats_trait_xq=args.stats_trait_xq,
                stats_trait_yq=args.stats_trait_yq
            )
            savefig(fig, "TMRCA_vs_BETA_byTrait", dpi=args.dpi)

    # ---- BY TRAIT: TMRCA vs |BETA|
    if not dfc.empty:
        dfc_abs = dfc.copy()
        dfc_abs["absBETA"] = dfc_abs[c_beta].abs()
        dfc_abs_f, counts_abt = filter_by_min_points(dfc_abs, c_trait,
                                                     xcol, "absBETA",
                                                     args.min_n_trait)
        if exclude_traits:
            dfc_abs_f = dfc_abs_f[~dfc_abs_f[c_trait].astype(str).isin(exclude_traits)]
        if args.min_n_trait > 0 and fname_suffix == "":
            counts_abt.to_csv(Path(args.outdir,
                                   "counts_TMRCA_vs_absBETA_byTrait.tsv"),
                              sep="\t", header=["n"])
        if not dfc_abs_f.empty:
            fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                          strip_pos=args.stats_position,
                                          strip_height_in=args.stats_panel_height)
            scatter_by_trait(
                dfc_abs_f, xcol=xcol, ycol="absBETA", trait_col=c_trait, ax=ax,
                title=f"{x_label} vs |BETA| (by trait; min_n={args.min_n_trait}){label_suffix}",
                regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                reg_filter=reg_filter,
                per_trait_stats=args.per_trait_stats,
                stats_traits=stats_traits,
                highlight_traits=highlight_traits,
                alpha_default=trait_alpha_default,
                alpha_highlight=trait_alpha_highlight,
                trait_colors=trait_colors,
                stats_outside=args.stats_outside, stats_position=args.stats_position,
                stats_margin=args.stats_margin,
                stats_panel=args.stats_panel, stats_panel_ax=ax_strip,
                xlabel_override=x_label,
                stats_trait_xq=args.stats_trait_xq,
                stats_trait_yq=args.stats_trait_yq
            )
            savefig(fig, "TMRCA_vs_absBETA_byTrait", dpi=args.dpi)

    # ---- BY TRAIT: TMRCA vs MAF
    dfc_maf = annotated.dropna(subset=[xcol, c_maf]).copy()
    if not dfc_maf.empty:
        dfc_maf_f, counts_maf = filter_by_min_points(dfc_maf, c_trait,
                                                     xcol, c_maf,
                                                     args.min_n_trait)
        if exclude_traits:
            dfc_maf_f = dfc_maf_f[~dfc_maf_f[c_trait].astype(str).isin(exclude_traits)]
        if args.min_n_trait > 0 and fname_suffix == "":
            counts_maf.to_csv(Path(args.outdir,
                                   "counts_TMRCA_vs_MAF_byTrait.tsv"),
                              sep="\t", header=["n"])
        if not dfc_maf_f.empty:
            fig, ax, ax_strip = make_axes(base_size, stats_strip=args.stats_panel,
                                          strip_pos=args.stats_position,
                                          strip_height_in=args.stats_panel_height)
            scatter_by_trait(
                dfc_maf_f, xcol=xcol, ycol=c_maf, trait_col=c_trait, ax=ax,
                title=f"{x_label} vs MAF (by trait; min_n={args.min_n_trait}){label_suffix}",
                regline=args.regline, show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                reg_filter=reg_filter,
                per_trait_stats=args.per_trait_stats,
                stats_traits=stats_traits,
                highlight_traits=highlight_traits,
                alpha_default=trait_alpha_default,
                alpha_highlight=trait_alpha_highlight,
                trait_colors=trait_colors,
                stats_outside=args.stats_outside, stats_position=args.stats_position,
                stats_margin=args.stats_margin,
                stats_panel=args.stats_panel, stats_panel_ax=ax_strip,
                xlabel_override=x_label,
                stats_trait_xq=args.stats_trait_xq,
                stats_trait_yq=args.stats_trait_yq
            )
            savefig(fig, "TMRCA_vs_MAF_byTrait", dpi=args.dpi)

# ============================ main ============================

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Annotate GWAS with TMRCA & summarize/plot "
            "(v13: keeps all advanced options + log10(TMRCA) option)."
        )
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

    # Percentile thresholds (high TMRCA)
    ap.add_argument("--top-percentiles", default="5",
                    help="Comma-separated list of percentages (e.g., '5,10,25') for top-X%% calculations.")

    # Data filtering for separate plot sets
    ap.add_argument("--tmrca-filter", type=float, default=None,
                    help="If set, also generate a second set of plots using only rows with TMRCA > cutoff (raw scale).")

    # Regression-only filter (all points shown, lines/stats only above this TMRCA)
    ap.add_argument("--tmrca-regfilter", type=float, default=None,
                    help="If set, regression lines/stats use only points with TMRCA > this threshold. "
                         "Interpreted in the same scale as x-axis (raw or log10).")

    # Trait-specific appearance
    ap.add_argument("--highlight-traits", default="HT",
                    help="Comma-separated trait names to highlight (higher alpha); e.g. 'HT,C13'.")
    ap.add_argument("--trait-alpha-default", type=float, default=0.25,
                    help="Alpha for non-highlight traits in colored trait plots.")
    ap.add_argument("--trait-alpha-highlight", type=float, default=0.9,
                    help="Alpha for highlighted traits in colored trait plots.")
    ap.add_argument("--trait-colors", default="",
                    help="Comma-separated mapping 'HT:red,C13:green,WWD:blue' for colored trait plots.")
    ap.add_argument("--stats-traits", default="",
                    help="Comma-separated trait names for which per-trait regression stats text is drawn on-plot "
                         "(by-trait plots only).")

    # Exclude traits from colored-by-trait plots entirely
    ap.add_argument("--exclude-traits", default="",
                    help="Comma-separated list of traits to exclude from colored-by-trait plots (no points/lines).")

    # Per-trait stats text position (quantiles of x/y within trait)
    ap.add_argument("--stats-trait-xq", type=float, default=80.0,
                    help="X-axis quantile (0-100) used for positioning per-trait stats text.")
    ap.add_argument("--stats-trait-yq", type=float, default=95.0,
                    help="Y-axis quantile (0-100) used for positioning per-trait stats text.")

    # New: log10(TMRCA) option
    ap.add_argument("--tmrca-log10", action="store_true",
                    help="If set, plot x-axis as log10(TMRCA) and compute regression/stats on log10(TMRCA). "
                         "Rows with TMRCA <= 0 are dropped from plots/regression in this mode.")

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

    for c in [c_pos, c_beta, c_maf]:
        gwas[c] = pd.to_numeric(gwas[c], errors="coerce")
    for c in [t_start, t_end, t_mean, t_lo, t_hi]:
        tmrca[c] = pd.to_numeric(tmrca[c], errors="coerce")
    if t_len:
        tmrca[t_len] = pd.to_numeric(tmrca[t_len], errors="coerce")

    gwas = gwas.dropna(subset=[c_scaf, c_pos])

    keep_cols = [t_scaf, t_start, t_end, t_mean, t_lo, t_hi] + ([t_len] if t_len else [])
    annotated = merge_asof_per_chrom(
        gwas, tmrca[keep_cols].copy(),
        scaf_col=c_scaf, pos_col=c_pos,
        t_scaf_col=t_scaf, t_start_col=t_start, t_end_col=t_end
    )

    rename_map = {t_scaf:"tmrca_scaffold", t_start:"start_tmrca", t_end:"end_tmrca",
                  t_mean:"mean_tmrca", t_lo:"Low_CI", t_hi:"UP_CI"}
    if t_len:
        rename_map[t_len] = "seg_length"
    annotated = annotated.rename(columns=rename_map)

    # ===================== x-axis TMRCA variable: raw or log10 =====================
    if args.tmrca_log10:
        annotated["tmrca_x"] = np.where(
            annotated["mean_tmrca"] > 0,
            np.log10(annotated["mean_tmrca"]),
            np.nan
        )
        x_label = "TMRCA (log10 scale)"
    else:
        annotated["tmrca_x"] = annotated["mean_tmrca"]
        x_label = "TMRCA"

    xcol = "tmrca_x"

    # ===================== high-TMRCA: percent of SNPs in top-X% =====================
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

    # ===================== low-TMRCA: below global mean TMRCA =====================
    global_mean_tmrca = float(np.nanmean(tmrca[t_mean].to_numpy()))
    overall_low, by_trait_low, annotated = summarize_snps_below_mean(
        annotated, trait_col=c_trait, mean_col="mean_tmrca",
        global_mean=global_mean_tmrca
    )

    low_overall_path = Path(args.outdir, "snps_below_global_mean_tmrca_overall.tsv")
    low_bytrait_path = Path(args.outdir, "snps_below_global_mean_tmrca_by_trait.tsv")
    below_only_path  = Path(args.outdir, "gwas_with_tmrca_below_global_mean.tsv")

    overall_low.to_csv(low_overall_path, sep="\t", index=False)
    by_trait_low.to_csv(low_bytrait_path, sep="\t", index=False)
    annotated[annotated["below_global_mean_tmrca"]].to_csv(below_only_path, sep="\t", index=False)

    # ===================== main outputs on full data =====================
    out_annot = Path(args.outdir) / "gwas_with_tmrca.tsv"
    annotated.to_csv(out_annot, sep="\t", index=False)

    trait_summ = trait_summary(annotated, trait_col=c_trait, tmrca_col="mean_tmrca")
    out_summ = Path(args.outdir) / "tmrca_trait_summary.tsv"
    trait_summ.to_csv(out_summ, sep="\t", index=False)

    # trait option parsing
    highlight_traits = parse_trait_list(args.highlight_traits)
    stats_traits     = parse_trait_list(args.stats_traits)
    trait_colors     = parse_trait_colors(args.trait_colors)
    exclude_traits   = parse_trait_list(args.exclude_traits)

    # ===================== plots: full data =====================
    make_all_plots(
        annotated, c_trait, c_beta, c_maf, xcol=xcol, x_label=x_label, args=args,
        label_suffix="", fname_suffix="",
        reg_filter=args.tmrca_regfilter,
        highlight_traits=highlight_traits,
        trait_alpha_default=args.trait_alpha_default,
        trait_alpha_highlight=args.trait_alpha_highlight,
        trait_colors=trait_colors,
        stats_traits=stats_traits,
        exclude_traits=exclude_traits
    )

    # ===================== plots: filtered by raw TMRCA cutoff (optional) =====================
    if args.tmrca_filter is not None:
        filt_val = args.tmrca_filter
        ann_filt = annotated[annotated["mean_tmrca"] > filt_val].copy()
        if not ann_filt.empty:
            suffix_label = f" (TMRCA>{filt_val})"
            suffix_fname = f"_TMRCAgt{str(filt_val).replace('.','p')}"
            make_all_plots(
                ann_filt, c_trait, c_beta, c_maf, xcol=xcol, x_label=x_label, args=args,
                label_suffix=suffix_label, fname_suffix=suffix_fname,
                reg_filter=args.tmrca_regfilter,
                highlight_traits=highlight_traits,
                trait_alpha_default=args.trait_alpha_default,
                trait_alpha_highlight=args.trait_alpha_highlight,
                trait_colors=trait_colors,
                stats_traits=stats_traits,
                exclude_traits=exclude_traits
            )
        else:
            print(f"[WARN] No rows with TMRCA > {filt_val}; filtered plots skipped.")

    print("[OK] Wrote:")
    print(f"  Annotated GWAS:                 {out_annot}")
    print(f"  Top-TMRCA overall:              {overall_path}")
    print(f"  Top-TMRCA by trait:             {bytrait_path}")
    print(f"  Low-TMRCA overall (< mean):     {low_overall_path}")
    print(f"  Low-TMRCA by trait (< mean):    {low_bytrait_path}")
    print(f"  SNPs below mean TMRCA (subset): {below_only_path}")
    print(f"  Trait TMRCA summary:            {out_summ}")
    print(f"  figures & counts under {args.outdir}/")

if __name__ == "__main__":
    main()
