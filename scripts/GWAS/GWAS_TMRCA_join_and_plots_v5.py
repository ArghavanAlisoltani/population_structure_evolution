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

python GWAS_TMRCA_join_and_plots_v5.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025/To_merge_with_TMRCA_MTAG_SIgnificant_sumstat.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_MTAG_trait_tmrca_results_v5" \
  --allpoints-color "#000000" \
  --regline --show-lm --show-r --show-p \
  --stats-outside --stats-position right \
  --min-n-trait 5 

cd /Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files

python GWAS_TMRCA_join_and_plots_v5.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files/combined_pval_maf_hits/combined_hits_full_merge_with_tmrca.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_single_trait_tmrca_results_v5" \
  --allpoints-color "#000000" \
  --regline --show-lm --show-r --show-p \
  --stats-outside --stats-position top \
  --min-n-trait 5 
  
'''
# Optional: Pearson r and p via SciPy if available
try:
    from scipy import stats as sstats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

# ------------------------- utils -------------------------

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
    For each scaffold:
      - asof-merge SNPs to nearest segment start at/left of position
      - keep rows only if position ∈ [start, end]
      - CRITICAL: if not inside, set ALL TMRCA columns to NaN (not just start/end)
    Coerces join keys to float64 and sorts per scaffold.
    """
    out = []
    for scaf, gsub in gwas.groupby(scaf_col, sort=False):
        tsub = tmrca[tmrca[t_scaf_col] == scaf]
        gsub = gsub.copy()
        if tsub.empty:
            out.append(gsub)
            continue

        # match key dtypes and sort
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

        # position truly inside interval?
        inside = merged[pos_col].between(merged[t_start_col], merged[t_end_col], inclusive="both")

        # Clear *all* TMRCA fields for non-inside matches
        tmrca_cols = [c for c in tsub.columns if c != t_scaf_col]
        merged.loc[~inside, tmrca_cols] = np.nan

        out.append(merged)

    return pd.concat(out, ignore_index=True) if out else gwas.copy()

def trait_summary(df, trait_col, tmrca_col):
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
    # Fallback: r only
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

def place_outside_text(fig, ax, text, position="top", color="black"):
    if not text:
        return
    box = ax.get_position()  # figure-relative bbox
    if position == "top":
        # add some headroom
        plt.subplots_adjust(top=0.86)
        x = (box.x0 + box.x1) / 2.0
        y = box.y1 + 0.03
        ha, va = "center", "bottom"
    elif position == "bottom":
        plt.subplots_adjust(bottom=0.20)
        x = (box.x0 + box.x1) / 2.0
        y = box.y0 - 0.12
        ha, va = "center", "top"
    else:  # right
        plt.subplots_adjust(right=0.80)
        x = box.x1 + 0.02
        y = (box.y0 + box.y1) / 2.0
        ha, va = "left", "center"

    fig.text(x, y, text, color=color, fontsize=9, ha=ha, va=va,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.75))

def scatter_basic(x, y, ax, title, xlabel, ylabel,
                  color="black", alpha=0.6, s=10,
                  regline=False, show_lm=False, show_r=False, show_p=False,
                  stats_outside=False, stats_position="top"):
    ax.scatter(x, y, s=s, c=color, alpha=alpha, edgecolors="none")
    if regline:
        add_reg_line(ax, x, y, color=color)

    if show_lm or show_r or show_p:
        txt = format_stats_text(x, y, show_lm, show_r, show_p)
        if stats_outside:
            place_outside_text(ax.figure, ax, txt, position=stats_position, color=color)
        else:
            # in-plot fallback (not requested, but keep support)
            ax.text(0.02, 0.98, txt, transform=ax.transAxes,
                    ha="left", va="top", fontsize=9, color=color,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.75))

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.15)

def scatter_by_trait(df, xcol, ycol, trait_col, ax, title,
                     regline=False, show_lm=False, show_r=False, show_p=False,
                     per_trait_stats=False, stats_outside=False, stats_position="top"):
    traits = df[trait_col].astype(str).fillna("NA").unique().tolist()
    cmap = cm.get_cmap("tab20", max(10, len(traits)))

    # Overall stats across filtered subset
    if (show_lm or show_r or show_p) and not per_trait_stats and not df.empty:
        txt = format_stats_text(df[xcol].to_numpy(), df[ycol].to_numpy(), show_lm, show_r, show_p)
        if stats_outside:
            place_outside_text(ax.figure, ax, txt, position=stats_position, color="black")
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

        # Per-trait stats go in legend label (keeps plot uncluttered)
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
    """Keep only traits with >= min_n finite x & y points (only for colored plots)."""
    if min_n <= 0 or df.empty:
        return df, pd.Series(dtype=int)
    dvalid = df[[trait_col, xcol, ycol]].dropna()
    counts = dvalid.groupby(trait_col, observed=True).size().sort_values(ascending=False)
    keep = counts[counts >= min_n].index
    return df[df[trait_col].isin(keep)], counts

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Annotate GWAS with TMRCA & summarize/plot (v5: NaN-on-no-overlap, all-points color, outside stats, trait min-N)"
    )
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--tmrca", required=True)
    ap.add_argument("--outdir", default="gwas_tmrca_out")

    # column names
    ap.add_argument("--trait-col", default="Significant_Trait")
    ap.add_argument("--scaf-col", default="Scaffold")
    ap.add_argument("--pos-col",  default="position", type=str)
    ap.add_argument("--beta-col", default="BETA")
    ap.add_argument("--maf-col",  default="MAF")

    # plotting toggles + stats overlays
    ap.add_argument("--regline", action="store_true", help="Draw OLS regression line(s)")
    ap.add_argument("--show-lm", action="store_true", help="Show line formula y = mx + b")
    ap.add_argument("--show-r", action="store_true", help="Show Pearson r")
    ap.add_argument("--show-p", action="store_true", help="Show two-tailed p-value for Pearson r")
    ap.add_argument("--per-trait-stats", action="store_true",
                    help="For colored plots, append stats per trait inside the legend")

    # NEW: print stats outside the axes so dots aren’t hidden
    ap.add_argument("--stats-outside", action="store_true",
                    help="Place stats outside the plot area instead of inside")
    ap.add_argument("--stats-position", default="top", choices=["top","bottom","right"],
                    help="Where to place outside stats text")

    # NEW: color for black 'all-points' plots only
    ap.add_argument("--allpoints-color", default="#000000",
                    help="Color for all-points (non-trait-colored) figures")

    # NEW: minimum points per trait for colored-by-trait plots only
    ap.add_argument("--min-n-trait", type=int, default=0,
                    help="Exclude traits with < N points in colored-by-trait plots")

    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # ---- load data
    gwas = pd.read_csv(args.gwas, sep="\t", dtype={})
    tmrca = pd.read_csv(args.tmrca, sep="\t", dtype={})

    # normalize column names
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

    # numeric coercions (float to tolerate NaN)
    for c in [c_pos, c_beta, c_maf]: gwas[c] = pd.to_numeric(gwas[c], errors="coerce")
    for c in [t_start, t_end, t_mean, t_lo, t_hi]: tmrca[c] = pd.to_numeric(tmrca[c], errors="coerce")
    if t_len: tmrca[t_len] = pd.to_numeric(tmrca[t_len], errors="coerce")

    gwas = gwas.dropna(subset=[c_scaf, c_pos])

    # interval join
    keep_cols = [t_scaf, t_start, t_end, t_mean, t_lo, t_hi] + ([t_len] if t_len else [])
    annotated = merge_asof_per_chrom(
        gwas, tmrca[keep_cols].copy(),
        scaf_col=c_scaf, pos_col=c_pos,
        t_scaf_col=t_scaf, t_start_col=t_start, t_end_col=t_end
    )

    # rename TMRCA columns in output
    rename_map = {t_scaf:"tmrca_scaffold", t_start:"start_tmrca", t_end:"end_tmrca",
                  t_mean:"mean_tmrca", t_lo:"Low_CI", t_hi:"UP_CI"}
    if t_len: rename_map[t_len] = "seg_length"
    annotated = annotated.rename(columns=rename_map)

    # write full annotated table
    out_annot = Path(args.outdir) / "gwas_with_tmrca.tsv"
    annotated.to_csv(out_annot, sep="\t", index=False)

    # trait summary
    trait_summ = trait_summary(annotated, trait_col=c_trait, tmrca_col="mean_tmrca")
    out_summ = Path(args.outdir) / "tmrca_trait_summary.tsv"
    trait_summ.to_csv(out_summ, sep="\t", index=False)

    def savefig(fig, name):
        fig.tight_layout()
        fig.savefig(Path(args.outdir, f"{name}.png"), dpi=args.dpi)
        fig.savefig(Path(args.outdir, f"{name}.pdf"))

    # ----------------------- PLOTS (ALL-POINTS; one color) -----------------------
    dfp = annotated.dropna(subset=["mean_tmrca", c_beta]).copy()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfp[c_beta].to_numpy(), dfp["mean_tmrca"].to_numpy(), ax,
                  title="BETA vs mean TMRCA", xlabel="BETA", ylabel="mean TMRCA",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "BETA_vs_TMRCA_all")

    dfp["absBETA"] = dfp[c_beta].abs()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfp["absBETA"].to_numpy(), dfp["mean_tmrca"].to_numpy(), ax,
                  title="|BETA| vs mean TMRCA", xlabel="|BETA|", ylabel="mean TMRCA",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "absBETA_vs_TMRCA_all")

    dfm = annotated.dropna(subset=["mean_tmrca", c_maf]).copy()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfm[c_maf].to_numpy(), dfm["mean_tmrca"].to_numpy(), ax,
                  title="MAF vs mean TMRCA", xlabel="MAF", ylabel="mean TMRCA",
                  color=args.allpoints_color, regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "MAF_vs_TMRCA_all")

    # ----------------------- PLOTS (BY-TRAIT with min-N filter) -----------------------
    # BETA vs mean TMRCA (colored)
    dfc = annotated.dropna(subset=["mean_tmrca", c_beta]).copy()
    dfc_bt, counts_bt = filter_by_min_points(dfc, c_trait, c_beta, "mean_tmrca", args.min_n_trait)
    if args.min_n_trait > 0:
        counts_bt.to_csv(Path(args.outdir, "counts_BETA_vs_TMRCA_byTrait.tsv"), sep="\t", header=["n"])
    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfc_bt, xcol=c_beta, ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title=f"BETA vs mean TMRCA (by trait; min_n={args.min_n_trait})",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "BETA_vs_TMRCA_byTrait")

    # |BETA| vs mean TMRCA (colored)
    dfc_abs = dfc.copy()
    dfc_abs["absBETA"] = dfc_abs[c_beta].abs()
    dfc_abs_f, counts_abt = filter_by_min_points(dfc_abs, c_trait, "absBETA", "mean_tmrca", args.min_n_trait)
    if args.min_n_trait > 0:
        counts_abt.to_csv(Path(args.outdir, "counts_absBETA_vs_TMRCA_byTrait.tsv"), sep="\t", header=["n"])
    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfc_abs_f, xcol="absBETA", ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title=f"|BETA| vs mean TMRCA (by trait; min_n={args.min_n_trait})",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "absBETA_vs_TMRCA_byTrait")

    # MAF vs mean TMRCA (colored)
    dfc_maf = annotated.dropna(subset=["mean_tmrca", c_maf]).copy()
    dfc_maf_f, counts_maf = filter_by_min_points(dfc_maf, c_trait, c_maf, "mean_tmrca", args.min_n_trait)
    if args.min_n_trait > 0:
        counts_maf.to_csv(Path(args.outdir, "counts_MAF_vs_TMRCA_byTrait.tsv"), sep="\t", header=["n"])
    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfc_maf_f, xcol=c_maf, ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title=f"MAF vs mean TMRCA (by trait; min_n={args.min_n_trait})",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats,
                     stats_outside=args.stats_outside, stats_position=args.stats_position)
    savefig(fig, "MAF_vs_TMRCA_byTrait")

    print(f"[OK] Wrote:\n  {out_annot}\n  {out_summ}\n  figures & counts under {args.outdir}/")

if __name__ == "__main__":
    main()
