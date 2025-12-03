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

python GWAS_TMRCA_join_and_plots_v3.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/MTAG_NOV_2025/To_merge_with_TMRCA_MTAG_SIgnificant_sumstat.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_MTAG_trait_tmrca_results_v3" \
  --regline --show-lm --show-r --show-p

cd /Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files

python GWAS_TMRCA_join_and_plots_v3.py \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/sumstat_files/combined_pval_maf_hits/combined_hits_full_merge_with_tmrca.txt" \
  --tmrca "/Users/aria/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_as_scaffold1.tsv" \
  --outdir "gwas_single_trait_tmrca_results_v3" \
  --regline --show-lm --show-r --show-p
'''

# Optional: Pearson r, p-value via SciPy if available
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
    For each scaffold, asof-merge positions onto segment starts, then
    mask rows where position > end_tmrca. Coerces join keys to float64.
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
        merged.loc[~inside, [t_start_col, t_end_col]] = np.nan
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
    # Fallback: r only (no SciPy)
    r = np.corrcoef(x[msk], y[msk])[0,1]
    return float(r), None

def add_reg_line(ax, x, y, color="black", lw=1.2, alpha=0.8):
    m, b = ols_line(x, y)
    if m is None:
        return None, None
    xs = np.linspace(np.nanmin(x[np.isfinite(x)]), np.nanmax(x[np.isfinite(x)]), 200)
    ax.plot(xs, m*xs + b, color=color, linewidth=lw, alpha=alpha)
    return m, b

def annotate_stats(ax, x, y, show_lm, show_r, show_p, color="black", loc="ul"):
    """
    loc: 'ul'|'ur'|'ll'|'lr' (upper/lower left/right)
    """
    lines = []
    m, b = ols_line(x, y)
    r, p = pearson_r_p(x, y)

    if show_lm and (m is not None):
        lines.append(f"y = {m:.3g}x + {b:.3g}")
    if show_r and (r is not None):
        lines.append(f"r = {r:.3f}")
    if show_p:
        if p is not None:
            lines.append(f"p = {p:.2e}")
        else:
            lines.append("p = (requires SciPy)")

    if not lines: 
        return

    txt = "\n".join(lines)
    ha = "left" if loc.endswith("l") else "right"
    va = "top"  if loc.startswith("u") else "bottom"
    x0 = 0.02 if ha=="left" else 0.98
    y0 = 0.98 if va=="top"  else 0.02
    ax.text(x0, y0, txt, transform=ax.transAxes, ha=ha, va=va,
            fontsize=9, color=color,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.75))

def scatter_basic(x, y, ax, title, xlabel, ylabel,
                  color="black", alpha=0.6, s=10,
                  regline=False, show_lm=False, show_r=False, show_p=False,
                  annot_loc="ul"):
    ax.scatter(x, y, s=s, c=color, alpha=alpha, edgecolors="none")
    if regline:
        add_reg_line(ax, x, y, color=color)
    if any([show_lm, show_r, show_p]):
        annotate_stats(ax, x, y, show_lm, show_r, show_p, color=color, loc=annot_loc)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.15)

def scatter_by_trait(df, xcol, ycol, trait_col, ax, title,
                     regline=False, show_lm=False, show_r=False, show_p=False,
                     per_trait_stats=False, annot_loc="ul"):
    traits = df[trait_col].astype(str).fillna("NA").unique().tolist()
    cmap = cm.get_cmap("tab20", max(10, len(traits)))

    # Overall annotation (across all points), if requested and *not* per-trait
    if (show_lm or show_r or show_p) and not per_trait_stats:
        annotate_stats(ax, df[xcol].to_numpy(), df[ycol].to_numpy(),
                       show_lm, show_r, show_p, color="black", loc=annot_loc)

    for i, t in enumerate(traits):
        sub = df[df[trait_col].astype(str).fillna("NA") == t]
        col = cmap(i)
        ax.scatter(sub[xcol], sub[ycol], s=10, alpha=0.7, color=col, label=t, edgecolors="none")
        if regline:
            add_reg_line(ax, sub[xcol].to_numpy(), sub[ycol].to_numpy(), color=col)
        if per_trait_stats and (show_lm or show_r or show_p):
            # put stats inside legend label
            xv = sub[xcol].to_numpy(); yv = sub[ycol].to_numpy()
            m, b = ols_line(xv, yv)
            r, p = pearson_r_p(xv, yv)
            bits = []
            if show_lm and (m is not None): bits.append(f"m={m:.2g}")
            if show_r and (r is not None):  bits.append(f"r={r:.2f}")
            if show_p:
                bits.append(f"p={p:.1e}" if p is not None else "p=?")
            if bits:
                ax.scatter([], [], color=col, label=f"{t} ({', '.join(bits)})")

    ax.set_title(title)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.legend(markerscale=1.5, fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.15)

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser(description="Annotate GWAS with TMRCA & summarize/plot (v3 with stats overlays)")
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--tmrca", required=True)
    ap.add_argument("--outdir", default="gwas_tmrca_out")

    # column names
    ap.add_argument("--trait-col", default="Significant_Trait")
    ap.add_argument("--scaf-col", default="Scaffold")
    ap.add_argument("--pos-col",  default="position", type=str)
    ap.add_argument("--beta-col", default="BETA")
    ap.add_argument("--maf-col",  default="MAF")

    # plotting toggles
    ap.add_argument("--regline", action="store_true", help="Draw OLS regression line(s)")
    ap.add_argument("--show-lm", action="store_true", help="Show line formula y = mx + b")
    ap.add_argument("--show-r", action="store_true", help="Show Pearson r")
    ap.add_argument("--show-p", action="store_true", help="Show two-tailed p-value for Pearson r")
    ap.add_argument("--per-trait-stats", action="store_true",
                    help="For by-trait plots, append stats per trait inside the legend")
    ap.add_argument("--annot-loc", default="ul", choices=["ul","ur","ll","lr"],
                    help="Where to place the stats box for 'all-points' plots")

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

    # numeric coercions (keep as float to tolerate NaN)
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

    # helper to save figs
    def savefig(fig, name):
        fig.tight_layout()
        fig.savefig(Path(args.outdir, f"{name}.png"), dpi=args.dpi)
        fig.savefig(Path(args.outdir, f"{name}.pdf"))

    # ----------------------- PLOTS -----------------------
    # BETA vs mean TMRCA (all points, single color)
    dfp = annotated.dropna(subset=["mean_tmrca", c_beta]).copy()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfp[c_beta].to_numpy(), dfp["mean_tmrca"].to_numpy(), ax,
                  title="BETA vs mean TMRCA", xlabel="BETA", ylabel="mean TMRCA",
                  color="black", regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  annot_loc=args.annot_loc)
    savefig(fig, "BETA_vs_TMRCA_all")

    # BETA vs mean TMRCA colored by trait
    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfp, xcol=c_beta, ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title="BETA vs mean TMRCA (by trait)",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats, annot_loc=args.annot_loc)
    savefig(fig, "BETA_vs_TMRCA_byTrait")

    # |BETA| vs mean TMRCA
    dfp["absBETA"] = dfp[c_beta].abs()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfp["absBETA"].to_numpy(), dfp["mean_tmrca"].to_numpy(), ax,
                  title="|BETA| vs mean TMRCA", xlabel="|BETA|", ylabel="mean TMRCA",
                  color="black", regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  annot_loc=args.annot_loc)
    savefig(fig, "absBETA_vs_TMRCA_all")

    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfp, xcol="absBETA", ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title="|BETA| vs mean TMRCA (by trait)",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats, annot_loc=args.annot_loc)
    savefig(fig, "absBETA_vs_TMRCA_byTrait")

    # MAF vs mean TMRCA
    dfm = annotated.dropna(subset=["mean_tmrca", c_maf]).copy()
    fig, ax = plt.subplots(figsize=(7,5))
    scatter_basic(dfm[c_maf].to_numpy(), dfm["mean_tmrca"].to_numpy(), ax,
                  title="MAF vs mean TMRCA", xlabel="MAF", ylabel="mean TMRCA",
                  color="black", regline=args.regline,
                  show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                  annot_loc=args.annot_loc)
    savefig(fig, "MAF_vs_TMRCA_all")

    fig, ax = plt.subplots(figsize=(8,5.5))
    scatter_by_trait(dfm, xcol=c_maf, ycol="mean_tmrca",
                     trait_col=c_trait, ax=ax,
                     title="MAF vs mean TMRCA (by trait)",
                     regline=args.regline,
                     show_lm=args.show_lm, show_r=args.show_r, show_p=args.show_p,
                     per_trait_stats=args.per_trait_stats, annot_loc=args.annot_loc)
    savefig(fig, "MAF_vs_TMRCA_byTrait")

    print(f"[OK] Wrote\n  {out_annot}\n  {out_summ}\n  figures under {args.outdir}/")

if __name__ == "__main__":
    main()

