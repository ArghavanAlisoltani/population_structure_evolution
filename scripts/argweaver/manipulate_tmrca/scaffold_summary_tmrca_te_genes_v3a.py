#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Run example

'''
import argparse, os, re
from typing import List
import pandas as pd
import numpy as np

# ---------- helpers ----------

def load_table(path, default_sep="\t"):
    path = os.path.expanduser(path)
    try:
        return pd.read_csv(path, sep=default_sep, low_memory=False)
    except Exception:
        return pd.read_csv(path, delim_whitespace=True, low_memory=False)

def require(df: pd.DataFrame, cols: List[str], label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {label}: {missing}\nAvailable: {list(df.columns)}")

def sanitize(col: str) -> str:
    col = re.sub(r"\s+", "_", str(col))
    col = re.sub(r"[^0-9A-Za-z_]+", "_", col)
    return col

def intervals_union_bp(df: pd.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pd.Series:
    """
    Per-scaffold union coverage length (bp), using 1-based inclusive lengths: end - start + 1.
    Overlaps/adjacencies merged and counted once.
    Returns a Series indexed by scaffold with 'union_bp' values.
    """
    out = {}
    if df.empty:
        return pd.Series(out, name="union_bp")
    for scf, sub in df.groupby(chrom_col, sort=False):
        if sub.empty:
            out[str(scf)] = 0
            continue
        s = pd.to_numeric(sub[start_col], errors="coerce").dropna().astype(np.int64).to_numpy()
        e = pd.to_numeric(sub[end_col],   errors="coerce").dropna().astype(np.int64).to_numpy()
        if len(s) == 0 or len(e) == 0:
            out[str(scf)] = 0
            continue
        arr = np.stack([s, e], axis=1)
        arr = arr[np.argsort(arr[:,0])]
        total = 0
        cur_s, cur_e = None, None
        for a, b in arr:
            if a > b: a, b = b, a
            if cur_s is None:
                cur_s, cur_e = int(a), int(b)
            else:
                if a <= cur_e + 1:  # overlap or touching
                    if b > cur_e:
                        cur_e = int(b)
                else:
                    total += (cur_e - cur_s + 1)
                    cur_s, cur_e = int(a), int(b)
        if cur_s is not None:
            total += (cur_e - cur_s + 1)
        out[str(scf)] = int(total)
    return pd.Series(out, name="union_bp")

def q(x, p):
    try:
        return float(np.nanquantile(x, p))
    except Exception:
        return np.nan

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description=("Per-scaffold summary across TMRCA, genes, and TEs, including coverage, "
                     "TE-class breakdown, and top-percentile TMRCA metrics. "
                     "Also outputs an extra combined row for scaffold_1a + scaffold_1b as scaffold_1.")
    )
    # Inputs
    ap.add_argument("--tmrca",
        default="~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025/all_tmrca_corrected_position.tsv",
        help="TMRCA segments TSV (CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    ap.add_argument("--mrna",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt",
        help="mRNA annotation TSV (new_seqid, new_start, new_end, mrna_id)")
    ap.add_argument("--te",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv",
        help="TE annotation TSV (seqid, sequence_ontology, start, end)")

    # Output
    ap.add_argument("--out",
        default="scaffold_summary_v3.tsv",
        help="Output TSV path")

    # Column overrides if needed
    ap.add_argument("--tmrca-chrom", default="CHROM")
    ap.add_argument("--tmrca-start", default="start_tmrca")
    ap.add_argument("--tmrca-end",   default="end_tmrca")
    ap.add_argument("--tmrca-mean",  default="mean_tmrca")

    ap.add_argument("--mrna-chrom",  default="new_seqid")
    ap.add_argument("--mrna-start",  default="new_start")
    ap.add_argument("--mrna-end",    default="new_end")
    ap.add_argument("--mrna-id",     default="mrna_id")

    ap.add_argument("--te-chrom",    default="seqid")
    ap.add_argument("--te-start",    default="start")
    ap.add_argument("--te-end",      default="end")
    ap.add_argument("--te-class",    default="sequence_ontology")

    # Top-percentile threshold for "high TMRCA" (global)
    ap.add_argument("--top-percentile", type=float, default=25.0,
                    help="Top percentile for high TMRCA threshold (e.g., 25 => top quartile).")

    args = ap.parse_args()

    # ---------- Load data ----------
    tdf = load_table(args.tmrca)
    require(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], "TMRCA file")
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    tdf[args.tmrca_start] = pd.to_numeric(tdf[args.tmrca_start], errors="coerce")
    tdf[args.tmrca_end]   = pd.to_numeric(tdf[args.tmrca_end],   errors="coerce")
    tdf[args.tmrca_mean]  = pd.to_numeric(tdf[args.tmrca_mean],  errors="coerce")
    tdf["seg_len"]        = (tdf[args.tmrca_end] - tdf[args.tmrca_start] + 1).clip(lower=0)

    gdf = load_table(args.mrna)
    require(gdf, [args.mrna_chrom, args.mrna_start, args.mrna_end], "mRNA file")
    gdf[args.mrna_chrom] = gdf[args.mrna_chrom].astype(str)
    gdf[args.mrna_start] = pd.to_numeric(gdf[args.mrna_start], errors="coerce")
    gdf[args.mrna_end]   = pd.to_numeric(gdf[args.mrna_end],   errors="coerce")
    if args.mrna_id not in gdf.columns:
        gdf[args.mrna_id] = np.nan
    else:
        gdf[args.mrna_id] = gdf[args.mrna_id].astype(str)

    tte = load_table(args.te)
    require(tte, [args.te_chrom, args.te_start, args.te_end], "TE file")
    if args.te_class not in tte.columns:
        tte[args.te_class] = "TE"
    tte[args.te_chrom] = tte[args.te_chrom].astype(str)
    tte[args.te_start] = pd.to_numeric(tte[args.te_start], errors="coerce")
    tte[args.te_end]   = pd.to_numeric(tte[args.te_end],   errors="coerce")
    tte[args.te_class] = tte[args.te_class].astype(str)

    # ---------- Universe of scaffolds (all inputs) ----------
    scaff_all = pd.Index(
        sorted(set(tdf[args.tmrca_chrom].unique())
             | set(gdf[args.mrna_chrom].unique())
             | set(tte[args.te_chrom].unique())),
        dtype="object"
    )
    base = pd.DataFrame({"scaffold": scaff_all})

    # ---------- TMRCA per-scaffold stats ----------
    tmrca_stats = (
        tdf.groupby(args.tmrca_chrom)[args.tmrca_mean]
           .agg(n_tmrca_segments="size",
                tmrca_min="min",
                tmrca_q1=lambda s: q(s, 0.25),
                tmrca_median="median",
                tmrca_q3=lambda s: q(s, 0.75),
                tmrca_max="max",
                tmrca_mean="mean")
           .reset_index()
           .rename(columns={args.tmrca_chrom: "scaffold"})
    )

    # TMRCA SD (robust creation)
    tmrca_sd_df = (
        tdf.groupby(args.tmrca_chrom)[args.tmrca_mean]
           .agg(lambda s: float(np.nanstd(s, ddof=1)) if s.notna().sum() >= 2 else np.nan)
           .reset_index()
           .rename(columns={args.tmrca_chrom: "scaffold", args.tmrca_mean: "tmrca_sd"})
    )

    # TMRCA coverage (union bp) & seg length stats
    tmrca_cov = intervals_union_bp(
        tdf.rename(columns={args.tmrca_chrom:"scaffold",
                            args.tmrca_start:"start",
                            args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    seg_len_stats = (tdf.groupby(args.tmrca_chrom)["seg_len"]
                        .agg(median_seg_len="median", mean_seg_len="mean")
                        .reset_index().rename(columns={args.tmrca_chrom:"scaffold"}))

    # Global top-percentile threshold
    top_p = float(args.top_percentile)
    top_thr = tdf[args.tmrca_mean].quantile(1.0 - top_p/100.0)
    tdf["is_top"] = tdf[args.tmrca_mean] >= top_thr

    top_count = (tdf.groupby([args.tmrca_chrom, "is_top"])
                    .size().unstack(fill_value=0).reset_index()
                    .rename(columns={args.tmrca_chrom:"scaffold"}))
    if True not in top_count.columns:  top_count[True]  = 0
    if False not in top_count.columns: top_count[False] = 0
    top_count = top_count.rename(columns={True:"n_top_tmrca_segments", False:"n_non_top_tmrca_segments"})

    top_cov = intervals_union_bp(
        tdf.loc[tdf["is_top"]].rename(columns={args.tmrca_chrom:"scaffold",
                                               args.tmrca_start:"start",
                                               args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("top_tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- Genes per scaffold ----------
    gene_counts = (gdf.groupby(args.mrna_chrom)[args.mrna_id]
                      .nunique(dropna=True)
                      .reset_index()
                      .rename(columns={args.mrna_chrom:"scaffold", args.mrna_id:"n_genes"}))

    gene_cov = intervals_union_bp(
        gdf.rename(columns={args.mrna_chrom:"scaffold",
                            args.mrna_start:"start",
                            args.mrna_end:"end"}),
        "scaffold","start","end"
    ).rename("gene_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- TEs per scaffold ----------
    te_tot = (tte.groupby(args.te_chrom).size()
                 .reset_index(name="n_TEs_total")
                 .rename(columns={args.te_chrom:"scaffold"}))

    n_te_classes = (tte.groupby(args.te_chrom)[args.te_class]
                       .nunique()
                       .reset_index(name="n_TE_classes")
                       .rename(columns={args.te_chrom:"scaffold"}))

    te_class_counts = (tte.groupby([args.te_chrom, args.te_class])
                          .size()
                          .reset_index(name="count")
                          .rename(columns={args.te_chrom:"scaffold", args.te_class:"TE_class"}))
    if len(te_class_counts) > 0:
        te_class_counts["TE_col"] = te_class_counts["TE_class"].map(lambda x: "TEclass_" + sanitize(x))
        te_wide = (te_class_counts.pivot_table(index="scaffold", columns="TE_col",
                                               values="count", fill_value=0, aggfunc="sum")
                                 .reset_index())
    else:
        te_wide = pd.DataFrame(columns=["scaffold"])

    te_cov = intervals_union_bp(
        tte.rename(columns={args.te_chrom:"scaffold",
                            args.te_start:"start",
                            args.te_end:"end"}),
        "scaffold","start","end"
    ).rename("te_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- Merge all pieces on ALL scaffolds ----------
    out = (base
           .merge(tmrca_stats, on="scaffold", how="left")
           .merge(tmrca_sd_df, on="scaffold", how="left")
           .merge(tmrca_cov,   on="scaffold", how="left")
           .merge(seg_len_stats, on="scaffold", how="left")
           .merge(top_count,   on="scaffold", how="left")
           .merge(top_cov,     on="scaffold", how="left")
           .merge(gene_counts, on="scaffold", how="left")
           .merge(gene_cov,    on="scaffold", how="left")
           .merge(te_tot,      on="scaffold", how="left")
           .merge(n_te_classes,on="scaffold", how="left")
           .merge(te_cov,      on="scaffold", how="left")
           .merge(te_wide,     on="scaffold", how="left")
    )

    # Ensure tmrca_sd column exists
    if "tmrca_sd" not in out.columns:
        out["tmrca_sd"] = np.nan

    # Fill count/coverage NaNs with 0; keep TMRCA stats as floats
    zero_int_cols = ["n_tmrca_segments","tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
                     "top_tmrca_coverage_bp","n_genes","gene_coverage_bp","n_TEs_total","n_TE_classes","te_coverage_bp"]
    for c in zero_int_cols:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)

    # Derived ratios (robust to missing columns)
    gene_cov_series = out["gene_coverage_bp"] if "gene_coverage_bp" in out else pd.Series(np.nan, index=out.index)
    te_cov_series   = out["te_coverage_bp"]   if "te_coverage_bp"   in out else pd.Series(np.nan, index=out.index)
    tm_cov_series   = out["tmrca_coverage_bp"] if "tmrca_coverage_bp" in out else pd.Series(np.nan, index=out.index)
    top_cov_series  = out["top_tmrca_coverage_bp"] if "top_tmrca_coverage_bp" in out else pd.Series(np.nan, index=out.index)

    out["te_gene_coverage_ratio"] = np.where(gene_cov_series > 0, te_cov_series / gene_cov_series, np.nan)
    denom_cov = gene_cov_series + te_cov_series
    out["te_fraction_of_covered_bp"] = np.where(denom_cov > 0, te_cov_series / denom_cov, np.nan)
    out["fraction_top_tmrca_bp"] = np.where(tm_cov_series > 0, top_cov_series / tm_cov_series, np.nan)

    # tmrca CV
    out["tmrca_cv"] = np.where(out["tmrca_mean"].notna() & (out["tmrca_mean"] != 0),
                               out["tmrca_sd"] / out["tmrca_mean"], np.nan)

    # ---------- Also output a combined row for scaffold_1a + scaffold_1b = scaffold_1 ----------
    def add_combined_row(out_df: pd.DataFrame) -> pd.DataFrame:
        a, b, new = "scaffold_1a", "scaffold_1b", "scaffold_1"
        if (a not in out_df["scaffold"].values) and (b not in out_df["scaffold"].values):
            return out_df

        row_a = out_df[out_df["scaffold"]==a].iloc[0] if a in out_df["scaffold"].values else None
        row_b = out_df[out_df["scaffold"]==b].iloc[0] if b in out_df["scaffold"].values else None

        def sget(r, col, default=np.nan):
            return r[col] if (r is not None and col in r.index and pd.notna(r[col])) else default

        combined = {"scaffold": new}
        # sum integer counts/coverages
        for col in ["n_tmrca_segments","tmrca_coverage_bp","n_top_tmrca_segments",
                    "n_non_top_tmrca_segments","top_tmrca_coverage_bp",
                    "n_genes","gene_coverage_bp","n_TEs_total","n_TE_classes","te_coverage_bp"]:
            va = sget(row_a, col, 0); vb = sget(row_b, col, 0)
            if pd.isna(va): va = 0
            if pd.isna(vb): vb = 0
            combined[col] = int(va) + int(vb)

        # pooled TMRCA stats: take weighted (by counts) for mean; min/max and quantiles from pooled NaN-robust concat
        # For simplicity, recompute from original TMRCA table (fewer edge cases)
        # (If either a or b missing in tdf, this still works.)
        mask = tdf[args.tmrca_chrom].isin([a,b])
        tm_vals = tdf.loc[mask, args.tmrca_mean].dropna().to_numpy()
        if tm_vals.size > 0:
            combined.update({
                "tmrca_min": float(np.nanmin(tm_vals)),
                "tmrca_q1":  q(tm_vals, 0.25),
                "tmrca_median": float(np.nanmedian(tm_vals)),
                "tmrca_q3":  q(tm_vals, 0.75),
                "tmrca_max": float(np.nanmax(tm_vals)),
                "tmrca_mean": float(np.nanmean(tm_vals)),
                "tmrca_sd": float(np.nanstd(tm_vals, ddof=1)) if tm_vals.size>=2 else np.nan,
                "tmrca_cv": (float(np.nanstd(tm_vals, ddof=1))/float(np.nanmean(tm_vals))) if (tm_vals.size>=2 and np.nanmean(tm_vals)!=0) else np.nan
            })
        else:
            for c in ["tmrca_min","tmrca_q1","tmrca_median","tmrca_q3","tmrca_max","tmrca_mean","tmrca_sd","tmrca_cv"]:
                combined[c] = np.nan

        # seg length stats (pooled)
        seg_vals = tdf.loc[mask, "seg_len"].dropna().to_numpy()
        combined["median_seg_len"] = float(np.nanmedian(seg_vals)) if seg_vals.size>0 else np.nan
        combined["mean_seg_len"]   = float(np.nanmean(seg_vals))   if seg_vals.size>0 else np.nan

        # derived ratios
        combined["te_gene_coverage_ratio"] = (combined["te_coverage_bp"]/combined["gene_coverage_bp"]
                                              if combined["gene_coverage_bp"]>0 else np.nan)
        denom = combined["gene_coverage_bp"] + combined["te_coverage_bp"]
        combined["te_fraction_of_covered_bp"] = (combined["te_coverage_bp"]/denom) if denom>0 else np.nan
        combined["fraction_top_tmrca_bp"] = (combined["top_tmrca_coverage_bp"]/combined["tmrca_coverage_bp"]
                                             if combined["tmrca_coverage_bp"]>0 else np.nan)

        # TE-class wide columns: sum the per-class counts
        tecols = [c for c in out_df.columns if c.startswith("TEclass_")]
        for c in tecols:
            va = sget(row_a, c, 0); vb = sget(row_b, c, 0)
            if pd.isna(va): va = 0
            if pd.isna(vb): vb = 0
            combined[c] = int(va) + int(vb)

        new_row = pd.DataFrame([combined])
        cols = list(out_df.columns)
        for c in new_row.columns:
            if c not in cols: cols.append(c)
        new_row = new_row.reindex(columns=cols)
        out_ext = pd.concat([out_df.reindex(columns=cols), new_row.reindex(columns=cols)], ignore_index=True)
        return out_ext

    out = add_combined_row(out)

    # ---------- Column ordering ----------
    fixed = ["scaffold",
             "n_tmrca_segments","tmrca_min","tmrca_q1","tmrca_median","tmrca_q3","tmrca_max","tmrca_mean",
             "tmrca_sd","tmrca_cv","median_seg_len","mean_seg_len",
             "tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
             "top_tmrca_coverage_bp","fraction_top_tmrca_bp",
             "n_genes","gene_coverage_bp",
             "n_TEs_total","n_TE_classes","te_coverage_bp",
             "te_gene_coverage_ratio","te_fraction_of_covered_bp"]
    te_cols = sorted([c for c in out.columns if c.startswith("TEclass_")])
    other_cols = [c for c in out.columns if c not in fixed + te_cols]
    out = out[fixed + te_cols + other_cols]

    # ---------- Write ----------
    out_path = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    meta_path = re.sub(r"\.tsv$", ".meta.txt", out_path)
    with open(meta_path, "w") as fh:
        fh.write(f"Global top-percentile for high TMRCA: {args.top_percentile}%\n")
        fh.write(f"Threshold (mean_tmrca >=): {tdf[args.tmrca_mean].quantile(1.0 - float(args.top_percentile)/100.0)}\n")
        fh.write("Coverage metrics use union of intervals (1-based inclusive lengths).\n")
        fh.write("scaffold_1 row is the combined summary of scaffold_1a and scaffold_1b.\n")
    print(f"Wrote: {out_path}")
    print(f"Wrote: {meta_path}")

if __name__ == "__main__":
    main()

