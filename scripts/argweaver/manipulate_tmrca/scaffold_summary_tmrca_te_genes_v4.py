#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Example Run


'''
import argparse, os, re
from typing import List, Tuple
import pandas as pd
import numpy as np

# ------------------------- helpers -------------------------

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

def q(x, p):
    try:
        return float(np.nanquantile(x, p))
    except Exception:
        return np.nan

def intervals_union_bp(df: pd.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pd.Series:
    """
    Per-scaffold union coverage length (bp), 1-based inclusive: end - start + 1.
    Overlaps/adjacencies merged once.
    """
    out = {}
    if df.empty:
        return pd.Series(out, name="union_bp")
    for scf, sub in df.groupby(chrom_col, sort=False):
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
                if a <= cur_e + 1:
                    if b > cur_e:
                        cur_e = int(b)
                else:
                    total += (cur_e - cur_s + 1)
                    cur_s, cur_e = int(a), int(b)
        if cur_s is not None:
            total += (cur_e - cur_s + 1)
        out[str(scf)] = int(total)
    return pd.Series(out, name="union_bp")

# ------------------------- TE hierarchy inference -------------------------

def infer_te_hierarchy(term: str) -> Tuple[str, str, str]:
    """
    Map a sequence_ontology string to (Class, Order, Superfamily).
    Rules are literature-based heuristics; adjust as needed for your annotation set.
    """
    if term is None or pd.isna(term):
        return ("Unclassified", "Unknown", "Unclassified")
    t = term.lower()

    # Quick helpers
    def any_in(s, keys): return any(k in s for k in keys)

    # --- Class I (Retrotransposons) ---
    if "retrotransposon" in t or any_in(t, ["line", "sine", "penelope", "dirs", "bel", "pao"]):
        # Orders & superfamilies
        if "ltr" in t or any_in(t, ["gypsy","copia","bel","pao"]):
            order = "LTR"
            if "gypsy" in t: sf = "Gypsy"
            elif "copia" in t: sf = "Copia"
            elif any_in(t, ["bel","pao","bel_pao","bel-pao"]): sf = "Bel-Pao"
            else: sf = "LTR_other"
            return ("Class I (Retrotransposon)", order, sf)
        if "line" in t or any_in(t, ["l1","jockey","cr1","rte","tx1","rnd"]):
            order = "LINE"
            if "l1" in t: sf = "L1"
            elif "rte" in t: sf = "RTE"
            elif "cr1" in t: sf = "CR1"
            elif "jockey" in t: sf = "Jockey"
            else: sf = "LINE_other"
            return ("Class I (Retrotransposon)", order, sf)
        if "sine" in t:
            return ("Class I (Retrotransposon)", "SINE", "SINE_other")
        if "penelope" in t:
            return ("Class I (Retrotransposon)", "Penelope", "Penelope")
        if "dirs" in t or "tyrosine_recombinase" in t:
            return ("Class I (Retrotransposon)", "DIRS", "DIRS")
        # fallthrough generic retro
        return ("Class I (Retrotransposon)", "Unknown", "Retro_other")

    # --- Class II (DNA transposons) ---
    if "transposon" in t or any_in(t, ["tir","helitron","maverick","polinton","crypton","is_element"]):
        if "helitron" in t:
            return ("Class II (DNA transposon)", "Helitron", "Helitron")
        if "crypton" in t:
            return ("Class II (DNA transposon)", "Crypton", "Crypton")
        if "maverick" in t or "polinton" in t:
            return ("Class II (DNA transposon)", "Maverick/Polinton", "Maverick/Polinton")
        if "tir" in t or "tn" in t or "transposon" in t:
            order = "TIR"
            if any_in(t, ["hat", "hobo", "ac"]): sf = "hAT"
            elif any_in(t, ["cacta","enspm","en-spm"]): sf = "CACTA/En-Spm"
            elif any_in(t, ["mutator","mudr"]): sf = "Mutator"
            elif any_in(t, ["harbinger","pif"]): sf = "PIF/Harbinger"
            elif any_in(t, ["tc1","mariner"]): sf = "Tc1/Mariner"
            elif "piggybac" in t: sf = "PiggyBac"
            else: sf = "TIR_other"
            return ("Class II (DNA transposon)", order, sf)
        return ("Class II (DNA transposon)", "Unknown", "DNA_other")

    # Unclassified / fragments / generic repeats
    if any_in(t, ["repeat_fragment","repeat", "unknown", "unclassified"]):
        return ("Unclassified", "Unknown", "Unclassified")

    # Default fallback
    return ("Unclassified", "Unknown", "Unclassified")

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser(
        description=("V4: Per-scaffold summary across TMRCA, genes, and TEs, "
                     "with TE Class→Order→Superfamily mapping and robust handling of scaffold_1 vs 1a/1b.")
    )
    # Inputs
    ap.add_argument("--tmrca",
        default="~/Desktop/OSU_projects/conifers/LP/ARGweaver/Nov_18_2025/all_tmrca_corrected_position.tsv",
        help="TMRCA TSV (CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    ap.add_argument("--mrna",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt",
        help="mRNA TSV (new_seqid, new_start, new_end, mrna_id)")
    ap.add_argument("--te",
        default="~/Desktop/OSU_projects/conifers/LP/EDTA/panTE/col8_readable.TEanno.cds.scaf1ab_v1.tsv",
        help="TE TSV (seqid, sequence_ontology, start, end, strand)")
    ap.add_argument("--out",
        default="scaffold_summary_v4.tsv",
        help="Output summary TSV")
    ap.add_argument("--te-hierarchy-out",
        default="TE_hierarchy_v4.tsv",
        help="Output TSV for TE hierarchy mapping")

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

    # High-TMRCA threshold
    ap.add_argument("--top-percentile", type=float, default=25.0,
                    help="Top percentile for high TMRCA (global).")

    args = ap.parse_args()

    # ---------- Load inputs ----------
    tdf = load_table(args.tmrca)
    require(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], "TMRCA file")
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    for c in (args.tmrca_start, args.tmrca_end, args.tmrca_mean):
        tdf[c] = pd.to_numeric(tdf[c], errors="coerce")
    tdf["seg_len"] = (tdf[args.tmrca_end] - tdf[args.tmrca_start] + 1).clip(lower=0)

    gdf = load_table(args.mrna)
    require(gdf, [args.mrna_chrom, args.mrna_start, args.mrna_end], "mRNA file")
    gdf[args.mrna_chrom] = gdf[args.mrna_chrom].astype(str)
    for c in (args.mrna_start, args.mrna_end):
        gdf[c] = pd.to_numeric(gdf[c], errors="coerce")
    if args.mrna_id not in gdf.columns:
        gdf[args.mrna_id] = np.nan
    else:
        gdf[args.mrna_id] = gdf[args.mrna_id].astype(str)

    tte = load_table(args.te)
    require(tte, [args.te_chrom, args.te_class, args.te_start, args.te_end], "TE file")
    tte[args.te_chrom] = tte[args.te_chrom].astype(str)
    tte[args.te_class] = tte[args.te_class].astype(str)
    for c in (args.te_start, args.te_end):
        tte[c] = pd.to_numeric(tte[c], errors="coerce")

    # ---------- TE hierarchy mapping ----------
    uniq_terms = pd.Series(sorted(tte[args.te_class].dropna().unique()), name="sequence_ontology")
    mapped = uniq_terms.to_frame()
    mapped[["Class","Order","Superfamily"]] = mapped["sequence_ontology"].apply(
        lambda s: pd.Series(infer_te_hierarchy(s))
    )
    # Save hierarchy mapping (so you can edit/tune)
    os.makedirs(os.path.dirname(os.path.expanduser(args.te_hierarchy_out)) or ".", exist_ok=True)
    mapped.to_csv(os.path.expanduser(args.te_hierarchy_out), sep="\t", index=False)

    # Attach hierarchy to TE rows
    tte = tte.merge(mapped, left_on=args.te_class, right_on="sequence_ontology", how="left")

    # ---------- Universe of scaffolds across all inputs ----------
    scaff_all = pd.Index(
        sorted(set(tdf[args.tmrca_chrom].dropna().unique())
             | set(gdf[args.mrna_chrom].dropna().unique())
             | set(tte[args.te_chrom].dropna().unique())),
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
           .rename(columns={args.tmrca_chrom:"scaffold"})
    )
    tmrca_sd_df = (
        tdf.groupby(args.tmrca_chrom)[args.tmrca_mean]
           .agg(lambda s: float(np.nanstd(s, ddof=1)) if s.notna().sum() >= 2 else np.nan)
           .reset_index()
           .rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_mean:"tmrca_sd"})
    )
    tmrca_cov = intervals_union_bp(
        tdf.rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_start:"start", args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})
    seg_len_stats = (tdf.groupby(args.tmrca_chrom)["seg_len"]
                        .agg(median_seg_len="median", mean_seg_len="mean")
                        .reset_index().rename(columns={args.tmrca_chrom:"scaffold"}))

    # High-TMRCA metrics
    top_thr = tdf[args.tmrca_mean].quantile(1.0 - float(args.top_percentile)/100.0)
    tdf["is_top"] = tdf[args.tmrca_mean] >= top_thr
    top_count = (tdf.groupby([args.tmrca_chrom,"is_top"]).size()
                    .unstack(fill_value=0).reset_index()
                    .rename(columns={args.tmrca_chrom:"scaffold"}))
    if True not in top_count.columns:  top_count[True]  = 0
    if False not in top_count.columns: top_count[False] = 0
    top_count = top_count.rename(columns={True:"n_top_tmrca_segments", False:"n_non_top_tmrca_segments"})
    top_cov = intervals_union_bp(
        tdf.loc[tdf["is_top"]].rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_start:"start", args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("top_tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- Genes ----------
    gene_counts = (gdf.groupby(args.mrna_chrom)[args.mrna_id]
                      .nunique(dropna=True)
                      .reset_index()
                      .rename(columns={args.mrna_chrom:"scaffold", args.mrna_id:"n_genes"}))
    gene_cov = intervals_union_bp(
        gdf.rename(columns={args.mrna_chrom:"scaffold", args.mrna_start:"start", args.mrna_end:"end"}),
        "scaffold","start","end"
    ).rename("gene_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- TEs ----------
    te_tot = (tte.groupby(args.te_chrom).size()
                 .reset_index(name="n_TEs_total")
                 .rename(columns={args.te_chrom:"scaffold"}))
    n_te_classes = (tte.groupby(args.te_chrom)[args.te_class]
                       .nunique().reset_index(name="n_TE_classes")
                       .rename(columns={args.te_chrom:"scaffold"}))
    te_cov = intervals_union_bp(
        tte.rename(columns={args.te_chrom:"scaffold", args.te_start:"start", args.te_end:"end"}),
        "scaffold","start","end"
    ).rename("te_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # per TE "sequence_ontology" class (wide, as in v3)
    te_class_counts = (tte.groupby([args.te_chrom, args.te_class])
                          .size().reset_index(name="count")
                          .rename(columns={args.te_chrom:"scaffold", args.te_class:"TE_class"}))
    if len(te_class_counts) > 0:
        te_class_counts["TE_col"] = te_class_counts["TE_class"].map(lambda x: "TEclass_" + sanitize(x))
        te_wide = (te_class_counts.pivot_table(index="scaffold", columns="TE_col", values="count",
                                               fill_value=0, aggfunc="sum").reset_index())
    else:
        te_wide = pd.DataFrame(columns=["scaffold"])

    # per TE *superfamily* counts (wide)
    te_sf_counts = (tte.groupby([args.te_chrom, "Superfamily"]).size()
                       .reset_index(name="count")
                       .rename(columns={args.te_chrom:"scaffold"}))
    if len(te_sf_counts) > 0:
        te_sf_counts["SF_col"] = te_sf_counts["Superfamily"].map(lambda x: "TEsf_" + sanitize(x))
        te_sf_wide = (te_sf_counts.pivot_table(index="scaffold", columns="SF_col", values="count",
                                               fill_value=0, aggfunc="sum").reset_index())
    else:
        te_sf_wide = pd.DataFrame(columns=["scaffold"])
    n_te_superfamilies = (tte.groupby(args.te_chrom)["Superfamily"].nunique()
                             .reset_index(name="n_TE_superfamilies")
                             .rename(columns={args.te_chrom:"scaffold"}))

    # 'repeat_fragment' diagnostics (optional but handy)
    repf_counts = (tte[tte[args.te_class].str.contains("repeat_fragment", case=False, na=False)]
                      .groupby(args.te_chrom).size()
                      .reset_index(name="repeat_fragment_count")
                      .rename(columns={args.te_chrom:"scaffold"}))
    repf_cov = intervals_union_bp(
        tte[tte[args.te_class].str.contains("repeat_fragment", case=False, na=False)]
           .rename(columns={args.te_chrom:"scaffold", args.te_start:"start", args.te_end:"end"}),
        "scaffold","start","end"
    ).rename("repeat_fragment_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- Merge everything on ALL scaffolds ----------
    out = (base
           .merge(tmrca_stats, on="scaffold", how="left")
           .merge(tmrca_sd_df, on="scaffold", how="left")
           .merge(tmrca_cov, on="scaffold", how="left")
           .merge(seg_len_stats, on="scaffold", how="left")
           .merge(top_count, on="scaffold", how="left")
           .merge(top_cov, on="scaffold", how="left")
           .merge(gene_counts, on="scaffold", how="left")
           .merge(gene_cov, on="scaffold", how="left")
           .merge(te_tot, on="scaffold", how="left")
           .merge(n_te_classes, on="scaffold", how="left")
           .merge(n_te_superfamilies, on="scaffold", how="left")
           .merge(te_cov, on="scaffold", how="left")
           .merge(repf_counts, on="scaffold", how="left")
           .merge(repf_cov, on="scaffold", how="left")
           .merge(te_wide, on="scaffold", how="left")
           .merge(te_sf_wide, on="scaffold", how="left")
    )

    # Fill counts/coverage
    for c in ["n_tmrca_segments","tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
              "top_tmrca_coverage_bp","n_genes","gene_coverage_bp","n_TEs_total",
              "n_TE_classes","n_TE_superfamilies","te_coverage_bp",
              "repeat_fragment_count","repeat_fragment_coverage_bp"]:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)
    if "tmrca_sd" not in out.columns:
        out["tmrca_sd"] = np.nan

    # Derived ratios/fractions
    def frac(a, b): 
        return np.where(b>0, a/b, np.nan)
    out["fraction_top_tmrca_bp"] = frac(out.get("top_tmrca_coverage_bp",0), out.get("tmrca_coverage_bp",0))
    out["te_gene_coverage_ratio"] = frac(out.get("te_coverage_bp",0), out.get("gene_coverage_bp",0))
    cov_sum = out.get("gene_coverage_bp",0) + out.get("te_coverage_bp",0)
    out["te_fraction_of_covered_bp"] = frac(out.get("te_coverage_bp",0), cov_sum)
    out["tmrca_cv"] = np.where(out["tmrca_mean"].notna() & (out["tmrca_mean"] != 0),
                               out["tmrca_sd"] / out["tmrca_mean"], np.nan)

    # ---------- Robust "scaffold_1" combination logic ----------
    def pooled_row_for(names: List[str]) -> pd.Series:
        sub = out[out["scaffold"].isin(names)].copy()
        if sub.empty:
            return None
        # numeric sums for counts/coverages
        sum_cols = ["n_tmrca_segments","tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
                    "top_tmrca_coverage_bp","n_genes","gene_coverage_bp","n_TEs_total","n_TE_classes",
                    "n_TE_superfamilies","te_coverage_bp","repeat_fragment_count","repeat_fragment_coverage_bp"]
        combined = {c:int(sub[c].fillna(0).sum()) for c in sum_cols if c in sub.columns}
        # pooled TMRCA distro from raw tdf
        mask = tdf[args.tmrca_chrom].isin(names)
        tm_vals = tdf.loc[mask, args.tmrca_mean].dropna().to_numpy()
        seg_vals = tdf.loc[mask, "seg_len"].dropna().to_numpy()
        combined.update({
            "tmrca_min": float(np.nanmin(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_q1":  q(tm_vals, 0.25) if tm_vals.size>0 else np.nan,
            "tmrca_median": float(np.nanmedian(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_q3":  q(tm_vals, 0.75) if tm_vals.size>0 else np.nan,
            "tmrca_max": float(np.nanmax(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_mean": float(np.nanmean(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_sd": float(np.nanstd(tm_vals, ddof=1)) if tm_vals.size>=2 else np.nan,
            "tmrca_cv": (float(np.nanstd(tm_vals, ddof=1))/float(np.nanmean(tm_vals))
                         if tm_vals.size>=2 and np.nanmean(tm_vals)!=0 else np.nan),
            "median_seg_len": float(np.nanmedian(seg_vals)) if seg_vals.size>0 else np.nan,
            "mean_seg_len":   float(np.nanmean(seg_vals))   if seg_vals.size>0 else np.nan,
        })
        # derived
        combined["fraction_top_tmrca_bp"] = (combined["top_tmrca_coverage_bp"]/combined["tmrca_coverage_bp"]
                                             if combined.get("tmrca_coverage_bp",0)>0 else np.nan)
        combined["te_gene_coverage_ratio"] = (combined["te_coverage_bp"]/combined["gene_coverage_bp"]
                                              if combined.get("gene_coverage_bp",0)>0 else np.nan)
        denom = combined.get("gene_coverage_bp",0) + combined.get("te_coverage_bp",0)
        combined["te_fraction_of_covered_bp"] = (combined.get("te_coverage_bp",0)/denom if denom>0 else np.nan)
        # sum wide TE class & superfamily columns
        for pref in ["TEclass_","TEsf_"]:
            cols = [c for c in out.columns if c.startswith(pref)]
            for c in cols:
                combined[c] = int(sub[c].fillna(0).sum()) if c in sub.columns else 0
        return pd.Series(combined)

    have_1 = "scaffold_1" in out["scaffold"].values
    have_1a = "scaffold_1a" in out["scaffold"].values
    have_1b = "scaffold_1b" in out["scaffold"].values
    if have_1a or have_1b:
        names_to_pool = [n for n in ["scaffold_1","scaffold_1a","scaffold_1b"] if n in out["scaffold"].values]
        pooled = pooled_row_for(names_to_pool)
        if pooled is not None:
            pooled["scaffold"] = "scaffold_1"
            # replace or append scaffold_1 row with pooled stats
            out = out[ out["scaffold"] != "scaffold_1" ]
            out = pd.concat([out, pooled.to_frame().T], ignore_index=True)

    # ---------- Final column order ----------
    fixed = ["scaffold",
             "n_tmrca_segments","tmrca_min","tmrca_q1","tmrca_median","tmrca_q3","tmrca_max","tmrca_mean",
             "tmrca_sd","tmrca_cv","median_seg_len","mean_seg_len",
             "tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
             "top_tmrca_coverage_bp","fraction_top_tmrca_bp",
             "n_genes","gene_coverage_bp",
             "n_TEs_total","n_TE_classes","n_TE_superfamilies","te_coverage_bp",
             "repeat_fragment_count","repeat_fragment_coverage_bp",
             "te_gene_coverage_ratio","te_fraction_of_covered_bp"]
    te_class_cols = sorted([c for c in out.columns if c.startswith("TEclass_")])
    te_sf_cols    = sorted([c for c in out.columns if c.startswith("TEsf_")])
    other_cols    = [c for c in out.columns if c not in fixed + te_class_cols + te_sf_cols]
    out = out[fixed + te_sf_cols + te_class_cols + other_cols]

    # ---------- Write outputs ----------
    out_path = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    hier_path = os.path.expanduser(args.te_hierarchy_out)
    print(f"Wrote summary:  {out_path}")
    print(f"Wrote hierarchy:{hier_path}")
    print(f"High-TMRCA threshold (mean_tmrca >=): {top_thr}")

if __name__ == "__main__":
    main()

