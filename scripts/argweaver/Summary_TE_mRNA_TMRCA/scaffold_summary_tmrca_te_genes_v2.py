#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, pandas as pd, numpy as np
from typing import Dict, List

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
    Return per-scaffold union coverage length (bp), treating coordinates as 1-based inclusive:
    length for an interval = end - start + 1. Overlaps are merged and counted once.
    """
    out = {}
    for scf, sub in df.groupby(chrom_col, sort=False):
        s = pd.to_numeric(sub[start_col], errors="coerce").dropna().astype(np.int64).to_numpy()
        e = pd.to_numeric(sub[end_col], errors="coerce").dropna().astype(np.int64).to_numpy()
        if len(s) == 0 or len(e) == 0:
            out[str(scf)] = 0
            continue
        arr = np.stack([s, e], axis=1)
        arr = arr[np.argsort(arr[:,0])]
        total = 0
        cur_s, cur_e = None, None
        for a, b in arr:
            if a > b:  # guard
                a, b = b, a
            if cur_s is None:
                cur_s, cur_e = int(a), int(b)
            else:
                if a <= cur_e + 1:  # overlapping or directly adjacent
                    if b > cur_e:
                        cur_e = int(b)
                else:
                    total += (cur_e - cur_s + 1)
                    cur_s, cur_e = int(a), int(b)
        if cur_s is not None:
            total += (cur_e - cur_s + 1)
        out[str(scf)] = int(total)
    return pd.Series(out, name="union_bp")

def main():
    ap = argparse.ArgumentParser(
        description="Per-scaffold summary of TMRCA, genes, TE counts and coverage (union bp)."
    )
    # Inputs (defaults set to your paths)
    ap.add_argument("--tmrca",
        default="~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025/all_tmrca_corrected_position.tsv",
        help="TMRCA segments TSV (CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    ap.add_argument("--mrna",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt",
        help="mRNA annotation TSV (new_seqid, new_start, new_end, mrna_id)")
    ap.add_argument("--te",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv",
        help="TE annotation TSV (seqid, sequence_ontology, start, end)")
    ap.add_argument("--out",
        default="scaffold_summary_v2.tsv",
        help="Output TSV")

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

    args = ap.parse_args()

    # ---------- TMRCA ----------
    tdf = load_table(args.tmrca)
    require(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], "TMRCA file")
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    tdf[args.tmrca_mean]  = pd.to_numeric(tdf[args.tmrca_mean], errors="coerce")
    # Per-scaffold stats
    def q(x, p): 
        try: return float(np.nanquantile(x, p))
        except Exception: return np.nan
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

    # ---------- Genes ----------
    gdf = load_table(args.mrna)
    require(gdf, [args.mrna_chrom, args.mrna_start, args.mrna_end], "mRNA file")
    gdf[args.mrna_chrom] = gdf[args.mrna_chrom].astype(str)
    # unique gene (mrna_id) count per scaffold (fallback to row count if mrna_id missing)
    if args.mrna_id in gdf.columns:
        gdf[args.mrna_id] = gdf[args.mrna_id].astype(str)
        genes_cnt = (gdf.groupby(args.mrna_chrom)[args.mrna_id]
                        .nunique(dropna=True)
                        .reset_index()
                        .rename(columns={args.mrna_chrom:"scaffold", args.mrna_id:"n_genes"}))
    else:
        genes_cnt = (gdf.groupby(args.mrna_chrom)
                        .size()
                        .reset_index(name="n_genes")
                        .rename(columns={args.mrna_chrom:"scaffold"}))
    # Gene union coverage (bp)
    gcov_series = intervals_union_bp(
        gdf.rename(columns={args.mrna_chrom:"scaffold",
                            args.mrna_start:"start",
                            args.mrna_end:"end"}),
        "scaffold","start","end"
    ).rename("gene_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- TEs ----------
    tte = load_table(args.te)
    require(tte, [args.te_chrom, args.te_start, args.te_end], "TE file")
    if args.te_class not in tte.columns:
        tte[args.te_class] = "TE"
    tte[args.te_chrom] = tte[args.te_chrom].astype(str)
    tte[args.te_class] = tte[args.te_class].astype(str)

    # total TE count per scaffold
    te_tot = (tte.groupby(args.te_chrom)
                 .size()
                 .reset_index(name="n_TEs_total")
                 .rename(columns={args.te_chrom:"scaffold"}))

    # per-class counts, wide
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

    # TE union coverage (bp)
    tecov_series = intervals_union_bp(
        tte.rename(columns={args.te_chrom:"scaffold",
                            args.te_start:"start",
                            args.te_end:"end"}),
        "scaffold","start","end"
    ).rename("te_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    # ---------- Merge everything ----------
    out = (tmrca_stats
           .merge(genes_cnt, on="scaffold", how="left")
           .merge(gcov_series, on="scaffold", how="left")
           .merge(te_tot, on="scaffold", how="left")
           .merge(tecov_series, on="scaffold", how="left")
           .merge(te_wide, on="scaffold", how="left"))

    # fill count/coverage NaNs
    for c in ["n_genes","gene_coverage_bp","n_TEs_total","te_coverage_bp"]:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)

    # order columns nicely
    fixed = ["scaffold",
             "n_tmrca_segments", "tmrca_min", "tmrca_q1", "tmrca_median",
             "tmrca_q3", "tmrca_max", "tmrca_mean",
             "n_genes", "gene_coverage_bp",
             "n_TEs_total", "te_coverage_bp"]
    te_cols = sorted([c for c in out.columns if c.startswith("TEclass_")])
    cols = fixed + [c for c in te_cols if c not in fixed]
    # keep any extras at the end
    extras = [c for c in out.columns if c not in cols]
    out = out[cols + extras]

    # write
    out_path = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")
    print(f"Rows: {len(out)}  Cols: {len(out.columns)}")

if __name__ == "__main__":
    main()

