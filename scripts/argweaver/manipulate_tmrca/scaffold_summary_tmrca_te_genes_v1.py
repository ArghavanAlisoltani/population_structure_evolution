#!/usr/bin/env python3
'''
Run example
python scaffold_summary_tmrca_te_genes.py \
  --tmrca  "~/Desktop/OSU_projects/conifers/LP/ARGweaver/Nov_18_2025/all_tmrca_corrected_position.tsv" \
  --mrna   "~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt" \
  --te     "~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv" \
  --out    "~/Desktop/OSU_projects/conifers/LP/ARGweaver/Nov_18_2025/scaffold_summary.tsv"


'''
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np

def load_table(path, default_sep="\t"):
    path = os.path.expanduser(path)
    try:
        return pd.read_csv(path, sep=default_sep, low_memory=False)
    except Exception:
        return pd.read_csv(path, delim_whitespace=True, low_memory=False)

def require(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {label}: {missing}\nAvailable: {list(df.columns)}")

def sanitize(col):
    # make safe column names for TE classes
    col = re.sub(r"\s+", "_", str(col))
    col = re.sub(r"[^0-9A-Za-z_]+", "_", col)
    return col

def main():
    ap = argparse.ArgumentParser(
        description="Per-scaffold summary of TMRCA, genes, and TE classes."
    )
    ap.add_argument("--tmrca",
        default="~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025/all_tmrca_corrected_position.tsv",
        help="TMRCA segments TSV (columns: CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    ap.add_argument("--mrna",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt",
        help="mRNA annotation TSV (columns: new_seqid, new_start, new_end, mrna_id)")
    ap.add_argument("--te",
        default="~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv",
        help="TE annotation TSV (columns: seqid, sequence_ontology, start, end)")
    ap.add_argument("--out",
        default="scaffold_summary.tsv",
        help="Output TSV path")
    # optional overrides if headers differ
    ap.add_argument("--tmrca-chrom", default="CHROM")
    ap.add_argument("--tmrca-mean",  default="mean_tmrca")
    ap.add_argument("--mrna-chrom",  default="new_seqid")
    ap.add_argument("--mrna-id",     default="mrna_id")
    ap.add_argument("--te-chrom",    default="seqid")
    ap.add_argument("--te-class",    default="sequence_ontology")
    args = ap.parse_args()

    # --- TMRCA summary per scaffold ---
    tdf = load_table(args.tmrca)
    require(tdf, [args.tmrca_chrom, args.tmrca_mean], "TMRCA file")
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    tdf[args.tmrca_mean]  = pd.to_numeric(tdf[args.tmrca_mean], errors="coerce")
    tmrca_sum = (
        tdf.groupby(args.tmrca_chrom, as_index=False)
           .agg(n_tmrca_segments=(args.tmrca_mean, "size"),
                median_mean_tmrca=(args.tmrca_mean, "median"),
                max_mean_tmrca=(args.tmrca_mean, "max"))
           .rename(columns={args.tmrca_chrom: "scaffold"})
    )

    # --- gene counts per scaffold (unique mrna_id) ---
    gdf = load_table(args.mrna)
    require(gdf, [args.mrna_chrom], "mRNA file")
    gdf[args.mrna_chrom] = gdf[args.mrna_chrom].astype(str)
    if args.mrna_id in gdf.columns:
        gdf[args.mrna_id] = gdf[args.mrna_id].astype(str)
        gene_counts = (
            gdf.groupby(args.mrna_chrom, as_index=False)[args.mrna_id]
               .nunique(dropna=True)
               .rename(columns={args.mrna_chrom: "scaffold", args.mrna_id: "n_genes"})
        )
    else:
        # fallback: count rows per scaffold if mrna_id column missing
        gene_counts = (
            gdf.groupby(args.mrna_chrom, as_index=False)
               .size()
               .rename(columns={args.mrna_chrom: "scaffold", "size": "n_genes"})
        )

    # --- TE totals and per-class counts per scaffold ---
    tdf_te = load_table(args.te)
    require(tdf_te, [args.te_chrom], "TE file")
    # class column may be absent; treat as "TE"
    if args.te_class not in tdf_te.columns:
        tdf_te[args.te_class] = "TE"

    tdf_te[args.te_chrom] = tdf_te[args.te_chrom].astype(str)
    tdf_te[args.te_class] = tdf_te[args.te_class].astype(str)

    # total TEs per scaffold
    te_tot = (
        tdf_te.groupby(args.te_chrom, as_index=False)
              .size()
              .rename(columns={args.te_chrom: "scaffold", "size": "n_TEs_total"})
    )

    # per-class wide table
    te_class_counts = (
        tdf_te.groupby([args.te_chrom, args.te_class], as_index=False)
              .size()
              .rename(columns={args.te_chrom: "scaffold",
                               args.te_class: "TE_class",
                               "size": "count"})
    )
    # sanitize class names for columns
    te_class_counts["TE_col"] = te_class_counts["TE_class"].map(lambda x: "TEclass_" + sanitize(x))
    te_wide = (
        te_class_counts.pivot_table(index="scaffold", columns="TE_col", values="count", fill_value=0, aggfunc="sum")
                       .reset_index()
    )

    # --- merge all parts ---
    out = tmrca_sum.merge(gene_counts, on="scaffold", how="left") \
                   .merge(te_tot, on="scaffold", how="left") \
                   .merge(te_wide, on="scaffold", how="left")

    # fill missing counts with 0; keep medians/max as float
    count_cols = ["n_genes", "n_TEs_total"] + [c for c in out.columns if c.startswith("TEclass_")]
    for c in count_cols:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)

    # nice ordering
    fixed = ["scaffold", "n_tmrca_segments", "median_mean_tmrca", "max_mean_tmrca", "n_genes", "n_TEs_total"]
    te_cols_sorted = sorted([c for c in out.columns if c.startswith("TEclass_")])
    out = out[fixed + te_cols_sorted]

    # write TSV
    out_path = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")
    print(f"Columns: {list(out.columns)}")

if __name__ == "__main__":
    main()
