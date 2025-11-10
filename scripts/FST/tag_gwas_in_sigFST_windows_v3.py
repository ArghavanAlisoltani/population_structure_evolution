#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Report mRNA IDs and TE classes for significant FST windows that overlap GWAS SNPs.

Inputs
------
--gwas   : GWAS table (tab-delimited) with: scaffold, position, trait, SNP
--wins   : Window table (tab-delimited) with: CHROM, WIN_START, WIN_END,
           q_poi_3vs4/q_poi_3vs5/q_poi_4vs5 and TE_* count columns
--ann    : mRNA annotation table (tab-delimited) with: new_seqid, new_start, new_end,
           and an mRNA ID column (auto-detected; override with --id_col if needed)
--q      : q-value threshold for significant FST windows (default 0.05)
--out_prefix : prefix for outputs

Outputs
-------
<out_prefix>.SNP_level_sigFST_windows_mRNA_TE.tsv
  columns: scaffold, position, trait, SNP, comparison,
           CHROM, WIN_START, WIN_END, gene_ids, gene_count, TE_present

<out_prefix>.Window_level_sigFST_summary.tsv
  columns: comparison, CHROM, WIN_START, WIN_END, gene_ids, gene_count,
           TE_present, n_GWAS_in_window
"""

import argparse
import sys
import pandas as pd
import numpy as np

def norm_scaf_str(s: pd.Series) -> pd.Series:
    # normalize scaffold labels; also fix common "scafold" typo
    return s.astype(str).str.strip().str.lower().str.replace("scafold_", "scaffold_", regex=False)

def detect_id_column(df: pd.DataFrame, user_col: str | None = None) -> str:
    if user_col and user_col in df.columns:
        return user_col
    # heuristic: prefer columns containing both 'mrna' and 'id', or 'gene' and 'id'
    cand = [c for c in df.columns if ("mrna" in c.lower() and "id" in c.lower())]
    if cand:
        return cand[0]
    cand = [c for c in df.columns if ("gene" in c.lower() and "id" in c.lower())]
    if cand:
        return cand[0]
    # last resort: any column named exactly 'mRNA_ID'/'mrna_id'/'gene_ID' variants
    for c in ["mRNA_ID","mrna_id","gene_ID","gene_id","ID","id"]:
        if c in df.columns:
            return c
    sys.exit("Could not detect mRNA ID column in the annotation; pass --id_col explicitly.")

def list_te_classes_present(win_row: pd.Series, te_class_cols: list[str]) -> str:
    if not te_class_cols:
        return ""
    vals = pd.to_numeric(win_row[te_class_cols], errors="coerce").fillna(0)
    present = [c for c, v in zip(te_class_cols, vals.values) if v > 0]
    # Return class names without the "TE_" prefix for readability
    present_clean = [p[3:] if p.startswith("TE_") else p for p in present]
    return ";".join(present_clean)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True, help="GWAS file (tab): scaffold, position, trait, SNP")
    ap.add_argument("--wins", required=True, help="Windows file (tab): CHROM, WIN_START, WIN_END, q_poi_* and TE_* counts")
    ap.add_argument("--ann",  required=True, help="mRNA annotation (tab): new_seqid, new_start, new_end, and an mRNA ID column")
    ap.add_argument("--id_col", default=None, help="Name of mRNA ID column in --ann (auto-detected if omitted)")
    ap.add_argument("--q", type=float, default=0.05, help="q-value threshold for sig FST windows (default 0.05)")
    ap.add_argument("--out_prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    # ---- Load inputs
    gw = pd.read_csv(args.gwas, sep="\t")
    wins = pd.read_csv(args.wins, sep="\t")
    ann = pd.read_csv(args.ann, sep="\t")

    # ---- Basic checks
    for col in ["scaffold","position"]:
        if col not in gw.columns:
            sys.exit(f"GWAS file must contain column: {col}")
    for col in ["CHROM","WIN_START","WIN_END"]:
        if col not in wins.columns:
            sys.exit(f"Windows file must contain column: {col}")
    for col in ["new_seqid","new_start","new_end"]:
        if col not in ann.columns:
            sys.exit("Annotation must contain columns: new_seqid, new_start, new_end")

    # ---- Normalize coords / types
    gw["scaffold_norm"] = norm_scaf_str(gw["scaffold"])
    wins["CHROM_norm"]  = norm_scaf_str(wins["CHROM"])
    ann["new_seqid_norm"] = norm_scaf_str(ann["new_seqid"])

    gw["position"]   = pd.to_numeric(gw["position"], errors="coerce").astype("Int64")
    wins["WIN_START"] = pd.to_numeric(wins["WIN_START"], errors="coerce").astype("Int64")
    wins["WIN_END"]   = pd.to_numeric(wins["WIN_END"], errors="coerce").astype("Int64")
    ann["new_start"]  = pd.to_numeric(ann["new_start"], errors="coerce").astype("Int64")
    ann["new_end"]    = pd.to_numeric(ann["new_end"], errors="coerce").astype("Int64")

    # ---- Determine comparisons present
    qcols_expected = ["q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"]
    qcols_present = [c for c in qcols_expected if c in wins.columns]
    if not qcols_present:
        sys.exit("No q_poi_* columns found in windows file (expected one or more of q_poi_3vs4, q_poi_3vs5, q_poi_4vs5).")

    # ---- TE class columns (counts)
    te_class_cols = [c for c in wins.columns
                     if c.startswith("TE_") and not c.startswith("TEdens_")]
    # Drop generic totals/densities if you want only classes
    drop_exact = set(["TE_count","TE_density"])
    te_class_cols = [c for c in te_class_cols if c not in drop_exact]

    # ---- mRNA ID column
    id_col = detect_id_column(ann, args.id_col)

    # ---- Build fast per-scaffold indexes for genes and GWAS
    # Sort annotation per scaffold by start for quick slicing
    ann = ann.sort_values(["new_seqid_norm","new_start","new_end"])
    g_by_scaf = {sc: sub[["new_start","new_end",id_col]].to_numpy()
                 for sc, sub in ann.groupby("new_seqid_norm", sort=False)}
    gw_by_scaf = {sc: sub for sc, sub in gw.groupby("scaffold_norm", sort=False)}

    snp_rows = []     # SNP-level output rows
    win_rows = []     # window-level summary rows (unique windows per comparison)

    # ---- Iterate comparisons
    for qc in qcols_present:
        comp = qc.replace("q_poi_", "")  # e.g. '3vs4'
        sig = wins[pd.to_numeric(wins[qc], errors="coerce") <= args.q].copy()
        if sig.empty:
            continue

        # Precompute TE_present string for each window
        if te_class_cols:
            sig["TE_present"] = sig.apply(lambda r: list_te_classes_present(r, te_class_cols), axis=1)
        else:
            sig["TE_present"] = ""

        # For each scaffold, map SNPs to overlapping windows
        for scaf, gw_sub in gw_by_scaf.items():
            subW = sig[sig["CHROM_norm"] == scaf]
            if subW.empty:
                continue

            # windows arrays
            w_st = subW["WIN_START"].astype(int).to_numpy()
            w_en = subW["WIN_END"].astype(int).to_numpy()
            # Optional: speed by sorting windows (already likely sorted)
            w_te = subW["TE_present"].to_numpy()
            # gene list cache per window index
            gene_cache = {}

            # Build gene array for this scaffold once
            genes_arr = g_by_scaf.get(scaf, None)  # columns: new_start, new_end, id_col
            # iterate SNPs
            for _, row in gw_sub.iterrows():
                p = int(row["position"])
                # mask windows containing p
                mask = (w_st <= p) & (p <= w_en)
                if not mask.any():
                    continue

                idxs = np.where(mask)[0]
                for k in idxs:
                    w_start = int(w_st[k]); w_end = int(w_en[k])
                    # get gene IDs overlapping this window (cache to avoid recompute)
                    if k not in gene_cache:
                        if genes_arr is None or genes_arr.size == 0:
                            gene_ids = []
                        else:
                            gs = genes_arr
                            gmask = (gs[:,0] <= w_end) & (gs[:,1] >= w_start)  # overlap
                            if gmask.any():
                                # ids can be any dtype; cast to str and join
                                gene_ids = [str(x) for x in gs[gmask][:,2].tolist()]
                            else:
                                gene_ids = []
                        gene_cache[k] = ";".join(gene_ids), len(gene_ids)
                    genes_str, genes_n = gene_cache[k]
                    te_str = str(w_te[k]) if isinstance(w_te[k], str) else ""

                    # SNP-level row
                    snp_rows.append({
                        "scaffold": row["scaffold"],
                        "position": int(row["position"]),
                        "trait": row.get("trait", ""),
                        "SNP": row.get("SNP", ""),
                        "comparison": comp,
                        "CHROM": str(subW.iloc[k]["CHROM"]),
                        "WIN_START": w_start,
                        "WIN_END": w_end,
                        "gene_ids": genes_str,
                        "gene_count": genes_n,
                        "TE_present": te_str
                    })

            # Window-level rows (unique windows on this scaffold)
            # Count GWAS SNPs in each window
            if not gw_sub.empty:
                snp_pos = gw_sub["position"].astype(int).to_numpy()
                counts = []
                genes_for_win = []
                for k in range(len(subW)):
                    n_in = int(((snp_pos >= w_st[k]) & (snp_pos <= w_en[k])).sum())
                    counts.append(n_in)
                    # ensure gene cache filled for summary as well
                    if k not in gene_cache:
                        if genes_arr is None or genes_arr.size == 0:
                            gene_ids = []
                        else:
                            gs = genes_arr
                            gmask = (gs[:,0] <= w_en[k]) & (gs[:,1] >= w_st[k])
                            if gmask.any():
                                gene_ids = [str(x) for x in gs[gmask][:,2].tolist()]
                            else:
                                gene_ids = []
                        gene_cache[k] = ";".join(gene_ids), len(gene_ids)
                    genes_str, genes_n = gene_cache[k]
                    win_rows.append({
                        "comparison": comp,
                        "CHROM": str(subW.iloc[k]["CHROM"]),
                        "WIN_START": int(w_st[k]),
                        "WIN_END": int(w_en[k]),
                        "gene_ids": genes_str,
                        "gene_count": genes_n,
                        "TE_present": str(w_te[k]) if isinstance(w_te[k], str) else "",
                        "n_GWAS_in_window": int(counts[-1])
                    })

    # ---- Build and write outputs
    snp_df = pd.DataFrame(snp_rows,
                          columns=["scaffold","position","trait","SNP","comparison",
                                   "CHROM","WIN_START","WIN_END","gene_ids","gene_count","TE_present"])
    win_df = pd.DataFrame(win_rows,
                          columns=["comparison","CHROM","WIN_START","WIN_END","gene_ids","gene_count","TE_present","n_GWAS_in_window"])

    out1 = f"{args.out_prefix}.SNP_level_sigFST_windows_mRNA_TE.tsv"
    out2 = f"{args.out_prefix}.Window_level_sigFST_summary.tsv"
    snp_df.to_csv(out1, sep="\t", index=False)
    win_df.to_csv(out2, sep="\t", index=False)

    print(f"Wrote SNP-level:    {out1}  (rows: {len(snp_df)})")
    print(f"Wrote window-level: {out2}  (rows: {len(win_df)})")

if __name__ == "__main__":
    main()

