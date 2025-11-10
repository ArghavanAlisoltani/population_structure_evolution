#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
v4: Tag GWAS SNPs in significant FST windows (per comparison), attach window bounds, TE presence,
overlapping mRNA IDs, AND run GO enrichment for genes in sig windows that also overlap GWAS
vs the entire mRNA set.

Inputs
------
--gwas       : GWAS table (tab) with columns: scaffold, position, [trait], [SNP]
--wins       : windows table (tab) with: CHROM, WIN_START, WIN_END, q_poi_3vs4/3vs5/4vs5,
               TE_* class count columns (optional), and (optionally) gene overlap info (not required)
--ann        : mRNA annotation (tab) with: new_seqid, new_start, new_end, mRNA ID column, and a GO column
--id_col     : (optional) name of mRNA ID column in --ann (auto-detected if omitted)
--go_col     : (optional) name of GO column in --ann (auto-detected if omitted)
--go_sep     : (optional) GO term separator in --ann, default auto-detect among [';', '|', ',', ' ']
--q          : q-value threshold for FST significance on q_poi_* (default 0.05)
--out_prefix : output prefix for generated files

Outputs
-------
<out_prefix>.GWAS_with_FSTwindow_flags_v2.tsv
<out_prefix>.sigFST_windows_summary_v2.tsv
<out_prefix>.GO_enrichment_<comp>.tsv
<out_prefix>.GO_enrichment_any.tsv
"""

import argparse
import sys
import pandas as pd
import numpy as np
from math import comb

# ---------- Helpers ----------
def norm_scaf(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().str.replace("scafold_", "scaffold_", regex=False)

def detect_id_column(df: pd.DataFrame, user_col: str | None) -> str:
    if user_col and user_col in df.columns:
        return user_col
    cands = [c for c in df.columns if ("mrna" in c.lower() and "id" in c.lower())]
    if cands: return cands[0]
    cands = [c for c in df.columns if ("gene" in c.lower() and "id" in c.lower())]
    if cands: return cands[0]
    for c in ["mRNA_ID","mrna_id","gene_ID","gene_id","ID","id"]:
        if c in df.columns: return c
    sys.exit("Cannot detect mRNA ID column in annotation. Pass --id_col.")

def detect_go_column(df: pd.DataFrame, user_col: str | None) -> str:
    if user_col and user_col in df.columns:
        return user_col
    # prefer columns that look like GO lists
    cands = [c for c in df.columns if c.lower() in ("go","go_terms","go_bp","go_biological_process")]
    if cands: return cands[0]
    # any column containing 'go'
    cands = [c for c in df.columns if "go" in c.lower()]
    if cands: return cands[0]
    sys.exit("Cannot detect GO column in annotation. Pass --go_col.")

def autodetect_sep(sample: str) -> str:
    for sep in [';', '|', ',', ' ']:
        if sep in sample:
            return sep
    return ';'

def list_te_classes_present(win_row: pd.Series, te_class_cols: list[str]) -> str:
    if not te_class_cols:
        return ""
    vals = pd.to_numeric(win_row[te_class_cols], errors="coerce").fillna(0)
    present = [c for c, v in zip(te_class_cols, vals.values) if v > 0]
    # strip "TE_" prefix
    return ";".join([p[3:] if p.startswith("TE_") else p for p in present])

# One-sided Fisher (over-representation) via hypergeometric tail
def hypergeom_pmf(k, K, N, n):
    # K successes in population N; draw n; prob of k successes
    if k < 0 or k > n or k > K or n > N or K > N: return 0.0
    return comb(K, k) * comb(N - K, n - k) / comb(N, n)

def fisher_overrep_p(a, b, c, d):
    # contingency:
    # [a b]
    # [c d]
    N = a + b + c + d
    n = a + b         # foreground size
    K = a + c         # background successes
    # P(X >= a)
    p = 0.0
    kmax = min(n, K)
    for k in range(a, kmax + 1):
        p += hypergeom_pmf(k, K, N, n)
    return min(1.0, p)

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(m, dtype=float)
    ranked[order] = pvals[order] * m / (np.arange(m) + 1)
    # cumulative minimum from the right
    for i in range(m - 2, -1, -1):
        ranked[order[i]] = min(ranked[order[i]], ranked[order[i + 1]])
    ranked = np.clip(ranked, 0, 1)
    return ranked

# ---------- Core ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--wins", required=True)
    ap.add_argument("--ann",  required=True)
    ap.add_argument("--id_col", default=None)
    ap.add_argument("--go_col", default=None)
    ap.add_argument("--go_sep", default=None, help="GO separator; auto-detected if omitted")
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    # Load inputs
    gw = pd.read_csv(args.gwas, sep="\t")
    wins = pd.read_csv(args.wins, sep="\t")
    ann = pd.read_csv(args.ann, sep="\t")

    # Basic checks
    for col in ["scaffold","position"]:
        if col not in gw.columns:
            sys.exit(f"GWAS must contain: scaffold, position (missing {col})")
    for col in ["CHROM","WIN_START","WIN_END"]:
        if col not in wins.columns:
            sys.exit("Windows must contain: CHROM, WIN_START, WIN_END")
    for col in ["new_seqid","new_start","new_end"]:
        if col not in ann.columns:
            sys.exit("Annotation must contain: new_seqid, new_start, new_end")

    # Normalize scaffolds, ensure ints
    gw["scaffold_norm"] = norm_scaf(gw["scaffold"])
    wins["CHROM_norm"]  = norm_scaf(wins["CHROM"])
    ann["new_seqid_norm"] = norm_scaf(ann["new_seqid"])

    gw["position"]    = pd.to_numeric(gw["position"], errors="coerce").astype("Int64")
    wins["WIN_START"] = pd.to_numeric(wins["WIN_START"], errors="coerce").astype("Int64")
    wins["WIN_END"]   = pd.to_numeric(wins["WIN_END"], errors="coerce").astype("Int64")
    ann["new_start"]  = pd.to_numeric(ann["new_start"], errors="coerce").astype("Int64")
    ann["new_end"]    = pd.to_numeric(ann["new_end"], errors="coerce").astype("Int64")

    # Detect columns
    id_col = detect_id_column(ann, args.id_col)
    go_col = detect_go_column(ann, args.go_col)

    # GO separator
    if args.go_sep:
        go_sep = args.go_sep
    else:
        # peek a non-null value
        sample = ann[go_col].dropna().astype(str)
        go_sep = autodetect_sep(sample.iloc[0]) if not sample.empty else ';'

    # Build gene->GO mapping
    ann_go = ann[[id_col, go_col]].copy()
    ann_go[go_col] = ann_go[go_col].astype(str)
    # expand into lists
    ann_go_list = {}
    for gid, terms in ann_go.values:
        if pd.isna(terms) or str(terms).strip().lower() in ("", "na", "nan"):
            ann_go_list[str(gid)] = set()
        else:
            parts = [t.strip() for t in str(terms).split(go_sep) if t.strip() != ""]
            ann_go_list[str(gid)] = set(parts)
    background_genes = set(str(x) for x in ann[id_col].astype(str).values)

    # Which comparisons exist
    qcols_present = [c for c in ["q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"] if c in wins.columns]
    if not qcols_present:
        sys.exit("No q_poi_* columns found in windows file.")

    # TE class columns
    te_class_cols = [c for c in wins.columns if c.startswith("TE_") and not c.startswith("TEdens_")]
    te_class_cols = [c for c in te_class_cols if c not in {"TE_count","TE_density"}]

    # Sort and index
    ann_sorted = ann.sort_values(["new_seqid_norm","new_start","new_end"])
    genes_by_scaf = {
        sc: sub[["new_start","new_end",id_col]].to_numpy()
        for sc, sub in ann_sorted.groupby("new_seqid_norm", sort=False)
    }
    gw_by_scaf = {sc: sub for sc, sub in gw.groupby("scaffold_norm", sort=False)}

    # Collect outputs (as in v2)
    snp_rows = []
    win_rows = []

    # Foreground gene sets per comparison and pooled
    fg_genes_per_comp = { }
    fg_genes_any = set()

    for qc in qcols_present:
        comp = qc.replace("q_poi_", "")  # e.g., 3vs4
        sig = wins[pd.to_numeric(wins[qc], errors="coerce") <= args.q].copy()
        if sig.empty:
            fg_genes_per_comp[comp] = set()
            continue

        # Precompute TE presence string
        if te_class_cols:
            sig["TE_present"] = sig.apply(lambda r: list_te_classes_present(r, te_class_cols), axis=1)
        else:
            sig["TE_present"] = ""

        comp_fg_genes = set()

        # Per scaffold
        for scaf, gw_sub in gw_by_scaf.items():
            subW = sig[sig["CHROM_norm"] == scaf]
            if subW.empty:
                continue
            w_st = subW["WIN_START"].astype(int).to_numpy()
            w_en = subW["WIN_END"].astype(int).to_numpy()
            w_te = subW["TE_present"].to_numpy()
            gene_cache = {}  # window index -> (gene_str, gene_n, gene_set)
            genes_arr = genes_by_scaf.get(scaf, None)

            # Iterate SNPs
            for _, row in gw_sub.iterrows():
                p = int(row["position"])
                mask = (w_st <= p) & (p <= w_en)
                if not mask.any():
                    continue
                idxs = np.where(mask)[0]
                for k in idxs:
                    w_start = int(w_st[k]); w_end = int(w_en[k])

                    if k not in gene_cache:
                        if genes_arr is None or genes_arr.size == 0:
                            gset = set()
                        else:
                            gs = genes_arr
                            gmask = (gs[:,0] <= w_end) & (gs[:,1] >= w_start)
                            gset = set(str(x) for x in gs[gmask][:,2].tolist()) if gmask.any() else set()
                        gene_cache[k] = (";".join(sorted(gset)), len(gset), gset)
                    genes_str, genes_n, gset = gene_cache[k]
                    comp_fg_genes |= gset
                    te_str = str(w_te[k]) if isinstance(w_te[k], str) else ""

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

            # window-level counts
            if not gw_sub.empty:
                snp_pos = gw_sub["position"].astype(int).to_numpy()
                for k in range(len(subW)):
                    n_in = int(((snp_pos >= w_st[k]) & (snp_pos <= w_en[k])).sum())
                    if k not in gene_cache:
                        if genes_arr is None or genes_arr.size == 0:
                            gset = set()
                        else:
                            gs = genes_arr
                            gmask = (gs[:,0] <= w_en[k]) & (gs[:,1] >= w_st[k])
                            gset = set(str(x) for x in gs[gmask][:,2].tolist()) if gmask.any() else set()
                        gene_cache[k] = (";".join(sorted(gset)), len(gset), gset)
                    genes_str, genes_n, gset = gene_cache[k]
                    win_rows.append({
                        "comparison": comp,
                        "CHROM": str(subW.iloc[k]["CHROM"]),
                        "WIN_START": int(w_st[k]),
                        "WIN_END": int(w_en[k]),
                        "gene_ids": genes_str,
                        "gene_count": genes_n,
                        "TE_present": str(w_te[k]) if isinstance(w_te[k], str) else "",
                        "n_GWAS_in_window": n_in
                    })

        fg_genes_per_comp[comp] = comp_fg_genes
        fg_genes_any |= comp_fg_genes

    # ---------- Write SNP-level and window-level (same as v2) ----------
    snp_df = pd.DataFrame(
        snp_rows,
        columns=["scaffold","position","trait","SNP","comparison",
                 "CHROM","WIN_START","WIN_END","gene_ids","gene_count","TE_present"]
    )
    win_df = pd.DataFrame(
        win_rows,
        columns=["comparison","CHROM","WIN_START","WIN_END","gene_ids","gene_count","TE_present","n_GWAS_in_window"]
    )

    snp_out = f"{args.out_prefix}.GWAS_with_FSTwindow_flags_v2.tsv"
    win_out = f"{args.out_prefix}.sigFST_windows_summary_v2.tsv"
    snp_df.to_csv(snp_out, sep="\t", index=False)
    win_df.to_csv(win_out, sep="\t", index=False)
    print(f"Wrote SNP-level:    {snp_out}  (rows: {len(snp_df)})")
    print(f"Wrote window-level: {win_out}  (rows: {len(win_df)})")

    # ---------- GO enrichment ----------
    def one_enrichment(fg_genes: set[str], label: str):
        fg = set(str(x) for x in fg_genes if x in background_genes)
        bg = set(background_genes)

        fg_size = len(fg)
        bg_size = len(bg)
        if fg_size == 0 or bg_size == 0:
            return pd.DataFrame(columns=[
                "comparison","GO","fg_count","bg_count","fg_size","bg_size",
                "p_over","q_over","enrichment","fg_example_genes"
            ])

        # Build counts per GO (unique per gene)
        # Background GO counts
        go_bg_counts = {}
        for gid in bg:
            for go in ann_go_list.get(gid, set()):
                go_bg_counts[go] = go_bg_counts.get(go, 0) + 1
        # Foreground GO counts
        go_fg_counts = {}
        for gid in fg:
            for go in ann_go_list.get(gid, set()):
                go_fg_counts[go] = go_fg_counts.get(go, 0) + 1

        # Union of all GO terms observed in background
        all_gos = set(go_bg_counts.keys())
        if not all_gos:
            return pd.DataFrame(columns=[
                "comparison","GO","fg_count","bg_count","fg_size","bg_size",
                "p_over","q_over","enrichment","fg_example_genes"
            ])

        records = []
        pvals = []
        # optional: a small gene list per GO (examples)
        # build reverse index for foreground genes
        go_to_fg_genes = {}
        for gid in fg:
            for go in ann_go_list.get(gid, set()):
                go_to_fg_genes.setdefault(go, set()).add(gid)

        for go in all_gos:
            a = go_fg_counts.get(go, 0)
            c = go_bg_counts.get(go, 0)
            b = fg_size - a
            d = bg_size - c
            p = fisher_overrep_p(a, b, c, d)
            pvals.append(p)
            enr = (a / fg_size) / (c / bg_size) if c > 0 and a > 0 else 0.0
            ex_genes = ";".join(list(sorted(go_to_fg_genes.get(go, set())))[:15])
            records.append([label, go, a, c, fg_size, bg_size, p, None, enr, ex_genes])

        res = pd.DataFrame(records, columns=[
            "comparison","GO","fg_count","bg_count","fg_size","bg_size",
            "p_over","q_over","enrichment","fg_example_genes"
        ])
        res["q_over"] = bh_fdr(res["p_over"].to_numpy())
        res = res.sort_values(["q_over","p_over"], ascending=[True, True]).reset_index(drop=True)
        return res

    # Per comparison
    for comp, fgset in fg_genes_per_comp.items():
        res = one_enrichment(fgset, comp)
        outp = f"{args.out_prefix}.GO_enrichment_{comp}.tsv"
        res.to_csv(outp, sep="\t", index=False)
        print(f"Wrote GO enrichment: {outp} (rows: {len(res)})")

    # Any comparison pooled
    res_any = one_enrichment(fg_genes_any, "any")
    out_any = f"{args.out_prefix}.GO_enrichment_any.tsv"
    res_any.to_csv(out_any, sep="\t", index=False)
    print(f"Wrote GO enrichment: {out_any} (rows: {len(res_any)})")

if __name__ == "__main__":
    main()

