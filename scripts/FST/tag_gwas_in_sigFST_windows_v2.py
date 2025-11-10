#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tag each GWAS SNP with whether it falls inside significant FST windows,
and also annotate the matched windows' boundaries and whether any TE or gene
overlaps those windows.

Inputs:
  --gwas  : GWAS table (tab-delimited) with columns: scaffold, position, trait, SNP
  --wins  : Window table (tab-delimited) with columns: CHROM, WIN_START, WIN_END, q_poi_3vs4/3vs5/4vs5, and TE/gene info
  --q     : FDR/Q threshold for significance on q_poi_* (default 0.05)
  --out_prefix : prefix for outputs

Outputs:
  <out_prefix>.GWAS_with_FSTwindow_flags_v2.tsv
     - one row per GWAS SNP
     - per comparison: sigFST_* (bool), n windows overlapped, concatenated bounds,
       and booleans for TE_any and GENE_any aggregated across matched windows

  <out_prefix>.sigFST_windows_summary_v2.tsv
     - one row per significant window per comparison
     - window boundaries, TE_any, GENE_any, count of GWAS SNPs inside
"""

import argparse
import sys
import pandas as pd
import numpy as np

def norm_scaf(x: pd.Series) -> pd.Series:
    return x.astype(str).str.strip().str.lower()

def detect_gene_any(df: pd.DataFrame) -> pd.Series:
    """
    Return a boolean Series indicating if any gene overlaps the window.
    Tries common columns: 'gene_count', 'mRNA_count', 'genes_count', etc.
    Falls back to 'gene_ids' non-empty if present; else False.
    """
    cand_counts = [c for c in df.columns
                   if c.lower() in ("gene_count","genes_count","mrna_count","mrna_counts","n_genes","gene_n")]
    if cand_counts:
        cc = pd.to_numeric(df[cand_counts[0]], errors="coerce").fillna(0)
        return cc > 0

    # fallback: gene_ids style
    if "gene_ids" in df.columns:
        s = df["gene_ids"].astype(str).str.strip().str.lower()
        bad = s.isna() | (s == "") | (s == "na") | (s == "nan")
        return ~bad
    # nothing available
    return pd.Series(False, index=df.index)

def detect_te_any(df: pd.DataFrame) -> pd.Series:
    """
    Return a boolean Series indicating if any TE overlaps the window.
    Prefer 'TE_count' if present; else any 'TE_*' count column (not TEdens_*).
    """
    if "TE_count" in df.columns:
        cc = pd.to_numeric(df["TE_count"], errors="coerce").fillna(0)
        return cc > 0

    te_cols = [c for c in df.columns
               if c.startswith("TE_") and not c.startswith("TEdens_")]
    if te_cols:
        sub = df[te_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        return (sub.sum(axis=1) > 0)
    # nothing available
    return pd.Series(False, index=df.index)

def collect_bounds_str(starts: np.ndarray, ends: np.ndarray) -> str:
    """Join window bounds as 'start-end' separated by ';'."""
    if starts.size == 0:
        return ""
    return ";".join([f"{int(s)}-{int(e)}" for s, e in zip(starts, ends)])

def annotate_one_comparison(gw: pd.DataFrame, sigW: pd.DataFrame, label_suffix: str):
    """
    For a single comparison:
      - gw: GWAS dataframe with ['scaffold_norm','position']
      - sigW: significant windows with columns ['CHROM_norm','WIN_START','WIN_END','TE_any','GENE_any']
      - label_suffix: e.g. '3vs4'
    Adds columns to gw in place:
      sigFST_<suffix> (bool)
      FSTwin_<suffix>_n (int)
      FSTwin_<suffix>_bounds (str; 'start-end;start-end;...')
      FSTwin_<suffix>_TE_any (bool)
      FSTwin_<suffix>_GENE_any (bool)
    Returns array of indices (0..nGWAS-1) that were inside any window.
    """
    lab = label_suffix
    flag_col   = f"sigFST_{lab}"
    n_col      = f"FSTwin_{lab}_n"
    b_col      = f"FSTwin_{lab}_bounds"
    te_col     = f"FSTwin_{lab}_TE_any"
    gene_col   = f"FSTwin_{lab}_GENE_any"

    # init
    gw[flag_col] = False
    gw[n_col]    = 0
    gw[b_col]    = ""
    gw[te_col]   = False
    gw[gene_col] = False

    # nothing to do
    if sigW.empty:
        return np.array([], dtype=int)

    # group by scaffold for speed
    hit_indices = []

    for chrom, subW in sigW.groupby("CHROM_norm"):
        # GWAS indices on this scaffold
        idx = np.where(gw["scaffold_norm"].values == chrom)[0]
        if idx.size == 0:
            continue

        pos = gw.loc[idx, "position"].astype(int).values
        starts = subW["WIN_START"].astype(int).values
        ends   = subW["WIN_END"].astype(int).values
        te_any = subW["TE_any"].values.astype(bool)
        ge_any = subW["GENE_any"].values.astype(bool)

        # For each SNP pos, find all windows that contain it
        for j, p in enumerate(pos):
            mask = (starts <= p) & (p <= ends)
            if not mask.any():
                continue
            hit_indices.append(idx[j])

            # Update aggregates
            gw.at[idx[j], flag_col] = True
            gw.at[idx[j], n_col] = int(mask.sum())

            # bounds list
            bstr = collect_bounds_str(starts[mask], ends[mask])
            gw.at[idx[j], b_col] = bstr

            # OR across overlapped windows
            gw.at[idx[j], te_col] = bool(te_any[mask].any())
            gw.at[idx[j], gene_col] = bool(ge_any[mask].any())

    return np.unique(hit_indices)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True, help="GWAS file (tab), columns: scaffold, position, trait, SNP")
    ap.add_argument("--wins", required=True, help="Windows file (tab), columns: CHROM, WIN_START, WIN_END, q_poi_*")
    ap.add_argument("--q", type=float, default=0.05, help="q-value threshold for significant FST windows (default 0.05)")
    ap.add_argument("--out_prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    # Load
    gw = pd.read_csv(args.gwas, sep="\t")
    wins = pd.read_csv(args.wins, sep="\t")

    # Basic checks
    required_gw = {"scaffold","position"}
    if not required_gw.issubset(gw.columns):
        sys.exit(f"GWAS file must contain columns: {required_gw}")

    required_w = {"CHROM","WIN_START","WIN_END"}
    if not required_w.issubset(wins.columns):
        sys.exit(f"Windows file must contain columns: {required_w}")

    # Normalize scaffolds
    gw["scaffold_norm"] = norm_scaf(gw["scaffold"])
    wins["CHROM_norm"]  = norm_scaf(wins["CHROM"])

    # Ensure ints
    gw["position"]      = pd.to_numeric(gw["position"], errors="coerce").astype("Int64")
    wins["WIN_START"]   = pd.to_numeric(wins["WIN_START"], errors="coerce").astype("Int64")
    wins["WIN_END"]     = pd.to_numeric(wins["WIN_END"], errors="coerce").astype("Int64")

    # Determine available comparisons
    qcols_order = ["q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"]
    qcols_present = [c for c in qcols_order if c in wins.columns]
    if not qcols_present:
        sys.exit("No q_poi_* columns found. Expected one or more of: q_poi_3vs4, q_poi_3vs5, q_poi_4vs5")

    # Precompute TE_any and GENE_any on the full windows table, then subset
    wins = wins.copy()
    wins["TE_any"]   = detect_te_any(wins)
    wins["GENE_any"] = detect_gene_any(wins)

    # Collect window-level summary rows
    win_summary_rows = []

    # Per comparison annotate
    for qc in qcols_present:
        suffix = qc.replace("q_poi_", "")      # e.g., '3vs4'
        # Significant windows
        sigW = wins.loc[pd.to_numeric(wins[qc], errors="coerce") <= args.q,
                        ["CHROM","CHROM_norm","WIN_START","WIN_END","TE_any","GENE_any"]].dropna(subset=["WIN_START","WIN_END"])

        # Tag GWAS
        annotate_one_comparison(gw, sigW, suffix)

        # Window-level summary: count how many GWAS SNPs fall in each sig window
        if not sigW.empty:
            # Count overlaps per window (vectorized-ish)
            counts = []
            for chrom, subW in sigW.groupby("CHROM_norm", sort=False):
                idx = np.where(gw["scaffold_norm"].values == chrom)[0]
                if idx.size == 0:
                    counts.extend([0]*len(subW))
                    continue
                pos = gw.loc[idx, "position"].astype(int).values
                s = subW["WIN_START"].astype(int).values
                e = subW["WIN_END"].astype(int).values
                n_in = []
                for k in range(len(subW)):
                    n_in.append(int(((pos >= s[k]) & (pos <= e[k])).sum()))
                counts.extend(n_in)

            tmp = sigW.copy()
            tmp["comparison"] = suffix
            tmp["n_GWAS_in_window"] = counts
            win_summary_rows.append(tmp[["comparison","CHROM","WIN_START","WIN_END","TE_any","GENE_any","n_GWAS_in_window"]])

    # sigFST_any across comparisons
    sig_cols = [c for c in gw.columns if c.startswith("sigFST_")]
    gw["sigFST_any"] = gw[sig_cols].any(axis=1) if sig_cols else False

    # Arrange columns: core + flags first
    # Collect the per-comparison columns in a nice order
    percomp_cols = []
    for qc in qcols_present:
        suf = qc.replace("q_poi_", "")
        percomp_cols.extend([
            f"sigFST_{suf}",
            f"FSTwin_{suf}_n",
            f"FSTwin_{suf}_bounds",
            f"FSTwin_{suf}_TE_any",
            f"FSTwin_{suf}_GENE_any",
        ])

    core = ["scaffold","position"]
    if "trait" in gw.columns: core.append("trait")
    if "SNP" in gw.columns: core.append("SNP")

    ordered = core + percomp_cols + ["sigFST_any"]
    rest = [c for c in gw.columns if c not in ordered and c != "scaffold_norm"]
    gw_out = gw[ordered + rest].copy()

    # Write SNP-level output
    snp_out = f"{args.out_prefix}.GWAS_with_FSTwindow_flags_v2.tsv"
    gw_out.to_csv(snp_out, sep="\t", index=False)

    # Write window-level output
    if win_summary_rows:
        win_df = pd.concat(win_summary_rows, ignore_index=True)
    else:
        win_df = pd.DataFrame(columns=["comparison","CHROM","WIN_START","WIN_END","TE_any","GENE_any","n_GWAS_in_window"])
    win_out = f"{args.out_prefix}.sigFST_windows_summary_v2.tsv"
    win_df.to_csv(win_out, sep="\t", index=False)

    print(f"Wrote SNP-level: {snp_out}")
    print(f"Wrote window-level: {win_out}")

if __name__ == "__main__":
    main()

