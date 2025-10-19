#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Augment FST windows with MTAG GWAS hits and summarize MTAG overlap among FST-significant windows.

Adds columns:
  - MTAG_GWAS_n_hits
  - MTAG_GWAS_positions
  - MTAG_GWAS_traits
  - MTAG_GWAS_SNPs
  - MTAG_GWAS_has_hit

Creates per-comparison summaries (3vs4, 3vs5, 4vs5) counting overlap among windows
with q_poi_<comp> <= --q_threshold, and writes detail tables.

Usage example (paths adapted below).


python add_mtag_gwas_to_windows_and_summaries_v2.py \
  -w fst_windows_with_TEcounts__augmented.tsv \
  -g "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/Aria_MTAG_summary.txt" \
  -o fst_windows__WITH_MTAG.tsv \
  --q_threshold 0.05 \
  --win_chrom CHROM --win_start WIN_START --win_end WIN_END \
  --gwas_pos position --gwas_scaf scaffold --gwas_trait trait --gwas_snp SNP \
  --summary_out MTAG_FST_sig_overlap_summary.tsv \
  --details_dir MTAG_FST_sig_details

"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

COMPS = ["3vs4", "3vs5", "4vs5"]

def to_num(s):
    return pd.to_numeric(s, errors="coerce")

def add_mtag_to_windows(win, gwas, win_chrom, win_start, win_end,
                        gwas_scaf, gwas_pos, gwas_trait, gwas_snp):
    # Ensure numeric positions
    win = win.copy()
    win[win_start] = to_num(win[win_start])
    win[win_end]   = to_num(win[win_end])
    gwas = gwas.copy()
    gwas[gwas_pos] = to_num(gwas[gwas_pos])

    # Drop missing coords
    win = win.dropna(subset=[win_start, win_end])
    gwas = gwas.dropna(subset=[gwas_pos])

    # Cast to int64 for comparisons
    win[win_start] = win[win_start].astype(np.int64)
    win[win_end]   = win[win_end].astype(np.int64)
    gwas[gwas_pos] = gwas[gwas_pos].astype(np.int64)

    # Pre-create outputs
    n = len(win)
    nhits  = np.zeros(n, dtype=np.int64)
    pos_s  = np.array([""]*n, dtype=object)
    trt_s  = np.array([""]*n, dtype=object)
    snp_s  = np.array([""]*n, dtype=object)
    hashit = np.zeros(n, dtype=np.int64)

    # Index GWAS by scaffold
    g_groups = {k: v for k, v in gwas.groupby(gwas_scaf, sort=False)}

    # Process windows per scaffold
    for chrom, wsub in win.groupby(win_chrom, sort=False):
        if chrom not in g_groups:
            continue
        gsub = g_groups[chrom].sort_values(gwas_pos)
        gpos = gsub[gwas_pos].to_numpy()
        gtr  = gsub[gwas_trait].astype(str).to_numpy()
        gsnp = gsub[gwas_snp].astype(str).to_numpy()

        wsub = wsub.sort_values(win_start)
        widx = wsub.index.to_numpy()
        ws   = wsub[win_start].to_numpy()
        we   = wsub[win_end].to_numpy()

        j = 0
        m = len(gpos)
        for k in range(len(widx)):
            S = int(ws[k]); E = int(we[k])
            while j < m and gpos[j] < S:
                j += 1
            jj = j
            hits = []
            while jj < m and gpos[jj] <= E:
                hits.append(jj); jj += 1
            if hits:
                nhits[widx[k]] = len(hits)
                pos_set   = sorted({int(gpos[t]) for t in hits})
                trait_set = sorted({gtr[t]       for t in hits})
                snp_set   = sorted({gsnp[t]      for t in hits})
                pos_s[widx[k]] = ";".join(map(str, pos_set))
                trt_s[widx[k]] = ";".join(trait_set)
                snp_s[widx[k]] = ";".join(snp_set)
                hashit[widx[k]] = 1

    # Drop any previous MTAG_* columns to avoid duplicates
    for c in list(win.columns):
        if c.startswith("MTAG_GWAS_"):
            del win[c]

    win["MTAG_GWAS_n_hits"]    = nhits
    win["MTAG_GWAS_positions"] = pos_s
    win["MTAG_GWAS_traits"]    = trt_s
    win["MTAG_GWAS_SNPs"]      = snp_s
    win["MTAG_GWAS_has_hit"]   = hashit

    return win

def per_comparison_summaries(win_aug, q_threshold, win_chrom, comp_dirs=True):
    """
    Build a tidy summary per comparison and (optionally) write detail tables.
    Returns a summary dataframe.
    """
    rows = []
    for comp in COMPS:
        qcol = f"q_poi_{comp}"
        if qcol not in win_aug.columns:
            raise SystemExit(f"Missing column: {qcol}")

        sub = win_aug.copy()
        sub["__is_sig"] = to_num(sub[qcol]) <= q_threshold
        sig = sub.loc[sub["__is_sig"]].copy()
        n_sig = len(sig)
        if n_sig == 0:
            rows.append(dict(
                comparison=comp,
                q_threshold=q_threshold,
                n_sig_windows=0,
                n_sig_with_MTAG=0,
                prop_sig_with_MTAG=0.0,
                total_MTAG_hits_in_sig=0,
                unique_MTAG_positions_in_sig=0,
                unique_MTAG_SNPs_in_sig=0,
                unique_MTAG_traits_in_sig=0
            ))
            continue

        with_hit = sig.loc[sig["MTAG_GWAS_has_hit"] == 1]
        n_with = int((sig["MTAG_GWAS_has_hit"] == 1).sum())
        prop   = (n_with / n_sig) if n_sig else 0.0

        # aggregate uniques across sig windows
        all_pos = set()
        all_snp = set()
        all_trt = set()
        total_hits = int(sig["MTAG_GWAS_n_hits"].sum())
        for _, r in with_hit.iterrows():
            if isinstance(r["MTAG_GWAS_positions"], str) and r["MTAG_GWAS_positions"]:
                all_pos.update(r["MTAG_GWAS_positions"].split(";"))
            if isinstance(r["MTAG_GWAS_SNPs"], str) and r["MTAG_GWAS_SNPs"]:
                all_snp.update(r["MTAG_GWAS_SNPs"].split(";"))
            if isinstance(r["MTAG_GWAS_traits"], str) and r["MTAG_GWAS_traits"]:
                all_trt.update(r["MTAG_GWAS_traits"].split(";"))

        rows.append(dict(
            comparison=comp,
            q_threshold=q_threshold,
            n_sig_windows=int(n_sig),
            n_sig_with_MTAG=int(n_with),
            prop_sig_with_MTAG=round(prop, 6),
            total_MTAG_hits_in_sig=total_hits,
            unique_MTAG_positions_in_sig=len(all_pos - {""}),
            unique_MTAG_SNPs_in_sig=len(all_snp - {""}),
            unique_MTAG_traits_in_sig=len(all_trt - {""})
        ))

    return pd.DataFrame(rows).sort_values("comparison")

def write_per_comp_details(win_aug, outdir, q_threshold, win_chrom):
    """
    For each comparison, write detail tables of FST-significant windows that have MTAG hits.
    Also write a scaffold-level summary.
    """
    for comp in COMPS:
        comp_dir = Path(outdir, comp)
        comp_dir.mkdir(parents=True, exist_ok=True)
        qcol = f"q_poi_{comp}"

        sub = win_aug.copy()
        sub["__is_sig"] = to_num(sub[qcol]) <= q_threshold
        sig = sub.loc[sub["__is_sig"]].copy()

        # Detail of sig windows with hits
        detail = sig.loc[sig["MTAG_GWAS_has_hit"] == 1, [
            win_chrom, "WIN_START", "WIN_END", qcol,
            "MTAG_GWAS_n_hits", "MTAG_GWAS_positions", "MTAG_GWAS_traits", "MTAG_GWAS_SNPs"
        ]].sort_values([win_chrom, "WIN_START", "WIN_END"])
        detail.to_csv(comp_dir / f"fst_sig_windows_with_MTAG_{comp}.tsv", sep="\t", index=False)

        # Per-scaffold counts within sig windows
        scaf_summary = (sig
            .assign(has_hit=(sig["MTAG_GWAS_has_hit"] == 1).astype(int))
            .groupby(win_chrom, as_index=False)
            .agg(n_sig_windows=("__is_sig", "sum"),
                 n_sig_with_MTAG=("has_hit", "sum"),
                 total_MTAG_hits_in_sig=("MTAG_GWAS_n_hits", "sum"))
            .sort_values("n_sig_with_MTAG", ascending=False))
        scaf_summary.to_csv(comp_dir / f"per_scaffold_MTAG_in_FSTsig_{comp}.tsv", sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser(description="Add MTAG GWAS columns and summarize overlap within FST-significant windows.")
    ap.add_argument("-w", "--windows", required=True, help="Windows TSV (e.g., fst_windows_with_TEcounts__augmented.tsv)")
    ap.add_argument("-g", "--gwas", required=True, help="MTAG GWAS file (tab-delimited)")
    ap.add_argument("-o", "--out", required=True, help="Output augmented TSV with MTAG_* columns")
    ap.add_argument("--q_threshold", type=float, default=0.05, help="Poisson q threshold for FST significance (default 0.05)")
    # Windows columns
    ap.add_argument("--win_chrom", default="CHROM")
    ap.add_argument("--win_start", default="WIN_START")
    ap.add_argument("--win_end",   default="WIN_END")
    # GWAS columns
    ap.add_argument("--gwas_pos",   default="position")
    ap.add_argument("--gwas_scaf",  default="scaffold")
    ap.add_argument("--gwas_trait", default="trait")
    ap.add_argument("--gwas_snp",   default="SNP")
    # Extra outputs
    ap.add_argument("--summary_out", default=None, help="Optional: path for the per-comparison summary TSV")
    ap.add_argument("--details_dir", default=None, help="Optional: directory to write per-comparison detail tables")
    args = ap.parse_args()

    # Load
    win  = pd.read_csv(args.windows, sep="\t", low_memory=False)
    gwas = pd.read_csv(args.gwas,    sep="\t", low_memory=False)

    # Sanity checks
    for c in [args.win_chrom, args.win_start, args.win_end]:
        if c not in win.columns:
            raise SystemExit(f"Windows missing column: {c}")
    for c in [args.gwas_scaf, args.gwas_pos, args.gwas_trait, args.gwas_snp]:
        if c not in gwas.columns:
            raise SystemExit(f"GWAS missing column: {c}")
    for comp in COMPS:
        if f"q_poi_{comp}" not in win.columns:
            raise SystemExit(f"Windows missing column: q_poi_{comp}")

    # Augment windows with MTAG columns
    win_aug = add_mtag_to_windows(
        win, gwas,
        win_chrom=args.win_chrom, win_start=args.win_start, win_end=args.win_end,
        gwas_scaf=args.gwas_scaf, gwas_pos=args.gwas_pos,
        gwas_trait=args.gwas_trait, gwas_snp=args.gwas_snp
    )

    # Write augmented windows
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    win_aug.to_csv(args.out, sep="\t", index=False)
    print(f"[done] Augmented windows written: {args.out}")

    # Build and write per-comparison summary
    summary = per_comparison_summaries(win_aug, args.q_threshold, args.win_chrom)
    if args.summary_out:
        Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary_out, sep="\t", index=False)
        print(f"[done] Per-comparison summary written: {args.summary_out}")
    else:
        # If no path given, also print to stdout
        print("\n=== MTAG overlap among FST-significant windows ===")
        print(summary.to_string(index=False))

    # Write per-comparison details (optional)
    if args.details_dir:
        ddir = Path(args.details_dir); ddir.mkdir(parents=True, exist_ok=True)
        write_per_comp_details(win_aug, ddir, args.q_threshold, args.win_chrom)
        print(f"[done] Per-comparison detail tables written under: {ddir}")

if __name__ == "__main__":
    main()

