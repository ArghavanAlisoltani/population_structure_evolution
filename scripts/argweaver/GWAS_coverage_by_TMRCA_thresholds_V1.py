#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS_coverage_by_TMRCA_thresholds_V1.py
---------------------------------------
Augments tmrca_with_mrna.tsv with multiple TMRCA percentile flags (top 5,10,15,20,25)
and computes, for each threshold, the number and percentage of GWAS positions
that fall inside segments with mean_tmrca in the top X% of the genome.
Also (re)adds per-segment GWAS annotations (gwas_n, gwas_traits, gwas_snp_ids, gwas_any).

---
USAGE
python GWAS_coverage_by_TMRCA_thresholds_V1.py \
  --tmrca fisher_mrna_out_v2/tmrca_with_mrna.tsv \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/Aria_MTAG_summary.txt" \
  --outdir tmrca_gwas_multi_out
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import bisect

def load_tmrca(path, chrom_col, start_col, end_col, mean_col):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    for c in [chrom_col, start_col, end_col, mean_col]:
        if c not in df.columns:
            raise ValueError(f"Required TMRCA column '{c}' not found in {path}")
    df[start_col] = pd.to_numeric(df[start_col], errors='coerce').astype('Int64')
    df[end_col]   = pd.to_numeric(df[end_col], errors='coerce').astype('Int64')
    df = df.dropna(subset=[start_col, end_col]).copy()
    df[start_col] = df[start_col].astype(int)
    df[end_col]   = df[end_col].astype(int)
    if (df[end_col] < df[start_col]).any():
        raise ValueError("Found segments with end < start in TMRCA file.")
    return df

def load_gwas(path, scaffold_col, pos_col, trait_col, snp_col):
    try:
        df = pd.read_csv(path, sep='\t', low_memory=False)
    except Exception:
        df = pd.read_csv(path, delim_whitespace=True, low_memory=False, header=0)
    lower = {c.lower(): c for c in df.columns}
    def pick(name): return lower.get(name.lower(), name)
    scaffold_col = pick(scaffold_col); pos_col = pick(pos_col); trait_col = pick(trait_col); snp_col = pick(snp_col)
    for c in [scaffold_col, pos_col, trait_col, snp_col]:
        if c not in df.columns:
            raise ValueError(f"GWAS column '{c}' not found; have: {list(df.columns)}")
    df[pos_col] = pd.to_numeric(df[pos_col], errors='coerce').astype('Int64')
    df = df.dropna(subset=[pos_col]).copy()
    df[pos_col] = df[pos_col].astype(int)
    df[scaffold_col] = df[scaffold_col].astype(str)
    df[trait_col] = df[trait_col].astype(str)
    df[snp_col] = df[snp_col].astype(str)
    return df.rename(columns={scaffold_col:'g_chrom', pos_col:'g_pos', trait_col:'g_trait', snp_col:'g_snp'})[['g_chrom','g_pos','g_trait','g_snp']]

def build_gwas_index(gdf):
    idx = {}
    for chrom, sub in gdf.groupby('g_chrom'):
        sub = sub.sort_values('g_pos').reset_index(drop=True)
        idx[chrom] = {'pos': sub['g_pos'].to_numpy(), 'trait': sub['g_trait'].to_numpy(), 'snp': sub['g_snp'].to_numpy()}
    return idx

def annotate_segments_with_gwas(tdf, gidx, chrom_col, start_col, end_col):
    gwas_n = []; gwas_traits = []; gwas_snps = []
    for ch, st, en in zip(tdf[chrom_col].astype(str).to_numpy(), tdf[start_col].to_numpy(), tdf[end_col].to_numpy()):
        if ch in gidx:
            arr_pos = gidx[ch]['pos']; arr_trait = gidx[ch]['trait']; arr_snp = gidx[ch]['snp']
            L = bisect.bisect_left(arr_pos, int(st)); R = bisect.bisect_right(arr_pos, int(en))
            if L < R:
                traits = arr_trait[L:R].tolist(); snps = arr_snp[L:R].tolist()
                seen=set(); tlist=[]; [tlist.append(t) for t in traits if not (t in seen or seen.add(t))]
                seen2=set(); slist=[]; [slist.append(s) for s in snps if not (s in seen2 or seen2.add(s))]
                gwas_n.append(R-L); gwas_traits.append(';'.join(tlist)); gwas_snps.append(';'.join(slist))
            else:
                gwas_n.append(0); gwas_traits.append(''); gwas_snps.append('')
        else:
            gwas_n.append(0); gwas_traits.append(''); gwas_snps.append('')
    tdf['gwas_n'] = gwas_n; tdf['gwas_traits'] = gwas_traits; tdf['gwas_snp_ids'] = gwas_snps; tdf['gwas_any'] = tdf['gwas_n'] > 0
    return tdf

def compute_threshold_flags(tdf, mean_col, percents):
    for p in percents:
        cutoff = np.percentile(tdf[mean_col].to_numpy(), 100 - p)
        tdf[f'is_top_p{p}'] = tdf[mean_col] >= cutoff
    return tdf

def percent_coverage_by_threshold(tdf, gidx, chrom_col, start_col, end_col, percents):
    total_positions = sum(len(gidx[ch]['pos']) for ch in gidx)
    rows = []; trait_rows = []
    for p in percents:
        n_overlap = 0
        per_trait = {}
        for ch, arr in gidx.items():
            pos = arr['pos']; trait = arr['trait']
            sub = tdf[(tdf[chrom_col].astype(str) == ch) & (tdf[f'is_top_p{p}'])][[start_col, end_col]]
            if sub.empty or pos.size == 0:
                continue
            sub = sub.sort_values(start_col).to_numpy()
            overlapped_idx = set()
            for st, en in sub:
                L = bisect.bisect_left(pos, int(st)); R = bisect.bisect_right(pos, int(en))
                if L < R:
                    overlapped_idx.update(range(L, R))
            n_overlap += len(overlapped_idx)
            for idx in overlapped_idx:
                tr = trait[idx]; per_trait[tr] = per_trait.get(tr, 0) + 1
        pct = (100.0 * n_overlap / total_positions) if total_positions > 0 else 0.0
        rows.append({'threshold_percent': p, 'n_gwas_total': total_positions, 'n_overlap': n_overlap, 'pct_overlap': pct})
        for tr, cnt in sorted(per_trait.items()):
            trait_rows.append({'threshold_percent': p, 'trait': tr, 'n_overlap': cnt})
    return pd.DataFrame(rows), pd.DataFrame(trait_rows) if trait_rows else pd.DataFrame(columns=['threshold_percent','trait','n_overlap'])

def main():
    ap = argparse.ArgumentParser(description="GWAS coverage vs TMRCA thresholds + per-segment GWAS annotation.")
    ap.add_argument('--tmrca', required=True)
    ap.add_argument('--gwas', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--chrom-col', default='CHROM')
    ap.add_argument('--start-col', default='start_tmrca')
    ap.add_argument('--end-col', default='end_tmrca')
    ap.add_argument('--mean-col', default='mean_tmrca')
    ap.add_argument('--thresholds', default='5,10,15,20,25')
    ap.add_argument('--gwas-scaffold-col', default='scaffold')
    ap.add_argument('--gwas-pos-col', default='position')
    ap.add_argument('--gwas-trait-col', default='trait')
    ap.add_argument('--gwas-snp-col', default='SNP')
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tdf = load_tmrca(args.tmrca, args.chrom_col, args.start_col, args.end_col, args.mean_col)
    gdf = load_gwas(args.gwas, args.gwas_scaffold_col, args.gwas_pos_col, args.gwas_trait_col, args.gwas_snp_col)
    gidx = build_gwas_index(gdf)

    percents = [int(x.strip()) for x in args.thresholds.split(',') if x.strip()!='']
    percents = sorted(set(percents))

    tdf = annotate_segments_with_gwas(tdf, gidx, args.chrom_col, args.start_col, args.end_col)
    tdf = compute_threshold_flags(tdf, args.mean_col, percents)
    tdf.to_csv(outdir/'tmrca_with_mrna_gwas_multi_thr.tsv', sep='\t', index=False)

    cov_df, trait_df = percent_coverage_by_threshold(tdf, gidx, args.chrom_col, args.start_col, args.end_col, percents)
    cov_df.to_csv(outdir/'gwas_coverage_by_tmrca_thresholds.tsv', sep='\t', index=False)
    if not trait_df.empty:
        trait_df.to_csv(outdir/'gwas_coverage_by_threshold_and_trait.tsv', sep='\t', index=False)

    with open(outdir/'summary_multi_thresholds.txt', 'w') as f:
        for _, r in cov_df.iterrows():
            f.write(f"Top {int(r['threshold_percent']):>2d}%: {int(r['n_overlap'])}/{int(r['n_gwas_total'])} GWAS positions in top segments ({r['pct_overlap']:.2f}%).\n")

if __name__ == '__main__':
    main()
