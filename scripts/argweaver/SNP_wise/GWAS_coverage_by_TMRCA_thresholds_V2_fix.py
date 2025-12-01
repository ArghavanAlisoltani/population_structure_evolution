#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS_coverage_by_TMRCA_thresholds_V2_fix.py
- Fixes KeyError: 'trait' when merging per-trait totals by renaming g_trait -> trait.
- Otherwise identical to V2.


python GWAS_coverage_by_TMRCA_thresholds_V2_fix.py \
  --tmrca all_tmrca_corrected_position.tsv \
  --gwas "/Users/aria/Desktop/OSU_projects/conifers/LP/Soms_shared/Aria_MTAG_summary.txt" \
  --outdir tmrca_gwas_multi_out_v2

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
    df[scaffold_col] = df[scaffold_col].astype(str).str.strip()
    df[trait_col] = df[trait_col].astype(str).str.strip()
    df[snp_col] = df[snp_col].astype(str).str.strip()
    # ensure uniqueness per (chrom,pos,trait)
    df = df.drop_duplicates(subset=[scaffold_col, pos_col, trait_col]).reset_index(drop=True)
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
                seen=set(); tlist=[]
                for t in traits:
                    if t not in seen: seen.add(t); tlist.append(t)
                seen2=set(); slist=[]
                for s in snps:
                    if s not in seen2: seen2.add(s); slist.append(s)
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

def percent_coverage_by_threshold(tdf, gidx, gdf, chrom_col, start_col, end_col, percents):
    total_positions = gdf.shape[0]
    # per-trait totals, with a proper 'trait' column for merging
    trait_totals = (gdf.groupby('g_trait').size()
                    .rename('n_trait_total')
                    .reset_index()
                    .rename(columns={'g_trait':'trait'}))

    rows = []
    trait_rows_counts = []

    for p in percents:
        n_overlap = 0
        trait_over_idx = {}

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
                tr = trait[idx]
                trait_over_idx.setdefault(tr, set()).add(idx)

        pct = (100.0 * n_overlap / total_positions) if total_positions > 0 else 0.0
        rows.append({'threshold_percent': p, 'n_gwas_total': total_positions, 'n_overlap': n_overlap, 'pct_overlap': pct})

        for tr, idxset in sorted(trait_over_idx.items()):
            trait_rows_counts.append({'threshold_percent': p, 'trait': tr, 'n_overlap': len(idxset)})

    trait_counts_df = pd.DataFrame(trait_rows_counts)
    if not trait_counts_df.empty:
        trait_cov = trait_counts_df.merge(trait_totals, on='trait', how='left')
        trait_cov['pct_overlap'] = (100.0 * trait_cov['n_overlap'] / trait_cov['n_trait_total'].replace(0, np.nan)).fillna(0.0)
    else:
        trait_cov = pd.DataFrame(columns=['threshold_percent','trait','n_overlap','n_trait_total','pct_overlap'])

    return pd.DataFrame(rows), trait_counts_df, trait_cov

def main():
    ap = argparse.ArgumentParser(description="GWAS coverage vs TMRCA thresholds (+ per-trait counts & percentages).")
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
    percents = sorted(set(int(x.strip()) for x in args.thresholds.split(',') if x.strip()!=''))

    tdf = annotate_segments_with_gwas(tdf, gidx, args.chrom_col, args.start_col, args.end_col)
    tdf = compute_threshold_flags(tdf, args.mean_col, percents)
    tdf.to_csv(outdir/'tmrca_with_mrna_gwas_multi_thr.tsv', sep='\t', index=False)

    overall_cov_df, trait_counts_df, trait_cov_df = percent_coverage_by_threshold(
        tdf, gidx, gdf, args.chrom_col, args.start_col, args.end_col, percents
    )
    overall_cov_df.to_csv(outdir/'gwas_coverage_by_tmrca_thresholds.tsv', sep='\t', index=False)
    trait_counts_df.to_csv(outdir/'gwas_coverage_by_threshold_and_trait.tsv', sep='\t', index=False)
    trait_cov_df.to_csv(outdir/'gwas_trait_coverage_by_threshold.tsv', sep='\t', index=False)

    if not trait_cov_df.empty:
        wide = trait_cov_df.pivot(index='trait', columns='threshold_percent', values='pct_overlap').fillna(0.0)
        wide = wide.sort_index(axis=1)
        wide.to_csv(outdir/'gwas_trait_coverage_by_threshold_wide.tsv', sep='\t')

    with open(outdir/'summary_multi_thresholds.txt', 'w') as f:
        for _, r in overall_cov_df.iterrows():
            f.write(f"Top {int(r['threshold_percent']):>2d}%: {int(r['n_overlap'])}/{int(r['n_gwas_total'])} GWAS positions in top segments ({r['pct_overlap']:.2f}%).\n")
        if not trait_cov_df.empty:
            f.write("\nPer-trait percentages (Top 5% if available):\n")
            p5 = 5 if 5 in set(overall_cov_df['threshold_percent']) else min(percents)
            sub = trait_cov_df[trait_cov_df['threshold_percent'] == p5]
            for _, r in sub.sort_values(['pct_overlap','trait'], ascending=[False, True]).iterrows():
                f.write(f"  {r['trait']}: {int(r['n_overlap'])}/{int(r['n_trait_total'])} ({r['pct_overlap']:.2f}%)\n")


if __name__ == '__main__':
    main()
