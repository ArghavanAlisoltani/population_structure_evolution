#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate_GWAS_with_TMRCA_V2.py
------------------------------
Attach mean TMRCA to each GWAS row and label 'top'/'background'.
NEW: per-trait stats (mean, median, SD) over non-missing mean TMRCA.




run example
outdir="$HOME/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/gwas_tmrca_annotations"
mkdir -p "$outdir"

python Annotate_GWAS_with_TMRCA_V2.py \
  --tmrca  $HOME/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_corrected_position.tsv \
  --gwas   $HOME/Desktop/OSU_projects/conifers/LP/Soms_shared/Aria_MTAG_summary.txt \
  --gwas-scaffold-col scaffold \
  --gwas-pos-col position \
  --gwas-trait-col trait \
  --gwas-snp-col SNP \
  --top-percentile 5 \
  --out   "$outdir/gwas_mtag_with_tmrca.tsv"

"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import bisect

def load_table(path, sep='\t'):
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except Exception:
        return pd.read_csv(path, delim_whitespace=True, low_memory=False)

def require_cols(df, cols, where):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {where}: {missing}. Available: {list(df.columns)}")

def build_interval_index(df, chrom_col, start_col, end_col, payload_cols):
    idx = {}
    for chrom, sub in df.groupby(chrom_col):
        sub = sub.sort_values([start_col, end_col]).reset_index(drop=True)
        entry = {
            'starts': sub[start_col].to_numpy(dtype=np.int64),
            'ends':   sub[end_col].to_numpy(dtype=np.int64),
            'payload': {c: sub[c].to_numpy() for c in payload_cols}
        }
        idx[str(chrom)] = entry
    return idx

def point_query_first(entry, pos):
    starts = entry['starts']; ends = entry['ends']
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and ends[i] >= pos:
        return i
    j = i + 1
    if 0 <= j < len(starts) and starts[j] <= pos <= ends[j]:
        return j
    return -1

def main():
    ap = argparse.ArgumentParser(description="Annotate GWAS/MTAG positions with mean TMRCA and top/background; per-trait stats.")
    ap.add_argument('--tmrca', required=True)
    ap.add_argument('--gwas', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--top-percentile', type=float, default=5.0)

    ap.add_argument('--tmrca-chrom', default='CHROM')
    ap.add_argument('--tmrca-start', default='start_tmrca')
    ap.add_argument('--tmrca-end',   default='end_tmrca')
    ap.add_argument('--tmrca-mean',  default='mean_tmrca')

    ap.add_argument('--gwas-scaffold-col', default='scaffold')
    ap.add_argument('--gwas-pos-col',      default='position')
    ap.add_argument('--gwas-trait-col',    default=None)
    ap.add_argument('--gwas-snp-col',      default=None)

    args = ap.parse_args()

    tdf = load_table(args.tmrca)
    require_cols(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], 'TMRCA')
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    tdf[args.tmrca_start] = pd.to_numeric(tdf[args.tmrca_start], errors='coerce').astype('Int64')
    tdf[args.tmrca_end]   = pd.to_numeric(tdf[args.tmrca_end], errors='coerce').astype('Int64')
    tdf = tdf.dropna(subset=[args.tmrca_start, args.tmrca_end]).copy()
    tdf[args.tmrca_start] = tdf[args.tmrca_start].astype(int)
    tdf[args.tmrca_end]   = tdf[args.tmrca_end].astype(int)

    gdf = load_table(args.gwas)
    require_cols(gdf, [args.gwas_scaffold_col, args.gwas_pos_col], 'GWAS')
    gdf[args.gwas_scaffold_col] = gdf[args.gwas_scaffold_col].astype(str)
    gdf[args.gwas_pos_col]      = pd.to_numeric(gdf[args.gwas_pos_col], errors='coerce').astype('Int64')
    gdf = gdf.dropna(subset=[args.gwas_pos_col]).copy()
    gdf[args.gwas_pos_col]      = gdf[args.gwas_pos_col].astype(int)

    tmrca_idx = build_interval_index(tdf, args.tmrca_chrom, args.tmrca_start, args.tmrca_end, [args.tmrca_mean])
    cutoff = np.percentile(tdf[args.tmrca_mean].to_numpy(), 100 - args.top_percentile)

    mean_vals = []; groups = []
    for ch, pos in zip(gdf[args.gwas_scaffold_col].to_numpy(), gdf[args.gwas_pos_col].to_numpy()):
        chs = str(ch)
        mt = np.nan; grp = 'NA'
        entry = tmrca_idx.get(chs)
        if entry is not None:
            i = point_query_first(entry, int(pos))
            if i >= 0:
                mt = float(entry['payload'][args.tmrca_mean][i])
                grp = 'top' if mt >= cutoff else 'background'
        mean_vals.append(mt)
        groups.append(grp)

    out = gdf.copy()
    out['mean_tmrca_at_pos'] = mean_vals
    out['tmrca_group'] = groups

    outpath = Path(args.out); outpath.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outpath, sep='\t', index=False)

    # Summary
    n = len(out)
    n_has = int(np.isfinite(out['mean_tmrca_at_pos']).sum())
    n_top = int((out['tmrca_group'] == 'top').sum())
    with open(outpath.with_suffix('.summary.txt'), 'w') as f:
        f.write(f"GWAS rows: {n}\n")
        f.write(f"Annotated with segment mean TMRCA: {n_has} ({(n_has/n*100 if n else 0):.2f}%)\n")
        f.write(f"Top {args.top_percentile:.1f}% TMRCA: {n_top} ({(n_top/n*100 if n else 0):.2f}%)\n")

    # Per-trait stats
    if args.gwas_trait_col and args.gwas_trait_col in out.columns:
        trait_col = args.gwas_trait_col
        non_na = out[np.isfinite(out['mean_tmrca_at_pos'])].copy()
        if not non_na.empty:
            stats = (non_na
                     .groupby(trait_col)['mean_tmrca_at_pos']
                     .agg(n_with_tmrca='size',
                          mean_tmrca='mean',
                          median_tmrca='median',
                          sd_tmrca='std')
                     .reset_index())
            stats['sd_tmrca'] = stats['sd_tmrca'].fillna(0.0)
            stats_path = outpath.with_suffix('.trait_stats.tsv')
            stats.to_csv(stats_path, sep='\t', index=False)

            with open(outpath.with_suffix('.summary.txt'), 'a') as f:
                f.write("\nPer-trait mean/median/SD (n with TMRCA):\n")
                for _, r in stats.iterrows():
                    f.write(f"  {r[trait_col]}: mean={r['mean_tmrca']:.3f}, median={r['median_tmrca']:.3f}, SD={r['sd_tmrca']:.3f} (n={int(r['n_with_tmrca'])})\n")
        else:
            with open(outpath.with_suffix('.summary.txt'), 'a') as f:
                f.write("\nPer-trait stats: no rows with finite mean_tmrca_at_pos.\n")

if __name__ == '__main__':
    main()
