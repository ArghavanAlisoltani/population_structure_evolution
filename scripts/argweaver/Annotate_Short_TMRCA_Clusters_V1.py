#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate_Short_TMRCA_Clusters_V1.py
-----------------------------------
Identify and annotate *short* TMRCA segments that are *nearby* on the same scaffold
and assign them the same cluster ID.

Config (CLI):
  --short-max        (default 100) short if length < short-max (use --short-criterion le to make <=)
  --near-gap         (default 200) maximum allowed gap between consecutive short segments on a scaffold
  --chrom-col        (default CHROM)
  --start-col        (default start_tmrca)
  --end-col          (default end_tmrca)

Outputs:
  - tmrca_with_short_clusters.tsv  (per-segment annotations)
  - short_clusters.tsv             (cluster summary)
  - summary.txt


  ---
  Run example
  
  python Annotate_Short_TMRCA_Clusters_V1.py \
  --tmrca annotated_tmrca_4_GPT_13columns.tsv \
  --short-max 100 \
  --short-criterion lt \
  --near-gap 200 \
  --outdir short_tmrca_clusters_out

"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def load_tmrca(path, chrom_col, start_col, end_col):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    for c in (chrom_col, start_col, end_col):
        if c not in df.columns:
            raise ValueError(f"Required column '{c}' not found in {path}")
    df[chrom_col] = df[chrom_col].astype(str)
    df[start_col] = pd.to_numeric(df[start_col], errors='coerce').astype('Int64')
    df[end_col]   = pd.to_numeric(df[end_col], errors='coerce').astype('Int64')
    df = df.dropna(subset=[start_col, end_col]).copy()
    df[start_col] = df[start_col].astype(int)
    df[end_col]   = df[end_col].astype(int)
    if (df[end_col] < df[start_col]).any():
        raise ValueError("Found rows with end < start in TMRCA file.")
    return df

def annotate_short_clusters(df, chrom_col, start_col, end_col, short_max, short_criterion, near_gap):
    df = df.sort_values([chrom_col, start_col, end_col]).reset_index(drop=True).copy()
    df['seg_length'] = (df[end_col] - df[start_col] + 1).astype(int)

    if short_criterion == 'lt':
        is_short = df['seg_length'] < short_max
    else:
        is_short = df['seg_length'] <= short_max
    df['is_short'] = is_short

    cluster_id = np.full(len(df), np.nan, dtype=float)
    cluster_size = np.full(len(df), np.nan, dtype=float)
    cluster_start = np.full(len(df), np.nan, dtype=float)
    cluster_end = np.full(len(df), np.nan, dtype=float)
    cluster_span = np.full(len(df), np.nan, dtype=float)
    gap_to_prev = np.full(len(df), np.nan, dtype=float)

    next_cluster = 0

    for chrom, sub_idx in df.groupby(chrom_col, sort=False).groups.items():
        idxs = list(sub_idx)
        prev_short_i = None
        current_cluster = None
        current_min = None
        current_max = None
        current_members = []

        def flush_cluster():
            nonlocal current_cluster, current_min, current_max, current_members
            if current_cluster is None or len(current_members) == 0:
                current_cluster = None; current_min = None; current_max = None; current_members = []
                return
            size = len(current_members)
            span = (current_max - current_min + 1) if (current_min is not None and current_max is not None) else np.nan
            for i in current_members:
                cluster_id[i] = current_cluster
                cluster_size[i] = size
                cluster_start[i] = current_min
                cluster_end[i]   = current_max
                cluster_span[i]  = span
            current_cluster = None; current_min = None; current_max = None; current_members = []

        for i in idxs:
            if not is_short.iat[i]:
                prev_short_i = None
                flush_cluster()
                continue

            st = int(df.at[i, start_col]); en = int(df.at[i, end_col])
            if prev_short_i is None:
                current_cluster = next_cluster; next_cluster += 1
                current_min = st; current_max = en
                current_members = [i]
                prev_short_i = i
                gap_to_prev[i] = np.nan
            else:
                prev_en = int(df.at[prev_short_i, end_col])
                gap = st - prev_en - 1
                gap_to_prev[i] = gap
                if gap <= near_gap:
                    current_members.append(i)
                    prev_short_i = i
                    current_min = min(current_min, st)
                    current_max = max(current_max, en)
                else:
                    flush_cluster()
                    current_cluster = next_cluster; next_cluster += 1
                    current_min = st; current_max = en
                    current_members = [i]
                    prev_short_i = i

        flush_cluster()

    df['gap_to_prev_short_on_scaffold'] = gap_to_prev
    df['short_cluster_id'] = cluster_id
    df['short_cluster_size'] = cluster_size
    df['short_cluster_start'] = cluster_start
    df['short_cluster_end'] = cluster_end
    df['short_cluster_span'] = cluster_span
    return df

def summarize_clusters(df, chrom_col, start_col, end_col):
    x = df.dropna(subset=['short_cluster_id']).copy()
    if x.empty:
        return pd.DataFrame(columns=['short_cluster_id','scaffold','short_cluster_start','short_cluster_end',
                                     'short_cluster_span','n_short_segments','total_short_bases',
                                     'mean_gap_between_members','max_gap_between_members'])
    x['short_cluster_id'] = x['short_cluster_id'].astype(int)
    x_sorted = x.sort_values([chrom_col, 'short_cluster_id', start_col])

    # collect gaps inside each cluster (ignore the first member per cluster)
    def collect_gaps(s):
        vals = s.values.tolist()
        if len(vals) <= 1:
            return []
        # first is NaN (or gap to an earlier short outside cluster), keep from second onward
        return [v for v in vals[1:] if pd.notna(v)]

    per_cluster_gaps = x_sorted.groupby('short_cluster_id')['gap_to_prev_short_on_scaffold']                                .apply(collect_gaps).rename('gaps_list').reset_index()

    g = x.groupby('short_cluster_id').agg(
        scaffold=(chrom_col, 'first'),
        short_cluster_start=('short_cluster_start', 'min'),
        short_cluster_end=('short_cluster_end', 'max'),
        short_cluster_span=('short_cluster_span', 'max'),
        n_short_segments=('short_cluster_id', 'size'),
        total_short_bases=('seg_length', 'sum')
    ).reset_index()

    g = g.merge(per_cluster_gaps, on='short_cluster_id', how='left')
    def safe_mean(v): 
        v = [vv for vv in v if pd.notna(vv)]
        return float(np.mean(v)) if len(v) else np.nan
    def safe_max(v):
        v = [vv for vv in v if pd.notna(vv)]
        return float(np.max(v)) if len(v) else np.nan
    g['mean_gap_between_members'] = g['gaps_list'].apply(safe_mean)
    g['max_gap_between_members']  = g['gaps_list'].apply(safe_max)
    g = g.drop(columns=['gaps_list']).sort_values(['scaffold','short_cluster_start']).reset_index(drop=True)
    return g

def main():
    ap = argparse.ArgumentParser(description="Annotate clusters of nearby short TMRCA segments on the same scaffold.")
    ap.add_argument('--tmrca', required=True, help='Input TMRCA TSV')
    ap.add_argument('--outdir', required=True, help='Output directory')
    ap.add_argument('--chrom-col', default='CHROM')
    ap.add_argument('--start-col', default='start_tmrca')
    ap.add_argument('--end-col', default='end_tmrca')
    ap.add_argument('--short-max', type=int, default=100, help='Short segment threshold (bp)')
    ap.add_argument('--short-criterion', choices=['lt','le'], default='lt', help='Use length < short-max (lt) or <= (le)')
    ap.add_argument('--near-gap', type=int, default=200, help='Max gap (bp) between short segments to keep them in the same cluster')
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = load_tmrca(args.tmrca, args.chrom_col, args.start_col, args.end_col)

    ann = annotate_short_clusters(df, args.chrom_col, args.start_col, args.end_col,
                                  args.short_max, args.short_criterion, args.near_gap)
    ann.to_csv(outdir / 'tmrca_with_short_clusters.tsv', sep='\t', index=False)

    clusters = summarize_clusters(ann, args.chrom_col, args.start_col, args.end_col)
    clusters.to_csv(outdir / 'short_clusters.tsv', sep='\t', index=False)

    with open(outdir / 'summary.txt', 'w') as f:
        n_short = int(ann['is_short'].sum())
        n_clusters = int(clusters.shape[0])
        f.write(f"Short criterion: length {'<' if args.short_criterion=='lt' else '<='} {args.short_max} bp\n")
        f.write(f"Nearby gap threshold: <= {args.near_gap} bp on same scaffold\n")
        f.write(f"Total segments: {ann.shape[0]}\n")
        f.write(f"Short segments: {n_short}\n")
        f.write(f"Short clusters: {n_clusters}\n")

if __name__ == '__main__':
    main()
