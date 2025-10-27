#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate_SNPs_with_TMRCA_TE_mRNA_V1.py
--------------------------------------
Annotate SNPs with TMRCA, TE, and mRNA.
See header for details and CLI options.

--- 
Run Examples

python Annotate_SNPs_with_TMRCA_TE_mRNA_V1.py \
  --tmrca annotated_tmrca_4_GPT_13columns.tsv \
  --snps  "/Users/aria/Desktop/OSU_projects/conifers/LP/vcf_v1/positions_split_poly_s100_scaffolds.tsv" \
  --te    "/Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv" \
  --mrna  "/Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt" \
  --top-percentile 5 \
  --out   snps_annotated_tmrca_te_mrna.tsv


  python Annotate_SNPs_with_TMRCA_TE_mRNA_V1.py \
  --tmrca  "~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025/all_tmrca_corrected_position.tsv" \
  --snps   "/Users/aria/Desktop/OSU_projects/conifers/LP/vcf_v1/positions_split_poly_s100_scaffolds.tsv" \
  --te     "/Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv" \
  --mrna   "/Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt" \
  --te-chrom seqid --te-start start --te-end end --te-class sequence_ontology \
  --top-percentile 5 \
  --out   "~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025/snp_annotations/snps_annotated_tmrca_te_mrna.tsv"

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

def point_query_all(entry, pos):
    starts = entry['starts']; ends = entry['ends']
    i = bisect.bisect_left(starts, pos)
    hits = []
    k = i - 1
    while k >= 0 and ends[k] >= pos:
        if starts[k] <= pos <= ends[k]:
            hits.append(k)
        k -= 1
    k = i
    n = len(starts)
    while k < n and starts[k] <= pos:
        if ends[k] >= pos:
            hits.append(k)
        k += 1
    return sorted(set(hits))

def main():
    ap = argparse.ArgumentParser(description="Annotate SNPs with TMRCA, TE, and mRNA.")
    ap.add_argument('--tmrca', required=True)
    ap.add_argument('--snps', required=True)
    ap.add_argument('--te', default=None)
    ap.add_argument('--mrna', default=None)
    ap.add_argument('--out', required=True)
    ap.add_argument('--top-percentile', type=float, default=5.0)

    ap.add_argument('--tmrca-chrom', default='CHROM')
    ap.add_argument('--tmrca-start', default='start_tmrca')
    ap.add_argument('--tmrca-end',   default='end_tmrca')
    ap.add_argument('--tmrca-mean',  default='mean_tmrca')

    ap.add_argument('--snp-chrom', default='CHROM')
    ap.add_argument('--snp-pos',   default='POS')

    ap.add_argument('--te-chrom', default='CHROM')
    ap.add_argument('--te-start', default='start')
    ap.add_argument('--te-end',   default='end')
    ap.add_argument('--te-class', default='TE_class')

    ap.add_argument('--mrna-id',   default='mrna_id')
    ap.add_argument('--mrna-seq',  default='new_seqid')
    ap.add_argument('--mrna-start',default='new_start')
    ap.add_argument('--mrna-end',  default='new_end')

    args = ap.parse_args()

    tdf = load_table(args.tmrca)
    require_cols(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], 'TMRCA')
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    tdf[args.tmrca_start] = pd.to_numeric(tdf[args.tmrca_start], errors='coerce').astype('Int64')
    tdf[args.tmrca_end]   = pd.to_numeric(tdf[args.tmrca_end], errors='coerce').astype('Int64')
    tdf = tdf.dropna(subset=[args.tmrca_start, args.tmrca_end]).copy()
    tdf[args.tmrca_start] = tdf[args.tmrca_start].astype(int)
    tdf[args.tmrca_end]   = tdf[args.tmrca_end].astype(int)

    sdf = load_table(args.snps)
    require_cols(sdf, [args.snp_chrom, args.snp_pos], 'SNPs')
    sdf[args.snp_chrom] = sdf[args.snp_chrom].astype(str)
    sdf[args.snp_pos]   = pd.to_numeric(sdf[args.snp_pos], errors='coerce').astype('Int64')
    sdf = sdf.dropna(subset=[args.snp_pos]).copy()
    sdf[args.snp_pos]   = sdf[args.snp_pos].astype(int)

    te_idx = None
    if args.te:
        tedf = load_table(args.te)
        require_cols(tedf, [args.te_chrom, args.te_start, args.te_end], 'TE')
        if args.te_class not in tedf.columns:
            tedf[args.te_class] = ""
        tedf[args.te_chrom] = tedf[args.te_chrom].astype(str)
        tedf[args.te_start] = pd.to_numeric(tedf[args.te_start], errors='coerce').astype('Int64')
        tedf[args.te_end]   = pd.to_numeric(tedf[args.te_end], errors='coerce').astype('Int64')
        tedf = tedf.dropna(subset=[args.te_start, args.te_end]).copy()
        tedf[args.te_start] = tedf[args.te_start].astype(int)
        tedf[args.te_end]   = tedf[args.te_end].astype(int)
        te_idx = build_interval_index(tedf, args.te_chrom, args.te_start, args.te_end, [args.te_class])

    mrna_idx = None
    if args.mrna:
        mdf = load_table(args.mrna)
        require_cols(mdf, [args.mrna_id, args.mrna_seq, args.mrna_start, args.mrna_end], 'mRNA')
        mdf[args.mrna_seq]   = mdf[args.mrna_seq].astype(str)
        mdf[args.mrna_start] = pd.to_numeric(mdf[args.mrna_start], errors='coerce').astype('Int64')
        mdf[args.mrna_end]   = pd.to_numeric(mdf[args.mrna_end], errors='coerce').astype('Int64')
        mdf = mdf.dropna(subset=[args.mrna_start, args.mrna_end]).copy()
        mdf[args.mrna_start] = mdf[args.mrna_start].astype(int)
        mdf[args.mrna_end]   = mdf[args.mrna_end].astype(int)
        mrna_idx = build_interval_index(mdf, args.mrna_seq, args.mrna_start, args.mrna_end, [args.mrna_id])

    tmrca_idx = build_interval_index(tdf, args.tmrca_chrom, args.tmrca_start, args.tmrca_end, [args.tmrca_mean])

    cutoff = np.percentile(tdf[args.tmrca_mean].to_numpy(), 100 - args.top_percentile)

    snp_tmrca = []; snp_group = []; snp_in_te = []; snp_te_class = []; snp_mrna_ids = []

    for ch, pos in zip(sdf[args.snp_chrom].astype(str).to_numpy(), sdf[args.snp_pos].to_numpy()):
        chs = str(ch)

        mt = np.nan; grp = "NA"
        entry = tmrca_idx.get(chs)
        if entry is not None:
            i = point_query_first(entry, int(pos))
            if i >= 0:
                mt = float(entry['payload'][args.tmrca_mean][i])
                grp = "top" if mt >= cutoff else "background"
        snp_tmrca.append(mt); snp_group.append(grp)

        in_te = False; te_class_val = ""
        if te_idx is not None and chs in te_idx:
            hits = point_query_all(te_idx[chs], int(pos))
            if hits:
                in_te = True
                classes = [te_idx[chs]['payload'][args.te_class][h] for h in hits]
                classes = [str(x) for x in classes if pd.notna(x)]
                seen=set(); uniq=[]
                for c in classes:
                    if c not in seen:
                        seen.add(c); uniq.append(c)
                te_class_val = ";".join(uniq)
        snp_in_te.append(in_te); snp_te_class.append(te_class_val)

        mrna_ids_val = ""
        if mrna_idx is not None and chs in mrna_idx:
            hits = point_query_all(mrna_idx[chs], int(pos))
            if hits:
                ids = [mrna_idx[chs]['payload'][args.mrna_id][h] for h in hits]
                ids = [str(x) for x in ids if pd.notna(x)]
                seen=set(); uniq=[]
                for c in ids:
                    if c not in seen:
                        seen.add(c); uniq.append(c)
                mrna_ids_val = ";".join(uniq)
        snp_mrna_ids.append(mrna_ids_val)

    out = sdf.copy()
    out['mean_tmrca_at_pos'] = snp_tmrca
    out['tmrca_group'] = snp_group
    out['in_TE'] = snp_in_te
    out['TE_class'] = snp_te_class
    out['mrna_id'] = snp_mrna_ids

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)

    n = out.shape[0]
    n_top = int((out['tmrca_group'] == 'top').sum())
    n_te  = int(out['in_TE'].sum()) if 'in_TE' in out.columns else 0
    with open(Path(args.out).with_suffix('.summary.txt'), 'w') as f:
        f.write(f"SNPs total: {n}\n")
        f.write(f"Top {args.top_percentile:.1f}% TMRCA SNPs: {n_top} ({(n_top/n*100.0 if n else 0):.2f}%)\n")
        if args.te:
            f.write(f"SNPs in TE: {n_te} ({(n_te/n*100.0 if n else 0):.2f}%)\n")

if __name__ == '__main__':
    main()
