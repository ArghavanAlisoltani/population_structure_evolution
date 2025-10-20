#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Length_Matched_Fisher_mRNA_V5.py
--------------------------------
Same logic as V4 (inclusive lengths, short-segment cap to 1 mRNA, regional deduplication),
plus robust length matching with bins so results do not go empty when exact lengths
are sparse (common at the region level).

Key additions
- --len-match-mode {exact,bin}    (default: bin)
- --len-bin-width INT             (default: 200)
- For mode=bin, lengths are grouped into bins and Fisher's test is run per bin.
  min-len-freq applies to bin counts.

Outputs
- tmrca_with_mrna.tsv, regions.tsv        (unchanged)
- fisher_overall_mRNA_segments.tsv        (unchanged)
- fisher_overall_mRNA_regions.tsv         (unchanged)
- If mode=exact:   fisher_by_length_mRNA_segments.tsv / fisher_by_length_mRNA_regions.tsv
- If mode=bin:     fisher_by_lenbin_mRNA_segments.tsv / fisher_by_lenbin_mRNA_regions.tsv
- length_frequencies_mRNA_segments.tsv and length_frequencies_mRNA_regions.tsv always reflect the chosen key (length or bin).

Run example
python Length_Matched_Fisher_mRNA_V5.py \
  --tmrca annotated_tmrca_4_GPT_13columns.tsv \
  --mrna /Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt \
  --top-percentile 5 \
  --min-len-freq 50 \
  --cluster-gap 400 \
  --short-min 2 --short-max 400 \
  --len-match-mode bin --len-bin-width 200 \
  --outdir fisher_mrna_out_v5
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from bisect import bisect_left
from scipy.stats import fisher_exact

def ok_int(x):
    try:
        return int(x)
    except Exception:
        try:
            f = float(x)
            if not np.isfinite(f) or f < 0:
                return None
            return int(round(f))
        except Exception:
            return None

def load_tmdf(path, chrom_col, start_col, end_col, mean_col):
    df = pd.read_csv(path, sep='\t', dtype={chrom_col: str}, low_memory=False)
    for c in [chrom_col, start_col, end_col, mean_col]:
        if c not in df.columns:
            raise ValueError(f"Required TMRCA column '{c}' not found in {path}")
    df[start_col] = df[start_col].apply(ok_int)
    df[end_col]   = df[end_col].apply(ok_int)
    df = df.dropna(subset=[start_col, end_col])
    df = df[df[end_col] >= df[start_col]]
    return df

def load_mrna(path, id_col, seq_col, start_col, end_col):
    df = pd.read_csv(path, sep='\t', dtype={seq_col: str}, low_memory=False)
    for c in [id_col, seq_col, start_col, end_col]:
        if c not in df.columns:
            raise ValueError(f"Required mRNA column '{c}' not found in {path}")
    df[start_col] = df[start_col].apply(ok_int)
    df[end_col]   = df[end_col].apply(ok_int)
    df = df.dropna(subset=[start_col, end_col])
    df = df[df[end_col] >= df[start_col]]
    return df[[id_col, seq_col, start_col, end_col]].rename(
        columns={id_col:"mrna_id", seq_col:"mrna_seqid", start_col:"mrna_start", end_col:"mrna_end"})

def build_chr_index(mrna):
    idx = {}
    for chrom, sub in mrna.groupby("mrna_seqid"):
        sub = sub.sort_values("mrna_start").reset_index(drop=True)
        idx[chrom] = {
            "starts": sub["mrna_start"].to_numpy(),
            "ends":   sub["mrna_end"].to_numpy(),
            "ids":    sub["mrna_id"].astype(str).to_numpy()
        }
    return idx

def find_overlap_indices(chr_idx, seg_start, seg_end):
    starts = chr_idx["starts"]; ends = chr_idx["ends"]
    n = len(starts)
    hits = []
    i = bisect_left(starts, seg_start)
    j = i-1
    while j >= 0 and ends[j] >= seg_start:
        if not (ends[j] < seg_start or starts[j] > seg_end):
            hits.append(j)
        j -= 1
    k = i
    while k < n and starts[k] <= seg_end:
        if not (ends[k] < seg_start or starts[k] > seg_end):
            hits.append(k)
        k += 1
    # unique preserve order
    seen=set(); out=[]
    for h in hits:
        if h not in seen:
            seen.add(h); out.append(h)
    return out

def choose_nearest_one(chr_idx, seg_start, seg_end, candidate_idxs):
    if not candidate_idxs:
        return []
    seg_mid = (seg_start + seg_end) / 2.0
    starts = chr_idx["starts"]; ends = chr_idx["ends"]; ids = chr_idx["ids"]
    dists = []
    for i in candidate_idxs:
        g_mid = (starts[i] + ends[i]) / 2.0
        dists.append((abs(g_mid - seg_mid), i))
    dists.sort(key=lambda x: (x[0], x[1]))
    best_idx = dists[0][1]
    return [ids[best_idx]]

def find_overlaps(chr_idx, seg_start, seg_end):
    idxs = find_overlap_indices(chr_idx, seg_start, seg_end)
    if not idxs:
        return []
    ids = chr_idx["ids"][idxs].tolist()
    seen=set(); out=[]
    for x in ids:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int); ranks[order] = np.arange(1, n+1)
    fdr = pvals * n / ranks
    fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
    out = np.empty_like(pvals, dtype=float); out[order] = np.clip(fdr_sorted, 0, 1)
    return out

def safe_fisher(a,b,c,d):
    try:
        odds, p = fisher_exact([[a,b],[c,d]], alternative='two-sided')
    except Exception:
        odds, p = (np.nan, 1.0)
    rr = ((a+0.5)/(a+b+1.0)) / ((c+0.5)/(c+d+1.0))
    return odds, rr, p

def cluster_segments(df, chrom_col, start_col, end_col, cluster_gap):
    df = df.sort_values([chrom_col, start_col, end_col]).reset_index(drop=True).copy()
    region_id = []
    current_region = -1
    prev_chrom = None
    prev_end = None
    for _, row in df.iterrows():
        ch = row[chrom_col]; st = row[start_col]; en = row[end_col]
        if prev_chrom != ch:
            current_region += 1
            region_id.append(current_region)
            prev_chrom = ch
            prev_end = en
            continue
        gap = st - prev_end - 1
        if gap <= cluster_gap:
            region_id.append(current_region)
            prev_end = max(prev_end, en)
        else:
            current_region += 1
            region_id.append(current_region)
            prev_end = en
    df["region_id"] = region_id
    return df

def summarize_regions(tmdf, chrom_col, start_col, end_col, mean_col):
    g = tmdf.groupby("region_id", sort=True)
    out = g.agg(
        region_chrom=(chrom_col, "first"),
        region_start=(start_col, "min"),
        region_end=(end_col, "max"),
        n_segments=(start_col, "size"),
        region_max_mean_tmrca=(mean_col, "max"),
        region_mean_mean_tmrca=(mean_col, "mean")
    ).reset_index()
    out["region_length"] = (out["region_end"] - out["region_start"] + 1).astype(int)
    return out

def add_segment_overlaps(tmdf, chr_index, chrom_col, start_col, end_col, short_min, short_max):
    raw_counts = []; raw_ids = []
    capped_counts = []; capped_ids = []
    seg_len = (tmdf[end_col] - tmdf[start_col] + 1).astype(int).to_numpy()
    for ch, st, en, L in zip(tmdf[chrom_col].astype(str).to_numpy(),
                              tmdf[start_col].to_numpy(),
                              tmdf[end_col].to_numpy(),
                              seg_len):
        if ch in chr_index:
            idxs = find_overlap_indices(chr_index[ch], int(st), int(en))
            if idxs:
                ids = chr_index[ch]["ids"][idxs].tolist()
                seen=set(); uniq_ids=[]
                for x in ids:
                    if x not in seen:
                        seen.add(x); uniq_ids.append(x)
                raw_counts.append(len(uniq_ids)); raw_ids.append(";".join(uniq_ids))
                if short_min <= L <= short_max:
                    chosen = choose_nearest_one(chr_index[ch], int(st), int(en), idxs)
                    capped_counts.append(len(chosen)); capped_ids.append(";".join(chosen))
                else:
                    capped_counts.append(len(uniq_ids)); capped_ids.append(";".join(uniq_ids))
            else:
                raw_counts.append(0); raw_ids.append("")
                capped_counts.append(0); capped_ids.append("")
        else:
            raw_counts.append(0); raw_ids.append("")
            capped_counts.append(0); capped_ids.append("")
    tmdf["overlap_mrna_n_raw"]  = raw_counts
    tmdf["overlap_mrna_ids_raw"] = raw_ids
    tmdf["overlap_mrna_n"]      = capped_counts
    tmdf["overlap_mrna_ids"]    = capped_ids
    return tmdf

def run_length_matched_fisher(table, length_key, is_top_col, hit_col, min_len_freq, out_path):
    lf = table[length_key].value_counts().sort_index().rename_axis(length_key).reset_index(name="count")
    lf.to_csv(out_path.parent / ("length_frequencies_" + out_path.name.replace("fisher_by_", "").replace(".tsv","") + ".tsv"),
              sep='\t', index=False)
    keep = lf.loc[lf["count"] >= min_len_freq, length_key].to_numpy()
    rows = []
    for L in keep:
        sub = table.loc[table[length_key] == L]
        hit = sub[hit_col] > 0
        top = sub[is_top_col].astype(bool)
        a = int(((top) & (hit)).sum()); b = int(((top) & (~hit)).sum())
        c = int(((~top) & (hit)).sum()); d = int(((~top) & (~hit)).sum())
        odds, rr, p = safe_fisher(a,b,c,d)
        rows.append({length_key:L, "N_total":int(len(sub)), "N_top":int(top.sum()), "N_bkg":int((~top).sum()),
                     "TOP_hit":a, "TOP_miss":b, "BKG_hit":c, "BKG_miss":d,
                     "TOP_rate":a/(a+b+1e-15), "BKG_rate":c/(c+d+1e-15),
                     "odds_ratio":odds, "risk_ratio":rr, "fisher_p":p})
    res = pd.DataFrame(rows)
    if not res.empty:
        res["fdr_bh"] = benjamini_hochberg(res["fisher_p"].to_numpy())
    res.to_csv(out_path, sep='\t', index=False)

def main():
    ap = argparse.ArgumentParser(description="mRNA overlaps with TMRCA; short-segment cap; regional dedup; length-matched Fisher with bins.")
    ap.add_argument("--tmrca", required=True)
    ap.add_argument("--mrna", required=True)
    ap.add_argument("--top-percentile", type=float, default=5.0)
    ap.add_argument("--min-len-freq", type=int, default=50, help="Minimum number of items per length/len-bin to test")
    ap.add_argument("--cluster-gap", type=int, default=400)
    ap.add_argument("--short-min", type=int, default=2)
    ap.add_argument("--short-max", type=int, default=400)
    ap.add_argument("--chrom-col", default="CHROM")
    ap.add_argument("--start-col", default="start_tmrca")
    ap.add_argument("--end-col", default="end_tmrca")
    ap.add_argument("--mean-col", default="mean_tmrca")
    ap.add_argument("--mrna-id-col", default="mrna_id")
    ap.add_argument("--mrna-seq-col", default="new_seqid")
    ap.add_argument("--mrna-start-col", default="new_start")
    ap.add_argument("--mrna-end-col", default="new_end")
    ap.add_argument("--len-match-mode", choices=["exact","bin"], default="bin")
    ap.add_argument("--len-bin-width", type=int, default=200)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    tmdf = load_tmdf(args.tmrca, args.chrom_col, args.start_col, args.end_col, args.mean_col)
    mdf  = load_mrna(args.mrna, args.mrna_id_col, args.mrna_seq_col, args.mrna_start_col, args.mrna_end_col)

    tmdf["seg_length"] = (tmdf[args.end_col] - tmdf[args.start_col] + 1).astype(int)
    cutoff = np.percentile(tmdf[args.mean_col].to_numpy(), 100 - args.top_percentile)
    tmdf["is_top"] = tmdf[args.mean_col] >= cutoff

    chr_index = build_chr_index(mdf)
    tmdf = add_segment_overlaps(tmdf, chr_index, args.chrom_col, args.start_col, args.end_col, args.short_min, args.short_max)

    tmdf = cluster_segments(tmdf, args.chrom_col, args.start_col, args.end_col, args.cluster_gap)
    regions = summarize_regions(tmdf, args.chrom_col, args.start_col, args.end_col, args.mean_col)
    regions = regions.merge(tmdf.groupby("region_id")["is_top"].any().rename("region_is_top").reset_index(),
                            on="region_id", how="left")

    reg_n = []; reg_ids = []
    for _, r in regions.iterrows():
        ch = r["region_chrom"]; st = int(r["region_start"]); en = int(r["region_end"])
        if ch in chr_index:
            ids = find_overlaps(chr_index[ch], st, en)
        else:
            ids = []
        reg_n.append(len(ids)); reg_ids.append(";".join(ids) if ids else "")
    regions["region_overlap_mrna_n"]   = reg_n
    regions["region_overlap_mrna_ids"] = reg_ids

    tmdf = tmdf.merge(regions[["region_id","region_overlap_mrna_n","region_overlap_mrna_ids","region_length"]],
                      on="region_id", how="left")

    tmdf_out = tmdf[[args.chrom_col, args.start_col, args.end_col, args.mean_col, "seg_length", "is_top",
                     "overlap_mrna_n_raw","overlap_mrna_ids_raw",
                     "overlap_mrna_n","overlap_mrna_ids",
                     "region_id","region_overlap_mrna_n","region_overlap_mrna_ids","region_length"]]
    tmdf_out.to_csv(outdir/"tmrca_with_mrna.tsv", sep='\t', index=False)
    regions.to_csv(outdir/"regions.tsv", sep='\t', index=False)

    # pooled Fishers
    hit_seg = (tmdf_out["overlap_mrna_n"] > 0); top_seg = tmdf_out["is_top"].astype(bool)
    a = int(((top_seg) & (hit_seg)).sum()); b = int(((top_seg) & (~hit_seg)).sum())
    c = int(((~top_seg) & (hit_seg)).sum()); d = int(((~top_seg) & (~hit_seg)).sum())
    odds, rr, p = safe_fisher(a,b,c,d)
    pd.DataFrame([{"level":"segment","category":"Any_mRNA","distance_bp":0,
                   "TOP_hit":a,"TOP_miss":b,"BKG_hit":c,"BKG_miss":d,
                   "TOP_rate":a/(a+b+1e-15), "BKG_rate":c/(c+d+1e-15),
                   "odds_ratio":odds, "risk_ratio":rr, "fisher_p":p, "fdr_bh":p}])        .to_csv(outdir/"fisher_overall_mRNA_segments.tsv", sep='\t', index=False)

    hit_reg = (regions["region_overlap_mrna_n"] > 0); top_reg = regions["region_is_top"].astype(bool)
    a = int(((top_reg) & (hit_reg)).sum()); b = int(((top_reg) & (~hit_reg)).sum())
    c = int(((~top_reg) & (hit_reg)).sum()); d = int(((~top_reg) & (~hit_reg)).sum())
    odds, rr, p = safe_fisher(a,b,c,d)
    pd.DataFrame([{"level":"region","category":"Any_mRNA","distance_bp":0,
                   "TOP_hit":a,"TOP_miss":b,"BKG_hit":c,"BKG_miss":d,
                   "TOP_rate":a/(a+b+1e-15), "BKG_rate":c/(c+d+1e-15),
                   "odds_ratio":odds, "risk_ratio":rr, "fisher_p":p, "fdr_bh":p}])        .to_csv(outdir/"fisher_overall_mRNA_regions.tsv", sep='\t', index=False)

    # length matching mode
    if args.len_match_mode == "exact":
        seg_key = "seg_length"; reg_key = "region_length"
        seg_out = outdir/"fisher_by_length_mRNA_segments.tsv"
        reg_out = outdir/"fisher_by_length_mRNA_regions.tsv"
        tmdf_out["length_key"] = tmdf_out[seg_key]
        regions["length_key"]  = regions[reg_key]
    else:
        def to_bin(x, w):
            start = (int(x)-1)//w*w + 1
            end   = start + w - 1
            return f"[{start}-{end}]"
        tmdf_out["seg_len_bin"]   = tmdf_out["seg_length"].apply(lambda x: to_bin(x, args.len_bin_width))
        regions["region_len_bin"] = regions["region_length"].apply(lambda x: to_bin(x, args.len_bin_width))
        seg_key = "seg_len_bin"; reg_key = "region_len_bin"
        seg_out = outdir/"fisher_by_lenbin_mRNA_segments.tsv"
        reg_out = outdir/"fisher_by_lenbin_mRNA_regions.tsv"
        tmdf_out["length_key"] = tmdf_out[seg_key]
        regions["length_key"]  = regions[reg_key]

    run_length_matched_fisher(tmdf_out, seg_key, "is_top", "overlap_mrna_n", args.min_len_freq, seg_out)
    run_length_matched_fisher(regions, reg_key, "region_is_top", "region_overlap_mrna_n", args.min_len_freq, reg_out)

    # length frequencies written above per out file; also emit generic ones for convenience
    tmdf_out["length_key"].value_counts().sort_index().rename_axis("length_key").reset_index(name="count")        .to_csv(outdir/"length_frequencies_mRNA_segments.tsv", sep='\t', index=False)
    regions["length_key"].value_counts().sort_index().rename_axis("length_key").reset_index(name="count")        .to_csv(outdir/"length_frequencies_mRNA_regions.tsv", sep='\t', index=False)

    with open(outdir/"summary_mRNA.txt", "w") as f:
        f.write("mRNA overlaps (0 bp) with inclusive lengths, short-segment cap, regional deduplication\n")
        f.write(f"Top percentile = {args.top_percentile:.1f}% (cutoff on {args.mean_col})\n")
        f.write(f"Cluster gap (bp) = {args.cluster_gap}\n")
        f.write(f"Short-segment cap range (bp): [{args.short_min}, {args.short_max}]\n")
        f.write(f"Length-match mode = {args.len_match_mode}; bin width = {args.len_bin_width}\n")
        f.write(f"Segments total: {tmdf.shape[0]:,}; Regions total: {regions.shape[0]:,}\n")

if __name__ == "__main__":
    try:
        import scipy  # noqa: F401
    except Exception:
        print("ERROR: This script requires scipy (for Fisher's exact test). Install via: pip install scipy", file=sys.stderr)
        raise SystemExit(1)
    main()
