#!/usr/bin/env python3
"""
TMRCA_bins_genes_TEs_SNPs_V1.py

- Compute average TMRCA per gene (mRNA) from overlapping TMRCA segments
- Compute average TMRCA per TE element (any TE), and subsets for Copia & Gypsy
- Annotate SNPs with TMRCA at their genomic position (if needed)
- Define TMRCA bins (automatic log10 or user-defined edges)
- For each bin, compute:
    * number of TMRCA segments
    * number of SNPs
    * number of genes (mRNAs)
    * number of TE elements (all)
    * number of Copia elements
    * number of Gypsy elements
- Write all per-bin frequencies into one TSV file

Assumes TMRCA segments form a non-overlapping, ordered partition of the genome.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------

def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)


def pick_col(cols, candidates, required=True, label="column"):
    """Pick first matching column from candidates (case-insensitive)."""
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    if required:
        raise SystemExit(f"ERROR: Could not find {label}. Looked for: {candidates}")
    return None


def merge_asof_per_chrom_points(pts_df, seg_df, pt_scaf, pt_pos,
                                seg_scaf, seg_start, seg_end, seg_mean):
    """
    Annotate point positions (SNPs) with TMRCA segments via a per-scaffold merge_asof.

    For each scaffold:
      - merge pts on position vs segment start ("backward")
      - then keep only matches where pos in [start, end]
      - if not inside, set mean_tmrca_at_pos = NaN
    """
    out = []
    for scaf, pts_sub in pts_df.groupby(pt_scaf, sort=False):
        seg_sub = seg_df[seg_df[seg_scaf] == scaf]
        pts_sub = pts_sub.copy()
        if seg_sub.empty:
            pts_sub["mean_tmrca_at_pos"] = np.nan
            out.append(pts_sub)
            continue

        pts_sub[pt_pos] = pd.to_numeric(pts_sub[pt_pos], errors="coerce").astype("float64")
        seg_sub = seg_sub.copy()
        seg_sub[seg_start] = pd.to_numeric(seg_sub[seg_start], errors="coerce").astype("float64")
        seg_sub[seg_end]   = pd.to_numeric(seg_sub[seg_end], errors="coerce").astype("float64")

        pts_sub.sort_values(pt_pos, inplace=True)
        seg_sub.sort_values(seg_start, inplace=True)

        # keep only needed columns in segments
        seg_keep = seg_sub[[seg_start, seg_end, seg_mean]].copy()

        merged = pd.merge_asof(
            pts_sub,
            seg_keep,
            left_on=pt_pos, right_on=seg_start,
            direction="backward",
            allow_exact_matches=True
        )

        inside = merged[pt_pos].between(merged[seg_start],
                                        merged[seg_end],
                                        inclusive="both")
        merged.loc[~inside, seg_mean] = np.nan
        merged.rename(columns={seg_mean: "mean_tmrca_at_pos"}, inplace=True)

        # drop helper columns seg_start, seg_end to avoid confusion
        merged.drop(columns=[seg_start, seg_end], inplace=True)

        out.append(merged)

    return pd.concat(out, ignore_index=True) if out else pts_df.copy()


def build_segments_by_scaffold(seg_df, scaf_col, start_col, end_col, mean_col):
    """
    Build a dict: scaffold -> (starts, ends, means) numpy arrays, sorted by start.
    """
    seg_dict = {}
    for scaf, sub in seg_df.groupby(scaf_col, sort=False):
        sub = sub.copy()
        sub[start_col] = pd.to_numeric(sub[start_col], errors="coerce")
        sub[end_col]   = pd.to_numeric(sub[end_col], errors="coerce")
        sub[mean_col]  = pd.to_numeric(sub[mean_col], errors="coerce")
        sub = sub.sort_values(start_col)
        starts = sub[start_col].to_numpy()
        ends   = sub[end_col].to_numpy()
        means  = sub[mean_col].to_numpy()
        seg_dict[str(scaf)] = (starts, ends, means)
    return seg_dict


def assign_interval_tmrca(interval_df,
                          iv_scaf, iv_start, iv_end,
                          seg_dict,
                          out_mean_col="tmrca_avg",
                          out_nseg_col="n_segments_overlapped"):
    """
    For each interval in interval_df, find all TMRCA segments on same scaffold
    whose [start,end] overlaps [iv_start,iv_end], and assign the mean of their
    mean_tmrca to out_mean_col. Also record how many segments were used.
    """
    tm_list = []
    nseg_list = []

    for _, row in interval_df.iterrows():
        scaf = str(row[iv_scaf])
        gstart = float(row[iv_start])
        gend   = float(row[iv_end])

        if scaf not in seg_dict:
            tm_list.append(np.nan)
            nseg_list.append(0)
            continue

        seg_starts, seg_ends, seg_means = seg_dict[scaf]

        # segments assumed non-overlapping & ordered
        left = np.searchsorted(seg_ends, gstart, side="left")
        right = np.searchsorted(seg_starts, gend, side="right") - 1

        if left <= right and right >= 0 and left < len(seg_means):
            tm_val = float(np.nanmean(seg_means[left:right+1]))
            nseg   = int(right - left + 1)
        else:
            tm_val = np.nan
            nseg   = 0

        tm_list.append(tm_val)
        nseg_list.append(nseg)

    interval_df[out_mean_col] = tm_list
    interval_df[out_nseg_col] = nseg_list
    return interval_df


def build_tmrca_bins(tmrca_values, bins_str=None, nbins=8):
    """
    Build TMRCA bin edges (linear scale, but can be log10 spaced).

    - If bins_str is provided (comma-separated numbers), use those as edges
      (and extend to cover full data if needed).
    - Otherwise, automatically define ~log10-spaced bins across min..max.
    """
    vals = pd.to_numeric(tmrca_values, errors="coerce").dropna()
    vals = vals[vals > 0]
    if vals.empty:
        raise SystemExit("ERROR: No positive TMRCA values to define bins.")

    vals_min = float(vals.min())
    vals_max = float(vals.max())

    if bins_str:
        edges = [float(x.strip()) for x in str(bins_str).split(",") if x.strip() != ""]
        edges = sorted(edges)

        if edges[0] > vals_min:
            edges = [vals_min] + edges
        if edges[-1] < vals_max:
            edges = edges + [vals_max]

        return np.array(edges, dtype=float)

    # automatic log10 bins
    log_min = float(np.floor(np.log10(vals_min)))
    log_max = float(np.ceil(np.log10(vals_max)))
    nbins = max(1, int(nbins))
    edges = np.logspace(log_min, log_max, nbins + 1)
    # make sure lower edge not above min
    edges[0] = min(edges[0], vals_min)
    return edges


def count_in_bins(values, bins_edges):
    """
    Count how many finite values fall into each bin (using left-closed, right-open).
    """
    series = pd.to_numeric(values, errors="coerce")
    cats = pd.cut(series, bins=bins_edges, right=False, include_lowest=True)
    bin_order = cats.cat.categories
    counts = cats.value_counts().reindex(bin_order, fill_value=0)
    return counts, bin_order


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=(
            "Compute average TMRCA per gene and TE, then bin TMRCA values and "
            "compute frequencies for segments, SNPs, genes, TEs, Copia, Gypsy."
        )
    )

    # core inputs
    p.add_argument("--tmrca", required=True,
                   help="TMRCA segments TSV (CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    p.add_argument("--mrna", required=True,
                   help="mRNA annotation TSV (must have scaffold, start, end, gene ID).")
    p.add_argument("--te", required=True,
                   help="TE annotation TSV (must have scaffold, start, end, sequence_ontology or similar).")
    p.add_argument("--snps", required=True,
                   help="SNP positions TSV (must have scaffold & position).")

    # column names: TMRCA segments
    p.add_argument("--tmrca-scaf-col", default="CHROM")
    p.add_argument("--tmrca-start-col", default="start_tmrca")
    p.add_argument("--tmrca-end-col", default="end_tmrca")
    p.add_argument("--tmrca-mean-col", default="mean_tmrca")

    # column names: mRNA
    p.add_argument("--mrna-scaf-col", default="new_seqid")
    p.add_argument("--mrna-start-col", default="new_start")
    p.add_argument("--mrna-end-col", default="new_end")
    p.add_argument("--mrna-id-col", default="mrna_id")

    # column names: TE
    p.add_argument("--te-scaf-col", default="seqid")
    p.add_argument("--te-start-col", default="start")
    p.add_argument("--te-end-col", default="end")
    p.add_argument("--te-class-col", default="sequence_ontology")

    # column names: SNPs
    p.add_argument("--snp-scaf-col", default="CHROM")
    p.add_argument("--snp-pos-col", default="POS")
    p.add_argument("--snp-tmrca-col", default=None,
                   help="If provided and exists, use this TMRCA column for SNPs instead of re-annotating.")

    # binning options
    p.add_argument("--bins", default=None,
                   help="Optional comma-separated bin edges (e.g., '0,50,500,5000,50000,200000'). "
                        "If omitted, bins are built automatically (log10).")
    p.add_argument("--nbins", type=int, default=8,
                   help="Number of automatic bins if --bins is not given (default: 8).")

    # output
    p.add_argument("--outdir", default="tmrca_bins_out_v1",
                   help="Output directory.")

    args = p.parse_args()
    ensure_dir(args.outdir)

    # -----------------------------------------------------------------
    # 1. Load TMRCA segments
    # -----------------------------------------------------------------
    tmrca = pd.read_csv(args.tmrca, sep="\t", dtype={})
    tcols = tmrca.columns.tolist()

    t_scaf = pick_col(tcols, [args.tmrca_scaf_col], label="TMRCA scaffold")
    t_start = pick_col(tcols, [args.tmrca_start_col], label="TMRCA start")
    t_end   = pick_col(tcols, [args.tmrca_end_col], label="TMRCA end")
    t_mean  = pick_col(tcols, [args.tmrca_mean_col], label="TMRCA mean")

    tmrca[t_scaf] = tmrca[t_scaf].astype(str)
    tmrca[t_start] = pd.to_numeric(tmrca[t_start], errors="coerce")
    tmrca[t_end]   = pd.to_numeric(tmrca[t_end], errors="coerce")
    tmrca[t_mean]  = pd.to_numeric(tmrca[t_mean], errors="coerce")

    # -----------------------------------------------------------------
    # 2. Build segment dict for interval-based TMRCA assignment
    # -----------------------------------------------------------------
    seg_dict = build_segments_by_scaffold(tmrca, t_scaf, t_start, t_end, t_mean)

    # -----------------------------------------------------------------
    # 3. Load mRNA annotation and assign average TMRCA per gene
    # -----------------------------------------------------------------
    mrna = pd.read_csv(args.mrna, sep="\t", dtype={})
    mcols = mrna.columns.tolist()

    m_scaf = pick_col(mcols, [args.mrna_scaf_col], label="mRNA scaffold")
    m_start = pick_col(mcols, [args.mrna_start_col], label="mRNA start")
    m_end   = pick_col(mcols, [args.mrna_end_col], label="mRNA end")
    m_id    = pick_col(mcols, [args.mrna_id_col], label="mRNA ID")

    mrna[m_scaf] = mrna[m_scaf].astype(str)
    mrna[m_start] = pd.to_numeric(mrna[m_start], errors="coerce")
    mrna[m_end]   = pd.to_numeric(mrna[m_end], errors="coerce")

    mrna_tm = assign_interval_tmrca(
        mrna.copy(),
        iv_scaf=m_scaf, iv_start=m_start, iv_end=m_end,
        seg_dict=seg_dict,
        out_mean_col="tmrca_avg",
        out_nseg_col="n_segments_overlapped"
    )

    gene_out = mrna_tm[[m_id, m_scaf, m_start, m_end, "tmrca_avg", "n_segments_overlapped"]].copy()
    gene_out.to_csv(Path(args.outdir, "gene_tmrca.tsv"), sep="\t", index=False)

    # -----------------------------------------------------------------
    # 4. Load TE annotation and assign average TMRCA per TE element
    # -----------------------------------------------------------------
    te = pd.read_csv(args.te, sep="\t", dtype={})
    tecols = te.columns.tolist()

    te_scaf = pick_col(tecols, [args.te_scaf_col], label="TE scaffold")
    te_start = pick_col(tecols, [args.te_start_col], label="TE start")
    te_end   = pick_col(tecols, [args.te_end_col], label="TE end")
    te_class = pick_col(tecols, [args.te_class_col], label="TE class / sequence_ontology")

    te[te_scaf] = te[te_scaf].astype(str)
    te[te_start] = pd.to_numeric(te[te_start], errors="coerce")
    te[te_end]   = pd.to_numeric(te[te_end], errors="coerce")

    # give each TE an ID row index if not present
    te["te_id"] = np.arange(1, len(te) + 1)

    te_tm = assign_interval_tmrca(
        te.copy(),
        iv_scaf=te_scaf, iv_start=te_start, iv_end=te_end,
        seg_dict=seg_dict,
        out_mean_col="tmrca_avg",
        out_nseg_col="n_segments_overlapped"
    )

    te_out_all = te_tm[["te_id", te_scaf, te_start, te_end, te_class, "tmrca_avg", "n_segments_overlapped"]].copy()
    te_out_all.to_csv(Path(args.outdir, "te_tmrca_all.tsv"), sep="\t", index=False)

    # subsets: Copia & Gypsy based on substring in te_class
    te_tm["te_class_string"] = te_tm[te_class].astype(str)
    te_copia = te_tm[te_tm["te_class_string"].str.contains("Copia", case=False, na=False)].copy()
    te_gypsy = te_tm[te_tm["te_class_string"].str.contains("Gypsy", case=False, na=False)].copy()

    te_copia_out = te_copia[["te_id", te_scaf, te_start, te_end, te_class, "tmrca_avg", "n_segments_overlapped"]].copy()
    te_gypsy_out = te_gypsy[["te_id", te_scaf, te_start, te_end, te_class, "tmrca_avg", "n_segments_overlapped"]].copy()

    te_copia_out.to_csv(Path(args.outdir, "te_tmrca_copia.tsv"), sep="\t", index=False)
    te_gypsy_out.to_csv(Path(args.outdir, "te_tmrca_gypsy.tsv"), sep="\t", index=False)

    # -----------------------------------------------------------------
    # 5. Load SNP positions and annotate with TMRCA if needed
    # -----------------------------------------------------------------
    snps = pd.read_csv(args.snps, sep="\t", dtype={})
    scols = snps.columns.tolist()

    s_scaf = pick_col(scols, [args.snp_scaf_col], label="SNP scaffold")
    s_pos  = pick_col(scols, [args.snp_pos_col], label="SNP position")

    snps[s_scaf] = snps[s_scaf].astype(str)
    snps[s_pos]  = pd.to_numeric(snps[s_pos], errors="coerce")

    # if user has given a column name for TMRCA and it exists, use it
    if args.snp_tmrca_col and args.snp_tmrca_col in snps.columns:
        snps["snp_tmrca"] = pd.to_numeric(snps[args.snp_tmrca_col], errors="coerce")
    else:
        # annotate from segments
        snps_anno = merge_asof_per_chrom_points(
            pts_df=snps,
            seg_df=tmrca[[t_scaf, t_start, t_end, t_mean]].copy(),
            pt_scaf=s_scaf, pt_pos=s_pos,
            seg_scaf=t_scaf, seg_start=t_start, seg_end=t_end, seg_mean=t_mean
        )
        snps = snps_anno
        snps["snp_tmrca"] = pd.to_numeric(snps["mean_tmrca_at_pos"], errors="coerce")

    snps_out = snps.copy()
    snps_out.to_csv(Path(args.outdir, "snps_with_tmrca.tsv"), sep="\t", index=False)

    # -----------------------------------------------------------------
    # 6. Define TMRCA bins based on original segment TMRCA values
    # -----------------------------------------------------------------
    bin_edges = build_tmrca_bins(tmrca[t_mean], bins_str=args.bins, nbins=args.nbins)

    # -----------------------------------------------------------------
    # 7. Count frequencies per bin for segments, SNPs, genes, TEs, Copia, Gypsy
    # -----------------------------------------------------------------
    seg_counts, bin_cats = count_in_bins(tmrca[t_mean], bin_edges)

    snp_counts, _ = count_in_bins(snps["snp_tmrca"], bin_edges)
    gene_counts, _ = count_in_bins(gene_out["tmrca_avg"], bin_edges)
    te_all_counts, _ = count_in_bins(te_out_all["tmrca_avg"], bin_edges)
    te_copia_counts, _ = count_in_bins(te_copia_out["tmrca_avg"], bin_edges)
    te_gypsy_counts, _ = count_in_bins(te_gypsy_out["tmrca_avg"], bin_edges)

    total_seg = int(seg_counts.sum())
    total_snps = int(snp_counts.sum())
    total_genes = int(gene_counts.sum())
    total_te_all = int(te_all_counts.sum())
    total_te_copia = int(te_copia_counts.sum())
    total_te_gypsy = int(te_gypsy_counts.sum())

    rows = []
    for i, cat in enumerate(bin_cats):
        left = float(cat.left)
        right = float(cat.right)
        label = f"[{left:.3g},{right:.3g})"

        n_seg = int(seg_counts.iloc[i])
        n_snp = int(snp_counts.iloc[i])
        n_gene = int(gene_counts.iloc[i])
        n_te_all = int(te_all_counts.iloc[i])
        n_te_copia = int(te_copia_counts.iloc[i])
        n_te_gypsy = int(te_gypsy_counts.iloc[i])

        rows.append({
            "bin_index": i,
            "bin_left": left,
            "bin_right": right,
            "bin_label": label,
            "n_segments": n_seg,
            "n_snps": n_snp,
            "n_genes": n_gene,
            "n_TEs_all": n_te_all,
            "n_TEs_Copia": n_te_copia,
            "n_TEs_Gypsy": n_te_gypsy,
            "frac_segments": n_seg / total_seg if total_seg > 0 else np.nan,
            "frac_snps": n_snp / total_snps if total_snps > 0 else np.nan,
            "frac_genes": n_gene / total_genes if total_genes > 0 else np.nan,
            "frac_TEs_all": n_te_all / total_te_all if total_te_all > 0 else np.nan,
            "frac_TEs_Copia": n_te_copia / total_te_copia if total_te_copia > 0 else np.nan,
            "frac_TEs_Gypsy": n_te_gypsy / total_te_gypsy if total_te_gypsy > 0 else np.nan,
        })

    bins_df = pd.DataFrame(rows)
    bins_df.to_csv(Path(args.outdir, "tmrca_bins_frequencies.tsv"), sep="\t", index=False)

    print("[OK] Wrote:")
    print(f"  gene TMRCA per mRNA:     {Path(args.outdir, 'gene_tmrca.tsv')}")
    print(f"  TE TMRCA (all):          {Path(args.outdir, 'te_tmrca_all.tsv')}")
    print(f"  TE TMRCA (Copia only):   {Path(args.outdir, 'te_tmrca_copia.tsv')}")
    print(f"  TE TMRCA (Gypsy only):   {Path(args.outdir, 'te_tmrca_gypsy.tsv')}")
    print(f"  SNPs with TMRCA:         {Path(args.outdir, 'snps_with_tmrca.tsv')}")
    print(f"  Bin frequencies:         {Path(args.outdir, 'tmrca_bins_frequencies.tsv')}")

if __name__ == "__main__":
    main()

