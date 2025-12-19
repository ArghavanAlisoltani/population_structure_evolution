#!/usr/bin/env python3
"""
TMRCA_bins_genes_TEs_SNPs_V4.py

- Reads TMRCA segments, SNP-level TMRCA, gene annotation, and TE annotation.
- Computes:
    * one TMRCA per gene (mean of overlapping segments)
    * one TMRCA per TE (mean of overlapping segments)
    * identifies Copia and Gypsy TE subsets
- Bins TMRCA values into user-specified or automatic bins
- Computes frequencies of:
    * segments
    * SNPs
    * genes
    * all TEs
    * Copia TEs
    * Gypsy TEs
  across TMRCA bins.
- Outputs a single TSV: tmrca_bins_frequencies.tsv

New in V4:
- --manual-bins: manual bin edges (comma-separated)
- --dedup-overlapping-tes: collapse overlapping TE intervals by keeping only
  the longest TE per overlapping cluster (per scaffold).
- uses "gene" terminology instead of "mRNA".
"""

import argparse
import os
from typing import List, Optional

import numpy as np
import pandas as pd


# ------------------------- helpers -------------------------


def parse_manual_bins(bin_str: str) -> List[float]:
    edges = [float(x.strip()) for x in bin_str.split(",") if x.strip() != ""]
    if len(edges) < 2:
        raise ValueError("Need at least two values for --manual-bins.")
    if sorted(edges) != edges:
        raise ValueError("Bin edges must be sorted ascending.")
    return edges


def make_auto_bins(tmrca_vals: np.ndarray, n_bins: int = 100) -> List[float]:
    """Make log-spaced-ish bins between min>0 and max."""
    x = tmrca_vals[np.isfinite(tmrca_vals) & (tmrca_vals > 0)]
    if x.size == 0:
        raise ValueError("No positive finite TMRCA values to define bins.")
    lo = x.min()
    hi = x.max()
    # log-spaced edges
    edges = np.logspace(np.log10(lo), np.log10(hi), num=n_bins + 1)
    return list(edges)


def dedup_overlapping_tes(
    df: pd.DataFrame,
    chrom_col: str = "seqid",
    start_col: str = "start",
    end_col: str = "end",
) -> pd.DataFrame:
    """
    Deduplicate overlapping TE intervals.

    For each scaffold/chrom, sort by start, end. Walk through and collapse
    overlapping blocks; in each overlapping block, keep only the TE with the
    largest length (end - start + 1).
    """
    out_rows = []
    for chrom, sub in df.sort_values([chrom_col, start_col, end_col]).groupby(chrom_col):
        sub = sub.reset_index(drop=True)
        if sub.empty:
            continue
        current_idx = 0
        current_start = sub.loc[0, start_col]
        current_end = sub.loc[0, end_col]
        current_len = current_end - current_start + 1

        chosen_idx = 0  # index within sub

        for i in range(1, len(sub)):
            s = sub.loc[i, start_col]
            e = sub.loc[i, end_col]
            length = e - s + 1

            if s <= current_end:  # overlap
                # extend cluster end
                if e > current_end:
                    current_end = e
                # if this TE is longer, pick it as representative
                if length > current_len:
                    current_len = length
                    chosen_idx = i
            else:
                # flush chosen TE from previous cluster
                out_rows.append(sub.loc[chosen_idx].to_dict())
                # start new cluster
                current_idx = i
                current_start = s
                current_end = e
                current_len = length
                chosen_idx = i

        # flush last cluster
        out_rows.append(sub.loc[chosen_idx].to_dict())

    return pd.DataFrame(out_rows)


def compute_entity_tmrca(
    tmrca_df: pd.DataFrame,
    entity_df: pd.DataFrame,
    chrom_col_tmrca: str,
    start_col_tmrca: str,
    end_col_tmrca: str,
    mean_col_tmrca: str,
    chrom_col_ent: str,
    start_col_ent: str,
    end_col_ent: str,
    id_col_ent: str,
) -> pd.DataFrame:
    """
    For each entity (gene or TE), compute mean of mean_tmrca for overlapping
    segments.

    Returns a DataFrame with columns [id_col_ent, 'tmrca'].
    Entities with no overlapping segments get tmrca = NaN.
    """
    records = []

    # Ensure proper sorting for efficient overlap scanning
    tmrca_sorted = tmrca_df.sort_values([chrom_col_tmrca, start_col_tmrca])

    for chrom, ent_sub in entity_df.groupby(chrom_col_ent):
        t_sub = tmrca_sorted[tmrca_sorted[chrom_col_tmrca] == chrom]
        if t_sub.empty:
            # no segments on this chrom -> all NaN
            for _, row in ent_sub.iterrows():
                records.append(
                    {
                        id_col_ent: row[id_col_ent],
                        "tmrca": np.nan,
                    }
                )
            continue

        t_starts = t_sub[start_col_tmrca].to_numpy()
        t_ends = t_sub[end_col_tmrca].to_numpy()
        t_mean = t_sub[mean_col_tmrca].to_numpy()

        for _, row in ent_sub.iterrows():
            s = int(row[start_col_ent])
            e = int(row[end_col_ent])

            # segments with start < e
            right = np.searchsorted(t_starts, e, side="left")
            cand_idx = np.arange(right)
            if cand_idx.size == 0:
                tmrca_val = np.nan
            else:
                overlap_mask = t_ends[cand_idx] > s
                vals = t_mean[cand_idx[overlap_mask]]
                tmrca_val = np.nanmean(vals) if vals.size > 0 else np.nan

            records.append(
                {
                    id_col_ent: row[id_col_ent],
                    "tmrca": tmrca_val,
                }
            )

    return pd.DataFrame(records)


def assign_bins(values: np.ndarray, bin_edges: List[float]) -> pd.Series:
    """Return categorical bin labels for each value (NaN for out-of-range)."""
    bins = pd.IntervalIndex.from_breaks(bin_edges, closed="right")
    return pd.cut(values, bins=bins)


def fraction_per_bin(bin_series: pd.Series) -> pd.DataFrame:
    """
    Given a series of Interval (bins) with possible NaN, return:
    index: bins
    columns: count, fraction
    """
    counts = bin_series.value_counts(dropna=True).sort_index()
    total = counts.sum()
    frac = counts / total if total > 0 else counts * np.nan
    out = pd.DataFrame({"count": counts, "fraction": frac})
    out.index.name = "bin"
    return out


# ------------------------- main -------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Summarize TMRCA bins for segments, SNPs, genes, and TEs (with optional manual bins and TE dedup)."
    )
    ap.add_argument(
        "--tmrca",
        required=True,
        help="TMRCA segments file (TSV) with columns CHROM,start_tmrca,end_tmrca,mean_tmrca,...",
    )
    ap.add_argument(
        "--snps",
        required=True,
        help="SNP-level TMRCA file (TSV) with a column giving per-SNP TMRCA.",
    )
    ap.add_argument(
        "--snp-tmrca-col",
        default="mean_tmrca_at_pos",
        help="Column name in SNP file containing TMRCA per SNP (default: mean_tmrca_at_pos)",
    )
    ap.add_argument(
        "--genes",
        required=True,
        help="Gene annotation file (TSV), e.g. 1ab_mRNA_ID_merged_interproscan.txt",
    )
    ap.add_argument(
        "--tes",
        required=True,
        help="TE annotation file (TSV), e.g. col8_readable.TEanno.cds.scaf1ab_v1.tsv",
    )
    ap.add_argument(
        "--dedup-overlapping-tes",
        action="store_true",
        help="Collapse overlapping TEs on each scaffold and keep only the longest TE per overlapping block.",
    )
    ap.add_argument(
        "--manual-bins",
        type=str,
        default=None,
        help="Comma-separated list of manual TMRCA bin edges, e.g. '0,50,100,200,...,140000'. If not provided, automatic log-spaced bins are used.",
    )
    ap.add_argument(
        "--n-bins",
        type=int,
        default=100,
        help="Number of automatic bins if --manual-bins is not used (default: 100).",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Output directory (will be created if not exists).",
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---------------- read TMRCA segments ----------------
    print("Reading TMRCA segments:", args.tmrca)
    tmrca = pd.read_csv(args.tmrca, sep="\t")
    # Basic checks
    required_tmrca_cols = ["CHROM", "start_tmrca", "end_tmrca", "mean_tmrca"]
    for c in required_tmrca_cols:
        if c not in tmrca.columns:
            raise ValueError(f"TMRCA file missing column: {c}")

    # ---------------- read SNP TMRCA ----------------
    print("Reading SNP TMRCA:", args.snps)
    snps = pd.read_csv(args.snps, sep="\t")
    if args.snp_tmrca_col not in snps.columns:
        raise ValueError(
            f"SNP file missing TMRCA column '{args.snp_tmrca_col}'. Columns: {snps.columns.tolist()}"
        )
    snp_tmrca_vals = snps[args.snp_tmrca_col].astype(float).to_numpy()

    # ---------------- read genes ----------------
    print("Reading genes:", args.genes)
    genes = pd.read_csv(args.genes, sep="\t")
    required_gene_cols = ["new_seqid", "new_start", "new_end", "mRNA_ID"]
    for c in required_gene_cols:
        if c not in genes.columns:
            raise ValueError(f"Gene file missing column: {c}")
    # treat mRNA_ID as gene_id
    genes = genes.rename(columns={"mRNA_ID": "gene_id"})

    # ---------------- read TEs ----------------
    print("Reading TEs:", args.tes)
    tes = pd.read_csv(args.tes, sep="\t")
    required_te_cols = ["seqid", "sequence_ontology", "start", "end"]
    for c in required_te_cols:
        if c not in tes.columns:
            raise ValueError(f"TE file missing column: {c}")

    if args.dedup_overlapping_tes:
        print("Deduplicating overlapping TEs (keeping longest per overlapping block)...")
        n_before = len(tes)
        tes = dedup_overlapping_tes(tes)
        n_after = len(tes)
        print(f"  TEs before dedup: {n_before:,}; after dedup: {n_after:,}")

    # ---------------- compute gene TMRCA ----------------
    print("Computing per-gene TMRCA (mean of overlapping segments)...")
    gene_tmrca_df = compute_entity_tmrca(
        tmrca_df=tmrca,
        entity_df=genes,
        chrom_col_tmrca="CHROM",
        start_col_tmrca="start_tmrca",
        end_col_tmrca="end_tmrca",
        mean_col_tmrca="mean_tmrca",
        chrom_col_ent="new_seqid",
        start_col_ent="new_start",
        end_col_ent="new_end",
        id_col_ent="gene_id",
    )

    # ---------------- compute TE TMRCA ----------------
    print("Computing per-TE TMRCA (mean of overlapping segments)...")
    # give each TE an id
    tes = tes.copy()
    tes["te_id"] = np.arange(len(tes))
    te_tmrca_df = compute_entity_tmrca(
        tmrca_df=tmrca,
        entity_df=tes,
        chrom_col_tmrca="CHROM",
        start_col_tmrca="start_tmrca",
        end_col_tmrca="end_tmrca",
        mean_col_tmrca="mean_tmrca",
        chrom_col_ent="seqid",
        start_col_ent="start",
        end_col_ent="end",
        id_col_ent="te_id",
    )
    tes = tes.merge(te_tmrca_df, on="te_id", how="left")

    # Identify Copia and Gypsy subsets
    tes["is_copia"] = tes["sequence_ontology"].str.contains("Copia", case=False, na=False)
    tes["is_gypsy"] = tes["sequence_ontology"].str.contains("Gypsy", case=False, na=False)

    # ---------------- define bins ----------------
    print("Defining TMRCA bins...")
    seg_tmrca_vals = tmrca["mean_tmrca"].astype(float).to_numpy()
    if args.manual_bins is not None:
        bin_edges = parse_manual_bins(args.manual_bins)
        print("  Using manual bins with edges:", bin_edges)
    else:
        bin_edges = make_auto_bins(seg_tmrca_vals, n_bins=args.n_bins)
        print(f"  Using {args.n_bins} automatic log-spaced bins "
              f"from {min(bin_edges):.3g} to {max(bin_edges):.3g}")

    # ---------------- assign bins ----------------
    print("Assigning bins to segments, SNPs, genes, and TEs...")

    seg_bins = assign_bins(seg_tmrca_vals, bin_edges)
    snp_bins = assign_bins(snp_tmrca_vals, bin_edges)

    gene_bins = assign_bins(gene_tmrca_df["tmrca"].to_numpy(dtype=float), bin_edges)

    te_bins_all = assign_bins(tes["tmrca"].to_numpy(dtype=float), bin_edges)
    te_bins_copia = assign_bins(
        tes.loc[tes["is_copia"], "tmrca"].to_numpy(dtype=float), bin_edges
    )
    te_bins_gypsy = assign_bins(
        tes.loc[tes["is_gypsy"], "tmrca"].to_numpy(dtype=float), bin_edges
    )

    # ---------------- frequencies ----------------
    seg_freq = fraction_per_bin(seg_bins).rename(
        columns={"count": "segments_count", "fraction": "segments_fraction"}
    )
    snp_freq = fraction_per_bin(snp_bins).rename(
        columns={"count": "snps_count", "fraction": "snps_fraction"}
    )
    gene_freq = fraction_per_bin(gene_bins).rename(
        columns={"count": "genes_count", "fraction": "genes_fraction"}
    )
    te_freq_all = fraction_per_bin(te_bins_all).rename(
        columns={"count": "tes_count", "fraction": "tes_fraction"}
    )
    te_freq_copia = fraction_per_bin(te_bins_copia).rename(
        columns={"count": "copia_count", "fraction": "copia_fraction"}
    )
    te_freq_gypsy = fraction_per_bin(te_bins_gypsy).rename(
        columns={"count": "gypsy_count", "fraction": "gypsy_fraction"}
    )

    # join on bin (Interval)
    freq_all = (
        seg_freq.join(snp_freq, how="outer")
        .join(gene_freq, how="outer")
        .join(te_freq_all, how="outer")
        .join(te_freq_copia, how="outer")
        .join(te_freq_gypsy, how="outer")
        .sort_index()
    )

    # unpack Interval to numeric bounds & a label
    freq_all = freq_all.reset_index()
    freq_all["bin_left"] = freq_all["bin"].apply(lambda x: x.left if pd.notna(x) else np.nan)
    freq_all["bin_right"] = freq_all["bin"].apply(lambda x: x.right if pd.notna(x) else np.nan)
    freq_all["bin_label"] = freq_all["bin"].astype(str)

    # re-order columns
    col_order = [
        "bin_label",
        "bin_left",
        "bin_right",
        "segments_count",
        "segments_fraction",
        "snps_count",
        "snps_fraction",
        "genes_count",
        "genes_fraction",
        "tes_count",
        "tes_fraction",
        "copia_count",
        "copia_fraction",
        "gypsy_count",
        "gypsy_fraction",
    ]
    # keep only those that exist
    col_order = [c for c in col_order if c in freq_all.columns]
    freq_all = freq_all[col_order]

    out_tsv = os.path.join(args.outdir, "tmrca_bins_frequencies.tsv")
    freq_all.to_csv(out_tsv, sep="\t", index=False)
    print("Wrote:", out_tsv)
    print("Done.")


if __name__ == "__main__":
    main()

