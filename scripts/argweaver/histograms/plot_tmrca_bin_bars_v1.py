#!/usr/bin/env python3
"""
plot_tmrca_bin_bars_v1.py

Make overlapping bar charts from tmrca_bins_frequencies.tsv:

1) SNPs vs (mRNA + TE)
2) SNPs vs (Copia + Gypsy)
3) SNPs vs mRNA
4) SNPs vs TE
5) segments vs (mRNA + TE)
6) segments vs (Copia + Gypsy)
7) segments vs mRNA
8) segments vs TE
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_figsize(s):
    parts = str(s).split(",")
    if len(parts) != 2:
        return (10, 6)
    return (float(parts[0]), float(parts[1]))


def plot_combo(df, x, bin_labels,
               left_col, right_cols,
               left_label, right_labels,
               color_left, colors_right,
               alpha_left, alphas_right,
               fontsize, figsize, outpath,
               title="", left_ylabel="", right_ylabel=""):
    """
    Generic plotting helper:
    - left_col goes to ax1 (left y-axis)
    - right_cols (1 or 2) go to ax2 (right y-axis)
    - Bars overlap at same x, but with different colors/alphas.
    """
    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()

    # left bars
    y_left = df[left_col].to_numpy()
    ax1.bar(x, y_left,
            color=color_left,
            alpha=alpha_left,
            width=0.8,
            label=left_label)

    # right bars (one or two series)
    for i, (col, c, a, lab) in enumerate(zip(right_cols, colors_right, alphas_right, right_labels)):
        if col is None:
            continue
        y_right = df[col].to_numpy()
        ax2.bar(x, y_right,
                color=c,
                alpha=a,
                width=0.5,
                label=lab)

    # x-axis labels (maybe sparse if many bins)
    ax1.set_xticks(x)
    if len(x) > 40:
        step = max(1, len(x) // 40)
        ax1.set_xticks(x[::step])
        ax1.set_xticklabels(bin_labels[::step], rotation=90, fontsize=fontsize*0.7)
    else:
        ax1.set_xticklabels(bin_labels, rotation=90, fontsize=fontsize*0.7)

    ax1.set_ylabel(left_ylabel, fontsize=fontsize)
    ax2.set_ylabel(right_ylabel, fontsize=fontsize)
    ax1.set_xlabel("TMRCA bin", fontsize=fontsize)
    ax1.set_title(title, fontsize=fontsize+2)

    # Combine legends from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = labels1 + labels2
    ax1.legend(handles, labels, fontsize=fontsize*0.8, loc="upper right")

    ax1.tick_params(axis='both', labelsize=fontsize*0.8)
    ax2.tick_params(axis='y', labelsize=fontsize*0.8)

    fig.tight_layout()
    outpath = Path(outpath)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Make overlapping bar plots from tmrca_bins_frequencies.tsv"
    )
    ap.add_argument("--bins-file", required=True,
                    help="tmrca_bins_frequencies.tsv from TMRCA_bins_genes_TEs_SNPs_V3.py")
    ap.add_argument("--out-prefix", default="tmrca_bins_plots",
                    help="Prefix for output image files (PNG).")

    # appearance options
    ap.add_argument("--figsize", default="12,6",
                    help="Figure size as 'width,height' (default: 12,6).")
    ap.add_argument("--fontsize", type=float, default=10.0,
                    help="Base font size (default: 10).")

    # colors
    ap.add_argument("--color-snp", default="grey",
                    help="Color for SNP bars (default: grey).")
    ap.add_argument("--color-seg", default="grey",
                    help="Color for segment bars (default: grey).")
    ap.add_argument("--color-mrna", default="blue",
                    help="Color for mRNA bars (default: blue).")
    ap.add_argument("--color-te", default="red",
                    help="Color for TE bars (default: red).")
    ap.add_argument("--color-copia", default="pink",
                    help="Color for Copia bars (default: pink).")
    ap.add_argument("--color-gypsy", default="green",
                    help="Color for Gypsy bars (default: green).")

    # transparency
    ap.add_argument("--alpha-left", type=float, default=0.6,
                    help="Alpha (transparency) for left-axis bars (default: 0.6).")
    ap.add_argument("--alpha-right1", type=float, default=0.5,
                    help="Alpha for first right-axis series (default: 0.5).")
    ap.add_argument("--alpha-right2", type=float, default=0.5,
                    help="Alpha for second right-axis series (default: 0.5).")

    args = ap.parse_args()

    df = pd.read_csv(args.bins_file, sep="\t")
    # Expect columns: bin_label, frac_segments, frac_snps, frac_genes,
    # frac_TEs_all, frac_TEs_Copia, frac_TEs_Gypsy
    required_cols = [
        "bin_label",
        "frac_segments",
        "frac_snps",
        "frac_genes",
        "frac_TEs_all",
        "frac_TEs_Copia",
        "frac_TEs_Gypsy",
    ]
    for c in required_cols:
        if c not in df.columns:
            raise SystemExit(f"ERROR: Column '{c}' not found in {args.bins_file}")

    x = np.arange(len(df))
    bin_labels = df["bin_label"].astype(str).tolist()
    figsize = parse_figsize(args.figsize)
    fs = args.fontsize

    # 1) SNPs vs (mRNA + TE)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_snps",
        right_cols=["frac_genes", "frac_TEs_all"],
        left_label="SNPs",
        right_labels=["mRNA", "TE (all)"],
        color_left=args.color_snp,
        colors_right=[args.color_mrna, args.color_te],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1, args.alpha_right2],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_1_snps_vs_mrna_te.png",
        title="SNPs vs mRNA & TE across TMRCA bins",
        left_ylabel="Fraction of SNPs",
        right_ylabel="Fraction of mRNA / TE"
    )

    # 2) SNPs vs (Copia + Gypsy)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_snps",
        right_cols=["frac_TEs_Copia", "frac_TEs_Gypsy"],
        left_label="SNPs",
        right_labels=["Copia", "Gypsy"],
        color_left=args.color_snp,
        colors_right=[args.color_copia, args.color_gypsy],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1, args.alpha_right2],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_2_snps_vs_copia_gypsy.png",
        title="SNPs vs Copia & Gypsy across TMRCA bins",
        left_ylabel="Fraction of SNPs",
        right_ylabel="Fraction of Copia / Gypsy"
    )

    # 3) SNPs vs mRNA
    plot_combo(
        df, x, bin_labels,
        left_col="frac_snps",
        right_cols=["frac_genes"],
        left_label="SNPs",
        right_labels=["mRNA"],
        color_left=args.color_snp,
        colors_right=[args.color_mrna],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_3_snps_vs_mrna.png",
        title="SNPs vs mRNA across TMRCA bins",
        left_ylabel="Fraction of SNPs",
        right_ylabel="Fraction of mRNA"
    )

    # 4) SNPs vs TE (all)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_snps",
        right_cols=["frac_TEs_all"],
        left_label="SNPs",
        right_labels=["TE (all)"],
        color_left=args.color_snp,
        colors_right=[args.color_te],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_4_snps_vs_te.png",
        title="SNPs vs TE (all) across TMRCA bins",
        left_ylabel="Fraction of SNPs",
        right_ylabel="Fraction of TE (all)"
    )

    # 5) segments vs (mRNA + TE)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_segments",
        right_cols=["frac_genes", "frac_TEs_all"],
        left_label="Segments",
        right_labels=["mRNA", "TE (all)"],
        color_left=args.color_seg,
        colors_right=[args.color_mrna, args.color_te],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1, args.alpha_right2],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_5_segments_vs_mrna_te.png",
        title="Segments vs mRNA & TE across TMRCA bins",
        left_ylabel="Fraction of segments",
        right_ylabel="Fraction of mRNA / TE"
    )

    # 6) segments vs (Copia + Gypsy)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_segments",
        right_cols=["frac_TEs_Copia", "frac_TEs_Gypsy"],
        left_label="Segments",
        right_labels=["Copia", "Gypsy"],
        color_left=args.color_seg,
        colors_right=[args.color_copia, args.color_gypsy],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1, args.alpha_right2],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_6_segments_vs_copia_gypsy.png",
        title="Segments vs Copia & Gypsy across TMRCA bins",
        left_ylabel="Fraction of segments",
        right_ylabel="Fraction of Copia / Gypsy"
    )

    # 7) segments vs mRNA
    plot_combo(
        df, x, bin_labels,
        left_col="frac_segments",
        right_cols=["frac_genes"],
        left_label="Segments",
        right_labels=["mRNA"],
        color_left=args.color_seg,
        colors_right=[args.color_mrna],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_7_segments_vs_mrna.png",
        title="Segments vs mRNA across TMRCA bins",
        left_ylabel="Fraction of segments",
        right_ylabel="Fraction of mRNA"
    )

    # 8) segments vs TE (all)
    plot_combo(
        df, x, bin_labels,
        left_col="frac_segments",
        right_cols=["frac_TEs_all"],
        left_label="Segments",
        right_labels=["TE (all)"],
        color_left=args.color_seg,
        colors_right=[args.color_te],
        alpha_left=args.alpha_left,
        alphas_right=[args.alpha_right1],
        fontsize=fs,
        figsize=figsize,
        outpath=f"{args.out_prefix}_8_segments_vs_te.png",
        title="Segments vs TE (all) across TMRCA bins",
        left_ylabel="Fraction of segments",
        right_ylabel="Fraction of TE (all)"
    )

    print("Done. Figures written with prefix:", args.out_prefix)


if __name__ == "__main__":
    main()

