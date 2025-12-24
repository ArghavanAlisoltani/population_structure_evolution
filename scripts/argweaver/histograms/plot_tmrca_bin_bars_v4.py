#!/usr/bin/env python3
"""
Plot TMRCA-bin bar and area plots with dual y-axes.

v4:
- Same logic as v3 / original v1
- Area plots are smooth curves (no step)
- Optional rolling smoothing for area plots (--smooth-area-window)
- Only one y-label per axis side (for broken-axes case, labels only on lower axes)
- Option to use upper-bin values for x-axis labels (e.g. bin_right)
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_ybreak(s):
    """Parse 'low,high' into (low, high) or None."""
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None
    parts = s.split(",")
    if len(parts) != 2:
        raise ValueError(f"Could not parse ybreak '{s}', expected 'low,high'")
    low, high = [float(p.strip()) for p in parts]
    if not (low < high):
        raise ValueError(f"ybreak low must be < high, got {low}, {high}")
    return low, high


def add_y_break_marks(ax_upper, ax_lower, d=0.02):
    """
    Draw diagonal 'break' marks between two stacked axes.
    Purely cosmetic.
    """
    kwargs = dict(transform=ax_upper.transAxes, color='k', clip_on=False, linewidth=0.7)
    # bottom of upper axis
    ax_upper.plot((-d, +d), (-d, +d), **kwargs)
    ax_upper.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    # top of lower axis
    kwargs.update(transform=ax_lower.transAxes)
    ax_lower.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_lower.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)


def smooth_series(y, window):
    """Rolling mean smoothing (centered) for area plots."""
    if window is None or window <= 1:
        return y
    s = pd.Series(y)
    return s.rolling(window=window, center=True, min_periods=1).mean().to_numpy()


def dual_axis_plot(
    df,
    bin_col,
    labels,
    main_col,
    sec_cols,
    sec_types,
    out_prefix,
    suffix,
    main_label,
    sec_label,
    colors,
    alphas,
    widths,
    x_gap=0.01,
    ybreak_main=None,
    make_area=False,
    base_height=5.0,
    smooth_area_window=1,
    legend_loc="upper right",
):
    """
    Create bar (and optionally area) plots with dual y-axis and optional y-break.

    colors: dict type -> color (snp, seg, gene, te, copia, gypsy)
    alphas: dict type -> alpha
    widths: dict type -> width
    legend_loc:
        - standard matplotlib strings, e.g. 'upper right', 'upper left', 'best'
        - 'right-out' to place legend outside on the right
        - 'none' to hide legend
    """

    # X positions / labels
    x = np.arange(len(df))

    # Infer main type (only used for picking default color/alpha/width)
    main_type = "snp" if "snp" in main_col.lower() else (
        "seg" if "seg" in main_col.lower() else "snp"
    )
    main_color = colors[main_type]
    main_alpha = alphas[main_type]
    main_width = widths[main_type]

    # y values
    y_main = df[main_col].values
    max_main = np.nanmax(y_main) if len(y_main) else 1.0

    y_sec_all = []
    for sc in sec_cols:
        if sc in df.columns:
            y_sec_all.append(df[sc].values)
    max_sec = np.nanmax(np.concatenate(y_sec_all)) if y_sec_all else 1.0

    # Precompute smoothed arrays for area plots
    y_main_area = smooth_series(y_main, smooth_area_window)
    y_sec_area = {}
    for sc in sec_cols:
        if sc in df.columns:
            y_sec_area[sc] = smooth_series(df[sc].values, smooth_area_window)

    def _make_axes_no_break():
        fig, ax_main = plt.subplots(figsize=(12, base_height))
        ax_sec = ax_main.twinx()
        return fig, ax_main, ax_sec

    def _make_axes_with_break(low, high):
        fig, (ax_upper_main, ax_lower_main) = plt.subplots(
            2, 1,
            sharex=True,
            figsize=(12, base_height * 1.4),
            gridspec_kw={"height_ratios": [1, 2]}
        )
        ax_upper_sec = ax_upper_main.twinx()
        ax_lower_sec = ax_lower_main.twinx()

        # Set y-limits for main axis
        ax_lower_main.set_ylim(0, low)
        ax_upper_main.set_ylim(high, max_main)

        # Set y-limits for secondary axis (same numeric range; fractions)
        ax_lower_sec.set_ylim(0, low)
        ax_upper_sec.set_ylim(high, max_sec)

        ax_upper_main.tick_params(labelbottom=False)
        add_y_break_marks(ax_upper_main, ax_lower_main)

        return fig, (ax_upper_main, ax_lower_main), (ax_upper_sec, ax_lower_sec)

    def _plot_bars(ax_main, ax_sec):
        # main bars
        ax_main.bar(
            x,
            df[main_col],
            width=main_width,
            color=main_color,
            alpha=main_alpha,
            zorder=2,
        )

        # secondary bars
        for sc, typ in zip(sec_cols, sec_types):
            if sc not in df.columns:
                continue
            ax_sec.bar(
                x,
                df[sc],
                width=widths[typ],
                color=colors[typ],
                alpha=alphas[typ],
                zorder=1,
            )

        # Bring main axis to front
        ax_main.set_zorder(ax_sec.get_zorder() + 1)
        ax_main.patch.set_visible(False)

    def _plot_area(ax_main, ax_sec):
        # main area + smooth line
        ax_main.fill_between(
            x,
            0,
            y_main_area,
            color=main_color,
            alpha=main_alpha,
        )
        ax_main.plot(
            x,
            y_main_area,
            color=main_color,
            alpha=min(1.0, main_alpha + 0.2),
            linewidth=2,
            label=main_label,
        )

        # secondary areas + smooth lines
        for sc, typ in zip(sec_cols, sec_types):
            if sc not in df.columns:
                continue
            y_vals = y_sec_area.get(sc, df[sc].values)
            ax_sec.fill_between(
                x,
                0,
                y_vals,
                color=colors[typ],
                alpha=alphas[typ],
            )
            ax_sec.plot(
                x,
                y_vals,
                color=colors[typ],
                alpha=min(1.0, alphas[typ] + 0.2),
                linewidth=2,
                label=typ,
            )

        ax_main.set_zorder(ax_sec.get_zorder() + 1)
        ax_main.patch.set_visible(False)

    def _add_legend(fig, ax_main_for_legend, ax_sec_for_legend):
        if legend_loc is None:
            return
        loc_str = str(legend_loc).lower()
        if loc_str == "none":
            return

        handles_main, labels_main = ax_main_for_legend.get_legend_handles_labels()
        handles_sec, labels_sec = ax_sec_for_legend.get_legend_handles_labels()
        if not (handles_main or handles_sec):
            return

        if loc_str == "right-out":
            fig.legend(
                handles_main + handles_sec,
                labels_main + labels_sec,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
                fontsize=9,
            )
        else:
            fig.legend(
                handles_main + handles_sec,
                labels_main + labels_sec,
                loc=legend_loc,
                fontsize=9,
            )

    # ---------- BAR PLOT ----------
    if ybreak_main is None:
        fig, ax_main, ax_sec = _make_axes_no_break()
        _plot_bars(ax_main, ax_sec)

        # Labels and ticks – single set
        ax_main.set_ylabel(main_label)
        ax_sec.set_ylabel(sec_label)
        ax_main.set_xticks(x)
        ax_main.set_xticklabels(labels, rotation=90)
        ax_main.set_xlabel("TMRCA bin")
        ax_main.margins(x=x_gap)
        ax_sec.margins(x=x_gap)

        ax_main.set_title(suffix.replace("_", " ").lstrip("_"))
        _add_legend(fig, ax_main, ax_sec)

        fig.tight_layout()
        outfile = f"{out_prefix}{suffix}_bar.png"
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        low, high = ybreak_main
        fig, (ax_upper_main, ax_lower_main), (ax_upper_sec, ax_lower_sec) = _make_axes_with_break(low, high)

        _plot_bars(ax_upper_main, ax_upper_sec)
        _plot_bars(ax_lower_main, ax_lower_sec)

        # Only LOWER axes get labels & ticks
        ax_lower_main.set_ylabel(main_label)
        ax_lower_sec.set_ylabel(sec_label)
        ax_lower_main.set_xticks(x)
        ax_lower_main.set_xticklabels(labels, rotation=90)
        ax_lower_main.set_xlabel("TMRCA bin")

        # Remove y-axis labels from upper axes
        ax_upper_main.set_ylabel("")
        ax_upper_sec.set_ylabel("")
        ax_upper_main.tick_params(axis='y', labelleft=False)
        ax_upper_sec.tick_params(axis='y', labelright=False)

        ax_lower_main.margins(x=x_gap)
        ax_lower_sec.margins(x=x_gap)

        ax_upper_main.set_title(suffix.replace("_", " ").lstrip("_"))
        _add_legend(fig, ax_upper_main, ax_upper_sec)

        fig.tight_layout()
        outfile = f"{out_prefix}{suffix}_bar.png"
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)

    # ---------- AREA PLOT ----------
    if make_area:
        if ybreak_main is None:
            fig, ax_main, ax_sec = _make_axes_no_break()
            _plot_area(ax_main, ax_sec)

            ax_main.set_ylabel(main_label)
            ax_sec.set_ylabel(sec_label)
            ax_main.set_xticks(x)
            ax_main.set_xticklabels(labels, rotation=90)
            ax_main.set_xlabel("TMRCA bin")
            ax_main.margins(x=x_gap)
            ax_sec.margins(x=x_gap)

            ax_main.set_title(suffix.replace("_", " ").lstrip("_") + " (area)")
            _add_legend(fig, ax_main, ax_sec)

            fig.tight_layout()
            outfile = f"{out_prefix}{suffix}_area.png"
            fig.savefig(outfile, dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            low, high = ybreak_main
            fig, (ax_upper_main, ax_lower_main), (ax_upper_sec, ax_lower_sec) = _make_axes_with_break(low, high)

            _plot_area(ax_upper_main, ax_upper_sec)
            _plot_area(ax_lower_main, ax_lower_sec)

            # Labels only on lower axes
            ax_lower_main.set_ylabel(main_label)
            ax_lower_sec.set_ylabel(sec_label)
            ax_lower_main.set_xticks(x)
            ax_lower_main.set_xticklabels(labels, rotation=90)
            ax_lower_main.set_xlabel("TMRCA bin")

            ax_upper_main.set_ylabel("")
            ax_upper_sec.set_ylabel("")
            ax_upper_main.tick_params(axis='y', labelleft=False)
            ax_upper_sec.tick_params(axis='y', labelright=False)

            ax_lower_main.margins(x=x_gap)
            ax_lower_sec.margins(x=x_gap)

            ax_upper_main.set_title(suffix.replace("_", " ").lstrip("_") + " (area)")
            _add_legend(fig, ax_upper_main, ax_upper_sec)

            fig.tight_layout()
            outfile = f"{out_prefix}{suffix}_area.png"
            fig.savefig(outfile, dpi=300, bbox_inches="tight")
            plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Plot TMRCA bin bar + optional area charts with dual y-axes."
    )
    ap.add_argument("--input", "-i", required=True,
                    help="Input tmrca_bins_frequencies.tsv")
    ap.add_argument("--out-prefix", "-o", required=True,
                    help="Output prefix (no extension, e.g. 'tmrca_bins_plots_')")

    # Columns
    ap.add_argument("--col-bin", default="TMRCA_bin",
                    help="Bin label column [default: TMRCA_bin]")
    ap.add_argument("--col-bin-upper", default="",
                    help="Optional column giving upper bin edges for x-axis labels (e.g. bin_right)")
    ap.add_argument("--col-snp-frac", default="frac_snps",
                    help="SNP fraction column [default: frac_snps]")
    ap.add_argument("--col-seg-frac", default="frac_segments",
                    help="Segment fraction column [default: frac_segments]")
    ap.add_argument("--col-gene-frac", default="frac_genes",
                    help="Gene fraction column [default: frac_genes]")
    ap.add_argument("--col-te-frac", default="frac_te_all",
                    help="TE(all) fraction column [default: frac_te_all]")
    ap.add_argument("--col-copia-frac", default="frac_copia",
                    help="Copia fraction column [default: frac_copia]")
    ap.add_argument("--col-gypsy-frac", default="frac_gypsy",
                    help="Gypsy fraction column [default: frac_gypsy]")

    # Visual parameters
    ap.add_argument("--x-gap", type=float, default=0.01,
                    help="Fractional x-axis margin [default: 0.01]")
    ap.add_argument("--base-height", type=float, default=5.0,
                    help="Base figure height for non-broken plots [default: 5.0]")

    # Colors
    ap.add_argument("--col-snp", default="#BDBDBD",
                    help="Color for SNP bars [default: #BDBDBD]")
    ap.add_argument("--col-seg", default="#BDBDBD",
                    help="Color for segment bars [default: #BDBDBD]")
    ap.add_argument("--col-gene", default="#377EB8",
                    help="Color for gene bars [default: #377EB8]")
    ap.add_argument("--col-te", default="#E41A1C",
                    help="Color for TE(all) bars [default: #E41A1C]")
    ap.add_argument("--col-copia", default="#F781BF",
                    help="Color for Copia bars [default: #F781BF]")
    ap.add_argument("--col-gypsy", default="#4DAF4A",
                    help="Color for Gypsy bars [default: #4DAF4A]")

    # Alphas
    ap.add_argument("--alpha-snp", type=float, default=0.7,
                    help="Alpha for SNP bars [default: 0.7]")
    ap.add_argument("--alpha-seg", type=float, default=0.7,
                    help="Alpha for segment bars [default: 0.7]")
    ap.add_argument("--alpha-gene", type=float, default=0.7,
                    help="Alpha for gene bars [default: 0.7]")
    ap.add_argument("--alpha-te", type=float, default=0.4,
                    help="Alpha for TE(all) bars [default: 0.4]")
    ap.add_argument("--alpha-copia", type=float, default=0.4,
                    help="Alpha for Copia bars [default: 0.4]")
    ap.add_argument("--alpha-gypsy", type=float, default=0.4,
                    help="Alpha for Gypsy bars [default: 0.4]")

    # Bar widths
    ap.add_argument("--width-snp", type=float, default=0.8,
                    help="Bar width for SNP bars [default: 0.8]")
    ap.add_argument("--width-seg", type=float, default=0.8,
                    help="Bar width for segment bars [default: 0.8]")
    ap.add_argument("--width-gene", type=float, default=0.8,
                    help="Bar width for gene bars [default: 0.8]")
    ap.add_argument("--width-te", type=float, default=0.8,
                    help="Bar width for TE(all) bars [default: 0.8]")
    ap.add_argument("--width-copia", type=float, default=0.8,
                    help="Bar width for Copia bars [default: 0.8]")
    ap.add_argument("--width-gypsy", type=float, default=0.8,
                    help="Bar width for Gypsy bars [default: 0.8]")

    # Y-axis break
    ap.add_argument("--ybreak-main", default="",
                    help="Optional main y-axis break as 'low,high' [default: none]")

    # Area plots
    ap.add_argument("--make-area", action="store_true",
                    help="Also make area plots (in addition to bar plots)")

    # Smoothing for area plots
    ap.add_argument("--smooth-area-window", type=int, default=1,
                    help="Window size (in bins) for smoothing area plots; 1 = no smoothing [default: 1]")

    # Legend position
    ap.add_argument(
        "--legend-loc",
        default="upper right",
        help=(
            "Legend location; examples: 'upper right', 'upper left', 'lower right', "
            "'best', 'right-out', or 'none' to hide [default: 'upper right']"
        )
    )

    args = ap.parse_args()

    # Load data
    df = pd.read_csv(args.input, sep="\t", header=0)

    # Ensure required columns exist
    required_cols = [
        args.col_bin, args.col_snp_frac, args.col_seg_frac,
        args.col_gene_frac, args.col_te_frac,
        args.col_copia_frac, args.col_gypsy_frac
    ]
    for c in required_cols:
        if c not in df.columns:
            raise SystemExit(f"Required column '{c}' not found in input.")

    # Preserve file order for bins
    df = df.copy()
    df[args.col_bin] = pd.Categorical(
        df[args.col_bin],
        categories=df[args.col_bin],
        ordered=True,
    )

    # Decide x-axis labels:
    # - If upper-bin column given, use that (e.g. bin_right: 50, 200, 400...)
    # - Otherwise, use bin label column (e.g. "[0,50)")
    if args.col_bin_upper and args.col_bin_upper in df.columns:
        labels = df[args.col_bin_upper].astype(str).tolist()
    else:
        labels = df[args.col_bin].astype(str).tolist()

    colors = {
        "snp": args.col_snp,
        "seg": args.col_seg,
        "gene": args.col_gene,
        "te": args.col_te,
        "copia": args.col_copia,
        "gypsy": args.col_gypsy,
    }
    alphas = {
        "snp": args.alpha_snp,
        "seg": args.alpha_seg,
        "gene": args.alpha_gene,
        "te": args.alpha_te,
        "copia": args.alpha_copia,
        "gypsy": args.alpha_gypsy,
    }
    widths = {
        "snp": args.width_snp,
        "seg": args.width_seg,
        "gene": args.width_gene,
        "te": args.width_te,
        "copia": args.width_copia,
        "gypsy": args.width_gypsy,
    }

    ybreak_main = parse_ybreak(args.ybreak_main) if args.ybreak_main else None

    bin_col = args.col_bin
    snp_col = args.col_snp_frac
    seg_col = args.col_seg_frac
    gene_col = args.col_gene_frac
    te_col = args.col_te_frac
    copia_col = args.col_copia_frac
    gypsy_col = args.col_gypsy_frac

    # 1) SNPs vs genes + TE
    dual_axis_plot(
        df, bin_col, labels,
        main_col=snp_col,
        sec_cols=[gene_col, te_col],
        sec_types=["gene", "te"],
        out_prefix=args.out_prefix,
        suffix="_1_snps_vs_genes_te",
        main_label="Fraction of SNPs",
        sec_label="Fraction of genes / TE",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 2) SNPs vs Copia + Gypsy
    dual_axis_plot(
        df, bin_col, labels,
        main_col=snp_col,
        sec_cols=[copia_col, gypsy_col],
        sec_types=["copia", "gypsy"],
        out_prefix=args.out_prefix,
        suffix="_2_snps_vs_copia_gypsy",
        main_label="Fraction of SNPs",
        sec_label="Fraction of Copia / Gypsy",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 3) SNPs vs genes only
    dual_axis_plot(
        df, bin_col, labels,
        main_col=snp_col,
        sec_cols=[gene_col],
        sec_types=["gene"],
        out_prefix=args.out_prefix,
        suffix="_3_snps_vs_genes",
        main_label="Fraction of SNPs",
        sec_label="Fraction of genes",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 4) SNPs vs TE only
    dual_axis_plot(
        df, bin_col, labels,
        main_col=snp_col,
        sec_cols=[te_col],
        sec_types=["te"],
        out_prefix=args.out_prefix,
        suffix="_4_snps_vs_te",
        main_label="Fraction of SNPs",
        sec_label="Fraction of TE (all)",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 5) Segments vs genes + TE
    dual_axis_plot(
        df, bin_col, labels,
        main_col=seg_col,
        sec_cols=[gene_col, te_col],
        sec_types=["gene", "te"],
        out_prefix=args.out_prefix,
        suffix="_5_segments_vs_genes_te",
        main_label="Fraction of segments",
        sec_label="Fraction of genes / TE",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 6) Segments vs Copia + Gypsy
    dual_axis_plot(
        df, bin_col, labels,
        main_col=seg_col,
        sec_cols=[copia_col, gypsy_col],
        sec_types=["copia", "gypsy"],
        out_prefix=args.out_prefix,
        suffix="_6_segments_vs_copia_gypsy",
        main_label="Fraction of segments",
        sec_label="Fraction of Copia / Gypsy",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 7) Segments vs genes only
    dual_axis_plot(
        df, bin_col, labels,
        main_col=seg_col,
        sec_cols=[gene_col],
        sec_types=["gene"],
        out_prefix=args.out_prefix,
        suffix="_7_segments_vs_genes",
        main_label="Fraction of segments",
        sec_label="Fraction of genes",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 8) Segments vs TE only
    dual_axis_plot(
        df, bin_col, labels,
        main_col=seg_col,
        sec_cols=[te_col],
        sec_types=["te"],
        out_prefix=args.out_prefix,
        suffix="_8_segments_vs_te",
        main_label="Fraction of segments",
        sec_label="Fraction of TE (all)",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )

    # 9) SNPs vs genes + Copia + Gypsy
    dual_axis_plot(
        df, bin_col, labels,
        main_col=snp_col,
        sec_cols=[gene_col, copia_col, gypsy_col],
        sec_types=["gene", "copia", "gypsy"],
        out_prefix=args.out_prefix,
        suffix="_9_snps_vs_genes_copia_gypsy",
        main_label="Fraction of SNPs",
        sec_label="Fraction of genes / Copia / Gypsy",
        colors=colors,
        alphas=alphas,
        widths=widths,
        x_gap=args.x_gap,
        ybreak_main=ybreak_main,
        make_area=args.make_area,
        base_height=args.base_height,
        smooth_area_window=args.smooth_area_window,
        legend_loc=args.legend_loc,
    )


if __name__ == "__main__":
    main()

