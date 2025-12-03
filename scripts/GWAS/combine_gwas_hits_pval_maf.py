#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def get_trait_name(filename: str) -> str:
    """
    Extract trait name from filename.
    Examples:
      'C13_sumstat.txt' -> 'C13'
      'DBH30_sumstat.txt' -> 'DBH30'
    """
    base = Path(filename).name
    if base.endswith("_sumstat.txt"):
        return base[:-len("_sumstat.txt")]
    elif base.endswith(".txt"):
        return base[:-4]
    else:
        return base


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine GWAS hits across traits using P-value and MAF thresholds. "
            "Outputs two combined files in a new folder: "
            "1) simple format (position, trait, scaffold, SNP), "
            "2) full format with all original columns + trait."
        )
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="Input summary statistic files (e.g. C13_sumstat.txt DBH30_sumstat.txt ...)",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=1e-6,
        help="Maximum P-value to call a SNP significant (default: 1e-6, i.e. P <= 1e-6).",
    )
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help=(
            "Minimum MAF to keep a SNP (default: 0.01, i.e. MAF >= 0.01). "
            "Set lower if you want to include rarer variants."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="combined_pval_maf_hits",
        help="Output directory to write combined files (default: combined_pval_maf_hits).",
    )

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    simple_hits_list = []
    full_hits_list = []

    for fname in args.files:
        path = Path(fname)
        trait = get_trait_name(path.name)

        print(f"\nProcessing {path} (trait = {trait})")

        # Read file (space/tab-delimited)
        df = pd.read_csv(path, delim_whitespace=True)

        required_cols = {
            "SNP", "CHR", "BP", "A1", "A2", "MAF",
            "BETA", "SE", "Z", "P", "N"
        }
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"File {fname} is missing required columns: {missing}")

        # Drop rows with missing MAF or P
        df = df[df["MAF"].notna() & df["P"].notna()].copy()

        # Apply MAF (minimum) and P-value thresholds
        mask = (df["MAF"] >= args.maf_threshold) & (df["P"] <= args.p_threshold)
        hits = df.loc[mask].copy()

        print(f"  SNPs passing MAF >= {args.maf_threshold} and P <= {args.p_threshold}: {hits.shape[0]}")

        if hits.empty:
            continue

        # Add trait column
        hits["trait"] = trait

        # Save for full combined file (all original columns + trait)
        full_hits_list.append(hits)

        # Build simple format: position, trait, scaffold, SNP
        # Parse scaffold from SNP column: split at last underscore
        # e.g. 'SSCAFFOLD_1_648169' -> scaffold = 'SSCAFFOLD_1', pos_str = '648169'
        snp_split = hits["SNP"].astype(str).str.rsplit("_", n=1, expand=True)
        scaffold = snp_split[0]

        position = hits["BP"].astype(int)

        simple = pd.DataFrame({
            "position": position,
            "trait": trait,
            "scaffold": scaffold,
            "SNP": hits["SNP"].astype(str),
        })
        simple_hits_list.append(simple)

    # If no hits at all, stop here
    if not simple_hits_list:
        print("\nNo SNPs passed the MAF and P-value filters in any file.")
        print(f"No combined files written to {out_dir}.")
        return

    # Combine and sort simple hits
    combined_simple = pd.concat(simple_hits_list, ignore_index=True)
    combined_simple.sort_values(
        by=["trait", "scaffold", "position"], inplace=True
    )

    # Combine and sort full hits
    combined_full = pd.concat(full_hits_list, ignore_index=True)
    # A sensible sort: by trait, CHR, BP
    if "CHR" in combined_full.columns and "BP" in combined_full.columns:
        combined_full.sort_values(
            by=["trait", "CHR", "BP"], inplace=True
        )
    else:
        combined_full.sort_values(
            by=["trait"], inplace=True
        )

    # Write outputs in the new folder
    simple_path = out_dir / "combined_hits_simple.tsv"
    full_path = out_dir / "combined_hits_full.tsv"

    combined_simple.to_csv(simple_path, sep="\t", index=False)
    combined_full.to_csv(full_path, sep="\t", index=False)

    print(f"\nWrote simple combined hits to: {simple_path}")
    print(f"Wrote full combined hits (all columns + trait) to: {full_path}")
    print(f"Total SNPs in simple combined hits: {combined_simple.shape[0]}")
    print(f"Total SNPs in full combined hits:   {combined_full.shape[0]}")


if __name__ == "__main__":
    main()

