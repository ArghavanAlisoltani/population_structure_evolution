#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import pandas as pd


def get_trait_name(filename: str) -> str:
    """
    Extract trait name from a file name like 'C13_sumstat.txt'
    -> 'C13'
    """
    base = Path(filename).name
    if base.endswith("_sumstat.txt"):
        return base[:-len("_sumstat.txt")]
    elif base.endswith(".txt"):
        return base[:-4]
    else:
        # fallback
        return base


def main():
    parser = argparse.ArgumentParser(
        description="Combine significant SNPs from multiple GWAS sumstat files."
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="Input summary statistic files (e.g. C13_sumstat.txt DBH30_sumstat.txt ...)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Maximum FDR to call a SNP significant (default: 0.10)",
    )
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help="Minimum MAF to keep a SNP (default: 0.01)",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="combined_hits.tsv",
        help="Output file name (default: combined_hits.tsv)",
    )

    args = parser.parse_args()

    all_hits = []

    for fname in args.files:
        trait = get_trait_name(fname)

        print(f"Processing {fname} (trait = {trait})")

        # Read the sumstat file; works with space or tab delimiters
        df = pd.read_csv(fname, delim_whitespace=True)

        # Ensure required columns exist
        required_cols = {"SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "Z", "P", "N", "fdr"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"File {fname} is missing columns: {missing}")

        # Drop rows with missing MAF or fdr
        df = df[df["MAF"].notna() & df["fdr"].notna()].copy()

        # Apply filters: FDR and MAF
        hit_mask = (df["fdr"] <= args.fdr_threshold) & (df["MAF"] >= args.maf_threshold)
        hits = df.loc[hit_mask].copy()

        if hits.empty:
            print(f"  No SNPs passed filters in {fname}")
            continue

        # Parse scaffold from SNP column by splitting at the last underscore
        # e.g. 'SSCAFFOLD_1_648169' -> scaffold = 'SSCAFFOLD_1', pos_str = '648169'
        snp_split = hits["SNP"].str.rsplit("_", n=1, expand=True)
        scaffold = snp_split[0]
        # position we will take from BP, but you could also compare with snp_split[1]
        position = hits["BP"].astype(int)

        hits["position"] = position
        hits["scaffold"] = scaffold
        hits["trait"] = trait

        # Select only desired columns
        subset = hits[["position", "trait", "scaffold", "SNP"]].copy()
        all_hits.append(subset)

    if not all_hits:
        print("No SNPs passed the filters in any file. No output written.")
        return

    combined = pd.concat(all_hits, ignore_index=True)

    # Sort by trait / scaffold / position to keep it tidy
    combined.sort_values(by=["trait", "scaffold", "position"], inplace=True)

    # Write output
    combined.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote combined hits to: {args.out}")
    print(f"Total SNPs: {combined.shape[0]}")


if __name__ == "__main__":
    main()

