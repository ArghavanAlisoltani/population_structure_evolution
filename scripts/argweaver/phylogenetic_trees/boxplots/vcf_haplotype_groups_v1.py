#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
from collections import Counter, defaultdict
'''
Run example
python3 vcf_haplotype_groups.py \
  --vcf merged.vcf.gz \
  --scaffold scaffold_4 \
  --positions 983057685,983057688,983057694,983057707,983057714,983057718,983057724 \
  --outdir hap_scaffold4

'''
DEFAULT_SITES = [
    ("scaffold_4", 983057685),
    ("scaffold_4", 983057688),
    ("scaffold_4", 983057694),
    ("scaffold_4", 983057707),
    ("scaffold_4", 983057714),
    ("scaffold_4", 983057718),
    ("scaffold_4", 983057724),
]

def run(cmd, check=True, text=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    if check and p.returncode != 0:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(cmd)}\n\nSTDERR:\n{p.stderr}"
        )
    return p.stdout

def parse_sites(args):
    sites = []
    if args.sites_file:
        # Accept: "CHROM  POS" (tab/space), optionally with header lines starting with '#'
        with open(args.sites_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"Bad sites_file line (need CHROM POS): {line}")
                chrom = parts[0]
                pos = int(parts[1])
                sites.append((chrom, pos))
    elif args.positions:
        if not args.scaffold:
            raise ValueError("--positions requires --scaffold (or use --sites_file).")
        pos_list = [int(x.strip()) for x in args.positions.split(",") if x.strip()]
        sites = [(args.scaffold, p) for p in pos_list]
    else:
        sites = list(DEFAULT_SITES)

    if len(sites) == 0:
        raise ValueError("No sites provided.")
    return sites

def allele_to_base(ref, alt, allele_idx):
    """
    allele_idx: int (0=REF, 1..=ALT index), or None for missing.
    Returns a single-character base if possible; otherwise 'N'.
    """
    if allele_idx is None:
        return "N"

    if allele_idx == 0:
        a = ref
    else:
        alts = alt.split(",") if alt else []
        j = allele_idx - 1
        if j < 0 or j >= len(alts):
            return "N"
        a = alts[j]

    # This pipeline is intended for SNP hap strings (1 bp per site).
    # If indel/multibase allele occurs, return 'N' to keep fixed length.
    if len(a) != 1:
        return "N"
    return a

def parse_gt(gt):
    """
    Return (a1, a2, phased_flag) where a1/a2 are ints or None if missing.
    Handles '0|1', '0/1', '.|.', './.' etc.
    """
    gt = gt.strip()
    if gt in (".", "./.", ".|.", "./", ".|"):
        return (None, None, False)

    phased = "|" in gt
    sep = "|" if phased else ("/" if "/" in gt else None)
    if sep is None:
        # Unexpected GT; treat as missing
        return (None, None, False)

    parts = gt.split(sep)
    if len(parts) != 2:
        return (None, None, phased)

    def to_int(x):
        x = x.strip()
        if x == "." or x == "":
            return None
        try:
            return int(x)
        except ValueError:
            return None

    return (to_int(parts[0]), to_int(parts[1]), phased)

def main():
    ap = argparse.ArgumentParser(
        description="Build per-sample phased haplotype strings at selected SNP positions from a VCF, and group identical haplotypes."
    )
    ap.add_argument("--vcf", required=True, help="Input VCF (.vcf.gz) indexed with .tbi/.csi")
    ap.add_argument("--outdir", default="hap_out", help="Output directory")
    ap.add_argument("--scaffold", default=None, help="Scaffold/contig ID (used with --positions)")
    ap.add_argument("--positions", default=None, help="Comma-separated 1-based positions (used with --scaffold)")
    ap.add_argument("--sites_file", default=None, help="File with CHROM POS per line (recommended for multi-scaffold)")
    ap.add_argument("--bcftools", default="bcftools", help="Path to bcftools executable")
    ap.add_argument("--missing_char", default="N", help="Character for missing/invalid allele (default N)")
    ap.add_argument("--drop_any_missing", action="store_true",
                    help="If set, drop any haplotype containing missing_char from grouping (still output per-sample fasta).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Sites (keep user order!)
    sites = parse_sites(args)
    site_to_index = {(c, p): i for i, (c, p) in enumerate(sites)}
    m = len(sites)

    # 2) Sample list (order must match bcftools query output)
    try:
        sample_list = run([args.bcftools, "query", "-l", args.vcf]).strip().splitlines()
    except Exception as e:
        print(f"ERROR: failed to get sample list using bcftools query -l\n{e}", file=sys.stderr)
        sys.exit(2)

    if not sample_list:
        print("ERROR: no samples found in VCF.", file=sys.stderr)
        sys.exit(2)

    n = len(sample_list)

    # 3) Prepare per-sample hap arrays (fixed length m)
    hap1 = {s: [args.missing_char] * m for s in sample_list}
    hap2 = {s: [args.missing_char] * m for s in sample_list}

    # 4) Query only the requested sites
    regions_path = os.path.join(args.outdir, "sites.regions.txt")
    with open(regions_path, "w") as f:
        for chrom, pos in sites:
            f.write(f"{chrom}:{pos}-{pos}\n")

    fmt = r"%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"
    try:
        raw = run([args.bcftools, "query", "-R", regions_path, "-f", fmt, args.vcf])
    except Exception as e:
        print(f"ERROR: bcftools query failed.\n{e}", file=sys.stderr)
        sys.exit(2)

    found = set()
    loci_info = { (c,p): None for (c,p) in sites }

    for line in raw.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[2]
        alt = parts[3]
        gts = parts[4:]

        key = (chrom, pos)
        if key not in site_to_index:
            continue

        idx = site_to_index[key]
        found.add(key)
        loci_info[key] = (ref, alt)

        # Ensure genotype column count matches sample list; if not, skip safely
        if len(gts) != n:
            print(
                f"WARNING: site {chrom}:{pos} has {len(gts)} GT columns but VCF has {n} samples; skipping this site.",
                file=sys.stderr,
            )
            loci_info[key] = None
            continue

        # Fill per-sample hap bases
        for i, s in enumerate(sample_list):
            a1, a2, _phased = parse_gt(gts[i])
            b1 = allele_to_base(ref, alt, a1)
            b2 = allele_to_base(ref, alt, a2)
            if b1 != args.missing_char:
                hap1[s][idx] = b1
            else:
                hap1[s][idx] = args.missing_char
            if b2 != args.missing_char:
                hap2[s][idx] = b2
            else:
                hap2[s][idx] = args.missing_char

    # 5) Write loci_used.tsv (in requested order) and report missing sites
    loci_tsv = os.path.join(args.outdir, "loci_used.tsv")
    missing_sites = []
    with open(loci_tsv, "w") as out:
        out.write("index\tCHROM\tPOS\tREF\tALT\tstatus\n")
        for i, (chrom, pos) in enumerate(sites):
            info = loci_info.get((chrom, pos))
            if info is None:
                out.write(f"{i+1}\t{chrom}\t{pos}\t.\t.\tMISSING_IN_VCF_OR_BAD_COLUMNS\n")
                missing_sites.append((chrom, pos))
            else:
                ref, alt = info
                out.write(f"{i+1}\t{chrom}\t{pos}\t{ref}\t{alt}\tOK\n")

    # 6) Write per-sample FASTA (two haplotypes each)
    fasta_per_sample = os.path.join(args.outdir, "haplotypes_per_sample.fasta")
    with open(fasta_per_sample, "w") as out:
        for s in sample_list:
            s1 = "".join(hap1[s])
            s2 = "".join(hap2[s])
            out.write(f">{s}_1\n{s1}\n")
            out.write(f">{s}_2\n{s2}\n")

    # 7) Group haplotypes by identical sequence
    hap_records = []  # (sample, hap_idx, seq)
    for s in sample_list:
        s1 = "".join(hap1[s])
        s2 = "".join(hap2[s])
        hap_records.append((s, 1, s1))
        hap_records.append((s, 2, s2))

    # Optionally drop any haplotype with missing_char from grouping
    if args.drop_any_missing:
        hap_records_for_group = [r for r in hap_records if args.missing_char not in r[2]]
    else:
        hap_records_for_group = hap_records

    seq_counts = Counter([seq for (_s, _h, seq) in hap_records_for_group])
    seq_samples = defaultdict(set)
    for s, _h, seq in hap_records_for_group:
        seq_samples[seq].add(s)

    # Stable ID assignment: by count desc, then sequence lexicographically
    seq_sorted = sorted(seq_counts.items(), key=lambda x: (-x[1], x[0]))
    seq_to_hid = {}
    for k, (seq, _cnt) in enumerate(seq_sorted, start=1):
        seq_to_hid[seq] = f"H{k}"

    # 8) Write haplotype_groups.tsv
    groups_tsv = os.path.join(args.outdir, "haplotype_groups.tsv")
    denom = sum(seq_counts.values()) if seq_counts else 0
    with open(groups_tsv, "w") as out:
        out.write("hap_id\tsequence\tn_haplotypes\tn_samples\tfrequency\n")
        for seq, cnt in seq_sorted:
            hid = seq_to_hid[seq]
            ns = len(seq_samples[seq])
            freq = (cnt / denom) if denom else 0.0
            out.write(f"{hid}\t{seq}\t{cnt}\t{ns}\t{freq:.6f}\n")

    # 9) Write unique haplotypes FASTA
    fasta_unique = os.path.join(args.outdir, "haplotypes_unique.fasta")
    with open(fasta_unique, "w") as out:
        for seq, cnt in seq_sorted:
            hid = seq_to_hid[seq]
            out.write(f">{hid} count={cnt}\n{seq}\n")

    # 10) Write per-sample hap calls table
    calls_tsv = os.path.join(args.outdir, "sample_haplotype_calls.tsv")
    with open(calls_tsv, "w") as out:
        out.write("sample\thap1_id\thap2_id\thap1_seq\thap2_seq\tgenotype_label\n")
        for s in sample_list:
            s1 = "".join(hap1[s])
            s2 = "".join(hap2[s])
