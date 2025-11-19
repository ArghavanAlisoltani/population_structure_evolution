#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, sys, csv
from collections import defaultdict
from typing import List, Tuple, Dict
import pandas as pd
import numpy as np

# ------------------------- helpers -------------------------

def ensure_outdir(path: str):
    path = os.path.expanduser(path)
    os.makedirs(path, exist_ok=True)
    return path

def load_table(path, default_sep="\t"):
    path = os.path.expanduser(path)
    try:
        return pd.read_csv(path, sep=default_sep, low_memory=False)
    except Exception:
        return pd.read_csv(path, delim_whitespace=True, low_memory=False)

def require(df: pd.DataFrame, cols: List[str], label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {label}: {missing}\nAvailable: {list(df.columns)}")

def sanitize(col: str) -> str:
    col = re.sub(r"\s+", "_", str(col))
    col = re.sub(r"[^0-9A-Za-z_]+", "_", col)
    return col

def q(x, p):
    try:
        return float(np.nanquantile(x, p))
    except Exception:
        return np.nan

def intervals_union_bp_df(df: pd.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pd.Series:
    """
    Compute union coverage per scaffold for a *moderate* sized df (TMRCA/mRNA).
    1-based inclusive lengths: end-start+1.
    """
    out = {}
    if df.empty: 
        return pd.Series(out, name="union_bp")
    for scf, sub in df.groupby(chrom_col, sort=False):
        s = pd.to_numeric(sub[start_col], errors="coerce").dropna().astype(np.int64).to_numpy()
        e = pd.to_numeric(sub[end_col],   errors="coerce").dropna().astype(np.int64).to_numpy()
        if len(s)==0 or len(e)==0:
            out[str(scf)] = 0
            continue
        arr = np.stack([s, e], axis=1)
        arr = arr[np.argsort(arr[:,0])]
        total = 0
        cur_s, cur_e = None, None
        for a, b in arr:
            if a > b: a, b = b, a
            if cur_s is None:
                cur_s, cur_e = int(a), int(b)
            else:
                if a <= cur_e + 1:
                    if b > cur_e: cur_e = int(b)
                else:
                    total += (cur_e - cur_s + 1)
                    cur_s, cur_e = int(a), int(b)
        if cur_s is not None:
            total += (cur_e - cur_s + 1)
        out[str(scf)] = int(total)
    return pd.Series(out, name="union_bp")

# ------------------------- TE hierarchy inference -------------------------

def infer_te_hierarchy(term: str) -> Tuple[str, str, str]:
    """
    Map a sequence_ontology string to (Class, Order, Superfamily).
    Literature-based heuristics; adjust if needed for your annotation scheme.
    """
    if term is None or (isinstance(term, float) and np.isnan(term)) or term == "":
        return ("Unclassified", "Unknown", "Unclassified")
    t = term.lower()
    def any_in(s, keys): return any(k in s for k in keys)

    # Class I
    if "retrotransposon" in t or any_in(t, ["line","sine","penelope","dirs","bel","pao"]):
        if "ltr" in t or any_in(t, ["gypsy","copia","bel","pao"]):
            order = "LTR"
            if "gypsy" in t: sf = "Gypsy"
            elif "copia" in t: sf = "Copia"
            elif any_in(t, ["bel","pao","bel_pao","bel-pao"]): sf = "Bel-Pao"
            else: sf = "LTR_other"
            return ("Class I (Retrotransposon)", order, sf)
        if "line" in t or any_in(t, ["l1","jockey","cr1","rte","tx1"]):
            order = "LINE"
            if "l1" in t: sf = "L1"
            elif "rte" in t: sf = "RTE"
            elif "cr1" in t: sf = "CR1"
            elif "jockey" in t: sf = "Jockey"
            else: sf = "LINE_other"
            return ("Class I (Retrotransposon)", order, sf)
        if "sine" in t:
            return ("Class I (Retrotransposon)", "SINE", "SINE_other")
        if "penelope" in t:
            return ("Class I (Retrotransposon)", "Penelope", "Penelope")
        if "dirs" in t or "tyrosine_recombinase" in t:
            return ("Class I (Retrotransposon)", "DIRS", "DIRS")
        return ("Class I (Retrotransposon)", "Unknown", "Retro_other")

    # Class II
    if "transposon" in t or any_in(t, ["tir","helitron","maverick","polinton","crypton","is_element"]):
        if "helitron" in t:
            return ("Class II (DNA transposon)", "Helitron", "Helitron")
        if "crypton" in t:
            return ("Class II (DNA transposon)", "Crypton", "Crypton")
        if "maverick" in t or "polinton" in t:
            return ("Class II (DNA transposon)", "Maverick/Polinton", "Maverick/Polinton")
        if "tir" in t or "transposon" in t or any_in(t, ["tc1","mariner","hat","cacta","enspm","mudr","mutator","pif","harbinger","piggybac"]):
            order = "TIR"
            if any_in(t, ["hat","hobo","ac"]): sf = "hAT"
            elif any_in(t, ["cacta","enspm","en-spm"]): sf = "CACTA/En-Spm"
            elif any_in(t, ["mutator","mudr"]): sf = "Mutator"
            elif any_in(t, ["harbinger","pif"]): sf = "PIF/Harbinger"
            elif any_in(t, ["tc1","mariner"]): sf = "Tc1/Mariner"
            elif "piggybac" in t: sf = "PiggyBac"
            else: sf = "TIR_other"
            return ("Class II (DNA transposon)", order, sf)
        return ("Class II (DNA transposon)", "Unknown", "DNA_other")

    if any_in(t, ["repeat_fragment","repeat","unknown","unclassified"]):
        return ("Unclassified", "Unknown", "Unclassified")
    return ("Unclassified", "Unknown", "Unclassified")

# ------------------------- TE streaming summary -------------------------

def summarize_te_streaming(te_path: str, te_chrom: str, te_class_col: str, te_start: str, te_end: str,
                           hierarchy_out: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, bool]:
    """
    Stream a very large TE TSV and compute:
      - per-scaffold totals (n_TEs_total, te_coverage_bp [exact if sorted], n_TE_classes, n_TE_superfamilies)
      - per-scaffold counts by 'sequence_ontology' (class) -> wide
      - per-scaffold counts by Superfamily -> wide
      - Write TE hierarchy mapping file
    Returns (te_summary_df, te_class_wide_df, te_sf_wide_df, coverage_exact_bool)
    """

    te_path = os.path.expanduser(te_path)
    # tracking
    n_total = 0
    class_per_scf_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    sf_per_scf_counts:    Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    class_set_per_scf:    Dict[str, set] = defaultdict(set)
    sf_set_per_scf:       Dict[str, set] = defaultdict(set)
    total_counts_per_scf: Dict[str, int] = defaultdict(int)

    # union coverage streaming (requires sorted by scaffold then start)
    coverage_exact = True
    last_scf = None
    last_start = -1
    cur_s = None
    cur_e = None
    te_union_bp: Dict[str, int] = defaultdict(int)

    # collect unique sequence_ontology terms to write mapping
    unique_terms = set()

    # read header
    with open(te_path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for req in (te_chrom, te_class_col, te_start, te_end):
            if req not in reader.fieldnames:
                raise ValueError(f"TE file missing column: {req}. Found: {reader.fieldnames}")

        for row in reader:
            scf = row[te_chrom]
            try:
                s = int(float(row[te_start]))
                e = int(float(row[te_end]))
            except Exception:
                continue
            if e < s: s, e = e, s

            cl = row[te_class_col] if row[te_class_col] is not None else ""
            unique_terms.add(cl)

            # infer hierarchy
            klass, order, sf = infer_te_hierarchy(cl)

            # counts
            total_counts_per_scf[scf] += 1
            class_per_scf_counts[scf][cl] += 1
            sf_per_scf_counts[scf][sf] += 1
            class_set_per_scf[scf].add(cl)
            sf_set_per_scf[scf].add(sf)
            n_total += 1

            # coverage (exact if sorted)
            if last_scf is None:
                last_scf, last_start = scf, s
                cur_s, cur_e = s, e
            else:
                if scf == last_scf:
                    if s < last_start:
                        coverage_exact = False
                    last_start = s
                    if coverage_exact:
                        if s <= cur_e + 1:
                            if e > cur_e: cur_e = e
                        else:
                            te_union_bp[scf] += (cur_e - cur_s + 1)
                            cur_s, cur_e = s, e
                else:
                    # flush previous scaffold
                    if coverage_exact and cur_s is not None:
                        te_union_bp[last_scf] += (cur_e - cur_s + 1)
                    # switch scaffold
                    last_scf, last_start = scf, s
                    cur_s, cur_e = s, e

        # flush last scaffold window
        if coverage_exact and last_scf is not None and cur_s is not None:
            te_union_bp[last_scf] += (cur_e - cur_s + 1)

    # If unsorted, approximate coverage = sum of lengths
    if not coverage_exact:
        te_union_bp = defaultdict(int)  # reset
        # second pass: sum lengths (still streaming and fast)
        with open(te_path, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                scf = row[te_chrom]
                try:
                    s = int(float(row[te_start])); e = int(float(row[te_end]))
                except Exception:
                    continue
                if e < s: s, e = e, s
                te_union_bp[scf] += (e - s + 1)

    # Build hierarchy mapping file
    mapped = sorted(unique_terms)
    map_rows = []
    for term in mapped:
        k, o, sf = infer_te_hierarchy(term)
        map_rows.append({"sequence_ontology": term, "Class": k, "Order": o, "Superfamily": sf})
    hier_df = pd.DataFrame(map_rows)
    hier_df.to_csv(os.path.expanduser(hierarchy_out), sep="\t", index=False)

    # Build per-scaffold summary
    scaffolds = sorted(set(total_counts_per_scf.keys()) | set(te_union_bp.keys()))
    recs = []
    for scf in scaffolds:
        n_tot = total_counts_per_scf.get(scf, 0)
        n_classes = len(class_set_per_scf.get(scf, set()))
        n_sfs     = len(sf_set_per_scf.get(scf, set()))
        cov_bp    = int(te_union_bp.get(scf, 0))
        recs.append({
            "scaffold": scf,
            "n_TEs_total": n_tot,
            "n_TE_classes": n_classes,
            "n_TE_superfamilies": n_sfs,
            "te_coverage_bp": cov_bp
        })
    te_summary_df = pd.DataFrame(recs)

    # Build wide class counts
    if scaffolds:
        # classes
        rows = []
        for scf, d in class_per_scf_counts.items():
            for cl, c in d.items():
                rows.append({"scaffold": scf, "TE_class": cl, "count": c})
        te_class_long = pd.DataFrame(rows)
        if not te_class_long.empty:
            te_class_long["TEclass_col"] = te_class_long["TE_class"].map(lambda x: "TEclass_"+sanitize(x))
            te_class_wide = (te_class_long
                             .pivot_table(index="scaffold", columns="TEclass_col", values="count",
                                          fill_value=0, aggfunc="sum")
                             .reset_index())
        else:
            te_class_wide = pd.DataFrame(columns=["scaffold"])
        # superfamilies
        rows = []
        for scf, d in sf_per_scf_counts.items():
            for sf, c in d.items():
                rows.append({"scaffold": scf, "Superfamily": sf, "count": c})
        te_sf_long = pd.DataFrame(rows)
        if not te_sf_long.empty:
            te_sf_long["TEsf_col"] = te_sf_long["Superfamily"].map(lambda x: "TEsf_"+sanitize(x))
            te_sf_wide = (te_sf_long
                          .pivot_table(index="scaffold", columns="TEsf_col", values="count",
                                       fill_value=0, aggfunc="sum")
                          .reset_index())
        else:
            te_sf_wide = pd.DataFrame(columns=["scaffold"])
    else:
        te_class_wide = pd.DataFrame(columns=["scaffold"])
        te_sf_wide    = pd.DataFrame(columns=["scaffold"])

    return te_summary_df, te_class_wide, te_sf_wide, coverage_exact

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser(
        description=("V5: Fast per-scaffold summaries for TMRCA, genes, and massive TE tables by streaming; "
                     "then merged into a final scaffold summary. Writes all outputs to a new folder.")
    )
    # Inputs
    ap.add_argument("--tmrca",
        required=True,
        help="TMRCA TSV (columns: CHROM, start_tmrca, end_tmrca, mean_tmrca, ...)")
    ap.add_argument("--mrna",
        required=True,
        help="mRNA TSV (columns: new_seqid, new_start, new_end, mrna_id)")
    ap.add_argument("--te",
        required=True,
        help="TE TSV (columns: seqid, sequence_ontology, start, end[, strand])")

    # Output folder
    ap.add_argument("--outdir",
        required=True,
        help="Output folder to write all intermediate and final summaries")

    # Column overrides
    ap.add_argument("--tmrca-chrom", default="CHROM")
    ap.add_argument("--tmrca-start", default="start_tmrca")
    ap.add_argument("--tmrca-end",   default="end_tmrca")
    ap.add_argument("--tmrca-mean",  default="mean_tmrca")

    ap.add_argument("--mrna-chrom",  default="new_seqid")
    ap.add_argument("--mrna-start",  default="new_start")
    ap.add_argument("--mrna-end",    default="new_end")
    ap.add_argument("--mrna-id",     default="mrna_id")

    ap.add_argument("--te-chrom",    default="seqid")
    ap.add_argument("--te-class",    default="sequence_ontology")
    ap.add_argument("--te-start",    default="start")
    ap.add_argument("--te-end",      default="end")

    # High-TMRCA threshold
    ap.add_argument("--top-percentile", type=float, default=25.0,
                    help="Top percentile for high TMRCA (global).")

    args = ap.parse_args()
    outdir = ensure_outdir(args.outdir)

    # ---------- TMRCA (independent) ----------
    tdf = load_table(args.tmrca)
    require(tdf, [args.tmrca_chrom, args.tmrca_start, args.tmrca_end, args.tmrca_mean], "TMRCA")
    tdf[args.tmrca_chrom] = tdf[args.tmrca_chrom].astype(str)
    for c in (args.tmrca_start, args.tmrca_end, args.tmrca_mean):
        tdf[c] = pd.to_numeric(tdf[c], errors="coerce")
    tdf["seg_len"] = (tdf[args.tmrca_end] - tdf[args.tmrca_start] + 1).clip(lower=0)

    top_thr = tdf[args.tmrca_mean].quantile(1.0 - float(args.top_percentile)/100.0)
    tdf["is_top"] = tdf[args.tmrca_mean] >= top_thr

    tmrca_stats = (
        tdf.groupby(args.tmrca_chrom)[args.tmrca_mean]
           .agg(n_tmrca_segments="size",
                tmrca_min="min",
                tmrca_q1=lambda s: q(s, 0.25),
                tmrca_median="median",
                tmrca_q3=lambda s: q(s, 0.75),
                tmrca_max="max",
                tmrca_mean="mean")
           .reset_index()
           .rename(columns={args.tmrca_chrom:"scaffold"})
    )
    tmrca_sd = (tdf.groupby(args.tmrca_chrom)[args.tmrca_mean]
                  .agg(lambda s: float(np.nanstd(s, ddof=1)) if s.notna().sum()>=2 else np.nan)
                  .reset_index().rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_mean:"tmrca_sd"}))
    tmrca_cov = intervals_union_bp_df(
        tdf.rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_start:"start", args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})
    seg_lens = (tdf.groupby(args.tmrca_chrom)["seg_len"]
                  .agg(median_seg_len="median", mean_seg_len="mean")
                  .reset_index().rename(columns={args.tmrca_chrom:"scaffold"}))
    top_count = (tdf.groupby([args.tmrca_chrom,"is_top"]).size()
                   .unstack(fill_value=0).reset_index()
                   .rename(columns={args.tmrca_chrom:"scaffold", True:"n_top_tmrca_segments", False:"n_non_top_tmrca_segments"}))
    top_cov = intervals_union_bp_df(
        tdf.loc[tdf["is_top"]].rename(columns={args.tmrca_chrom:"scaffold", args.tmrca_start:"start", args.tmrca_end:"end"}),
        "scaffold","start","end"
    ).rename("top_tmrca_coverage_bp").reset_index().rename(columns={"index":"scaffold"})

    tmrca_summary = (tmrca_stats.merge(tmrca_sd, on="scaffold", how="left")
                                 .merge(tmrca_cov, on="scaffold", how="left")
                                 .merge(seg_lens, on="scaffold", how="left")
                                 .merge(top_count, on="scaffold", how="left")
                                 .merge(top_cov, on="scaffold", how="left"))
    tmrca_summary["tmrca_cv"] = np.where(tmrca_summary["tmrca_mean"].notna() & (tmrca_summary["tmrca_mean"]!=0),
                                         tmrca_summary["tmrca_sd"]/tmrca_summary["tmrca_mean"], np.nan)
    tmrca_summary_path = os.path.join(outdir, "tmrca_summary.tsv")
    tmrca_summary.to_csv(tmrca_summary_path, sep="\t", index=False)

    # ---------- mRNA (independent) ----------
    gdf = load_table(args.mrna)
    require(gdf, [args.mrna_chrom, args.mrna_start, args.mrna_end], "mRNA")
    gdf[args.mrna_chrom] = gdf[args.mrna_chrom].astype(str)
    for c in (args.mrna_start, args.mrna_end):
        gdf[c] = pd.to_numeric(gdf[c], errors="coerce")
    if args.mrna_id not in gdf.columns:
        gdf[args.mrna_id] = np.nan
    else:
        gdf[args.mrna_id] = gdf[args.mrna_id].astype(str)

    gene_counts = (gdf.groupby(args.mrna_chrom)[args.mrna_id]
                     .nunique(dropna=True).reset_index()
                     .rename(columns={args.mrna_chrom:"scaffold", args.mrna_id:"n_genes"}))
    gene_cov = intervals_union_bp_df(
        gdf.rename(columns={args.mrna_chrom:"scaffold", args.mrna_start:"start", args.mrna_end:"end"}),
        "scaffold","start","end"
    ).rename("gene_coverage_bp").reset_index().rename(columns={"index":"scaffold"})
    mrna_summary = gene_counts.merge(gene_cov, on="scaffold", how="outer")
    mrna_summary_path = os.path.join(outdir, "mrna_summary.tsv")
    mrna_summary.to_csv(mrna_summary_path, sep="\t", index=False)

    # ---------- TE (independent, streaming) ----------
    te_hierarchy_path = os.path.join(outdir, "TE_hierarchy_v5.tsv")
    te_summary_df, te_class_wide, te_sf_wide, te_cov_exact = summarize_te_streaming(
        args.te, args.te_chrom, args.te_class, args.te_start, args.te_end, te_hierarchy_path
    )
    te_summary_path = os.path.join(outdir, "te_summary.tsv")
    te_class_wide_path = os.path.join(outdir, "te_class_wide.tsv")
    te_sf_wide_path = os.path.join(outdir, "te_superfamily_wide.tsv")
    te_summary_df.to_csv(te_summary_path, sep="\t", index=False)
    te_class_wide.to_csv(te_class_wide_path, sep="\t", index=False)
    te_sf_wide.to_csv(te_sf_wide_path, sep="\t", index=False)

    # ---------- Merge by scaffold (outer) ----------
    scaff_all = pd.Index(sorted(set(tmrca_summary["scaffold"].dropna().astype(str).unique())
                             | set(mrna_summary["scaffold"].dropna().astype(str).unique())
                             | set(te_summary_df["scaffold"].dropna().astype(str).unique())), dtype="object")
    base = pd.DataFrame({"scaffold": scaff_all})

    out = (base
           .merge(tmrca_summary, on="scaffold", how="left")
           .merge(mrna_summary,  on="scaffold", how="left")
           .merge(te_summary_df, on="scaffold", how="left")
           .merge(te_sf_wide,    on="scaffold", how="left")
           .merge(te_class_wide, on="scaffold", how="left"))

    # Safe integer fills
    for c in ["n_tmrca_segments","tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
              "top_tmrca_coverage_bp","n_genes","gene_coverage_bp","n_TEs_total",
              "n_TE_classes","n_TE_superfamilies","te_coverage_bp"]:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)

    # Derived fractions
    def frac(a, b): return np.where(b>0, a/b, np.nan)
    if "tmrca_coverage_bp" in out.columns and "top_tmrca_coverage_bp" in out.columns:
        out["fraction_top_tmrca_bp"] = frac(out["top_tmrca_coverage_bp"], out["tmrca_coverage_bp"])
    if "gene_coverage_bp" in out.columns and "te_coverage_bp" in out.columns:
        out["te_gene_coverage_ratio"] = frac(out["te_coverage_bp"], out["gene_coverage_bp"])
        cov_sum = out["gene_coverage_bp"] + out["te_coverage_bp"]
        out["te_fraction_of_covered_bp"] = frac(out["te_coverage_bp"], cov_sum)

    # ---------- Robust scaffold_1 pooling (works with 1, 1a, 1b in any combo) ----------
    have = set(out["scaffold"].astype(str))
    targets = [n for n in ["scaffold_1","scaffold_1a","scaffold_1b"] if n in have]
    if targets:
        sub = out[out["scaffold"].isin(targets)].copy()
        pooled = {"scaffold": "scaffold_1"}
        # sum integer columns
        int_cols = [c for c in out.columns if re.search(r"^n_|_bp$", c)]
        for c in int_cols:
            pooled[c] = int(sub[c].fillna(0).sum())
        # pooled TMRCA stats recomputed from raw TMRCA table for accuracy
        mask = tdf[args.tmrca_chrom].isin(targets)
        tm_vals = tdf.loc[mask, args.tmrca_mean].dropna().to_numpy()
        seg_vals = tdf.loc[mask, "seg_len"].dropna().to_numpy()
        pooled.update({
            "tmrca_min": float(np.nanmin(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_q1":  q(tm_vals, 0.25) if tm_vals.size>0 else np.nan,
            "tmrca_median": float(np.nanmedian(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_q3":  q(tm_vals, 0.75) if tm_vals.size>0 else np.nan,
            "tmrca_max": float(np.nanmax(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_mean": float(np.nanmean(tm_vals)) if tm_vals.size>0 else np.nan,
            "tmrca_sd": float(np.nanstd(tm_vals, ddof=1)) if tm_vals.size>=2 else np.nan,
            "tmrca_cv": (float(np.nanstd(tm_vals, ddof=1))/float(np.nanmean(tm_vals))
                         if tm_vals.size>=2 and np.nanmean(tm_vals)!=0 else np.nan),
            "median_seg_len": float(np.nanmedian(seg_vals)) if seg_vals.size>0 else np.nan,
            "mean_seg_len":   float(np.nanmean(seg_vals))   if seg_vals.size>0 else np.nan,
        })
        # derived
        if pooled.get("tmrca_coverage_bp",0) > 0 and "top_tmrca_coverage_bp" in pooled:
            pooled["fraction_top_tmrca_bp"] = pooled["top_tmrca_coverage_bp"]/pooled["tmrca_coverage_bp"]
        if pooled.get("gene_coverage_bp",0) > 0:
            pooled["te_gene_coverage_ratio"] = pooled.get("te_coverage_bp",0) / pooled["gene_coverage_bp"]
        denom = pooled.get("gene_coverage_bp",0) + pooled.get("te_coverage_bp",0)
        pooled["te_fraction_of_covered_bp"] = (pooled.get("te_coverage_bp",0)/denom) if denom>0 else np.nan

        # sum all wide TE columns too
        wide_cols = [c for c in out.columns if c.startswith("TEsf_") or c.startswith("TEclass_")]
        for c in wide_cols:
            pooled[c] = int(sub[c].fillna(0).sum()) if c in sub.columns else 0

        # drop any existing 'scaffold_1' then append pooled
        out = out[out["scaffold"] != "scaffold_1"]
        out = pd.concat([out, pd.DataFrame([pooled])], ignore_index=True)

    # ---------- Order columns & write ----------
    fixed = ["scaffold",
             "n_tmrca_segments","tmrca_min","tmrca_q1","tmrca_median","tmrca_q3","tmrca_max","tmrca_mean",
             "tmrca_sd","tmrca_cv","median_seg_len","mean_seg_len",
             "tmrca_coverage_bp","n_top_tmrca_segments","n_non_top_tmrca_segments",
             "top_tmrca_coverage_bp","fraction_top_tmrca_bp",
             "n_genes","gene_coverage_bp",
             "n_TEs_total","n_TE_classes","n_TE_superfamilies","te_coverage_bp",
             "te_gene_coverage_ratio","te_fraction_of_covered_bp"]
    te_sf_cols    = sorted([c for c in out.columns if c.startswith("TEsf_")])
    te_class_cols = sorted([c for c in out.columns if c.startswith("TEclass_")])
    other_cols    = [c for c in out.columns if c not in fixed + te_sf_cols + te_class_cols]
    final_cols = fixed + te_sf_cols + te_class_cols + other_cols
    out = out.reindex(columns=final_cols)

    final_path = os.path.join(outdir, "scaffold_summary_v5.tsv")
    out.to_csv(final_path, sep="\t", index=False)

    # meta
    with open(os.path.join(outdir, "README_v5.txt"), "w") as fh:
        fh.write("V5 per-scaffold summaries\n")
        fh.write(f"High-TMRCA threshold (mean_tmrca >=): {top_thr}\n")
        fh.write(f"TE coverage exact: {te_cov_exact}\n")
        if not te_cov_exact:
            fh.write("NOTE: TE file not sorted by seqid,start. Exact union coverage disabled.\n")
            fh.write("      To enable exact union coverage, pre-sort TE file, e.g.:\n")
            fh.write("      LC_ALL=C sort -t$'\\t' -k1,1 -k3,3n your_TE.tsv > your_TE.sorted.tsv\n")

    print(f"[OK] Wrote: {tmrca_summary_path}")
    print(f"[OK] Wrote: {mrna_summary_path}")
    print(f"[OK] Wrote: {te_summary_path}")
    print(f"[OK] Wrote: {te_class_wide_path}")
    print(f"[OK] Wrote: {te_sf_wide_path}")
    print(f"[OK] Wrote: {final_path}")
    print(f"[OK] TE hierarchy: {os.path.join(outdir, 'TE_hierarchy_v5.tsv')}")
    if not te_cov_exact:
        print("** TE coverage used approximate sum(lengths) because TE file was not sorted by seqid,start.")
        print("   Pre-sort to get exact union coverage (see README_v5.txt).")

if __name__ == "__main__":
    main()

