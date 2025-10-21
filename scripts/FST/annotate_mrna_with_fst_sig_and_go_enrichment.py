#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Annotate mRNAs with per-comparison FST significance (overlap with Poisson-q significant windows),
split into sig vs not_sig files, and run in-house GO enrichment (sig vs not_sig) per comparison.

Inputs
------
1) mRNA annotation TSV (after your scaffold_1 a/b fix), must have:
   - new_seqid, new_start, new_end
   - a gene ID column (default: mrna_id)
   - a GO column (default: ipr_go_terms) containing pipe-separated GO IDs (e.g. 'GO:0004601|GO:0006979|...')

2) FST windows TSV:
   - CHROM, WIN_START, WIN_END
   - q_poi_3vs4, q_poi_3vs5, q_poi_4vs5

Outputs (per -o OUTDIR)
-----------------------
- mrna_with_fst_overlap.tsv             # original mRNA + FST_3vs4/FST_3vs5/FST_4vs5 ('sig'/'not_sig')
- <comp>/
    mrna_sig_<comp>.tsv                 # full annotation rows where FST_<comp> == 'sig'
    mrna_not_sig_<comp>.tsv             # full annotation rows where FST_<comp> == 'not_sig'
    mrna_sig_<comp>.GO.tsv              # ONE-COLUMN (GO only) table from sig set
    mrna_not_sig_<comp>.GO.tsv          # ONE-COLUMN (GO only) table from not_sig set
    go_enrichment_<comp>.tsv            # Fisher (sig vs not_sig), BH-FDR q-values

Notes
-----
- Overlap logic: a gene is 'sig' for a comparison if it overlaps ANY window with q_poi_<comp> <= threshold.
  Otherwise it's 'not_sig'.
- GO enrichment is presence/absence per gene (not weighted by multiple GO repeats on the same gene).


Run example
-----
python annotate_mrna_with_fst_sig_and_go_enrichment.py \
  -a /Users/aria/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/1ab_mRNA_ID_merged_interproscan.txt \
  -w fst_windows_with_TEcounts__augmented.tsv \
  -o mrna_fst_sig_go_v1 \
  --q_threshold 0.05 \
  --go_col ipr_go_terms \
  --id_col mrna_id \
  --ann_chr new_seqid --ann_start new_start --ann_end new_end \
  --win_chr CHROM --win_start WIN_START --win_end WIN_END

"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

COMPS = ["3vs4", "3vs5", "4vs5"]

# ---------- utils ----------

def to_num(a):
    return pd.to_numeric(a, errors="coerce")

def clean_coord_series(s):
    # robust numeric coercion (handles stray text, NA tokens, commas)
    s = s.astype(str).str.strip()
    s = s.replace({"": np.nan, "NA": np.nan, "NaN": np.nan, "None": np.nan, ".": np.nan})
    s = s.str.replace(",", "", regex=False)
    # grab first integer substring if needed
    s = s.str.extract(r"(-?\d+)", expand=False)
    return pd.to_numeric(s, errors="coerce")

def interval_overlap_annotate_genes(mrna_df, win_df, q_col, thresh,
                                    mrna_chr="new_seqid", mrna_start="new_start", mrna_end="new_end",
                                    win_chr="CHROM", win_start="WIN_START", win_end="WIN_END"):
    """
    Return a pandas Series with values 'sig' or 'not_sig' for each mRNA,
    based on overlap with windows where q_col <= thresh.
    """
    label = np.array(["not_sig"] * len(mrna_df), dtype=object)

    # filter significant windows for this comparison
    sig_win = win_df.loc[to_num(win_df[q_col]) <= thresh, [win_chr, win_start, win_end]].copy()
    if sig_win.empty:
        return pd.Series(label, index=mrna_df.index, name="tmp")

    # per-chromosome sweep
    # pre-group windows
    for chrom, sub_mrna in mrna_df.groupby(mrna_chr, sort=False):
        sub_win = sig_win.loc[sig_win[win_chr] == chrom]
        if sub_win.empty:
            continue

        # windows arrays
        ws = clean_coord_series(sub_win[win_start]).astype("Int64").dropna()
        we = clean_coord_series(sub_win[win_end]).astype("Int64").dropna()
        if ws.empty or we.empty:
            continue
        w = pd.DataFrame({"s": ws.astype(np.int64).values, "e": we.astype(np.int64).values})
        w = w.sort_values("s").to_numpy()  # shape (n,2)

        # gene arrays with indices
        gs = clean_coord_series(sub_mrna[mrna_start]).astype("Int64")
        ge = clean_coord_series(sub_mrna[mrna_end]).astype("Int64")
        valid = (~gs.isna()) & (~ge.isna())
        if not valid.any():
            continue
        idx = sub_mrna.index[valid].to_numpy()
        gs = gs[valid].astype(np.int64).to_numpy()
        ge = ge[valid].astype(np.int64).to_numpy()

        # two-pointer approach: for each gene, advance window pointer until e >= gene.start
        j = 0
        nW = w.shape[0]
        for k in range(len(idx)):
            S = int(gs[k]); E = int(ge[k])
            # move window pointer forward while window ends before gene starts
            while j < nW and w[j,1] < S:
                j += 1
            jj = j
            overlapped = False
            # check from current pointer while window starts <= gene.end
            while jj < nW and w[jj,0] <= E:
                # overlap exists if w[jj].e >= S and w[jj].s <= E  (both guaranteed by loops)
                overlapped = True
                break
            if overlapped:
                label[np.where(mrna_df.index == idx[k])[0][0]] = "sig"

    return pd.Series(label, index=mrna_df.index, name="tmp")

def extract_go_set_per_gene(go_series, go_col_name="ipr_go_terms"):
    """
    Returns:
      - dict: gene_index -> set of GO IDs (GO:XXXXXXX)
      - set: all GO IDs observed
    """
    all_terms = set()
    per_gene = {}
    vals = go_series.fillna("").astype(str).values
    for i, s in enumerate(vals):
        if not s or s == "-" or s.strip() == "":
            per_gene[i] = set()
            continue
        parts = [p.strip() for p in s.split("|") if p.strip()]
        gos = {p for p in parts if p.startswith("GO:")}
        per_gene[i] = gos
        all_terms.update(gos)
    return per_gene, all_terms

def bh_fdr(pvals):
    """Benjamini–Hochberg FDR (returns array of q-values in same order)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n+1)
    q = p * n / ranks
    # monotone
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q[ranks.argsort()] = q_sorted  # map back
    return np.clip(q, 0, 1)

def go_enrichment(sig_df, bg_df, go_col="ipr_go_terms"):
    """
    Fisher exact (alternative='greater') on presence/absence per gene for each GO term.
    Returns a dataframe with p and BH q.
    """
    # Build GO presence sets
    sig_map, sig_all = extract_go_set_per_gene(sig_df[go_col], go_col)
    bg_map, bg_all = extract_go_set_per_gene(bg_df[go_col], go_col)
    universe = sorted(sig_all.union(bg_all))
    n_sig = len(sig_map)
    n_bg  = len(bg_map)
    rows = []
    for go in universe:
        a = sum(1 for i in sig_map if go in sig_map[i])
        c = sum(1 for i in bg_map  if go in bg_map[i])
        b = n_sig - a
        d = n_bg  - c
        # Skip degenerate (no information)
        if (a + c) == 0:
            continue
        table = [[a, b], [c, d]]
        try:
            odds, p = fisher_exact(table, alternative="greater")
        except Exception:
            odds, p = np.nan, 1.0
        rows.append((go, a, n_sig, c, n_bg, odds, p))
    if not rows:
        return pd.DataFrame(columns=["GO","sig_with_GO","sig_total","bg_with_GO","bg_total","odds_ratio","p_value","q_value"])
    out = pd.DataFrame(rows, columns=["GO","sig_with_GO","sig_total","bg_with_GO","bg_total","odds_ratio","p_value"])
    out["q_value"] = bh_fdr(out["p_value"].values)
    out = out.sort_values(["q_value","p_value","odds_ratio"], ascending=[True, True, False])
    return out

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Annotate mRNAs with FST significant overlap per comparison and run GO enrichment (sig vs not_sig).")
    ap.add_argument("-a", "--annotation", required=True, help="mRNA annotation TSV (with new_seqid/new_start/new_end)")
    ap.add_argument("-w", "--windows",     required=True, help="FST windows TSV")
    ap.add_argument("-o", "--outdir",      required=True, help="Output directory")
    ap.add_argument("--q_threshold", type=float, default=0.05, help="Poisson q threshold to call windows significant (default 0.05)")
    ap.add_argument("--go_col", default="ipr_go_terms", help="GO column name in annotation (default: ipr_go_terms)")
    ap.add_argument("--id_col", default="mrna_id", help="mRNA ID column name in annotation (default: mrna_id)")
    # column names (override if your files differ)
    ap.add_argument("--ann_chr", default="new_seqid")
    ap.add_argument("--ann_start", default="new_start")
    ap.add_argument("--ann_end", default="new_end")
    ap.add_argument("--win_chr", default="CHROM")
    ap.add_argument("--win_start", default="WIN_START")
    ap.add_argument("--win_end", default="WIN_END")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Load
    mrna = pd.read_csv(args.annotation, sep="\t", low_memory=False)
    win  = pd.read_csv(args.windows, sep="\t", low_memory=False)

    # Sanity
    for c in [args.ann_chr, args.ann_start, args.ann_end, args.id_col, args.go_col]:
        if c not in mrna.columns:
            raise SystemExit(f"Annotation missing column: {c}")
    for c in [args.win_chr, args.win_start, args.win_end, "q_poi_3vs4", "q_poi_3vs5", "q_poi_4vs5"]:
        if c not in win.columns:
            raise SystemExit(f"Windows missing column: {c}")

    # Per-comparison gene labels
    for comp in COMPS:
        colname = f"FST_{comp}"
        qcol    = f"q_poi_{comp}"
        lab = interval_overlap_annotate_genes(
            mrna_df=mrna, win_df=win, q_col=qcol, thresh=args.q_threshold,
            mrna_chr=args.ann_chr, mrna_start=args.ann_start, mrna_end=args.ann_end,
            win_chr=args.win_chr, win_start=args.win_start, win_end=args.win_end
        )
        mrna[colname] = lab.values

    # Write combined annotation with labels
    combined_path = Path(outdir, "mrna_with_fst_overlap.tsv")
    mrna.to_csv(combined_path, sep="\t", index=False)

    # Per-comparison splits + GO-only + enrichment
    for comp in COMPS:
        subdir = Path(outdir, comp); subdir.mkdir(parents=True, exist_ok=True)
        colname = f"FST_{comp}"

        sig_df = mrna.loc[mrna[colname] == "sig"].copy()
        bg_df  = mrna.loc[mrna[colname] == "not_sig"].copy()

        sig_df.to_csv(subdir / f"mrna_sig_{comp}.tsv", sep="\t", index=False)
        bg_df.to_csv(subdir / f"mrna_not_sig_{comp}.tsv", sep="\t", index=False)

        # GO-only (one column)
        sig_go = sig_df[[args.go_col]].rename(columns={args.go_col: "GO"})
        bg_go  = bg_df[[args.go_col]].rename(columns={args.go_col: "GO"})
        sig_go.to_csv(subdir / f"mrna_sig_{comp}.GO.tsv", sep="\t", index=False)
        bg_go.to_csv(subdir / f"mrna_not_sig_{comp}.GO.tsv", sep="\t", index=False)

        # Enrichment
        enr = go_enrichment(sig_df, bg_df, go_col=args.go_col)
        enr.to_csv(subdir / f"go_enrichment_{comp}.tsv", sep="\t", index=False)

    # README
    with open(Path(outdir, "README.txt"), "w") as fh:
        fh.write(f"""mRNA FST overlap + GO enrichment

Inputs:
  annotation: {args.annotation}
  windows:    {args.windows}

Params:
  q_threshold: {args.q_threshold}
  GO column:   {args.go_col}
  ID column:   {args.id_col}
  Overlap columns: mRNA({args.ann_chr},{args.ann_start},{args.ann_end}) vs Windows({args.win_chr},{args.win_start},{args.win_end})

Outputs:
  mrna_with_fst_overlap.tsv   # adds FST_3vs4, FST_3vs5, FST_4vs5 = 'sig'/'not_sig'
  <comp>/mrna_sig_<comp>.tsv and mrna_not_sig_<comp>.tsv
  <comp>/mrna_sig_<comp>.GO.tsv and mrna_not_sig_<comp>.GO.tsv  # GO column only
  <comp>/go_enrichment_<comp>.tsv  # Fisher (greater), BH-FDR q-values
""")
if __name__ == "__main__":
    main()

