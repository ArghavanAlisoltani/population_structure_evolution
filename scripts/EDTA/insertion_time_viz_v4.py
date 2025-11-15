#!/usr/bin/env python3
"""
Summarize & visualize LTR_retriever *.list (pass.list) files.

Header example (tab-delimited; note the single header cell "5_TSD 3_TSD"):
#LTR_loc	Category	Motif	TSD	5_TSD 3_TSD	Internal	Identity	Strand	SuperFamily	TE_type	Insertion_Time

Usage:
  python insertion_time_viz_v4.py <input.pass.list> <outdir/>

Outputs:
  <outdir>/summary_by_superfamily.tsv
  <outdir>/summary_overall.tsv
  <outdir>/hist_insertion_time.png
  <outdir>/box_insertion_time_by_superfamily.png
  <outdir>/scatter_identity_vs_insertion_time.png
  <outdir>/clean_passlist.tsv
"""

import os, re, sys, pathlib
import numpy as np
import pandas as pd
import matplotlib

# headless-safe plotting
if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ensure_dir(d):
    pathlib.Path(d).mkdir(parents=True, exist_ok=True)

def strip_prefix(val, prefixes=('motif:', 'TSD:', 'tsd:', 'IN:', 'in:')):
    s = str(val)
    for p in prefixes:
        if s.startswith(p):
            return s[len(p):]
    return s

def find_col(df, targets):
    """Return actual df column matching any target after normalization."""
    normmap = {}
    for c in df.columns:
        key = str(c).strip().lstrip('#').replace(' ', '_').lower()
        normmap[key] = c
    for t in targets:
        key = str(t).strip().lstrip('#').replace(' ', '_').lower()
        if key in normmap:
            return normmap[key]
    return None

def rename_unnamed_after(df, colname, newname):
    """If the next column is Unnamed:*, rename it to newname."""
    try:
        idx = df.columns.get_loc(colname)
    except Exception:
        return
    if idx + 1 < len(df.columns):
        nxt = df.columns[idx + 1]
        if str(nxt).startswith("Unnamed"):
            df.rename(columns={nxt: newname}, inplace=True)

def main(infile, outdir):
    ensure_dir(outdir)

    # ---- read ----
    df = pd.read_csv(infile, sep='\t', dtype=str, engine='python')

    # ---- detect location column ----
    loc_col = find_col(df, ['LTR_loc', '#LTR_loc'])
    if not loc_col:
        sys.stderr.write(f"[ERROR] Could not find LTR location column (#LTR_loc/LTR_loc). Seen: {list(df.columns)}\n")
        sys.exit(2)

    # ---- handle header oddity: '5_TSD 3_TSD' ----
    merged = find_col(df, ['5_TSD 3_TSD'])
    has_5  = find_col(df, ['5_TSD'])
    has_3  = find_col(df, ['3_TSD'])

    if merged and (not has_5 or not has_3):
        two = df[merged].astype(str).str.split(r'\s+', n=1, expand=True)
        df['5_TSD'] = two[0]
        if hasattr(two, "columns") and (1 in list(two.columns)):
            df['3_TSD'] = two[1]
        else:
            rename_unnamed_after(df, merged, '3_TSD')
            if '3_TSD' not in df.columns:
                df['3_TSD'] = pd.NA
        df.drop(columns=[merged], inplace=True)
    else:
        if has_5 and not has_3:
            rename_unnamed_after(df, has_5, '3_TSD')

    # ---- clean simple text fields ----
    for col in ['Motif','TSD','5_TSD','3_TSD','Internal','SuperFamily','TE_type','Category','Strand']:
        if col in df.columns:
            df[col] = df[col].astype(str).map(strip_prefix)

    # ---- parse coordinates via regex extract (robust) ----
    # loc like: contig:a..b (a may be > b)
    extracted = df[loc_col].astype(str).str.extract(r'^(?P<Contig>[^:]+):(?P<a>\d+)\.\.(?P<b>\d+)$', expand=True)
    # numeric
    for c in ['a','b']:
        extracted[c] = pd.to_numeric(extracted[c], errors='coerce')

    # compute Start/End/Mid (order-agnostic)
    start = extracted[['a','b']].min(axis=1)
    end   = extracted[['a','b']].max(axis=1)
    mid   = ((start + end) // 2).astype('Int64')

    df['Contig'] = extracted['Contig']
    df['Start']  = start.astype('Int64')
    df['End']    = end.astype('Int64')
    df['Mid']    = mid

    # ---- parse Internal range if present ----
    if 'Internal' in df.columns:
        ir = df['Internal'].astype(str).str.extract(r'^(?:IN:)?(?P<IN_Start>\d+)\.\.(?P<IN_End>\d+)$', expand=True)
        df['IN_Start']  = pd.to_numeric(ir['IN_Start'], errors='coerce').astype('Int64')
        df['IN_End']    = pd.to_numeric(ir['IN_End'], errors='coerce').astype('Int64')
        df['IN_Length'] = (df['IN_End'] - df['IN_Start'] + 1).astype('Int64')

    # ---- numeric conversions ----
    if 'Identity' in df.columns:
        df['Identity'] = pd.to_numeric(df['Identity'], errors='coerce')
        df['identity_frac'] = df['Identity'].where(df['Identity'] <= 1.0, df['Identity'] / 100.0)
        df['identity_pct']  = df['identity_frac'] * 100.0

    if 'Insertion_Time' in df.columns:
        df['Insertion_Time'] = pd.to_numeric(df['Insertion_Time'], errors='coerce')  # typically years
        df['Insertion_Myr']  = df['Insertion_Time'] / 1e6

    # element span
    df['Length'] = (df['End'] - df['Start'] + 1).astype('Int64')

    # ---- summaries ----
    if 'SuperFamily' in df.columns:
        by_fam = (df.groupby('SuperFamily', dropna=False)
                    .agg(n=('SuperFamily','size'),
                         median_Myr=('Insertion_Myr','median'),
                         q1_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x, 25)),
                         q3_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x, 75)),
                         median_identity=('identity_pct','median'))
                    .reset_index()
                    .sort_values(['n','median_Myr'], ascending=[False, True]))
    else:
        by_fam = pd.DataFrame(columns=['SuperFamily','n','median_Myr','q1_Myr','q3_Myr','median_identity'])

    by_fam.to_csv(os.path.join(outdir, 'summary_by_superfamily.tsv'), sep='\t', index=False)

    overall = {
        'n_total': int(len(df)),
        'n_pass': int((df['Category'] == 'pass').sum()) if 'Category' in df.columns else np.nan,
        'median_Myr': float(df['Insertion_Myr'].median()) if 'Insertion_Myr' in df.columns else np.nan,
        'q1_Myr': float(np.nanpercentile(df['Insertion_Myr'], 25)) if 'Insertion_Myr' in df.columns else np.nan,
        'q3_Myr': float(np.nanpercentile(df['Insertion_Myr'], 75)) if 'Insertion_Myr' in df.columns else np.nan,
        'median_identity_pct': float(df['identity_pct'].median()) if 'identity_pct' in df.columns else np.nan,
        'contigs': int(df['Contig'].nunique())
    }
    pd.DataFrame([overall]).to_csv(os.path.join(outdir, 'summary_overall.tsv'), sep='\t', index=False)

    # ---- plots ----
    if 'Insertion_Myr' in df.columns:
        plt.figure(figsize=(7, 4.5))
        df['Insertion_Myr'].dropna().plot(kind='hist', bins=50)
        plt.xlabel('Insertion time (Myr)'); plt.ylabel('Count')
        plt.title('LTR insertion time distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'hist_insertion_time.png'), dpi=180)
        plt.close()

    if 'SuperFamily' in df.columns and df['SuperFamily'].nunique() > 1 and 'Insertion_Myr' in df.columns:
        order = by_fam['SuperFamily'].tolist()
        data = [df.loc[df['SuperFamily'] == fam, 'Insertion_Myr'].dropna() for fam in order]
        plt.figure(figsize=(8, 4.5))
        plt.boxplot(data, labels=order, showfliers=False)
        plt.ylabel('Insertion time (Myr)')
        plt.title('Insertion time by SuperFamily')
        plt.xticks(rotation=30, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'box_insertion_time_by_superfamily.png'), dpi=180)
        plt.close()

    if 'identity_pct' in df.columns and 'Insertion_Myr' in df.columns:
        plt.figure(figsize=(6.5, 4.5))
        plt.scatter(df['identity_pct'], df['Insertion_Myr'], s=8, alpha=0.6)
        plt.xlabel('LTR–LTR identity (%)')
        plt.ylabel('Insertion time (Myr)')
        plt.title('Identity vs insertion time')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'scatter_identity_vs_insertion_time.png'), dpi=180)
        plt.close()

    # ---- cleaned table ----
    keep_cols = [
        'Contig','Start','End','Mid','Strand','SuperFamily','TE_type',
        'Category','Motif','TSD','5_TSD','3_TSD',
        'IN_Start','IN_End','IN_Length',
        'identity_frac','identity_pct','Insertion_Time','Insertion_Myr','Length',
        loc_col
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df[keep_cols].to_csv(os.path.join(outdir, 'clean_passlist.tsv'), sep='\t', index=False)

    print(f"[DONE] Wrote outputs to: {outdir}")

if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        sys.stderr.write("Usage: python insertion_time_viz_v4.py <input.pass.list> <outdir/>\n")
        sys.exit(1)
    infile = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) == 3 else "passlist_out"
    main(infile, outdir)

