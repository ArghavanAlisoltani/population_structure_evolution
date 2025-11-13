#!/usr/bin/env python3
"""
Summarize & visualize LTR_retriever *.list (pass.list) files.

Expected columns (TSV):
#LTR_loc  Category  Motif  TSD  5_TSD  3_TSD  Internal  Identity  Strand  SuperFamily  TE_type  Insertion_Time

Example run:
  python summarize_passlist.py input.pass.list outdir/

Outputs:
  outdir/summary_by_superfamily.tsv
  outdir/summary_overall.tsv
  outdir/hist_insertion_time.png
  outdir/box_insertion_time_by_superfamily.png
  outdir/scatter_identity_vs_insertion_time.png
"""

import os, re, sys, math, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- helpers ----------
loc_re = re.compile(r'^(?P<contig>[^:]+):(?P<a>\d+)\.\.(?P<b>\d+)$')
rng_re = re.compile(r'^(?:IN:)?(?P<a>\d+)\.\.(?P<b>\d+)$', re.IGNORECASE)

def parse_loc(s):
    m = loc_re.match(str(s))
    if not m:
        return pd.Series([None, np.nan, np.nan, np.nan])
    a = int(m.group('a')); b = int(m.group('b'))
    start = min(a,b); end = max(a,b)
    mid = (start + end) // 2
    return pd.Series([m.group('contig'), start, end, mid])

def parse_range(s):
    m = rng_re.match(str(s))
    if not m:
        return pd.Series([np.nan, np.nan, np.nan])
    a = int(m.group('a')); b = int(m.group('b'))
    start = min(a,b); end = max(a,b)
    length = end - start + 1
    return pd.Series([start, end, length])

def strip_prefix(val, prefixes=('motif:', 'TSD:', 'tsd:', 'IN:', 'in:')):
    s = str(val)
    for p in prefixes:
        if s.startswith(p):
            return s[len(p):]
    return s

def ensure_dir(d):
    pathlib.Path(d).mkdir(parents=True, exist_ok=True)

# ---------- main ----------
def main(infile, outdir):
    ensure_dir(outdir)

    # read tsv
    df = pd.read_csv(infile, sep='\t', dtype=str)
    # fix odd header tokens
    df = df.rename(columns={
        '#LTR_loc': 'LTR_loc',
        '5_TSD 3_TSD': '5_TSD_3_TSD'  # safety if header came in merged; we’ll split below if needed
    })

    # if a merged '5_TSD 3_TSD' slipped in as one column, try to split by whitespace
    if '5_TSD_3_TSD' in df.columns and '5_TSD' not in df.columns and '3_TSD' not in df.columns:
        # attempt split into two new columns from the merged field
        two = df['5_TSD_3_TSD'].str.split(r'\s+', n=1, expand=True)
        df['5_TSD'] = two[0]
        df['3_TSD'] = two[1]
        df.drop(columns=['5_TSD_3_TSD'], inplace=True)

    # basic cleaning
    for col in ['Motif','TSD','5_TSD','3_TSD','Internal','SuperFamily','TE_type','Category','Strand']:
        if col in df.columns:
            df[col] = df[col].astype(str).map(strip_prefix)

    # parse location columns
    loc_cols = df['LTR_loc'].map(parse_loc)
    loc_cols.columns = ['Contig','Start','End','Mid']
    df = pd.concat([df, loc_cols], axis=1)

    if 'Internal' in df.columns:
        intr = df['Internal'].map(parse_range)
        intr.columns = ['IN_Start','IN_End','IN_Length']
        df = pd.concat([df, intr], axis=1)

    # numeric columns
    for col in ['Identity', 'Insertion_Time']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Identity: file uses fraction (e.g., 0.9339). Convert to percent as well.
    if 'Identity' in df.columns:
        df['identity_frac'] = df['Identity'].where(df['Identity'] <= 1.0, df['Identity']/100.0)
        df['identity_pct']  = df['identity_frac'] * 100.0

    # Insertion time: LTR_retriever default units ~ years (with μ≈1.3e-8 y^-1)
    if 'Insertion_Time' in df.columns:
        df['Insertion_Myr'] = df['Insertion_Time'] / 1e6

    # length of full LTR-RT span
    df['Length'] = (df['End'] - df['Start'] + 1).astype('Int64')

    # ---- summaries ----
    # by SuperFamily
    by_fam = (df
              .groupby('SuperFamily', dropna=False)
              .agg(n=('SuperFamily','size'),
                   median_Myr=('Insertion_Myr','median'),
                   q1_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x, 25)),
                   q3_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x, 75)),
                   median_identity=('identity_pct','median'))
              .reset_index()
              .sort_values(['n','median_Myr'], ascending=[False, True]))

    by_fam.to_csv(os.path.join(outdir, 'summary_by_superfamily.tsv'), sep='\t', index=False)

    # overall summary
    overall = {
        'n_total': len(df),
        'n_pass': int((df.get('Category','') == 'pass').sum()) if 'Category' in df.columns else np.nan,
        'median_Myr': float(df['Insertion_Myr'].median()) if 'Insertion_Myr' in df.columns else np.nan,
        'q1_Myr': float(np.nanpercentile(df['Insertion_Myr'], 25)) if 'Insertion_Myr' in df.columns else np.nan,
        'q3_Myr': float(np.nanpercentile(df['Insertion_Myr'], 75)) if 'Insertion_Myr' in df.columns else np.nan,
        'median_identity_pct': float(df['identity_pct'].median()) if 'identity_pct' in df.columns else np.nan,
        'contigs': int(df['Contig'].nunique())
    }
    pd.DataFrame([overall]).to_csv(os.path.join(outdir, 'summary_overall.tsv'), sep='\t', index=False)

    # ---- plots ----
    plt.figure(figsize=(7,4.5))
    df['Insertion_Myr'].dropna().plot(kind='hist', bins=50)
    plt.xlabel('Insertion time (Myr)')
    plt.ylabel('Count')
    plt.title('LTR insertion time distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'hist_insertion_time.png'), dpi=180)
    plt.close()

    if 'SuperFamily' in df.columns and df['SuperFamily'].nunique() > 1:
        plt.figure(figsize=(8,4.5))
        order = by_fam['SuperFamily'].tolist()
        data = [df.loc[df['SuperFamily']==fam, 'Insertion_Myr'].dropna() for fam in order]
        plt.boxplot(data, labels=order, showfliers=False)
        plt.ylabel('Insertion time (Myr)')
        plt.title('Insertion time by SuperFamily')
        plt.xticks(rotation=30, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir,'box_insertion_time_by_superfamily.png'), dpi=180)
        plt.close()

    if 'identity_pct' in df.columns and 'Insertion_Myr' in df.columns:
        plt.figure(figsize=(6.5,4.5))
        plt.scatter(df['identity_pct'], df['Insertion_Myr'], s=8, alpha=0.6)
        plt.xlabel('LTR–LTR identity (%)')
        plt.ylabel('Insertion time (Myr)')
        plt.title('Identity vs insertion time')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir,'scatter_identity_vs_insertion_time.png'), dpi=180)
        plt.close()

    # write a cleaned table (useful downstream)
    keep_cols = ['Contig','Start','End','Mid','Strand','SuperFamily','TE_type',
                 'Category','Motif','TSD','5_TSD','3_TSD','IN_Start','IN_End','IN_Length',
                 'identity_frac','identity_pct','Insertion_Time','Insertion_Myr','Length','LTR_loc']
    keep_cols = [c for c in keep_cols if c in df.columns]
    df[keep_cols].to_csv(os.path.join(outdir, 'clean_passlist.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python summarize_passlist.py <input.pass.list> <outdir/>\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
