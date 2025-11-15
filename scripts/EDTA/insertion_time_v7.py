#!/usr/bin/env python3
"""
insertion_time_v7.py
Summarize & visualize LTR_retriever *.list AND rescale insertion times, plus:
- Copia vs Gypsy boxplots (original + rescaled)
- Summary stats for Copia/Gypsy overall and by scaffold

Usage:
  python insertion_time_v7.py <input_or_concat.tsv> <outdir/> \
      --rate-new 2.2e-9 --rate-old 1.3e-8 --gen-years 50
"""

import os, re, sys, argparse, pathlib
import numpy as np
import pandas as pd
import matplotlib
if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ensure_dir(d): pathlib.Path(d).mkdir(parents=True, exist_ok=True)

def strip_prefix(val, prefixes=('motif:', 'TSD:', 'tsd:', 'IN:', 'in:')):
    s = str(val)
    for p in prefixes:
        if s.startswith(p): return s[len(p):]
    return s

def find_col(df, targets):
    normmap = {}
    for c in df.columns:
        key = str(c).strip().lstrip('#').replace(' ', '_').lower()
        normmap[key] = c
    for t in targets:
        key = str(t).strip().lstrip('#').replace(' ', '_').lower()
        if key in normmap: return normmap[key]
    return None

def rename_unnamed_after(df, colname, newname):
    try: idx = df.columns.get_loc(colname)
    except Exception: return
    if idx + 1 < len(df.columns):
        nxt = df.columns[idx + 1]
        if str(nxt).startswith("Unnamed"):
            df.rename(columns={nxt: newname}, inplace=True)

def safe_name(s): return re.sub(r'[^A-Za-z0-9._-]+', '_', str(s).strip())

def normalize_family(x):
    s = str(x).strip().lower()
    if "gypsy" in s: return "Gypsy"
    if "copia" in s: return "Copia"
    return None

def add_rescaled_columns(df, rate_old, rate_new, gen_years):
    if 'Insertion_Time' not in df.columns: return df
    factor = rate_old / rate_new
    it = pd.to_numeric(df['Insertion_Time'], errors='coerce')
    df['Insertion_Time_rescaled_years'] = it * factor
    df['Insertion_Time_rescaled_Myr']   = df['Insertion_Time_rescaled_years'] / 1e6
    df['Age_gen_rescaled']              = df['Insertion_Time_rescaled_years'] / float(gen_years)
    return df

def plot_hist(values, outpath, title, xlabel='Insertion time (Myr)', bins=50):
    vals = pd.Series(values).dropna().values
    if vals.size == 0: return
    plt.figure(figsize=(7, 4.5))
    plt.hist(vals, bins=bins)
    plt.xlabel(xlabel); plt.ylabel('Count'); plt.title(title)
    plt.tight_layout(); plt.savefig(outpath, dpi=180); plt.close()

def familia_boxplot(df, col, outpath, title):
    # col is 'Insertion_Myr' or 'Insertion_Time_rescaled_Myr'
    tmp = df[['Family', col]].dropna()
    tmp = tmp[tmp['Family'].isin(['Copia','Gypsy'])]
    if tmp.empty: return
    order = ['Copia','Gypsy']
    data = [tmp.loc[tmp['Family']==fam, col].values for fam in order if fam in tmp['Family'].unique()]
    labels = [fam for fam in order if fam in tmp['Family'].unique()]
    if len(data)==0: return
    plt.figure(figsize=(6.5, 4.5))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.ylabel('Insertion time (Myr)'); plt.title(title)
    plt.tight_layout(); plt.savefig(outpath, dpi=180); plt.close()

def summarize_family(df, col):
    # overall Copia/Gypsy
    sub = df[df['Family'].isin(['Copia','Gypsy'])]
    if sub.empty:
        return pd.DataFrame(columns=['Family','n',f'median_{col}',f'q1_{col}',f'q3_{col}',f'mean_{col}'])
    return (sub.groupby('Family', dropna=False)[col]
            .agg(n='count',
                 median=lambda x: np.nanmedian(x),
                 q1=lambda x: np.nanpercentile(x,25),
                 q3=lambda x: np.nanpercentile(x,75),
                 mean=lambda x: np.nanmean(x))
            .reset_index()
            .rename(columns={'median':f'median_{col}','q1':f'q1_{col}','q3':f'q3_{col}','mean':f'mean_{col}'}))

def summarize_family_by_scaffold(df, col):
    sub = df[df['Family'].isin(['Copia','Gypsy'])]
    if sub.empty:
        return pd.DataFrame(columns=['Contig','Family','n',f'median_{col}',f'q1_{col}',f'q3_{col}',f'mean_{col}'])
    g = sub.groupby(['Contig','Family'], dropna=False)[col]
    out = g.agg(n='count',
                median=lambda x: np.nanmedian(x),
                q1=lambda x: np.nanpercentile(x,25),
                q3=lambda x: np.nanpercentile(x,75),
                mean=lambda x: np.nanmean(x)).reset_index()
    return out.rename(columns={'median':f'median_{col}','q1':f'q1_{col}','q3':f'q3_{col}','mean':f'mean_{col}'})

def main(infile, outdir, rate_old, rate_new, gen_years):
    ensure_dir(outdir)
    df = pd.read_csv(infile, sep='\t', dtype=str, engine='python')

    loc_col = find_col(df, ['LTR_loc', '#LTR_loc'])
    if not loc_col:
        sys.stderr.write(f"[ERROR] Could not find LTR location column (#LTR_loc/LTR_loc). Seen: {list(df.columns)}\n")
        sys.exit(2)

    merged = find_col(df, ['5_TSD 3_TSD'])
    has_5  = find_col(df, ['5_TSD'])
    has_3  = find_col(df, ['3_TSD'])
    if merged and (not has_5 or not has_3):
        two = df[merged].astype(str).str.split(r'\s+', n=1, expand=True)
        df['5_TSD'] = two[0]
        if hasattr(two, "columns") and (1 in list(two.columns)): df['3_TSD'] = two[1]
        else:
            rename_unnamed_after(df, merged, '3_TSD')
            if '3_TSD' not in df.columns: df['3_TSD'] = pd.NA
        df.drop(columns=[merged], inplace=True)
    else:
        if has_5 and not has_3: rename_unnamed_after(df, has_5, '3_TSD')

    for col in ['Motif','TSD','5_TSD','3_TSD','Internal','SuperFamily','TE_type','Category','Strand']:
        if col in df.columns: df[col] = df[col].astype(str).map(strip_prefix)

    extracted = df[loc_col].astype(str).str.extract(r'^(?P<Contig>[^:]+):(?P<a>\d+)\.\.(?P<b>\d+)$', expand=True)
    for c in ['a','b']: extracted[c] = pd.to_numeric(extracted[c], errors='coerce')
    start = extracted[['a','b']].min(axis=1); end = extracted[['a','b']].max(axis=1)
    mid = ((start + end) // 2).astype('Int64')
    df['Contig'] = extracted['Contig']; df['Start'] = start.astype('Int64'); df['End'] = end.astype('Int64'); df['Mid'] = mid

    if 'Internal' in df.columns:
        ir = df['Internal'].astype(str).str.extract(r'^(?:IN:)?(?P<IN_Start>\d+)\.\.(?P<IN_End>\d+)$', expand=True)
        df['IN_Start'] = pd.to_numeric(ir['IN_Start'], errors='coerce').astype('Int64')
        df['IN_End']   = pd.to_numeric(ir['IN_End'], errors='coerce').astype('Int64')
        df['IN_Length']= (df['IN_End'] - df['IN_Start'] + 1).astype('Int64')

    if 'Identity' in df.columns:
        df['Identity'] = pd.to_numeric(df['Identity'], errors='coerce')
        df['identity_frac'] = df['Identity'].where(df['Identity'] <= 1.0, df['Identity']/100.0)
        df['identity_pct']  = df['identity_frac'] * 100.0

    if 'Insertion_Time' in df.columns:
        df['Insertion_Time'] = pd.to_numeric(df['Insertion_Time'], errors='coerce')
        df['Insertion_Myr']  = df['Insertion_Time'] / 1e6

    df['Length'] = (df['End'] - df['Start'] + 1).astype('Int64')

    # Normalized family tag for Copia/Gypsy analyses
    if 'SuperFamily' in df.columns:
        df['Family'] = df['SuperFamily'].map(normalize_family)
    else:
        df['Family'] = None

    # RESCALE
    df = add_rescaled_columns(df, rate_old=rate_old, rate_new=rate_new, gen_years=gen_years)

    # ---- Existing summaries ----
    if 'SuperFamily' in df.columns:
        by_fam = (df.groupby('SuperFamily', dropna=False)
                    .agg(n=('SuperFamily','size'),
                         median_Myr=('Insertion_Myr','median'),
                         median_Myr_rescaled=('Insertion_Time_rescaled_Myr','median'),
                         q1_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x,25)),
                         q3_Myr =('Insertion_Myr', lambda x: np.nanpercentile(x,75)),
                         median_identity=('identity_pct','median'))
                    .reset_index()
                    .sort_values(['n','median_Myr'], ascending=[False, True]))
    else:
        by_fam = pd.DataFrame(columns=['SuperFamily','n','median_Myr','median_Myr_rescaled','q1_Myr','q3_Myr','median_identity'])
    by_fam.to_csv(os.path.join(outdir, 'summary_by_superfamily.tsv'), sep='\t', index=False)

    overall = {
        'n_total': int(len(df)),
        'n_pass': int((df['Category']=='pass').sum()) if 'Category' in df.columns else np.nan,
        'median_Myr': float(df['Insertion_Myr'].median()) if 'Insertion_Myr' in df.columns else np.nan,
        'median_Myr_rescaled': float(df['Insertion_Time_rescaled_Myr'].median()) if 'Insertion_Time_rescaled_Myr' in df.columns else np.nan,
        'q1_Myr': float(np.nanpercentile(df['Insertion_Myr'],25)) if 'Insertion_Myr' in df.columns else np.nan,
        'q3_Myr': float(np.nanpercentile(df['Insertion_Myr'],75)) if 'Insertion_Myr' in df.columns else np.nan,
        'median_identity_pct': float(df['identity_pct'].median()) if 'identity_pct' in df.columns else np.nan,
        'contigs': int(df['Contig'].nunique()),
        'rate_old': rate_old, 'rate_new': rate_new, 'gen_years': gen_years
    }
    pd.DataFrame([overall]).to_csv(os.path.join(outdir, 'summary_overall.tsv'), sep='\t', index=False)

    # ---- Overall plots ----
    if 'Insertion_Myr' in df.columns:
        plot_hist(df['Insertion_Myr'], os.path.join(outdir,'hist_insertion_time.png'),
                  'LTR insertion time (original rate)', 'Insertion time (Myr)', bins=50)
    if 'Insertion_Time_rescaled_Myr' in df.columns:
        plot_hist(df['Insertion_Time_rescaled_Myr'], os.path.join(outdir,'hist_insertion_time_rescaled.png'),
                  'LTR insertion time (rescaled)', 'Insertion time (Myr, rescaled)', bins=50)

    if 'SuperFamily' in df.columns and df['SuperFamily'].nunique()>1 and 'Insertion_Myr' in df.columns:
        order = by_fam['SuperFamily'].tolist()
        data = [df.loc[df['SuperFamily']==fam, 'Insertion_Myr'].dropna() for fam in order]
        if any(len(x)>0 for x in data):
            plt.figure(figsize=(8, 4.5))
            plt.boxplot(data, labels=order, showfliers=False)
            plt.ylabel('Insertion time (Myr)'); plt.title('Insertion time by SuperFamily (original rate)')
            plt.xticks(rotation=30, ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outdir, 'box_insertion_time_by_superfamily.png'), dpi=180); plt.close()

    if 'identity_pct' in df.columns and 'Insertion_Myr' in df.columns:
        plt.figure(figsize=(6.5, 4.5))
        plt.scatter(df['identity_pct'], df['Insertion_Myr'], s=8, alpha=0.6)
        plt.xlabel('LTR–LTR identity (%)'); plt.ylabel('Insertion time (Myr, original rate)')
        plt.title('Identity vs insertion time'); plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'scatter_identity_vs_insertion_time.png'), dpi=180); plt.close()

    # ---- NEW: Copia vs Gypsy boxplots (overall) ----
    if 'Family' in df.columns:
        if 'Insertion_Myr' in df.columns:
            familia_boxplot(df, 'Insertion_Myr',
                            os.path.join(outdir,'box_copia_vs_gypsy.png'),
                            'Copia vs Gypsy (original rate)')
        if 'Insertion_Time_rescaled_Myr' in df.columns:
            familia_boxplot(df, 'Insertion_Time_rescaled_Myr',
                            os.path.join(outdir,'box_copia_vs_gypsy_rescaled.png'),
                            'Copia vs Gypsy (rescaled)')

    # ---- NEW: summary stats Copia/Gypsy overall and by scaffold (original + rescaled) ----
    overall_o = summarize_family(df.assign(**{'Insertion_Myr':pd.to_numeric(df.get('Insertion_Myr'), errors='coerce')}), 'Insertion_Myr') \
                if 'Insertion_Myr' in df.columns else pd.DataFrame()
    overall_r = summarize_family(df.assign(**{'Insertion_Time_rescaled_Myr':pd.to_numeric(df.get('Insertion_Time_rescaled_Myr'), errors='coerce')}), 'Insertion_Time_rescaled_Myr') \
                if 'Insertion_Time_rescaled_Myr' in df.columns else pd.DataFrame()
    # merge (outer) on Family
    if not overall_o.empty or not overall_r.empty:
        all_overall = pd.merge(overall_o, overall_r, on='Family', how='outer')
        all_overall.to_csv(os.path.join(outdir, 'summary_copia_gypsy_overall.tsv'), sep='\t', index=False)

    bysc_o = summarize_family_by_scaffold(df.assign(**{'Insertion_Myr':pd.to_numeric(df.get('Insertion_Myr'), errors='coerce')}), 'Insertion_Myr') \
             if 'Insertion_Myr' in df.columns else pd.DataFrame()
    bysc_r = summarize_family_by_scaffold(df.assign(**{'Insertion_Time_rescaled_Myr':pd.to_numeric(df.get('Insertion_Time_rescaled_Myr'), errors='coerce')}), 'Insertion_Time_rescaled_Myr') \
             if 'Insertion_Time_rescaled_Myr' in df.columns else pd.DataFrame()
    if not bysc_o.empty or not bysc_r.empty:
        all_bysc = pd.merge(bysc_o, bysc_r, on=['Contig','Family'], how='outer')
        all_bysc.to_csv(os.path.join(outdir, 'summary_copia_gypsy_by_scaffold.tsv'), sep='\t', index=False)

    # ---- Per-class histograms (original + rescaled) ----
    if 'SuperFamily' in df.columns:
        class_dir = os.path.join(outdir, 'class_hists'); ensure_dir(class_dir)
        all_vals_o = df['Insertion_Myr'].dropna().values if 'Insertion_Myr' in df.columns else np.array([])
        bins_o = np.linspace(np.nanmin(all_vals_o), np.nanmax(all_vals_o), 51) if all_vals_o.size>0 else np.linspace(0.0,1.0,51)
        all_vals_r = df['Insertion_Time_rescaled_Myr'].dropna().values if 'Insertion_Time_rescaled_Myr' in df.columns else np.array([])
        bins_r = np.linspace(np.nanmin(all_vals_r), np.nanmax(all_vals_r), 51) if all_vals_r.size>0 else np.linspace(0.0,1.0,51)

        for fam in sorted(df['SuperFamily'].dropna().unique()):
            tag = safe_name(fam)
            if 'Insertion_Myr' in df.columns:
                arr = df.loc[df['SuperFamily']==fam, 'Insertion_Myr'].dropna().values
                if arr.size:
                    plt.figure(figsize=(7,4.5)); plt.hist(arr, bins=bins_o)
                    plt.xlabel('Insertion time (Myr)'); plt.ylabel('Count')
                    plt.title(f'Insertion time — {fam} (original, n={arr.size})')
                    plt.tight_layout(); plt.savefig(os.path.join(class_dir, f"hist_insertion_time_{tag}.png"), dpi=180); plt.close()
                    counts, edges = np.histogram(arr, bins=bins_o)
                    mids = (edges[:-1]+edges[1:])/2.0
                    pd.DataFrame({'bin_left':edges[:-1],'bin_right':edges[1:],'bin_mid_Myr':mids,'count':counts}) \
                        .to_csv(os.path.join(class_dir, f"hist_counts_{tag}.tsv"), sep='\t', index=False)
            if 'Insertion_Time_rescaled_Myr' in df.columns:
                arrR = df.loc[df['SuperFamily']==fam, 'Insertion_Time_rescaled_Myr'].dropna().values
                if arrR.size:
                    plt.figure(figsize=(7,4.5)); plt.hist(arrR, bins=bins_r)
                    plt.xlabel('Insertion time (Myr, rescaled)'); plt.ylabel('Count')
                    plt.title(f'Insertion time — {fam} (rescaled, n={arrR.size})')
                    plt.tight_layout(); plt.savefig(os.path.join(class_dir, f"hist_insertion_time_rescaled_{tag}.png"), dpi=180); plt.close()

    # ---- Cleaned table ----
    keep_cols = ['Contig','Start','End','Mid','Strand','SuperFamily','Family','TE_type',
                 'Category','Motif','TSD','5_TSD','3_TSD','IN_Start','IN_End','IN_Length',
                 'identity_frac','identity_pct','Insertion_Time','Insertion_Myr',
                 'Insertion_Time_rescaled_years','Insertion_Time_rescaled_Myr','Age_gen_rescaled',
                 'Length', loc_col]
    # keep SourceFile if present (from concatenation)
    if 'SourceFile' in df.columns: keep_cols.append('SourceFile')
    keep_cols = [c for c in keep_cols if c in df.columns]
    # numeric cast for output clarity
    for c in ['Insertion_Myr','Insertion_Time_rescaled_Myr','identity_pct','Length','Age_gen_rescaled']:
        if c in df.columns: df[c] = pd.to_numeric(df[c], errors='coerce')
    df[keep_cols].to_csv(os.path.join(outdir, 'clean_passlist.tsv'), sep='\t', index=False)

    print(f"[DONE] Outputs -> {outdir}")
    print(f"   Used rate_old={rate_old:.3e}, rate_new={rate_new:.3e}, gen_years={gen_years}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Summarize/visualize LTR pass.list, rescale ages, and add Copia vs Gypsy stats/plots.")
    p.add_argument("infile", help="Input .pass.list or concatenated TSV")
    p.add_argument("outdir", nargs="?", default="passlist_out", help="Output directory")
    p.add_argument("--rate-old", type=float, default=1.3e-8, help="Old per-year rate used originally")
    p.add_argument("--rate-new", type=float, default=2.2e-9, help="New per-year rate to rescale to")
    p.add_argument("--gen-years", type=float, default=50.0, help="Years per generation for Age_gen_rescaled")
    args = p.parse_args()
    main(args.infile, args.outdir, rate_old=args.rate_old, rate_new=args.rate_new, gen_years=args.gen_years)

