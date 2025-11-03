# save as: tag_gwas_in_sigFST_windows.py
# run:
# python tag_gwas_in_sigFST_windows.py Aria_MTAG_summary.txt fst_windows_with_TEcounts__augmented.tsv 0.05 GWAS_with_FSTwindow_flags.tsv
import sys, pandas as pd, numpy as np

gwas_path, wins_path = sys.argv[1], sys.argv[2]
qthr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
out_path = sys.argv[4] if len(sys.argv) > 4 else "GWAS_with_FSTwindow_flags.tsv"

gw = pd.read_csv(gwas_path, sep="\t")
w  = pd.read_csv(wins_path, sep="\t")

gw["scaffold_norm"] = gw["scaffold"].str.strip().str.lower()
w["CHROM_norm"]     = w["CHROM"].astype(str).str.strip().str.lower()

qcols = [c for c in ["q_poi_3vs4","q_poi_3vs5","q_poi_4vs5"] if c in w.columns]
for qc in qcols:
    lab = "sigFST_" + qc.split("_",2)[-1]  # e.g., 3vs4
    sig = w.loc[w[qc] <= qthr, ["CHROM","CHROM_norm","WIN_START","WIN_END"]].copy()
    if sig.empty:
        gw[lab] = False
        continue
    sig["WIN_START"] = sig["WIN_START"].astype(int)
    sig["WIN_END"]   = sig["WIN_END"].astype(int)

    # Mark SNPs inside any significant window (per scaffold)
    flag = np.zeros(len(gw), dtype=bool)
    for chrom, sub in sig.groupby("CHROM_norm"):
        idx = np.where(gw["scaffold_norm"].values == chrom)[0]
        if idx.size == 0: continue
        pos = gw.loc[idx, "position"].astype(int).values
        s = sub[["WIN_START","WIN_END"]].values
        # For each SNP position, check if it lies in any [start, end] (vectorized enough at this scale)
        for j, p in enumerate(pos):
            if np.any((s[:,0] <= p) & (p <= s[:,1])):
                flag[idx[j]] = True
    gw[lab] = flag

# Any comparison?
sig_cols = [c for c in gw.columns if c.startswith("sigFST_")]
gw["sigFST_any"] = gw[sig_cols].any(axis=1) if sig_cols else False

# Write
order = ["scaffold","position","trait","SNP"] + sig_cols + ["sigFST_any"]
other = [c for c in gw.columns if c not in order and c != "scaffold_norm"]
gw[order + other].to_csv(out_path, sep="\t", index=False)
print(f"Wrote: {out_path}")

