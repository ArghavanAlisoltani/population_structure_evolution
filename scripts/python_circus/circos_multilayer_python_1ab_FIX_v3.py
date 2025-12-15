# circos_multilayer_python_1ab_FIX_v3.py
# - FIX: use sector_map[name] instead of circos.sectors[name]
# - Robust LG numeric ordering
# - Carries meta (e.g., FDR) through 1→1a/1b split safely

import re
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pycirclize import Circos

# ===============================
# CONFIG
# ===============================
DATA = {
    "tmrca":      "all_tmrca_corrected_position.tsv",
    "snp":        "positions_All_scaffold_poly_s100.txt",
    "te":         "col8_readable.TEanno.cds.scaf1ab.tsv",
    "mrna":       "1ab_mRNA_ID_merged_interproscan.txt",
    "gwas_orig":  "combined_hits_fdr_0.1_MAF_0.01.tsv",
    "trait_gwas": "FC_sumstat.txt",
    "lg_map":     "canLG.v1_Aria_version.txt",
}

CHR1, CHR1A, CHR1B = "scaffold_1", "scaffold_1a", "scaffold_1b"
GAP_END  = 1418382258.0
OFFSET   = 1427634029.0

SHOW_SNP_DENSITY = True
SNP_WINDOW_BP    = 1_000_000

SHOW_TE          = False
TE_MAJOR_TO_SHOW = {"LTR", "DNA", "Helitron", "LINE", "SINE"}
TE_COLORS = {
    "LTR:Gypsy":  "#1b9e77",
    "LTR:Copia":  "#d95f02",
    "DNA":        "#7570b3",
    "Helitron":   "#e7298a",
    "LINE":       "#66a61e",
    "SINE":       "#e6ab02",
}
SHOW_MRNA        = False

SHOW_TMRCA       = False
TMRCA_WIN_BP     = 500_000
TMRCA_STEP_BP    = None
TMRCA_OFFSET_BP  = 1.0

SHOW_TRAIT_GWAS  = True
FDR_CUTOFF       = 0.10
MAF_CUTOFF       = 0.01
TRAIT_POINT_SIZE = 4.0

SHOW_ORIG_GWAS   = False
ORIG_GWAS_SIZE   = 2.5

OUT_JPEG = "circos_multilayer_python_1ab.jpg"
OUT_PDF  = "circos_multilayer_python_1ab.pdf"
DPI      = 300
FIG_W, FIG_H = 10, 10

# ===============================
# Helpers
# ===============================
def load_tsv(path: str) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists(): return None
    return pd.read_csv(p, sep="\t", dtype=str, low_memory=False)

def normalize_scaffold_val(x: str) -> str:
    if pd.isna(x): return np.nan
    s = str(x).strip().lower()
    m = re.match(r"^s{0,2}scaffold_([0-9]+[ab]?)(_.*)?$", s)
    if m: return f"scaffold_{m.group(1)}"
    if re.match(r"^scaffold_[0-9]+[ab]?$", s): return s
    m = re.match(r"^chr?([0-9]+[ab]?)$", s)
    if m: return f"scaffold_{m.group(1)}"
    if re.match(r"^[0-9]+[ab]?$", s): return f"scaffold_{s}"
    return s

def parse_te_major_sub(s: str) -> tuple[str|None, str|None]:
    if pd.isna(s): return (None, None)
    t = str(s).lower()
    if "helitron" in t: return ("Helitron", "Helitron")
    if re.search(r"(^|[^a-z])line", t): return ("LINE", "LINE")
    if re.search(r"(^|[^a-z])sine", t): return ("SINE", "SINE")
    if re.search(r"(^|[^a-z])(dna|tir|hat|tc.?mar|mutator|cacta|enspm|harbinger|mariner|pogo)", t): return ("DNA", "DNA")
    if "ltr" in t:
        if "gypsy" in t or "ty3" in t: return ("LTR","Gypsy")
        if "copia" in t or "ty1" in t: return ("LTR","Copia")
        return ("LTR", None)
    return (None, None)

# ===============================
# Loaders (standardize columns)
# ===============================
def load_tmrca_std(path: str) -> pd.DataFrame | None:
    df = load_tsv(path);  
    if df is None: return None
    req = ["CHROM","start_tmrca","end_tmrca","mean_tmrca","Low_CI","UP_CI"]
    if not set(req).issubset(df.columns):
        raise ValueError("TMRCA file missing: " + ", ".join(set(req)-set(df.columns)))
    out = pd.DataFrame({
        "CHROM": df["CHROM"].map(normalize_scaffold_val),
        "START": pd.to_numeric(df["start_tmrca"], errors="coerce") + TMRCA_OFFSET_BP,
        "END":   pd.to_numeric(df["end_tmrca"],   errors="coerce") + TMRCA_OFFSET_BP,
        "mean":  pd.to_numeric(df["mean_tmrca"],  errors="coerce"),
        "lo":    pd.to_numeric(df["Low_CI"],      errors="coerce"),
        "hi":    pd.to_numeric(df["UP_CI"],       errors="coerce"),
    }).dropna(subset=["CHROM","START","END"])
    return out[out["END"] >= out["START"]]

def load_snp_std(path: str) -> pd.DataFrame | None:
    df = load_tsv(path); 
    if df is None: return None
    cols = {c.lower(): c for c in df.columns}
    chrom_col = next((cols.get(k) for k in ["chrom","chr","scaffold"]), None)
    pos_col   = next((cols.get(k) for k in ["pos","bp","position"]), None)
    if not chrom_col or not pos_col:
        raise ValueError("SNP file must have CHROM/CHR/scaffold and POS/BP/position.")
    return pd.DataFrame({
        "CHROM": df[chrom_col].map(normalize_scaffold_val),
        "POS":   pd.to_numeric(df[pos_col], errors="coerce"),
    }).dropna(subset=["CHROM","POS"])

def load_te_std(path: str) -> pd.DataFrame | None:
    df = load_tsv(path); 
    if df is None: return None
    cols = {c.lower(): c for c in df.columns}
    chr_col = next((cols.get(k) for k in ["seqid","scaffold","chrom","chr","gff_seqid"]), None)
    st_col  = next((cols.get(k) for k in ["start","gff_start"]), None)
    en_col  = next((cols.get(k) for k in ["end","gff_end"]), None)
    cls_col = next((cols.get(k) for k in ["sequence_ontology","class","te_class","superfamily"]), None)
    if not chr_col or not st_col or not en_col:
        raise ValueError("TE file must include scaffold/seqid and start/end.")
    out = pd.DataFrame({
        "CHROM": df[chr_col].map(normalize_scaffold_val),
        "START": pd.to_numeric(df[st_col], errors="coerce"),
        "END":   pd.to_numeric(df[en_col], errors="coerce"),
        "CLS":   df[cls_col] if cls_col else None
    }).dropna(subset=["CHROM","START","END"])
    out = out[out["END"] >= out["START"]]
    major, sub = [], []
    for val in (out["CLS"] if "CLS" in out.columns else [None]*len(out)):
        M,S = parse_te_major_sub(val)
        major.append(M); sub.append(S)
    out["MAJOR"] = major
    out["SUB"]   = sub
    return out

def load_mrna_std(path: str) -> pd.DataFrame | None:
    df = load_tsv(path); 
    if df is None: return None
    cols = {c.lower(): c for c in df.columns}
    chr_col = next((cols.get(k) for k in ["new_seqid","gff_seqid","seqid","chrom","scaffold"]), None)
    st_col  = next((cols.get(k) for k in ["new_start","gff_start","start"]), None)
    en_col  = next((cols.get(k) for k in ["new_end","gff_end","end"]), None)
    if not chr_col or not st_col or not en_col:
        raise ValueError("mRNA file must have new_seqid/new_start/new_end (or gff_*).")
    out = pd.DataFrame({
        "CHROM": df[chr_col].map(normalize_scaffold_val),
        "START": pd.to_numeric(df[st_col], errors="coerce"),
        "END":   pd.to_numeric(df[en_col], errors="coerce"),
    }).dropna(subset=["CHROM","START","END"])
    return out[out["END"] >= out["START"]]

def load_gwas_orig_std(path: str) -> pd.DataFrame | None:
    df = load_tsv(path); 
    if df is None: return None
    cols = {c.lower(): c for c in df.columns}
    sc_col  = next((cols.get(k) for k in ["scaffold","chr","chrom"]), None)
    pos_col = next((cols.get(k) for k in ["position","pos","bp"]), None)
    tr_col  = next((cols.get(k) for k in ["trait"]), None)
    if not sc_col or not pos_col:
        raise ValueError("Original GWAS must include scaffold/CHR and position/POS/BP.")
    out = pd.DataFrame({
        "CHROM": df[sc_col].map(normalize_scaffold_val),
        "POS":   pd.to_numeric(df[pos_col], errors="coerce"),
        "TRAIT": df[tr_col] if tr_col else ""
    }).dropna(subset=["CHROM","POS"])
    return out

def load_trait_gwas_std(path: str, trait_name: str | None = None) -> pd.DataFrame | None:
    df = load_tsv(path); 
    if df is None: return None
    low = {c.lower(): c for c in df.columns}
    if "scaffold" in low:
        sc = df[low["scaffold"]]
    elif "chr" in low:
        sc = df[low["chr"]]
    elif "chrom" in low:
        sc = df[low["chrom"]]
    elif "snp" in low:
        snp = df[low["snp"]].astype(str).str.lower()
        sc = snp.str.replace(r"^s{0,2}scaffold_([0-9]+[ab]?)(_.*)?$", r"scaffold_\1", regex=True)
    else:
        raise ValueError("Trait GWAS must contain scaffold/chr/chrom or snp.")
    pos_col = low.get("bp") or low.get("pos") or low.get("position")
    if not pos_col:
        raise ValueError("Trait GWAS must contain BP/POS/position.")
    fdr_col = low.get("fdr") or low.get("p")
    if not fdr_col:
        raise ValueError("Trait GWAS must contain 'fdr' or 'p'.")
    maf_col = low.get("maf")

    out = pd.DataFrame({
        "CHROM": pd.Series(sc).map(normalize_scaffold_val),
        "POS":   pd.to_numeric(df[pos_col], errors="coerce"),
        "FDR":   pd.to_numeric(df[fdr_col], errors="coerce"),
        "MAF":   pd.to_numeric(df[maf_col], errors="coerce") if maf_col else np.nan,
        "TRAIT": trait_name or Path(path).name
    }).dropna(subset=["CHROM","POS","FDR"])
    out["FDR"] = out["FDR"].clip(lower=np.finfo(float).tiny, upper=1.0)
    return out

def load_lg_map(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ["Candidate_LG_ID","Scaffold_ID","Strand","Scaffold_length"]:
        if col not in df.columns:
            raise ValueError("LG map missing required columns.")
    out = pd.DataFrame({
        "LG":      df["Candidate_LG_ID"].astype(str),
        "CHROM":   df["Scaffold_ID"].map(normalize_scaffold_val),
        "Strand":  df["Strand"].astype(str),
        "length":  pd.to_numeric(df["Scaffold_length"], errors="coerce"),
    }).dropna(subset=["CHROM","length"])
    out["order_in_LG"] = out.groupby("LG").cumcount() + 1
    out = out.sort_values(["LG","order_in_LG"], key=lambda s: pd.to_numeric(s, errors="coerce"))
    out["lg_scaf_start"] = out.groupby("LG")["length"].cumsum().shift(fill_value=0)
    out["lg_scaf_end"]   = out["lg_scaf_start"] + out["length"]
    return out

# ===============================
# 1a/1b merge/split (with meta)
# ===============================
def merge_1ab_points_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    orig = out["CHROM"].astype(str)
    out.loc[orig == CHR1B, "POS"] = out.loc[orig == CHR1B, "POS"].astype(float) + (OFFSET - 1.0)
    out.loc[orig.isin([CHR1A, CHR1B, CHR1]), "CHROM"] = CHR1
    return out

def split_1_to_1ab_points_with_meta(df: pd.DataFrame, meta_cols: list[str] | None = None) -> pd.DataFrame:
    if meta_cols is None: meta_cols = []
    df = df.copy()
    rows = []
    m_a = (df["CHROM"] == CHR1) & (df["POS"].astype(float) <= GAP_END)
    if m_a.any():
        a = df.loc[m_a, ["CHROM","POS"] + meta_cols].copy()
        a["CHROM"] = CHR1A
        rows.append(a)
    m_b = (df["CHROM"] == CHR1) & (df["POS"].astype(float) >= OFFSET)
    if m_b.any():
        b = df.loc[m_b, ["CHROM","POS"] + meta_cols].copy()
        b["CHROM"] = CHR1B
        b["POS"] = b["POS"].astype(float) - OFFSET + 1.0
        rows.append(b)
    rest = df.loc[~df["CHROM"].isin([CHR1]), ["CHROM","POS"] + meta_cols]
    if len(rest): rows.append(rest)
    if rows: return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=["CHROM","POS"] + meta_cols)

def map_points_to_LG_with_meta(df_points_std: pd.DataFrame, scaff_info: pd.DataFrame, meta_cols: list[str] | None = None) -> pd.DataFrame:
    if meta_cols is None: meta_cols = []
    if df_points_std is None or df_points_std.empty:
        cols = ["LG","pseudo_pos"] + meta_cols
        return pd.DataFrame({c: pd.Series(dtype=float if c!="LG" else str) for c in cols})
    df = df_points_std.copy()
    df["CHROM"] = df["CHROM"].map(normalize_scaffold_val)
    df["POS"]   = pd.to_numeric(df["POS"], errors="coerce")
    df = df.dropna(subset=["CHROM","POS"])
    merged = merge_1ab_points_df(df[["CHROM","POS"] + meta_cols])
    split  = split_1_to_1ab_points_with_meta(merged, meta_cols=meta_cols)
    if split.empty:
        cols = ["LG","pseudo_pos"] + meta_cols
        return pd.DataFrame({c: pd.Series(dtype=float if c!="LG" else str) for c in cols})
    m = split.merge(scaff_info, on="CHROM", how="inner")
    pseudo = np.where(
        m["Strand"] == "+",
        m["lg_scaf_start"] + m["POS"].astype(float),
        m["lg_scaf_start"] + (m["length"].astype(float) - m["POS"].astype(float) + 1.0),
    )
    out = pd.DataFrame({"LG": m["LG"].astype(str), "pseudo_pos": pseudo})
    for k in meta_cols:
        out[k] = m[k].to_numpy()
    return out

def merge_1ab_intervals_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    orig = out["CHROM"].astype(str)
    out.loc[orig == CHR1B, ["START","END"]] = out.loc[orig == CHR1B, ["START","END"]].astype(float) + (OFFSET - 1.0)
    out.loc[orig.isin([CHR1A, CHR1B, CHR1]), "CHROM"] = CHR1
    return out

def split_1_to_1ab_intervals_with_meta(df: pd.DataFrame, meta_cols: list[str] | None = None) -> pd.DataFrame:
    if meta_cols is None: meta_cols = []
    df = df.copy()
    rows = []
    mask1 = (df["CHROM"] == CHR1) & (df["END"].astype(float) >= 1.0) & (df["START"].astype(float) <= GAP_END)
    if mask1.any():
        a = df.loc[mask1, ["CHROM","START","END"] + meta_cols].copy()
        a["CHROM"] = CHR1A
        a["START"] = np.maximum(a["START"].astype(float), 1.0)
        a["END"]   = np.minimum(a["END"].astype(float), GAP_END)
        rows.append(a)
    mask2 = (df["CHROM"] == CHR1) & (df["END"].astype(float) >= OFFSET)
    if mask2.any():
        b = df.loc[mask2, ["CHROM","START","END"] + meta_cols].copy()
        b["CHROM"] = CHR1B
        b["START"] = np.maximum(b["START"].astype(float), OFFSET) - OFFSET + 1.0
        b["END"]   = b["END"].astype(float) - OFFSET + 1.0
        rows.append(b)
    rest = df.loc[~df["CHROM"].isin([CHR1]) , ["CHROM","START","END"] + meta_cols]
    if len(rest): rows.append(rest)
    if rows:
        out = pd.concat(rows, ignore_index=True)
        out = out[out["END"].astype(float) >= out["START"].astype(float)]
        return out
    return pd.DataFrame(columns=["CHROM","START","END"] + meta_cols)

def map_intervals_to_LG_with_meta(df_intervals: pd.DataFrame, scaff_info: pd.DataFrame, meta_cols: list[str] | None = None) -> pd.DataFrame:
    if meta_cols is None: meta_cols = []
    if df_intervals is None or df_intervals.empty:
        cols = ["LG","pseudo_start","pseudo_end"] + meta_cols
        return pd.DataFrame({c: pd.Series(dtype=float if c!="LG" else str) for c in cols})
    df = df_intervals.copy()
    df["CHROM"] = df["CHROM"].map(normalize_scaffold_val)
    df["START"] = pd.to_numeric(df["START"], errors="coerce")
    df["END"]   = pd.to_numeric(df["END"],   errors="coerce")
    df = df.dropna(subset=["CHROM","START","END"])
    df = df[df["END"] >= df["START"]]
    merged = merge_1ab_intervals_df(df[["CHROM","START","END"] + meta_cols])
    split  = split_1_to_1ab_intervals_with_meta(merged, meta_cols=meta_cols)
    if split.empty:
        cols = ["LG","pseudo_start","pseudo_end"] + meta_cols
        return pd.DataFrame({c: pd.Series(dtype=float if c!="LG" else str) for c in cols})
    m = split.merge(scaff_info, on="CHROM", how="inner")
    ps = np.where(
        m["Strand"] == "+",
        m["lg_scaf_start"] + m["START"].astype(float),
        m["lg_scaf_start"] + (m["length"].astype(float) - m["END"].astype(float) + 1.0),
    )
    pe = np.where(
        m["Strand"] == "+",
        m["lg_scaf_start"] + m["END"].astype(float),
        m["lg_scaf_start"] + (m["length"].astype(float) - m["START"].astype(float) + 1.0),
    )
    ps, pe = np.minimum(ps, pe), np.maximum(ps, pe)
    out = pd.DataFrame({"LG": m["LG"].astype(str), "pseudo_start": ps, "pseudo_end": pe})
    for k in meta_cols:
        out[k] = m[k].to_numpy()
    return out

# ===============================
# Windows & smoothing
# ===============================
def map_points_to_LG(df_points_std: pd.DataFrame, scaff_info: pd.DataFrame) -> pd.DataFrame:
    return map_points_to_LG_with_meta(df_points_std, scaff_info, meta_cols=[])[["LG","pseudo_pos"]]

def build_snp_windows_LG(snp_std: pd.DataFrame, scaff_info: pd.DataFrame, window_size=1_000_000) -> pd.DataFrame:
    pts = map_points_to_LG(snp_std, scaff_info)
    if pts.empty:
        return pd.DataFrame(columns=["LG","start","end","n_snps"])
    lg_lengths = scaff_info.groupby("LG", as_index=False)["lg_scaf_end"].max().rename(columns={"lg_scaf_end":"lg_length"})
    rows = []
    for lg, sub in pts.groupby("LG"):
        L = float(lg_lengths.loc[lg_lengths["LG"]==lg, "lg_length"].values[0])
        idx = np.floor(sub["pseudo_pos"].values / window_size).astype(int)
        unique, counts = np.unique(idx, return_counts=True)
        nwin = int(math.floor(L / window_size) + 1)
        cmap = dict(zip(unique, counts))
        for w in range(nwin):
            start = w * window_size
            end   = min((w+1)*window_size, L)
            n = int(cmap.get(w, 0))
            rows.append((lg, start, end, n))
    return pd.DataFrame(rows, columns=["LG","start","end","n_snps"])

def smooth_tmrca_LG(tm_std: pd.DataFrame, scaff_info: pd.DataFrame, win_bp=500_000, step_bp=None) -> pd.DataFrame:
    if tm_std is None or tm_std.empty:
        return pd.DataFrame(columns=["LG","win_start","win_end","mid","mean_s","lo_s","hi_s"])
    if step_bp is None:
        step_bp = max(1, int(win_bp//5))
    segs_map = map_intervals_to_LG_with_meta(tm_std[["CHROM","START","END","mean","lo","hi"]], scaff_info, meta_cols=["mean","lo","hi"])
    if segs_map.empty:
        return pd.DataFrame(columns=["LG","win_start","win_end","mid","mean_s","lo_s","hi_s"])
    lg_lengths = scaff_info.groupby("LG", as_index=False)["lg_scaf_end"].max().rename(columns={"lg_scaf_end":"lg_length"})
    out = []
    for lg, g in segs_map.groupby("LG"):
        L = float(lg_lengths.loc[lg_lengths["LG"]==lg, "lg_length"].values[0])
        starts = np.arange(0, max(0, L - win_bp + 1), step_bp, dtype=float)
        ends   = np.minimum(starts + win_bp, L)
        s = g["pseudo_start"].values.astype(float)
        e = g["pseudo_end"].values.astype(float)
        mn = g["mean"].values.astype(float)
        lo = g["lo"].values.astype(float)
        hi = g["hi"].values.astype(float)
        for ws, we in zip(starts, ends):
            ol = np.minimum(we, e) - np.maximum(ws, s)
            mask = ol > 0
            if not np.any(mask):
                continue
            w = ol[mask]
            mean_s = np.sum(mn[mask] * w) / np.sum(w)
            lo_s   = np.sum(lo[mask] * w) / np.sum(w)
            hi_s   = np.sum(hi[mask] * w) / np.sum(w)
            out.append((lg, ws, we, (ws+we)/2.0, mean_s, lo_s, hi_s))
    return pd.DataFrame(out, columns=["LG","win_start","win_end","mid","mean_s","lo_s","hi_s"])

# ===============================
# Load inputs
# ===============================
tmrca_df   = load_tmrca_std(DATA["tmrca"])
snp_df     = load_snp_std(DATA["snp"])
te_df      = load_te_std(DATA["te"])
mrna_df    = load_mrna_std(DATA["mrna"])
gwas_orig  = load_gwas_orig_std(DATA["gwas_orig"])
trait_df   = load_trait_gwas_std(DATA["trait_gwas"])
scaff_info = load_lg_map(DATA["lg_map"])

if tmrca_df is None or snp_df is None or scaff_info is None or tmrca_df.empty or snp_df.empty or scaff_info.empty:
    raise SystemExit("Required files missing or empty (TMRCA, SNP, LG map).")

# LG order: strictly numeric where possible
def _numkey(s):
    try: return (0, int(float(s)))
    except: return (1, str(s))
lg_order = sorted(scaff_info["LG"].astype(str).unique().tolist(), key=_numkey)

# Make ordered sectors dict
lg_lengths = scaff_info.groupby("LG", as_index=False)["lg_scaf_end"].max().rename(columns={"lg_scaf_end":"lg_length"})
sectors = {lg: float(lg_lengths.loc[lg_lengths["LG"]==lg, "lg_length"].values[0]) for lg in lg_order}

# Precompute layers
snp_win = build_snp_windows_LG(snp_df, scaff_info, window_size=SNP_WINDOW_BP) if SHOW_SNP_DENSITY else pd.DataFrame()

present_te_labels = set()
te_map = pd.DataFrame()
if SHOW_TE and te_df is not None and not te_df.empty:
    te_keep = te_df.copy()
    te_keep = te_keep[te_keep["MAJOR"].isin(TE_MAJOR_TO_SHOW)]
    te_keep["DISPLAY"] = np.where(
        (te_keep["MAJOR"]=="LTR") & (te_keep["SUB"].isin(["Gypsy","Copia"])),
        te_keep["MAJOR"] + ":" + te_keep["SUB"],
        te_keep["MAJOR"]
    )
    te_keep = te_keep.dropna(subset=["DISPLAY"])
    if not te_keep.empty:
        te_map = map_intervals_to_LG_with_meta(te_keep[["CHROM","START","END","DISPLAY"]], scaff_info, meta_cols=["DISPLAY"])
        present_te_labels = set(te_map["DISPLAY"].unique().tolist())

mrna_map = pd.DataFrame()
if SHOW_MRNA and mrna_df is not None and not mrna_df.empty:
    mrna_map = map_intervals_to_LG_with_meta(mrna_df[["CHROM","START","END"]], scaff_info, meta_cols=[])

tm_smooth = smooth_tmrca_LG(tmrca_df, scaff_info, win_bp=TMRCA_WIN_BP, step_bp=TMRCA_STEP_BP) if SHOW_TMRCA else pd.DataFrame()

trait_map = pd.DataFrame()
if SHOW_TRAIT_GWAS and trait_df is not None and not trait_df.empty:
    td = trait_df.copy()
    td = td[(td["MAF"].isna()) | (td["MAF"] >= MAF_CUTOFF)]
    pts_with_meta = td[["CHROM","POS","FDR"]].copy()
    trait_map = map_points_to_LG_with_meta(pts_with_meta, scaff_info, meta_cols=["FDR"])
    trait_map["FDR"] = pd.to_numeric(trait_map["FDR"], errors="coerce").clip(lower=np.finfo(float).tiny, upper=1.0)
    trait_map = trait_map.dropna(subset=["FDR"])

orig_map = pd.DataFrame()
if SHOW_ORIG_GWAS and gwas_orig is not None and not gwas_orig.empty:
    orig_map = map_points_to_LG(gwas_orig[["CHROM","POS"]], scaff_info)

# ===============================
# Build Circos
# ===============================
circos = Circos(sectors, space=4)  # gap in degrees
# Make name→sector map (FIX)
sector_map = {s.name: s for s in circos.sectors}

# Outer thin axis + LG labels
for s in circos.sectors:
    s.axis(fc="none", ec="black", lw=1.0, alpha=0.6)
    s.text(f"LG{s.name}", r=105, size=10)

# 1) SNP density
if SHOW_SNP_DENSITY and not snp_win.empty:
    max_n = max(1, int(snp_win["n_snps"].max()))
    cmap = LinearSegmentedColormap.from_list("snp_heat",
        ["#ffffff","#006400","#ffa500","#ff8c00","#ff0000","#8b0000"], N=256)
    norm = Normalize(vmin=0, vmax=max_n)
    for lg in lg_order:
        lg_sub = snp_win[snp_win["LG"] == lg]
        if lg_sub.empty: continue
        tr = sector_map[lg].add_track((80, 90))
        tr.axis()
        for _, row in lg_sub.iterrows():
            x1, x2 = float(row["start"]), float(row["end"])
            color = cmap(norm(int(row["n_snps"])))
            tr.rect(x1, x2, ec=None, color=color)

# 2) TE track
if SHOW_TE and not te_map.empty:
    for lg in lg_order:
        lg_sub = te_map[te_map["LG"] == lg]
        if lg_sub.empty: continue
        tr = sector_map[lg].add_track((75, 80))
        tr.axis()
        for _, row in lg_sub.iterrows():
            x1, x2 = float(row["pseudo_start"]), float(row["pseudo_end"])
            disp = str(row["DISPLAY"])
            color = TE_COLORS.get(disp, "#999999")
            tr.rect(x1, x2, ec=None, color=color)

# 3) mRNA track
if SHOW_MRNA and not mrna_map.empty:
    for lg in lg_order:
        lg_sub = mrna_map[mrna_map["LG"] == lg]
        if lg_sub.empty: continue
        tr = sector_map[lg].add_track((70, 75))
        tr.axis()
        for _, row in lg_sub.iterrows():
            tr.rect(float(row["pseudo_start"]), float(row["pseudo_end"]), ec=None, color="#F17CB0")

# 4) TMRCA (CI + mean)
if SHOW_TMRCA and not tm_smooth.empty:
    for lg in lg_order:
        sub = tm_smooth[tm_smooth["LG"] == lg]
        if sub.empty: continue
        ylo = float(np.nanmin(sub[["lo_s","hi_s"]].values))
        yhi = float(np.nanmax(sub[["lo_s","hi_s"]].values))
        if not np.isfinite(ylo) or ylo == yhi:
            ylo, yhi = 1e-6, 10.0
        tr = sector_map[lg].add_track((55, 70), r_pad_ratio=0.05)
        tr.axis()
        x = sub["mid"].values.astype(float)
        lo = sub["lo_s"].values.astype(float)
        hi = sub["hi_s"].values.astype(float)
        mn = sub["mean_s"].values.astype(float)
        tr.fill_between(x, lo, y2=hi, fc="#9ecae1", ec="none", alpha=0.35)
        tr.line(x, mn, color="#08519c", lw=1.8)

# 5) Trait GWAS (−log10 FDR)
if SHOW_TRAIT_GWAS and not trait_map.empty:
    trait_map = trait_map.copy()
    trait_map["neglog10"] = -np.log10(trait_map["FDR"].astype(float))
    ymax = float(np.nanmax(trait_map["neglog10"].values))
    if not np.isfinite(ymax) or ymax <= 0: ymax = 1.0
    for lg in lg_order:
        sub = trait_map[trait_map["LG"] == lg]
        if sub.empty: continue
        x = sub["pseudo_pos"].values.astype(float)
        y = sub["neglog10"].values.astype(float)
        a = np.where(sub["FDR"].astype(float) <= FDR_CUTOFF, 1.0, 0.5)
        colors = [(0.415,0.318,0.639,alpha) for alpha in a]
        tr = sector_map[lg].add_track((45, 55), r_pad_ratio=0.05)
        tr.axis()
        tr.scatter(x, y, s=TRAIT_POINT_SIZE, color=colors)

# 6) Original GWAS (red)
if SHOW_ORIG_GWAS and not orig_map.empty:
    for lg in lg_order:
        sub = orig_map[orig_map["LG"] == lg]
        if sub.empty: continue
        x = sub["pseudo_pos"].values.astype(float)
        y = np.full_like(x, 0.5, dtype=float)
        tr = sector_map[lg].add_track((42, 45))
        tr.axis()
        tr.scatter(x, y, s=ORIG_GWAS_SIZE, color="red")

# Title & TE legend
fig = circos.plotfig()
fig.suptitle("Circos (Python): SNP density, TE majors (Gypsy/Copia), mRNA, TMRCA, GWAS", y=0.98, fontsize=12)

present_te_labels = [k for k in ["LTR:Gypsy","LTR:Copia","DNA","Helitron","LINE","SINE"]
                     if k in present_te_labels]
if present_te_labels:
    handles, labels = [], []
    for key in present_te_labels:
        patch = plt.Line2D([0],[0], marker="s", color="w",
                           markerfacecolor=TE_COLORS.get(key,"#999999"),
                           markersize=10, linestyle="None")
        handles.append(patch); labels.append(key)
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 6),
               frameon=False, bbox_to_anchor=(0.5, 0.02))

# Save
if OUT_JPEG: fig.savefig(OUT_JPEG, dpi=DPI, bbox_inches="tight")
if OUT_PDF:  fig.savefig(OUT_PDF, dpi=DPI, bbox_inches="tight")
print(f"Saved: {OUT_JPEG} {('and ' + OUT_PDF) if OUT_PDF else ''}")

