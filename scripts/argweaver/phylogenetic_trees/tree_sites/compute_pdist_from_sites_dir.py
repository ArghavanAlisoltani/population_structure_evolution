#!/usr/bin/env python3
import argparse, glob, os
import numpy as np
import multiprocessing as mp

MISSING_BYTES = np.frombuffer(b"Nn-.?Xx", dtype=np.uint8)  # treat these as missing

def parse_names_from_sites(path: str):
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("#NAMES\t"):
                return line.rstrip("\n").split("\t")[1:]
    raise ValueError(f"No #NAMES line found in: {path}")

def update_counts(valid, mismatch, allele_strings, n):
    joined = "".join(allele_strings).encode("ascii", errors="ignore")
    arr = np.frombuffer(joined, dtype=np.uint8)
    if arr.size != len(allele_strings) * n:
        raise ValueError("Unexpected size while converting allele strings to matrix.")
    codes = arr.reshape((len(allele_strings), n)).copy()

    # set missing chars to 0
    mask_missing = np.isin(codes, MISSING_BYTES)
    codes[mask_missing] = 0

    V = (codes != 0).astype(np.uint8)
    valid_chunk = V.T @ V

    # matches for each observed state in this chunk (excluding 0)
    uniq = np.unique(codes)
    uniq = uniq[uniq != 0]
    matches = np.zeros((n, n), dtype=np.int64)
    for code in uniq:
        M = (codes == code).astype(np.uint8)
        matches += (M.T @ M)

    mismatch_chunk = valid_chunk - matches
    valid += valid_chunk
    mismatch += mismatch_chunk
    return valid, mismatch

def process_files(file_list, names_ref, chunk_sites=5000, max_sites_per_file=0):
    n = len(names_ref)
    valid = np.zeros((n, n), dtype=np.int64)
    mismatch = np.zeros((n, n), dtype=np.int64)

    for path in file_list:
        names = parse_names_from_sites(path)
        if names != names_ref:
            raise ValueError(f"#NAMES mismatch in {path}")

        chunk = []
        seen = 0
        with open(path, "rt") as f:
            for line in f:
                if not line or line[0] == "#":
                    continue
                line = line.rstrip("\n")
                if not line:
                    continue
                try:
                    _, s = line.split("\t", 1)
                except ValueError:
                    continue
                if len(s) != n:
                    raise ValueError(f"Allele string length != n in {path}: {len(s)} vs {n}")
                chunk.append(s)
                seen += 1
                if max_sites_per_file > 0 and seen >= max_sites_per_file:
                    break
                if len(chunk) >= chunk_sites:
                    valid, mismatch = update_counts(valid, mismatch, chunk, n)
                    chunk = []
        if chunk:
            valid, mismatch = update_counts(valid, mismatch, chunk, n)

    return valid, mismatch

def worker(args_tuple):
    files, names_ref, chunk_sites, max_sites_per_file = args_tuple
    return process_files(files, names_ref, chunk_sites, max_sites_per_file)

def chunk_list(items, n_chunks):
    n_chunks = max(1, n_chunks)
    chunks = [[] for _ in range(n_chunks)]
    for i, x in enumerate(items):
        chunks[i % n_chunks].append(x)
    return [c for c in chunks if c]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites_dir", required=True)
    ap.add_argument("--pattern", default="*.sites")
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--outprefix", default="ALL_SITES")
    ap.add_argument("--n_jobs", type=int, default=4)
    ap.add_argument("--chunk_sites", type=int, default=5000)
    ap.add_argument("--max_sites_per_file", type=int, default=0)
    ap.add_argument("--mp_start", default="spawn", choices=["spawn","fork","forkserver"])
    args = ap.parse_args()

    if args.recursive:
        files = glob.glob(os.path.join(args.sites_dir, "**", args.pattern), recursive=True)
    else:
        files = glob.glob(os.path.join(args.sites_dir, args.pattern))
    files = sorted(files)
    if not files:
        raise SystemExit("No .sites files found.")

    names_ref = parse_names_from_sites(files[0])
    n = len(names_ref)
    print(f"Found {len(files)} files; n_samples={n}")

    batches = chunk_list(files, args.n_jobs)
    work = [(b, names_ref, args.chunk_sites, args.max_sites_per_file) for b in batches]

    if args.n_jobs > 1 and len(work) > 1:
        ctx = mp.get_context(args.mp_start)
        with ctx.Pool(args.n_jobs) as pool:
            results = pool.map(worker, work)
    else:
        results = [worker(work[0])]

    valid = np.zeros((n, n), dtype=np.int64)
    mismatch = np.zeros((n, n), dtype=np.int64)
    for v, m in results:
        valid += v
        mismatch += m

    # p-distance (robust for 0/1 or A/C/G/T, etc.)
    valid_f = valid.astype(np.float64)
    mismatch_f = mismatch.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(valid_f > 0, mismatch_f / valid_f, np.nan)

    # IMPORTANT: do NOT set valid==0 pairs to 0 (that makes fake “identical” samples).
    # Fill missing overlaps with the mean of observed distances (conservative).
    mean_p = np.nanmean(p)
    p = np.where(np.isnan(p), mean_p, p)
    np.fill_diagonal(p, 0.0)

    dist_path = args.outprefix + ".pdist.tsv"
    ovl_path  = args.outprefix + ".overlap_sites.tsv"

    with open(dist_path, "wt") as f:
        f.write("\t" + "\t".join(names_ref) + "\n")
        for i, nm in enumerate(names_ref):
            f.write(nm + "\t" + "\t".join(f"{p[i,j]:.8f}" for j in range(n)) + "\n")

    with open(ovl_path, "wt") as f:
        f.write("\t" + "\t".join(names_ref) + "\n")
        for i, nm in enumerate(names_ref):
            f.write(nm + "\t" + "\t".join(str(int(valid[i,j])) for j in range(n)) + "\n")

    print("Wrote:")
    print(" ", dist_path)
    print(" ", ovl_path)

if __name__ == "__main__":
    main()

