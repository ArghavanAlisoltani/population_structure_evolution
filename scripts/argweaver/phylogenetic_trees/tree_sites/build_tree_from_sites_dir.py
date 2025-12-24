#!/usr/bin/env python3
import argparse, glob, os
import numpy as np
import multiprocessing as mp

# ----------------------------
# Parsing
# ----------------------------
def parse_names_from_sites(path: str):
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("#NAMES\t"):
                parts = line.rstrip("\n").split("\t")
                return parts[1:]
    raise ValueError(f"No #NAMES line found in: {path}")

def lut_table():
    lut = np.zeros(256, dtype=np.uint8)
    for ch in (b"A", b"a"): lut[ch[0]] = 1
    for ch in (b"C", b"c"): lut[ch[0]] = 2
    for ch in (b"G", b"g"): lut[ch[0]] = 3
    for ch in (b"T", b"t"): lut[ch[0]] = 4
    # all others => 0 (missing)
    return lut

LUT = lut_table()

# ----------------------------
# Core counting (pairwise valid/mismatch) without building full alignment
# ----------------------------
def update_counts(valid, mismatch, allele_strings, n):
    # Convert chunk of strings into (m_sites x n_samples) code matrix
    joined = "".join(allele_strings).encode("ascii", errors="ignore")
    arr = np.frombuffer(joined, dtype=np.uint8)
    if arr.size != len(allele_strings) * n:
        raise ValueError("Unexpected encoding/size while converting allele strings to matrix.")
    codes = LUT[arr].reshape((len(allele_strings), n))

    V = (codes != 0).astype(np.uint8)        # valid allele (A/C/G/T)
    valid_chunk = V.T @ V                    # n x n counts of jointly valid sites

    matches = np.zeros((n, n), dtype=np.int64)
    for code in (1, 2, 3, 4):                # A,C,G,T
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
            raise ValueError(
                f"#NAMES mismatch in {path}\n"
                f"Expected first file names (n={len(names_ref)}), got (n={len(names)})"
            )

        chunk = []
        seen_sites = 0

        with open(path, "rt") as f:
            for line in f:
                if not line or line[0] == "#":
                    continue
                line = line.rstrip("\n")
                if not line:
                    continue

                try:
                    _, s = line.split("\t", 1)   # pos \t allele_string
                except ValueError:
                    continue

                if len(s) != n:
                    raise ValueError(f"Allele string length != n in {path}. Got {len(s)} expected {n}")

                chunk.append(s)
                seen_sites += 1
                if max_sites_per_file > 0 and seen_sites >= max_sites_per_file:
                    break

                if len(chunk) >= chunk_sites:
                    valid, mismatch = update_counts(valid, mismatch, chunk, n)
                    chunk = []

        if chunk:
            valid, mismatch = update_counts(valid, mismatch, chunk, n)

    return valid, mismatch

# ----------------------------
# JC69 distances
# ----------------------------
def jc69_distance(valid, mismatch):
    valid = valid.astype(np.float64)
    mismatch = mismatch.astype(np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(valid > 0, mismatch / valid, np.nan)

    p = np.clip(p, 0.0, 0.749999)  # avoid log(<=0)
    d = -0.75 * np.log(1.0 - (4.0/3.0) * p)

    d[np.isnan(d)] = 0.0
    np.fill_diagonal(d, 0.0)
    return d

# ----------------------------
# Neighbor-Joining (NJ) to Newick
# ----------------------------
class Node:
    __slots__ = ("name", "children")
    def __init__(self, name):
        self.name = name
        self.children = []  # list of (Node, branch_length)

def to_newick(node):
    if not node.children:
        return node.name
    parts = []
    for child, bl in node.children:
        bl = max(0.0, float(bl))
        parts.append(f"{to_newick(child)}:{bl:.8f}")
    return "(" + ",".join(parts) + ")"

def neighbor_joining(dist_mat, names):
    D = dist_mat.copy()
    labels = [Node(n) for n in names]
    active = list(range(len(names)))
    next_id = 1

    while len(active) > 2:
        m = len(active)
        r = np.array([D[i, active].sum() for i in active], dtype=np.float64)

        Qmin = None
        pair = None
        for a in range(m):
            i = active[a]
            for b in range(a+1, m):
                j = active[b]
                q = (m - 2) * D[i, j] - r[a] - r[b]
                if Qmin is None or q < Qmin:
                    Qmin = q
                    pair = (a, b)

        a, b = pair
        i = active[a]
        j = active[b]

        delta = (r[a] - r[b]) / (m - 2)
        li = 0.5 * (D[i, j] + delta)
        lj = 0.5 * (D[i, j] - delta)

        u = Node(f"U{next_id}")
        next_id += 1
        u.children.append((labels[i], li))
        u.children.append((labels[j], lj))

        new_row = {}
        for k in active:
            if k in (i, j):
                continue
            duk = 0.5 * (D[i, k] + D[j, k] - D[i, j])
            new_row[k] = duk

        labels.append(u)
        u_idx = len(labels) - 1

        D2 = np.zeros((D.shape[0] + 1, D.shape[1] + 1), dtype=np.float64)
        D2[:-1, :-1] = D
        D = D2
        for k, duk in new_row.items():
            D[u_idx, k] = duk
            D[k, u_idx] = duk

        active = [x for x in active if x not in (i, j)]
        active.append(u_idx)

    i, j = active[0], active[1]
    root = Node("ROOT")
    root.children.append((labels[i], D[i, j] / 2.0))
    root.children.append((labels[j], D[i, j] / 2.0))
    return root

# ----------------------------
# Multiprocessing worker (TOP-LEVEL => picklable)
# ----------------------------
def worker_process(args_tuple):
    file_batch, names_ref, chunk_sites, max_sites_per_file = args_tuple
    return process_files(file_batch, names_ref, chunk_sites=chunk_sites, max_sites_per_file=max_sites_per_file)

def chunk_list(items, n_chunks):
    n_chunks = max(1, n_chunks)
    chunks = [[] for _ in range(n_chunks)]
    for idx, x in enumerate(items):
        chunks[idx % n_chunks].append(x)
    return [c for c in chunks if c]

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites_dir", required=True)
    ap.add_argument("--pattern", default="*.sites")
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--outprefix", default="all_sites_tree")
    ap.add_argument("--n_jobs", type=int, default=4)
    ap.add_argument("--chunk_sites", type=int, default=5000)
    ap.add_argument("--max_sites_per_file", type=int, default=0)
    ap.add_argument("--mp_start", default="spawn", choices=["spawn","fork","forkserver"],
                    help="Mac default is spawn. Try fork on some systems for speed.")
    args = ap.parse_args()

    if args.recursive:
        files = glob.glob(os.path.join(args.sites_dir, "**", args.pattern), recursive=True)
    else:
        files = glob.glob(os.path.join(args.sites_dir, args.pattern))
    files = sorted(files)
    if not files:
        raise SystemExit("No .sites files found with given dir/pattern.")

    names_ref = parse_names_from_sites(files[0])
    n = len(names_ref)
    print(f"Found {len(files)} files; n_samples={n}")

    batches = chunk_list(files, args.n_jobs)
    work = [(b, names_ref, args.chunk_sites, args.max_sites_per_file) for b in batches]

    if args.n_jobs > 1 and len(work) > 1:
        ctx = mp.get_context(args.mp_start)
        with ctx.Pool(processes=args.n_jobs) as pool:
            results = pool.map(worker_process, work)
    else:
        results = [worker_process(work[0])]

    valid = np.zeros((n, n), dtype=np.int64)
    mismatch = np.zeros((n, n), dtype=np.int64)
    for v, m in results:
        valid += v
        mismatch += m

    dist = jc69_distance(valid, mismatch)
    root = neighbor_joining(dist, names_ref)
    newick = to_newick(root) + ";\n"

    nwk_path = args.outprefix + ".nwk"
    with open(nwk_path, "wt") as f:
        f.write(newick)

    dist_path = args.outprefix + ".dist.tsv"
    with open(dist_path, "wt") as f:
        f.write("\t" + "\t".join(names_ref) + "\n")
        for i, nm in enumerate(names_ref):
            row = "\t".join(f"{dist[i,j]:.8f}" for j in range(n))
            f.write(nm + "\t" + row + "\n")

    print("Wrote:")
    print(" ", nwk_path)
    print(" ", dist_path)

if __name__ == "__main__":
    main()

