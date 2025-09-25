# make_balanced_windows.R
# Input: 2-column table of variant positions with a header (e.g., CHROM, POS)
# Output: TSV with columns: Scaffold  start  end  nsnp

setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/balanced_100MB_window/")
# -------- user settings --------
infile      <- "~/Desktop/OSU_projects/conifers/LP/vcf_v1/positions_split_poly_s100_scaffolds.tsv"  # change if needed
outfile     <- "balanced_scaffold_windows_100Mb.tsv"
target_size <- 1e8      # ~100,000,000 bp per window
w_len       <- 0.5      # weight for length balance (0..1)
w_snp       <- 0.5      # weight for SNP-count balance (0..1)
# --------------------------------

stopifnot(abs(w_len + w_snp - 1) < 1e-9)

# --- read input (robust to column names) ---
suppressWarnings({
  d <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
})
if (!all(c("CHROM","POS") %in% names(d))) {
  # fallbacks: try common alternatives
  cn <- names(d)
  sc_col <- if ("Scaffold" %in% cn) "Scaffold" else cn[1]
  pos_col <- if ("Position" %in% cn) "Position" else if ("POS" %in% cn) "POS" else cn[2]
  d <- setNames(d[, c(sc_col, pos_col)], c("CHROM","POS"))
}
d$POS <- as.numeric(d$POS)
d <- d[is.finite(d$POS), ]
stopifnot(nrow(d) > 0)

# --- helpers ---
# Candidate boundaries: start(1), midpoints between consecutive SNPs, end(L=max POS).
# Also track cumulative SNP count to the LEFT of each boundary.
build_candidates <- function(pos_vec, L) {
  p <- sort(unique(pos_vec))
  n <- length(p)
  if (n <= 1) {
    # No or one SNP: only start and end boundaries; cum_snp = c(0, n)
    return(list(pos = c(1L, L), cum_snp = c(0L, n)))
  }
  mids <- floor((p[-n] + p[-1]) / 2)
  cand_pos <- c(1L, mids, L)
  # cum_snp: number of SNPs strictly to the left of each boundary
  cum_snp <- c(0L, seq_len(n - 1), n)
  list(pos = as.numeric(cand_pos), cum_snp = as.numeric(cum_snp))
}

# Greedy multi-objective cutpoint picker.
# Chooses k-1 increasing indices into candidate arrays to balance length and SNP count.
pick_cuts <- function(cand_pos, cum_snp, L, N, k, w_len, w_snp) {
  if (k <= 1) return(integer(0))
  cuts <- integer(k - 1)
  min_idx <- 2L                               # cannot cut at first boundary (start)
  max_idx_global <- length(cand_pos) - 1L     # cannot cut at last boundary (end)
  
  for (j in seq_len(k - 1)) {
    t_len <- j * L / k                        # target length for boundary j
    t_snp <- j * N / k                        # target SNP count for boundary j
    
    remain <- (k - 1) - (j - 1)
    max_idx <- max_idx_global - (remain - 1)  # leave room for remaining cuts
    
    if (min_idx > max_idx) {
      # no space left; bail
      return(integer(0))
    }
    idx <- seq.int(min_idx, max_idx)
    
    # normalized costs
    cost_len <- abs(cand_pos[idx] - t_len) / max(L, 1)
    cost_snp <- if (N > 0) abs(cum_snp[idx] - t_snp) / N else 0
    cost <- w_len * cost_len + w_snp * cost_snp
    
    best_local <- idx[which.min(cost)]
    cuts[j] <- best_local
    min_idx <- best_local + 1L
  }
  cuts
}

# Count SNPs per window using cumulative counts at boundaries
window_nsnp <- function(cum_snp, cut_indices) {
  all_idx <- c(1L, cut_indices, length(cum_snp))
  k <- length(all_idx) - 1L
  ns <- integer(k)
  for (i in seq_len(k)) {
    left_idx <- all_idx[i]
    right_idx <- all_idx[i + 1]
    ns[i] <- cum_snp[right_idx] - cum_snp[left_idx]
  }
  ns
}

# Build windows for a single scaffold
make_windows_one_scaffold <- function(scaf, pos_vec, target_size, w_len, w_snp) {
  p <- sort(unique(pos_vec))
  L <- max(p, na.rm = TRUE)
  N <- length(p)
  
  k <- max(1L, as.integer(ceiling(L / target_size)))   # number of windows
  
  cand <- build_candidates(p, L)
  cand_pos <- cand$pos
  cum_snp <- cand$cum_snp
  
  if (k == 1) {
    return(data.frame(Scaffold = scaf, start = 1L, end = L, nsnp = N, stringsAsFactors = FALSE))
  }
  
  cut_idx <- pick_cuts(cand_pos, cum_snp, L, N, k, w_len, w_snp)
  
  # Fallback if greedy failed (rare): snap equal-length targets to nearest candidates, enforcing order
  if (length(cut_idx) != (k - 1)) {
    targets <- (seq_len(k - 1) * L) / k
    cut_idx <- integer(k - 1)
    last <- 1L
    for (j in seq_len(k - 1)) {
      # ensure increasing indices and leave room for remaining cuts
      max_allow <- length(cand_pos) - (k - j)
      idx <- seq.int(last + 1L, max_allow)
      pick <- idx[which.min(abs(cand_pos[idx] - targets[j]))]
      cut_idx[j] <- pick
      last <- pick
    }
  }
  
  bounds <- c(1L, cand_pos[cut_idx], L)
  ns <- window_nsnp(cum_snp, cut_idx)
  
  data.frame(
    Scaffold = scaf,
    start = bounds[-length(bounds)],
    end   = bounds[-1],
    nsnp  = ns,
    stringsAsFactors = FALSE
  )
}

# --- run per scaffold ---
by_scaf <- split(d$POS, d$CHROM)
out_list <- vector("list", length(by_scaf))
nm <- names(by_scaf)
for (i in seq_along(by_scaf)) {
  out_list[[i]] <- make_windows_one_scaffold(nm[i], by_scaf[[i]], target_size, w_len, w_snp)
}
out <- do.call(rbind, out_list)
out$length<-(out$end-out$start)+1
out$ratio<-(out$nsnp/out$length)*1000
boxplot(out$ratio)
scaff1aa<-out[out$Scaffold=="scaffold_1aa",]
boxplot(scaff1aa$ratio)
# --- write TSV ---
write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat(sprintf("Wrote %d windows across %d scaffolds to %s\n", nrow(out), length(by_scaf), outfile))
