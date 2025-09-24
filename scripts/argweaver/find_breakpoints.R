# find_breakpoints.R
# Input: 2-column TSV from your VCF (#CHROM, POS), header present
setwd("~/path/to/folder/")
# ======= User settings =======
infile        <- "positions_poly_s100_All_1a1b_renamed.txt"  # 2 cols: CHROM, POS (header required)
qc_out        <- "breakpoints.qc.tsv"                        # detailed, per-scaffold QC table
bp_out        <- "breakpoints.txt"                           # 2 cols for awk: scaffold<TAB>breakpoint
LEN_THRESHOLD <- 1e9                                         # only split scaffolds longer than this
MAX_IMBALANCE <- 0.30                                        # require |len_a - len_b| / total_len <= this
MIN_GAP       <- 1                                           # minimum inter-SNP gap size (bp) to consider

# Optional: true scaffold lengths (override max POS). Set to NULL to skip.
# File format: 2 cols, tab-separated: scaffold <TAB> length
lengths_file <- NULL  # e.g., "scaffold_lengths.txt" or keep NULL to skip

# ======= Load input =======
suppressWarnings({
  d <- read.table(infile, header = TRUE, sep = "\t",
                  comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
})
# normalize expected columns
stopifnot(all(c("CHROM","POS") %in% names(d)))
d$POS <- as.numeric(d$POS)

# Optional true lengths
true_len <- NULL
if (!is.null(lengths_file)) {
  tl <- read.table(lengths_file, header = FALSE, sep = "\t",
                   col.names = c("scaffold","length"),
                   stringsAsFactors = FALSE)
  tl$length <- as.numeric(tl$length)
  true_len <- setNames(tl$length, tl$scaffold)
}

# ======= Core selector (Option A) =======
# Enforce a maximum imbalance; among those meeting it, pick the largest gap.
# If none meet it, fallback to the most balanced (smallest imbalance) gap.
pick_breakpoint_optionA <- function(xpos, total_len, max_imbalance = 0.30, min_gap = 1) {
  xpos <- sort(unique(xpos))
  if (length(xpos) < 2) {
    return(list(ok = FALSE, reason = "insufficient_snps"))
  }
  
  left  <- xpos[-length(xpos)]
  right <- xpos[-1]
  gaps  <- right - left
  
  # Filter by minimum gap size
  keep <- which(gaps >= min_gap)
  if (length(keep) == 0) {
    return(list(ok = FALSE, reason = "no_gap_ge_MIN_GAP"))
  }
  left  <- left[keep]; right <- right[keep]; gaps <- gaps[keep]
  
  # Candidate midpoints at gap mid, integer midpoint
  mids <- left + floor(gaps / 2)
  
  # Total length baseline
  L <- if (is.finite(total_len) && total_len > 0) total_len else max(xpos, na.rm = TRUE)
  
  # New scaffold lengths (1-based original coordinates):
  #   a: [1 .. mid-1],  b: [mid .. L] (b will be shifted so mid -> 1)
  len_a <- pmax(mids - 1, 0)
  len_b <- pmax(L - mids + 1, 0)
  
  # Imbalance metric (smaller is better)
  imb <- abs(len_a - len_b) / L
  
  # Keep candidates that meet the balance constraint
  ok_idx <- which(imb <= max_imbalance)
  
  if (length(ok_idx) > 0) {
    # Among balanced candidates, pick the largest gap
    best <- ok_idx[ which.max(gaps[ok_idx]) ]
    reason <- "meets_balance_max_gap"
    satisfied <- TRUE
  } else {
    # Fallback: pick the smallest imbalance overall
    best <- which.min(imb)
    reason <- "fallback_most_balanced"
    satisfied <- FALSE
  }
  
  list(
    ok         = TRUE,
    reason     = reason,
    satisfied  = satisfied,
    breakpoint = as.numeric(mids[best]),
    left_snp   = as.numeric(left[best]),
    right_snp  = as.numeric(right[best]),
    gap        = as.numeric(gaps[best]),
    dist_left  = as.numeric(mids[best] - left[best]),
    dist_right = as.numeric(right[best] - mids[best]),
    len_a      = as.numeric(len_a[best]),
    len_b      = as.numeric(len_b[best]),
    imbalance  = as.numeric(imb[best])
  )
}

# ======= Run per scaffold (only long ones) =======
by_scaf  <- split(d$POS, d$CHROM)
is_long  <- vapply(by_scaf, function(v) max(v, na.rm = TRUE) > LEN_THRESHOLD, logical(1))
targets  <- names(by_scaf)[is_long]

qc_rows <- lapply(targets, function(s) {
  xpos <- by_scaf[[s]]
  
  # Prefer true length if provided; else estimate from max observed POS
  total_len <- if (!is.null(true_len) && !is.na(true_len[s])) true_len[s] else max(xpos, na.rm = TRUE)
  
  res <- pick_breakpoint_optionA(xpos, total_len,
                                 max_imbalance = MAX_IMBALANCE,
                                 min_gap = MIN_GAP)
  
  if (!isTRUE(res$ok)) {
    data.frame(
      scaffold    = s,
      n_snps      = length(unique(xpos)),
      total_len   = total_len,
      max_pos     = max(xpos, na.rm = TRUE),
      breakpoint  = NA_real_,
      left_snp    = NA_real_,
      right_snp   = NA_real_,
      gap         = NA_real_,
      dist_left   = NA_real_,
      dist_right  = NA_real_,
      len_a       = NA_real_,
      len_b       = NA_real_,
      imbalance   = NA_real_,
      satisfied   = FALSE,
      reason      = res$reason,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      scaffold    = s,
      n_snps      = length(unique(xpos)),
      total_len   = total_len,
      max_pos     = max(xpos, na.rm = TRUE),
      breakpoint  = res$breakpoint,
      left_snp    = res$left_snp,
      right_snp   = res$right_snp,
      gap         = res$gap,
      dist_left   = res$dist_left,
      dist_right  = res$dist_right,
      len_a       = res$len_a,
      len_b       = res$len_b,
      imbalance   = res$imbalance,
      satisfied   = res$satisfied,
      reason      = res$reason,
      stringsAsFactors = FALSE
    )
  }
})

qc <- do.call(rbind, qc_rows)

# ======= Write outputs =======
write.table(qc, file = qc_out, sep = "\t", quote = FALSE, row.names = FALSE)

bp <- qc[!is.na(qc$breakpoint), c("scaffold","breakpoint")]
write.table(bp, file = bp_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(qc, file = bp_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf(
  "Wrote %d scaffolds to %s (QC) and %s (breakpoints). %d satisfied the balance constraint (MAX_IMBALANCE=%.2f).",
  nrow(qc), qc_out, bp_out, sum(qc$satisfied, na.rm = TRUE), MAX_IMBALANCE
))
