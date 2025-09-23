library(data.table)
library(readxl)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
setwd("~/Desktop/OSU_projects/conifers/LP/FSTcalc/provenances/100_greedy/fst_per_site")

# === Inputs ===
# Vector of vcftools per-site FST files (one per comparison)
files <- c(
  "PROC3_vs_PROC4.weir.fst",
  "PROC3_vs_PROC5.weir.fst",
  "PROC4_vs_PROC5.weir.fst"
)

# --- Libraries ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --- Helper to read and bin one file ---
read_and_bin <- function(path) {
  dt <- fread(path)
  # vcftools column is typically 'WEIR_AND_COCKERHAM_FST'
  fst_col <- grep("WEIR.*FST", names(dt), value = TRUE)
  if (length(fst_col) != 1) stop("Could not find FST column in: ", path)
  dt[, FST := as.numeric(get(fst_col))]
  # treat negatives as ~0
  dt[!is.na(FST) & FST < 0, FST := 0]
  # define bins
  breaks <- c(0, 0.05, 0.1, 0.15, 0.25, 1)
  labels <- c("<0.05", "0.05–0.1", "0.1–0.15", "0.15-0.25", "0.25-1")
  dt[, bin := cut(FST, breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)]
  # summarize
  out <- dt[, .N, by = bin][CJ(bin = labels, unique = TRUE), on = "bin"][is.na(N), N := 0L]
  out[, total := sum(N)]
  out[, pct := 100 * N / total]
  out[, comparison := sub("\\.weir.*$", "", basename(path))]
  out
}

# --- Combine all comparisons ---
summ <- rbindlist(lapply(files, read_and_bin))
#fwrite(summ, "fst_bin_summary_R.csv")

# --- Plot: grouped bars (percent per bin by comparison) ---
gg <- ggplot(summ, aes(x = comparison, y = pct, fill = bin)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(x = "Comparison", y = "% SNPs", fill = "FST bin",
       title = "Per-SNP FST distribution (binned)") +
  coord_flip() +
  scale_fill_manual(values = c("lightblue","pink","khaki","lightgreen","plum"))+
  theme_bw(base_size = 12)
gg
ggsave("fst_bin_summary_procs345_R.png", gg, width = 8, height = 5, dpi = 150)


# ============================
# Pairwise FST -> lower-tri table
# Updated for files like:
#   /mnt/data/PROC1_vs_PROC2.weir.fst
# ============================
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)

# Point to your folder with the attached files
fst_dir <- "."

# Pick up the attached pairwise files
files <- list.files(
  fst_dir,
  pattern = "^PROC[0-9]+_vs_PROC[0-9]+\\.weir\\.fst$",
  full.names = TRUE
)

stopifnot(length(files) > 0)

# Parse provenance codes from filenames like "PROC1_vs_PROC2.weir.fst"
parse_pair <- function(fname) {
  m <- str_match(basename(fname), "^(PROC[0-9]+)_vs_(PROC[0-9]+)\\.weir\\.fst$")
  if (is.na(m[1,1])) return(NULL)
  tibble(pop1 = m[1,2], pop2 = m[1,3])
}

# Compute a single FST per pair from one file (per-site OR windowed)
compute_one <- function(path) {
  nm <- parse_pair(path)
  if (is.null(nm)) return(NULL)
  
  # Try reading flexibly (tsv/space-delimited)
  df <- tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) tryCatch(read_delim(path, delim = " ", show_col_types = FALSE),
                                 error = function(e2) NULL)
  )
  if (is.null(df)) return(NULL)
  
  cols <- tolower(names(df))
  
  # Windowed?
  is_windowed <- any(str_detect(cols, "weighted_fst")) &&
    any(str_detect(cols, "n_variants|n_snps|n_sites"))
  
  if (is_windowed) {
    # Use WEIGHTED_FST and weight by the number of variants in the window
    wcol <- names(df)[which(str_detect(cols, "weighted_fst"))[1]]
    ncol <- names(df)[which(str_detect(cols, "n_variants|n_snps|n_sites"))[1]]
    tmp <- df %>%
      transmute(
        wfst = !! rlang::sym(wcol),
        n    = suppressWarnings(as.numeric(!! rlang::sym(ncol)))
      ) %>%
      filter(is.finite(wfst), is.finite(n), n > 0)
    # Truncate negative windows to 0 (common, stabilizing)
    tmp$wfst <- pmax(tmp$wfst, 0)
    val <- weighted.mean(tmp$wfst, w = tmp$n)
  } else {
    # Per-site file: find the FST column (usually WEIR_AND_COCKERHAM_FST)
    fstcol_idx <- which(str_detect(cols, "fst"))
    stopifnot(length(fstcol_idx) >= 1)
    fstcol <- names(df)[fstcol_idx[1]]
    
    tmp <- df %>%
      transmute(fst = suppressWarnings(as.numeric(!! rlang::sym(fstcol)))) %>%
      filter(is.finite(fst))
    # Truncate negative per-site estimates to 0 for the summary
    tmp$fst <- pmax(tmp$fst, 0)
    val <- mean(tmp$fst)
  }
  
  tibble(pop1 = nm$pop1, pop2 = nm$pop2, fst = as.numeric(val))
}

pairwise <- map_dfr(files, compute_one)

# Choose display order (PROC1..PROC5). Adjust if you have more/other codes.
prov_order <- sort(unique(c(pairwise$pop1, pairwise$pop2)))
# Example if you want a fixed order:
# prov_order <- paste0("PROC", 1:5)

# Helper: build a lower-triangular table like your example
make_lower_tri <- function(pw, prov_order, digits = 3, blank_upper = TRUE, diag_blank = TRUE) {
  pw2 <- pw %>%
    mutate(pair = paste(pmin(pop1, pop2), pmax(pop1, pop2), sep = " vs ")) %>%
    group_by(pair) %>% slice(1) %>% ungroup() %>%
    transmute(pop1 = pmin(pop1, pop2),
              pop2 = pmax(pop1, pop2),
              fst  = fst)
  
  M <- matrix(NA_real_, nrow = length(prov_order), ncol = length(prov_order),
              dimnames = list(prov_order, prov_order))
  
  for (i in seq_len(nrow(pw2))) {
    a <- pw2$pop1[i]; b <- pw2$pop2[i]; v <- pw2$fst[i]
    if (a %in% prov_order && b %in% prov_order) {
      M[a, b] <- v
      M[b, a] <- v
    }
  }
  
  if (diag_blank) diag(M) <- NA_real_
  if (blank_upper) M[upper.tri(M, diag = FALSE)] <- NA_real_
  
  M_round <- format(round(M, digits), nsmall = digits)
  M_round[is.na(M)] <- ""
  M_round
}

lower_table <- make_lower_tri(pairwise, prov_order, digits = 3)

# ----- Optional: replace PROC codes with human-friendly names -----
# Fill this mapping if you want pretty row/col names:
# label_map <- c(
#   PROC1 = "Deer Mtn.",
#   PROC2 = "Inverness",
#   PROC3 = "Judy Cr.",
#   PROC4 = "Swan Hills",
#   PROC5 = "Virginia H."
# )
# if (exists("label_map")) {
#   # rename the dimnames of the lower-tri table
#   rn <- rownames(lower_table); cn <- colnames(lower_table)
#   rownames(lower_table) <- ifelse(rn %in% names(label_map), label_map[rn], rn)
#   colnames(lower_table) <- ifelse(cn %in% names(label_map), label_map[cn], cn)
# }

# Show the table in the console
lower_table

# Save outputs
out_pairs <- file.path(fst_dir, "pairwise_fst_values.csv")
out_lower <- file.path(fst_dir, "pairwise_fst_lower_triangle.csv")
readr::write_csv(pairwise, out_pairs)
# For the lower-tri, write as a simple CSV matrix
write.table(lower_table, file = out_lower, sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

message("Saved:\n  ", out_pairs, "\n  ", out_lower)

# ============================
# Pairwise global FST + block-jackknife CIs
# ============================
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)

fst_dir <- "."   # folder with PROCx_vs_PROCy files

# --- Helper: robust read (tab or space) ---
read_any <- function(path) {
  tryCatch(
    suppressMessages(read_tsv(path, show_col_types = FALSE, progress = FALSE)),
    error = function(e) suppressMessages(read_delim(path, delim = " ", show_col_types = FALSE))
  )
}

# --- Helper: detect common columns case-insensitively ---
pick_col <- function(df, pattern) {
  cols <- names(df)
  hit <- cols[grepl(pattern, tolower(cols))]
  if (length(hit)) hit[1] else NA_character_
}

# --- Parse pair names from "PROC1_vs_PROC2..." ---
parse_pair <- function(fname) {
  m <- str_match(basename(fname), "^(.+?)_vs_(.+?)\\.")
  if (is.na(m[1,1])) return(NULL)
  tibble(pop1 = m[1,2], pop2 = m[1,3])
}

# --- Core: compute point estimate + block-jackknife CI for 1 file ---
compute_jackknife_one <- function(path,
                                  block_by = c("chrom", "snps"),
                                  snps_per_block = 1000L) {
  block_by <- match.arg(block_by)
  nm <- parse_pair(path)
  if (is.null(nm)) return(NULL)
  
  df <- read_any(path)
  if (is.null(df) || !nrow(df)) return(NULL)
  
  # Column detection
  chrom_col <- pick_col(df, "chrom|chr|scaffold|contig")
  pos_col   <- pick_col(df, "pos|position|bin_start|start")
  
  wfst_col  <- pick_col(df, "weighted_fst")
  n_col     <- pick_col(df, "n_variants|n_snps|n_sites")
  fst_col   <- {
    cand <- names(df)[grepl("fst", tolower(names(df)))]
    cand <- cand[!grepl("weighted", tolower(cand))]
    if (length(cand)) cand[1] else NA_character_
  }
  
  # Identify "mode"
  is_windowed <- !is.na(wfst_col) && !is.na(n_col)
  
  if (is_windowed) {
    # ---- Windowed mode (preferred) ----
    if (is.na(chrom_col)) stop("Need a CHROM/contig column for block jackknife.")
    dat <- df %>%
      transmute(
        CHROM = .data[[chrom_col]],
        value = as.numeric(.data[[wfst_col]]),
        n     = as.numeric(.data[[n_col]])
      ) %>%
      filter(is.finite(value), is.finite(n), n > 0)
    
    # truncate negative window estimates to 0 (stabilizes summaries)
    dat$value <- pmax(dat$value, 0)
    
    # blocks
    if (block_by == "chrom") {
      dat$block <- as.character(dat$CHROM)
    } else {
      # fixed-SNP blocks: order by CHROM then assume rows are roughly ordered
      dat <- dat %>% arrange(CHROM)
      idx <- rep(seq_len(ceiling(nrow(dat) / snps_per_block)),
                 each = snps_per_block, length.out = nrow(dat))
      dat$block <- paste0("B", idx)
    }
    
    # full estimate: weighted mean by the #variants in window
    theta_full <- with(dat, weighted.mean(value, w = n))
    
    # jackknife over blocks
    blocks <- unique(dat$block)
    thetas <- vapply(blocks, function(b) {
      tmp <- dat %>% filter(block != b)
      with(tmp, weighted.mean(value, w = n))
    }, numeric(1))
    
  } else {
    # ---- Per-site mode (approximation) ----
    if (is.na(fst_col)) stop("Could not find an FST column.")
    if (is.na(chrom_col)) stop("Need a CHROM/contig column for block jackknife.")
    
    dat <- df %>%
      transmute(
        CHROM = .data[[chrom_col]],
        value = suppressWarnings(as.numeric(.data[[fst_col]]))
      ) %>%
      filter(is.finite(value))
    
    # common practice: set negative per-site estimates to 0 before summarising
    dat$value <- pmax(dat$value, 0)
    
    # blocks
    if (block_by == "chrom") {
      dat$block <- as.character(dat$CHROM)
    } else {
      dat <- dat %>% arrange(CHROM)
      idx <- rep(seq_len(ceiling(nrow(dat) / snps_per_block)),
                 each = snps_per_block, length.out = nrow(dat))
      dat$block <- paste0("B", idx)
    }
    
    # full estimate: mean of (truncated) per-site FST (approximation)
    theta_full <- mean(dat$value)
    
    # jackknife over blocks: recompute the mean after dropping each block
    blocks <- unique(dat$block)
    thetas <- vapply(blocks, function(b) {
      mean(dat$value[dat$block != b])
    }, numeric(1))
  }
  
  B <- length(blocks)
  if (B < 5) warning("Few blocks (", B, "); CIs may be unstable. Consider larger/fewer blocks.")
  
  theta_bar <- mean(thetas)
  # Jackknife variance estimator
  var_jack <- (B - 1) / B * sum((thetas - theta_bar)^2)
  se_jack  <- sqrt(var_jack)
  ci95     <- c(theta_full - 1.96 * se_jack, theta_full + 1.96 * se_jack)
  
  tibble(
    pop1 = nm$pop1,
    pop2 = nm$pop2,
    mode = if (is_windowed) "window_weighted" else "per_site_mean",
    n_blocks = B,
    est = theta_full,
    se = se_jack,
    ci_low = max(0, ci95[1]),
    ci_high = min(1, ci95[2])
  )
}

# ---- Run over all pairwise files in folder ----
files <- list.files(fst_dir, pattern = "_vs_.*\\.weir(\\.windowed)?\\.fst$", full.names = TRUE)
stopifnot(length(files) > 0)

res <- map_dfr(files, compute_jackknife_one, block_by = "chrom", snps_per_block = 1000)

# Pretty lower-tri table from point estimates
prov_order <- sort(unique(c(res$pop1, res$pop2)))
make_lower_tri <- function(pw, order, digits = 3) {
  M <- matrix(NA_real_, length(order), length(order), dimnames = list(order, order))
  for (i in seq_len(nrow(pw))) {
    a <- pw$pop1[i]; b <- pw$pop2[i]; v <- pw$est[i]
    if (a %in% order && b %in% order) { M[a, b] <- v; M[b, a] <- v }
  }
  diag(M) <- NA; M[upper.tri(M)] <- NA
  out <- format(round(M, digits), nsmall = digits); out[is.na(M)] <- ""
  out
}

lower_table <- make_lower_tri(res, prov_order, digits = 3)

# Save results
write_csv(res, file.path(fst_dir, "pairwise_fst_jackknife_results.csv"))
write.table(lower_table, file = file.path(fst_dir, "pairwise_fst_lower_triangle.csv"),
            sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

print(res)
lower_table

##########################################################
# visualize
##########################################################
# Ensure a consistent order for axes
procs <- sort(unique(c(res$pop1, res$pop2)))

# ---- 2) Heatmap (symmetric matrix with lower-tri shown) ----
# Build full symmetric matrix of estimates
M <- matrix(NA_real_, length(procs), length(procs),
            dimnames = list(procs, procs))
for (i in seq_len(nrow(res))) {
  a <- res$pop1[i]; b <- res$pop2[i]; v <- res$est[i]
  M[a,b] <- v; M[b,a] <- v
}
diag(M) <- NA
# Mask upper triangle to emulate lower-tri view
M[upper.tri(M)] <- NA

heat_df <- as.data.frame(as.table(M)) 
names(heat_df)<-c("pop1", "pop2", "est") 
heat_df<-heat_df[!is.na(heat_df$est),]

p_heat <- ggplot(heat_df, aes(pop2, pop1, fill = est)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", est)), color = "white", size = 8) +
  scale_fill_gradient(low = "#2c7bb6", high = "#d7191c", na.value = "grey90") +
  labs(title = "Pairwise FST (lower-tri view)",
       x = "", y = "", fill = "FST") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_heat)

# ---- 3) CI “forest” plot (sorted by estimate) ----
res_plot <- res %>%
  mutate(pair = paste(pmin(pop1,pop2), pmax(pop1,pop2), sep = " vs ")) %>%
  arrange(est) %>%
  mutate(pair = factor(pair, levels = pair))

p_forest <- ggplot(res_plot, aes(x = est, y = pair)) +
  geom_point() +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15) +
  geom_vline(xintercept = c(0.05, 0.15, 0.25), linetype = "dashed", alpha = 0.4) +
  labs(title = "Global pairwise FST with 95% block-jackknife CIs",
       x = "FST", y = "Pair") +
  theme_minimal(base_size = 13)
print(p_forest)

# ---- 4) MDS map from the FST matrix (treat FST as a distance) ----
# Fill missing with 0 on diagonal, ensure symmetry
M_full <- matrix(0, length(procs), length(procs),
                 dimnames = list(procs, procs))
for (i in seq_len(nrow(res))) {
  a <- res$pop1[i]; b <- res$pop2[i]; v <- res$est[i]
  M_full[a,b] <- v; M_full[b,a] <- v
}
# Classical MDS (2D)
mds <- cmdscale(as.dist(M_full), k = 2, eig = TRUE)
mds_df <- tibble(proc = rownames(mds$points),
                 Dim1 = mds$points[,1],
                 Dim2 = mds$points[,2])

p_mds <- ggplot(mds_df, aes(Dim1, Dim2, label = proc)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel() +
  labs(title = "MDS of pairwise FST distances",
       x = "Dim 1", y = "Dim 2") +
  theme_minimal(base_size = 12)
print(p_mds)
