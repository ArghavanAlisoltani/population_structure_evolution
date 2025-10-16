# Collapse overlapping transcripts: keep the longest (per contig, with optional keys)
# Requires: data.table (fast) + dplyr (optional tidy niceties)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- CONFIG ---------------------------------------------------------------

in_file  <- "input_annotations.tsv"   # TSV/CSV (auto delim if fread)
out_file <- "annotations_dedup.tsv"

# Include additional grouping constraints for "overlap" definition:
# If TRUE, overlaps are collapsed only when these also match.
require_same_orientation <- FALSE   # set TRUE if you want + vs - kept separate
require_same_nr_id       <- FALSE   # set TRUE to merge only when NR.ID is identical

# Tie-breakers when lengths are identical:
prefer_higher_coverage <- TRUE

# ---- READ -----------------------------------------------------------------

dt <- data.table::fread(in_file)
# Ensure required columns
req <- c("contig","start_pos","end_pos","transcript_id","length","coverage")
missing <- setdiff(req, names(dt))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse=", "))

# Coerce numeric + fix reversed coords if any
dt[, start_pos := as.numeric(start_pos)]
dt[, end_pos   := as.numeric(end_pos)]
dt[, length    := as.numeric(length)]
if ("coverage" %in% names(dt)) dt[, coverage := as.numeric(coverage)] else dt[, coverage := NA_real_]

# Normalize: start <= end
swap_idx <- dt$start_pos > dt$end_pos
if (any(swap_idx)) {
  tmp <- dt$start_pos[swap_idx]
  dt$start_pos[swap_idx] <- dt$end_pos[swap_idx]
  dt$end_pos[swap_idx]   <- tmp
}

# ---- DEFINE CLUSTERING KEYS ----------------------------------------------

# Base key: contig
cluster_keys <- c("contig")
if (require_same_orientation && "orientation" %in% names(dt)) {
  cluster_keys <- c(cluster_keys, "orientation")
}
if (require_same_nr_id && "NR.ID" %in% names(dt)) {
  cluster_keys <- c(cluster_keys, "NR.ID")
}

# ---- BUILD OVERLAP CLUSTERS (linear sweep per key block) ------------------

setkeyv(dt, c(cluster_keys, "start_pos", "end_pos"))

# Function to assign cluster IDs to overlapping intervals within a block
assign_overlap_clusters <- function(block) {
  # block is a data.table sorted by start_pos, end_pos
  n <- nrow(block)
  if (n == 0L) return(integer(0))
  cluster_id <- integer(n)
  current_id <- 1L
  current_max_end <- block$end_pos[1]
  cluster_id[1] <- current_id

  for (i in 2:n) {
    s <- block$start_pos[i]
    e <- block$end_pos[i]
    # If starts before or at current_max_end => overlaps current cluster
    if (!is.na(s) && !is.na(current_max_end) && s <= current_max_end) {
      cluster_id[i] <- current_id
      # extend the cluster span if needed
      if (!is.na(e) && e > current_max_end) current_max_end <- e
    } else {
      # starts after current cluster ends => new cluster
      current_id <- current_id + 1L
      cluster_id[i] <- current_id
      current_max_end <- e
    }
  }
  cluster_id
}

# Apply per group of cluster_keys
dt[, .cluster_id := assign_overlap_clusters(.SD), by = cluster_keys]

# ---- SELECT WINNER PER CLUSTER -------------------------------------------

# Order by selection criteria:
# 1) longest length (desc)
# 2) higher coverage (desc) if available/desired
# 3) earlier start_pos (asc)
# 4) keep first (stable)
ord_expr <- list(-length, start_pos)  # base
if (prefer_higher_coverage && "coverage" %in% names(dt)) {
  ord_expr <- append(list(-coverage), ord_expr, after = 1)
}

# Rank rows within each cluster and take the top
dt[, .rank := frankv(list(!!!ord_expr), ties.method = "min"), by = c(cluster_keys, ".cluster_id")]
winners <- dt[.rank == 1]

# ---- OPTIONAL: keep a summary of what got merged (e.g., merged transcript_ids) ----
# Build a summary per cluster to attach back if you want transparency
cluster_summary <- dt[, .(
  merged_transcripts = paste(unique(transcript_id), collapse = ";"),
  n_collapsed = .N
), by = c(cluster_keys, ".cluster_id")]

winners <- merge(
  winners,
  cluster_summary,
  by = c(cluster_keys, ".cluster_id"),
  all.x = TRUE
)

# ---- CLEAN & WRITE --------------------------------------------------------

# You can drop helper columns if you like:
winners[, c(".cluster_id", ".rank") := NULL]

# Keep original column order where possible and append summaries at the end
orig_cols <- names(dt)
orig_cols <- setdiff(orig_cols, c(".cluster_id", ".rank"))  # remove helpers
keep_order <- c(orig_cols, setdiff(names(winners), orig_cols))

setcolorder(winners, keep_order)

fwrite(winners, out_file, sep = "\t", na = "")

cat(sprintf(
  "Collapsed %d -> %d rows (%d removed). Grouped by: %s. Overlap criterion: intervals on same group.\n",
  nrow(dt), nrow(winners), nrow(dt) - nrow(winners), paste(cluster_keys, collapse = ", ")
))
