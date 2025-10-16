# Collapse overlapping intervals per contig (optionally also by orientation / NR.ID)
# Keep the longest; tie-breakers: higher coverage, then earlier start.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- CONFIG ---------------------------------------------------------------
setwd("~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly")

in_file <- "Aria_curated_annotation.tsv"   # for TSV/CSV (auto-delimiter detection)
out_file <- "annotations_dedup.tsv"



# Collapse only when these ALSO match (flip to TRUE if desired)
require_same_orientation <- FALSE
require_same_nr_id       <- FALSE

# ---- READ ----
dt <- fread(in_file)

# sanity
req <- c("contig","start_pos","end_pos","transcript_id","length")
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

# numeric + normalize coordinates
numify <- function(x) suppressWarnings(as.numeric(x))
dt[, `:=`(
  start_pos = numify(start_pos),
  end_pos   = numify(end_pos),
  length    = numify(length)
)]
if (!"coverage" %in% names(dt)) dt[, coverage := NA_real_] else dt[, coverage := numify(coverage)]

swap <- dt$start_pos > dt$end_pos
if (any(swap)) dt[swap, `:=`(start_pos = end_pos, end_pos = start_pos)]

# ---- GROUPING KEYS ----
cluster_keys <- c("contig")
if (require_same_orientation && "orientation" %in% names(dt)) cluster_keys <- c(cluster_keys, "orientation")
if (require_same_nr_id && "NR.ID" %in% names(dt))            cluster_keys <- c(cluster_keys, "NR.ID")

# ---- SORT (use setorderv, not tidyverse !!!) ----
ord_cols <- c(cluster_keys, "start_pos", "end_pos")
# make sure they all exist (defensive)
stopifnot(all(ord_cols %in% names(dt)))
setorderv(dt, ord_cols)   # ascending by default

# ---- ASSIGN OVERLAP CLUSTERS (length .N guaranteed) ----
dt[, .cluster_id := {
  if (.N == 0L) integer(0) else {
    out <- integer(.N)
    cid <- 1L
    cur_end <- end_pos[1L]
    out[1L] <- cid
    if (.N > 1L) {
      for (i in 2L:.N) {
        if (!is.na(start_pos[i]) && !is.na(cur_end) && start_pos[i] <= cur_end) {
          if (!is.na(end_pos[i]) && end_pos[i] > cur_end) cur_end <- end_pos[i]
        } else {
          cid <- cid + 1L
          cur_end <- end_pos[i]
        }
        out[i] <- cid
      }
    }
    out
  }
}, by = cluster_keys]

# ---- PICK WINNER PER CLUSTER ----
# Longest, then higher coverage, then earlier start
dt[, .sel_rank := frankv(list(-length, -coverage, start_pos), ties.method = "min"),
   by = c(cluster_keys, ".cluster_id")]

winners <- dt[.sel_rank == 1]

# Optional transparency about what got merged
merged_summary <- dt[, .(
  merged_transcripts = paste(unique(transcript_id), collapse = ";"),
  n_collapsed = .N
), by = c(cluster_keys, ".cluster_id")]

winners <- merge(winners, merged_summary,
                 by = c(cluster_keys, ".cluster_id"), all.x = TRUE)

# Cleanup & write
winners[, c(".cluster_id", ".sel_rank") := NULL]
orig <- names(dt); orig <- setdiff(orig, c(".cluster_id", ".sel_rank"))
setcolorder(winners, c(orig, setdiff(names(winners), orig)))
fwrite(winners, out_file, sep = "\t", na = "")

cat(sprintf("Collapsed %d -> %d rows. Grouped by: %s\n",
            nrow(dt), nrow(winners), paste(cluster_keys, collapse = ", ")))
