# ============================
# Merge per-site FST files by position, clean, filter <0.15, annotate
# ============================
library(tidyverse)
library(readr)
library(readxl)
library(janitor)
library(data.table)
setwd("~/Desktop/OSU_projects/conifers/LP/FSTcalc/provenances/100_greedy/fst_per_site")

fst_dir <- "."

# ---- 1) Collect the per-site FST files ----
fst_files <- list.files(fst_dir, pattern = "^PROC.*_vs_.*\\.weir\\.fst$", full.names = TRUE)
stopifnot(length(fst_files) > 0)

# Helper: robust column picker
pick_col <- function(df, pattern) {
  cols <- names(df)
  hit <- cols[grepl(pattern, tolower(cols))]
  if (length(hit)) hit[1] else NA_character_
}

# Helper: read one FST file and return CHROM, POS, <pairname>_FST
read_fst_one <- function(path) {
  df <- suppressMessages(
    tryCatch(read_tsv(path, show_col_types = FALSE, progress = FALSE),
             error = function(e) read_delim(path, delim = " ", show_col_types = FALSE))
  )
  stopifnot(nrow(df) > 0)
  
  # infer columns
  chrom_col <- pick_col(df, "chrom|chr|scaffold|contig")
  pos_col   <- pick_col(df, "^pos$|position|site")
  # choose the non-weighted per-site FST column
  fst_col <- names(df)[grepl("fst", tolower(names(df))) & !grepl("weighted", tolower(names(df)))]
  fst_col <- fst_col[1]
  
  if (is.na(chrom_col) || is.na(pos_col) || is.na(fst_col)) {
    stop(paste("Could not find CHROM/POS/FST in", basename(path)))
  }
  
  # pair name from filename: "PROCx_vs_PROCy"
  base <- basename(path)
  pair <- sub("\\.weir\\.fst$", "", base)
  
  tibble(
    CHROM = df[[chrom_col]],
    POS   = as.integer(df[[pos_col]]),
    !!pair := suppressWarnings(as.numeric(df[[fst_col]]))
  )
}

message("Reading ", length(fst_files), " FST files …")
fst_list <- lapply(fst_files, read_fst_one)

# ---- 2) Merge all pairs by CHROM & POS (full join) ----
fst_wide <- reduce(fst_list, full_join, by = c("CHROM", "POS")) %>%
  arrange(CHROM, POS)

# Identify FST columns
fst_cols <- setdiff(names(fst_wide), c("CHROM","POS"))

# ---- 3) Clean values: NA/NaN/Inf -> 0; negatives -> 0 ----
clean_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  x
}
fst_clean <- fst_wide %>%
  mutate(across(all_of(fst_cols), clean_num))

# ---- 4) Rowwise max and filter < 0.15 ----
fst_clean <- fst_clean %>%
  mutate(row_max = do.call(pmax, c(across(all_of(fst_cols)), list(na.rm = TRUE))))

fst_lt005 <- fst_clean %>% filter(row_max > 0.15)

# Save intermediate/filtered tables
write_csv(fst_clean, file.path(fst_dir, "fst_per_site_merged_clean_v1.csv"))
write_csv(fst_lt005, file.path(fst_dir, "fst_per_site_merged_rowmax_lt0.15_v1.csv"))

message("Kept ", nrow(fst_lt005), " positions with row_max > 0.15 out of ", nrow(fst_clean))

