# De-duplicate by (contig, start_pos, end_pos, orientation, gene_id)
# collapsing rows that only differ by tissue_source.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(readxl)
  library(tidyr)
  library(purrr)
})

# ---- CONFIG ---------------------------------------------------------------

# Set one of these and leave the other NULL
in_file_tsv <- "input_annotations.tsv"   # for TSV/CSV (auto-delimiter detection)
in_file_xlsx <- NULL                     # e.g., "Example_Gene_V1_header_gpt_input.xlsx"

out_file <- "annotations_dedup.tsv"

# If TRUE, try to normalize transcript IDs down to a gene-level ID
# by removing common isoform suffixes (e.g., .t1, -RA, _i2).
strip_isoform_suffix <- TRUE

# Prefer rows with higher coverage; then longer length; then first occurrence
prefer_by <- c("coverage", "length")

# ---- READ ----------------------------------------------------------------

read_table_flex <- function(in_file_tsv, in_file_xlsx) {
  if (!is.null(in_file_xlsx)) {
    # First (or only) sheet by default
    df <- readxl::read_excel(in_file_xlsx)
    df <- as.data.frame(df)
  } else if (!is.null(in_file_tsv)) {
    # Try to guess delimiter (tab vs comma)
    # readr::read_delim will auto-detect when delim = NULL (readr >= 2.0)
    df <- readr::read_delim(in_file_tsv, delim = NULL, guess_max = 100000, show_col_types = FALSE)
  } else {
    stop("Provide either in_file_tsv or in_file_xlsx.")
  }
  df
}

df <- read_table_flex(in_file_tsv, in_file_xlsx)

# ---- VALIDATE REQUIRED COLUMNS -------------------------------------------

required_cols <- c(
  "contig","tissue_source","start_pos","end_pos","orientation","transcript_id",
  "length","coverage","Unique_ID"
)

missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
}

# ---- HELPERS --------------------------------------------------------------

# Try to produce a gene-level ID from transcript_id
derive_gene_id <- function(x, strip = TRUE) {
  x <- as.character(x)
  if (!strip) return(x)
  # remove common isoform suffixes (edit/add patterns as needed)
  # Examples matched: ".t1", ".T2", "-RA", "-RB", "_i1", "_iso2"
  x %>%
    str_replace("[._-][tT]\\d+$", "") %>%      # .t1 / -t2 / _T3
    str_replace("[-_][Rr][A-Za-z0-9]+$", "") %>% # -RA / -RB / _Ra
    str_replace("[-_][iI][0-9]+$", "") %>%     # _i1 / -i2
    str_replace("[._-]iso\\d+$", "")           # .iso1 / -iso2
}

# Ensure numeric types for ranking
to_num <- function(v) suppressWarnings(as.numeric(v))

# ---- MAIN ----------------------------------------------------------------

df1 <- df %>%
  mutate(
    start_pos = to_num(start_pos),
    end_pos   = to_num(end_pos),
    coverage  = to_num(coverage),
    length    = to_num(length),
    gene_id   = derive_gene_id(transcript_id, strip = strip_isoform_suffix)
  )

# Redundancy key: same region + same gene
group_key <- c("contig","start_pos","end_pos","orientation","gene_id")

# Aggregate tissues per redundant group (preserve info)
tissues_by_group <- df1 %>%
  group_by(across(all_of(group_key))) %>%
  summarize(
    tissues_merged = paste(sort(unique(tissue_source)), collapse = ";"),
    .groups = "drop"
  )

# Choose a single representative per redundant group
order_cols <- intersect(prefer_by, names(df1))
if (length(order_cols) == 0) order_cols <- character(0)

dedup <- df1 %>%
  group_by(across(all_of(group_key))) %>%
  {
    # build arrange() dynamically: descending by each preference
    if (length(order_cols) > 0) {
      arrange(., across(all_of(order_cols), ~ desc(.x)), .by_group = TRUE)
    } else {
      .
    }
  } %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  # attach merged tissue list (useful when we drop duplicates)
  left_join(tissues_by_group, by = group_key) %>%
  # keep original columns first, then add helpers at the end
  select(all_of(names(df)), gene_id, tissues_merged)

# ---- WRITE ----------------------------------------------------------------

readr::write_tsv(dedup, out_file, na = "")

message(sprintf(
  "Done. Input rows: %d | Output rows: %d | Collapsed by: %s",
  nrow(df), nrow(dedup), paste(group_key, collapse = ", ")
))
