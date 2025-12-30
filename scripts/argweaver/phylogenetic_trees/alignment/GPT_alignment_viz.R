############################################################
## ARGweaver .sites -> selected positions alignment plot
## + FASTA output with sample names
##
## Run in RStudio:
## source("sites_selected_positions_alignment.R")
############################################################
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/Dec_16_2025")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# User inputs
# -----------------------------
sites_file <- "sites_for_tree_0_rep/outargs_scaffold_4_900000001_1050000000.0.sites"  # <-- change
outdir <- "sites_selected_positions_out"                          # <-- change if needed
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# positions to extract (order matters!)
selected_positions <- c(
  983057685,
  983057688,
  983057694,
  983057707,
  983057714,
  983057718,
  983057724
)

# output prefix for files
prefix <- sprintf("selected_sites_%s_pos%d_%d",
                  gsub("\\.sites$", "", basename(sites_file)),
                  min(selected_positions), max(selected_positions))

# -----------------------------
# Helper: read .sites
# -----------------------------
read_argweaver_sites <- function(path) {
  # Read all lines fast
  lines <- readLines(path, warn = FALSE)
  
  # Find NAMES line
  nm_i <- grep("^#NAMES\\b", lines)
  if (length(nm_i) != 1) stop("Could not find a single #NAMES line in: ", path)
  names_line <- lines[nm_i]
  samples <- strsplit(sub("^#NAMES\\s+", "", names_line), "\\s+")[[1]]
  
  # Find REGION line (optional but useful)
  reg_i <- grep("^#REGION\\b", lines)
  region <- NA_character_
  if (length(reg_i) >= 1) region <- lines[reg_i[1]]
  
  # Data lines = those that start with a number (POS)
  dat_i <- grep("^[0-9]+\\s+", lines)
  if (length(dat_i) == 0) stop("No site rows found (POS SEQ...) in: ", path)
  
  # Parse into data.table: pos + seqstring
  # Some seq strings can be long; keep as character
  dt <- rbindlist(lapply(lines[dat_i], function(x) {
    # split into first token (pos) and the rest (sequence)
    # sequence has no spaces in normal .sites; still safe:
    parts <- strsplit(x, "\\s+")[[1]]
    list(pos = as.integer(parts[1]), seq = parts[2])
  }))
  
  # Sanity check: sequence length == number of samples
  n <- length(samples)
  bad <- dt[nchar(seq) != n]
  if (nrow(bad) > 0) {
    stop(sprintf(
      "Found %d site rows where seq length != #samples (%d). Example pos=%s length=%d",
      nrow(bad), n, bad$pos[1], nchar(bad$seq[1])
    ))
  }
  
  list(samples = samples, region = region, dt = dt)
}

# -----------------------------
# Helper: build sample x position matrix
# -----------------------------
extract_positions_matrix <- function(sites_obj, positions, fill_missing = "N") {
  samples <- sites_obj$samples
  dt <- sites_obj$dt
  
  # Keep only requested positions that exist
  dt_sub <- dt[pos %in% positions]
  found_pos <- sort(unique(dt_sub$pos))
  missing_pos <- setdiff(positions, found_pos)
  
  # Create a matrix [samples x positions] in the *requested order*
  mat <- matrix(fill_missing, nrow = length(samples), ncol = length(positions))
  rownames(mat) <- samples
  colnames(mat) <- as.character(positions)
  
  # Fill from dt_sub
  # For each position, split seqstring into bases
  for (p in found_pos) {
    s <- dt_sub[pos == p][1, seq]
    bases <- strsplit(s, "")[[1]]
    mat[, as.character(p)] <- bases
  }
  
  list(mat = mat, found = found_pos, missing = missing_pos)
}

# -----------------------------
# Helper: write FASTA
# -----------------------------
write_fasta_from_matrix <- function(mat, fasta_path) {
  con <- file(fasta_path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_len(nrow(mat))) {
    nm <- rownames(mat)[i]
    seq <- paste0(mat[i, ], collapse = "")
    writeLines(paste0(">", nm), con)
    writeLines(seq, con)
  }
}

# -----------------------------
# Helper: plot alignment tiles
# -----------------------------
plot_alignment_tiles <- function(mat, title = NULL, tile_height = 0.9, base_font = 11) {
  # Long format
  df <- as.data.table(mat, keep.rownames = "sample")
  df_long <- melt(df, id.vars = "sample", variable.name = "position", value.name = "base")
  
  # Make sample order match the tree/ARG order (as in #NAMES line)
  df_long[, sample := factor(sample, levels = rev(rownames(mat)))]  # reverse so first sample is top
  
  # position as factor to preserve order
  df_long[, position := factor(position, levels = colnames(mat))]
  
  ggplot(df_long, aes(x = position, y = sample, fill = base)) +
    geom_tile(height = tile_height) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    labs(title = title, x = "Position", y = NULL, fill = "Base") +
    theme_bw(base_size = base_font) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

# -----------------------------
# MAIN
# -----------------------------
sites_obj <- read_argweaver_sites(sites_file)

res <- extract_positions_matrix(sites_obj, selected_positions, fill_missing = "N")
mat <- res$mat

# Report missing positions
if (length(res$missing) > 0) {
  message("WARNING: These requested positions were not found in the .sites file and were filled with 'N':")
  message(paste(res$missing, collapse = ", "))
} else {
  message("All requested positions were found in the .sites file.")
}

# Write TSV matrix (samples x positions)
tsv_path <- file.path(outdir, paste0(prefix, ".matrix.tsv"))
mat_dt <- data.table(sample = rownames(mat), mat)
fwrite(mat_dt, tsv_path, sep = "\t")

# Write FASTA
fasta_path <- file.path(outdir, paste0(prefix, ".fasta"))
write_fasta_from_matrix(mat, fasta_path)

# Plot alignment
plot_title <- paste0(
  basename(sites_file),
  "\nSelected positions: ", paste(selected_positions, collapse = ", ")
)
p <- plot_alignment_tiles(mat, title = plot_title, tile_height = 0.9, base_font = 11)
p
# Save figure
png_path <- file.path(outdir, paste0(prefix, ".alignment.png"))
pdf_path <- file.path(outdir, paste0(prefix, ".alignment.pdf"))

ggsave(png_path, p, width = 10, height = max(6, 0.03 * nrow(mat) + 2), dpi = 300)
ggsave(pdf_path, p, width = 10, height = max(6, 0.03 * nrow(mat) + 2))

message("Wrote:")
message("  ", tsv_path)
message("  ", fasta_path)
message("  ", png_path)
message("  ", pdf_path)
