# dotplots_scaffolds.R visualize and save plots for each scaffold. Please find the example of runs at the end of the code.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# Helper to make safe filenames
sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# Main function
make_scaffold_dotplots <- function(
    infile,
    out_dir = "plots_out",
    bin_size = NULL,          # e.g., 1e5 for 100 kb. Use NULL for unbinned.
    units = c("Mb","kb","bp"),# controls x-axis units/labels
    multi_pdf = TRUE,         # TRUE = one multi-page PDF; FALSE = one image per scaffold
    device = c("pdf","png"),  # used when multi_pdf=FALSE
    width = 10, height = 3
) {
  units <- match.arg(units)
  device <- match.arg(device)
  
  # unit scale for x-axis
  scale_factor <- switch(units, Mb = 1e-6, kb = 1e-3, bp = 1)
  x_lab <- paste0("Position (", units, ")")
  
  # Read data (expects CHROM and POS)
  df <- readr::read_tsv(infile, col_types = cols())
  if (!all(c("CHROM","POS") %in% names(df))) {
    stop("Input must contain columns named 'CHROM' and 'POS'.")
  }
  df <- df %>% mutate(POS = as.numeric(POS))
  
  # Output dir
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  scaffolds <- unique(df$CHROM)
  
  # If unbinned: points along y=0 (with slight vertical jitter so dots are visible)
  if (is.null(bin_size)) {
    message("Plotting UNBINNED dot plots by scaffold...")
    
    # Multi-page PDF option
    if (multi_pdf) {
      pdf(file = file.path(out_dir, "scaffold_dots_unbinned.pdf"),
          width = width, height = height, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
    }
    
    for (sc in scaffolds) {
      d <- df %>% filter(CHROM == sc) %>%
        mutate(xu = POS * scale_factor)
      
      p <- ggplot(d, aes(x = xu)) +
        geom_point(aes(y = 0), size = 0.6, alpha = 0.6) +
        # Add a subtle rug to emphasize dense regions
        geom_rug(sides = "b", alpha = 0.25) +
        scale_y_continuous(NULL, breaks = NULL) +
        labs(title = sc, x = x_lab, y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_blank()
        )
      
      if (multi_pdf) {
        print(p)
      } else {
        out_file <- file.path(out_dir, paste0("dots_unbinned_", sanitize_filename(sc), ".", device))
        ggsave(out_file, p, width = width, height = height, dpi = 300)
      }
    }
    
  } else {
    # BINNED: non-overlapping windows of size bin_size (bp)
    message("Plotting BINNED (non-overlapping) dot plots by scaffold...")
    if (!is.numeric(bin_size) || bin_size <= 0) {
      stop("bin_size must be a positive number of base pairs (e.g., 1e5 for 100 kb).")
    }
    
    # Precompute binned counts per scaffold
    df_binned <- df %>%
      mutate(bin = floor((POS - 1) / bin_size)) %>%
      group_by(CHROM, bin) %>%
      summarise(
        n_snp = n(),
        start_bp = bin * bin_size + 1,
        end_bp   = (bin + 1) * bin_size,
        mid_bp   = (start_bp + end_bp) / 2,
        .groups = "drop"
      ) %>%
      mutate(mid_u = mid_bp * scale_factor)
    
    if (nrow(df_binned) == 0) {
      stop("After binning, no data remained—check bin_size and input.")
    }
    
    # Multi-page PDF option
    if (multi_pdf) {
      pdf(file = file.path(out_dir, paste0("scaffold_dots_binned_", format(bin_size, scientific = FALSE), "bp.pdf")),
          width = width, height = height, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
    }
    
    for (sc in scaffolds) {
      d <- df_binned %>% filter(CHROM == sc)
      
      p <- ggplot(d, aes(x = mid_u, y = n_snp)) +
        # Lollipop style: a stem up to the dot
        geom_segment(aes(xend = mid_u, y = 0, yend = n_snp), linewidth = 0.3, alpha = 0.6) +
        geom_point(size = 1.2, alpha = 0.9) +
        labs(
          title = paste0(sc, "  (bin = ", format(bin_size, scientific = FALSE), " bp)"),
          x = x_lab, y = "SNPs per bin"
        ) +
        theme_minimal(base_size = 12)
      
      if (multi_pdf) {
        print(p)
      } else {
        out_file <- file.path(
          out_dir,
          paste0("dots_binned_", format(bin_size, scientific = FALSE), "bp_", sanitize_filename(sc), ".", device)
        )
        ggsave(out_file, p, width = width, height = height, dpi = 300)
      }
    }
  }
  
  message("Done. Files written to: ", normalizePath(out_dir))
}

# ============================
# Example usage
# ============================
setwd("~/Desktop/OSU_projects/conifers/LP/vcf_v1")
#1) Unbinned, one multi-page PDF (one page per scaffold)
make_scaffold_dotplots("positions_split_poly_s100_scaffolds.tsv",
                       out_dir = "plots_unbinned", bin_size = NULL,
                       units = "Mb", multi_pdf = TRUE)

#2) Binned at 100 kb, one multi-page PDF (one page per scaffold)
make_scaffold_dotplots("positions_split_poly_s100_scaffolds.tsv",
                       out_dir = "plots_binned_100kb", bin_size = 1e5,
                       units = "Mb", multi_pdf = TRUE)

# 3) If you prefer separate PNGs per scaffold instead of a single PDF:
make_scaffold_dotplots("positions_split_poly_s100_scaffolds.tsv",
                       out_dir = "plots_unbinned_png", bin_size = NULL,
                       units = "Mb", multi_pdf = FALSE, device = "png")

make_scaffold_dotplots("positions_split_poly_s100_scaffolds.tsv",
                       out_dir = "plots_binned_500kb_png", bin_size = 5e5,
                       units = "Mb", multi_pdf = FALSE, device = "png")

