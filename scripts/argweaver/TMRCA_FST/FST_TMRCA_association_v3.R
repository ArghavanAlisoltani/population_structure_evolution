############################################################
##  FST (per SNP) + TMRCA segments merge + regression plots
##  Run in RStudio (no terminal needed)
############################################################
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/MAF_FST_TMRCA")

## ---------------------------
## 0) Packages (install if needed)
## ---------------------------
pkgs <- c("data.table", "ggplot2")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(data.table)
library(ggplot2)

## ---------------------------
## 1) User inputs (EDIT THESE PATHS)
## ---------------------------
fst_file   <- "FST_per_site_v1.tsv"
tmrca_file <- "tmrca_v1.tsv"

out_merged_tsv <- "FST_plus_meanTMRCA_v3_a.tsv"
out_plot_dir   <- "fst_tmrca_plots_spearman_v3_filtered_tmrca200"
dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)

## NEW: filtered outputs directory (keeps original merged TSV unchanged)
out_filtered_dir <- file.path(out_plot_dir, "filtered_tables")
dir.create(out_filtered_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------
## 2) Reading files
## ---------------------------
fst  <- fread(fst_file)
tmrca <- fread(tmrca_file)

## ---------------------------
## 3) Basic column checks + standardization
## ---------------------------
stopifnot(all(c("Scaffold", "POS") %in% names(fst)))
stopifnot(all(c("CHROM", "start_tmrca", "end_tmrca", "mean_tmrca") %in% names(tmrca)))

fst[, POS := as.numeric(POS)]
tmrca[, start_tmrca := as.numeric(start_tmrca)]
tmrca[, end_tmrca   := as.numeric(end_tmrca)]
tmrca[, mean_tmrca  := as.numeric(mean_tmrca)]

fst_cols <- intersect(c("PROC3_vs_PROC4", "PROC3_vs_PROC5", "PROC4_vs_PROC5"), names(fst))
for (cc in fst_cols) fst[, (cc) := as.numeric(get(cc))]

## ---------------------------
## 4) Interval join: map each SNP position to the segment it falls in
##    Uses [start, end) via tiny epsilon subtraction from end.
## ---------------------------
eps <- 1e-6

tmrca2 <- copy(tmrca)
tmrca2[, `:=`(
  seg_start = start_tmrca,
  seg_end   = end_tmrca - eps
)]

fst2 <- copy(fst)
fst2[, `:=`(
  snp_start = POS,
  snp_end   = POS
)]

setkey(tmrca2, CHROM, seg_start, seg_end)
setkey(fst2,  Scaffold, snp_start, snp_end)

merged <- foverlaps(
  x = fst2,
  y = tmrca2,
  by.x = c("Scaffold", "snp_start", "snp_end"),
  by.y = c("CHROM", "seg_start", "seg_end"),
  type = "within",
  nomatch = NA,
  mult = "first"
)

merged[, c("snp_start", "snp_end") := NULL]

if ("mean_tmrca" %in% names(merged)) {
  setcolorder(
    merged,
    c("Scaffold", "POS", "mean_tmrca",
      setdiff(names(merged), c("Scaffold", "POS", "mean_tmrca")))
  )
}

## ---------------------------
## 5) Write merged TSV (UNCHANGED)
## ---------------------------
fwrite(merged, out_merged_tsv, sep = "\t", quote = FALSE, na = "NA")
message("Wrote merged file: ", out_merged_tsv)

## ---------------------------
## 6) Filtering options (UPDATED)
##    These filters affect ONLY plotting/regression + filtered TSV outputs.
##    The original merged TSV remains unchanged.
##
## NOTE:
##  - TMRCA filters apply to mean_tmrca values assigned to SNPs after merge.
## ---------------------------
filter_opts <- list(
  # Keep only SNPs that have a matched segment (mean_tmrca not NA)
  keep_only_nonmissing_tmrca = TRUE,
  
  # TMRCA threshold options:
  # - set tmrca_min to 0 to enforce mean_tmrca > 0 (or >= 0 depending on inclusive flags)
  # - set to NULL to disable that bound
  tmrca_min = 0,       # e.g. 0, 100, 1000, NULL
  tmrca_max = NULL,    # e.g. 50000, NULL
  
  # Control inclusivity of bounds:
  tmrca_min_inclusive = FALSE, # FALSE => mean_tmrca > tmrca_min ; TRUE => >=
  tmrca_max_inclusive = TRUE,  # FALSE => mean_tmrca < tmrca_max ; TRUE => <=
  
  # FST threshold options (applied to RAW FST column values before any transform)
  fst_min_threshold = NULL,    # e.g. 0.05 ; NULL disables
  fst_min_inclusive = TRUE,    # TRUE => >= ; FALSE => >
  
  # Write filtered TSV used for each plot/regression
  write_filtered_tsv = TRUE
)

## ---------------------------
## 7) Plot controls (EDIT THESE)
## ---------------------------
plot_opts <- list(
  # Points
  point_color = "#2C7FB8",
  point_alpha = 0.35,
  point_size  = 1.2,
  
  # Regression line + CI ribbon
  line_color  = "#D95F0E",
  line_size   = 1.0,
  show_ci     = TRUE,
  ci_level    = 0.95,
  ci_alpha    = 0.18,
  
  # Text / theme
  base_family      = "Helvetica",
  base_size        = 12,
  title_size       = 14,
  axis_title_size  = 12,
  axis_text_size   = 10,
  
  # Annotation position
  ann_x      = Inf,
  ann_y      = Inf,
  ann_hjust  = 1.05,
  ann_vjust  = 1.20,
  ann_size   = 4.2,
  
  # Axis labels
  xlab = "mean TMRCA",
  ylab = "FST",
  
  # Optional transforms (affect modeling + axes)
  log10_x = FALSE,
  log10_y = FALSE,
  
  # Show regression formula on plot (optional)
  show_lm_formula = TRUE,
  formula_digits  = 4,
  
  # Output
  width_in  = 6.5,
  height_in = 5.2,
  dpi       = 300
)

## ---------------------------
## 8) Helper: apply tmrca filtering cleanly
## ---------------------------
apply_tmrca_filters <- function(dd,
                                tmrca_min = NULL, tmrca_max = NULL,
                                tmrca_min_inclusive = FALSE,
                                tmrca_max_inclusive = TRUE) {
  
  if (!is.null(tmrca_min)) {
    if (isTRUE(tmrca_min_inclusive)) {
      dd <- dd[mean_tmrca >= tmrca_min]
    } else {
      dd <- dd[mean_tmrca >  tmrca_min]
    }
  }
  
  if (!is.null(tmrca_max)) {
    if (isTRUE(tmrca_max_inclusive)) {
      dd <- dd[mean_tmrca <= tmrca_max]
    } else {
      dd <- dd[mean_tmrca <  tmrca_max]
    }
  }
  
  dd
}

## ---------------------------
## 9) Plotting function (UPDATED)
## ---------------------------
make_assoc_plot <- function(df, y_col, out_png, out_filtered_tsv = NULL,
                            opts = plot_opts,
                            cor_method = "spearman",
                            # filtering controls
                            keep_only_nonmissing_tmrca = TRUE,
                            tmrca_min = NULL, tmrca_max = NULL,
                            tmrca_min_inclusive = FALSE, tmrca_max_inclusive = TRUE,
                            fst_min_threshold = NULL, fst_min_inclusive = TRUE,
                            write_filtered_tsv = TRUE) {
  
  stopifnot("mean_tmrca" %in% names(df))
  stopifnot(y_col %in% names(df))
  
  dd <- copy(df)
  
  ## ---- Filtering BEFORE plot/regression ----
  if (isTRUE(keep_only_nonmissing_tmrca)) {
    dd <- dd[!is.na(mean_tmrca)]
  }
  
  # Keep only finite values for tmrca and this FST column
  dd <- dd[is.finite(mean_tmrca) & is.finite(get(y_col))]
  
  # Apply TMRCA thresholds
  dd <- apply_tmrca_filters(
    dd,
    tmrca_min = tmrca_min,
    tmrca_max = tmrca_max,
    tmrca_min_inclusive = tmrca_min_inclusive,
    tmrca_max_inclusive = tmrca_max_inclusive
  )
  
  # Apply FST threshold (raw FST values)
  if (!is.null(fst_min_threshold)) {
    if (isTRUE(fst_min_inclusive)) {
      dd <- dd[get(y_col) >= fst_min_threshold]
    } else {
      dd <- dd[get(y_col) >  fst_min_threshold]
    }
  }
  
  ## Write filtered TSV (new file), keeping original columns only
  if (isTRUE(write_filtered_tsv)) {
    if (is.null(out_filtered_tsv)) {
      stop("write_filtered_tsv=TRUE but out_filtered_tsv is NULL. Provide a path.")
    }
    fwrite(dd, out_filtered_tsv, sep = "\t", quote = FALSE, na = "NA")
    message("Wrote filtered TSV (used for plot/regression): ", out_filtered_tsv)
  }
  
  ## If too few points, plot points only
  if (nrow(dd) < 3) {
    warning(sprintf("Too few points after filtering for %s (n=%d). Plot will have points only.",
                    y_col, nrow(dd)))
    
    p <- ggplot(dd, aes(x = mean_tmrca, y = get(y_col))) +
      geom_point(alpha = opts$point_alpha, size = opts$point_size, color = opts$point_color) +
      labs(title = paste0(y_col, " vs mean_tmrca (n<3 after filtering)"),
           x = opts$xlab, y = opts$ylab) +
      theme_classic(base_size = opts$base_size, base_family = opts$base_family)
    
    print(p)
    ggsave(out_png, plot = p, width = opts$width_in, height = opts$height_in, dpi = opts$dpi)
    return(invisible(p))
  }
  
  ## ---- Build x/y used for modeling/plot (possibly transformed) ----
  dd_plot <- copy(dd)
  dd_plot[, x := mean_tmrca]
  dd_plot[, y := get(y_col)]
  
  if (isTRUE(opts$log10_x)) {
    dd_plot <- dd_plot[x > 0]
    dd_plot[, x := log10(x)]
  }
  if (isTRUE(opts$log10_y)) {
    dd_plot <- dd_plot[y > 0]
    dd_plot[, y := log10(y)]
  }
  
  if (nrow(dd_plot) < 3) {
    warning(sprintf("Too few points after log transform(s) for %s (n=%d). Plot will have points only.",
                    y_col, nrow(dd_plot)))
    
    p <- ggplot(dd_plot, aes(x = x, y = y)) +
      geom_point(alpha = opts$point_alpha, size = opts$point_size, color = opts$point_color) +
      labs(title = paste0(y_col, " vs mean_tmrca (n<3 after transforms)"),
           x = if (isTRUE(opts$log10_x)) paste0("log10(", opts$xlab, ")") else opts$xlab,
           y = if (isTRUE(opts$log10_y)) paste0("log10(", opts$ylab, ")") else opts$ylab) +
      theme_classic(base_size = opts$base_size, base_family = opts$base_family)
    
    print(p)
    ggsave(out_png, plot = p, width = opts$width_in, height = opts$height_in, dpi = opts$dpi)
    return(invisible(p))
  }
  
  ## ---- Correlation ----
  ct <- cor.test(dd_plot$x, dd_plot$y, method = cor_method)
  r_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  ## ---- Linear regression ----
  fit <- lm(y ~ x, data = dd_plot)
  r2 <- summary(fit)$r.squared
  
  ## Optional regression formula (in plotted x/y scale; after transforms if enabled)
  formula_line <- ""
  if (isTRUE(opts$show_lm_formula)) {
    b0 <- coef(fit)[1]
    b1 <- coef(fit)[2]
    d  <- opts$formula_digits
    sign_str <- ifelse(b1 >= 0, "+", "-")
    formula_line <- sprintf("\nLM: y = %.*f %s %.*f x", d, b0, sign_str, d, abs(b1))
  }
  
  ## Filter note for annotation
  filter_note <- ""
  if (isTRUE(keep_only_nonmissing_tmrca)) filter_note <- paste0(filter_note, "TMRCA!=NA; ")
  
  if (!is.null(tmrca_min)) {
    op <- if (isTRUE(tmrca_min_inclusive)) ">=" else ">"
    filter_note <- paste0(filter_note, "TMRCA", op, tmrca_min, "; ")
  }
  if (!is.null(tmrca_max)) {
    op <- if (isTRUE(tmrca_max_inclusive)) "<=" else "<"
    filter_note <- paste0(filter_note, "TMRCA", op, tmrca_max, "; ")
  }
  
  if (!is.null(fst_min_threshold)) {
    op <- if (isTRUE(fst_min_inclusive)) ">=" else ">"
    filter_note <- paste0(filter_note, "FST", op, fst_min_threshold, "; ")
  }
  
  if (nzchar(filter_note)) filter_note <- paste0("\nFilters: ", filter_note)
  
  ann_text <- sprintf(
    "%s r = %.3f\np = %.2g\nR² = %.3f\nn = %d%s%s",
    tools::toTitleCase(cor_method), r_val, p_val, r2, nrow(dd_plot),
    formula_line, filter_note
  )
  
  ## ---- Plot ----
  p <- ggplot(dd_plot, aes(x = x, y = y)) +
    geom_point(alpha = opts$point_alpha, size = opts$point_size, color = opts$point_color) +
    geom_smooth(
      method = "lm",
      se     = isTRUE(opts$show_ci),
      level  = opts$ci_level,
      linewidth = opts$line_size,
      color  = opts$line_color,
      alpha  = opts$ci_alpha
    ) +
    labs(
      title = paste0(y_col, " vs mean_tmrca"),
      x     = if (isTRUE(opts$log10_x)) paste0("log10(", opts$xlab, ")") else opts$xlab,
      y     = if (isTRUE(opts$log10_y)) paste0("log10(", opts$ylab, ")") else opts$ylab
    ) +
    annotate(
      "text",
      x = opts$ann_x, y = opts$ann_y,
      label = ann_text,
      hjust = opts$ann_hjust, vjust = opts$ann_vjust,
      size  = opts$ann_size,
      family = opts$base_family
    ) +
    theme_classic(base_size = opts$base_size, base_family = opts$base_family) +
    theme(
      plot.title   = element_text(size = opts$title_size, face = "bold"),
      axis.title   = element_text(size = opts$axis_title_size),
      axis.text    = element_text(size = opts$axis_text_size)
    )
  
  print(p)
  
  ggsave(out_png, plot = p, width = opts$width_in, height = opts$height_in, dpi = opts$dpi)
  
  invisible(p)
}

## ---------------------------
## 10) Make the three plots requested
## ---------------------------
targets <- c("PROC3_vs_PROC4", "PROC3_vs_PROC5", "PROC4_vs_PROC5")
targets <- targets[targets %in% names(merged)]

if (length(targets) == 0) {
  stop("None of the target FST columns were found in the merged data.")
}

for (col in targets) {
  
  ## Plot output
  out_png <- file.path(out_plot_dir, paste0(col, "_vs_meanTMRCA.png"))
  
  ## Filtered TSV output filename with tags
  fst_tag <- if (is.null(filter_opts$fst_min_threshold)) "noFstFilter" else paste0("fstGE", filter_opts$fst_min_threshold)
  tmrca_tag1 <- if (isTRUE(filter_opts$keep_only_nonmissing_tmrca)) "tmrcaNonMissing" else "tmrcaAll"
  
  tmrca_tag2 <- "noTmrcaThresh"
  if (!is.null(filter_opts$tmrca_min) || !is.null(filter_opts$tmrca_max)) {
    lo_op <- if (isTRUE(filter_opts$tmrca_min_inclusive)) "GE" else "GT"
    hi_op <- if (isTRUE(filter_opts$tmrca_max_inclusive)) "LE" else "LT"
    lo <- if (is.null(filter_opts$tmrca_min)) "NA" else paste0(lo_op, filter_opts$tmrca_min)
    hi <- if (is.null(filter_opts$tmrca_max)) "NA" else paste0(hi_op, filter_opts$tmrca_max)
    tmrca_tag2 <- paste0("tmrca_", lo, "_", hi)
  }
  
  out_filtered_tsv <- file.path(out_filtered_dir, paste0(col, "_filtered_", tmrca_tag1, "_", tmrca_tag2, "_", fst_tag, ".tsv"))
  
  make_assoc_plot(
    merged,
    y_col = col,
    out_png = out_png,
    out_filtered_tsv = out_filtered_tsv,
    opts = plot_opts,
    cor_method = "spearman",
    
    keep_only_nonmissing_tmrca = filter_opts$keep_only_nonmissing_tmrca,
    tmrca_min = 200,
    tmrca_max = filter_opts$tmrca_max,
    tmrca_min_inclusive = filter_opts$tmrca_min_inclusive,
    tmrca_max_inclusive = filter_opts$tmrca_max_inclusive,
    
    fst_min_threshold = 0.01,
    fst_min_inclusive = filter_opts$fst_min_inclusive,
    
    write_filtered_tsv = filter_opts$write_filtered_tsv
  )
  
  message("Saved plot: ", out_png)
}

############################################################
## End of script
############################################################
