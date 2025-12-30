############################################################
##  Combine TMRCA + FST beyond regression:
##   Option 1: Joint outlier/quadrant scan
##   Option 3: Stratified (bin TMRCA, compare FST distributions)
##  Run in RStudio
############################################################
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/MAF_FST_TMRCA")

## ---------------------------
## 0) Packages
## ---------------------------
pkgs <- c("data.table", "ggplot2")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(data.table)
library(ggplot2)

## Optional (nice labels that don’t overlap). If not installed, code still works.
use_ggrepel <- TRUE
if (use_ggrepel && !"ggrepel" %in% rownames(installed.packages())) {
  install.packages("ggrepel")
}
if (use_ggrepel) library(ggrepel)

## ---------------------------
## 1) Inputs
## ---------------------------
merged_file <- "FST_plus_meanTMRCA.tsv"  # produced by your merge script

out_dir <- "tmrca_fst_joint_analyses"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- fread(merged_file)

## expected columns:
stopifnot(all(c("Scaffold","POS","mean_tmrca") %in% names(df)))

contrasts <- c("PROC3_vs_PROC4", "PROC3_vs_PROC5", "PROC4_vs_PROC5")
contrasts <- contrasts[contrasts %in% names(df)]
if (length(contrasts) == 0) stop("No contrast columns found in merged file.")

## Ensure numeric types
df[, POS := as.numeric(POS)]
df[, mean_tmrca := as.numeric(mean_tmrca)]
for (cc in contrasts) df[, (cc) := as.numeric(get(cc))]

## ---------------------------
## 2) Global filtering options (applies before BOTH option 1 and 3)
##    These filters remove points from plots + downstream calculations.
## ---------------------------
global_filter <- list(
  keep_only_nonmissing_tmrca = TRUE,   # keep only rows where mean_tmrca is present
  tmrca_min = 0,                       # e.g. 0 to enforce >0 (use inclusive flag below)
  tmrca_max = NULL,                    # e.g. 50000
  tmrca_min_inclusive = FALSE,         # FALSE => > tmrca_min ; TRUE => >=
  tmrca_max_inclusive = TRUE,          # FALSE => < tmrca_max ; TRUE => <=
  
  fst_min = NULL,                      # e.g. 0.00 or 0.05
  fst_max = NULL,                      # e.g. 0.99
  fst_min_inclusive = TRUE,
  fst_max_inclusive = TRUE
)

apply_global_filters <- function(dd, ycol, gf) {
  if (isTRUE(gf$keep_only_nonmissing_tmrca)) dd <- dd[!is.na(mean_tmrca)]
  dd <- dd[is.finite(mean_tmrca) & is.finite(get(ycol))]
  
  # TMRCA bounds
  if (!is.null(gf$tmrca_min)) {
    if (isTRUE(gf$tmrca_min_inclusive)) dd <- dd[mean_tmrca >= gf$tmrca_min]
    else                               dd <- dd[mean_tmrca >  gf$tmrca_min]
  }
  if (!is.null(gf$tmrca_max)) {
    if (isTRUE(gf$tmrca_max_inclusive)) dd <- dd[mean_tmrca <= gf$tmrca_max]
    else                               dd <- dd[mean_tmrca <  gf$tmrca_max]
  }
  
  # FST bounds (on the chosen contrast column)
  if (!is.null(gf$fst_min)) {
    if (isTRUE(gf$fst_min_inclusive)) dd <- dd[get(ycol) >= gf$fst_min]
    else                             dd <- dd[get(ycol) >  gf$fst_min]
  }
  if (!is.null(gf$fst_max)) {
    if (isTRUE(gf$fst_max_inclusive)) dd <- dd[get(ycol) <= gf$fst_max]
    else                             dd <- dd[get(ycol) <  gf$fst_max]
  }
  
  dd
}

## ---------------------------
## 3) Option 1: Joint outlier / Quadrant scan
## ---------------------------
opt1 <- list(
  fst_top_pct   = 0.01,  # "high FST" = top 1% (change to 0.05 for top 5%)
  tmrca_top_pct = 0.05,  # "high TMRCA" = top 5%
  tmrca_low_pct = 0.05,  # "low  TMRCA" = bottom 5% (young)
  
  label_top_n   = 0,    # label top N by joint score
  label_size    = 3.2,
  
  # plot aesthetics
  point_alpha   = 0.35,
  point_size    = 1.1,
  base_family   = "Helvetica",
  base_size     = 12
)

run_option1 <- function(dd, ycol, out_prefix, opt1) {
  ## build working columns
  x <- dd$mean_tmrca
  y <- dd[[ycol]]
  
  # Thresholds from empirical quantiles
  fst_thr_high   <- as.numeric(quantile(y, probs = 1 - opt1$fst_top_pct,   na.rm = TRUE, names = FALSE))
  tmrca_thr_high <- as.numeric(quantile(x, probs = 1 - opt1$tmrca_top_pct, na.rm = TRUE, names = FALSE))
  tmrca_thr_low  <- as.numeric(quantile(x, probs = opt1$tmrca_low_pct,     na.rm = TRUE, names = FALSE))
  
  dd2 <- copy(dd)
  dd2[, fst := get(ycol)]
  dd2[, tmrca := mean_tmrca]
  
  # Quadrant labels using high-FST threshold + (high vs low) TMRCA
  dd2[, fst_class := ifelse(fst >= fst_thr_high, "HighFST", "LowFST")]
  dd2[, tmrca_class := fifelse(tmrca >= tmrca_thr_high, "HighTMRCA",
                               fifelse(tmrca <= tmrca_thr_low,  "LowTMRCA", "MidTMRCA"))]
  
  # A simple, robust joint score for ranking:
  # rank product favors points that are extreme in BOTH dimensions.
  dd2[, rank_fst   := frank(-fst, ties.method = "average")]   # high fst => small rank
  dd2[, rank_tmrca := frank(-tmrca, ties.method = "average")] # high tmrca => small rank
  dd2[, joint_rankprod := rank_fst * rank_tmrca]
  
  # Candidates of interest (choose whichever category you want):
  # 1) Old + differentiated: HighFST & HighTMRCA
  cand_old_diff <- dd2[fst >= fst_thr_high & tmrca >= tmrca_thr_high]
  # 2) Young + differentiated: HighFST & LowTMRCA
  cand_young_diff <- dd2[fst >= fst_thr_high & tmrca <= tmrca_thr_low]
  
  # Write out tables
  fwrite(dd2, paste0(out_prefix, "_quadrants_all.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(cand_old_diff,  paste0(out_prefix, "_candidates_HighFST_HighTMRCA.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(cand_young_diff,paste0(out_prefix, "_candidates_HighFST_LowTMRCA.tsv"),  sep = "\t", quote = FALSE, na = "NA")
  
  # Labels: top N by joint score among HighFST candidates (feel free to change)
  lab <- dd2[fst >= fst_thr_high][order(joint_rankprod)][1:min(opt1$label_top_n, .N)]
  lab[, label := paste0(Scaffold, ":", POS)]
  
  # Scatter plot with thresholds
  p_scatter <- ggplot(dd2, aes(x = tmrca, y = fst)) +
    geom_point(aes(color = tmrca_class), alpha = opt1$point_alpha, size = opt1$point_size) +
    geom_vline(xintercept = tmrca_thr_high, linetype = "dashed") +
    geom_vline(xintercept = tmrca_thr_low,  linetype = "dashed") +
    geom_hline(yintercept = fst_thr_high,   linetype = "dashed") +
    labs(
      title = paste0("Option 1: Quadrant scan — ", ycol),
      subtitle = paste0(
        "HighFST = top ", opt1$fst_top_pct*100, "% (≥ ", signif(fst_thr_high, 4), "); ",
        "HighTMRCA = top ", opt1$tmrca_top_pct*100, "% (≥ ", signif(tmrca_thr_high, 4), "); ",
        "LowTMRCA = bottom ", opt1$tmrca_low_pct*100, "% (≤ ", signif(tmrca_thr_low, 4), ")"
      ),
      x = "mean_tmrca",
      y = ycol,
      color = "TMRCA class"
    ) +
    theme_classic(base_size = opt1$base_size, base_family = opt1$base_family)
  
  # Add labels (optional)
  if (use_ggrepel) {
    p_scatter <- p_scatter +
      ggrepel::geom_text_repel(
        data = lab,
        aes(x = tmrca, y = fst, label = label),
        size = opt1$label_size,
        max.overlaps = Inf
      )
  } else {
    p_scatter <- p_scatter +
      geom_text(
        data = lab,
        aes(x = tmrca, y = fst, label = label),
        size = opt1$label_size,
        vjust = -0.5
      )
  }
  
  # Quadrant counts barplot (how many SNPs in HighFST & each tmrca_class)
  q_counts <- dd2[, .N, by = .(fst_class, tmrca_class)]
  p_bar <- ggplot(q_counts, aes(x = tmrca_class, y = N, fill = fst_class)) +
    geom_col(position = "dodge") +
    labs(
      title = paste0("Option 1: Quadrant counts — ", ycol),
      x = "TMRCA class",
      y = "Number of SNPs"
    ) +
    theme_classic(base_size = opt1$base_size, base_family = opt1$base_family)
  
  # Save plots
  ggsave(paste0(out_prefix, "_opt1_scatter_quadrants.png"), p_scatter, width = 7, height = 5.5, dpi = 300)
  ggsave(paste0(out_prefix, "_opt1_quadrant_counts.png"),   p_bar,     width = 6.5, height = 4.8, dpi = 300)
  
  list(
    thresholds = list(fst_high = fst_thr_high, tmrca_high = tmrca_thr_high, tmrca_low = tmrca_thr_low),
    scatter = p_scatter,
    bar = p_bar
  )
}

## ---------------------------
## 4) Option 3: Stratify by TMRCA bins, compare FST distributions
## ---------------------------
opt3 <- list(
  n_bins = 10,          # number of TMRCA bins (quantile-based)
  show_jitter = TRUE,
  jitter_alpha = 0.15,
  jitter_size = 0.7,
  
  base_family = "Helvetica",
  base_size = 12
)

run_option3 <- function(dd, ycol, out_prefix, opt3) {
  dd2 <- copy(dd)
  dd2[, fst := get(ycol)]
  dd2[, tmrca := mean_tmrca]
  
  # Quantile bins (approximately equal counts per bin)
  probs <- seq(0, 1, length.out = opt3$n_bins + 1)
  brks <- unique(as.numeric(quantile(dd2$tmrca, probs = probs, na.rm = TRUE)))
  if (length(brks) < 3) stop("Not enough unique TMRCA values to bin meaningfully.")
  
  dd2[, tmrca_bin := cut(tmrca, breaks = brks, include.lowest = TRUE, right = TRUE)]
  dd2[, tmrca_bin := factor(tmrca_bin, levels = levels(tmrca_bin), ordered = TRUE)]
  
  # Nonparametric test across bins (robust if not normal)
  kw <- kruskal.test(fst ~ tmrca_bin, data = dd2)
  
  # Summary per bin for a trend plot
  sum_bin <- dd2[, .(
    n = .N,
    tmrca_median = median(tmrca, na.rm = TRUE),
    fst_median   = median(fst, na.rm = TRUE),
    fst_mean     = mean(fst, na.rm = TRUE)
  ), by = tmrca_bin]
  
  fwrite(sum_bin, paste0(out_prefix, "_opt3_bin_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
  
  # Violin + boxplot
  p_violin <- ggplot(dd2, aes(x = tmrca_bin, y = fst)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    labs(
      title = paste0("Option 3: FST across TMRCA bins — ", ycol),
      subtitle = paste0("Kruskal–Wallis p = ", signif(kw$p.value, 3), " (n_bins=", opt3$n_bins, ")"),
      x = "TMRCA bin (quantiles)",
      y = ycol
    ) +
    theme_classic(base_size = opt3$base_size, base_family = opt3$base_family) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (isTRUE(opt3$show_jitter)) {
    p_violin <- p_violin +
      geom_jitter(alpha = opt3$jitter_alpha, size = opt3$jitter_size, width = 0.15)
  }
  
  # Trend plot (median FST vs median TMRCA per bin)
  p_trend <- ggplot(sum_bin, aes(x = tmrca_median, y = fst_median)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste0("Option 3: Trend by bin medians — ", ycol),
      x = "Median TMRCA (per bin)",
      y = "Median FST (per bin)"
    ) +
    theme_classic(base_size = opt3$base_size, base_family = opt3$base_family)
  
  ggsave(paste0(out_prefix, "_opt3_violin_box.png"), p_violin, width = 8, height = 5.2, dpi = 300)
  ggsave(paste0(out_prefix, "_opt3_trend_medians.png"), p_trend, width = 6.8, height = 5.0, dpi = 300)
  
  list(kruskal = kw, violin = p_violin, trend = p_trend, bin_summary = sum_bin)
}

## ---------------------------
## 5) Run both options for each contrast
## ---------------------------
for (ycol in contrasts) {
  
  message("Running joint analyses for: ", ycol)
  
  dd <- copy(df)
  dd <- apply_global_filters(dd, ycol, global_filter)
  
  if (nrow(dd) < 50) {
    warning("Very few rows after filtering for ", ycol, " (n=", nrow(dd), "). Results may be unstable.")
  }
  
  out_prefix <- file.path(out_dir, ycol)
  
  ## Option 1
  res1 <- run_option1(dd, ycol, out_prefix, opt1)
  print(res1$scatter)
  print(res1$bar)
  
  ## Option 3
  res3 <- run_option3(dd, ycol, out_prefix, opt3)
  print(res3$violin)
  print(res3$trend)
}

message("Done. Outputs written to: ", out_dir)


############################################################
## chunk 2
############################################################

library(data.table)
library(ggplot2)

# Optional for mixed model test
pkgs2 <- c("lme4", "lmerTest")
to_install2 <- pkgs2[!pkgs2 %in% rownames(installed.packages())]
if (length(to_install2) > 0) install.packages(to_install2)
library(lme4)
library(lmerTest)

# ---- Input ----
merged_file <- "FST_plus_meanTMRCA.tsv"
dt <- fread(merged_file)

contrasts <- c("PROC3_vs_PROC4", "PROC3_vs_PROC5", "PROC4_vs_PROC5")
contrasts <- contrasts[contrasts %in% names(dt)]

# ---- Controls ----
multi_opts <- list(
  tmrca_log10 = FALSE,      # TRUE to use log10(mean_tmrca)
  tmrca_min = 0,            # e.g. 0 to enforce >0
  tmrca_min_inclusive = FALSE,
  
  fst_min = NULL,           # e.g. 0.05
  fst_min_inclusive = TRUE,
  fst_max = 0.99,           # e.g. 0.99 (set NULL to disable)
  fst_max_inclusive = TRUE,
  
  point_alpha = 0.15,
  point_size  = 0.6,
  smooth_method = "lm",     # "lm" or "gam"
  smooth_se = FALSE,
  
  base_family = "Helvetica",
  base_size   = 12,
  
  out_png = "FST_vs_TMRCA_3contrasts_facets.png",
  width_in = 10,
  height_in = 4.2,
  dpi = 300
)

# ---- Prep / long format ----
stopifnot(all(c("Scaffold","POS","mean_tmrca") %in% names(dt)))
dt[, POS := as.numeric(POS)]
dt[, mean_tmrca := as.numeric(mean_tmrca)]
for (cc in contrasts) dt[, (cc) := as.numeric(get(cc))]

dt[, SNP_ID := paste(Scaffold, POS, sep=":")]

long <- melt(
  dt,
  id.vars = c("SNP_ID","Scaffold","POS","mean_tmrca"),
  measure.vars = contrasts,
  variable.name = "contrast",
  value.name = "FST"
)

# ---- Filters (applied to what you plot + analyze) ----
long <- long[is.finite(mean_tmrca) & is.finite(FST)]

# TMRCA min
if (!is.null(multi_opts$tmrca_min)) {
  if (isTRUE(multi_opts$tmrca_min_inclusive)) long <- long[mean_tmrca >= multi_opts$tmrca_min]
  else                                       long <- long[mean_tmrca >  multi_opts$tmrca_min]
}

# FST min/max
if (!is.null(multi_opts$fst_min)) {
  if (isTRUE(multi_opts$fst_min_inclusive)) long <- long[FST >= multi_opts$fst_min]
  else                                     long <- long[FST >  multi_opts$fst_min]
}
if (!is.null(multi_opts$fst_max)) {
  if (isTRUE(multi_opts$fst_max_inclusive)) long <- long[FST <= multi_opts$fst_max]
  else                                     long <- long[FST <  multi_opts$fst_max]
}

# optional log10 TMRCA for plot/model
if (isTRUE(multi_opts$tmrca_log10)) {
  long <- long[mean_tmrca > 0]
  long[, tmrca_x := log10(mean_tmrca)]
  xlab <- "log10(mean_tmrca)"
} else {
  long[, tmrca_x := mean_tmrca]
  xlab <- "mean_tmrca"
}

# ---- Per-contrast correlation text for facets ----
stats <- long[, {
  ct <- cor.test(tmrca_x, FST, method = "spearman")
  list(rho = unname(ct$estimate), p = ct$p.value, n = .N)
}, by = contrast]

stats[, label := sprintf("Spearman rho=%.3f\np=%.2g\nn=%d", rho, p, n)]

# ---- One-figure visualization: 3 facets ----
p <- ggplot(long, aes(x = tmrca_x, y = FST)) +
  geom_point(alpha = multi_opts$point_alpha, size = multi_opts$point_size) +
  geom_smooth(method = multi_opts$smooth_method, se = multi_opts$smooth_se) +
  facet_wrap(~contrast, nrow = 1, scales = "free_y") +
  geom_text(
    data = stats,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.10,
    size = 3.6,
    family = multi_opts$base_family
  ) +
  labs(
    title = "FST vs TMRCA across contrasts (faceted)",
    x = xlab,
    y = "FST"
  ) +
  theme_classic(base_size = multi_opts$base_size, base_family = multi_opts$base_family)

print(p)
ggsave(multi_opts$out_png, p, width = multi_opts$width_in, height = multi_opts$height_in, dpi = multi_opts$dpi)

# ---- Statistical comparison across contrasts (does slope differ by contrast?) ----
# Mixed model accounts for the fact that each SNP contributes 3 measurements (paired).
# Interaction term tests whether the TMRCA-FST relationship differs among contrasts.
m <- lmer(FST ~ tmrca_x * contrast + (1|SNP_ID), data = long)
print(summary(m))
print(anova(m))   # look at tmrca_x:contrast p-value (difference in slopes across contrasts)

