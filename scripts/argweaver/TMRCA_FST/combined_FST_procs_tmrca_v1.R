############################################################
## Compare 3 FST contrasts vs TMRCA (High-FST only > 0.15)
## - Scatter (faceted) + barplots (composition + enrichment)
## Run in RStudio
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

## Optional for non-overlapping labels
use_ggrepel <- TRUE
if (use_ggrepel && !"ggrepel" %in% rownames(installed.packages())) {
  install.packages("ggrepel")
}
if (use_ggrepel) library(ggrepel)

## ---------------------------
## 1) Inputs / Outputs
## ---------------------------
merged_file <- "FST_plus_meanTMRCA.tsv"  # produced by your merge script

out_dir <- "compare_3contrasts_highFST_tmrcaClasses"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------
## 2) User controls (EDIT THESE)
## ---------------------------
opts <- list(
  ## High-FST definition (fixed cutoff)
  high_fst_cutoff = 0.15,        # keep only SNPs with FST > 0.15 for each contrast
  high_fst_inclusive = FALSE,    # FALSE: > 0.15 ; TRUE: >= 0.15
  
  ## TMRCA class thresholds (quantile-based)
  tmrca_low_pct  = 0.05,         # bottom 5%  => LowTMRCA
  tmrca_high_pct = 0.05,         # top 5%     => HighTMRCA
  
  ## Basic filtering of TMRCA values
  keep_only_nonmissing_tmrca = TRUE,
  tmrca_min = 0,                # set NULL to disable
  tmrca_min_inclusive = FALSE,  # FALSE: >0 ; TRUE: >=0
  
  ## Plot aesthetics
  point_alpha = 0.30,
  point_size  = 1.0,
  base_family = "Helvetica",
  base_size   = 12,
  
  ## Optional labeling of top hits (High-FST only)
  label_top_n = 0,              # set 0 to disable labels; e.g. 15 to label top 15 per contrast
  label_drop_scaffold_prefix = TRUE,  # remove "scaffold_" from labels
  label_show_pos_only = FALSE,         # if TRUE label="POS" only; else "scaffold:POS"
  label_segment_color = NA,            # NA removes connector lines; set e.g. "grey40" to show
  
  ## Colors for classes
  tmrca_colors = c(HighTMRCA = "#F8766D", MidTMRCA = "#619CFF", LowTMRCA = "#00BA38"),
  
  ## Output image sizes
  dpi = 300
)

## ---------------------------
## 3) Load data
## ---------------------------
dt <- fread(merged_file)

stopifnot(all(c("Scaffold","POS","mean_tmrca") %in% names(dt)))

contrasts <- c("PROC3_vs_PROC4", "PROC3_vs_PROC5", "PROC4_vs_PROC5")
contrasts <- contrasts[contrasts %in% names(dt)]
if (length(contrasts) == 0) stop("No PROC*_vs_PROC* columns found in merged file.")

## Ensure numeric
dt[, POS := as.numeric(POS)]
dt[, mean_tmrca := as.numeric(mean_tmrca)]
for (cc in contrasts) dt[, (cc) := as.numeric(get(cc))]

## unique SNP id (so we can track rows after melting)
dt[, SNP_ID := paste(Scaffold, POS, sep=":")]

## ---------------------------
## 4) Long format (one row per SNP per contrast)
## ---------------------------
long <- melt(
  dt,
  id.vars = c("SNP_ID","Scaffold","POS","mean_tmrca"),
  measure.vars = contrasts,
  variable.name = "contrast",
  value.name = "FST"
)

## ---------------------------
## 5) Global filtering (applies to everything below)
## ---------------------------
if (isTRUE(opts$keep_only_nonmissing_tmrca)) {
  long <- long[!is.na(mean_tmrca)]
}
long <- long[is.finite(mean_tmrca) & is.finite(FST)]

if (!is.null(opts$tmrca_min)) {
  if (isTRUE(opts$tmrca_min_inclusive)) long <- long[mean_tmrca >= opts$tmrca_min]
  else                                  long <- long[mean_tmrca >  opts$tmrca_min]
}

## ---------------------------
## 6) Define TMRCA classes ONCE (using all SNPs with valid TMRCA)
##     so contrasts are comparable (same thresholds).
## ---------------------------
tmrca_low_thr  <- as.numeric(quantile(long$mean_tmrca, probs = opts$tmrca_low_pct,  na.rm = TRUE, names = FALSE))
tmrca_high_thr <- as.numeric(quantile(long$mean_tmrca, probs = 1 - opts$tmrca_high_pct, na.rm = TRUE, names = FALSE))

long[, tmrca_class := fifelse(mean_tmrca >= tmrca_high_thr, "HighTMRCA",
                              fifelse(mean_tmrca <= tmrca_low_thr,  "LowTMRCA", "MidTMRCA"))]
long[, tmrca_class := factor(tmrca_class, levels = c("LowTMRCA","MidTMRCA","HighTMRCA"))]

## ---------------------------
## 7) Keep HIGH-FST only (per contrast) for comparison
## ---------------------------
if (isTRUE(opts$high_fst_inclusive)) {
  high <- long[FST >= opts$high_fst_cutoff]
} else {
  high <- long[FST >  opts$high_fst_cutoff]
}

## Write out the high-FST-only long table
fwrite(high, file.path(out_dir, "HighFST_only_long.tsv"), sep="\t", quote=FALSE, na="NA")

## If a contrast has zero high-FST points, keep going but warn
n_by_contrast <- high[, .N, by=contrast]
print(n_by_contrast)

## ---------------------------
## 8) (Figure A) Scatter: High-FST only, faceted by contrast, colored by TMRCA class
## ---------------------------
p_scatter <- ggplot(high, aes(x = mean_tmrca, y = FST, color = tmrca_class)) +
  geom_point(alpha = opts$point_alpha, size = opts$point_size) +
  geom_vline(xintercept = tmrca_low_thr,  linetype = "dashed") +
  geom_vline(xintercept = tmrca_high_thr, linetype = "dashed") +
  geom_hline(yintercept = opts$high_fst_cutoff, linetype = "dashed") +
  facet_wrap(~contrast, nrow = 1) +
  scale_color_manual(values = opts$tmrca_colors) +
  labs(
    title = paste0("High-FST only (", ifelse(opts$high_fst_inclusive, "≥", ">"), opts$high_fst_cutoff,
                   "): FST vs mean_tmrca across contrasts"),
    subtitle = paste0("TMRCA classes from global thresholds: Low ≤ ", signif(tmrca_low_thr, 5),
                      "; High ≥ ", signif(tmrca_high_thr, 5)),
    x = "mean_tmrca",
    y = "FST",
    color = "TMRCA class"
  ) +
  theme_classic(base_size = opts$base_size, base_family = opts$base_family)

## Optional labels: top N by FST within each contrast (High-FST only)
if (opts$label_top_n > 0) {
  lab <- high[order(-FST), head(.SD, opts$label_top_n), by=contrast]
  
  lab[, scaffold_clean := if (isTRUE(opts$label_drop_scaffold_prefix)) sub("^scaffold_", "", Scaffold) else Scaffold]
  lab[, label := if (isTRUE(opts$label_show_pos_only)) as.character(POS) else paste0(scaffold_clean, ":", POS)]
  
  if (use_ggrepel) {
    p_scatter <- p_scatter +
      ggrepel::geom_text_repel(
        data = lab,
        aes(label = label),
        size = 3.2,
        max.overlaps = Inf,
        segment.color = opts$label_segment_color
      )
  } else {
    p_scatter <- p_scatter +
      geom_text(data = lab, aes(label = label), size = 3.2, vjust = -0.5)
  }
}

print(p_scatter)
ggsave(file.path(out_dir, "FigA_scatter_highFST_facets.png"),
       p_scatter, width = 12, height = 4.8, dpi = opts$dpi)

## ---------------------------
## 9) (Figure B) Barplot: composition of TMRCA classes among High-FST SNPs (per contrast)
## ---------------------------
comp <- high[, .N, by=.(contrast, tmrca_class)]
comp[, prop := N / sum(N), by=contrast]

fwrite(comp, file.path(out_dir, "HighFST_tmrcaClass_composition.tsv"), sep="\t", quote=FALSE, na="NA")

p_bar_prop <- ggplot(comp, aes(x = contrast, y = prop, fill = tmrca_class)) +
  geom_col(position = "fill") +   # shows proportions (0..1)
  scale_fill_manual(values = opts$tmrca_colors) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  labs(
    title = paste0("TMRCA class composition among High-FST SNPs (", ifelse(opts$high_fst_inclusive, "≥", ">"), opts$high_fst_cutoff, ")"),
    x = "FST contrast",
    y = "Proportion of High-FST SNPs",
    fill = "TMRCA class"
  ) +
  theme_classic(base_size = opts$base_size, base_family = opts$base_family)

print(p_bar_prop)
ggsave(file.path(out_dir, "FigB_bar_tmrcaComposition_highFST.png"),
       p_bar_prop, width = 8.5, height = 5.0, dpi = opts$dpi)

## ---------------------------
## 10) (Optional Figure C) Enrichment: High-FST vs background for TMRCA classes
##     (Fold-enrichment = (HighFST proportion) / (Background proportion))
## ---------------------------
bg <- long[, .N, by=.(contrast, tmrca_class)]
bg[, bg_prop := N / sum(N), by=contrast]

hi <- comp[, .(contrast, tmrca_class, hi_N = N, hi_prop = prop)]

enr <- merge(hi, bg[, .(contrast, tmrca_class, bg_prop)], by=c("contrast","tmrca_class"), all.x=TRUE)
enr[, fold_enrichment := hi_prop / bg_prop]

fwrite(enr, file.path(out_dir, "HighFST_tmrcaClass_enrichment.tsv"), sep="\t", quote=FALSE, na="NA")

p_enr <- ggplot(enr, aes(x = tmrca_class, y = fold_enrichment, fill = tmrca_class)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_manual(values = opts$tmrca_colors, guide="none") +
  labs(
    title = "Fold-enrichment of TMRCA classes among High-FST SNPs (vs background)",
    subtitle = "Fold-enrichment = (High-FST class proportion) / (All SNPs proportion)",
    x = "TMRCA class",
    y = "Fold enrichment"
  ) +
  theme_classic(base_size = opts$base_size, base_family = opts$base_family)

print(p_enr)
ggsave(file.path(out_dir, "FigC_foldEnrichment_highFST_vs_background.png"),
       p_enr, width = 12, height = 4.8, dpi = opts$dpi)

message("Done. Outputs in: ", out_dir)

##########################################################################################
#statistics
##########################################################################################

library(data.table)

merged_file <- "FST_plus_meanTMRCA.tsv"
dt <- fread(merged_file)

contrasts <- c("PROC3_vs_PROC4","PROC3_vs_PROC5","PROC4_vs_PROC5")
contrasts <- contrasts[contrasts %in% names(dt)]

dt[, POS := as.numeric(POS)]
dt[, mean_tmrca := as.numeric(mean_tmrca)]
for (cc in contrasts) dt[, (cc) := as.numeric(get(cc))]
dt[, SNP_ID := paste(Scaffold, POS, sep=":")]

# --- Make long ---
long <- melt(
  dt,
  id.vars = c("SNP_ID","Scaffold","POS","mean_tmrca"),
  measure.vars = contrasts,
  variable.name = "contrast",
  value.name = "FST"
)

# Basic filters (match what you plotted)
long <- long[is.finite(mean_tmrca) & is.finite(FST)]
long <- long[mean_tmrca > 0]  # if you used >0

# Define global TMRCA classes ONCE (same thresholds across contrasts)
tmrca_low_pct  <- 0.05
tmrca_high_pct <- 0.05
tmrca_low_thr  <- as.numeric(quantile(long$mean_tmrca, probs = tmrca_low_pct, na.rm=TRUE))
tmrca_high_thr <- as.numeric(quantile(long$mean_tmrca, probs = 1-tmrca_high_pct, na.rm=TRUE))

long[, tmrca_class := fifelse(mean_tmrca >= tmrca_high_thr, "HighTMRCA",
                              fifelse(mean_tmrca <= tmrca_low_thr,  "LowTMRCA", "MidTMRCA"))]
long[, tmrca_class := factor(tmrca_class, levels=c("LowTMRCA","MidTMRCA","HighTMRCA"))]

# High-FST indicator
fst_cut <- 0.15
long[, HighFST := (FST > fst_cut)]

# -------------------------
# Test A: 2x3 (HighFST x tmrca_class) per contrast
# -------------------------
res_2x3 <- long[, {
  tab <- table(HighFST, tmrca_class)
  # chi-square is usually fine with large counts
  chi <- suppressWarnings(chisq.test(tab))
  # effect size: Cramer's V
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  cramerV <- sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
  list(
    chisq = as.numeric(chi$statistic),
    df = as.numeric(chi$parameter),
    p = chi$p.value,
    cramerV = cramerV
  )
}, by=contrast]

res_2x3[, p_adj_BH := p.adjust(p, method="BH")]
print(res_2x3)

# -------------------------
# Test B: Focused 2x2 (HighTMRCA vs not) per contrast with odds ratio
# -------------------------
res_2x2 <- long[, {
  tab <- table(HighFST, HighTMRCA = (tmrca_class=="HighTMRCA"))
  ft <- fisher.test(tab)  # gives OR and CI
  list(
    OR = unname(ft$estimate),
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p = ft$p.value
  )
}, by=contrast]

res_2x2[, p_adj_BH := p.adjust(p, method="BH")]
print(res_2x2)

# -------------------------
# Test C: Logistic regression across all contrasts
# (Does tmrca_class effect differ by contrast?)
# -------------------------
fit_glm <- glm(HighFST ~ tmrca_class * contrast, data=long, family=binomial())
print(summary(fit_glm))

# Key term(s): tmrca_class:contrast
# If interaction not significant -> similar enrichment pattern across contrasts
