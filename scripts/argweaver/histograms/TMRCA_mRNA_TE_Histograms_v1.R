# ================================
# 1. Setup
# ================================
library(tidyverse)
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025")
# ---- User-tunable options ----

# Bins for mean_tmrca_at_pos (generations)
tmrca_breaks <- c(0, 8, 20, 100, 1000, 2000, 5000, 10000, 30000, 60000, Inf)

# Optional custom labels (same length as breaks-1)
tmrca_labels <- c(
  "0–49",
  "50–99",
  "100–499",
  "500–999",
  "1000–1999",
  "2000–4999",
  "5000–9999",
  "10000–29999",
  "30000–59999",
  "60000+"
)

# Colors (change as you like)
hist_fill_color   <- "grey70"
mrna_line_color   <- "red3"
te_line_color     <- "blue3"

# Input and output files
in_file  <- "snp_annotations_top_5_percentile/no_NA_snps_annotated_tmrca_te_mrna.tsv"
out_tsv  <- "snps_tmrca_bins_summary.tsv"
out_hist_counts      <- "tmrca_hist_counts.png"
out_hist_mrna_combo  <- "tmrca_hist_mrna_combo.png"
out_hist_te_combo    <- "tmrca_hist_te_combo.png"

# ================================
# 2. Read data
# ================================
snps <- readr::read_tsv(in_file, show_col_types = FALSE)

# Make logical indicators for overlaps
snps <- snps %>%
  mutate(
    mrna_overlap = !is.na(mrna_id) & mrna_id != "",
    te_overlap   = as.logical(in_TE)  # in case it's 0/1 or "TRUE"/"FALSE"
  )

# ================================
# 3. Bin TMRCA and summarise
# ================================
snps_binned <- snps %>%
  mutate(
    tmrca_bin = cut(
      mean_tmrca_at_pos,
      breaks  = tmrca_breaks,
      labels  = tmrca_labels,
      include.lowest = TRUE,
      right   = FALSE  # [lower, upper)
    )
  ) %>%
  group_by(tmrca_bin) %>%
  summarise(
    n_snps   = n(),
    n_mrna   = sum(mrna_overlap, na.rm = TRUE),
    n_te     = sum(te_overlap,   na.rm = TRUE),
    prop_mrna = if_else(n_snps > 0, n_mrna / n_snps, NA_real_),
    prop_te   = if_else(n_snps > 0, n_te   / n_snps, NA_real_),
    .groups = "drop"
  )

# Write TSV with all bin stats
readr::write_tsv(snps_binned, out_tsv)

# For scaling the secondary y axes
max_count <- max(snps_binned$n_snps, na.rm = TRUE)

# ================================
# 4. Histogram of SNP counts per TMRCA bin
# ================================
p_counts <- ggplot(snps_binned, aes(x = tmrca_bin, y = n_snps)) +
  geom_col(fill = hist_fill_color) +
  labs(
    x = "Mean TMRCA bin (generations)",
    y = "Number of SNPs",
    title = "Distribution of SNPs across TMRCA bins"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(out_hist_counts, p_counts, width = 8, height = 5, dpi = 300)

# ================================
# 5. Histogram + mRNA line (2nd y-axis)
# ================================
p_mrna <- ggplot(snps_binned, aes(x = tmrca_bin)) +
  geom_col(aes(y = n_snps), fill = hist_fill_color) +
  geom_line(aes(y = prop_mrna * max_count, group = 1),
            color = mrna_line_color, linewidth = 1) +
  geom_point(aes(y = prop_mrna * max_count),
             color = mrna_line_color, size = 2) +
  scale_y_continuous(
    name = "Number of SNPs",
    sec.axis = sec_axis(
      ~ . / max_count,
      name = "Proportion of SNPs overlapping mRNA"
    )
  ) +
  labs(
    x = "Mean TMRCA bin (generations)",
    title = "SNP counts and mRNA overlap across TMRCA bins"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = mrna_line_color),
    axis.text.y.right  = element_text(color = mrna_line_color)
  )

ggsave(out_hist_mrna_combo, p_mrna, width = 8, height = 5, dpi = 300)

# ================================
# 6. Histogram + TE line (2nd y-axis)
# ================================
p_te <- ggplot(snps_binned, aes(x = tmrca_bin)) +
  geom_col(aes(y = n_snps), fill = hist_fill_color) +
  geom_line(aes(y = prop_te * max_count, group = 1),
            color = te_line_color, linewidth = 1) +
  geom_point(aes(y = prop_te * max_count),
             color = te_line_color, size = 2) +
  scale_y_continuous(
    name = "Number of SNPs",
    sec.axis = sec_axis(
      ~ . / max_count,
      name = "Proportion of SNPs in TEs"
    )
  ) +
  labs(
    x = "Mean TMRCA bin (generations)",
    title = "SNP counts and TE overlap across TMRCA bins"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = te_line_color),
    axis.text.y.right  = element_text(color = te_line_color)
  )

ggsave(out_hist_te_combo, p_te, width = 8, height = 5, dpi = 300)
