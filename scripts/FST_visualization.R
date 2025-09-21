# ============================================
# Plot 3 aligned per-site FST tracks by scaffold & region
# ============================================
library(tidyverse)
setwd("~/Desktop/OSU_projects/conifers/LP/FSTcalc/provenances/100_greedy/fst_per_site")
fstano <- data.frame(fread("All_fst_merged_with_anno.tsv",sep="\t") )
fstanooverlap<-fstano[fstano$overlaps_feature=="TRUE",]
fwrite(fstanooverlap,"fst_415_overlapped.tsv",sep="\t")
tbl<-data.frame(table(fstanooverlap$NR.ID))
# ---- Inputs you change ----
in_csv   <- "fst_per_site_merged_rowmax_lt0.15.csv"
scaffold <- "scaffold_1b"        # <- put your scaffold/contig name here (must match CHROM)
start_bp <- 867985277          # <- region start (bp)
end_bp   <- 867988573          # <- region end (b)

# Name the 3 comparisons exactly as they appear as column names in the CSV
# (e.g., "PROC1_vs_PROC3.weir.fst", etc.)
comparisons <- c(
  "PROC3_vs_PROC4",
  "PROC3_vs_PROC5",
  "PROC4_vs_PROC5"
)

# ---- Load & checks ----
df <- readr::read_csv(in_csv, show_col_types = FALSE)

stopifnot(all(c("CHROM","POS") %in% names(df)))
missing_cols <- setdiff(comparisons, names(df))
if (length(missing_cols)) {
  stop("These comparison columns are missing in the CSV: ",
       paste(missing_cols, collapse=", "))
}
if (length(comparisons) != 3) stop("Please provide exactly three comparison columns.")

# ---- Filter region & reshape to long ----
plot_df <- df %>%
  filter(CHROM == scaffold, POS >= start_bp, POS <= end_bp) %>%
  select(CHROM, POS, all_of(comparisons)) %>%
  pivot_longer(cols = all_of(comparisons),
               names_to = "comparison", values_to = "FST") %>%
  mutate(FST = ifelse(is.finite(FST), pmax(FST, 0), 0))  # guard: NA/neg -> 0

# Common y-limit across facets so they align visually
y_max <- max(plot_df$FST, na.rm = TRUE)
y_max <- ifelse(is.finite(y_max) && y_max > 0, y_max, 0.05)


#visualize per contig
library(tidyverse)
library(readxl)

# ---- 1) Load + filter annotation to the same scaffold/region ----
anno_xlsx <- "../../../../lodgepole_pine_assembly/Aria_curated_annotation_1ab.xlsx"

anno <- read_xlsx(anno_xlsx)

# Be flexible about column names (edit if yours differ)
seq_col   <- dplyr::first(intersect(names(anno), c("contig_1ab")))
start_col <- dplyr::first(intersect(names(anno), c("start1ab_pos")))
end_col   <- dplyr::first(intersect(names(anno), c("end1ab_pos")))
id_col    <- dplyr::first(intersect(names(anno), c("NR.ID")))

stopifnot(!is.na(seq_col), !is.na(start_col), !is.na(end_col))

genes_region <- anno %>%
  transmute(
    CHROM = .data[[seq_col]],
    START = as.numeric(.data[[start_col]]),
    END   = as.numeric(.data[[end_col]]),
    label = if (!is.null(id_col)) as.character(.data[[id_col]]) else NA_character_
  ) %>%
  filter(
    CHROM == scaffold,
    END   >= start_bp,
    START <= end_bp
  ) %>%
  mutate(
    x0 = pmax(START, start_bp),
    x1 = pmin(END,   end_bp)
  )

# Duplicate the same gene track for every facet level
facet_levels <- unique(plot_df$comparison)
ann_track <- genes_region %>%
  tidyr::crossing(comparison = facet_levels)

# ---- 2) Build the plot: add a padded y-limit and draw the gene track ----
y_pad <- y_max * 1.15   # extra headroom to draw the gene track
track_y <- y_max * 1.05 # vertical position of the track
lab_y   <- y_max * 1.10 # label position (optional)

p <- ggplot(plot_df, aes(POS, FST, alpha = FST >= 0.25)) +
  geom_point(size = 0.65) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0), guide = "none") +
  geom_hline(yintercept = c(0.15, 0.5), linetype = "dashed", alpha = 0.35) +
  # --- gene track: horizontal segments + (optional) labels ---
  geom_segment(data = ann_track,
               aes(x = x0, xend = x1, y = track_y, yend = track_y),
               inherit.aes = FALSE, linewidth = 1) +
  # Optional: labels above segments (comment out if too busy)
  geom_text(data = ann_track %>% mutate(lbl = ifelse(is.na(label), "", label)),
            aes(x = (x0 + x1)/2, y = lab_y, label = lbl),
            inherit.aes = FALSE, size = 2.8, vjust = 0, check_overlap = TRUE) +
  facet_wrap(~ comparison, ncol = 1, scales = "fixed") +
  coord_cartesian(xlim = c(start_bp, end_bp), ylim = c(0, y_pad)) +
  scale_x_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    title = paste0("Per-site FST on ", scaffold, " (",
                   scales::comma(start_bp), "–", scales::comma(end_bp), ")"),
    x = "Genomic position (bp)",
    y = "FST"
  ) +
  theme_minimal(base_size = 17) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p)


# ---- Optional: save to file ----
ggsave(paste0("XP_010267416.1_Per-site FST on ", scaffold, " (",scales::comma(start_bp), "–", scales::comma(end_bp), ").png"),
       p, width = 9, height = 7, dpi = 300)
