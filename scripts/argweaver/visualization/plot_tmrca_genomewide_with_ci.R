# plot_tmrca_genomewide_with_ci.R
# Genome-wide TMRCA line plot with confidence ribbons and lots of customization.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ========= USER OPTIONS =========
tmrca_file   <- "all_tmrca_corrected_position.tsv"  # path to your TMRCA table
out_png      <- "tmrca_genomewide.png"
out_pdf      <- "tmrca_genomewide.pdf"

win_bp       <- 1e6   # window size in bp for aggregation (e.g. 1 Mb)
ci_alpha     <- 0.25  # CI ribbon transparency
line_size    <- 0.7   # line thickness
font_size    <- 12    # base font size
ci_method    <- "mean"  # how to combine CIs per window: "mean" or "range"
# =================================

# --- helpers ---------------------------------------------------------------

# find a column name by a set of candidate regex patterns (case-insensitive)
pick_col <- function(nms, patterns) {
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("Could not find a required column matching: ", paste(patterns, collapse = " | "))
}

# order scaffolds by numeric part, then optional suffix (a<b<c)
order_scaffolds <- function(x) {
  # e.g., "scaffold_1a" -> num=1, suf="a"
  num <- suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", x)))
  suf <- sub(".*?[0-9]+(.*)", "\\1", x)
  suf[is.na(suf)] <- ""
  ord <- order(num, suf, x, na.last = TRUE)
  factor(x, levels = unique(x[ord]))
}

# build per-scaffold windows of fixed width
make_windows <- function(df_scaf, win_bp) {
  # df_scaf has columns: scaffold, scaffold_len
  df_scaf %>%
    rowwise() %>%
    do({
      n <- ceiling(.$scaffold_len / win_bp)
      tibble(scaffold = .$scaffold,
             win_start = seq(0, by = win_bp, length.out = n),
             win_end   = pmin(win_start + win_bp, .$scaffold_len))
    }) %>% ungroup()
}

# --- 1) Read & normalize columns ------------------------------------------
DT <- fread(tmrca_file)

nms <- names(DT)

col_scaf <- pick_col(nms, c("^chrom$", "^scaf", "^seqid$"))
col_start <- pick_col(nms, c("^start", "start_tmrca"))
col_end   <- pick_col(nms, c("^end", "end_tmrca"))
col_mean  <- pick_col(nms, c("^mean_tmrca$", "tmrca_mean"))
col_lo    <- pick_col(nms, c("^low(_|)ci$", "^lower_ci$", "low_ci", "lower[_ ]?ci"))
col_hi    <- pick_col(nms, c("^up(_|)ci$", "^upper_ci$", "up_ci", "upper[_ ]?ci"))

DT <- DT %>%
  rename(scaffold = !!col_scaf,
         seg_start = !!col_start,
         seg_end   = !!col_end,
         tmrca     = !!col_mean,
         lo_ci     = !!col_lo,
         hi_ci     = !!col_hi)

# Ensure numeric & 0-based to 1-based sanity
DT <- DT %>%
  mutate(seg_start = as.numeric(seg_start),
         seg_end   = as.numeric(seg_end),
         tmrca     = as.numeric(tmrca),
         lo_ci     = as.numeric(lo_ci),
         hi_ci     = as.numeric(hi_ci)) %>%
  filter(is.finite(seg_start), is.finite(seg_end), is.finite(tmrca)) %>%
  mutate(seg_len = pmax(0, seg_end - seg_start))  # non-negative

# --- 2) Scaffold lengths & ordering ---------------------------------------
scaf_len <- DT %>%
  group_by(scaffold) %>%
  summarise(scaffold_len = max(seg_end, na.rm = TRUE), .groups = "drop") %>%
  mutate(scaffold_ord = order_scaffolds(scaffold))

# cumulative offsets across genome (for plotting x in bp)
scaf_len <- scaf_len %>%
  arrange(scaffold_ord) %>%
  mutate(cum_offset = c(0, head(cumsum(scaffold_len), -1)))

# merge offsets back to segments
DT <- DT %>% inner_join(scaf_len, by = "scaffold") %>%
  mutate(seg_start_cum = cum_offset + seg_start,
         seg_end_cum   = cum_offset + seg_end)

# --- 3) Bin into fixed windows and compute weighted means ------------------
# Build windows per scaffold
WIN <- make_windows(scaf_len, win_bp)

# Weighted overlap within windows (data.table non-equi join)
DT_seg <- as.data.table(DT[, c("scaffold","seg_start","seg_end","tmrca","lo_ci","hi_ci","seg_len")])
setkey(DT_seg, scaffold, seg_start, seg_end)

WIN_dt <- as.data.table(WIN)
setnames(WIN_dt, c("scaffold","win_start","win_end"))
setkey(WIN_dt, scaffold, win_start, win_end)

# non-equi overlap
ov <- foverlaps(DT_seg, WIN_dt, by.x = c("scaffold","seg_start","seg_end"),
                by.y = c("scaffold","win_start","win_end"),
                nomatch = 0L)

# overlap length within the window
ov[, overlap := pmax(0, pmin(seg_end, win_end) - pmax(seg_start, win_start))]
ov <- ov[overlap > 0]

agg <- ov[, .(
  w_tmrca = sum(tmrca * overlap, na.rm = TRUE) / sum(overlap),
  # two options for CIs: mean of CIs across overlap, or range (min/max)
  lo = if (ci_method == "mean") sum(lo_ci * overlap, na.rm = TRUE) / sum(overlap) else min(lo_ci, na.rm = TRUE),
  hi = if (ci_method == "mean") sum(hi_ci * overlap, na.rm = TRUE) / sum(overlap) else max(hi_ci, na.rm = TRUE)
), by = .(scaffold, win_start, win_end)]

# add cumulative x (use window midpoint)
agg <- as_tibble(agg) %>%
  left_join(scaf_len[, c("scaffold","cum_offset")], by = "scaffold") %>%
  mutate(win_mid_cum = cum_offset + (win_start + win_end)/2)

# --- 4) Plot ---------------------------------------------------------------
# X label formatter in Gb
label_gb <- function(x) sprintf("%.2f", x/1e9)

g <- ggplot(agg, aes(x = win_mid_cum, y = w_tmrca)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey60", alpha = ci_alpha, colour = NA) +
  geom_line(size = line_size, colour = "steelblue") +
  scale_x_continuous(
    name = "Genome position (Gb)",
    labels = label_gb,
    breaks = pretty(agg$win_mid_cum, n = 10)
  ) +
  scale_y_continuous(name = "TMRCA (generations)", labels = label_comma()) +
  theme_minimal(base_size = font_size) +
  theme(panel.grid.minor = element_blank())

print(g)

ggsave(out_png, g, width = 12, height = 4, units = "in", dpi = 300)
ggsave(out_pdf, g, width = 12, height = 4, units = "in")
message("Saved: ", out_png, " and ", out_pdf)

# --- 5) Optional: write the binned data (for reuse) ------------------------
fwrite(agg, "tmrca_binned_windows.tsv", sep = "\t")
