# Genome-wide TMRCA line plot with CI ribbons, with strict scaffold ordering (1a, 1b, 2, ...)
# Input columns expected (case-insensitive ok):
# CHROM, start_tmrca, end_tmrca, mean_tmrca, Low_CI, UP_CI, seg_length

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ================== USER OPTIONS ==================
tmrca_file <- "all_tmrca_corrected_position.tsv"
out_png    <- "tmrca_genomewide_v2.png"
out_pdf    <- "tmrca_genomewide_v2.pdf"

win_bp        <- 1e6    # window size in bp
ci_alpha      <- 0.25   # CI ribbon alpha
line_size     <- 0.7
font_size     <- 12
ci_method     <- "mean" # "mean" or "range"
x_axis_mode   <- "gb"   # "gb" (numerical Gb axis) or "scaffold" (labels at scaffold midpoints)
# ==================================================

# ---- helpers ----
pick_col <- function(nms, patterns){
  for (p in patterns){
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("Could not find column matching: ", paste(patterns, collapse=" | "))
}

parse_scaffold <- function(x){
  # extract numeric and optional suffix after the digits
  num <- suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", x)))
  suf <- sub(".*?[0-9]+(.*)", "\\1", x)
  suf[is.na(suf)] <- ""
  data.frame(scaffold = x, sc_num = num, sc_suf = suf, stringsAsFactors = FALSE)
}

make_windows <- function(scaf_tbl, win_bp){
  scaf_tbl %>%
    rowwise() %>%
    do({
      n <- max(1, ceiling(.$scaffold_len / win_bp))
      tibble(scaffold = .$scaffold,
             win_start = seq(0, by = win_bp, length.out = n),
             win_end   = pmin(win_start + win_bp, .$scaffold_len))
    }) %>% ungroup()
}

# ---- read ----
DT <- fread(tmrca_file)
nms <- names(DT)

col_scaf <- pick_col(nms, c("^chrom$", "^scaf", "^seqid$"))
col_start <- pick_col(nms, c("^start", "start_tmrca"))
col_end   <- pick_col(nms, c("^end", "end_tmrca"))
col_mean  <- pick_col(nms, c("^mean_tmrca$", "tmrca_mean"))
col_lo    <- pick_col(nms, c("^low(_|)ci$", "^lower_ci$", "low_ci"))
col_hi    <- pick_col(nms, c("^up(_|)ci$", "^upper_ci$", "up_ci"))

DT <- DT %>%
  rename(scaffold = !!col_scaf,
         seg_start = !!col_start,
         seg_end   = !!col_end,
         tmrca     = !!col_mean,
         lo_ci     = !!col_lo,
         hi_ci     = !!col_hi) %>%
  mutate(seg_start = as.numeric(seg_start),
         seg_end   = as.numeric(seg_end),
         tmrca     = as.numeric(tmrca),
         lo_ci     = as.numeric(lo_ci),
         hi_ci     = as.numeric(hi_ci)) %>%
  filter(is.finite(seg_start), is.finite(seg_end), is.finite(tmrca)) %>%
  mutate(seg_len = pmax(0, seg_end - seg_start))

# ---- strict scaffold order: integer part then suffix ----
scaf_len <- DT %>%
  group_by(scaffold) %>%
  summarise(scaffold_len = max(seg_end, na.rm = TRUE), .groups = "drop")

ord <- parse_scaffold(scaf_len$scaffold) %>%
  right_join(scaf_len, by = "scaffold") %>%
  arrange(sc_num, sc_suf, scaffold)

ord$scaffold <- factor(ord$scaffold, levels = ord$scaffold)
ord <- ord %>%
  mutate(cum_offset = c(0, head(cumsum(scaffold_len), -1)),
         mid_cum    = cum_offset + scaffold_len/2)

# apply order to segments + cumulative coordinates
DT <- DT %>%
  mutate(scaffold = factor(scaffold, levels = levels(ord$scaffold))) %>%
  inner_join(ord[, c("scaffold","cum_offset")], by = "scaffold") %>%
  mutate(seg_start_cum = cum_offset + seg_start,
         seg_end_cum   = cum_offset + seg_end)

# ---- bin into fixed-width windows per scaffold ----
WIN <- make_windows(ord[, c("scaffold","scaffold_len")], win_bp)
DT_seg <- as.data.table(DT[, c("scaffold","seg_start","seg_end","tmrca","lo_ci","hi_ci")])
setkey(DT_seg, scaffold, seg_start, seg_end)
WIN_dt <- as.data.table(WIN)
setnames(WIN_dt, c("scaffold","win_start","win_end"))
setkey(WIN_dt, scaffold, win_start, win_end)

ov <- foverlaps(DT_seg, WIN_dt,
                by.x = c("scaffold","seg_start","seg_end"),
                by.y = c("scaffold","win_start","win_end"),
                nomatch = 0L)

ov[, overlap := pmax(0, pmin(seg_end, win_end) - pmax(seg_start, win_start))]
ov <- ov[overlap > 0]

agg <- ov[, .(
  w_tmrca = sum(tmrca * overlap) / sum(overlap),
  lo = if (ci_method == "mean") sum(lo_ci * overlap) / sum(overlap) else min(lo_ci),
  hi = if (ci_method == "mean") sum(hi_ci * overlap) / sum(overlap) else max(hi_ci)
), by = .(scaffold, win_start, win_end)]

agg <- as_tibble(agg) %>%
  left_join(ord[, c("scaffold","cum_offset","mid_cum","scaffold_len")], by = "scaffold") %>%
  mutate(win_mid_cum = cum_offset + (win_start + win_end)/2)

# ---- plot ----
label_gb <- function(x) sprintf("%.2f", x/1e9)

# scaffold boundary lines
vlines <- ord$cum_offset
# remove 0 if you don’t want a line at the origin
vlines <- vlines[vlines > 0]

p <- ggplot(agg, aes(x = win_mid_cum, y = w_tmrca)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey60", alpha = ci_alpha, colour = NA) +
  geom_line(size = line_size, colour = "black") +
  geom_vline(xintercept = vlines, colour = "grey85", size = 0.25)

if (x_axis_mode == "gb") {
  p <- p +
    scale_x_continuous(
      name = "Genome position (Gb)",
      labels = label_gb,
      breaks = pretty(agg$win_mid_cum, n = 10)
    )
} else {
  # x ticks at scaffold midpoints, labeled by scaffold id
  p <- p +
    scale_x_continuous(
      name = "Scaffolds (ordered: 1a, 1b, 2, ...)",
      breaks = ord$mid_cum,
      labels = as.character(ord$scaffold),
      expand = expansion(mult = c(0.001, 0.001))
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

p <- p +
  scale_y_continuous(name = "TMRCA (generations)", labels = label_comma()) +
  theme_minimal(base_size = font_size) +
  theme(panel.grid.minor = element_blank())

ggsave(out_png, p, width = 13, height = 4.5, units = "in", dpi = 300)
ggsave(out_pdf, p, width = 13, height = 4.5, units = "in")

message("Saved: ", out_png, " and ", out_pdf)

