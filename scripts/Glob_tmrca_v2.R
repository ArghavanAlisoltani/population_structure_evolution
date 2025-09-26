# =========================
# Complete R script: read multiple TMRCA files, window them,
# compute young/old flags, and plot with CI ribbons.
# =========================
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/all_tmrca/all_tmrca")

# install.packages(c("data.table","dplyr","ggplot2","scales"))
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# ---------- SETTINGS ----------
glob <- "outargs_scaffold_10_*.tmrca.txt"     # files to read
arg_start_is_zero_based <- TRUE
win <- 100000             # window size (bp)
tmrca_stat <- "mean"      # one of: "mean","mode","lo","hi"
bottom_q <- 0.01          # bottom 1% = "young"
top_q    <- 0.99          # top 1%    = "old"
out_png  <- "S10_tmrca_windowed_extremes.png"
out_pdf  <- "S10_tmrca_windowed_extremes.pdf"

# ---------- HELPERS ----------
read_tm <- function(f){
  dt <- fread(f, header = FALSE, integer64 = "numeric")
  if (ncol(dt) == 7) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_mode","tmrca_lo","tmrca_hi"))
  } else if (ncol(dt) == 6) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_lo","tmrca_hi"))
    dt[, tmrca_mode := NA_real_]
  } else stop(sprintf("Unexpected column count in %s", f))
  
  if (arg_start_is_zero_based) dt[, start := as.numeric(start) + 1] else dt[, start := as.numeric(start)]
  dt[, end := as.numeric(end)]
  dt[, ci_low  := pmin(tmrca_mean, tmrca_lo, tmrca_hi, na.rm = TRUE)]
  dt[, ci_high := pmax(tmrca_mean, tmrca_lo, tmrca_hi, na.rm = TRUE)]
  dt[]
}

split_to_windows <- function(D, win){
  D <- D[, .(chrom, start, end, y, ci_low, ci_high)]
  if (nrow(D) == 0L) return(data.table(chrom=character(), w_start=integer(), w_end=integer(),
                                       cov_sum=integer(), y_wmean=numeric(), ciw_wmean=numeric()))
  D[, wstart := (start - 1L) %/% win * win + 1L]
  D[, wend   := (end   - 1L) %/% win * win + win]
  
  # Expand each segment across overlapped windows
  out <- D[, .(w = seq.int(wstart, wend, by = win)), by = .(chrom, start, end, y, ci_low, ci_high)]
  out[, wx := pmax(w, start)]
  out[, wy := pmin(w + win - 1L, end)]
  out[, cov := pmax(0L, wy - wx + 1L)]
  out[cov > 0L, .(
    cov_sum   = sum(cov),
    y_wmean   = sum(y * cov) / sum(cov),
    ciw_wmean = sum((ci_high - ci_low) * cov) / sum(cov)
  ), by = .(chrom, w_start = w, w_end = w + win - 1L)]
}

# ---------- READ ALL TMRCA ----------
files <- Sys.glob(glob)
stopifnot(length(files) > 0)
tm <- rbindlist(lapply(files, read_tm), use.names = TRUE, fill = TRUE)

# choose center value
tm[, y := dplyr::case_when(
  tmrca_stat == "mean" ~ tmrca_mean,
  tmrca_stat == "mode" ~ ifelse(is.na(tmrca_mode), tmrca_mean, tmrca_mode),
  tmrca_stat == "lo"   ~ ci_low,
  tmrca_stat == "hi"   ~ ci_high,
  TRUE ~ tmrca_mean
)]

# ---------- WINDOWING ----------
winDT <- split_to_windows(tm, win)
stopifnot(all(c("chrom","w_start","w_end","y_wmean") %in% names(winDT)))

# add midpoints (Mb)
winDT <- winDT %>%
  mutate(mid_bp = (w_start + w_end) / 2,
         mid_mb = mid_bp / 1e6)

# ---------- FLAGS (young/mid/old) ----------
winDT <- winDT %>%
  group_by(chrom) %>%
  mutate(
    y_rank = percent_rank(y_wmean),  # 0..1
    flag   = case_when(
      !is.na(y_rank) & y_rank <= bottom_q ~ "young",
      !is.na(y_rank) & y_rank >= top_q    ~ "old",
      TRUE                                ~ "mid"
    )
  ) %>%
  ungroup()

# ---------- PLOT ----------
p <- ggplot(winDT, aes(x = mid_mb, y = y_wmean)) +
  geom_ribbon(aes(ymin = y_wmean - ciw_wmean/2,
                  ymax = y_wmean + ciw_wmean/2),
              alpha = 0.15) +
  geom_line(linewidth = 0.35) +
  geom_point(
    data = dplyr::filter(winDT, flag %in% c("NA","old")),
    aes(color = flag),
    size = 1.8,
    inherit.aes = TRUE
  ) +
  facet_wrap(~ chrom, scales = "free", ncol = 1) +
  scale_y_continuous(trans = "log10", labels = label_number()) +
  scale_color_manual(values = c(young = "#2C7BB6", old = "#D7191C")) +
  labs(
    x = "Position (Mb)",
    y = paste0("TMRCA (", tmrca_stat, ", log10)"),
    title = sprintf("Windowed TMRCA (win=%skb) with young/old extremes", win/1000)
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.minor = element_blank())

print(p)
ggsave(out_png, p, width = 18, height = 10, dpi = 300)
ggsave(out_pdf, p, width = 18, height = 10)
