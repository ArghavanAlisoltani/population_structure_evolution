# =========================
# ONE TMRCA file + ONE scaffold region + gene coordinates shown
# =========================
# install.packages(c("data.table","dplyr","ggplot2","scales","readxl","ggrepel"))
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(readxl)
library(ggrepel)

# ---------- USER SETTINGS ----------
tmrca_path <- "all_tmrca/All_60.tmrca.txt"     # path to one TMRCA file
ann_xlsx   <- "~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/Aria_curated_annotation_1ab.xlsx"                       # annotation workbook

contig   <- "scaffold_10"  # scaffold to visualize
x_start  <- 500000000     # start (bp)
x_end    <- 510000000      # end   (bp)

arg_start_is_zero_based <- TRUE
win <- 100000               # window size (bp)
tmrca_stat <- "mean"        # one of: "mean","mode","lo","hi"
bottom_q <- 0.01            # bottom 1% = "young"
top_q    <- 0.99            # top 1%    = "old"
logY <- TRUE                # log10 Y axis

label_genes   <- TRUE       # show gene IDs above lines
max_gene_rows <- Inf        # set to a number to limit labels if crowded

out_png <- paste(contig,"tmrca","region",x_start,x_end,".png",sep = "_")

# ---------- HELPERS ----------
read_tm_one <- function(f, zero_based = TRUE){
  dt <- fread(f, header = FALSE, integer64 = "numeric")
  if (ncol(dt) == 7) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_mode","tmrca_lo","tmrca_hi"))
  } else if (ncol(dt) == 6) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_lo","tmrca_hi"))
    dt[, tmrca_mode := NA_real_]
  } else stop(sprintf("Unexpected column count in %s", f))
  if (zero_based) dt[, start := as.numeric(start) + 1] else dt[, start := as.numeric(start)]
  dt[, end := as.numeric(end)]
  dt[, ci_low  := pmin(tmrca_mean, tmrca_lo, tmrca_hi, na.rm = TRUE)]
  dt[, ci_high := pmax(tmrca_mean, tmrca_lo, tmrca_hi, na.rm = TRUE)]
  dt[]
}

read_annotation_1ab <- function(xlsx, sheet = NULL, id_col = "NR.ID"){
  if (is.null(sheet)) sheet <- excel_sheets(xlsx)[1]
  ann <- read_excel(xlsx, sheet = sheet, guess_max = 200000)
  req <- c("contig_1ab","start1ab_pos","end1ab_pos")
  miss <- setdiff(req, names(ann))
  if (length(miss)) stop("Missing required columns in annotation: ", paste(miss, collapse = ", "))
  if (!id_col %in% names(ann)) ann[[id_col]] <- paste0("GENE_", seq_len(nrow(ann)))
  ann %>%
    transmute(
      chrom      = as.character(contig_1ab),
      gene_start = as.numeric(start1ab_pos),
      gene_end   = as.numeric(end1ab_pos),
      gene_id    = as.character(.data[[id_col]])
    ) %>%
    arrange(chrom, gene_start, gene_end)
}

split_to_windows_region <- function(D, win, x_start, x_end){
  if (nrow(D) == 0L) return(data.table(chrom=character(), w_start=integer(), w_end=integer(),
                                       cov_sum=integer(), y_wmean=numeric(), ciw_wmean=numeric()))
  D <- D[, .(chrom,
             start = pmax(start, x_start),
             end   = pmin(end,   x_end),
             y, ci_low, ci_high)]
  D <- D[end >= start]
  if (nrow(D) == 0L) return(data.table(chrom=character(), w_start=integer(), w_end=integer(),
                                       cov_sum=integer(), y_wmean=numeric(), ciw_wmean=numeric()))
  w0 <- ((x_start - 1L) %/% win) * win + 1L
  wN <- ((x_end   - 1L) %/% win) * win + 1L
  D[, `:=`(wstart = pmax(((start - 1L) %/% win) * win + 1L, w0),
           wend   = pmin(((end   - 1L) %/% win) * win + 1L, wN))]
  out <- D[, .(w = seq.int(wstart, wend, by = win)), by = .(chrom, start, end, y, ci_low, ci_high)]
  if (nrow(out) == 0L) return(data.table(chrom=character(), w_start=integer(), w_end=integer(),
                                         cov_sum=integer(), y_wmean=numeric(), ciw_wmean=numeric()))
  out[, wx := pmax(w, start)]
  out[, wy := pmin(w + win - 1L, end)]
  out[, cov := pmax(0L, wy - wx + 1L)]
  out[cov > 0L, .(
    cov_sum   = sum(cov),
    y_wmean   = sum(y * cov) / sum(cov),
    ciw_wmean = sum((ci_high - ci_low) * cov) / sum(cov)
  ), by = .(chrom, w_start = w, w_end = w + win - 1L)]
}

# ---------- READ ONE TMRCA & FILTER ----------
tm <- read_tm_one(tmrca_path, zero_based = arg_start_is_zero_based)
tm <- tm[chrom == contig & end >= x_start & start <= x_end]
stopifnot(nrow(tm) > 0L)

# choose center value
tm[, y := dplyr::case_when(
  tmrca_stat == "mean" ~ tmrca_mean,
  tmrca_stat == "mode" ~ ifelse(is.na(tmrca_mode), tmrca_mean, tmrca_mode),
  tmrca_stat == "lo"   ~ ci_low,
  tmrca_stat == "hi"   ~ ci_high,
  TRUE ~ tmrca_mean
)]

# ---------- WINDOWING ON SELECTED REGION ----------
winDT <- split_to_windows_region(tm, win, x_start, x_end)
stopifnot(nrow(winDT) > 0L)
winDT <- winDT %>% mutate(mid_bp = (w_start + w_end)/2, mid_mb = mid_bp/1e6)


# ---------- FLAGS (young/mid/old) ----------
winDT <- winDT %>%
  group_by(chrom) %>%
  mutate(
    y_rank = percent_rank(y_wmean),
    flag   = case_when(
      !is.na(y_rank) & y_rank <= bottom_q ~ "young",
      !is.na(y_rank) & y_rank >= top_q    ~ "old",
      TRUE                                ~ "mid"
    )
  ) %>% ungroup()

# ---------- GENES (lines ON TOP + labels ABOVE lines) ----------
genes <- read_annotation_1ab(ann_xlsx)
g <- genes %>%
  filter(chrom == contig, gene_end >= x_start, gene_start <= x_end) %>%
  mutate(
    xmin   = pmax(gene_start, x_start)/1e6,
    xmax   = pmin(gene_end,   x_end)/1e6,
    xmid   = (xmin + xmax)/2,
    gene_id = as.character(gene_id)
  )
if (is.finite(max_gene_rows)) g <- head(g, max_gene_rows)

# compute top y positions (respecting log scale)
y_top_data <- max(winDT$y_wmean + winDT$ciw_wmean/2, na.rm = TRUE)
y_top_line <- y_top_data * 1.05
y_top_text <- y_top_data * 1.15

# add constants to gene df for plotting
if (nrow(g) > 0) {
  g$y_line <- y_top_line
  g$y_text <- y_top_text
}

# ---------- PLOT ----------
p <- ggplot(winDT, aes(x = mid_mb, y = y_wmean)) +
  geom_ribbon(aes(ymin = y_wmean - ciw_wmean/2,
                  ymax = y_wmean + ciw_wmean/2),
              alpha = 0.15) +
  geom_line(linewidth = 0.5) +
  geom_point(
    data = dplyr::filter(winDT, flag %in% c("young","old")),
    aes(color = flag),
    size = 1.0,
    inherit.aes = TRUE
  ) +
  scale_color_manual(values = c(young = "#2C7BB6", old = "#D7191C")) +
  labs(
    x = sprintf("%s position (Mb)", contig),
    y = paste0("TMRCA (", tmrca_stat, if (logY) ", log10" else "", ")"),
    title = sprintf("%s:%s-%s  (win=%s kb)", contig,
                    format(x_start, big.mark=","), format(x_end, big.mark=","),
                    win/1000)
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# add gene lines & labels at the TOP
if (nrow(g) > 0) {
  p <- p +
    geom_segment(data = g,
                 aes(x = xmin, xend = xmax, y = y_line, yend = y_line),
                 linewidth = 1, inherit.aes = FALSE) +
    { if (label_genes)
      geom_text(data = g,
                aes(x = xmid, y = y_text, label = gene_id),
                size = 3, vjust = 0, inherit.aes = FALSE)
      else NULL } +
    expand_limits(y = y_top_text * 1.02)
}

if (logY) {
  p <- p + scale_y_continuous(trans = "log10", labels = label_number())
}
p <- p + scale_x_continuous(labels = label_number(accuracy = 1))

print(p)
ggsave(out_png, p, width = 12, height = 4, dpi = 300)
#ggsave(out_pdf, p, width = 12, height = 4)
