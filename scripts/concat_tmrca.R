library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(tools)
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/all_tmrca/all_tmrca/")

tmrcaAll<-data.frame(fread("All_tmrca_60_run.txt",header=F))
names(tmrcaAll)<-c("CHROM", "start", "end", "posterior_mean_TMRCA",
                "lower_2.5_percentile_TMRCA", "upper_97.5_percentile_TMRCA")
scaffods<-unique(tmrcaAll$CHROM)
summary(tmrcaAll$posterior_mean_TMRCA)
tmrcaAll_top<-tmrcaAll[tmrcaAll$posterior_mean_TMRCA>134021,]
fwrite(tmrcaAll_top,"top_tmrca_greater_than_134021.tsv", sep="\t")
scaffodstop<-unique(tmrcaAll_top$CHROM)

# =========================
# Violin plot of TMRCA (mean) per scaffold with vertical x-axis labels
# =========================

# install.packages(c("data.table","ggplot2","dplyr","scales"))
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

# ---- Option A: read multiple *.tmrca.txt files (ARGweaver format: chrom start end mean [mode] lo hi)
files <- Sys.glob("*.tmrca.txt")
stopifnot(length(files) > 0)
tm_list <- lapply(files, function(f) {
  dt <- fread(f, header = FALSE, integer64 = "numeric")
  if (ncol(dt) == 7) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_mode","tmrca_lo","tmrca_hi"))
  } else if (ncol(dt) == 6) {
    setnames(dt, c("chrom","start","end","tmrca_mean","tmrca_lo","tmrca_hi"))
  } else stop(sprintf("Unexpected column count in %s", f))
  dt[, .(chrom, tmrca_mean)]
})
tm <- rbindlist(tm_list, use.names = TRUE, fill = TRUE)

# ---- Option B: if you already have a combined table 'tm' with columns chrom, tmrca_mean, comment out Option A above ----
# tm <- fread("your_combined_tmrca.tsv")  # must contain 'chrom' and 'tmrca_mean' columns

# Reorder scaffolds by median TMRCA (optional for nicer ordering)
tm <- tm %>%
  mutate(scaffold = reorder(as.factor(chrom), tmrca_mean, FUN = median, na.rm = TRUE))

# Violin plot (log10 y), with mean (white dot) and median (black dot); x labels vertical
p <- ggplot(tm, aes(x = scaffold, y = tmrca_mean)) +
  geom_violin(fill = "grey85", color = "grey40", trim = FALSE) +
  stat_summary(fun = median, geom = "point", size = 1.4, color = "black") +
  stat_summary(fun = mean,   geom = "point", size = 2.0, color = "red", stroke = 0.6) +
  scale_y_continuous(trans = "log10", labels = label_number()) +
  labs(x = "Scaffold", y = "TMRCA (mean, log10 generations)", title = "TMRCA distribution by scaffold") +
  theme_bw(base_size = 17) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p)

# Save (adjust height if many scaffolds)
ggsave("tmrca_violin_per_scaffold.png", p, width = 12, height = 6, dpi = 300)
ggsave("tmrca_violin_per_scaffold.pdf", p, width = 12, height = 6)



# ================================
# Map each TMRCA segment to genes:
#  - list ALL overlapping genes
#  - if none overlap, list NEAREST gene(s) (ties) and distance
# Uses annotation columns: contig_1ab, start1ab_pos, end1ab_pos (and NR.ID if present)
# ================================

# install.packages(c("data.table","readxl","dplyr","GenomicRanges","stringr","tools"))
library(data.table)
library(readxl)
library(dplyr)
library(GenomicRanges)
library(stringr)
library(tools)

# -------- USER INPUTS --------
tmrca_path <- "All_tmrca_60_run.txt"  # your TMRCA table (all scaffolds)
ann_xlsx   <- "~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/Aria_curated_annotation_1ab.xlsx"   # annotation workbook

tmrca_zero_based <- TRUE       # set TRUE if tmrca start is 0-based half-open (ARGweaver default)
nearby_bp        <- 50000      # report nearest genes (ties) at this minimal distance; keep even if > nearby_bp
out_csv          <- "tmrca_to_genes_overlap_nearest.csv"

# -------- READ TMRCA (robust header/columns) --------
read_tmrca_generic <- function(path, zero_based = TRUE) {
  # fread auto-detects sep; if your file is truly TSV, it's fine
  dt <- fread(path, integer64 = "numeric")
  # Try to detect columns for chrom/start/end
  cn <- tolower(names(dt))
  guess_col <- function(cands, default) {
    i <- which(cn %in% cands)
    if (length(i)) names(dt)[i[1]] else default
  }
  col_chrom <- guess_col(c("chrom","chr","scaffold","contig"), "V1")
  col_start <- guess_col(c("start","pos_start","begin"),      "V2")
  col_end   <- guess_col(c("end","pos_end","stop"),           "V3")
  
  # Coerce and standardize
  dt[, chrom := as.character(get(col_chrom))]
  dt[, start := as.numeric(get(col_start))]
  dt[, end   := as.numeric(get(col_end))]
  
  if (zero_based) dt[, start := start + 1]  # 0-based -> 1-based inclusive
  # Keep original extra columns if present
  setcolorder(dt, c("chrom","start","end", setdiff(names(dt), c("chrom","start","end"))))
  dt[]
}

tm <- read_tmrca_generic(tmrca_path, zero_based = tmrca_zero_based)
stopifnot(all(c("chrom","start","end") %in% names(tm)))

# -------- READ ANNOTATION (required columns) --------
read_annotation_1ab <- function(xlsx, sheet = NULL, id_col = "NR.ID") {
  if (is.null(sheet)) sheet <- excel_sheets(xlsx)[1]
  ann <- read_excel(xlsx, sheet = sheet, guess_max = 200000)
  req <- c("contig_1ab","start1ab_pos","end1ab_pos")
  missing_req <- setdiff(req, names(ann))
  if (length(missing_req)) stop("Missing required columns in annotation: ", paste(missing_req, collapse = ", "))
  
  if (!id_col %in% names(ann)) {
    ann[[id_col]] <- paste0("GENE_", seq_len(nrow(ann)))
  }
  
  out <- ann |>
    transmute(
      chrom      = as.character(.data$contig_1ab),
      gene_start = as.numeric(.data$start1ab_pos),
      gene_end   = as.numeric(.data$end1ab_pos),
      gene_id    = as.character(.data[[id_col]])
    ) |>
    arrange(chrom, gene_start, gene_end)
  as.data.table(out)
}

genes <- read_annotation_1ab(ann_xlsx, id_col = "NR.ID")
stopifnot(all(c("chrom","gene_start","gene_end","gene_id") %in% names(genes)))

# -------- BUILD GRanges --------
tm[, tm_id := .I]  # keep original row id
tm_gr <- GRanges(seqnames = tm$chrom, ranges = IRanges(start = tm$start, end = tm$end), tm_id = tm$tm_id)

genes_gr <- GRanges(seqnames = genes$chrom, ranges = IRanges(start = genes$gene_start, end = genes$gene_end),
                    gene_id = genes$gene_id, gene_start = genes$gene_start, gene_end = genes$gene_end)

# -------- FIND OVERLAPS --------
ov <- findOverlaps(tm_gr, genes_gr, ignore.strand = TRUE)  # Hits(query, subject)

# Overlap table (may be empty)
ov_tab <- if (length(ov)) {
  data.table(
    tm_id    = mcols(tm_gr)$tm_id[queryHits(ov)],
    gene_id  = mcols(genes_gr)$gene_id[subjectHits(ov)],
    gene_chr = as.character(seqnames(genes_gr))[subjectHits(ov)],
    gene_start = mcols(genes_gr)$gene_start[subjectHits(ov)],
    gene_end   = mcols(genes_gr)$gene_end[subjectHits(ov)]
  )
} else data.table(tm_id = integer(), gene_id = character(), gene_chr=character(), gene_start = numeric(), gene_end = numeric())

# Collapse overlapping genes per TMRCA row
ov_collapsed <- ov_tab[, .(
  overlap_genes = paste(unique(gene_id), collapse = ";"),
  overlap_gene_n = uniqueN(gene_id)
), by = tm_id]

# -------- NEAREST genes for non-overlapping rows (ties kept) --------
no_ov_ids <- setdiff(tm$tm_id, ov_collapsed$tm_id)
nearest_tab <- data.table()

if (length(no_ov_ids)) {
  q_gr <- tm_gr[match(no_ov_ids, mcols(tm_gr)$tm_id)]
  # Get all nearest ties
  nh <- nearest(q_gr, genes_gr, select = "all", ignore.strand = TRUE)  # Hits with all nearest ties
  if (length(nh)) {
    q_idx <- queryHits(nh)
    s_idx <- subjectHits(nh)
    # distances
    dist_bp <- distance(q_gr[q_idx], genes_gr[s_idx], ignore.strand = TRUE)
    
    nearest_tab <- data.table(
      tm_id    = mcols(q_gr)$tm_id[q_idx],
      gene_id  = mcols(genes_gr)$gene_id[s_idx],
      gene_chr = as.character(seqnames(genes_gr))[s_idx],
      gene_start = mcols(genes_gr)$gene_start[s_idx],
      gene_end   = mcols(genes_gr)$gene_end[s_idx],
      nearest_dist_bp = as.numeric(dist_bp)
    )
    
    # collapse per tm_id: all tied nearest genes
    # also flag those within nearby_bp (if the minimal distance <= nearby_bp)
    nearest_collapsed <- nearest_tab[, .(
      nearest_genes        = paste(unique(gene_id), collapse = ";"),
      nearest_min_dist_bp  = min(nearest_dist_bp, na.rm = TRUE),
      nearest_genes_within = paste(unique(gene_id[nearest_dist_bp <= nearby_bp]), collapse = ";")
    ), by = tm_id]
    
    # ensure empty strings -> NA for neatness
    nearest_collapsed[nearest_genes_within == "", nearest_genes_within := NA_character_]
  } else {
    nearest_collapsed <- data.table(tm_id = no_ov_ids,
                                    nearest_genes = NA_character_,
                                    nearest_min_dist_bp = NA_real_,
                                    nearest_genes_within = NA_character_)
  }
} else {
  nearest_collapsed <- data.table(tm_id = integer(),
                                  nearest_genes = character(),
                                  nearest_min_dist_bp = numeric(),
                                  nearest_genes_within = character())
}

# -------- JOIN everything back to TMRCA rows --------
res <- as.data.table(tm)[, .(tm_id, chrom, start, end, everything = .SD), .SDcols = setdiff(names(tm), c("tm_id","chrom","start","end"))]

res <- merge(res, ov_collapsed, by = "tm_id", all.x = TRUE)
res <- merge(res, nearest_collapsed, by = "tm_id", all.x = TRUE)

# For convenience, also report a unified 'genes_reported' column:
res[, genes_reported := fifelse(!is.na(overlap_genes) & overlap_genes != "",
                                overlap_genes,
                                fifelse(!is.na(nearest_genes) & nearest_genes != "", nearest_genes, NA_character_))]

# Order columns nicely
ord <- c("tm_id","chrom","start","end",
         setdiff(names(res), c("tm_id","chrom","start","end")))
setcolorder(res, ord)

# -------- WRITE OUTPUT --------
fwrite(res, out_csv)

# --------- OPTIONAL: quick peek ----------
# print(res[1:20])
# table(is.na(res$overlap_genes))  # how many had overlaps vs. not

mytbl<-data.frame(table(res$chrom, res$overlap_genes))

  
  
