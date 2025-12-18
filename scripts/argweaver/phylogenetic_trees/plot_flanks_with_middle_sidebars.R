#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(data.table)
  library(patchwork)   # install.packages("patchwork")
})

# ------------------ CLI helpers ------------------
get_arg <- function(args, flag, default=NULL){
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if(length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}
as_logical <- function(x, default=FALSE){
  if (is.null(x)) return(default)
  tolower(x) %in% c("true","t","1","yes","y")
}

clean_id <- function(x){
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- gsub("^X", "", x)
  x
}
extract_digits <- function(x){
  x <- clean_id(x)
  has_digit <- grepl("\\d", x)
  y <- rep(NA_character_, length(x))
  y[has_digit] <- sub(".*?(\\d+).*", "\\1", x[has_digit])
  y
}
mode1 <- function(x){
  x <- x[!is.na(x)]
  if(length(x)==0) return(NA_character_)
  tb <- sort(table(x), decreasing=TRUE)
  names(tb)[1]
}

strip_nhx <- function(newick){
  gsub("\\[&&NHX:[^\\]]*\\]", "", newick, perl=TRUE)
}

map_tip_labels <- function(tr, sample_names){
  idx <- suppressWarnings(as.integer(tr$tip.label))
  if (any(is.na(idx))) stop("Tree tips are not numeric indices. Example tip: ", tr$tip.label[which(is.na(idx))[1]])
  if (any(idx < 0) || any(idx >= length(sample_names))) {
    bad <- idx[which(idx < 0 | idx >= length(sample_names))[1]]
    stop("Tip index out of range 0..", length(sample_names)-1, ". Example bad index: ", bad)
  }
  tr$tip.label <- sample_names[idx + 1]
  tr
}

# ------------------ .sites readers ------------------
read_sites <- function(sites_file) {
  lines <- readLines(sites_file)
  name_line <- lines[grepl("^#NAMES\\t", lines)][1]
  if (is.na(name_line)) stop("No #NAMES line found. Is this an ARGweaver .sites file?")
  nm <- strsplit(name_line, "\t", fixed=TRUE)[[1]][-1]
  n <- length(nm)

  site_lines <- lines[!grepl("^#", lines)]
  site_lines <- site_lines[nchar(site_lines) > 0]
  pos <- as.integer(sub("\\t.*$", "", site_lines))
  allele_str <- sub("^[0-9]+\\t", "", site_lines)

  bad <- which(nchar(allele_str) != n)
  if (length(bad) > 0) stop("Allele string length != n in .sites. First bad line: ", bad[1])

  list(names=nm, pos=pos, allele_str=allele_str)
}

extract_genotypes <- function(names, pos, allele_str, positions_numeric) {
  n <- length(names)
  idx <- match(positions_numeric, pos)
  if (anyNA(idx)) stop("Requested positions not found in .sites: ", paste(positions_numeric[is.na(idx)], collapse=", "))

  geno <- data.table(label = names)
  for (k in seq_along(idx)) {
    s <- allele_str[idx[k]]
    geno[[as.character(positions_numeric[k])]] <- substring(s, seq_len(n), seq_len(n))
  }
  geno
}

# ------------------ SMC scan + selected trees ------------------
scan_smc_tree_intervals <- function(smc_file){
  con <- file(smc_file, "r")
  on.exit(close(con), add=TRUE)

  sample_names <- NULL
  region <- list(chr=NA_character_, start=NA_integer_, end=NA_integer_)
  tree_i <- 0L
  meta <- vector("list", 10000)

  repeat{
    lines <- readLines(con, n=2000)
    if (length(lines) == 0) break
    for (ln in lines){
      if (startsWith(ln, "NAMES\t")) {
        f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
        sample_names <- f[-1]
      } else if (startsWith(ln, "REGION\t")) {
        f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
        region$chr   <- f[2]
        region$start <- as.integer(f[3])
        region$end   <- as.integer(f[4])
      } else if (startsWith(ln, "TREE\t")) {
        f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
        tree_i <- tree_i + 1L
        st <- as.integer(f[2]); en <- as.integer(f[3])
        meta[[tree_i]] <- list(tree_index=tree_i, start=st, end=en, length=(en - st + 1L))
      }
    }
  }
  if (is.null(sample_names)) stop("No NAMES line found in .smc")
  meta <- rbindlist(meta[seq_len(tree_i)])
  list(sample_names=sample_names, region=region, meta=meta)
}

pick_hit_plus_flanks <- function(meta, position, flank=1){
  hit <- meta[start <= position & end >= position, tree_index]
  if (length(hit) == 0) stop("Position not covered by any TREE segment: ", position)
  hit <- hit[1]
  idx <- seq(hit - flank, hit + flank)
  idx <- idx[idx >= 1 & idx <= nrow(meta)]
  idx
}

read_trees_by_index <- function(smc_file, sample_names, n_total, keep_indices, ladderize_tree=TRUE){
  keep_set <- rep(FALSE, n_total)
  keep_set[keep_indices] <- TRUE

  trees <- vector("list", length(keep_indices))
  names(trees) <- as.character(keep_indices)

  con <- file(smc_file, "r")
  on.exit(close(con), add=TRUE)

  tree_i <- 0L
  repeat{
    lines <- readLines(con, n=2000)
    if (length(lines) == 0) break
    for (ln in lines){
      if (!startsWith(ln, "TREE\t")) next
      tree_i <- tree_i + 1L
      if (!keep_set[tree_i]) next

      f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
      st <- as.integer(f[2]); en <- as.integer(f[3])
      newick <- strip_nhx(f[4])

      tr <- read.tree(text=newick)
      tr <- map_tip_labels(tr, sample_names)
      if (ladderize_tree) tr <- ladderize(tr)

      trees[[as.character(tree_i)]] <- list(tree=tr, start=st, end=en)
    }
  }

  missing <- keep_indices[!sapply(as.character(keep_indices), function(k) !is.null(trees[[k]]))]
  if (length(missing) > 0) stop("Did not find these TREE indices in file: ", paste(missing, collapse=", "))
  trees
}

# ------------------ Genome segment track ------------------
plot_segment_track <- function(meta_sub, focal_pos, x_unit=c("bp","Mb"), label_bp=TRUE){
  x_unit <- match.arg(x_unit)
  df <- copy(meta_sub)
  setorder(df, start)
  df[, is_hit := FALSE]
  df[tree_index == df$tree_index[ceiling(.N/2)], is_hit := TRUE]  # placeholder; we overwrite below

  bps <- df$start[-1]  # breakpoints

  if (x_unit == "Mb") {
    df[, `:=`(x1=start/1e6, x2=end/1e6)]
    bps2 <- bps/1e6
    focal2 <- focal_pos/1e6
    xlab <- "Genomic position (Mb)"
  } else {
    df[, `:=`(x1=start, x2=end)]
    bps2 <- bps
    focal2 <- focal_pos
    xlab <- "Genomic position (bp)"
  }

  p <- ggplot(df) +
    geom_rect(aes(xmin=x1, xmax=x2, ymin=0, ymax=1),
              fill="grey85", color="grey30", linewidth=0.3) +
    geom_vline(xintercept=bps2, color="red", linewidth=0.5) +
    geom_vline(xintercept=focal2, color="blue", linewidth=0.6) +
    scale_y_continuous(expand=c(0,0), breaks=NULL) +
    labs(x=xlab, y=NULL) +
    theme_minimal(base_size=10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(2, 10, 2, 10)
    )

  if (label_bp && length(bps2) > 0) {
    lab <- data.table(x=bps2)
    p <- p + geom_text(data=lab, aes(x=x, y=1.10, label="BP"),
                       size=3, color="red", vjust=0)
  }
  p
}

# ------------------ Sidebars aligned to middle tree ------------------
choose_key_mode <- function(tip_labels, ref_ids){
  tips_clean <- clean_id(tip_labels)
  tips_digits <- extract_digits(tip_labels)

  ref_clean <- clean_id(ref_ids)
  ref_digits <- extract_digits(ref_ids)

  exact_hits <- sum(tips_clean %in% ref_clean)
  digit_hits <- sum(!is.na(tips_digits) & (tips_digits %in% ref_digits))

  list(use_digits = (digit_hits > exact_hits),
       exact_hits = exact_hits,
       digit_hits = digit_hits)
}

# ------------------ Main ------------------
args <- commandArgs(trailingOnly=TRUE)

smc_file   <- get_arg(args, "--smc")
sites_file <- get_arg(args, "--sites")
pheno_file <- get_arg(args, "--pheno")

position   <- as.integer(get_arg(args, "--position", NA))
flank      <- as.integer(get_arg(args, "--flank", "1"))

trait_name <- get_arg(args, "--trait", "C13")
pos_string <- get_arg(args, "--allele_positions", NULL)  # default = --position

outdir     <- get_arg(args, "--outdir", "trees_out")
outprefix  <- get_arg(args, "--outprefix", "pos_flanks_mid_sidebars")

tree_scale <- tolower(get_arg(args, "--tree_scale", "distance"))  # distance|cladogram
ladderize_tr <- as_logical(get_arg(args, "--ladderize", "true"), TRUE)

show_tip_labels_mid <- as_logical(get_arg(args, "--tip_labels_mid", "false"), FALSE)
tip_size_mid <- as.numeric(get_arg(args, "--tip_size_mid", "1.3"))

# overall figure size
fig_w <- as.numeric(get_arg(args, "--width", "18"))
fig_h <- as.numeric(get_arg(args, "--height", "8"))
dpi   <- as.integer(get_arg(args, "--dpi", "300"))
format_out <- get_arg(args, "--format", "pdf")  # pdf|png|both

# widths: left tree | middle block | right tree
w_left  <- as.numeric(get_arg(args, "--w_left",  "4"))
w_mid   <- as.numeric(get_arg(args, "--w_mid",   "10"))
w_right <- as.numeric(get_arg(args, "--w_right", "4"))

# inside middle block: tree | proc | site | trait | allele
w_mid_tree  <- as.numeric(get_arg(args, "--w_mid_tree",  "5"))
w_proc      <- as.numeric(get_arg(args, "--w_proc",      "0.35"))
w_site      <- as.numeric(get_arg(args, "--w_site",      "0.35"))
w_trait     <- as.numeric(get_arg(args, "--w_trait",     "0.45"))
w_allele    <- as.numeric(get_arg(args, "--w_allele",    "0.9"))

x_unit   <- get_arg(args, "--x_unit", "bp")   # bp|Mb
bp_label <- as_logical(get_arg(args, "--bp_labels", "true"), TRUE)

if (is.null(smc_file) || is.null(sites_file) || is.null(pheno_file) || is.na(position)) {
  cat("Usage:\n",
      "  Rscript plot_flanks_with_middle_sidebars.R \\\n",
      "    --smc=FILE.smc --sites=FILE.sites --pheno=PHENO.txt \\\n",
      "    --position=983057685 --flank=1 --trait=C13 \\\n",
      "    --outdir=out --outprefix=scaffold4_pos983057685 --format=pdf\n", sep="")
  quit(status=1)
}

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

allele_positions <- if (is.null(pos_string) || pos_string == "") position else as.integer(strsplit(pos_string, ",")[[1]])

# --- scan SMC and read trees ---
scan <- scan_smc_tree_intervals(smc_file)
meta <- scan$meta
keep <- pick_hit_plus_flanks(meta, position=position, flank=flank)
meta_sub <- meta[tree_index %in% keep]
setorder(meta_sub, start)

# identify the hit tree index (middle of keep list by definition is hit, but safest by position)
hit_idx <- meta[start <= position & end >= position, tree_index][1]

trees <- read_trees_by_index(
  smc_file=smc_file,
  sample_names=scan$sample_names,
  n_total=nrow(meta),
  keep_indices=meta_sub$tree_index,
  ladderize_tree=ladderize_tr
)

# build 3 trees: left, middle(hit), right (if flank>1 you still get many; we show only immediate neighbors left/right)
# For visualization: choose left = segment just before hit, right = just after hit (if available)
ord <- meta_sub$tree_index
hit_pos_in_ord <- which(ord == hit_idx)
left_idx  <- if (hit_pos_in_ord > 1) ord[hit_pos_in_ord - 1] else NA_integer_
right_idx <- if (hit_pos_in_ord < length(ord)) ord[hit_pos_in_ord + 1] else NA_integer_

# Tree plotting helper
make_tree_plot <- function(tr, st, en, title_prefix="", show_labels=FALSE, tip_size=1.2){
  branch_mode <- if (tree_scale == "cladogram") "none" else "branch.length"
  seg_len <- en - st + 1L
  ttl <- sprintf("%sTREE #%s\n%d–%d (len=%s)", title_prefix, title_prefix, st, en, format(seg_len, big.mark=","))

  p <- ggtree(tr, branch.length=branch_mode, size=0.35) +
    ggtitle(sprintf("TREE %d–%d (len=%s)", st, en, format(seg_len, big.mark=","))) +
    theme_tree2() +
    theme(
      plot.title = element_text(size=9, hjust=0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  if (show_labels) p <- p + geom_tiplab(size=tip_size)
  p
}

# Left tree
if (!is.na(left_idx)) {
  L <- trees[[as.character(left_idx)]]
  p_left <- make_tree_plot(L$tree, L$start, L$end, show_labels=FALSE)
} else {
  p_left <- ggplot() + theme_void() + ggtitle("No upstream tree")
}

# Right tree
if (!is.na(right_idx)) {
  R <- trees[[as.character(right_idx)]]
  p_right <- make_tree_plot(R$tree, R$start, R$end, show_labels=FALSE)
} else {
  p_right <- ggplot() + theme_void() + ggtitle("No downstream tree")
}

# Middle tree (hit)
M <- trees[[as.character(hit_idx)]]
p_mid_tree <- make_tree_plot(M$tree, M$start, M$end, show_labels=show_tip_labels_mid, tip_size=tip_size_mid)

# Extract tip y coords from the middle tree plot
tip_df <- as.data.table(p_mid_tree$data[p_mid_tree$data$isTip, c("label","y")])
setnames(tip_df, c("label","y"))
tip_df[, label_clean := clean_id(label)]
tip_df[, label_digits := extract_digits(label)]

# --- phenotype sidebars (aligned to middle tips) ---
ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
stopifnot("codg" %in% names(ph), "proc" %in% names(ph), "site" %in% names(ph))
if (!(trait_name %in% names(ph))) stop("Trait not found in phenotype file: ", trait_name)

ph[, codg := clean_id(codg)]
ph[, codg_digits := extract_digits(codg)]

diag_ph <- choose_key_mode(tip_df$label, ph$codg)
cat("Phenotype ID matching:\n",
    "  exact hits:", diag_ph$exact_hits, "\n",
    "  digit hits:", diag_ph$digit_hits, "\n",
    "  using:", if (diag_ph$use_digits) "DIGITS" else "EXACT", "\n")

if (diag_ph$use_digits) {
  ph_agg <- ph[, .(
    proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
    site  = suppressWarnings(as.integer(mode1(as.character(site)))),
    trait = mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE)
  ), by=.(key = codg_digits)]
  tip_df[, key := label_digits]
} else {
  ph_agg <- ph[, .(
    proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
    site  = suppressWarnings(as.integer(mode1(as.character(site)))),
    trait = mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE)
  ), by=.(key = codg)]
  tip_df[, key := label_clean]
}

anno <- merge(tip_df, ph_agg, by="key", all.x=TRUE, sort=FALSE)

# Maps requested by you
proc_map <- c("1"="Deer Mtn","2"="Inverness River","3"="Judy Creek","4"="Swann Hills","5"="Virginia Hills")
site_map <- c("1"="JUDY","2"="VIRG","3"="SWAN","4"="TIME")

anno[, proc_lab := factor(proc_map[as.character(proc)],
                          levels=c("Deer Mtn","Inverness River","Judy Creek","Swann Hills","Virginia Hills"))]
anno[, site_lab := factor(site_map[as.character(site)],
                          levels=c("JUDY","VIRG","SWAN","TIME"))]

# --- allele heatmap (aligned to middle tips) ---
sx <- read_sites(sites_file)
geno_wide <- extract_genotypes(sx$names, sx$pos, sx$allele_str, allele_positions)
geno_long <- melt(as.data.table(geno_wide), id.vars="label", variable.name="position", value.name="allele")
geno_long[, position := factor(position, levels=as.character(allele_positions))]
geno_long[, label_clean := clean_id(label)]
geno_long[, label_digits := extract_digits(label)]

diag_sites <- choose_key_mode(tip_df$label, sx$names)
cat("Sites ID matching:\n",
    "  exact hits:", diag_sites$exact_hits, "\n",
    "  digit hits:", diag_sites$digit_hits, "\n",
    "  using:", if (diag_sites$use_digits) "DIGITS" else "EXACT", "\n")

if (diag_sites$use_digits) {
  geno_long[, key := label_digits]
  tip_keys <- tip_df[, .(key=label_digits, label, y)]
} else {
  geno_long[, key := label_clean]
  tip_keys <- tip_df[, .(key=label_clean, label, y)]
}

geno_key <- merge(geno_long[, .(key, position, allele)], tip_keys, by="key", all.x=TRUE, sort=FALSE)
geno_key <- geno_key[!is.na(label)]

# common y limits
yl <- range(tip_df$y, na.rm=TRUE)

# Sidebar plotting helpers
tile_col <- function(df, fill_var, legend_name, scale_fn){
  ggplot(df, aes(x=1, y=y, fill=.data[[fill_var]])) +
    geom_tile(width=0.95, height=0.9) +
    scale_y_continuous(limits=yl, expand=c(0,0)) +
    scale_fn +
    theme_void() +
    theme(plot.margin = margin(5, 2, 5, 2), legend.title=element_text(size=9), legend.text=element_text(size=8)) +
    guides(fill=guide_legend(title=legend_name))
}

p_proc <- ggplot(anno, aes(x=1, y=y, fill=proc_lab)) +
  geom_tile(width=0.95, height=0.9) +
  scale_fill_brewer(palette="Set2", na.value="grey92", name="Provenances") +
  scale_y_continuous(limits=yl, expand=c(0,0)) +
  theme_void() + theme(plot.margin = margin(5, 2, 5, 2))

p_site <- ggplot(anno, aes(x=1, y=y, fill=site_lab)) +
  geom_tile(width=0.95, height=0.9) +
  scale_fill_brewer(palette="Set3", na.value="grey92", name="Site") +
  scale_y_continuous(limits=yl, expand=c(0,0)) +
  theme_void() + theme(plot.margin = margin(5, 2, 5, 2))

p_trait <- ggplot(anno, aes(x=1, y=y, fill=trait)) +
  geom_tile(width=0.95, height=0.9) +
  scale_fill_viridis_c(na.value="grey92", name=trait_name) +
  scale_y_continuous(limits=yl, expand=c(0,0)) +
  theme_void() + theme(plot.margin = margin(5, 2, 5, 2))

allele_cols <- c("A"="#1b9e77","C"="#7570b3","G"="#d95f02","T"="#e7298a","N"="#bdbdbd","-"="#252525")
p_allele <- ggplot(geno_key, aes(x=position, y=y, fill=allele)) +
  geom_tile(height=0.9, color="grey85", linewidth=0.25) +
  scale_fill_manual(values=allele_cols, na.value="grey92", name="Allele") +
  scale_y_continuous(limits=yl, expand=c(0,0)) +
  theme_void() +
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
    plot.margin = margin(5, 2, 5, 2)
  )

# Middle block: tree + sidebars (no overlap, fixed widths)
mid_block <- p_mid_tree + p_proc + p_site + p_trait + p_allele +
  plot_layout(widths=c(w_mid_tree, w_proc, w_site, w_trait, w_allele), guides="collect") &
  theme(legend.position="right")

# Segment track for ALL selected segments (hit ± flank)
track <- plot_segment_track(meta_sub, focal_pos=position, x_unit=x_unit, label_bp=bp_label)

# Top row: left | mid_block | right (fixed layout)
top_row <- (p_left | mid_block | p_right) +
  plot_layout(widths=c(w_left, w_mid, w_right))

final <- top_row / track + plot_layout(heights=c(4, 1))

base <- file.path(outdir, outprefix)
tag  <- paste0(".pos", position, ".flank", flank)

if (format_out %in% c("pdf","both")) {
  ggsave(paste0(base, tag, ".pdf"), final, width=fig_w, height=fig_h)
}
if (format_out %in% c("png","both")) {
  ggsave(paste0(base, tag, ".png"), final, width=fig_w, height=fig_h, dpi=dpi)
}

cat("Wrote outputs under:", outdir, "\n")

