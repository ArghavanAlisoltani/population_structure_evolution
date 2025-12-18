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
  tree_i <- 0L
  meta <- vector("list", 10000)

  repeat{
    lines <- readLines(con, n=2000)
    if (length(lines) == 0) break
    for (ln in lines){
      if (startsWith(ln, "NAMES\t")) {
        f <- strsplit(ln, "\t", fixed=TRUE)[[1]]
        sample_names <- f[-1]
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
  list(sample_names=sample_names, meta=meta)
}

pick_hit_plus_flanks <- function(meta, position, flank=1){
  hit <- meta[start <= position & end >= position, tree_index]
  if (length(hit) == 0) stop("Position not covered by any TREE segment: ", position)
  hit <- hit[1]
  idx <- seq(hit - flank, hit + flank)
  idx <- idx[idx >= 1 & idx <= nrow(meta)]
  list(hit=hit, idx=idx)
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

# ------------------ Axis "break" via tail compression ------------------
apply_xbreak_compress <- function(p, mode=c("none","auto","manual"),
                                  break_at=NULL,
                                  shrink=0.03,
                                  trigger_ratio=6,
                                  q=0.90,
                                  draw_slash=TRUE){

  mode <- match.arg(mode)
  if (mode == "none") return(p)

  x0 <- p$data$x
  if (!is.numeric(x0) || length(x0)==0) return(p)

  # decide breakpoint
  if (mode == "manual") {
    if (is.null(break_at) || is.na(break_at)) return(p)
    bx <- break_at
  } else {
    # AUTO: only do it if max is a big outlier vs typical scale
    # use median of internal nodes as "typical"
    x_typ <- median(x0[!p$data$isTip], na.rm=TRUE)
    if (!is.finite(x_typ) || x_typ <= 0) x_typ <- median(x0, na.rm=TRUE)

    if (!is.finite(x_typ) || x_typ <= 0) return(p)

    if (max(x0, na.rm=TRUE) / x_typ < trigger_ratio) {
      return(p)  # no extreme tail; keep normal axis
    }
    bx <- as.numeric(quantile(x0, probs=q, na.rm=TRUE))
  }

  if (!is.finite(bx) || bx <= 0 || bx >= max(x0, na.rm=TRUE)) return(p)

  # piecewise compress tail
  tx <- function(x){
    ifelse(x <= bx, x, bx + (x - bx) * shrink)
  }

  p$data$x <- tx(x0)

  # keep axis labels meaningful: label original values at mapped positions
  max_x0 <- max(x0, na.rm=TRUE)
  breaks_orig <- unique(c(0, round(bx), round(max_x0)))
  breaks_disp <- tx(breaks_orig)

  # add a little visual "slash" near the x-axis
  if (draw_slash) {
    y0 <- min(p$data$y, na.rm=TRUE) - 0.8
    xr <- diff(range(p$data$x, na.rm=TRUE))
    dx <- 0.012 * xr
    dy <- 0.35
    xslash <- tx(bx)

    p <- p +
      annotate("segment", x=xslash - dx, xend=xslash + dx, y=y0 - dy, yend=y0 + dy, linewidth=0.5) +
      annotate("segment", x=xslash - dx*0.4, xend=xslash + dx*1.6, y=y0 - dy, yend=y0 + dy, linewidth=0.5) +
      coord_cartesian(clip="off") +
      theme(plot.margin = margin(5, 5, 18, 5))
  }

  p <- p + scale_x_continuous(breaks=breaks_disp, labels=breaks_orig)

  p
}

# ------------------ ID matching helper ------------------
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
tip_size_mid <- as.numeric(get_arg(args, "--tip_size_mid", "1.1"))

# x-break settings
xbreak_mode <- tolower(get_arg(args, "--xbreak", "none"))   # none|auto|manual
xbreak_at   <- as.numeric(get_arg(args, "--xbreak_at", NA)) # for manual
xbreak_shrink <- as.numeric(get_arg(args, "--xbreak_shrink", "0.03"))
xbreak_ratio  <- as.numeric(get_arg(args, "--xbreak_trigger", "6"))
xbreak_q      <- as.numeric(get_arg(args, "--xbreak_q", "0.90"))

# overall figure size
fig_w <- as.numeric(get_arg(args, "--width", "18"))
fig_h <- as.numeric(get_arg(args, "--height", "7"))
dpi   <- as.integer(get_arg(args, "--dpi", "300"))
format_out <- get_arg(args, "--format", "png")  # pdf|png|both

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

if (is.null(smc_file) || is.null(sites_file) || is.null(pheno_file) || is.na(position)) {
  cat("Usage:\n",
      "  Rscript plot_flanks_middle_sidebars_xbreak.R \\\n",
      "    --smc=FILE.smc --sites=FILE.sites --pheno=PHENO.txt \\\n",
      "    --position=983057685 --flank=1 --trait=C13 \\\n",
      "    --xbreak=auto --outdir=out --outprefix=scaffold4_pos983057685\n", sep="")
  quit(status=1)
}

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

allele_positions <- if (is.null(pos_string) || pos_string == "") position else as.integer(strsplit(pos_string, ",")[[1]])

# --- scan SMC and read trees ---
scan <- scan_smc_tree_intervals(smc_file)
meta <- scan$meta
sel <- pick_hit_plus_flanks(meta, position=position, flank=flank)
hit_idx <- sel$hit
keep <- sel$idx
meta_sub <- meta[tree_index %in% keep]
setorder(meta_sub, start)

trees <- read_trees_by_index(
  smc_file=smc_file,
  sample_names=scan$sample_names,
  n_total=nrow(meta),
  keep_indices=meta_sub$tree_index,
  ladderize_tree=ladderize_tr
)

# choose immediate neighbors around hit within meta_sub order
ord <- meta_sub$tree_index
hit_pos_in_ord <- which(ord == hit_idx)
left_idx  <- if (hit_pos_in_ord > 1) ord[hit_pos_in_ord - 1] else NA_integer_
right_idx <- if (hit_pos_in_ord < length(ord)) ord[hit_pos_in_ord + 1] else NA_integer_

# helper: build tree panel and optionally apply xbreak compression
make_tree_plot <- function(tr, st, en, show_labels=FALSE, tip_size=1.2){
  branch_mode <- if (tree_scale == "cladogram") "none" else "branch.length"
  seg_len <- en - st + 1L
  ttl <- sprintf("TREE %d–%d (len=%s)", st, en, format(seg_len, big.mark=","))

  p <- ggtree(tr, branch.length=branch_mode, size=0.35) +
    ggtitle(ttl) +
    theme_tree2() +
    theme(plot.title = element_text(size=9, hjust=0.5),
          plot.margin = margin(5, 5, 5, 5))

  if (show_labels) p <- p + geom_tiplab(size=tip_size)

  # Apply xbreak only if branch lengths are shown
  if (tree_scale != "cladogram") {
    p <- apply_xbreak_compress(
      p,
      mode = ifelse(xbreak_mode %in% c("none","auto","manual"), xbreak_mode, "none"),
      break_at = if (!is.na(xbreak_at)) xbreak_at else NULL,
      shrink = xbreak_shrink,
      trigger_ratio = xbreak_ratio,
      q = xbreak_q,
      draw_slash = TRUE
    )
  }

  p
}

# left tree
if (!is.na(left_idx)) {
  L <- trees[[as.character(left_idx)]]
  p_left <- make_tree_plot(L$tree, L$start, L$end, show_labels=FALSE)
} else {
  p_left <- ggplot() + theme_void() + ggtitle("No upstream tree")
}

# right tree
if (!is.na(right_idx)) {
  R <- trees[[as.character(right_idx)]]
  p_right <- make_tree_plot(R$tree, R$start, R$end, show_labels=FALSE)
} else {
  p_right <- ggplot() + theme_void() + ggtitle("No downstream tree")
}

# middle tree (hit)
M <- trees[[as.character(hit_idx)]]
p_mid_tree <- make_tree_plot(M$tree, M$start, M$end, show_labels=show_tip_labels_mid, tip_size=tip_size_mid)

# Extract tip y coords from the middle tree plot
tip_df <- as.data.table(p_mid_tree$data[p_mid_tree$data$isTip, c("label","y")])
setnames(tip_df, c("label","y"))
tip_df[, label_clean := clean_id(label)]
tip_df[, label_digits := extract_digits(label)]

# --- phenotype sidebars aligned to middle tips ---
ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
stopifnot("codg" %in% names(ph), "proc" %in% names(ph), "site" %in% names(ph))
if (!(trait_name %in% names(ph))) stop("Trait not found in phenotype file: ", trait_name)

ph[, codg := clean_id(codg)]
ph[, codg_digits := extract_digits(codg)]

diag_ph <- choose_key_mode(tip_df$label, ph$codg)
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

proc_map <- c("1"="Deer Mtn","2"="Inverness River","3"="Judy Creek","4"="Swann Hills","5"="Virginia Hills")
site_map <- c("1"="JUDY","2"="VIRG","3"="SWAN","4"="TIME")

anno[, proc_lab := factor(proc_map[as.character(proc)],
                          levels=c("Deer Mtn","Inverness River","Judy Creek","Swann Hills","Virginia Hills"))]
anno[, site_lab := factor(site_map[as.character(site)],
                          levels=c("JUDY","VIRG","SWAN","TIME"))]

# --- allele heatmap aligned to middle tips ---
sx <- read_sites(sites_file)
geno_wide <- extract_genotypes(sx$names, sx$pos, sx$allele_str, allele_positions)
geno_long <- melt(as.data.table(geno_wide), id.vars="label", variable.name="position", value.name="allele")
geno_long[, position := factor(position, levels=as.character(allele_positions))]
geno_long[, label_clean := clean_id(label)]
geno_long[, label_digits := extract_digits(label)]

diag_sites <- choose_key_mode(tip_df$label, sx$names)
if (diag_sites$use_digits) {
  geno_long[, key := label_digits]
  tip_keys <- tip_df[, .(key=label_digits, label, y)]
} else {
  geno_long[, key := label_clean]
  tip_keys <- tip_df[, .(key=label_clean, label, y)]
}
geno_key <- merge(geno_long[, .(key, position, allele)], tip_keys, by="key", all.x=TRUE, sort=FALSE)
geno_key <- geno_key[!is.na(label)]

yl <- range(tip_df$y, na.rm=TRUE)

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
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        plot.margin = margin(5, 2, 5, 2))

mid_block <- p_mid_tree + p_proc + p_site + p_trait + p_allele +
  plot_layout(widths=c(w_mid_tree, w_proc, w_site, w_trait, w_allele), guides="collect") &
  theme(legend.position="right")

top_row <- (p_left | mid_block | p_right) +
  plot_layout(widths=c(w_left, w_mid, w_right))

base <- file.path(outdir, outprefix)
tag  <- paste0(".pos", position, ".flank", flank, ".xbreak_", xbreak_mode)

if (format_out %in% c("pdf","both")) {
  ggsave(paste0(base, tag, ".pdf"), top_row, width=fig_w, height=fig_h)
}
if (format_out %in% c("png","both")) {
  ggsave(paste0(base, tag, ".png"), top_row, width=fig_w, height=fig_h, dpi=dpi)
}

cat("Wrote outputs under:", outdir, "\n")

