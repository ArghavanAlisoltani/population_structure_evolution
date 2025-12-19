#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(data.table)
  library(patchwork)
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

# ------------------ ID helpers ------------------
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
choose_key_mode <- function(tip_labels, ref_ids){
  tips_clean  <- clean_id(tip_labels)
  tips_digits <- extract_digits(tip_labels)
  ref_clean   <- clean_id(ref_ids)
  ref_digits  <- extract_digits(ref_ids)

  exact_hits <- sum(tips_clean %in% ref_clean)
  digit_hits <- sum(!is.na(tips_digits) & (tips_digits %in% ref_digits))

  list(use_digits = (digit_hits > exact_hits),
       exact_hits = exact_hits,
       digit_hits = digit_hits)
}

# ------------------ parse NHX from .smc Newick ------------------
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

# ------------------ x-axis as "time before present" ------------------
apply_time_before_present_labels <- function(p,
                                            show_title=TRUE,
                                            title="Generations before present (ARGweaver units)"){

  xmax <- max(p$data$x, na.rm=TRUE)

  br <- pretty(c(0, xmax), n=5)
  br <- sort(unique(c(br, 0, xmax)))
  br <- br[br >= 0 & br <= xmax]

  p <- p + scale_x_continuous(
    breaks = br,
    labels = function(v) format(round(xmax - v), big.mark=","),
    expand = c(0, 0)
  )

  if (show_title) {
    p <- p + xlab(title)
  } else {
    p <- p + xlab(NULL) + theme(axis.title.x = element_blank())
  }

  p
}

# add whitespace around ggplot OR patchwork object
add_panel_margin <- function(obj, panel_gap){
  m <- margin(5, panel_gap, 5, panel_gap)
  if (inherits(obj, "ggplot")) {
    obj + theme(plot.margin = m)
  } else {
    obj & theme(plot.margin = m)
  }
}

# ------------------ Main ------------------
args <- commandArgs(trailingOnly=TRUE)

smc_file   <- get_arg(args, "--smc")
sites_file <- get_arg(args, "--sites")
pheno_file <- get_arg(args, "--pheno")

position   <- as.integer(get_arg(args, "--position", NA))
flank      <- as.integer(get_arg(args, "--flank", "1"))

trait_name <- get_arg(args, "--trait", "C13")
pos_string <- get_arg(args, "--allele_positions", NULL)

outdir     <- get_arg(args, "--outdir", "trees_out")
outprefix  <- get_arg(args, "--outprefix", "pos_flanks_mid_sidebars")

tree_scale <- tolower(get_arg(args, "--tree_scale", "distance"))  # distance|cladogram
ladderize_tr <- as_logical(get_arg(args, "--ladderize", "true"), TRUE)

xlab_mode <- tolower(get_arg(args, "--xlab_mode", "time_before_present"))
if (!xlab_mode %in% c("distance","time_before_present")) xlab_mode <- "time_before_present"

show_tip_labels_hit <- as_logical(get_arg(args, "--tip_labels_hit", "true"), TRUE)
tip_size_hit <- as.numeric(get_arg(args, "--tip_size_hit", "0.9"))

show_tip_labels_flanks <- as_logical(get_arg(args, "--tip_labels_flanks", "false"), FALSE)
tip_size_flanks <- as.numeric(get_arg(args, "--tip_size_flanks", "0.6"))

# panel spacing (NEW)
panel_gap <- as.numeric(get_arg(args, "--panel_gap", "25"))  # points; try 10–80

# sidebar toggles (NEW)
show_sidebars <- as_logical(get_arg(args, "--show_sidebars", "true"), TRUE)
show_proc  <- as_logical(get_arg(args, "--show_proc",  "true"), TRUE)
show_site  <- as_logical(get_arg(args, "--show_site",  "true"), TRUE)
show_trait <- as_logical(get_arg(args, "--show_trait", "true"), TRUE)
show_allele<- as_logical(get_arg(args, "--show_allele","true"), TRUE)

if (!show_sidebars) {
  show_proc <- FALSE; show_site <- FALSE; show_trait <- FALSE; show_allele <- FALSE
}

need_pheno <- show_proc || show_site || show_trait
need_sites <- show_allele

# figure sizing
fig_w <- as.numeric(get_arg(args, "--width", "22"))
fig_h <- as.numeric(get_arg(args, "--height", "7"))
dpi   <- as.integer(get_arg(args, "--dpi", "300"))
format_out <- get_arg(args, "--format", "png")  # png|pdf|both

# panel widths
w_other   <- as.numeric(get_arg(args, "--w_other", "4"))
w_midblock <- as.numeric(get_arg(args, "--w_midblock", "10"))

# inside hit block: tree | proc | site | trait | allele
w_hit_tree  <- as.numeric(get_arg(args, "--w_hit_tree",  "5"))
w_proc      <- as.numeric(get_arg(args, "--w_proc",      "0.35"))
w_site      <- as.numeric(get_arg(args, "--w_site",      "0.35"))
w_trait     <- as.numeric(get_arg(args, "--w_trait",     "0.45"))
w_allele    <- as.numeric(get_arg(args, "--w_allele",    "0.9"))

# optional spacer inside middle block
mid_spacer <- as.numeric(get_arg(args, "--mid_spacer", "0"))  # set >0 to insert spacer column

if (is.null(smc_file) || is.na(position)) {
  cat("Usage:\n",
      "  Rscript plot_flanks_mid_optional_sidebars_gap_v3.R \\\n",
      "    --smc=FILE.smc --position=983057685 --flank=2 \\\n",
      "    [--sites=FILE.sites --pheno=PHENO.txt] \\\n",
      "    --trait=C13 --allele_positions=983057685 \\\n",
      "    --show_sidebars=true --show_proc=true --show_site=true --show_trait=true --show_allele=true \\\n",
      "    --panel_gap=25 --outdir=out --outprefix=plot --format=png\n", sep="")
  quit(status=1)
}
if (need_pheno && is.null(pheno_file)) stop("Sidebars need phenotype file: pass --pheno=PHENO.txt or set --show_proc=false --show_site=false --show_trait=false")
if (need_sites && is.null(sites_file)) stop("Allele sidebar needs .sites file: pass --sites=FILE.sites or set --show_allele=false")

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

allele_positions <- if (is.null(pos_string) || pos_string == "") position else as.integer(strsplit(pos_string, ",")[[1]])

# --- scan SMC + select trees ---
scan <- scan_smc_tree_intervals(smc_file)
meta <- scan$meta

sel <- pick_hit_plus_flanks(meta, position=position, flank=flank)
hit_idx <- sel$hit
meta_sub <- meta[tree_index %in% sel$idx]
setorder(meta_sub, start)

trees <- read_trees_by_index(
  smc_file=smc_file,
  sample_names=scan$sample_names,
  n_total=nrow(meta),
  keep_indices=meta_sub$tree_index,
  ladderize_tree=ladderize_tr
)

# --- tree plotting function ---
make_tree_plot <- function(tr, st, en, show_labels=FALSE, tip_size=0.8, show_x_title=TRUE){

  branch_mode <- if (tree_scale == "cladogram") "none" else "branch.length"
  seg_len <- en - st + 1L
  ttl <- sprintf("TREE %d–%d (len=%s)", st, en, format(seg_len, big.mark=","))

  p <- ggtree(tr, branch.length=branch_mode, size=0.35) +
    ggtitle(ttl) +
    theme_tree2() +
    theme(plot.title = element_text(size=9, hjust=0.5))

  if (show_labels) p <- p + geom_tiplab(size=tip_size)

  if (tree_scale != "cladogram") {
    if (xlab_mode == "time_before_present") {
      p <- apply_time_before_present_labels(
        p,
        show_title = show_x_title,
        title = "Generations before present (ARGweaver units)"
      )
    } else {
      p <- p + scale_x_continuous(expand=c(0,0))
      if (show_x_title) p <- p + xlab("Tree distance") else p <- p + xlab(NULL) + theme(axis.title.x=element_blank())
    }
  } else {
    if (!show_x_title) p <- p + xlab(NULL) + theme(axis.title.x=element_blank())
  }

  p
}

# --- HIT tree (middle) ---
hit_obj <- trees[[as.character(hit_idx)]]
p_hit_tree <- make_tree_plot(hit_obj$tree, hit_obj$start, hit_obj$end,
                             show_labels=show_tip_labels_hit, tip_size=tip_size_hit,
                             show_x_title=TRUE)

# if no sidebars requested, hit_block is just the tree
hit_block <- p_hit_tree

# ---- build sidebars only if needed ----
if (show_proc || show_site || show_trait || show_allele) {

  # tip y coordinates from hit tree
  tip_df <- as.data.table(p_hit_tree$data[p_hit_tree$data$isTip, c("label","y")])
  setnames(tip_df, c("label","y"))
  tip_df[, label_clean := clean_id(label)]
  tip_df[, label_digits := extract_digits(label)]
  yl <- range(tip_df$y, na.rm=TRUE)

  anno <- NULL
  p_proc <- p_site <- p_trait <- NULL

  if (need_pheno) {
    ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
    stopifnot("codg" %in% names(ph), "proc" %in% names(ph), "site" %in% names(ph))
    if (show_trait && !(trait_name %in% names(ph))) stop("Trait not found in phenotype file: ", trait_name)

    ph[, codg := clean_id(codg)]
    ph[, codg_digits := extract_digits(codg)]

    diag_ph <- choose_key_mode(tip_df$label, ph$codg)
    cat("Phenotype matching: exact=", diag_ph$exact_hits, " digits=", diag_ph$digit_hits,
        " using=", if (diag_ph$use_digits) "DIGITS" else "EXACT", "\n")

    if (diag_ph$use_digits) {
      ph_agg <- ph[, .(
        proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
        site  = suppressWarnings(as.integer(mode1(as.character(site)))),
        trait = if (show_trait) mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE) else NA_real_
      ), by=.(key = codg_digits)]
      tip_df[, key := label_digits]
    } else {
      ph_agg <- ph[, .(
        proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
        site  = suppressWarnings(as.integer(mode1(as.character(site)))),
        trait = if (show_trait) mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE) else NA_real_
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

    if (show_proc) {
      p_proc <- ggplot(anno, aes(x=1, y=y, fill=proc_lab)) +
        geom_tile(width=0.95, height=0.9) +
        scale_fill_brewer(palette="Set2", na.value="grey92", name="Provenances") +
        scale_y_continuous(limits=yl, expand=c(0,0)) +
        theme_void()
    }
    if (show_site) {
      p_site <- ggplot(anno, aes(x=1, y=y, fill=site_lab)) +
        geom_tile(width=0.95, height=0.9) +
        scale_fill_brewer(palette="Set3", na.value="grey92", name="Site") +
        scale_y_continuous(limits=yl, expand=c(0,0)) +
        theme_void()
    }
    if (show_trait) {
      p_trait <- ggplot(anno, aes(x=1, y=y, fill=trait)) +
        geom_tile(width=0.95, height=0.9) +
        scale_fill_viridis_c(na.value="grey92", name=trait_name) +
        scale_y_continuous(limits=yl, expand=c(0,0)) +
        theme_void()
    }
  }

  p_allele <- NULL
  if (need_sites) {
    sx <- read_sites(sites_file)
    geno_wide <- extract_genotypes(sx$names, sx$pos, sx$allele_str, allele_positions)
    geno_long <- melt(as.data.table(geno_wide), id.vars="label", variable.name="position", value.name="allele")
    geno_long[, position := factor(position, levels=as.character(allele_positions))]
    geno_long[, label_clean := clean_id(label)]
    geno_long[, label_digits := extract_digits(label)]

    diag_sites <- choose_key_mode(tip_df$label, sx$names)
    cat("Sites matching: exact=", diag_sites$exact_hits, " digits=", diag_sites$digit_hits,
        " using=", if (diag_sites$use_digits) "DIGITS" else "EXACT", "\n")

    if (diag_sites$use_digits) {
      geno_long[, key := label_digits]
      tip_keys <- tip_df[, .(key=label_digits, label, y)]
    } else {
      geno_long[, key := label_clean]
      tip_keys <- tip_df[, .(key=label_clean, label, y)]
    }
    geno_key <- merge(geno_long[, .(key, position, allele)], tip_keys, by="key", all.x=TRUE, sort=FALSE)
    geno_key <- geno_key[!is.na(label)]

    allele_cols <- c("A"="#1b9e77","C"="#7570b3","G"="#d95f02","T"="#e7298a","N"="#bdbdbd","-"="#252525")

    p_allele <- ggplot(geno_key, aes(x=position, y=y, fill=allele)) +
      geom_tile(height=0.9, color="grey85", linewidth=0.25) +
      scale_fill_manual(values=allele_cols, na.value="grey92", name="Allele") +
      scale_y_continuous(limits=yl, expand=c(0,0)) +
      theme_void() +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
  }

  # ---- assemble middle block dynamically (tree + optional sidebars) ----
  plot_list <- list(p_hit_tree)
  width_list <- c(w_hit_tree)

  if (mid_spacer > 0) {
    plot_list <- c(plot_list, list(plot_spacer()))
    width_list <- c(width_list, mid_spacer)
  }
  if (show_proc  && !is.null(p_proc))   { plot_list <- c(plot_list, list(p_proc));   width_list <- c(width_list, w_proc) }
  if (show_site  && !is.null(p_site))   { plot_list <- c(plot_list, list(p_site));   width_list <- c(width_list, w_site) }
  if (show_trait && !is.null(p_trait))  { plot_list <- c(plot_list, list(p_trait));  width_list <- c(width_list, w_trait) }
  if (show_allele&& !is.null(p_allele)) { plot_list <- c(plot_list, list(p_allele)); width_list <- c(width_list, w_allele) }

  hit_block <- wrap_plots(plot_list, nrow=1, widths=width_list, guides="collect") &
    theme(legend.position="right")
}

# --- Build ALL panels in genomic order (hit ± flank) ---
panel_indices <- meta_sub$tree_index

panels <- lapply(panel_indices, function(i){
  obj <- trees[[as.character(i)]]
  if (i == hit_idx) {
    hit_block
  } else {
    make_tree_plot(obj$tree, obj$start, obj$end,
                   show_labels=show_tip_labels_flanks, tip_size=tip_size_flanks,
                   show_x_title=FALSE)
  }
})

# NEW: add whitespace around each panel
panels <- lapply(panels, add_panel_margin, panel_gap=panel_gap)

widths_vec <- rep(w_other, length(panel_indices))
widths_vec[which(panel_indices == hit_idx)] <- w_midblock

final <- wrap_plots(panels, nrow=1, widths=widths_vec)

# --- Save ---
base <- file.path(outdir, outprefix)
tag <- paste0(".pos", position, ".flank", flank, ".xlab_", xlab_mode,
              ".gap", panel_gap,
              ".proc", show_proc, ".site", show_site, ".trait", show_trait, ".allele", show_allele)

if (format_out %in% c("pdf","both")) {
  ggsave(paste0(base, tag, ".pdf"), final, width=fig_w, height=fig_h)
}
if (format_out %in% c("png","both")) {
  ggsave(paste0(base, tag, ".png"), final, width=fig_w, height=fig_h, dpi=dpi)
}

cat("Wrote outputs under:", outdir, "\n")

