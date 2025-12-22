#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
})

# ============================================================
# Args
# ============================================================
get_arg <- function(flag, default=NULL){
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", gsub("\\-", "\\\\-", flag), "(=|$)"), args)
  if (length(hit) == 0) return(default)

  tok <- args[hit[1]]
  if (grepl("=", tok, fixed = TRUE)) {
    sub("^[^=]+=", "", tok)
  } else {
    if (hit[1] == length(args)) default else args[hit[1] + 1]
  }
}

as_bool <- function(x, default=FALSE){
  if (is.null(x) || length(x) == 0 || is.na(x)) return(default)
  if (is.logical(x)) return(x)
  xx <- tolower(trimws(as.character(x)))
  if (xx %in% c("true","t","1","yes","y")) return(TRUE)
  if (xx %in% c("false","f","0","no","n")) return(FALSE)
  default
}

parse_int_list <- function(x){
  if (is.null(x) || is.na(x) || x == "") return(integer())
  parts <- unlist(strsplit(x, "[,;\\s]+"))
  parts <- parts[parts != ""]
  suppressWarnings(as.integer(parts))
}

# ============================================================
# ID helpers
# ============================================================
clean_id <- function(x){
  x <- as.character(x)
  trimws(x)
}

strip_hap_suffix <- function(x){
  sub("([_./])([12])$", "", clean_id(x))
}

extract_digits_first <- function(x){
  x <- strip_hap_suffix(clean_id(x))
  m <- regexpr("[0-9]+", x)
  if (m[1] == -1) return(NA_integer_)
  suppressWarnings(as.integer(regmatches(x, m)))
}

# ============================================================
# File indexing: outargs_<scaffold>_<start>_<end>.<rep>.(smc|sites)
# ============================================================
parse_outargs_filename <- function(fname){
  b <- basename(fname)
  m <- regexec("^outargs_(.+?)_([0-9]+)_([0-9]+)\\.([0-9]+)\\.(smc|sites)$", b)
  mm <- regmatches(b, m)[[1]]
  if (length(mm) == 0) return(NULL)
  list(
    file = fname,
    scaffold_id = mm[2],
    start = as.numeric(mm[3]),
    end = as.numeric(mm[4]),
    rep = as.integer(mm[5]),
    ext = mm[6]
  )
}

index_dir_files <- function(dir, exts=c("smc","sites")){
  files <- unlist(lapply(exts, function(e) list.files(dir, pattern=paste0("\\.", e, "$"), full.names=TRUE)))
  recs <- lapply(files, parse_outargs_filename)
  recs <- recs[!vapply(recs, is.null, logical(1))]
  if (length(recs) == 0) return(data.table())
  rbindlist(recs)
}

pick_window_safe <- function(smc_index_dt, scaffold_q, position_q){
  dt_sc <- smc_index_dt[scaffold_id == scaffold_q]
  if (nrow(dt_sc) == 0) stop("No .smc files for scaffold ", scaffold_q)

  hits <- dt_sc[position_q >= start & position_q <= end]
  if (nrow(hits) == 0) {
    stop("No window covers position ", position_q, " on ", scaffold_q)
  }
  hits[, span := end - start]
  setorder(hits, span, start, rep)
  list(start=hits$start[1], end=hits$end[1])
}

# ============================================================
# Read .smc (NAMES + TREE segments)
# ============================================================
read_smc_meta <- function(smc_file){
  con <- file(smc_file, open="r")
  on.exit(close(con), add=TRUE)

  l1 <- readLines(con, n=1)
  if (length(l1) == 0 || !startsWith(l1, "NAMES")) stop("SMC missing NAMES line: ", smc_file)
  parts <- strsplit(l1, "\t", fixed=TRUE)[[1]]
  names_vec <- parts[-1]

  lines <- readLines(con)
  tree_lines <- lines[startsWith(lines, "TREE\t")]
  if (length(tree_lines) == 0) stop("No TREE lines found in SMC: ", smc_file)

  dt <- rbindlist(lapply(tree_lines, function(ln){
    p <- strsplit(ln, "\t", fixed=TRUE)[[1]]
    data.table(start=as.numeric(p[2]), end=as.numeric(p[3]), newick=p[4])
  }))
  dt[, idx := .I]
  list(names=names_vec, trees=dt)
}

relabel_tree_tips_from_names <- function(tr, names_vec){
  tl <- tr$tip.label
  if (is.null(tl) || length(tl) == 0) return(tr)
  is_num <- grepl("^[0-9]+$", tl)
  if (all(is_num)) {
    ii <- as.integer(tl) + 1L
    ok <- ii >= 1L & ii <= length(names_vec)
    if (all(ok)) tr$tip.label <- names_vec[ii]
  }
  tr
}

# ============================================================
# Read ARGweaver .sites (accepts NAMES or #NAMES; pos + allele_string)
# ============================================================
open_text_con <- function(path){
  if (grepl("\\.gz$", path, ignore.case=TRUE)) gzfile(path, "rt") else file(path, "r")
}

read_sites_argweaver <- function(sites_file){
  con <- open_text_con(sites_file)
  on.exit(close(con), add=TRUE)
  lines <- readLines(con, warn=FALSE)

  name_line <- lines[grepl("^#?NAMES[ \t]", lines)][1]
  if (is.na(name_line)) stop("Sites missing NAMES/#NAMES line: ", sites_file)

  nm_parts <- strsplit(name_line, "[\t ]+")[[1]]
  names_vec <- nm_parts[-1]
  n <- length(names_vec)
  if (n == 0) stop("Parsed 0 names from NAMES line in ", sites_file)

  dat <- lines[!grepl("^#", lines)]
  dat <- dat[grepl("^\\s*\\d", dat)]
  dat <- dat[nchar(dat) > 0]
  if (length(dat) == 0) stop("No site rows found in: ", sites_file)

  pos <- as.integer(sub("^\\s*([0-9]+).*$", "\\1", dat))
  allele_str <- sub("^\\s*[0-9]+[ \t]+", "", dat)

  bad <- which(nchar(allele_str) != n)
  if (length(bad) > 0) {
    stop("Allele string length != N at line ", bad[1],
         " (len=", nchar(allele_str[bad[1]]), " vs N=", n, ") in ", sites_file)
  }

  list(names=names_vec, pos=pos, allele_str=allele_str)
}

extract_genotypes_from_sites <- function(sites_obj, positions){
  positions <- as.integer(positions)
  idx <- match(positions, sites_obj$pos)
  if (anyNA(idx)) stop("Requested positions not found in .sites: ", paste(positions[is.na(idx)], collapse=", "))

  n <- length(sites_obj$names)
  out <- data.table(label = sites_obj$names)
  for (k in seq_along(idx)) {
    s <- sites_obj$allele_str[idx[k]]
    out[[as.character(positions[k])]] <- substring(s, 1:n, 1:n)
  }
  out
}

# ============================================================
# Phenotype annotation
# ============================================================
read_pheno_anno <- function(pheno_file, trait, tip_labels,
                           id_col="codg", proc_col="proc", site_col="site"){
  tips <- data.table(label = tip_labels)
  if (is.null(pheno_file) || pheno_file == "" || !file.exists(pheno_file)) {
    tips[, proc := NA_character_]
    tips[, site := NA_character_]
    tips[, trait_val := NA_real_]
    return(tips)
  }

  ph <- fread(pheno_file)
  if (!(id_col %in% names(ph))) stop("Pheno missing id_col=", id_col)
  if (!(proc_col %in% names(ph))) stop("Pheno missing proc_col=", proc_col)
  if (!(site_col %in% names(ph))) stop("Pheno missing site_col=", site_col)
  if (!(trait %in% names(ph))) stop("Pheno missing trait column=", trait)

  ph[, codg_clean := clean_id(get(id_col))]
  ph[, codg_base := strip_hap_suffix(codg_clean)]
  ph[, codg_digits := extract_digits_first(codg_base)]

  tips[, label_clean := clean_id(label)]
  tips[, label_base := strip_hap_suffix(label_clean)]
  tips[, label_digits := extract_digits_first(label_base)]

  m1 <- merge(tips, ph, by.x="label_base", by.y="codg_base", all.x=TRUE, sort=FALSE)
  exact_hits <- sum(!is.na(m1[[trait]]))

  digit_hits <- 0
  m2 <- NULL
  if (exact_hits == 0) {
    m2 <- merge(tips, ph, by.x="label_digits", by.y="codg_digits", all.x=TRUE, sort=FALSE)
    digit_hits <- sum(!is.na(m2[[trait]]))
  }

  using <- if (exact_hits > 0) "EXACT_BASE" else if (digit_hits > 0) "DIGITS_BASE" else "NONE"
  cat(sprintf("ID match: exact=%d digit=%d using=%s\n", exact_hits, digit_hits, using))

  out <- if (exact_hits > 0) m1 else if (!is.null(m2)) m2 else tips

  out_dt <- data.table(label = out$label)
  out_dt[, proc := as.character(out[[proc_col]])]
  out_dt[, site := as.character(out[[site_col]])]
  out_dt[, trait_val := suppressWarnings(as.numeric(out[[trait]]))]

  proc_map <- c("1"="Deer Mtn","2"="Inverness River","3"="Judy Creek","4"="Swann Hills","5"="Virginia Hills")
  site_map <- c("1"="JUDY","2"="VIRG","3"="SWAN","4"="TIME")

  if (any(!is.na(out_dt$proc)) && all(grepl("^\\d+$", out_dt$proc[!is.na(out_dt$proc)]))) {
    out_dt[, proc := proc_map[proc]]
  }
  if (any(!is.na(out_dt$site)) && all(grepl("^\\d+$", out_dt$site[!is.na(out_dt$site)]))) {
    out_dt[, site := site_map[site]]
  }

  out_dt
}

# ============================================================
# Plot helpers
# ============================================================
tree_height <- function(tr){
  max(ape::node.depth.edgelength(tr))
}

make_tree_plot <- function(tr, title,
                           show_axis_title=FALSE,
                           axis_title="Generations before present (ARGweaver units)",
                           font_title=11, font_axis_title=10, font_axis_text=9,
                           tip_label_size=2.2, show_tip_labels=FALSE,
                           linewidth_tree=0.35,
                           panel_margin_pt=6){

  tr <- ladderize(tr)
  h <- tree_height(tr)

  p <- ggtree(tr) +
    geom_tree(linewidth=linewidth_tree) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size=font_title, hjust=0.5),
      axis.title.x = if (show_axis_title) element_text(size=font_axis_title) else element_blank(),
      axis.text.x  = element_text(size=font_axis_text),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    )

  breaks <- sort(unique(c(0, pretty(c(0, h), n=5), h)))
  breaks <- breaks[breaks >= 0 & breaks <= h]
  p <- p +
    scale_x_continuous(
      breaks = breaks,
      labels = function(x) format(round(h - x), big.mark=","),
      expand = expansion(mult=c(0.02, 0.02))
    ) +
    xlab(if (show_axis_title) axis_title else NULL)

  if (show_tip_labels) {
    p <- p + geom_tiplab(size=tip_label_size)
  }
  p
}

get_tip_y_map <- function(tree_plot){
  dd <- tree_plot$data
  tip <- dd[dd$isTip, c("label","y")]
  tip <- as.data.table(tip)
  setnames(tip, c("label","y"))
  tip
}

make_discrete_sidebar <- function(dt, tipmap, value_col, legend_title,
                                  palette="Set2",
                                  font_legend=8,
                                  panel_margin_pt=6){
  x <- merge(tipmap, dt[, .(label, val=get(value_col))], by="label", all.x=TRUE, sort=FALSE)
  x[, x := 1]

  ggplot(x, aes(x=x, y=y, fill=val)) +
    geom_tile(width=0.95, height=0.95) +
    theme_void() +
    theme(
      legend.title = element_text(size=font_legend),
      legend.text  = element_text(size=font_legend),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    ) +
    scale_fill_brewer(palette=palette, na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0))
}

make_continuous_sidebar <- function(dt, tipmap, value_col, legend_title,
                                    font_legend=8,
                                    panel_margin_pt=6){
  x <- merge(tipmap, dt[, .(label, val=suppressWarnings(as.numeric(get(value_col))))], by="label", all.x=TRUE, sort=FALSE)
  x[, x := 1]

  ggplot(x, aes(x=x, y=y, fill=val)) +
    geom_tile(width=0.95, height=0.95) +
    theme_void() +
    theme(
      legend.title = element_text(size=font_legend),
      legend.text  = element_text(size=font_legend),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    ) +
    scale_fill_viridis_c(na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0))
}

make_allele_heatmap <- function(geno_long, tipmap, legend_title="Allele",
                                font_axis_text=8, font_legend=8,
                                panel_margin_pt=6){

  x <- merge(tipmap, geno_long, by="label", all.x=TRUE, sort=FALSE)
  x[, position := factor(position, levels=unique(position))]

  allele_cols <- c("A"="#1b9e77","C"="#7570b3","G"="#d95f02","T"="#e7298a","N"="#bdbdbd","-"="#252525")

  ggplot(x, aes(x=position, y=y, fill=allele)) +
    geom_tile(width=0.95, height=0.95, color="grey85", linewidth=0.15) +
    theme_minimal(base_size=font_axis_text) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size=font_legend),
      legend.text  = element_text(size=font_legend),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    ) +
    scale_fill_manual(values=allele_cols, na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0))
}

# ============================================================
# Plot one replicate
# ============================================================
plot_one_replicate <- function(
  smc_file, sites_file,
  position_center,
  allele_positions,
  pheno_file, trait,
  n_up=2, n_down=2,
  show_mid_sidebars=TRUE,
  show_tip_labels_mid=TRUE,
  show_proc=TRUE, show_site=TRUE, show_trait=TRUE, show_allele=TRUE,
  panel_gap=0.25,
  gap_before_allele=0.01,
  w_proc=0.10, w_site=0.10, w_trait=0.16, w_allele=0.16,
  w_tree_mid=1.0,
  w_mid_panel=1.0,
  sidebars_on_right=TRUE,
  allele_near_tree=TRUE,
  font_title=11, font_axis_text=9, font_axis_title=10,
  tip_label_size=2.2,
  linewidth_tree=0.35,
  axis_title="Generations before present (ARGweaver units)",
  axis_label_mid_only=TRUE,
  panel_margin_pt=6
){
  smc <- read_smc_meta(smc_file)
  trees <- smc$trees

  hit <- which(trees$start <= position_center & trees$end >= position_center)
  if (length(hit) == 0) stop("No TREE interval overlaps --position=", position_center, " in ", smc_file)
  hit_i <- hit[1]

  sel_idx <- seq(max(1, hit_i - n_up), min(nrow(trees), hit_i + n_down))
  trees_sel <- trees[idx %in% sel_idx]
  setorder(trees_sel, idx)

  tree_plots <- vector("list", nrow(trees_sel))
  is_mid <- trees_sel$idx == hit_i

  for (i in seq_len(nrow(trees_sel))){
    tr <- read.tree(text=trees_sel$newick[i])
    tr <- relabel_tree_tips_from_names(tr, smc$names)

    title <- sprintf("TREE %.0f–%.0f", trees_sel$start[i], trees_sel$end[i])
    tree_plots[[i]] <- make_tree_plot(
      tr,
      title=title,
      show_axis_title = if (axis_label_mid_only) is_mid[i] else TRUE,
      axis_title=axis_title,
      font_title=font_title,
      font_axis_title=font_axis_title,
      font_axis_text=font_axis_text,
      tip_label_size=tip_label_size,
      show_tip_labels = (is_mid[i] && show_tip_labels_mid),
      linewidth_tree=linewidth_tree,
      panel_margin_pt=panel_margin_pt
    )
  }

  # ---- middle composite (tree + sidebars) ----
  if (show_mid_sidebars) {
    mid_i <- which(is_mid)[1]
    mid_tree_plot <- tree_plots[[mid_i]]
    tipmap <- get_tip_y_map(mid_tree_plot)

    anno <- read_pheno_anno(pheno_file, trait, tipmap$label)

    # Build components
    proc_p <- if (show_proc)  make_discrete_sidebar(anno, tipmap, "proc", "Provenance",
                                                    palette="Set2", font_legend=font_axis_text, panel_margin_pt=panel_margin_pt) else NULL
    site_p <- if (show_site)  make_discrete_sidebar(anno, tipmap, "site", "Site",
                                                    palette="Set3", font_legend=font_axis_text, panel_margin_pt=panel_margin_pt) else NULL
    trait_p<- if (show_trait) make_continuous_sidebar(anno, tipmap, "trait_val", trait,
                                                      font_legend=font_axis_text, panel_margin_pt=panel_margin_pt) else NULL

    allele_p <- NULL
    if (show_allele && !is.null(sites_file) && sites_file != "" && file.exists(sites_file)) {
      sx <- read_sites_argweaver(sites_file)
      geno_wide <- extract_genotypes_from_sites(sx, allele_positions)
      geno_long <- melt(geno_wide, id.vars="label", variable.name="position", value.name="allele")
      geno_long[, label := as.character(label)]
      geno_long[, position := as.character(position)]
      allele_p <- make_allele_heatmap(geno_long, tipmap,
                                      legend_title="Allele",
                                      font_axis_text=font_axis_text,
                                      font_legend=font_axis_text,
                                      panel_margin_pt=panel_margin_pt)
    }

    # Order of sidebars block
    phenos <- list()
    phenos_w <- numeric()
    if (!is.null(proc_p))  { phenos <- c(phenos, list(proc_p));  phenos_w <- c(phenos_w, w_proc) }
    if (!is.null(site_p))  { phenos <- c(phenos, list(site_p));  phenos_w <- c(phenos_w, w_site) }
    if (!is.null(trait_p)) { phenos <- c(phenos, list(trait_p)); phenos_w <- c(phenos_w, w_trait) }

    # Assemble mid composite
    parts <- list()
    widths <- numeric()

    if (sidebars_on_right) {
      # tree first
      parts <- c(parts, list(mid_tree_plot))
      widths <- c(widths, w_tree_mid)

      # then sidebars (optionally allele nearest tree)
      if (!is.null(allele_p) && allele_near_tree) {
        parts <- c(parts, list(allele_p))
        widths <- c(widths, w_allele)

        if (gap_before_allele > 0 && length(phenos) > 0) {
          parts <- c(parts, list(plot_spacer()))
          widths <- c(widths, gap_before_allele)
        }

        if (length(phenos) > 0) {
          parts <- c(parts, phenos)
          widths <- c(widths, phenos_w)
        }
      } else {
        if (length(phenos) > 0) {
          parts <- c(parts, phenos)
          widths <- c(widths, phenos_w)
        }
        if (!is.null(allele_p)) {
          if (gap_before_allele > 0 && length(phenos) > 0) {
            parts <- c(parts, list(plot_spacer()))
            widths <- c(widths, gap_before_allele)
          }
          parts <- c(parts, list(allele_p))
          widths <- c(widths, w_allele)
        }
      }
    } else {
      # original behavior (sidebars left)
      if (length(phenos) > 0) {
        parts <- c(parts, phenos)
        widths <- c(widths, phenos_w)
      }
      if (!is.null(allele_p)) {
        if (gap_before_allele > 0 && length(phenos) > 0) {
          parts <- c(parts, list(plot_spacer()))
          widths <- c(widths, gap_before_allele)
        }
        parts <- c(parts, list(allele_p))
        widths <- c(widths, w_allele)
      }
      parts <- c(parts, list(mid_tree_plot))
      widths <- c(widths, w_tree_mid)
    }

    mid_comp <- wrap_plots(parts, nrow=1, widths=widths) + plot_layout(guides="collect")
    tree_plots[[mid_i]] <- mid_comp
  }

  # ---- Assemble flanks row ----
  panels <- list()
  widths <- numeric()

  for (i in seq_along(tree_plots)) {
    this_w <- if (is_mid[i]) w_mid_panel else 1.0
    panels <- c(panels, list(tree_plots[[i]]))
    widths <- c(widths, this_w)

    if (i < length(tree_plots) && panel_gap > 0) {
      panels <- c(panels, list(plot_spacer()))
      widths <- c(widths, panel_gap)
    }
  }

  wrap_plots(panels, nrow=1, widths=widths) + plot_layout(guides="collect")
}

# ============================================================
# MAIN
# ============================================================
main <- function(){
  smc_dir   <- get_arg("--smc_dir", "")
  sites_dir <- get_arg("--sites_dir", "")
  scaffold_q <- get_arg("--scaffold", get_arg("--scaffold_id", ""))
  position_q <- suppressWarnings(as.numeric(get_arg("--position", NA)))

  pheno_file <- get_arg("--pheno", "")
  trait <- get_arg("--trait", "C13")

  n_up   <- suppressWarnings(as.integer(get_arg("--n_up", NA)))
  n_down <- suppressWarnings(as.integer(get_arg("--n_down", NA)))
  n_flanks <- suppressWarnings(as.integer(get_arg("--n_flanks", 2)))
  if (is.na(n_up)) n_up <- n_flanks
  if (is.na(n_down)) n_down <- n_flanks

  reps <- parse_int_list(get_arg("--replicates", "0"))
  if (length(reps) == 0) reps <- c(0)

  outbase <- get_arg("--outbase", "plots")

  show_mid_sidebars <- as_bool(get_arg("--show_mid_sidebars", get_arg("--show_sidebars", "true")), TRUE)
  show_tip_labels_mid <- as_bool(get_arg("--show_tip_labels_mid", get_arg("--show_tip_labels_mid", "true")), TRUE)

  show_proc  <- as_bool(get_arg("--show_proc",  "true"), TRUE)
  show_site  <- as_bool(get_arg("--show_site",  "true"), TRUE)
  show_trait <- as_bool(get_arg("--show_trait", "true"), TRUE)
  show_allele<- as_bool(get_arg("--show_allele","true"), TRUE)

  sidebars_on_right <- as_bool(get_arg("--sidebars_on_right", "true"), TRUE)
  allele_near_tree  <- as_bool(get_arg("--allele_near_tree",  "true"), TRUE)

  panel_gap <- as.numeric(get_arg("--panel_gap", "0.25"))
  gap_before_allele <- as.numeric(get_arg("--gap_before_allele", "0.01"))

  w_proc  <- as.numeric(get_arg("--w_proc", "0.10"))
  w_site  <- as.numeric(get_arg("--w_site", "0.10"))
  w_trait <- as.numeric(get_arg("--w_trait", "0.16"))
  w_allele<- as.numeric(get_arg("--w_allele", "0.16"))
  w_tree_mid <- as.numeric(get_arg("--w_tree_mid", "1.0"))
  w_mid_panel <- as.numeric(get_arg("--w_mid_panel", "1.0"))

  font_title      <- as.numeric(get_arg("--font_title", get_arg("--title_size","11")))
  font_axis_text  <- as.numeric(get_arg("--font_axis_text", get_arg("--axis_text_size","9")))
  font_axis_title <- as.numeric(get_arg("--font_axis_title", get_arg("--axis_title_size","10")))
  tip_label_size  <- as.numeric(get_arg("--tip_label_size", get_arg("--tiplab_size","2.2")))
  linewidth_tree  <- as.numeric(get_arg("--linewidth_tree","0.35"))

  axis_title <- get_arg("--axis_title", "Generations before present (ARGweaver units)")
  axis_label_mid_only <- as_bool(get_arg("--axis_label_mid_only", "true"), TRUE)
  panel_margin_pt <- as.numeric(get_arg("--panel_margin_pt", "6"))

  pdf_out <- as_bool(get_arg("--pdf","true"), TRUE)
  png_out <- as_bool(get_arg("--png","true"), TRUE)
  fig_w <- as.numeric(get_arg("--width","22"))
  fig_h <- as.numeric(get_arg("--height","11"))
  dpi <- as.numeric(get_arg("--dpi","300"))

  pos_string <- get_arg("--positions", NA)
  allele_positions <- if (!is.na(pos_string)) parse_int_list(pos_string) else as.integer(position_q)
  if (length(allele_positions) == 0) allele_positions <- as.integer(position_q)

  if (smc_dir == "" || scaffold_q == "" || is.na(position_q)) {
    stop("AUTO MODE requires: --smc_dir, --scaffold, --position")
  }

  cat("AUTO MODE\n")

  smc_idx <- index_dir_files(smc_dir, exts="smc")
  if (nrow(smc_idx) == 0) stop("No .smc files found in --smc_dir=", smc_dir)

  sites_idx <- data.table()
  if (sites_dir != "" && dir.exists(sites_dir)) {
    sites_idx <- index_dir_files(sites_dir, exts="sites")
  }

  win <- pick_window_safe(smc_idx, scaffold_q, position_q)
  cat(sprintf("Selected window: %s %.0f-%.0f\n", scaffold_q, win$start, win$end))

  outdir <- file.path(outbase, paste0(scaffold_q, "_", sprintf("%.0f", position_q)))
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

  for (rep_id in reps) {
    smc_hit <- smc_idx[scaffold_id == scaffold_q & start == win$start & end == win$end & rep == rep_id]
    if (nrow(smc_hit) == 0) {
      warning("Rep ", rep_id, ": missing .smc for window ", win$start, "-", win$end, " (skipping)")
      next
    }
    smc_file_rep <- smc_hit$file[1]
    cat(sprintf("Rep %s: %s\n", rep_id, basename(smc_file_rep)))

    sites_file_rep <- ""
    if (nrow(sites_idx) > 0) {
      si <- sites_idx[scaffold_id == scaffold_q & start == win$start & end == win$end & rep == rep_id]
      if (nrow(si) > 0) sites_file_rep <- si$file[1]
    }
    if (sites_file_rep != "") cat("  sites: ", basename(sites_file_rep), "\n", sep="") else cat("  sites: NONE\n")

    p <- plot_one_replicate(
      smc_file=smc_file_rep,
      sites_file=sites_file_rep,
      position_center=position_q,
      allele_positions=allele_positions,
      pheno_file=pheno_file,
      trait=trait,
      n_up=n_up, n_down=n_down,
      show_mid_sidebars=show_mid_sidebars,
      show_tip_labels_mid=show_tip_labels_mid,
      show_proc=show_proc, show_site=show_site, show_trait=show_trait, show_allele=show_allele,
      panel_gap=panel_gap,
      gap_before_allele=gap_before_allele,
      w_proc=w_proc, w_site=w_site, w_trait=w_trait, w_allele=w_allele,
      w_tree_mid=w_tree_mid,
      w_mid_panel=w_mid_panel,
      sidebars_on_right=sidebars_on_right,
      allele_near_tree=allele_near_tree,
      font_title=font_title, font_axis_text=font_axis_text, font_axis_title=font_axis_title,
      tip_label_size=tip_label_size,
      linewidth_tree=linewidth_tree,
      axis_title=axis_title,
      axis_label_mid_only=axis_label_mid_only,
      panel_margin_pt=panel_margin_pt
    )

    outprefix <- file.path(outdir, paste0("rep", rep_id, "_", scaffold_q, "_", sprintf("%.0f", position_q)))
    if (pdf_out) ggsave(paste0(outprefix, ".pdf"), p, width=fig_w, height=fig_h, units="in", dpi=dpi, bg="white")
    if (png_out) ggsave(paste0(outprefix, ".png"), p, width=fig_w, height=fig_h, units="in", dpi=dpi, bg="white")
    cat("Wrote: ", outprefix, if (pdf_out) ".pdf" else "", if (png_out) ".png" else "", "\n", sep="")
  }

  cat("DONE. Output dir: ", outdir, "\n", sep="")
}

main()

