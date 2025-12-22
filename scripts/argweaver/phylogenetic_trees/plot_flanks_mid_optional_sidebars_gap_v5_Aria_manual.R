#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
})

# ---------------------------
# Helpers
# ---------------------------
get_arg <- function(args, flag, default=NULL){
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if(length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

as_bool <- function(x, default=FALSE){
  if (is.null(x)) return(default)
  x <- tolower(trimws(x))
  x %in% c("1","true","t","yes","y")
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

# ---------------------------
# Read ARGweaver .smc (TREE + NAMES)
# ---------------------------
read_smc_names <- function(smc_file){
  con <- file(smc_file, "r")
  on.exit(close(con), add=TRUE)
  repeat{
    ln <- readLines(con, n=1)
    if(length(ln)==0) stop("Could not find NAMES line in .smc")
    if(grepl("^(#)?NAMES\\t", ln)){
      parts <- strsplit(ln, "\t", fixed=TRUE)[[1]]
      nm <- parts[-1]
      return(nm)
    }
  }
}

read_smc_trees_table <- function(smc_file){
  # Reads only TREE lines; returns data.table(tree_i, start, end, newick)
  lines <- readLines(smc_file)
  tree_lines <- lines[grepl("^TREE\\t", lines)]
  if(length(tree_lines)==0) stop("No TREE lines found in .smc")

  starts <- integer(length(tree_lines))
  ends   <- integer(length(tree_lines))
  newick <- character(length(tree_lines))

  for(i in seq_along(tree_lines)){
    parts <- strsplit(tree_lines[i], "\t", fixed=TRUE)[[1]]
    # TREE  start  end  newick
    starts[i] <- as.integer(parts[2])
    ends[i]   <- as.integer(parts[3])
    newick[i] <- parts[4]
  }
  data.table(tree_i = seq_along(tree_lines),
             start = starts, end = ends, newick = newick)
}

relabel_tree_tips_from_names <- function(tr, names_vec){
  # ARGweaver often uses 0..n-1 as tip labels in the TREE Newick.
  # Map those indices to names_vec (0-based).
  tl <- tr$tip.label
  is_idx <- grepl("^\\d+$", tl)
  if(any(is_idx)){
    idx <- as.integer(tl[is_idx])
    # 0-based
    mapped <- rep(NA_character_, length(idx))
    ok <- idx >= 0 & idx < length(names_vec)
    mapped[ok] <- names_vec[idx[ok] + 1]
    # if something is out of range, keep original
    mapped[!ok] <- tl[is_idx][!ok]
    tl[is_idx] <- mapped
    tr$tip.label <- tl
  }
  tr
}

tree_height <- function(tr){
  # max root-to-tip distance (uses edge lengths)
  # node.depth.edgelength gives distance from root to each node
  max(ape::node.depth.edgelength(tr))
}

# ---------------------------
# Read .sites (alleles)
# ---------------------------
read_sites <- function(sites_file) {
  lines <- readLines(sites_file)
  name_line <- lines[grepl("^#?NAMES\\t", lines)][1]
  if (is.na(name_line)) stop("No NAMES/#NAMES line found. Is this an ARGweaver .sites file?")
  names <- strsplit(name_line, "\t", fixed=TRUE)[[1]][-1]
  n <- length(names)

  site_lines <- lines[!grepl("^#", lines)]
  site_lines <- site_lines[nchar(site_lines) > 0]
  pos <- as.integer(sub("\\t.*$", "", site_lines))
  allele_str <- sub("^[0-9]+\\t", "", site_lines)

  bad <- which(nchar(allele_str) != n)
  if (length(bad) > 0) stop("Allele string length != n in .sites. First bad line: ", bad[1])

  list(names=names, pos=pos, allele_str=allele_str)
}

extract_genotypes <- function(names, pos, allele_str, positions_numeric) {
  n <- length(names)
  idx <- match(positions_numeric, pos)
  if (anyNA(idx)) stop("Requested positions not found in .sites: ",
                       paste(positions_numeric[is.na(idx)], collapse=", "))

  geno <- data.table(label = names)
  for (k in seq_along(idx)) {
    s <- allele_str[idx[k]]
    geno[[as.character(positions_numeric[k])]] <- substring(s, seq_len(n), seq_len(n))
  }
  geno
}

# ---------------------------
# Phenotype join
# ---------------------------
read_pheno_anno <- function(pheno_file, trait_name, tip_labels, outprefix){
  ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
  stopifnot("codg" %in% names(ph), "proc" %in% names(ph), "site" %in% names(ph))
  if (!(trait_name %in% names(ph))) stop("Trait not found in pheno: ", trait_name)

  ph[, codg := clean_id(codg)]
  ph[, codg_digits := extract_digits(codg)]

  tips <- data.table(label = tip_labels)
  tips[, label_clean  := clean_id(label)]
  tips[, label_digits := extract_digits(label)]

  exact_hits <- sum(tips$label_clean %in% ph$codg)
  digit_hits <- sum(!is.na(tips$label_digits) & (tips$label_digits %in% ph$codg_digits))
  use_digits <- (digit_hits > exact_hits)

  cat("ID matching diagnostics\n")
  cat("  tips:", nrow(tips), "\n")
  cat("  exact matches (clean string):", exact_hits, "\n")
  cat("  digit matches (extract digits):", digit_hits, "\n")
  cat("  using key:", if (use_digits) "DIGITS" else "EXACT", "\n")

  if (!use_digits) {
    ph_agg <- ph[, .(
      proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
      site  = suppressWarnings(as.integer(mode1(as.character(site)))),
      mum   = mode1(as.character(mum)),
      trait = mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE)
    ), by=.(key = codg)]
    tips[, key := label_clean]
  } else {
    ph_agg <- ph[, .(
      proc  = suppressWarnings(as.integer(mode1(as.character(proc)))),
      site  = suppressWarnings(as.integer(mode1(as.character(site)))),
      mum   = mode1(as.character(mum)),
      trait = mean(suppressWarnings(as.numeric(get(trait_name))), na.rm=TRUE)
    ), by=.(key = codg_digits)]
    tips[, key := label_digits]
  }

  anno <- merge(tips, ph_agg, by="key", all.x=TRUE, sort=FALSE)

  unmatched <- anno[is.na(proc) & is.na(site) & is.na(trait), label]
  writeLines(unmatched, paste0(outprefix, ".unmatched_tree_tips.txt"))
  cat("  unmatched tips written to:", paste0(outprefix, ".unmatched_tree_tips.txt"), "\n")

  proc_map <- c("1"="Deer Mtn",
                "2"="Inverness River",
                "3"="Judy Creek",
                "4"="Swann Hills",
                "5"="Virginia Hills")
  site_map <- c("1"="JUDY",
                "2"="VIRG",
                "3"="SWAN",
                "4"="TIME")

  anno[, proc_lab := factor(proc_map[as.character(proc)],
                            levels=c("Deer Mtn","Inverness River","Judy Creek","Swann Hills","Virginia Hills"))]
  anno[, site_lab := factor(site_map[as.character(site)],
                            levels=c("JUDY","VIRG","SWAN","TIME"))]

  list(anno=anno, use_digits=use_digits)
}

# ---------------------------
# Plot pieces
# ---------------------------
make_tree_plot <- function(tr, title,
                           show_xlab=FALSE,
                           xlab_text="Generations before present (ARGweaver units)",
                           show_tip_labels=FALSE,
                           tip_label_size=2.0,
                           axis_text_size=8,
                           axis_title_size=9,
                           title_size=10,
                           branch_linewidth=0.4,
                           panel_margin_pt=6){

  tr <- ladderize(tr)

  h <- tree_height(tr)

  # Build base tree (root at 0; tips at h)
  p <- ggtree(tr, size=branch_linewidth) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size=title_size, hjust=0.5),
      axis.text.x = element_text(size=axis_text_size),
      axis.title.x = element_text(size=axis_title_size),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt),
      legend.position = "none"
    )

  # Keep geometry as-is; only relabel axis to show "time before present"
  # so that rightmost (tips) becomes 0, leftmost (root) becomes older.
  breaks <- sort(unique(c(0, pretty(c(0,h), n=4), h)))
  p <- p +
    scale_x_continuous(
      breaks = breaks,
      labels = function(x) format(round(h - x), big.mark=","),
      expand = expansion(mult=c(0.01, 0.01))
    )

  if(show_xlab){
    p <- p + xlab(xlab_text)
  } else {
    p <- p + xlab(NULL)
  }

  if(show_tip_labels){
    # show original tip IDs (haplotypes)
    p <- p +
      geom_tiplab(size=tip_label_size, align=FALSE) +
      coord_cartesian(clip="off") +
      theme(plot.margin = margin(panel_margin_pt, panel_margin_pt + 20, panel_margin_pt, panel_margin_pt))
  }

  p
}

make_tile_bar_discrete <- function(df, y_map, col_name, value_col,
                                  legend_title,
                                  palette="Set2",
                                  w=1,
                                  legend_title_size=9,
                                  legend_text_size=8,
                                  panel_margin_pt=6){
  dt <- merge(y_map, df[, .(label, val = get(value_col))], by="label", all.x=TRUE, sort=FALSE)
  dt[, x := 1L]

  ggplot(dt, aes(x=x, y=y, fill=val)) +
    geom_tile(width=0.95, height=0.95) +
    scale_fill_brewer(palette=palette, na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme_void() +
    theme(
      legend.title = element_text(size=legend_title_size),
      legend.text  = element_text(size=legend_text_size),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    )
}

make_tile_bar_continuous <- function(df, y_map, value_col,
                                     legend_title,
                                     legend_title_size=9,
                                     legend_text_size=8,
                                     panel_margin_pt=6){
  dt <- merge(y_map, df[, .(label, val = suppressWarnings(as.numeric(get(value_col))))], by="label", all.x=TRUE, sort=FALSE)
  dt[, x := 1L]

  ggplot(dt, aes(x=x, y=y, fill=val)) +
    geom_tile(width=0.95, height=0.95) +
    scale_fill_viridis_c(na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme_void() +
    theme(
      legend.title = element_text(size=legend_title_size),
      legend.text  = element_text(size=legend_text_size),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    )
}

make_allele_heatmap <- function(geno_long, y_map,
                                allele_colors,
                                legend_title="Allele",
                                axis_text_size=8,
                                legend_title_size=9,
                                legend_text_size=8,
                                panel_margin_pt=6){
  dt <- merge(y_map, geno_long, by="label", all.x=TRUE, sort=FALSE)
  dt[, position := factor(position, levels=unique(position))]

  ggplot(dt, aes(x=position, y=y, fill=allele)) +
    geom_tile(width=0.95, height=0.95, color="grey85", linewidth=0.15) +
    scale_fill_manual(values=allele_colors, na.value="grey92", name=legend_title) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_minimal(base_size = axis_text_size) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size=legend_title_size),
      legend.text  = element_text(size=legend_text_size),
      plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt)
    )
}

# ---------------------------
# Args
# ---------------------------
args <- commandArgs(trailingOnly=TRUE)

smc_file   <- get_arg(args, "--smc")
sites_file <- get_arg(args, "--sites", NULL)
pheno_file <- get_arg(args, "--pheno", NULL)

pos_center <- as.integer(get_arg(args, "--position"))
allele_pos_string <- get_arg(args, "--positions", if(!is.na(pos_center)) as.character(pos_center) else NULL)

trait_name <- get_arg(args, "--trait", "C13")
n_flanks   <- as.integer(get_arg(args, "--n_flanks", "1"))

outprefix  <- get_arg(args, "--outprefix", "flanks_mid_tree")
out_png    <- as_bool(get_arg(args, "--out_png", "true"), TRUE)
out_pdf    <- as_bool(get_arg(args, "--out_pdf", "false"), FALSE)

show_mid_sidebars <- as_bool(get_arg(args, "--show_mid_sidebars", "true"), TRUE)
show_proc  <- as_bool(get_arg(args, "--show_proc", "true"), TRUE)
show_site  <- as_bool(get_arg(args, "--show_site", "true"), TRUE)
show_trait <- as_bool(get_arg(args, "--show_trait", "true"), TRUE)
show_allele <- as_bool(get_arg(args, "--show_allele", "true"), TRUE)

show_tip_labels_mid <- as_bool(get_arg(args, "--show_tip_labels_mid", "true"), TRUE)
tip_label_size <- as.numeric(get_arg(args, "--tip_label_size", "2.0"))

# font controls
title_size <- as.numeric(get_arg(args, "--font_title", "10"))
axis_text_size <- as.numeric(get_arg(args, "--font_axis_text", "8"))
axis_title_size <- as.numeric(get_arg(args, "--font_axis_title", "9"))
legend_title_size <- as.numeric(get_arg(args, "--font_legend_title", "9"))
legend_text_size  <- as.numeric(get_arg(args, "--font_legend_text", "8"))

# spacing controls
panel_gap <- as.numeric(get_arg(args, "--panel_gap", "0.5"))          # outer gap between tree-panels
panel_margin_pt <- as.numeric(get_arg(args, "--panel_margin_pt", "6"))# whitespace around each panel (points)

sidebar_gap <- as.numeric(get_arg(args, "--sidebar_gap", "0.05"))     # gap between proc/site/trait
gap_before_allele <- as.numeric(get_arg(args, "--gap_before_allele", "0.05")) # reduce this to tighten allele

# width controls (relative)
w_tree  <- as.numeric(get_arg(args, "--w_tree", "4.5"))
w_proc  <- as.numeric(get_arg(args, "--w_proc", "0.6"))
w_site  <- as.numeric(get_arg(args, "--w_site", "0.6"))
w_trait <- as.numeric(get_arg(args, "--w_trait", "0.8"))
w_allele <- as.numeric(get_arg(args, "--w_allele", "1.2"))
w_gap_outer <- as.numeric(get_arg(args, "--w_gap_outer", as.character(panel_gap)))
w_gap_sidebar <- as.numeric(get_arg(args, "--w_gap_sidebar", as.character(sidebar_gap)))
w_gap_allele  <- as.numeric(get_arg(args, "--w_gap_allele", as.character(gap_before_allele)))

if (is.null(smc_file) || is.na(pos_center) || is.null(allele_pos_string)) {
  cat("Usage:\n",
      "  Rscript plot_flanks_mid_optional_sidebars_gap_v5.R \\\n",
      "    --smc=FILE.smc --position=983057685 --n_flanks=2 \\\n",
      "    --sites=FILE.sites --pheno=PHENO.txt --positions=983057685 --trait=C13 \\\n",
      "    --show_mid_sidebars=true --show_tip_labels_mid=true \\\n",
      "    --outprefix=OUT\n", sep="")
  quit(status=1)
}

allele_positions <- as.integer(strsplit(allele_pos_string, ",")[[1]])
if (any(is.na(allele_positions))) stop("Bad --positions (must be comma-separated integers).")

# ---------------------------
# Load SMC trees + pick center index + flanks
# ---------------------------
names_vec <- read_smc_names(smc_file)
tree_tbl  <- read_smc_trees_table(smc_file)

hit <- tree_tbl[start <= pos_center & end >= pos_center]
if (nrow(hit) == 0) stop("No TREE interval overlaps --position=", pos_center)
hit_i <- hit$tree_i[1]

sel_i <- seq(max(1, hit_i - n_flanks), min(nrow(tree_tbl), hit_i + n_flanks))
sel_tbl <- tree_tbl[tree_i %in% sel_i]
sel_tbl[, panel_role := ifelse(tree_i == hit_i, "MID", "FLANK")]

# ---------------------------
# Build tree plots (flanks + mid)
# ---------------------------
tree_plots <- list()

# If we will add sidebars, we need the mid tree y-map and annotations
mid_tree <- NULL
mid_tree_plot <- NULL
mid_comp <- NULL

allele_cols <- c("A"="#1b9e77","C"="#7570b3","G"="#d95f02","T"="#e7298a","N"="#bdbdbd","-"="#252525")

for (k in seq_len(nrow(sel_tbl))) {
  row <- sel_tbl[k]
  tr <- ape::read.tree(text=row$newick)
  tr <- relabel_tree_tips_from_names(tr, names_vec)

  title_txt <- sprintf("TREE %d–%d (len=%s)", row$start, row$end, format(row$end - row$start + 1, big.mark=","))

  is_mid <- (row$tree_i == hit_i)
  show_xlab_here <- is_mid  # only middle gets x-axis title

  p_tr <- make_tree_plot(
    tr, title=title_txt,
    show_xlab=show_xlab_here,
    xlab_text="Generations before present (ARGweaver units)",
    show_tip_labels=(is_mid && show_tip_labels_mid && (!show_mid_sidebars)), # if sidebars ON, labels often clutter; handled below
    tip_label_size=tip_label_size,
    axis_text_size=axis_text_size,
    axis_title_size=axis_title_size,
    title_size=title_size,
    panel_margin_pt=panel_margin_pt
  )

  if (is_mid) {
    mid_tree <- tr

    # We'll rebuild mid tree plot with tip labels allowed (works fine with sidebars because they are separate plots)
    mid_tree_plot <- make_tree_plot(
      tr, title=title_txt,
      show_xlab=TRUE,
      xlab_text="Generations before present (ARGweaver units)",
      show_tip_labels=show_tip_labels_mid,
      tip_label_size=tip_label_size,
      axis_text_size=axis_text_size,
      axis_title_size=axis_title_size,
      title_size=title_size,
      panel_margin_pt=panel_margin_pt
    )

    # add mid later (maybe with sidebars)
    tree_plots[[length(tree_plots)+1]] <- NULL
  } else {
    tree_plots[[length(tree_plots)+1]] <- p_tr
  }
}

# ---------------------------
# Build middle composite (optional sidebars)
# ---------------------------
# figure out where the mid panel goes in the final row
mid_pos_in_sel <- which(sel_tbl$tree_i == hit_i)

if (!is.null(mid_tree_plot)) {

  if (show_mid_sidebars) {
    # y-map from the mid tree plot (tip coordinates)
    tip_df <- mid_tree_plot$data[mid_tree_plot$data$isTip, c("label","y")]
    y_map <- data.table(label=as.character(tip_df$label), y=tip_df$y)
    y_map <- y_map[order(y)]  # stable order

    sidebar_list <- list()

    # Phenotype bars (optional)
    if (!is.null(pheno_file) && file.exists(pheno_file)) {
      ph <- read_pheno_anno(pheno_file, trait_name, y_map$label, outprefix)
      anno <- ph$anno

      if (show_proc) {
        p_proc <- make_tile_bar_discrete(
          df=anno, y_map=y_map,
          col_name="proc", value_col="proc_lab",
          legend_title="Provenances",
          palette="Set2",
          legend_title_size=legend_title_size,
          legend_text_size=legend_text_size,
          panel_margin_pt=panel_margin_pt
        )
        sidebar_list <- c(sidebar_list, list(p_proc))
      }

      if (show_site) {
        p_site <- make_tile_bar_discrete(
          df=anno, y_map=y_map,
          col_name="site", value_col="site_lab",
          legend_title="Site",
          palette="Set3",
          legend_title_size=legend_title_size,
          legend_text_size=legend_text_size,
          panel_margin_pt=panel_margin_pt
        )
        sidebar_list <- c(sidebar_list, list(p_site))
      }

      if (show_trait) {
        p_trait <- make_tile_bar_continuous(
          df=anno, y_map=y_map,
          value_col="trait",
          legend_title=trait_name,
          legend_title_size=legend_title_size,
          legend_text_size=legend_text_size,
          panel_margin_pt=panel_margin_pt
        )
        sidebar_list <- c(sidebar_list, list(p_trait))
      }
    }

    # Allele heatmap (optional)
    if (show_allele && !is.null(sites_file) && file.exists(sites_file)) {
      sx <- read_sites(sites_file)

      geno_wide <- extract_genotypes(sx$names, sx$pos, sx$allele_str, allele_positions)
      geno_long <- melt(geno_wide, id.vars="label", variable.name="position", value.name="allele")
      geno_long[, label := as.character(label)]
      geno_long[, position := factor(position, levels=as.character(allele_positions))]

      # Keep only labels present in tree (avoid mismatches)
      geno_long <- geno_long[label %in% y_map$label]

      p_allele <- make_allele_heatmap(
        geno_long=geno_long, y_map=y_map,
        allele_colors=allele_cols,
        legend_title="Allele",
        axis_text_size=axis_text_size,
        legend_title_size=legend_title_size,
        legend_text_size=legend_text_size,
        panel_margin_pt=panel_margin_pt
      )

      # Optionally insert a small gap before allele (user request)
      if (length(sidebar_list) > 0 && w_gap_allele > 0) {
        sidebar_list <- c(sidebar_list, list(plot_spacer()), list(p_allele))
      } else {
        sidebar_list <- c(sidebar_list, list(p_allele))
      }
    }

    # Build mid composite with controlled sidebar gaps
    # Insert spacers between phenotype bars (uniform small gap)
    assembled <- list(mid_tree_plot)
    widths <- c(w_tree)

    if (length(sidebar_list) > 0) {
      for (i in seq_along(sidebar_list)) {
        obj <- sidebar_list[[i]]
        if (inherits(obj, "gg")) {
          assembled <- c(assembled, list(obj))
          # choose width based on which sidebar it is
          # (rough matching: proc/site/trait/allele)
          widths <- c(widths,
                      ifelse(i==1, w_proc,
                             ifelse(i==2, w_site,
                                    ifelse(i==3, w_trait, w_allele))))
          # insert gap after each real sidebar except last
          if (i < length(sidebar_list)) {
            assembled <- c(assembled, list(plot_spacer()))
            # allele-gap uses w_gap_allele where that spacer was inserted as plot_spacer() already
            widths <- c(widths,
                        ifelse(i == length(sidebar_list)-1 && show_allele, w_gap_allele, w_gap_sidebar))
          }
        } else {
          # plot_spacer (for gap_before_allele)
          assembled <- c(assembled, list(obj))
          widths <- c(widths, w_gap_allele)
        }
      }
    }

    mid_comp <- wrap_plots(assembled, nrow=1, widths=widths) +
      plot_layout(guides="collect") &
      theme(legend.position="right")

  } else {
    mid_comp <- mid_tree_plot
  }
}

# ---------------------------
# Assemble the full row (with adjustable gaps between tree-panels)
# ---------------------------
outer_list <- list()
outer_widths <- c()

# iterate selected trees in order; replace mid with composite
for (k in seq_len(nrow(sel_tbl))) {
  is_mid <- (sel_tbl$tree_i[k] == hit_i)

  if (is_mid) {
    outer_list <- c(outer_list, list(mid_comp))
    outer_widths <- c(outer_widths, w_tree + w_proc + w_site + w_trait + w_allele + 1.2) # rough; patchwork will adapt
  } else {
    # find next non-null flank plot in tree_plots
    # tree_plots contains only flanks in order, with NULL placeholder for mid
    idx <- sum(sel_tbl$tree_i[1:k] != hit_i)
    outer_list <- c(outer_list, list(tree_plots[[idx]]))
    outer_widths <- c(outer_widths, w_tree)
  }

  # add outer spacer between tree panels (not after last)
  if (k < nrow(sel_tbl)) {
    outer_list <- c(outer_list, list(plot_spacer()))
    outer_widths <- c(outer_widths, w_gap_outer)
  }
}

final_plot <- wrap_plots(outer_list, nrow=1, widths=outer_widths) &
  theme(plot.margin = margin(panel_margin_pt, panel_margin_pt, panel_margin_pt, panel_margin_pt))

# ---------------------------
# Save
# ---------------------------
if (out_png) {
  ggsave(paste0(outprefix, ".png"), final_plot, width=16, height=14, dpi=300)
}
if (out_pdf) {
  ggsave(paste0(outprefix, ".pdf"), final_plot, width=16, height=14)
}

cat("Wrote:\n")
if (out_png) cat("  ", paste0(outprefix, ".png"), "\n", sep="")
if (out_pdf) cat("  ", paste0(outprefix, ".pdf"), "\n", sep="")
cat("  ", paste0(outprefix, ".unmatched_tree_tips.txt"), "\n", sep="")

