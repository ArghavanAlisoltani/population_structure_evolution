#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# CLI parser (flag=value)
# -----------------------------
get_arg <- function(args, flag, default=NULL){
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if(length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}
as_bool <- function(x, default=FALSE){
  if (is.null(x) || length(x)==0) return(default)
  tolower(trimws(as.character(x))) %in% c("true","t","1","yes","y")
}
split_csv <- function(x){
  if (is.null(x) || x=="") return(character())
  trimws(strsplit(x, ",")[[1]])
}

# -----------------------------
# ID helpers
# -----------------------------
clean_id <- function(x){
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x
}
extract_indiv <- function(hap){
  hap <- clean_id(hap)
  hap <- sub("([_\\.])(1|2)$", "", hap)  # drop _1/_2 or .1/.2 at end
  hap
}
extract_copy <- function(hap){
  hap <- clean_id(hap)
  if (grepl("([_\\.])1$", hap)) return(1L)
  if (grepl("([_\\.])2$", hap)) return(2L)
  return(NA_integer_)
}

# -----------------------------
# Parse filenames
# outargs_<scaffold>_<start>_<end>.<rep>.smc
# outargs_<scaffold>_<start>_<end>.<rep>.sites
# -----------------------------
parse_smc_name <- function(fn){
  b <- basename(fn)
  m <- regexec("^outargs_(.*)_([0-9]+)_([0-9]+)\\.([0-9]+)\\.smc$", b)
  r <- regmatches(b, m)[[1]]
  if (length(r)==0) return(NULL)
  list(scaffold=r[2], start=as.numeric(r[3]), end=as.numeric(r[4]), rep=r[5], base=b)
}

pick_window <- function(smc_dir, scaffold, position, recursive=TRUE){
  smc_files <- list.files(smc_dir, pattern="\\.smc$", full.names=TRUE, recursive=recursive)
  if (length(smc_files)==0) stop("No .smc files found in: ", smc_dir)

  parsed <- lapply(smc_files, parse_smc_name)
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed)==0) stop("No files matched naming: outargs_<scaffold>_<start>_<end>.<rep>.smc")

  dt <- rbindlist(lapply(seq_along(parsed), function(i){
    p <- parsed[[i]]
    data.table(file=smc_files[i], scaffold=p$scaffold, start=p$start, end=p$end, rep=p$rep)
  }))

  dt <- dt[scaffold == scaffold & start <= position & position <= end]
  if (nrow(dt)==0) stop("No window found containing ", scaffold, ":", position)

  dt[, width := end - start]
  setorder(dt, width, start, end)
  dt[1, .(scaffold, start, end)]
}

find_rep_file <- function(dir, scaffold, start, end, rep, ext=c("smc","sites"), recursive=TRUE){
  ext <- match.arg(ext)
  files <- list.files(dir, pattern=paste0("\\.", ext, "$"), full.names=TRUE, recursive=recursive)
  if (length(files)==0) return(NULL)

  if (ext=="smc"){
    ok <- vapply(files, function(f){
      p <- parse_smc_name(f)
      if (is.null(p)) return(FALSE)
      (p$scaffold == scaffold) && (p$start == start) && (p$end == end) && (p$rep == as.character(rep))
    }, logical(1))
  } else {
    ok <- vapply(files, function(f){
      b <- basename(f)
      m <- regexec("^outargs_(.*)_([0-9]+)_([0-9]+)\\.([0-9]+)\\.sites$", b)
      r <- regmatches(b, m)[[1]]
      if (length(r)==0) return(FALSE)
      (r[2]==scaffold) && (as.numeric(r[3])==start) && (as.numeric(r[4])==end) && (r[5]==as.character(rep))
    }, logical(1))
  }
  hits <- files[ok]
  if (length(hits)==0) return(NULL)
  hits[1]
}

# -----------------------------
# Read .smc TREE table
# -----------------------------
read_smc_trees <- function(smc_file){
  lines <- readLines(smc_file, warn=FALSE)
  keep <- grep("^TREE[ \t]", lines, value=TRUE)
  if (length(keep)==0) stop("No TREE lines found in .smc: ", smc_file)

  # Split robustly: tabs or spaces
  sp <- strsplit(keep, "\t", fixed=TRUE)
  bad <- which(vapply(sp, length, integer(1)) < 4)
  if (length(bad)>0){
    sp <- strsplit(keep, "\\s+")
  }

  start <- vapply(sp, function(x) as.numeric(x[2]), numeric(1))
  end   <- vapply(sp, function(x) as.numeric(x[3]), numeric(1))
  newick <- vapply(sp, function(x){
    # newick might contain tabs? unlikely; take remainder past 3 cols
    if (length(x) == 4) return(x[4])
    paste(x[4:length(x)], collapse=" ")
  }, character(1))

  dt <- data.table(tree_index=seq_along(start), start=start, end=end, newick=newick)
  # sort by start just in case
  setorder(dt, start, end)
  dt[, tree_index := seq_len(.N)]
  dt
}

pick_tree_index_for_position <- function(trees, position){
  # treat as inclusive intervals for safety
  hit <- trees[start <= position & position <= end]
  if (nrow(hit)==0){
    # fallback half-open
    hit <- trees[start <= position & position < end]
  }
  if (nrow(hit)==0) return(NA_integer_)
  hit$tree_index[1]
}

# -----------------------------
# Read .sites (ARGweaver format)
# -----------------------------
read_sites <- function(sites_file){
  lines <- readLines(sites_file, warn=FALSE)

  name_idx <- grep("^#NAMES\\s*\\t", lines)
  if (length(name_idx)==0) name_idx <- grep("^#NAMES", lines, ignore.case=TRUE)
  if (length(name_idx)==0) stop("Sites missing NAMES line: ", sites_file)

  name_line <- lines[name_idx[1]]
  names <- strsplit(name_line, "\t", fixed=TRUE)[[1]][-1]
  n <- length(names)
  if (n < 2) stop("Parsed <2 names from #NAMES in: ", sites_file)

  site_lines <- lines[!grepl("^#", lines)]
  site_lines <- site_lines[nchar(site_lines) > 0]

  pos <- as.numeric(sub("\\t.*$", "", site_lines))
  allele_str <- sub("^[0-9]+\\t", "", site_lines)

  bad <- which(nchar(allele_str) != n)
  if (length(bad)>0){
    stop("Allele string length != n in .sites. First bad line: ", bad[1],
         " (n=", n, ", len=", nchar(allele_str[bad[1]]), ")")
  }
  list(names=names, pos=pos, allele_str=allele_str)
}

# -----------------------------
# Site -> tree assignment
# -----------------------------
assign_sites_to_trees <- function(site_pos, trees){
  # trees assumed sorted by start
  # findInterval gives i where start[i] <= pos < start[i+1]
  idx <- findInterval(site_pos, trees$start)
  idx[idx==0] <- NA_integer_

  # drop anything beyond its end (handles inclusive vs half-open confusion)
  ok <- !is.na(idx) & site_pos <= trees$end[idx]
  idx[!ok] <- NA_integer_
  idx
}

# -----------------------------
# Build haplotype-long table for selected sites
# -----------------------------
sites_to_hap_long <- function(names, pos_vec, allele_vec){
  n <- length(names)
  alle_mat <- vapply(allele_vec, function(s) substring(s, seq_len(n), seq_len(n)), character(n))
  dt <- as.data.table(alle_mat)
  dt[, hap := names]
  dt_long <- melt(dt, id.vars="hap", variable.name="site_idx", value.name="allele")
  dt_long[, position := pos_vec[as.integer(sub("^V","", site_idx))]]
  dt_long[, site_idx := NULL]
  dt_long
}

hap_long_to_genotypes <- function(hap_long){
  hap_long[, indiv := extract_indiv(hap)]
  hap_long[, copy := vapply(hap, extract_copy, integer(1))]
  setorder(hap_long, indiv, position, copy, hap)

  geno <- hap_long[, .(
    allele1 = allele[1],
    allele2 = if (.N>=2) allele[2] else NA_character_,
    hap1 = hap[1],
    hap2 = if (.N>=2) hap[2] else NA_character_
  ), by=.(indiv, position)]

  geno[, genotype := {
    a1 <- allele1; a2 <- allele2
    ifelse(is.na(a1) | is.na(a2), NA_character_,
           paste0(sort(c(a1,a2)), collapse=""))
  }, by=.(indiv, position)]

  geno[, is_het := ifelse(is.na(genotype), NA,
                          substr(genotype,1,1) != substr(genotype,2,2))]

  geno
}

# -----------------------------
# Summaries
# -----------------------------
allele_counts_overall <- function(geno){
  a <- rbind(
    geno[!is.na(allele1), .(allele=allele1)],
    geno[!is.na(allele2), .(allele=allele2)]
  )
  a[, .(Count=.N), by=allele][order(-Count)]
}
genotype_counts_overall <- function(geno){
  geno[!is.na(genotype), .(Count=.N), by=genotype][order(-Count)]
}
het_hom_counts_overall <- function(geno){
  geno[!is.na(is_het), .(Count=.N), by=.(state=ifelse(is_het,"HET","HOM"))][order(state)]
}
counts_by_group <- function(geno, meta, group_col){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt <- dt[!is.na(genotype)]
  dt[, .(Count=.N), by=.(group=get(group_col), genotype)][order(group, -Count)]
}
allele_by_group <- function(geno, meta, group_col){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  a <- rbind(
    dt[!is.na(allele1), .(group=get(group_col), allele=allele1)],
    dt[!is.na(allele2), .(group=get(group_col), allele=allele2)]
  )
  a[, .(Count=.N), by=.(group, allele)][order(group, -Count)]
}
trait_by_allele_overall_unique_indiv <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  a <- rbind(
    dt[!is.na(allele1), .(indiv, allele=allele1, traitv)],
    dt[!is.na(allele2), .(indiv, allele=allele2, traitv)]
  )
  a <- a[!is.na(traitv)]
  a <- unique(a, by=c("indiv","allele"))
  a[, .(n=.N, mean=mean(traitv), median=median(traitv)), by=allele][order(-n)]
}
trait_by_genotype_overall_unique_indiv <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  dt <- dt[!is.na(traitv) & !is.na(genotype)]
  dt <- unique(dt[, .(indiv, genotype, traitv)], by=c("indiv","genotype"))
  dt[, .(n=.N, mean=mean(traitv), median=median(traitv)), by=genotype][order(-n)]
}

# -----------------------------
# MAIN
# -----------------------------
args <- commandArgs(trailingOnly=TRUE)

smc_dir    <- get_arg(args, "--smc_dir")
sites_dir  <- get_arg(args, "--sites_dir")
scaffold   <- get_arg(args, "--scaffold")
position   <- as.numeric(get_arg(args, "--position"))
pheno_file <- get_arg(args, "--pheno")
trait      <- get_arg(args, "--trait", "C13")
rep_str    <- get_arg(args, "--replicates", "0")
outbase    <- get_arg(args, "--outbase", "tree_contrib_stats")
recursive  <- as_bool(get_arg(args, "--recursive", "true"), TRUE)
n_near     <- as.integer(get_arg(args, "--n_near", "20"))

id_col   <- get_arg(args, "--id_col", "codg")
proc_col <- get_arg(args, "--proc_col", "proc")
site_col <- get_arg(args, "--site_col", "site")
mum_col  <- get_arg(args, "--mum_col", "mum")

if (is.null(smc_dir) || is.null(sites_dir) || is.null(scaffold) || is.na(position) || is.null(pheno_file)) {
  cat("Usage:\n",
      "  Rscript auto_tree_contrib_sites_stats_v2.R \\\n",
      "    --smc_dir=smc_files --sites_dir=sites_for_tree \\\n",
      "    --scaffold=scaffold_4 --position=983057685 \\\n",
      "    --pheno=PHENO_Charles_6_2025.txt --trait=C13 \\\n",
      "    --replicates=0,10,20,30,40,50 --outbase=tree_contrib_stats\n", sep="")
  quit(status=1)
}

reps <- split_csv(rep_str)
if (length(reps)==0) reps <- "0"

# Output dir
outdir0 <- file.path(outbase, paste0(scaffold, "_", format(position, scientific=FALSE, trim=TRUE)))
dir.create(outdir0, recursive=TRUE, showWarnings=FALSE)

# Auto window
win <- pick_window(smc_dir, scaffold, position, recursive=recursive)
cat("AUTO MODE\n")
cat("Selected window:", win$scaffold, win$start, "-", win$end, "\n")

# Phenotype
ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
stopifnot(id_col %in% names(ph))
ph[, indiv := clean_id(get(id_col))]

meta_cols <- intersect(c(proc_col, site_col, mum_col, trait), names(ph))
meta <- ph[, c("indiv", meta_cols), with=FALSE]

for (rep in reps) {
  rep_out <- file.path(outdir0, paste0("rep", rep))
  dir.create(rep_out, recursive=TRUE, showWarnings=FALSE)

  smc_file <- find_rep_file(smc_dir, win$scaffold, win$start, win$end, rep, ext="smc", recursive=recursive)
  sites_file <- find_rep_file(sites_dir, win$scaffold, win$start, win$end, rep, ext="sites", recursive=recursive)

  if (is.null(smc_file)) { cat("Rep", rep, ": missing .smc\n"); next }
  if (is.null(sites_file)) { cat("Rep", rep, ": missing .sites\n"); next }

  cat("Rep", rep, ":\n  smc:", basename(smc_file), "\n  sites:", basename(sites_file), "\n")

  trees <- read_smc_trees(smc_file)
  sel_tree <- pick_tree_index_for_position(trees, position)
  if (is.na(sel_tree)) stop("No TREE interval overlaps position in: ", smc_file)

  # Sites load
  sx <- read_sites(sites_file)
  site_pos <- as.numeric(sx$pos)
  allele_str <- sx$allele_str
  hap_names <- sx$names

  # Diagnostics: range + total
  sites_summary <- data.table(
    scaffold=scaffold,
    query_position=position,
    window_start=win$start, window_end=win$end,
    smc_file=smc_file,
    sites_file=sites_file,
    n_sites_total=length(site_pos),
    sites_min=min(site_pos),
    sites_max=max(site_pos),
    n_trees=nrow(trees),
    selected_tree_index=sel_tree,
    selected_tree_start=trees[tree_index==sel_tree]$start,
    selected_tree_end=trees[tree_index==sel_tree]$end
  )
  fwrite(sites_summary, file.path(rep_out, "rep_summary_diagnostics.tsv"), sep="\t")

  # Assign each site to a tree interval
  site_tree_idx <- assign_sites_to_trees(site_pos, trees)

  # tree site counts
  trees[, n_sites := 0L]
  tab <- as.data.table(table(site_tree_idx, useNA="no"))
  if (nrow(tab)>0) {
    tab[, site_tree_idx := as.integer(as.character(site_tree_idx))]
    setnames(tab, c("site_tree_idx","N"), c("tree_index","n_sites"))
    trees[tab, n_sites := i.n_sites, on="tree_index"]
  }
  trees[, is_selected := (tree_index == sel_tree)]

  # Write all trees + site counts + selected flag
  fwrite(trees[, .(tree_index, start, end, n_sites, is_selected)],
         file.path(rep_out, "trees_all_with_site_counts.tsv"), sep="\t")

  # Save selected tree newick + info
  sel_info <- trees[tree_index==sel_tree, .(tree_index, start, end, n_sites)]
  fwrite(sel_info, file.path(rep_out, "selected_tree_info.tsv"), sep="\t")
  writeLines(trees[tree_index==sel_tree]$newick, file.path(rep_out, "selected_tree.nwk"))

  # List the exact positions contributing to selected tree
  contrib_idx <- which(site_tree_idx == sel_tree)
  contrib_pos <- site_pos[contrib_idx]
  fwrite(data.table(position=contrib_pos),
         file.path(rep_out, "selected_tree_site_positions.tsv"), sep="\t")

  # Also: nearest sites around query position (helps explain 0-site cases)
  ord <- order(abs(site_pos - position))
  near <- ord[seq_len(min(n_near, length(ord)))]
  near_dt <- data.table(
    position=site_pos[near],
    dist_to_query=abs(site_pos[near]-position),
    tree_index=site_tree_idx[near]
  )
  near_dt <- merge(near_dt, trees[, .(tree_index, tree_start=start, tree_end=end)], by="tree_index", all.x=TRUE)
  setorder(near_dt, dist_to_query)
  fwrite(near_dt, file.path(rep_out, "nearest_sites_to_query_with_tree_interval.tsv"), sep="\t")

  # If no sites mapped to selected tree, stop after writing diagnostics (this may be REAL)
  if (length(contrib_idx) == 0) {
    cat("  NOTE: 0 sites mapped to selected tree interval. See diagnostics TSVs in: ", rep_out, "\n", sep="")
    next
  }

  # Build genotype tables for contributing sites
  hap_long <- sites_to_hap_long(hap_names, contrib_pos, allele_str[contrib_idx])
  geno <- hap_long_to_genotypes(hap_long)
  geno_meta <- merge(geno, meta, by="indiv", all.x=TRUE)

  fwrite(geno_meta[order(position, indiv)],
         file.path(rep_out, "genotypes_long_per_position.tsv"), sep="\t")

  fwrite(allele_counts_overall(geno), file.path(rep_out, "allele_overall.tsv"), sep="\t")
  fwrite(genotype_counts_overall(geno), file.path(rep_out, "genotype_overall.tsv"), sep="\t")
  fwrite(het_hom_counts_overall(geno), file.path(rep_out, "het_hom_overall.tsv"), sep="\t")

  # Group summaries if columns exist
  if (proc_col %in% names(meta)) {
    fwrite(counts_by_group(geno, meta, proc_col), file.path(rep_out, "genotype_by_provenance.tsv"), sep="\t")
    fwrite(allele_by_group(geno, meta, proc_col), file.path(rep_out, "allele_by_provenance.tsv"), sep="\t")
  }
  if (site_col %in% names(meta)) {
    fwrite(counts_by_group(geno, meta, site_col), file.path(rep_out, "genotype_by_site.tsv"), sep="\t")
    fwrite(allele_by_group(geno, meta, site_col), file.path(rep_out, "allele_by_site.tsv"), sep="\t")
  }
  if (mum_col %in% names(meta)) {
    fwrite(counts_by_group(geno, meta, mum_col), file.path(rep_out, "genotype_by_mum.tsv"), sep="\t")
    fwrite(allele_by_group(geno, meta, mum_col), file.path(rep_out, "allele_by_mum.tsv"), sep="\t")
  }

  # Trait summaries if present
  if (trait %in% names(meta)) {
    fwrite(trait_by_allele_overall_unique_indiv(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_allele_overall_unique_indiv.tsv")), sep="\t")
    fwrite(trait_by_genotype_overall_unique_indiv(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_genotype_overall_unique_indiv.tsv")), sep="\t")
  }

  cat("  Wrote outputs to: ", rep_out, "\n", sep="")
}

cat("DONE. Output base dir: ", outdir0, "\n", sep="")

