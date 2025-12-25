#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Small CLI parser (flag=value)
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
# ID helpers (haplotypes -> individual)
# -----------------------------
clean_id <- function(x){
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x
}

# Extract base individual ID from haplotype tip label
# Examples: "5001_1" -> "5001", "5001_2" -> "5001"
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

mode1 <- function(x){
  x <- x[!is.na(x)]
  if (length(x)==0) return(NA_character_)
  tb <- sort(table(x), decreasing=TRUE)
  names(tb)[1]
}

# -----------------------------
# Parse filenames and auto-pick window
# Expected: outargs_<scaffold>_<start>_<end>.<rep>.smc
# -----------------------------
parse_outargs_name <- function(fn){
  b <- basename(fn)
  m <- regexec("^outargs_(.*)_([0-9]+)_([0-9]+)\\.([0-9]+)\\.smc$", b)
  r <- regmatches(b, m)[[1]]
  if (length(r)==0) return(NULL)
  list(scaffold=r[2], start=as.numeric(r[3]), end=as.numeric(r[4]), rep=r[5], base=b)
}

pick_window <- function(smc_dir, scaffold, position, recursive=TRUE){
  smc_files <- list.files(smc_dir, pattern="\\.smc$", full.names=TRUE, recursive=recursive)
  if (length(smc_files)==0) stop("No .smc files found in: ", smc_dir)

  parsed <- lapply(smc_files, parse_outargs_name)
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed)==0) stop("No files matched expected naming: outargs_<scaffold>_<start>_<end>.<rep>.smc")

  dt <- rbindlist(lapply(seq_along(parsed), function(i){
    p <- parsed[[i]]
    data.table(file=smc_files[i], scaffold=p$scaffold, start=p$start, end=p$end, rep=p$rep, base=p$base)
  }))

  dt <- dt[scaffold == scaffold & start <= position & position <= end]
  if (nrow(dt)==0) stop("No window found containing ", scaffold, ":", position)

  dt[, width := end - start]
  setorder(dt, width, start, end)
  win <- dt[1, .(scaffold, start, end)]
  win
}

find_rep_file <- function(dir, scaffold, start, end, rep, ext=c("smc","sites"), recursive=TRUE){
  ext <- match.arg(ext)
  base <- sprintf("outargs_%s_%s_%s.%s.%s",
                  scaffold,
                  format(start, scientific=FALSE, trim=TRUE),
                  format(end, scientific=FALSE, trim=TRUE),
                  rep, ext)

  # However filenames often keep leading zeros. So do robust search by parsing again:
  files <- list.files(dir, pattern=paste0("\\.", ext, "$"), full.names=TRUE, recursive=recursive)
  if (length(files)==0) return(NULL)

  if (ext=="smc"){
    ok <- vapply(files, function(f){
      p <- parse_outargs_name(f)
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
# Read .smc and find TREE interval overlapping position
# -----------------------------
read_smc_trees <- function(smc_file){
  lines <- readLines(smc_file, warn=FALSE)
  keep <- grep("^TREE\\t", lines, value=TRUE)
  if (length(keep)==0) stop("No TREE lines found in .smc: ", smc_file)

  # TREE  <start> <end> <newick>
  parts <- tstrsplit(keep, "\t", fixed=TRUE)
  dt <- data.table(
    start=as.numeric(parts[[2]]),
    end=as.numeric(parts[[3]]),
    newick=parts[[4]]
  )
  dt
}

pick_tree_interval <- function(tree_dt, position){
  # inclusive start, exclusive end is typical for intervals
  hit <- tree_dt[start <= position & position < end]
  if (nrow(hit)==0){
    # fallback inclusive end
    hit <- tree_dt[start <= position & position <= end]
  }
  if (nrow(hit)==0) return(NULL)
  hit[1]
}

# -----------------------------
# Read .sites (ARGweaver format)
# -----------------------------
read_sites <- function(sites_file){
  lines <- readLines(sites_file, warn=FALSE)
  name_idx <- grep("^#NAMES\\s*\\t", lines)
  if (length(name_idx)==0){
    # try any case
    name_idx <- grep("^#NAMES", lines, ignore.case=TRUE)
  }
  if (length(name_idx)==0){
    stop("Sites missing NAMES line: ", sites_file)
  }
  name_line <- lines[name_idx[1]]
  names <- strsplit(name_line, "\t", fixed=TRUE)[[1]]
  names <- names[-1]
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

  data.table(pos=pos, allele_str=allele_str, n=n, names=list(names))
}

# Expand allele strings -> haplotype long table for selected positions
sites_to_hap_long <- function(names, pos_vec, allele_vec){
  n <- length(names)
  # Build matrix: rows = positions, cols = haplotypes
  # Efficient extraction: substring for each site
  alle_mat <- vapply(allele_vec, function(s) substring(s, seq_len(n), seq_len(n)), character(n))
  # alle_mat: n x m, convert
  dt <- as.data.table(alle_mat)
  dt[, hap := names]
  dt_long <- melt(dt, id.vars="hap", variable.name="site_idx", value.name="allele")
  dt_long[, position := pos_vec[as.integer(sub("^V","", site_idx))]]
  dt_long[, site_idx := NULL]
  dt_long
}

# Convert hap alleles -> individual genotypes per position
hap_long_to_genotypes <- function(hap_long){
  hap_long[, indiv := extract_indiv(hap)]
  hap_long[, copy := vapply(hap, extract_copy, integer(1))]
  # if copy not present, we still group by indiv and take first two haplotypes
  hap_long <- hap_long[order(indiv, position, copy, hap)]

  geno <- hap_long[, .(
    allele1 = allele[1],
    allele2 = if (.N>=2) allele[2] else NA_character_,
    hap1 = hap[1],
    hap2 = if (.N>=2) hap[2] else NA_character_
  ), by=.(indiv, position)]

  # genotype string
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
  # count alleles across diploid genotypes
  a <- rbind(
    geno[!is.na(allele1), .(allele=allele1)],
    geno[!is.na(allele2), .(allele=allele2)]
  )
  a[, .(Count=.N), by=allele][order(-Count)]
}

allele_counts_by_group <- function(geno, meta, group_col){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  a <- rbind(
    dt[!is.na(allele1), .(allele=allele1, group=get(group_col), position)],
    dt[!is.na(allele2), .(allele=allele2, group=get(group_col), position)]
  )
  a[, .(Count=.N), by=.(group, allele)][order(group, -Count)]
}

genotype_counts_overall <- function(geno){
  geno[!is.na(genotype), .(Count=.N), by=genotype][order(-Count)]
}

het_hom_counts <- function(geno){
  geno[!is.na(is_het), .(Count=.N), by=.(state=ifelse(is_het, "HET", "HOM"))][order(state)]
}

per_position_het_hom <- function(geno){
  geno[!is.na(is_het), .(Count=.N), by=.(position, state=ifelse(is_het, "HET", "HOM"))][order(position, state)]
}

trait_by_allele_per_position <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  # each allele observation (diploid)
  a <- rbind(
    dt[!is.na(allele1), .(position, indiv, allele=allele1, traitv)],
    dt[!is.na(allele2), .(position, indiv, allele=allele2, traitv)]
  )
  a <- a[!is.na(traitv)]
  a[, .(
    n = uniqueN(indiv),
    mean = mean(traitv),
    median = median(traitv)
  ), by=.(position, allele)][order(position, allele)]
}

trait_by_genotype_per_position <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  dt <- dt[!is.na(traitv) & !is.na(genotype)]
  dt[, .(
    n = uniqueN(indiv),
    mean = mean(traitv),
    median = median(traitv)
  ), by=.(position, genotype)][order(position, genotype)]
}

trait_by_allele_overall_unique_indiv <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  a <- rbind(
    dt[!is.na(allele1), .(indiv, allele=allele1, traitv)],
    dt[!is.na(allele2), .(indiv, allele=allele2, traitv)]
  )
  a <- a[!is.na(traitv)]
  # unique individuals carrying allele at any contributing site
  a <- unique(a, by=c("indiv","allele"))
  a[, .(n=.N, mean=mean(traitv), median=median(traitv)), by=allele][order(-n)]
}

trait_by_genotype_overall_unique_indiv <- function(geno, meta, trait){
  dt <- merge(geno, meta, by="indiv", all.x=TRUE)
  dt[, traitv := suppressWarnings(as.numeric(get(trait)))]
  dt <- dt[!is.na(traitv) & !is.na(genotype)]
  # unique individuals that ever exhibit a genotype across contributing sites
  dt <- unique(dt[, .(indiv, genotype, traitv)], by=c("indiv","genotype"))
  dt[, .(n=.N, mean=mean(traitv), median=median(traitv)), by=genotype][order(-n)]
}

# -----------------------------
# Main
# -----------------------------
args <- commandArgs(trailingOnly=TRUE)

smc_dir   <- get_arg(args, "--smc_dir")
sites_dir <- get_arg(args, "--sites_dir")
scaffold  <- get_arg(args, "--scaffold")
position  <- as.numeric(get_arg(args, "--position"))
pheno     <- get_arg(args, "--pheno")
trait     <- get_arg(args, "--trait", "C13")
rep_str   <- get_arg(args, "--replicates", "0")
outbase   <- get_arg(args, "--outbase", "tree_informative_sites_stats")
recursive <- as_bool(get_arg(args, "--recursive", "true"), TRUE)

id_col    <- get_arg(args, "--id_col", "codg")   # phenotype ID column
proc_col  <- get_arg(args, "--proc_col", "proc")
site_col  <- get_arg(args, "--site_col", "site")
mum_col   <- get_arg(args, "--mum_col",  "mum")

if (is.null(smc_dir) || is.null(sites_dir) || is.null(scaffold) || is.na(position) || is.null(pheno)) {
  cat("Usage:\n",
      "  Rscript auto_tree_informative_sites_stats_v1.R \\\n",
      "    --smc_dir=smc_files --sites_dir=sites_for_tree \\\n",
      "    --scaffold=scaffold_4 --position=983057685 \\\n",
      "    --pheno=PHENO_Charles_6_2025.txt --trait=C13 \\\n",
      "    --replicates=0,10,20,30,40,50 --outbase=tree_stats_out\n", sep="")
  quit(status=1)
}

reps <- split_csv(rep_str)
if (length(reps)==0) reps <- "0"

# output base dir
outdir0 <- file.path(outbase, paste0(scaffold, "_", format(position, scientific=FALSE, trim=TRUE)))
dir.create(outdir0, recursive=TRUE, showWarnings=FALSE)

# auto window
win <- pick_window(smc_dir, scaffold, position, recursive=recursive)
cat("AUTO MODE\n")
cat("Selected window:", win$scaffold, win$start, "-", win$end, "\n")

# phenotype load
ph <- fread(pheno, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
stopifnot(id_col %in% names(ph))
ph[, indiv := clean_id(get(id_col))]

# keep meta fields if present
meta_cols <- intersect(c(proc_col, site_col, mum_col, trait), names(ph))
meta <- ph[, c("indiv", meta_cols), with=FALSE]

# map proc/site labels (optional, keeps original numeric too)
proc_map <- c("1"="Deer Mtn","2"="Inverness River","3"="Judy Creek","4"="Swann Hills","5"="Virginia Hills")
site_map <- c("1"="JUDY","2"="VIRG","3"="SWAN","4"="TIME")

if (proc_col %in% names(meta)) {
  meta[, proc_lab := proc_map[as.character(get(proc_col))]]
}
if (site_col %in% names(meta)) {
  meta[, site_lab := site_map[as.character(get(site_col))]]
}

# process each replicate
for (rep in reps) {
  rep_out <- file.path(outdir0, paste0("rep", rep))
  dir.create(rep_out, recursive=TRUE, showWarnings=FALSE)

  smc_file <- find_rep_file(smc_dir, win$scaffold, win$start, win$end, rep, ext="smc", recursive=recursive)
  sites_file <- find_rep_file(sites_dir, win$scaffold, win$start, win$end, rep, ext="sites", recursive=recursive)

  if (is.null(smc_file)) {
    cat("Rep", rep, ": missing .smc file for this window\n")
    next
  }
  if (is.null(sites_file)) {
    cat("Rep", rep, ": missing .sites file for this window\n")
    next
  }

  cat("Rep", rep, ":\n  smc:", basename(smc_file), "\n  sites:", basename(sites_file), "\n")

  trees <- read_smc_trees(smc_file)
  hit <- pick_tree_interval(trees, position)
  if (is.null(hit)) stop("No TREE interval overlaps --position=", position, " in: ", smc_file)

  tree_start <- hit$start
  tree_end   <- hit$end

  # Read sites
  sx <- read_sites(sites_file)
  names_hap <- unlist(sx$names[[1]])
  sites_dt <- sx[, .(pos, allele_str)]

  # informative sites contributing to that local tree
  contrib <- sites_dt[pos >= tree_start & pos < tree_end]
  n_sites <- nrow(contrib)

  # write interval + counts
  interval_dt <- data.table(
    scaffold=scaffold,
    query_position=position,
    window_start=win$start, window_end=win$end,
    tree_start=tree_start, tree_end=tree_end,
    smc_file=smc_file,
    sites_file=sites_file,
    n_informative_sites=n_sites
  )
  fwrite(interval_dt, file.path(rep_out, "tree_interval_summary.tsv"), sep="\t")

  # write positions list
  fwrite(contrib[, .(position=pos)], file.path(rep_out, "informative_site_positions.tsv"), sep="\t")

  if (n_sites == 0) {
    cat("  No informative sites in this tree interval.\n")
    next
  }

  # build haplotype-long then genotypes
  hap_long <- sites_to_hap_long(names_hap, contrib$pos, contrib$allele_str)
  geno <- hap_long_to_genotypes(hap_long)

  # merge metadata (proc/site/mum/trait)
  geno_meta <- merge(geno, meta, by.x="indiv", by.y="indiv", all.x=TRUE)

  # ------------------- outputs: detailed tables -------------------
  fwrite(geno_meta[order(position, indiv)],
         file.path(rep_out, "genotypes_long_per_position.tsv"), sep="\t")

  # per-position hom/het
  fwrite(per_position_het_hom(geno),
         file.path(rep_out, "het_hom_counts_by_position.tsv"), sep="\t")

  # overall hom/het
  fwrite(het_hom_counts(geno),
         file.path(rep_out, "het_hom_counts_overall.tsv"), sep="\t")

  # genotype overall pooled across all contributing sites
  fwrite(genotype_counts_overall(geno),
         file.path(rep_out, "genotype_overall.tsv"), sep="\t")

  # allele overall pooled across all contributing sites
  fwrite(allele_counts_overall(geno),
         file.path(rep_out, "allele_overall.tsv"), sep="\t")

  # genotype counts by groups (pooled across all sites)
  if (proc_col %in% names(geno_meta)) {
    tmp <- geno_meta[!is.na(genotype), .(Count=.N), by=.(proc=get(proc_col), genotype)][order(proc, -Count)]
    fwrite(tmp, file.path(rep_out, "genotype_by_provenance.tsv"), sep="\t")

    tmp2 <- allele_counts_by_group(geno, meta, group_col=proc_col)
    fwrite(tmp2, file.path(rep_out, "allele_by_provenance.tsv"), sep="\t")
  }

  if (site_col %in% names(geno_meta)) {
    tmp <- geno_meta[!is.na(genotype), .(Count=.N), by=.(site=get(site_col), genotype)][order(site, -Count)]
    fwrite(tmp, file.path(rep_out, "genotype_by_site.tsv"), sep="\t")

    tmp2 <- allele_counts_by_group(geno, meta, group_col=site_col)
    fwrite(tmp2, file.path(rep_out, "allele_by_site.tsv"), sep="\t")
  }

  if (mum_col %in% names(geno_meta)) {
    tmp <- geno_meta[!is.na(genotype), .(Count=.N), by=.(mum=get(mum_col), genotype)][order(mum, -Count)]
    fwrite(tmp, file.path(rep_out, "genotype_by_mum.tsv"), sep="\t")

    tmp2 <- allele_counts_by_group(geno, meta, group_col=mum_col)
    fwrite(tmp2, file.path(rep_out, "allele_by_mum.tsv"), sep="\t")
  }

  # Trait summaries
  if (trait %in% names(meta)) {
    fwrite(trait_by_allele_per_position(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_allele_by_position.tsv")), sep="\t")

    fwrite(trait_by_genotype_per_position(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_genotype_by_position.tsv")), sep="\t")

    fwrite(trait_by_allele_overall_unique_indiv(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_allele_overall_unique_indiv.tsv")), sep="\t")

    fwrite(trait_by_genotype_overall_unique_indiv(geno, meta, trait),
           file.path(rep_out, paste0(trait, ".trait_by_genotype_overall_unique_indiv.tsv")), sep="\t")
  }

  cat("  Wrote TSVs to: ", rep_out, "\n", sep="")
}

cat("DONE. Output base dir: ", outdir0, "\n", sep="")

