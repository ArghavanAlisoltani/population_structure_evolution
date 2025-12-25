#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

get_arg <- function(args, flag, default=NULL){
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if(length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}
as_bool <- function(x, default=FALSE){
  if (is.null(x) || length(x)==0) return(default)
  tolower(trimws(as.character(x))) %in% c("true","t","1","yes","y")
}

# ---------- read ARGweaver .sites ----------
read_sites_one <- function(sites_file){
  lines <- readLines(sites_file, warn=FALSE)

  name_line <- lines[grepl("^#NAMES\\t", lines)]
  if (length(name_line) == 0) name_line <- lines[grepl("^#NAMES", lines)]
  if (length(name_line) == 0) stop("No #NAMES line in: ", sites_file)

  names <- strsplit(name_line[1], "\t", fixed=TRUE)[[1]][-1]
  n <- length(names)
  if (n < 2) stop("Parsed <2 haplotypes from #NAMES in: ", sites_file)

  site_lines <- lines[!grepl("^#", lines)]
  site_lines <- site_lines[nchar(site_lines) > 0]

  allele_str <- sub("^[0-9]+\\t", "", site_lines)
  bad <- which(nchar(allele_str) != n)
  if (length(bad) > 0){
    stop("Allele string length != n in ", sites_file,
         ". First bad line: ", bad[1],
         " (expected ", n, ", got ", nchar(allele_str[bad[1]]), ")")
  }

  list(names=names, allele_str=allele_str)
}

build_seqs_from_allele_str <- function(names, allele_str){
  n <- length(names)
  seqs <- vapply(seq_len(n), function(i){
    paste0(substr(allele_str, i, i), collapse="")
  }, character(1))
  names(seqs) <- names
  seqs
}

write_fasta <- function(seqs, file){
  con <- file(file, "w")
  on.exit(close(con), add=TRUE)
  for (nm in names(seqs)){
    writeLines(paste0(">", nm), con)
    writeLines(seqs[[nm]], con)
  }
}

write_phylip <- function(seqs, file){
  n <- length(seqs)
  L <- nchar(seqs[[1]])
  if (!all(nchar(seqs) == L)) stop("Sequences not equal length; cannot write PHYLIP.")
  con <- file(file, "w")
  on.exit(close(con), add=TRUE)
  writeLines(sprintf("%d %d", n, L), con)
  for (nm in names(seqs)){
    writeLines(paste(nm, seqs[[nm]]), con)
  }
}

# ---------- parse filename ----------
# expected:
# outargs_scaffold_1a_000000001_150000000.0.sites
# outargs_scaffold_2_000000001_150000000.0.sites
parse_sites_name <- function(path){
  b <- basename(path)
  m <- regexec("^outargs_scaffold_([^_]+)_([0-9]+)_([0-9]+)\\.([0-9]+)\\.sites$", b)
  r <- regmatches(b, m)[[1]]
  if (length(r)==0) return(NULL)
  scaff <- r[2]         # e.g. "1a" or "2" or "1b"
  start <- as.numeric(r[3])
  end   <- as.numeric(r[4])
  rep   <- as.integer(r[5])

  # split scaffold into numeric + optional letters
  m2 <- regexec("^([0-9]+)([A-Za-z]*)$", scaff)
  r2 <- regmatches(scaff, m2)[[1]]
  if (length(r2)==0) return(NULL)
  scaff_num <- as.integer(r2[2])
  scaff_suf <- tolower(r2[3])  # "" or "a" or "b"

  list(file=path, base=b, scaff=scaff, scaff_num=scaff_num, scaff_suf=scaff_suf,
       start=start, end=end, rep=rep)
}

# suffix ordering:
# a < b < c < ... < ""   (none last)  OR  "" first
suffix_rank <- function(suf, none_last=TRUE){
  suf <- tolower(suf)
  if (suf == "") return(if (none_last) 999L else -1L)
  # map a->1, b->2, ...
  chars <- strsplit(suf, "")[[1]]
  # handle multi-letter just in case
  val <- 0L
  for (ch in chars){
    if (ch < "a" || ch > "z") next
    val <- val * 27L + (utf8ToInt(ch) - utf8ToInt("a") + 1L)
  }
  val
}

# ---------------- main ----------------
args <- commandArgs(trailingOnly=TRUE)

sites_dir    <- get_arg(args, "--sites_dir")
outprefix    <- get_arg(args, "--outprefix", "ALL_rep0_concat")
recursive    <- as_bool(get_arg(args, "--recursive", "true"), TRUE)
none_last    <- as_bool(get_arg(args, "--none_last", "true"), TRUE)  # controls scaffold_1 vs 1a/1b
write_fa     <- as_bool(get_arg(args, "--fasta", "true"), TRUE)
write_phy    <- as_bool(get_arg(args, "--phylip", "true"), TRUE)

if (is.null(sites_dir)){
  cat("Usage:\n",
      "  Rscript concat_sites_rep0_scaffold_order.R \\\n",
      "    --sites_dir=sites_for_tree \\\n",
      "    --outprefix=ALL_rep0_concat \\\n",
      "    --recursive=true --none_last=true --fasta=true --phylip=true\n", sep="")
  quit(status=1)
}

# Only rep0 sites
files <- list.files(sites_dir, pattern="\\.0\\.sites$", full.names=TRUE, recursive=recursive)
if (length(files)==0) stop("No replicate-0 .sites files found (*.0.sites) in: ", sites_dir)

parsed <- lapply(files, parse_sites_name)
parsed <- parsed[!vapply(parsed, is.null, logical(1))]
if (length(parsed)==0) stop("No files matched expected naming outargs_scaffold_<ID>_<start>_<end>.0.sites")

dt <- rbindlist(lapply(parsed, as.data.table), fill=TRUE)
dt[, suf_rank := vapply(scaff_suf, suffix_rank, integer(1), none_last=none_last)]

# Sort: scaffold number -> suffix -> window start -> window end
setorder(dt, scaff_num, suf_rank, start, end)

cat("Files in concatenation order (first 12):\n")
print(dt[1:min(12,.N), .(base, scaff, start, end)])

# Read first file to set haplotype order
x1 <- read_sites_one(dt$file[1])
ref_names <- x1$names
ref_set <- sort(ref_names)

# init concat
seqs_concat <- setNames(rep("", length(ref_names)), ref_names)

for (f in dt$file){
  x <- read_sites_one(f)
  if (!identical(sort(x$names), ref_set)){
    stop("Haplotype set differs in file: ", f)
  }
  seqs <- build_seqs_from_allele_str(x$names, x$allele_str)

  # reorder if needed
  if (!identical(x$names, ref_names)){
    seqs <- seqs[ref_names]
  }

  for (nm in ref_names){
    seqs_concat[[nm]] <- paste0(seqs_concat[[nm]], seqs[[nm]])
  }
  cat("Appended ", basename(f), " (", nchar(seqs[[1]]), " sites)\n", sep="")
}

L <- nchar(seqs_concat[[1]])
cat("Final alignment length:", L, " sites\n")

if (write_fa){
  fa <- paste0(outprefix, ".fa")
  write_fasta(seqs_concat, fa)
  cat("Wrote:", fa, "\n")
}
if (write_phy){
  phy <- paste0(outprefix, ".phy")
  write_phylip(seqs_concat, phy)
  cat("Wrote:", phy, "\n")
}

# also write the file order table for reproducibility
fwrite(dt[, .(base, scaff, scaff_num, scaff_suf, start, end, rep)],
       paste0(outprefix, ".files_in_order.tsv"), sep="\t")
cat("Wrote:", paste0(outprefix, ".files_in_order.tsv"), "\n")

