#!/usr/bin/env Rscript

# Example usage:
# Rscript sites_tree_heatmap_with_pheno_v4.R \
#   --sites=sites_for_tree/outargs_scaffold_4_900000001_1050000000.50.sites \
#   --pheno=PHENO_Charles_6_2025.txt \
#   --positions=983057685,983057700 --trait=C13 \
#   --outprefix=scaffold4_pos983057685_heatmap --tree_scale=distance --hide_labels=true
suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(data.table)
})

# ---------- helpers ----------
get_arg <- function(args, flag, default=NULL){
  hit <- grep(paste0("^", flag, "="), args, value=TRUE)
  if(length(hit)==0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
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

# ---------- .sites ----------
read_sites <- function(sites_file) {
  lines <- readLines(sites_file)
  name_line <- lines[grepl("^#NAMES\\t", lines)][1]
  if (is.na(name_line)) stop("No #NAMES line found in .sites file.")
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

sites_to_alignment <- function(names, allele_str) {
  n <- length(names)
  seqs <- vapply(seq_len(n), function(i) paste0(substr(allele_str, i, i), collapse=""), character(1))
  names(seqs) <- names
  ape::as.DNAbin(strsplit(seqs, ""))
}

extract_genotypes <- function(names, pos, allele_str, positions_numeric) {
  n <- length(names)
  idx <- match(positions_numeric, pos)
  if (anyNA(idx)) stop("Requested positions not found: ", paste(positions_numeric[is.na(idx)], collapse=", "))

  geno <- data.table(label = names)
  for (k in seq_along(idx)) {
    s <- allele_str[idx[k]]
    geno[[as.character(positions_numeric[k])]] <- substring(s, seq_len(n), seq_len(n))
  }
  geno
}

# ---------- args ----------
args <- commandArgs(trailingOnly=TRUE)
sites_file  <- get_arg(args, "--sites")
pheno_file  <- get_arg(args, "--pheno")
pos_string  <- get_arg(args, "--positions")
trait_name  <- get_arg(args, "--trait", "C13")
outprefix   <- get_arg(args, "--outprefix", "sites_tree_pheno")
hide_labels <- tolower(get_arg(args, "--hide_labels", "true")) %in% c("true","t","1","yes","y")
tree_scale  <- tolower(get_arg(args, "--tree_scale", "distance"))  # distance|cladogram

if (is.null(sites_file) || is.null(pheno_file) || is.null(pos_string)) {
  cat("Usage:\n",
      "  Rscript sites_tree_heatmap_with_pheno_v4.R \\\n",
      "    --sites=FILE.sites --pheno=PHENO.txt --positions=983057685,983057700 \\\n",
      "    --trait=C13 --outprefix=out --tree_scale=distance --hide_labels=true\n", sep="")
  quit(status=1)
}
positions_numeric <- as.integer(strsplit(pos_string, ",")[[1]])
if (any(is.na(positions_numeric))) stop("Bad --positions (must be comma-separated integers).")

# ---------- 1) Tree from ONE .sites ----------
x <- read_sites(sites_file)
aln <- sites_to_alignment(x$names, x$allele_str)

d <- dist.dna(aln, model="JC69", pairwise.deletion=TRUE, as.matrix=FALSE)
tree <- bionj(d)
tree <- ladderize(tree)
write.tree(tree, file=paste0(outprefix, ".nwk"))

# ---------- 2) Alleles (long) ----------
geno_wide <- extract_genotypes(x$names, x$pos, x$allele_str, positions_numeric)
geno_long <- melt(geno_wide, id.vars="label", variable.name="position", value.name="allele")
geno_long[, position := factor(position, levels=as.character(positions_numeric))]
geno_long[, label_clean := clean_id(label)]
geno_long[, label_digits := extract_digits(label)]

# ---------- 3) Phenotype join (robust) ----------
ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE)
stopifnot("codg" %in% names(ph), "proc" %in% names(ph), "site" %in% names(ph), "mum" %in% names(ph))
if (!(trait_name %in% names(ph))) stop("Trait not found in pheno: ", trait_name)

ph[, codg := clean_id(codg)]
ph[, codg_digits := extract_digits(codg)]

tips <- data.table(label = tree$tip.label)
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

unmatched <- anno[is.na(proc) & is.na(site) & is.na(mum) & is.na(trait), label]
writeLines(unmatched, paste0(outprefix, ".unmatched_tree_tips.txt"))
cat("  unmatched tips written to:", paste0(outprefix, ".unmatched_tree_tips.txt"), "\n")

proc_map <- c("1"="Deer Mtn","2"="Inverness River","3"="Judy Creek","4"="Swann Hills","5"="Virginia Hills")
site_map <- c("1"="JUDY","2"="VIRG","3"="SWAN","4"="TIME")

anno[, proc_lab := factor(proc_map[as.character(proc)],
                          levels=c("Deer Mtn","Inverness River","Judy Creek","Swann Hills","Virginia Hills"))]
anno[, site_lab := factor(site_map[as.character(site)],
                          levels=c("Judy Creek","Virginia Hills","Swann Hills","TIME"))]

mum_levels <- unique(na.omit(anno$mum))
use_mum_numeric <- length(mum_levels) > 25
if (use_mum_numeric) anno[, mum_num := suppressWarnings(as.numeric(mum))] else anno[, mum_fac := factor(mum)]

# ---------- 4) Plot ----------
branch_mode <- if (tree_scale == "cladogram") "none" else "branch.length"
p <- ggtree(tree, branch.length=branch_mode)
if (!hide_labels) p <- p + geom_tiplab(size=2.6)

offset0 <- 0.3
w_bar   <- 0.02
gap     <- 0.01



# SITE
# p <- p +
#   ggnewscale::new_scale_fill() +
#   geom_fruit(data=anno, geom=geom_tile,
#              mapping=aes(y=label, x=1, fill=site_lab),
#              offset=offset0 + (w_bar+gap)*1, pwidth=w_bar) +
#   scale_fill_brewer(palette="Set3", na.value="grey92", name="Site")

# MUM
# if (use_mum_numeric) {
#   p <- p +
#     ggnewscale::new_scale_fill() +
#     geom_fruit(data=anno, geom=geom_tile,
#                mapping=aes(y=label, x=1, fill=mum_num),
#                offset=offset0 + (w_bar+gap)*2, pwidth=w_bar) +
#     scale_fill_viridis_c(na.value="grey92", name="Mum")
# } else {
#   p <- p +
#     ggnewscale::new_scale_fill() +
#     geom_fruit(data=anno, geom=geom_tile,
#                mapping=aes(y=label, x=1, fill=mum_fac),
#                offset=offset0 + (w_bar+gap)*2, pwidth=w_bar) +
#     scale_fill_hue(na.value="grey92", name="Mum")
# }


# ALLELES (FIXED JOIN: ensure there is exactly ONE 'label' column)
geno_key <- copy(geno_long)
geno_key[, key := if (!use_digits) label_clean else label_digits]
geno_key <- geno_key[, .(key, position, allele)]                       # DROP label before merge
geno_key <- merge(geno_key, tips[, .(key, label)], by="key", all.x=TRUE, sort=FALSE)
geno_key <- geno_key[!is.na(label)]

allele_cols <- c("A"="#1b9e77","C"="#7570b3","G"="#d95f02","T"="#e7298a","N"="#bdbdbd","-"="#252525")

p <- p +
  ggnewscale::new_scale_fill() +
  geom_fruit(data=geno_key, geom=geom_tile,
             mapping=aes(y=label, x=position, fill=allele),
             offset=5.1, pwidth=0.5, color="grey85") +
  scale_fill_manual(values=allele_cols, na.value="grey92", name="Allele") +
  theme(
    legend.position="right",
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
  )

# PROC
p <- p +
  ggnewscale::new_scale_fill() +
  geom_fruit(data=anno, geom=geom_tile,
             mapping=aes(y=label, x=1, fill=proc_lab),
             offset=0.002, pwidth=26.0000001) +
  scale_fill_brewer(palette="Set2", na.value="grey92", name="Provenances")

# TRAIT
p <- p +
  ggnewscale::new_scale_fill() +
  geom_fruit(data=anno, geom=geom_tile,
             mapping=aes(y=label, x=1, fill=trait),
             offset=0.003, pwidth=26) +
  scale_fill_viridis_c(na.value="grey92", name=trait_name)

#ggsave(paste0(outprefix, ".png"), p, width=10, height=12, dpi=300)
ggsave(paste0(outprefix, ".pdf"), p, width=8, height=16)

cat("Wrote:\n",
    "  ", paste0(outprefix, ".nwk"), "\n",
    "  ", paste0(outprefix, ".png"), "\n",
    "  ", paste0(outprefix, ".pdf"), "\n", sep="")

