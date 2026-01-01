############################################################
## Genotype counts per locus from VCF (AA/AB/BB/missing)
## Output TSV with per-variant genotype counts + het count
############################################################
setwd("~/Desktop/OSU_projects/conifers/LP/Files_for_Katherine/GBS_SNPs")

# install.packages("vcfR")
suppressPackageStartupMessages({
  library(vcfR)
})

# -------- user inputs ----------
vcf_file <- "Imputed_whole_panel_Esteban_Soms_shared.vcf.gz"  # .vcf or .vcf.gz
out_tsv  <- "genotype_counts_per_locus.tsv"

# Optional: analyze only one scaffold (set to NULL to keep all)
only_chrom <- NULL  # e.g. "scaffold_4"

# -------- load VCF ----------
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Optional filter to one chromosome
if (!is.null(only_chrom)) {
  fix <- getFIX(vcf)
  keep <- fix[, "CHROM"] == only_chrom
  vcf <- vcf[keep, ]
}

# -------- extract GT matrix ----------
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)  # variants x samples

# helper: recode GT to AA/AB/BB/NA using REF=0, ALT=1
recode_gt <- function(x) {
  # normalize phased/unphased
  x <- gsub("\\|", "/", x)
  # keep only diploid 0/0,0/1,1/0,1/1; everything else -> NA
  out <- rep(NA_character_, length(x))
  out[x == "0/0"] <- "AA"
  out[x == "0/1" | x == "1/0"] <- "AB"
  out[x == "1/1"] <- "BB"
  out
}

# Apply per row (variant); keep memory reasonable
gt_rec <- t(apply(gt, 1, recode_gt))  # variants x samples

# -------- counts per variant ----------
count_levels <- function(v) {
  tab <- table(factor(v, levels = c("AA","AB","BB")), useNA = "ifany")
  c(
    n_AA = unname(tab["AA"]),
    n_AB = unname(tab["AB"]),   # heterozygotes
    n_BB = unname(tab["BB"]),
    n_missing = sum(is.na(v)),
    n_called = sum(!is.na(v))
  )
}

counts <- t(apply(gt_rec, 1, count_levels))

# -------- variant metadata ----------
fix <- getFIX(vcf)
# FIX columns: CHROM POS ID REF ALT QUAL FILTER INFO
res <- data.frame(
  CHROM = fix[, "CHROM"],
  POS   = as.integer(fix[, "POS"]),
  ID    = fix[, "ID"],
  REF   = fix[, "REF"],
  ALT   = fix[, "ALT"],
  stringsAsFactors = FALSE
)

# Add counts
res$n_AA      <- counts[, "n_AA"]
res$n_AB      <- counts[, "n_AB"]      # heterozygotes (separate column)
res$n_BB      <- counts[, "n_BB"]
res$n_missing <- counts[, "n_missing"]
res$n_called  <- counts[, "n_called"]

# sanity: expected sample count
res$n_total_samples <- ncol(gt)

# -------- write TSV ----------
write.table(
  res, file = out_tsv, sep = "\t",
  quote = FALSE, row.names = FALSE
)

message("Wrote: ", out_tsv)
