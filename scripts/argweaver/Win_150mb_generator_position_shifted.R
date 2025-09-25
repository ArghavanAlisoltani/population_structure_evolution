setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/Win_150_1aa_1ab/")
# fixed_windows_150Mb.R
# Input: a 2-column table with header: CHROM, POS
# Output: TSV with columns: Scaffold  start  end  nsnp  length

infile <- "~/Desktop/OSU_projects/conifers/LP/vcf_v1/positions_split_poly_s100_scaffolds.tsv"  # change if needed
outfile <- "windows_150Mb.tsv"
W       <- 150000000L  # window size

# ---- read input (robust to column names) ----
suppressWarnings({
  d <- read.table(infile, header = TRUE, sep = "\t",
                  comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
})
if (!all(c("CHROM","POS") %in% names(d))) {
  cn <- names(d)
  sc_col  <- if ("Scaffold" %in% cn) "Scaffold" else cn[1]
  pos_col <- if ("Position" %in% cn) "Position" else if ("POS" %in% cn) "POS" else cn[2]
  d <- setNames(d[, c(sc_col, pos_col)], c("CHROM","POS"))
}
d$POS <- as.numeric(d$POS)
d <- d[is.finite(d$POS), ]
stopifnot(nrow(d) > 0)

# ---- make fixed windows per scaffold ----
make_windows <- function(scaf_pos) {
  L <- max(scaf_pos, na.rm = TRUE)
  k <- as.integer(ceiling(L / W))
  starts <- as.integer((0:(k-1)) * W + 1L)
  ends   <- pmin(as.integer((1:k) * W), L)
  # bin index for each SNP (1..k), inclusive on right
  bins <- pmin(pmax(ceiling(scaf_pos / W), 1L), k)
  nsnp <- tabulate(bins, nbins = k)
  data.frame(
    start  = starts,
    end    = ends,
    nsnp   = nsnp,
    length = ends - starts + 1L,
    check.names = FALSE
  )
}

by_scaf <- split(d$POS, d$CHROM)
out_list <- lapply(names(by_scaf), function(sc) {
  w <- make_windows(by_scaf[[sc]])
  cbind(Scaffold = sc, w, row.names = NULL)
})
out <- do.call(rbind, out_list)
out$ratio<-(out$nsnp/out$length)*1000
scaff1aa<-out[out$Scaffold=="scaffold_1aa",]
boxplot(scaff1aa$ratio)
tbl<-data.frame(table(out$Scaffold))
# ---- write TSV ----
write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat(sprintf("Wrote %d windows across %d scaffolds to %s\n",
            nrow(out), length(by_scaf), outfile))
