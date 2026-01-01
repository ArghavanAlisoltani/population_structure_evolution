#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(optparse)
})

# ----------------------------
# Helpers
# ----------------------------
cmd_exists <- function(x) nzchar(Sys.which(x))

stop2 <- function(...) stop(paste0(...), call. = FALSE)

split_csv <- function(x) {
  x <- gsub("\\s+", "", x)
  if (!nzchar(x)) return(character())
  unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
}

# Convert GT + REF/ALT into phased alleles
# Returns list(h1=char, h2=char) with bases (A/C/G/T/N)
gt_to_haps <- function(gt, ref, alt) {
  if (is.na(gt) || gt %in% c(".", "./.", ".|.")) return(list(h1="N", h2="N"))
  sep <- if (grepl("\\|", gt)) "|" else if (grepl("/", gt)) "/" else NA
  if (is.na(sep)) return(list(h1="N", h2="N"))
  a <- strsplit(gt, sep, fixed = TRUE)[[1]]
  if (length(a) != 2) return(list(h1="N", h2="N"))
  conv <- function(x) {
    if (x %in% c(".", "")) return("N")
    if (x == "0") return(ref)
    if (x == "1") return(alt)
    # if multi-allelic (shouldn’t happen after your filters), mark N
    return("N")
  }
  list(h1 = conv(a[1]), h2 = conv(a[2]))
}

# Ensure VCF is bgzip + tabix indexed (needed for fast random access)
ensure_bgzip_index <- function(vcf, outdir) {
  if (grepl("\\.vcf\\.gz$", vcf) && file.exists(paste0(vcf, ".tbi"))) {
    return(normalizePath(vcf))
  }
  if (!cmd_exists("bgzip") || !cmd_exists("tabix")) {
    stop2(
      "VCF is not bgzip+tabix indexed, and bgzip/tabix not found in PATH.\n",
      "Fix options:\n",
      "  1) bgzip -c input.vcf > input.vcf.gz; tabix -p vcf input.vcf.gz\n",
      "  2) Or ensure bcftools/htslib provides bgzip/tabix.\n"
    )
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  vcf_gz <- file.path(outdir, basename(vcf))
  if (!grepl("\\.gz$", vcf_gz)) vcf_gz <- paste0(vcf_gz, ".gz")

  message("Indexing VCF for fast queries:\n  ", vcf, "\n-> ", vcf_gz)
  system2("bgzip", c("-c", shQuote(vcf)), stdout = vcf_gz)
  system2("tabix", c("-p", "vcf", shQuote(vcf_gz)))

  if (!file.exists(paste0(vcf_gz, ".tbi"))) {
    stop2("Failed to create tabix index: ", vcf_gz, ".tbi")
  }
  normalizePath(vcf_gz)
}

# Read multi positions from VCF with bcftools query (fast)
read_positions_from_vcf <- function(vcf_gz, scaffold, positions, outdir) {
  if (!cmd_exists("bcftools")) stop2("bcftools not found in PATH.")

  # Create BED-like regions file: CHROM  start0  end (1-based pos -> start=pos-1)
  reg <- data.table(
    CHROM = scaffold,
    start0 = as.integer(positions) - 1L,
    end1   = as.integer(positions)
  )
  reg <- reg[start0 >= 0]
  reg_file <- file.path(outdir, "regions.bed")
  fwrite(reg, reg_file, sep = "\t", col.names = FALSE)

  # sample IDs in the same order bcftools prints genotypes
  sample_ids <- system2("bcftools", c("query", "-l", shQuote(vcf_gz)), stdout = TRUE)
  if (length(sample_ids) < 1) stop2("No samples found in VCF: ", vcf_gz)

  # Query: CHROM POS REF ALT then GTs
  fmt <- "%CHROM\t%POS\t%REF\t%ALT[\\t%GT]\\n"
  cmd <- c("query", "-R", shQuote(reg_file), "-f", fmt, shQuote(vcf_gz))
  lines <- system2("bcftools", cmd, stdout = TRUE)
  if (length(lines) < 1) stop2("No variants returned for requested positions. Check scaffold/positions.")

  # Parse into data.table
  dt <- fread(text = lines, sep = "\t", header = FALSE)
  if (ncol(dt) < 5) stop2("Unexpected bcftools query output; got ", ncol(dt), " columns.")

  setnames(dt, c("CHROM", "POS", "REF", "ALT", sample_ids))
  dt
}

# Build haplotype strings per sample across ordered positions
build_haplotypes <- function(dt_vcf) {
  # Order by POS
  setorder(dt_vcf, POS)

  sample_cols <- setdiff(names(dt_vcf), c("CHROM", "POS", "REF", "ALT"))
  if (length(sample_cols) < 1) stop2("No sample genotype columns found after parsing VCF.")

  # For each sample, iterate positions, build hap1/hap2 strings
  haps <- lapply(sample_cols, function(sid) {
    gts <- dt_vcf[[sid]]
    ref <- dt_vcf$REF
    alt <- dt_vcf$ALT

    h1 <- character(nrow(dt_vcf))
    h2 <- character(nrow(dt_vcf))
    for (i in seq_len(nrow(dt_vcf))) {
      hh <- gt_to_haps(gts[i], ref[i], alt[i])
      h1[i] <- hh$h1
      h2[i] <- hh$h2
    }
    list(sample_id = sid,
         hap1_seq  = paste0(h1, collapse = ""),
         hap2_seq  = paste0(h2, collapse = ""))
  })

  wide <- rbindlist(haps)
  wide[, zygosity := ifelse(hap1_seq == hap2_seq, "homozygous", "heterozygous")]

  # membership: heterozygotes contribute to BOTH hap groups; homozygotes once
  mem <- rbindlist(list(
    wide[, .(sample_id, zygosity, hap_group_seq = hap1_seq, which_hap="hap1")],
    wide[, .(sample_id, zygosity, hap_group_seq = hap2_seq, which_hap="hap2")]
  ), use.names = TRUE)
  mem <- unique(mem, by = c("sample_id", "hap_group_seq"))

  # map hap seq -> hap code
  hap_levels <- sort(unique(mem$hap_group_seq))
  hap_map <- data.table(
    hap_group_seq = hap_levels,
    hap_code = sprintf("H%02d", seq_along(hap_levels))
  )
  mem <- merge(mem, hap_map, by = "hap_group_seq", all.x = TRUE)

  # group counts by unique individuals
  hap_counts <- mem[, .(n_ind = uniqueN(sample_id)), by = .(hap_code, hap_group_seq)]
  setorder(hap_counts, -n_ind, hap_code)

  list(wide = wide, mem = mem, hap_map = hap_map, hap_counts = hap_counts, pos = dt_vcf$POS)
}

# Phenotype merge (ID matching: exact or digits)
merge_pheno <- function(mem, pheno_file, sep, id_col, traits) {
  ph <- fread(pheno_file, sep = sep)
  if (!(id_col %in% names(ph))) stop2("Phenotype file missing id_col: ", id_col)

  ph[, sample_id_raw := as.character(get(id_col))]
  ph[, sample_id_digits := str_extract(sample_id_raw, "\\d+")]
  ph[, sample_id := ifelse(!is.na(sample_id_digits) & nzchar(sample_id_digits), sample_id_digits, sample_id_raw)]

  keep_cols <- unique(c("sample_id", traits))
  missing_traits <- setdiff(traits, names(ph))
  if (length(missing_traits)) stop2("Phenotype file missing trait columns: ", paste(missing_traits, collapse=", "))

  ph2 <- ph[, ..keep_cols]
  # coerce traits to numeric where possible
  for (tr in traits) ph2[[tr]] <- suppressWarnings(as.numeric(ph2[[tr]]))

  merge(mem, ph2, by = "sample_id", all.x = TRUE)
}

# ----------------------------
# CLI
# ----------------------------
opt <- list(
  make_option("--vcf", type="character", help="Input VCF (.vcf or .vcf.gz). If not gz+tbi, will be bgzip+tabix indexed into outdir."),
  make_option("--pheno", type="character", help="Phenotype table (tsv)."),
  make_option("--pheno_sep", type="character", default="\t", help="Phenotype delimiter. Default: tab."),
  make_option("--id_col", type="character", default="codg", help="Phenotype ID column. Default: codg"),
  make_option("--scaffold", type="character", help="Scaffold/CHROM (e.g., scaffold_4)"),
  make_option("--positions", type="character", help="Comma-separated positions (e.g., 983057685,983057688,...)"),
  make_option("--traits", type="character", default="C13", help="Comma-separated trait columns to plot. First trait is used for the figure title/y-axis. Default: C13"),
  make_option("--outdir", type="character", default="vcf_hap_box_out", help="Output directory base."),
  make_option("--png", type="logical", default=TRUE, help="Write PNG. Default TRUE"),
  make_option("--pdf", type="logical", default=FALSE, help="Write PDF. Default FALSE"),
  make_option("--width", type="double", default=12, help="Plot width. Default 12"),
  make_option("--height", type="double", default=6, help="Plot height. Default 6"),
  make_option("--base_size", type="double", default=12, help="Theme base font size. Default 12"),
  make_option("--jitter_size", type="double", default=2.0, help="Jitter point size. Default 2.0"),
  make_option("--jitter_alpha", type="double", default=0.75, help="Jitter alpha. Default 0.75"),
  make_option("--box_alpha", type="double", default=0.35, help="Box alpha. Default 0.35"),
  make_option("--x_angle", type="double", default=45, help="X label angle. Default 45")
)

args <- parse_args(OptionParser(option_list = opt))

if (is.null(args$vcf) || is.null(args$pheno) || is.null(args$scaffold) || is.null(args$positions)) {
  stop2("Missing required args. Required: --vcf --pheno --scaffold --positions")
}

positions <- as.integer(split_csv(args$positions))
if (anyNA(positions) || length(positions) < 1) stop2("Could not parse --positions")

traits <- split_csv(args$traits)
if (length(traits) < 1) stop2("Could not parse --traits")

outdir <- file.path(args$outdir, paste0(args$scaffold, "_", min(positions), "_", max(positions)))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Main
# ----------------------------
message("Preparing VCF...")
vcf_gz <- ensure_bgzip_index(args$vcf, outdir)

message("Reading requested positions from VCF...")
dt_vcf <- read_positions_from_vcf(vcf_gz, args$scaffold, positions, outdir)

# write extracted genotype table (wide)
fwrite(dt_vcf, file.path(outdir, "vcf_selected_positions_raw.tsv"), sep = "\t")

message("Building haplotypes...")
hap <- build_haplotypes(dt_vcf)

# Outputs: hap tables
fwrite(hap$hap_map,    file.path(outdir, "hap_code_to_sequence.tsv"), sep = "\t")
fwrite(hap$hap_counts, file.path(outdir, "haplotype_counts_unique_individuals.tsv"), sep = "\t")
fwrite(hap$mem,        file.path(outdir, "sample_to_haplotype_groups.tsv"), sep = "\t")

# Merge phenotype
message("Merging phenotype...")
plot_dt <- merge_pheno(hap$mem, args$pheno, args$pheno_sep, args$id_col, traits)

# drop missing first trait for plotting
main_trait <- traits[1]
plot_dt <- plot_dt[!is.na(get(main_trait))]

# Order hap groups by n_ind (desc)
order_groups <- hap$hap_counts[order(-n_ind)]$hap_code
plot_dt[, hap_code := factor(hap_code, levels = order_groups)]

# Build labels with n
tmp <- hap$hap_counts[match(order_groups, hap$hap_counts$hap_code)]
x_lab_map <- setNames(
  paste0(tmp$hap_code, "(", tmp$hap_group_seq, ") (n=", tmp$n_ind, ")"),
  tmp$hap_code
)
leg_lab_map <- setNames(
  paste0(tmp$hap_code, "=", tmp$hap_group_seq),
  tmp$hap_code
)

# Color palette (auto-expand safely)
n_groups <- length(order_groups)
pal <- grDevices::hcl.colors(n_groups, palette = "Dark 3")
hap_cols <- setNames(pal, order_groups)

# Plot
p <- ggplot(plot_dt, aes(x = hap_code, y = .data[[main_trait]], fill = hap_code, color = hap_code)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = args$box_alpha, linewidth = 0.5) +
  geom_jitter(aes(shape = zygosity), width = 0.18, height = 0,
              size = args$jitter_size, alpha = args$jitter_alpha) +
  scale_x_discrete(labels = x_lab_map) +
  scale_fill_manual(values = hap_cols, breaks = order_groups, labels = leg_lab_map, drop = FALSE) +
  scale_color_manual(values = hap_cols, breaks = order_groups, labels = leg_lab_map, drop = FALSE) +
  scale_shape_manual(values = c("homozygous" = 16, "heterozygous" = 17)) +
  labs(
    x = "Haplotype group (sequence across selected positions; n unique individuals)",
    y = main_trait,
    title = paste0(main_trait, " by haplotype group (heterozygotes contribute to both haplotypes)"),
    subtitle = paste0(args$scaffold, " positions: ", paste(hap$pos, collapse = ", "))
  ) +
  theme_classic(base_size = args$base_size) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5),
    axis.title.x = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    axis.text.x  = element_text(angle = args$x_angle, hjust = 1, vjust = 1),
    legend.title = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Save plot data used
fwrite(plot_dt, file.path(outdir, paste0("plot_data_", main_trait, "_hapgroups.tsv")), sep = "\t")

if (isTRUE(args$png)) {
  ggsave(file.path(outdir, paste0("box_jitter_", main_trait, "_by_HAPLOTYPES_vcf.png")),
         p, width = args$width, height = args$height, dpi = 300)
}
if (isTRUE(args$pdf)) {
  ggsave(file.path(outdir, paste0("box_jitter_", main_trait, "_by_HAPLOTYPES_vcf.pdf")),
         p, width = args$width, height = args$height)
}

message("DONE. Output dir:\n  ", normalizePath(outdir))

