# install.packages("tidyverse")
library(tidyverse)

# site_tsv: your two-column file like shown (first col = FORMAT string starting with GT, second col = sample name)
# meta_tsv: optional TSV with columns sample, group (e.g., provenance)
plot_AF_for_site <- function(site_tsv,
                             meta_tsv = NULL,
                             outprefix = "site_AF",
                             ref_allele = "A",
                             alt_allele = "C") {

  # 1) read and keep only lines that look like genotype records
  df <- read_tsv(site_tsv, col_names = c("fmt","sample"), show_col_types = FALSE)
  df <- df %>% filter(str_detect(fmt, "^([0-9.][/|][0-9.]):"))

  # 2) extract GT (first token before ':') and map to ALT allele count (diploid)
  df <- df %>%
    mutate(GT = sub(":.*", "", fmt)) %>%
    mutate(alt_cnt = case_when(
      GT %in% c("0/0","0|0") ~ 0L,
      GT %in% c("0/1","1/0","0|1","1|0") ~ 1L,
      GT %in% c("1/1","1|1") ~ 2L,
      TRUE ~ NA_integer_
    )) %>%
    mutate(called = !is.na(alt_cnt))

  # 3) attach groups (if provided)
  if (!is.null(meta_tsv)) {
    meta <- read_tsv(meta_tsv, show_col_types = FALSE) %>%
      rename(sample = 1, group = 2) %>% # or set exact names if they differ
      mutate(group = as.character(group))
    df <- df %>% left_join(meta, by = "sample")
  }
  df <- df %>% mutate(group = if_else(is.na(group), "All", group))

  # 4) per-group AF & CI (ignore missing ./.)
  af <- df %>%
    filter(called) %>%
    group_by(group) %>%
    summarise(
      n_samples = n(),
      alt_sum   = sum(alt_cnt),
      af        = alt_sum / (2*n_samples),
      # 95% exact binomial CI on ALT allele count out of 2N chromosomes
      ci_lo     = binom.test(alt_sum, 2*n_samples)$conf.int[1],
      ci_hi     = binom.test(alt_sum, 2*n_samples)$conf.int[2],
      .groups = "drop"
    ) %>%
    arrange(desc(af))

  write_tsv(af, paste0(outprefix, ".allele_freq.tsv"))

  # 5) genotype composition per group (proportions)
  geno_comp <- df %>%
    filter(GT %in% c("0/0","0/1","1/1","0|0","0|1","1|0","1|1")) %>%
    mutate(GT = recode(GT, "0|0"="0/0","0|1"="0/1","1|0"="0/1","1|1"="1/1")) %>%
    count(group, GT) %>%
    group_by(group) %>%
    mutate(prop = n/sum(n)) %>%
    ungroup()

  write_tsv(geno_comp, paste0(outprefix, ".genotype_props.tsv"))

  # 6) plots
  g1 <- ggplot(af, aes(x = reorder(group, af), y = af)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = paste0("ALT freq (", alt_allele, ")"),
         title = "ALT allele frequency by group (95% CI)",
         subtitle = "Missing genotypes (./.) excluded") +
    theme_bw(base_size = 12)
  ggsave(paste0(outprefix, ".AF_by_group.png"), g1, width = 6, height = max(3, 0.35*nrow(af)+1), dpi = 300)

  g2 <- ggplot(geno_comp, aes(x = group, y = prop, fill = GT)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = NULL, y = "Genotype proportion", fill = "GT",
         title = "Genotype composition by group") +
    theme_bw(base_size = 12)
  ggsave(paste0(outprefix, ".genotype_composition.png"), g2, width = 6, height = max(3, 0.35*length(unique(geno_comp$group))+1), dpi = 300)

  invisible(list(af = af, geno = geno_comp))
}

# Example:
# plot_AF_for_site("site_scaf6_72516188.tsv",
#                  meta_tsv = "samples_to_provenance.tsv",   # two columns: sample \t group
#                  outprefix = "scaf6_72516188_AF",
#                  ref_allele = "A", alt_allele = "C")
