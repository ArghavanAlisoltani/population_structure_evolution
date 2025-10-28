# install.packages(c("tidyverse"))
library(tidyverse)

fisher_snps_tmrca <- function(
    snp_file,
    outdir = "fisher_snps_out",
    min_class_n = 25  # only test TE classes seen in >= this many SNPs (across top+bkg)
){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- load & basic cleanup ----------
  df <- read_tsv(snp_file, show_col_types = FALSE)
  
  # keep only SNPs that have a TMRCA label we care about
  df <- df %>%
    filter(tmrca_group %in% c("top","background")) %>%
    mutate(
      tmrca_group = factor(tmrca_group, levels = c("background","top")),
      in_TE = case_when(
        is.logical(in_TE) ~ in_TE,
        is.na(in_TE) ~ FALSE,
        TRUE ~ in_TE %in% c(TRUE, "TRUE", "true", "T", "1")
      ),
      mrna_present = !is.na(mrna_id) & nzchar(trimws(mrna_id)),
      TE_class = if_else(is.na(TE_class), "", TE_class)
    )
  
  # handy helper
  safe_fisher <- function(mat2x2){
    ft <- fisher.test(mat2x2)
    tibble(
      odds_ratio = if (!is.null(ft$estimate)) unname(ft$estimate) else NA_real_,
      ci_lo = if (!is.null(ft$conf.int)) ft$conf.int[1] else NA_real_,
      ci_hi = if (!is.null(ft$conf.int)) ft$conf.int[2] else NA_real_,
      p_value = ft$p.value
    )
  }
  
  # ---------- 1) Overall: in_TE ~ tmrca_group ----------
  # 2x2: rows=tmrca_group (bkg, top), cols=in_TE (FALSE, TRUE)
  tab_te <- table(df$tmrca_group, factor(df$in_TE, levels=c(FALSE, TRUE)))
  res_te <- safe_fisher(tab_te) %>%
    mutate(
      test = "in_TE_vs_tmrca_group",
      top_TE       = unname(tab_te["top","TRUE"]),
      top_nonTE    = unname(tab_te["top","FALSE"]),
      bkg_TE       = unname(tab_te["background","TRUE"]),
      bkg_nonTE    = unname(tab_te["background","FALSE"])
    ) %>%
    relocate(test)
  
  # ---------- 2) Overall: mRNA presence ~ tmrca_group ----------
  tab_mrna <- table(df$tmrca_group, factor(df$mrna_present, levels=c(FALSE, TRUE)))
  res_mrna <- safe_fisher(tab_mrna) %>%
    mutate(
      test = "mRNA_present_vs_tmrca_group",
      top_mrna      = unname(tab_mrna["top","TRUE"]),
      top_nonmrna   = unname(tab_mrna["top","FALSE"]),
      bkg_mrna      = unname(tab_mrna["background","TRUE"]),
      bkg_nonmrna   = unname(tab_mrna["background","FALSE"])
    ) %>%
    relocate(test)
  
  overall_out <- bind_rows(res_te, res_mrna)
  write_tsv(overall_out, file.path(outdir, "fisher_overall_te_mrna.tsv"))
  
  # ---------- quick proportion plots (overall) ----------
  # TE %
  prop_te <- df %>%
    group_by(tmrca_group) %>%
    summarise(n = n(), k = sum(in_TE), .groups="drop") %>%
    rowwise() %>%
    mutate(
      p = k/n,
      # exact binomial CI per group
      ci_lo = binom.test(k, n)$conf.int[1],
      ci_hi = binom.test(k, n)$conf.int[2]
    ) %>% ungroup()
  
  g_te <- ggplot(prop_te, aes(x = tmrca_group, y = p)) +
    geom_col(width = 0.55, fill = "grey60") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = "% SNPs in TE", title = "TE overlap by TMRCA group") +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, "prop_inTE_by_tmrca_group.png"), g_te, width = 5, height = 4, dpi = 300)
  
  # mRNA %
  prop_mrna <- df %>%
    group_by(tmrca_group) %>%
    summarise(n = n(), k = sum(mrna_present), .groups="drop") %>%
    rowwise() %>%
    mutate(
      p = k/n,
      ci_lo = binom.test(k, n)$conf.int[1],
      ci_hi = binom.test(k, n)$conf.int[2]
    ) %>% ungroup()
  
  g_mrna <- ggplot(prop_mrna, aes(x = tmrca_group, y = p)) +
    geom_col(width = 0.55, fill = "grey60") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = "% SNPs overlapping mRNA", title = "mRNA overlap by TMRCA group") +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, "prop_mRNA_by_tmrca_group.png"), g_mrna, width = 5, height = 4, dpi = 300)
  
  # ---------- 3) Per-TE-class Fisher: class present ~ tmrca_group ----------
  # If a SNP has multiple classes (e.g., 'Gypsy;Copia'), they were de-duplicated upstream;
  # here we split to long format so each class is tested.
  te_long <- df %>%
    filter(TE_class != "") %>%
    separate_rows(TE_class, sep = ";", convert = FALSE) %>%
    mutate(TE_class = trimws(TE_class))
  
  totals <- df %>% count(tmrca_group, name = "N_group")
  
  class_counts <- te_long %>%
    count(TE_class, tmrca_group, name = "present") %>%
    complete(TE_class, tmrca_group, fill = list(present = 0)) %>%
    left_join(totals, by = "tmrca_group") %>%
    mutate(absent = N_group - present)
  
  # keep only classes with enough signal
  keep_classes <- class_counts %>%
    group_by(TE_class) %>%
    summarise(total_present = sum(present), .groups="drop") %>%
    filter(total_present >= min_class_n) %>%
    pull(TE_class)
  
  class_results <- map_dfr(keep_classes, function(cls){
    sub <- class_counts %>% filter(TE_class == cls)
    # build 2x2: rows=group (bkg, top), cols=present/absent
    mat <- matrix(
      c(
        sub$present[sub$tmrca_group=="background"],
        sub$absent [sub$tmrca_group=="background"],
        sub$present[sub$tmrca_group=="top"],
        sub$absent [sub$tmrca_group=="top"]
      ),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("background","top"), c("present","absent"))
    )
    ft <- fisher.test(mat)
    tibble(
      TE_class = cls,
      top_present = mat["top","present"],     top_absent = mat["top","absent"],
      bkg_present = mat["background","present"], bkg_absent = mat["background","absent"],
      odds_ratio = if (!is.null(ft$estimate)) unname(ft$estimate) else NA_real_,
      ci_lo = if (!is.null(ft$conf.int)) ft$conf.int[1] else NA_real_,
      ci_hi = if (!is.null(ft$conf.int)) ft$conf.int[2] else NA_real_,
      p_value = ft$p.value
    )
  }) %>%
    mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  write_tsv(class_results, file.path(outdir, "fisher_by_TE_class.tsv"))
  
  # ---------- forest plot for TE classes (OR with 95% CI) ----------
  if (nrow(class_results) > 0) {
    # clip extreme ORs for log scale prettiness
    plot_df <- class_results %>%
      mutate(
        TE_class = fct_reorder(TE_class, odds_ratio),
        logOR = log2(odds_ratio),
        logLO = log2(ci_lo),
        logHI = log2(ci_hi)
      )
    
    g_forest <- ggplot(plot_df, aes(x = logOR, y = TE_class)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_errorbarh(aes(xmin = logLO, xmax = logHI), height = 0.15) +
      geom_point(size = 2) +
      labs(
        x = "log2(odds ratio)  (right = enrichment in TOP)",
        y = NULL,
        title = "TE-class enrichment in TOP vs BACKGROUND (Fisher exact)"
      ) +
      theme_bw(base_size = 12)
    ggsave(file.path(outdir, "forest_TE_classes_log2OR.png"), g_forest, width = 7, height = max(4, 0.3*nrow(plot_df)+2), dpi = 300)
  }
  
  message("Done. Files written to: ", normalizePath(outdir))
}

# ---------- run ----------
# fisher_snps_tmrca("snps_annotated_tmrca_te_mrna.tsv", outdir = "fisher_snps_out", min_class_n = 25)
