################################################################################
#   Figure 6 — Panel E: Protein-Phenotype Heatmap
#   Top hub proteins (baseline expression) correlated with Age, VL, LBM traits
################################################################################
#
# STAT AUDIT — Panel E (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Pearson r (Age, Baseline VL): computed via cor.test(). PASS.
#    ADDED: 95% CI from cor.test()$conf.int stored in ci_lo/ci_hi columns.
# 2. Partial r (Delta VL, Delta LBM): ppcor::pcor.test() controlling for
#    age group. ADDED: Fisher z-transform CI (k=1 covariate).
# 3. BH correction: applied across ALL protein-trait pairs in a single
#    p.adjust() call (n_top * 4 tests). PASS.
# 4. Stars (*/**/***)  based on BH-adjusted p-values (padj). PASS.
# 5. Note: heatmap tile text shows point estimate + stars; exported CSV
#    now includes ci_lo and ci_hi columns for all correlations.
# ---------------------------------------------------------------------------
#

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Select top 15 proteins per top-3 module by kME * |GS| score ----
top3_colors <- gsub("^ME", "", top3)

top_prots <- map_dfr(top3_colors, function(mod) {
  kme_c     <- paste0("kME", mod)
  mod_prots <- module_df$uniprot_id[module_df$module_color == mod]
  mod_prots <- intersect(mod_prots, rownames(kME_all))
  mod_prots <- intersect(mod_prots, rownames(gs_vl))
  tibble(
    uniprot_id = mod_prots,
    kME        = kME_all[mod_prots, kme_c],
    GS         = abs(gs_vl[mod_prots, 1]),
    module     = mod,
    score      = kME * GS
  ) %>%
    arrange(desc(score)) %>%
    head(15)
}) %>%
  left_join(module_df %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

# ---- Build correlation + p-value + CI matrices ----
n_top   <- nrow(top_prots)
n_subj  <- length(common_subj)
traits  <- c("Age", "Baseline VL", "Delta VL", "Delta LBM")
heat_mat  <- matrix(NA_real_, nrow = n_top, ncol = length(traits),
                    dimnames = list(top_prots$gene, traits))
pval_heat <- heat_mat
ci_lo_mat <- heat_mat
ci_hi_mat <- heat_mat

age_num_vec <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0
)
base_vl_vec <- pheno_wide$VL_Pre[match(common_subj, pheno_wide$subject_key)]

for (i in seq_len(n_top)) {
  prot      <- top_prots$uniprot_id[i]
  base_vals <- pre_expr[common_subj, prot]

  # Age — Pearson r with CI from cor.test()
  ct <- tryCatch(cor.test(base_vals, age_num_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Age"]  <- ct$estimate
    pval_heat[i, "Age"] <- ct$p.value
    ci_lo_mat[i, "Age"] <- ct$conf.int[1]
    ci_hi_mat[i, "Age"] <- ct$conf.int[2]
  }

  # Baseline VL — Pearson r with CI from cor.test()
  ct <- tryCatch(cor.test(base_vals, base_vl_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Baseline VL"]  <- ct$estimate
    pval_heat[i, "Baseline VL"] <- ct$p.value
    ci_lo_mat[i, "Baseline VL"] <- ct$conf.int[1]
    ci_hi_mat[i, "Baseline VL"] <- ct$conf.int[2]
  }

  # Delta VL — partial correlation controlling for age group
  pc <- tryCatch({
    res <- ppcor::pcor.test(base_vals, delta_vl_vec, age_num_vec, method = "pearson")
    ci  <- fisher_z_ci(res$estimate, n_subj, k = 1)
    list(estimate = res$estimate, p.value = res$p.value,
         ci_lo = ci[["lo"]], ci_hi = ci[["hi"]])
  }, error = function(e) NULL)
  if (!is.null(pc)) {
    heat_mat[i, "Delta VL"]  <- pc$estimate
    pval_heat[i, "Delta VL"] <- pc$p.value
    ci_lo_mat[i, "Delta VL"] <- pc$ci_lo
    ci_hi_mat[i, "Delta VL"] <- pc$ci_hi
  }

  # Delta LBM — partial correlation controlling for age group
  pc <- tryCatch({
    res <- ppcor::pcor.test(base_vals, delta_lbm_vec, age_num_vec, method = "pearson")
    ci  <- fisher_z_ci(res$estimate, n_subj, k = 1)
    list(estimate = res$estimate, p.value = res$p.value,
         ci_lo = ci[["lo"]], ci_hi = ci[["hi"]])
  }, error = function(e) NULL)
  if (!is.null(pc)) {
    heat_mat[i, "Delta LBM"]  <- pc$estimate
    pval_heat[i, "Delta LBM"] <- pc$p.value
    ci_lo_mat[i, "Delta LBM"] <- pc$ci_lo
    ci_hi_mat[i, "Delta LBM"] <- pc$ci_hi
  }
}

# ---- Assemble long-format heatmap data ----
# Apply BH correction across ALL correlation tests in the heatmap
pval_adj <- matrix(p.adjust(as.vector(pval_heat), method = "BH"),
                   nrow = nrow(pval_heat), dimnames = dimnames(pval_heat))

heat_long <- expand.grid(gene = rownames(heat_mat), trait = colnames(heat_mat),
                         stringsAsFactors = FALSE) %>%
  mutate(cor   = as.vector(heat_mat),
         ci_lo = as.vector(ci_lo_mat),
         ci_hi = as.vector(ci_hi_mat),
         pval  = as.vector(pval_heat),
         padj  = as.vector(pval_adj),
         stars = sig_stars(padj),
         module = top_prots$module[match(gene, top_prots$gene)])

# Capitalize module names for facet labels
heat_long$module <- str_to_title(heat_long$module)

gene_order <- top_prots %>%
  mutate(module_title = str_to_title(module)) %>%
  arrange(module_title, desc(kME)) %>%
  pull(gene)
heat_long$gene <- factor(heat_long$gene, levels = rev(gene_order))

# ---- Plot ----
pE <- ggplot(heat_long, aes(x = trait, y = gene, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%s", cor, stars)), size = 1.5, color = "black",
            na.rm = TRUE) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r",
                       na.value = "white") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  labs(title    = "E  Protein-Phenotype Associations",
       subtitle = sprintf("Top hub proteins (baseline vs traits) | BH-corrected across %d tests | Delta traits adjusted for age group",
                          sum(!is.na(pval_heat))),
       x = NULL, y = NULL) +
  THEME_PUB +
  LEGEND_THEME +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y  = element_text(size = 4.5),
        legend.key.width   = unit(15, "mm"),
        legend.key.height  = unit(2, "mm"),
        strip.text.y = element_text(angle = 0, size = 5.5))

# ---- Save ----
ggsave(file.path(RPT_DIR, "panel_E_heatmap.pdf"), pE,
       width = 250, height = 350, units = "mm")
ggsave(file.path(RPT_DIR, "panel_E_heatmap.png"), pE,
       width = 250, height = 350, units = "mm", dpi = 300)

write_csv(heat_long, file.path(DAT_DIR, "05_panel_E_heatmap.csv"))

cat("Panel E done\n")
