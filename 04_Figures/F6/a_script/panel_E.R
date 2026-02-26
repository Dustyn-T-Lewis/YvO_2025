################################################################################
#   Figure 6 — Panel E: Protein-Phenotype Heatmap
#   Top hub proteins (baseline expression) correlated with Age, VL, LBM traits
################################################################################

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

# ---- Build correlation + p-value matrices ----
n_top   <- nrow(top_prots)
traits  <- c("Age", "Baseline VL", "Delta VL", "Delta LBM")
heat_mat  <- matrix(NA_real_, nrow = n_top, ncol = length(traits),
                    dimnames = list(top_prots$gene, traits))
pval_heat <- heat_mat

age_num_vec <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0
)
base_vl_vec <- pheno_wide$VL_Pre[match(common_subj, pheno_wide$subject_key)]

for (i in seq_len(n_top)) {
  prot      <- top_prots$uniprot_id[i]
  base_vals <- pre_expr[common_subj, prot]

  # Age
  ct <- tryCatch(cor.test(base_vals, age_num_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Age"]  <- ct$estimate
    pval_heat[i, "Age"] <- ct$p.value
  }

  # Baseline VL
  ct <- tryCatch(cor.test(base_vals, base_vl_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Baseline VL"]  <- ct$estimate
    pval_heat[i, "Baseline VL"] <- ct$p.value
  }

  # Delta VL
  ct <- tryCatch(cor.test(base_vals, delta_vl_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Delta VL"]  <- ct$estimate
    pval_heat[i, "Delta VL"] <- ct$p.value
  }

  # Delta LBM
  ct <- tryCatch(cor.test(base_vals, delta_lbm_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Delta LBM"]  <- ct$estimate
    pval_heat[i, "Delta LBM"] <- ct$p.value
  }
}

# ---- Assemble long-format heatmap data ----
# Apply BH correction across ALL correlation tests in the heatmap
pval_adj <- matrix(p.adjust(as.vector(pval_heat), method = "BH"),
                   nrow = nrow(pval_heat), dimnames = dimnames(pval_heat))

heat_long <- expand.grid(gene = rownames(heat_mat), trait = colnames(heat_mat),
                         stringsAsFactors = FALSE) %>%
  mutate(cor   = as.vector(heat_mat),
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
       subtitle = sprintf("Top hub proteins (baseline vs traits) | BH-corrected across %d tests",
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

write_csv(heat_long, file.path(DAT_DIR, "fig6_panel_E_heatmap.csv"))

cat("Panel E done\n")
