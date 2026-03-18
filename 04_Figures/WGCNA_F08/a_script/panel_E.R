# Figure 6 — Panel E: Protein-Phenotype Heatmap
# Top hub proteins (baseline expression) vs Age, VL, LBM traits
#
# STAT AUDIT (2026-02-27):
# 1. Pearson r (Age, Baseline VL) + 95% CI via cor.test(). PASS.
# 2. Partial r (Delta VL, Delta LBM) controlling for age + Fisher z CI.
# 3. BH correction across all n_top*4 tests. PASS.
# 4. Stars based on padj. Exported CSV includes CIs.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ppcor)
})

RPT <- "04_Figures/WGCNA_F08/b_reports"
DAT <- "04_Figures/WGCNA_F08/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

module_df  <- read_csv(file.path(DAT, "module_df.csv"), show_col_types = FALSE)
kME_all    <- readRDS(file.path(DAT, "kME_all.rds"))
gs_vl      <- readRDS(file.path(DAT, "gs_vl.rds"))
pre_expr   <- readRDS(file.path(DAT, "pre_expr.rds"))
pheno_wide <- read_csv(file.path(DAT, "pheno_wide.csv"), show_col_types = FALSE)
subj_age   <- read_csv(file.path(DAT, "subj_age.csv"), show_col_types = FALSE)
shared     <- readRDS(file.path(DAT, "shared_objects.rds"))

top3          <- shared$top3
common_subj   <- shared$common_subj
delta_vl_vec  <- shared$delta_vl_vec
delta_lbm_vec <- shared$delta_lbm_vec

message("Panel E: protein-phenotype heatmap...")

PE_W <- 250  # panel width mm
PE_H <- 350  # panel height mm

txt_cell  <- scale_text(BASE_GENE, PE_W) * 0.5
txt_gene  <- scale_text(BASE_GENE, PE_W) * 0.7
txt_axis  <- scale_text(BASE_STAT, PE_W)
txt_strip <- scale_text(BASE_GENE, PE_W) * 0.85

# Top 15 proteins per top-3 module by kME * |GS| score
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

  ct <- tryCatch(cor.test(base_vals, age_num_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Age"]  <- ct$estimate
    pval_heat[i, "Age"] <- ct$p.value
    ci_lo_mat[i, "Age"] <- ct$conf.int[1]
    ci_hi_mat[i, "Age"] <- ct$conf.int[2]
  }

  ct <- tryCatch(cor.test(base_vals, base_vl_vec), error = function(e) NULL)
  if (!is.null(ct)) {
    heat_mat[i, "Baseline VL"]  <- ct$estimate
    pval_heat[i, "Baseline VL"] <- ct$p.value
    ci_lo_mat[i, "Baseline VL"] <- ct$conf.int[1]
    ci_hi_mat[i, "Baseline VL"] <- ct$conf.int[2]
  }

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

heat_long$module <- str_to_title(heat_long$module)

gene_order <- top_prots %>%
  mutate(module_title = str_to_title(module)) %>%
  arrange(module_title, desc(kME)) %>%
  pull(gene)
heat_long$gene <- factor(heat_long$gene, levels = rev(gene_order))

pE <- ggplot(heat_long, aes(x = trait, y = gene, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%s", cor, stars)), size = txt_cell,
            color = "black", na.rm = TRUE) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r",
                       na.value = "white") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  labs(title    = "Protein-Phenotype Associations",
       subtitle = sprintf("Top hub proteins (baseline vs traits) | BH-corrected across %d tests | Delta traits adjusted for age group",
                          sum(!is.na(pval_heat))),
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(legend.position  = "none",
        axis.text.x      = element_text(angle = 45, hjust = 1, size = txt_axis),
        axis.text.y      = element_text(size = txt_gene),
        strip.text.y     = element_text(angle = 0, size = txt_strip))

ggsave(file.path(RPT, "panel_E_heatmap.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_E_heatmap.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)

write_csv(heat_long, file.path(DAT, "05_panel_E_heatmap.csv"))

message("  Panel E saved")
