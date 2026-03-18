# Panel F — Delta pathway vs delta phenotype.
# Q5: Do training-induced pathway changes track phenotype gains?
# Method: Pearson correlation (delta LV vs delta VL/LBM), BH correction.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

DAT <- "04_Figures/PLIER_F09/c_data"
RPT <- "04_Figures/PLIER_F09/b_results"

# --- Load data ---
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
meta <- as.data.frame(dal$metadata) |>
  mutate(subject = sub("_(Pre|Post)$", "", Col_ID))

lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))
lv_cols <- grep("^LV\\d+$", names(lv_meta), value = TRUE)

# --- Compute subject-level deltas ---
lv_wide <- lv_meta |>
  select(subject, Group, Timepoint, all_of(lv_cols)) |>
  pivot_wider(names_from = Timepoint, values_from = all_of(lv_cols),
              names_sep = ".")

# Delta LV scores
delta_lv <- data.frame(subject = lv_wide$subject, Group = lv_wide$Group)
for (lv in lv_cols) {
  delta_lv[[lv]] <- lv_wide[[paste0(lv, ".Post")]] - lv_wide[[paste0(lv, ".Pre")]]
}

# Delta phenotype
pheno_wide <- lv_meta |>
  select(subject, Group, Timepoint, VL_thick_cm, DXA_LBM_kg) |>
  distinct() |>
  pivot_wider(names_from = Timepoint, values_from = c(VL_thick_cm, DXA_LBM_kg),
              names_sep = "_") |>
  mutate(delta_VL = VL_thick_cm_Post - VL_thick_cm_Pre,
         delta_LBM = DXA_LBM_kg_Post - DXA_LBM_kg_Pre)

delta_df <- delta_lv |>
  inner_join(pheno_wide |> select(subject, delta_VL, delta_LBM), by = "subject")

cat(sprintf("Delta correlations: n=%d subjects\n", nrow(delta_df)))

# --- Correlations ---
align_df <- read_csv(file.path(DAT, "05_lv_pathway_alignment_summary.csv"),
                     show_col_types = FALSE)
top_label <- align_df |>
  filter(FDR < 0.05) |>
  group_by(LV) |>
  slice_max(AUC, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(lv_name = paste0("LV", LV),
         pathway_label = clean_pathway_name(pathway))

cor_results <- lapply(lv_cols, function(lv) {
  res <- list()
  if (sum(!is.na(delta_df$delta_VL)) >= 5) {
    ct <- cor.test(delta_df[[lv]], delta_df$delta_VL, use = "complete.obs")
    res$vl <- tibble(LV = lv, phenotype = "delta_VL", r = ct$estimate,
                     pvalue = ct$p.value, n = sum(!is.na(delta_df$delta_VL)))
  }
  if (sum(!is.na(delta_df$delta_LBM)) >= 5) {
    ct <- cor.test(delta_df[[lv]], delta_df$delta_LBM, use = "complete.obs")
    res$lbm <- tibble(LV = lv, phenotype = "delta_LBM", r = ct$estimate,
                      pvalue = ct$p.value, n = sum(!is.na(delta_df$delta_LBM)))
  }
  bind_rows(res)
}) |> bind_rows() |>
  group_by(phenotype) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  ungroup() |>
  left_join(top_label |> select(lv_name, pathway_label),
            by = c("LV" = "lv_name")) |>
  mutate(label = ifelse(is.na(pathway_label), LV, pathway_label),
         stars = sig_stars(padj))

write_csv(cor_results, file.path(DAT, "16_panel_F_delta_vs_delta.csv"))

# --- Scatter plots for top associations ---
make_delta_scatter <- function(pheno_col, pheno_name, top_n = 4) {
  sub_res <- cor_results |>
    filter(phenotype == pheno_col) |>
    arrange(pvalue) |>
    slice_head(n = top_n)

  if (nrow(sub_res) == 0) return(NULL)

  plot_data <- delta_df |>
    select(subject, Group, all_of(sub_res$LV), !!sym(pheno_col)) |>
    pivot_longer(cols = all_of(sub_res$LV), names_to = "LV", values_to = "delta_lv") |>
    rename(delta_pheno = !!pheno_col) |>
    filter(!is.na(delta_pheno)) |>
    left_join(sub_res |> select(LV, label, r, padj, stars), by = "LV") |>
    mutate(facet_label = sprintf("%s (r=%.2f%s)", label, r, stars))

  ggplot(plot_data, aes(x = delta_lv, y = delta_pheno, color = Group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.5,
                fill = "grey80", alpha = 0.2) +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 2) +
    scale_color_manual(values = AGE_COLORS) +
    labs(title = sprintf("F  \u0394 Pathway vs \u0394 %s", pheno_name),
         subtitle = "Pearson r (Post - Pre), BH-adjusted",
         x = "\u0394 LV Score", y = paste0("\u0394 ", pheno_name)) +
    FIG_THEME +
    theme(legend.position = "bottom")
}

pF_vl <- make_delta_scatter("delta_VL", "VL Thickness (cm)")
if (!is.null(pF_vl)) {
  ggsave(file.path(RPT, "16_panel_F_delta_vs_delta_vl.pdf"), pF_vl,
         width = 200, height = 170, units = "mm")
  ggsave(file.path(RPT, "16_panel_F_delta_vs_delta_vl.png"), pF_vl,
         width = 200, height = 170, units = "mm", dpi = 300)
}

pF_lbm <- make_delta_scatter("delta_LBM", "LBM (kg)")
if (!is.null(pF_lbm)) {
  ggsave(file.path(RPT, "16_panel_F_delta_vs_delta_lbm.pdf"), pF_lbm,
         width = 200, height = 170, units = "mm")
  ggsave(file.path(RPT, "16_panel_F_delta_vs_delta_lbm.png"), pF_lbm,
         width = 200, height = 170, units = "mm", dpi = 300)
}

for (ph in c("delta_VL", "delta_LBM")) {
  n_sig <- sum(cor_results$phenotype == ph & cor_results$padj < 0.05, na.rm = TRUE)
  cat(sprintf("Panel F (%s): %d LVs with padj < 0.05\n", ph, n_sig))
}
