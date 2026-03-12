# Panel D — Baseline pathway activity predicts training response.
# Q3: Does baseline LV activity predict delta VL thickness / delta LBM?
# Method: Pearson correlation per LV, BH correction.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

DAT <- "04_Figures/F5_PLIER/c_data"
RPT <- "04_Figures/F5_PLIER/b_results"

# --- Load data ---
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
meta <- as.data.frame(dal$metadata) |>
  mutate(subject = sub("_(Pre|Post)$", "", Col_ID))

lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))

lv_cols <- grep("^LV\\d+$", names(lv_meta), value = TRUE)

# --- Compute delta phenotype per subject ---
pheno_wide <- lv_meta |>
  select(subject, Group, Timepoint, VL_thick_cm, DXA_LBM_kg) |>
  pivot_wider(names_from = Timepoint, values_from = c(VL_thick_cm, DXA_LBM_kg),
              names_sep = "_")

pheno_wide <- pheno_wide |>
  mutate(
    delta_VL = VL_thick_cm_Post - VL_thick_cm_Pre,
    delta_LBM = DXA_LBM_kg_Post - DXA_LBM_kg_Pre
  ) |>
  filter(!is.na(delta_VL) | !is.na(delta_LBM))

# Baseline LV scores
pre <- lv_meta |>
  filter(Timepoint == "Pre") |>
  mutate(subject = sub("_(Pre|Post)$", "", sample_id))

pred_df <- pre |>
  inner_join(pheno_wide |> select(subject, Group, delta_VL, delta_LBM),
             by = c("subject", "Group"))

cat(sprintf("Baseline prediction: n=%d subjects with phenotype data\n", nrow(pred_df)))

# --- Load pathway labels ---
align_df <- read_csv(file.path(DAT, "05_lv_pathway_alignment_summary.csv"),
                     show_col_types = FALSE)
top_label <- align_df |>
  filter(FDR < 0.05) |>
  group_by(LV) |>
  slice_max(AUC, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(lv_name = paste0("LV", LV),
         pathway_label = clean_pathway_name(pathway))

# --- Correlate baseline LV with delta phenotype ---
cor_results <- lapply(lv_cols, function(lv) {
  res_vl <- if (sum(!is.na(pred_df$delta_VL)) >= 5) {
    ct <- cor.test(pred_df[[lv]], pred_df$delta_VL, method = "pearson",
                   use = "complete.obs")
    tibble(LV = lv, phenotype = "delta_VL", r = ct$estimate,
           pvalue = ct$p.value, n = sum(!is.na(pred_df$delta_VL)))
  } else NULL

  res_lbm <- if (sum(!is.na(pred_df$delta_LBM)) >= 5) {
    ct <- cor.test(pred_df[[lv]], pred_df$delta_LBM, method = "pearson",
                   use = "complete.obs")
    tibble(LV = lv, phenotype = "delta_LBM", r = ct$estimate,
           pvalue = ct$p.value, n = sum(!is.na(pred_df$delta_LBM)))
  } else NULL

  bind_rows(res_vl, res_lbm)
}) |> bind_rows() |>
  group_by(phenotype) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  ungroup() |>
  left_join(top_label |> select(lv_name, pathway_label),
            by = c("LV" = "lv_name")) |>
  mutate(label = ifelse(is.na(pathway_label), LV, pathway_label),
         stars = sig_stars(padj))

write_csv(cor_results, file.path(DAT, "14_panel_D_baseline_prediction.csv"))

# --- Scatter plot helper ---
make_scatter <- function(pheno_col, pheno_name, top_n = 4) {
  sub_res <- cor_results |>
    filter(phenotype == pheno_col) |>
    arrange(pvalue) |>
    slice_head(n = top_n)

  if (nrow(sub_res) == 0) return(NULL)

  plot_data <- pred_df |>
    select(subject, Group, all_of(sub_res$LV), delta_pheno = !!sym(pheno_col)) |>
    pivot_longer(cols = all_of(sub_res$LV), names_to = "LV", values_to = "lv_score") |>
    filter(!is.na(delta_pheno)) |>
    left_join(sub_res |> select(LV, label, r, padj, stars), by = "LV") |>
    mutate(facet_label = sprintf("%s (r=%.2f%s)", label, r, stars))

  ggplot(plot_data, aes(x = lv_score, y = delta_pheno, color = Group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.5,
                fill = "grey80", alpha = 0.2) +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 2) +
    scale_color_manual(values = AGE_COLORS) +
    labs(title = sprintf("D  Baseline LV Predicts %s", pheno_name),
         subtitle = sprintf("Pearson r, top %d LVs by p-value", top_n),
         x = "Baseline LV Score", y = paste0("\u0394 ", pheno_name)) +
    FIG_THEME +
    theme(legend.position = "bottom")
}

# VL thickness
pD_vl <- make_scatter("delta_VL", "VL Thickness (cm)")
if (!is.null(pD_vl)) {
  ggsave(file.path(RPT, "14_panel_D_baseline_prediction_vl.pdf"), pD_vl,
         width = 200, height = 170, units = "mm")
  ggsave(file.path(RPT, "14_panel_D_baseline_prediction_vl.png"), pD_vl,
         width = 200, height = 170, units = "mm", dpi = 300)
}

# LBM
pD_lbm <- make_scatter("delta_LBM", "LBM (kg)")
if (!is.null(pD_lbm)) {
  ggsave(file.path(RPT, "14_panel_D_baseline_prediction_lbm.pdf"), pD_lbm,
         width = 200, height = 170, units = "mm")
  ggsave(file.path(RPT, "14_panel_D_baseline_prediction_lbm.png"), pD_lbm,
         width = 200, height = 170, units = "mm", dpi = 300)
}

# Summary
for (ph in c("delta_VL", "delta_LBM")) {
  sub <- cor_results |> filter(phenotype == ph, padj < 0.05)
  cat(sprintf("Panel D (%s): %d/%d LVs with padj < 0.05\n",
              ph, nrow(sub), sum(cor_results$phenotype == ph)))
}
