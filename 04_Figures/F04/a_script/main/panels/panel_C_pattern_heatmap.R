# F04 Panel C: Training Concordance Pattern Heatmap
# Config wrapper for shared/comparison_panels/panel_C_pattern_heatmap.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

ROW_H <- 0.078  # needed for key position calculations below

cfg <- list(
  fig_id     = "F04",
  contrast_x = "Training_Young",
  contrast_y = "Training_Old",
  title      = "Training Response Patterns",
  col_headers = c("Tr.(Y)", "Tr.(O)"),
  sort_col   = "logFC_Training_Young",

  rpt_png = "04_Figures/F04/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F04/b_reports/main/pdf/panels",
  dat     = "04_Figures/F04/c_data",

  classify_fn = function(dep_df) {
    dep_df %>%
      dplyr::filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
      dplyr::filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
             pi_score_Interaction < 0.05) %>%
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
          logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
          TRUE ~ "Discordant"
        ),
        sig_cat = dplyr::case_when(
          pi_score_Training_Young < 0.05 & pi_score_Training_Old < 0.05 ~ "Both",
          pi_score_Training_Young < 0.05 ~ "Tr.(Y)",
          pi_score_Training_Old < 0.05   ~ "Tr.(O)",
          pi_score_Interaction < 0.05    ~ "Inter.",
          TRUE ~ "NS"
        )
      )
  },

  QUAD_ORDER      = c("Concordant Up", "Concordant Down", "Discordant"),
  QUAD_COLORS     = c("Concordant Up" = "#D32F2F", "Concordant Down" = "#1976D2",
                      "Discordant" = "#388E3C"),
  QUAD_BG         = c("Concordant Up" = "#FFCDD2", "Concordant Down" = "#BBDEFB",
                      "Discordant" = "#C8E6C9", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Concordant Up" = "#8B0000", "Concordant Down" = "#0D47A1",
                      "Discordant" = "#1B5E20"),
  SIG_COLORS      = c("Both" = "#2E7D32", "Tr.(Y)" = "#E05A4E",
                      "Tr.(O)" = "#5DA5DA", "Inter." = "#7B5EA7", "NS" = "grey70"),

  display_labels = c(
    "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
    "Amino Acid & Cofactor Metabolism" = "AA & Cofactor Metab."
  ),

  col_header_colors = c(
    CONTRAST_COLORS["Training_Young"],   # #E05A4E for "Tr.(Y)"
    CONTRAST_COLORS["Training_Old"]      # #5DA5DA for "Tr.(O)"
  ),

  bar_scale            = 0.38,
  bar_ref_width        = 19,
  key_y_base           = ROW_H * 7.5,
  key_dy               = ROW_H * 2.2,
  key_x_sig            = NULL,
  protein_count_x_mult = 7.5,
  count_tick_y_label   = ROW_H * 2.3,
  count_tick_filter    = function(df) df,

  sig_cats       = c("Tr.(Y)", "Tr.(O)", "Both", "Inter."),
  sig_cat_labels = c("Sig Young", "Sig Old", "Sig Both", "Interaction")
)

source("04_Figures/shared/comparison_panels/panel_C_pattern_heatmap.R")
