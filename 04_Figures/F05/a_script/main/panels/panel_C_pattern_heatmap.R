# F05 Panel C: Aging Reversal Pattern Heatmap
# Config wrapper for shared/comparison_panels/panel_C_pattern_heatmap.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

ROW_H <- 0.078

cfg <- list(
  fig_id     = "F05",
  contrast_x = "Aging",
  contrast_y = "Training_Old",
  title      = "Aging Reversal Patterns",
  col_headers = c("Aging", "Tr.(O)"),
  sort_col   = "logFC_Aging",

  rpt_png = "04_Figures/F05/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F05/b_reports/main/pdf/panels",
  dat     = "04_Figures/F05/c_data",

  classify_fn = function(dep_df) {
    dep_df %>%
      dplyr::filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
      dplyr::filter(pi_score_Aging < 0.05 | pi_score_Training_Old < 0.05) %>%
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_Aging > 0 & logFC_Training_Old < 0 ~ "Reversed Up",
          logFC_Aging < 0 & logFC_Training_Old > 0 ~ "Reversed Down",
          TRUE ~ "Non-reversed"
        ),
        sig_cat = dplyr::case_when(
          pi_score_Aging < 0.05 & pi_score_Training_Old < 0.05 ~ "Both",
          pi_score_Aging < 0.05          ~ "Aging",
          pi_score_Training_Old < 0.05   ~ "Tr.(O)",
          TRUE ~ "NS"
        )
      )
  },

  QUAD_ORDER      = c("Reversed Up", "Reversed Down", "Non-reversed"),
  QUAD_COLORS     = c("Reversed Up" = "#D32F2F", "Reversed Down" = "#1976D2",
                      "Non-reversed" = "#388E3C"),
  QUAD_BG         = c("Reversed Up" = "#FFCDD2", "Reversed Down" = "#BBDEFB",
                      "Non-reversed" = "#C8E6C9", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Reversed Up" = "#8B0000", "Reversed Down" = "#0D47A1",
                      "Non-reversed" = "#1B5E20"),
  SIG_COLORS      = c("Both" = "#2E7D32", "Aging" = "#4CAF50",
                      "Tr.(O)" = "#5DA5DA", "NS" = "grey70"),

  display_labels = c(
    "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
    "Amino Acid & Cofactor Metabolism" = "AA & Cofactor\nMetab."
  ),

  col_header_colors = c(
    CONTRAST_COLORS["Aging"],            # #4CAF50 for "Aging"
    CONTRAST_COLORS["Training_Old"]      # #5DA5DA for "Tr.(O)"
  ),

  bar_scale            = 0.20,
  bar_ref_width        = 35,
  key_y_base           = ROW_H * 15.5,
  key_dy               = ROW_H * 3.8,
  key_x_sig            = NULL,
  protein_count_x_mult = 15,
  count_tick_y_label   = ROW_H * 2.6,
  count_tick_filter    = function(df) dplyr::filter(df, val != 15),

  sig_cats       = c("Tr.(O)", "Aging", "Both"),
  sig_cat_labels = c("Sig Training", "Sig Aging", "Sig Both")
)

source("04_Figures/shared/comparison_panels/panel_C_pattern_heatmap.R")
