# F05 Panel B: fGSEA NES Scatter — Aging Reversal
# Config wrapper for shared/comparison_panels/panel_B_nes_scatter.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id      = "F05",
  contrast_x  = "Aging",
  contrast_y  = "Training_Old",
  title       = "Pathway-Level Reversal (fGSEA)",
  axis_x_label = "Aging",
  axis_y_label = "Training Old",
  subtitle_metric = "reversed",
  subtitle_interpretation = "negative \u03c1 = training opposes aging effects",
  ref_slope   = -1,
  panel_w     = 146,
  label_border_size = 0.25,

  rpt_png = "04_Figures/F05/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F05/b_reports/main/pdf/panels",
  dat     = "04_Figures/F05/c_data",

  sig_colors     = SIG_COLORS_F3,
  sig_label_fill = SIG_LABEL_FILL_F3,
  sig_label_text = SIG_LABEL_TEXT_F3,
  sig_draw_order = c("Sig Training only", "Sig Aging only", "Sig Both"),

  quadrant_defs = list(
    sig_both_label = "Sig Both",
    sig_x_label    = "Sig Aging only",
    sig_y_label    = "Sig Training only",
    # Reversal: off-diagonal = blue (reversed), diagonal = red (exacerbated)
    bg_blue_1 = c(0, Inf, -Inf, 0),      # bottom-right (reversed)
    bg_blue_2 = c(-Inf, 0, 0, Inf),      # top-left (reversed)
    bg_red_1  = c(0, Inf, 0, Inf),       # top-right (exacerbated)
    bg_red_2  = c(-Inf, 0, -Inf, 0),     # bottom-left (exacerbated)
    label_tr = "Exacerbated",  color_tr = "#D6604D",
    label_tl = "Reversed",     color_tl = "#4393C3",
    label_bl = "Exacerbated",  color_bl = "#D6604D",
    label_br = "Reversed",     color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q2 + q4  # reversed = TL + BR
  ),

  display_overrides = c(
    "Unfolded Protein Response"         = "UPR",
    "Ribosome Biogenesis"               = "Ribo Bio",
    "Amino Acid Metabolism"             = "AA Metabolism",
    "Fatty Acid Metabolism"             = "FA Metabolism",
    "Fatty Acid Beta Oxidation"         = "FA Beta-Oxidation",
    "Plasma Membrane Protein Loc."      = "PM Protein Loc.",
    "Cytoplasmic Translation"           = "Cytoplasmic Transl.",
    "Mitochondrial Organization"        = "Mito Org.",
    "Precursor Metabolites & Energy"    = "Precursor Metab. & Energy",
    "Mitochondrial Transport"           = "Mito Transport",
    "Mitochondrial Protein Import"      = "Mito Protein Import",
    "Mitochondrial Protein Degradation" = "Mito Protein Deg.",
    "Extracellular Matrix Organization" = "ECM Org.",
    "Heme Metabolism"                   = "Heme Metab.",
    "Ketone Metabolism"                 = "Ketone Metab.",
    "Peroxisome"                        = "Peroxisome",
    "Muscle System"                     = "Muscle System"
  )
)

source("04_Figures/shared/comparison_panels/panel_B_nes_scatter.R")
