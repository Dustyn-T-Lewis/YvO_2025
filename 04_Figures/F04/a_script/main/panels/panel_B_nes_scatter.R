# F04 Panel B: fGSEA NES Scatter — Training Concordance
# Config wrapper for shared/comparison_panels/panel_B_nes_scatter.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id      = "F04",
  contrast_x  = "Training_Young",
  contrast_y  = "Training_Old",
  title       = "Pathway-Level Concordance (fGSEA)",
  axis_x_label = "Training Young",
  axis_y_label = "Training Old",
  subtitle_metric = "concordant",
  subtitle_interpretation = "positive \u03c1 = shared pathway regulation across age groups",
  ref_slope   = 1,
  panel_w     = 146,
  label_border_size = 0.25,

  rpt_png = "04_Figures/F04/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F04/b_reports/main/pdf/panels",
  dat     = "04_Figures/F04/c_data",

  sig_colors     = SIG_COLORS_F2,
  sig_label_fill = SIG_LABEL_FILL_F2,
  sig_label_text = SIG_LABEL_TEXT_F2,
  sig_draw_order = c("Sig Old only", "Sig Young only", "Sig Both", "Interaction"),

  quadrant_defs = list(
    sig_both_label = "Sig Both",
    sig_x_label    = "Sig Young only",
    sig_y_label    = "Sig Old only",
    # Concordance: diagonal = red (same direction), off-diagonal = blue
    bg_red_1  = c(0, Inf, 0, Inf),       # top-right
    bg_red_2  = c(-Inf, 0, -Inf, 0),     # bottom-left
    bg_blue_1 = c(0, Inf, -Inf, 0),      # bottom-right
    bg_blue_2 = c(-Inf, 0, 0, Inf),      # top-left
    label_tr = "Concordant Up",   color_tr = "#D6604D",
    label_tl = "Discordant",      color_tl = "#4393C3",
    label_bl = "Concordant Down", color_bl = "#D6604D",
    label_br = "Discordant",      color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q1 + q3  # concordant = TR + BL
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
    "Protein Folding"                   = "Protein Folding",
    "Spermatogenesis"                   = "Spermatogenesis",
    "Adipogenesis"                      = "Adipogenesis",
    "Heme Metabolism"                   = "Heme Metab.",
    "Ketone Metabolism"                 = "Ketone Metab.",
    "Detoxification"                    = "Detoxification"
  )
)

source("04_Figures/shared/comparison_panels/panel_B_nes_scatter.R")
