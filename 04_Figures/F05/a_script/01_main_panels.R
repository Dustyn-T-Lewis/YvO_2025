#!/usr/bin/env Rscript
# F05 Main — Aging Reversal (5-panel composite)
# A: Quadrant ORA  B: Pattern heatmap  C: fry barcode
# D: NES scatter   E: RRHO2
#
# 3-column geometry-aware layout (panels labeled in visual reading order):
#   Top row:    A (Quadrant ORA)  | B (Pattern heatmap, full-height)
#   Bottom row: C (fry barcode)   | D (NES scatter) | E (RRHO2)
#
# Panel B (pattern heatmap) loads AnnotationDbi (via go_slim_categories.R);
# the S4 select() masking is repaired inside go_slim_categories.R.

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

source(here::here("04_Figures", "shared", "style.R"))

BASE    <- here::here("04_Figures", "F05")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ======================================================================
# Source panels (A first, then D/C/E, then B last for stat-snapshot ordering)
# ======================================================================

message("=== F05 Composite: sourcing panels ===")

# Panel A — ORA scatter + flanking bars (complex, extracted helper)
source(here::here("04_Figures", "F05", "a_script", "_panel_A_ORA.R"))
n_total_A    <- nrow(scatter_df)
n_sig_A      <- n_sig
n_enrich_A   <- n_enrich
r_spear_A    <- r_spear

# Panel D — fGSEA NES scatter (config wrapper -> shared engine)
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
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
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
source(here::here("04_Figures", "shared", "comparison_panels", "panel_D_nes_scatter.R"))
n_pw_D       <- nrow(fgsea_wide)
n_sig_pw_D   <- n_total_sig
rho_D        <- as.numeric(nes_cor_all$estimate)
rho_lo_D     <- nes_ci_all[1]
rho_hi_D     <- nes_ci_all[2]
pw_rev_D     <- pw_rev_frac

# Panel C — fry rotation test (config wrapper -> shared engine)
ROW_H_fry <- 0.078
cfg <- list(
  fig_id           = "F05",
  contrast_source  = "Aging",
  contrast_test    = "Training_Old",
  set_prefix       = "aging",
  expected_up      = "Down",
  expected_down    = "Up",
  driving_up_sign  = "neg",   # aging-up set reversal = t_TO < 0
  driving_dn_sign  = "pos",   # aging-down set reversal = t_TO > 0
  has_circularity  = TRUE,
  up_title_fmt = "Aging-Up DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t",
  dn_title_fmt = "Aging-Down DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t%s",
  fig_color    = unname(CONTRAST_COLORS["Aging"]),  # #4CAF50
  stat_corner_up = "bottomleft",
  stat_corner_dn = "topright",
  ora_flank_up_label = "Reversed (Up\u2192Down)",
  ora_flank_dn_label = "Reversed (Down\u2192Up)",
  ora_supp_up_label  = "Aging-Up (Reversed)",
  ora_supp_dn_label  = "Aging-Down (Reversed)",
  ora_supp_title     = "Leading-Edge ORA: fry Driving Proteins (Reversal)",
  ora_supp_subtitle  = "Hypergeometric ORA on reversal-driving proteins | top 3 per set",
  label_map = c(
    "Amino Acid Catabolic Process"                           = "AA Catabol.",
    "Fatty Acid Catabolic Process"                           = "FA Catabol.",
    "Amino Acid Metabolic Process"                           = "AA Metabolic Process",
    "Establishment Or Maintenance Of Cell Polarity"          = "Cell Polarity",
    "Generation Of Precursor Metabolites And Energy"         = "Precursor Metab. & Energy",
    "Protein Localization To Plasma Membrane"                = "Plasma Membrane Protein Loc.",
    "Organic Acid Catabolic Process"                         = "Organic Acid Catabolism",
    "Membraneless Organelle Assembly"                        = "Membraneless Org. Assembly"
  ),
  force_inside_labels = c("AA Catabol.", "FA Catabol."),
  long_label_mode     = "truncate",
  title        = "fry Gene-Set Rotation Test: Aging Reversal",
  subtitle_fmt = "Rotation-based set test (exact GSEA analogue) | Circularity r = %.3f | dupCor = %.3f | n = %d proteins",
  panel_w     = 178,
  rpt_png     = PNL_PNG, rpt_pdf = PNL_PDF,
  rpt_sup_png = file.path(BASE, "b_reports", "supp", "png", "panels"),
  rpt_sup_pdf = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
  dat         = DAT
)
source(here::here("04_Figures", "shared", "comparison_panels", "panel_C_fry.R"))
cor_imp_C    <- cor_imp
n_all_C      <- n_all
circ_r_C     <- circ_r

# Panel E — RRHO2 (config wrapper -> shared engine)
cfg <- list(
  fig_id     = "F05",
  t_col_1    = "t_Aging",
  t_col_2    = "t_Training_Old",
  rrho_labels = c("Aging", "Training (Old)"),
  title       = "Threshold-Free Reversal (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm off-diagonal = training reverses aging | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Aging rank"~(Up %->% Down)),
  axis_label_2 = expression("Training (Old) rank"~(Up %->% Down)),
  quadrant_labels = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (Age\u2191 Tr\u2193)",
    DU = "Reversed (Age\u2193 Tr\u2191)"
  ),
  hotspot_export_names = list(
    UU = "Exacerbated_Up",
    DD = "Exacerbated_Down",
    UD = "Reversed_AgingUp_TrainingDown",
    DU = "Reversed_AgingDown_TrainingUp"
  ),
  ora_min_size = 10,
  ora_quadrant_names = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (Aging Up / Training Down)",
    DU = "Reversed (Aging Down / Training Up)"
  ),
  ora_grouped = list(
    file_1_quads     = c("ora_UD", "ora_DU"),    # reversal -> concordant file
    file_2_quads     = c("ora_UU", "ora_DD"),    # exacerbation -> discordant file
    note_if_empty_2  = "No pathways enriched in exacerbation quadrants (padj<0.05)"
  ),
  ora_colors = ORA_QUAD_COLORS_F3,
  summary_quadrant_names = list(
    UU = "Exacerbated_Up",                  UU_slug = "exacerbated_up",
    DD = "Exacerbated_Down",                DD_slug = "exacerbated_down",
    UD = "Reversed_AgingUp_TrainingDown",   UD_slug = "reversed_aging_up",
    DU = "Reversed_AgingDown_TrainingUp",   DU_slug = "reversed_aging_down"
  ),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  supp = list(
    rpt_png        = file.path(BASE, "b_reports", "supp", "png", "panels"),
    rpt_pdf        = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
    ora_bar_title  = "Enriched Pathways by Reversal Quadrant",
    ora_quad_order = c("Reversed (Aging Up / Training Down)",
                       "Reversed (Aging Down / Training Up)",
                       "Exacerbated Up", "Exacerbated Down"),
    ora_quad_short = c(
      "Reversed (Aging Up / Training Down)"  = "Reversed\n(Up \u2192 Down)",
      "Reversed (Aging Down / Training Up)"  = "Reversed\n(Down \u2192 Up)",
      "Exacerbated Up"                       = "Exacerbated Up",
      "Exacerbated Down"                     = "Exacerbated Down"
    )
  )
)
source(here::here("04_Figures", "shared", "comparison_panels", "panel_E_rrho2.R"))
n_shared_E   <- n_shared
max_rev_E    <- max(max_UD, max_DU)
n_rev_E      <- if (max_UD >= max_DU) n_UD else n_DU

# Panel B — Pattern heatmap (last: loads AnnotationDbi, clobbers select())
ROW_H <- 0.078
cfg <- list(
  fig_id     = "F05",
  contrast_x = "Aging",
  contrast_y = "Training_Old",
  title      = "Aging Reversal Patterns",
  col_headers = c("Aging", "Tr.(O)"),
  sort_col   = "logFC_Aging",
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  classify_fn = function(dep_df) {
    dep_df |>
      dplyr::filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) |>
      dplyr::filter(pi_score_Aging < 0.05 | pi_score_Training_Old < 0.05) |>
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
  QUAD_COLORS     = c("Reversed Up" = "#B2182B", "Reversed Down" = "#2166AC",
                      "Non-reversed" = "#1B7837"),
  QUAD_BG         = c("Reversed Up" = "#F4D9D2", "Reversed Down" = "#D5DEEF",
                      "Non-reversed" = "#C8E0CD", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Reversed Up" = "#67001F", "Reversed Down" = "#053061",
                      "Non-reversed" = "#00441B"),
  SIG_COLORS      = c("Both" = "#2E7D32", "Aging" = "#E05A4E",
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
source(here::here("04_Figures", "shared", "comparison_panels", "panel_B_pattern_heatmap.R"))

# Restore RPT paths (shared engines clobber RPT_PDF/RPT_PNG to panels subdir)
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

# ======================================================================
# Composite layout (3-column, geometry-aware)
# ======================================================================

COMP_W     <- 420
COMP_H     <- 310
PRINT_SCALE <- 380 / 178
TAG_SZ     <- round(10 * PRINT_SCALE * 0.85)
TTL_SZ     <- round(10 * PRINT_SCALE * 0.85)
SUB_SZ     <- round(7 * PRINT_SCALE * 0.85)

# Stat-snapshot title strings
ttl_A <- "Quadrant ORA (Reversal)"
sub_A <- sprintf("N = %d | %d DEPs (\u03a0) | %d enriched (FDR) | \u03c1 = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "fry: Reversal"
sub_C <- sprintf("n = %d | dupCor = %.3f", n_all_C, cor_imp_C)
ttl_D <- "Pathway Reversal"
sub_D <- sprintf("\u03c1 = %.2f | %.0f%% reversed",
                 rho_D, pw_rev_D * 100)
ttl_E <- "RRHO2 Reversal"
sub_E <- sprintf("%d genes | max %d", n_shared_E, n_rev_E)

# Patchwork layout
layout <- paste(
  "##############",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "##############",
  "##############",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  sep = "\n"
)

composite <- composite + plot_annotation(theme = theme(plot.margin = margin(-2.5, -1, -2.5, -1, "mm")))
pC_fry  <- pC_fry + plot_annotation(theme = theme(plot.margin = margin(3, 5, 0, 0, "mm")))
pD      <- pD + theme(plot.margin = margin(-2.8, 5, 2.8, -5, "mm"))
pE_heat <- pE_heat + theme(plot.margin = margin(-2.1, -0.2, 3.4, -3.5, "mm"),
                           axis.title = element_text(face = "bold", size = 8))
pB <- pB + coord_cartesian(xlim = c(-0.25, X_BAR_MAX + 1.75),
                           ylim = c(BAR_YMAX + ROW_H * 4.5, -ROW_H * 0.05),
                           expand = FALSE)

fig <- wrap_elements(full = composite) +
       wrap_elements(full = pB) +
       wrap_elements(full = pC_fry) +
       wrap_elements(full = pD) +
       wrap_elements(full = pE_heat) +
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6)))

# Panel tags, titles, subtitles via cowplot
X_A <- 0.005;  X_B <- 0.549;  X_C <- 0.012;  X_D <- 0.406;  X_E <- 0.693
X_TTL      <- 0.030
TAG_DY     <- -0.002
SUB_OFFSET <- 0.020

Y_A <- 0.984;  Y_B <- 0.984
Y_C <- 0.512;  Y_D <- 0.511;  Y_E <- 0.511

composite_final <- ggdraw(fig) +
  draw_label("A",    x = X_A,           y = Y_A - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A,  x = X_A + X_TTL,   y = Y_A,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A,  x = X_A + X_TTL,   y = Y_A - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",    x = X_B,           y = Y_B - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B,  x = X_B + X_TTL,   y = Y_B,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B,  x = X_B + X_TTL,   y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",    x = X_C,           y = Y_C - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C,  x = X_C + X_TTL,   y = Y_C,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C,  x = X_C + X_TTL,   y = Y_C - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",    x = X_D,           y = Y_D - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D,  x = X_D + X_TTL,   y = Y_D,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D,  x = X_D + X_TTL,   y = Y_D - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E",    x = X_E,           y = Y_E - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E,  x = X_E + X_TTL,   y = Y_E,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E,  x = X_E + X_TTL,   y = Y_E - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

ggsave(file.path(RPT_PDF, "MAIN_F05_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F05_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F05 composite (5-panel, 3-column layout) saved")
