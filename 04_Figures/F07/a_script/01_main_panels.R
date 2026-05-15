#!/usr/bin/env Rscript
# F07 Main — Phenotype Prediction (2-panel composite)
# A (top) : Top-12 per-module ROCs ranked by permutation p
# B (bot) : Module-Phenotype Coupling age-dependent 2x3 hero scatter grid
#
# Sources ALL data-producing scripts in dependency order (SUPP panels first
# to generate CSVs needed by panel_A and the xlsx), then builds the main
# composite + xlsx + cleanup.

setwd(rprojroot::find_rstudio_root_file())

library(patchwork)
library(cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")

BASE    <- "04_Figures/F07"
DAT     <- file.path(BASE, "c_data")

# Source all data-producing scripts in dependency order

message("=== F07 Composite: sourcing all panels ===")

# 1. SUPP module grid (writes c_data/module_grid/*.csv needed by panel_A)
source("04_Figures/F07/a_script/_supp_module_grid.R")

# 2. MAIN panel A (reads module_grid CSVs)
source("04_Figures/F07/a_script/_panel_A_auc_bars.R")

# 3. MAIN panel B (writes panel_B_full_screen_bh.csv)
source("04_Figures/F07/a_script/_panel_B_hero_grid.R")

# 4. SUPP panel B grid (reads panel_B_full_screen_bh.csv)
source("04_Figures/F07/a_script/_supp_panel_B_grid.R")

# 5. SUPP prepare ROC data (writes classifier_pilot_*.csv)
source("04_Figures/F07/a_script/_supp_prepare_roc.R")

# 6. SUPP ROC panel (reads classifier_pilot_*.csv)
source("04_Figures/F07/a_script/_supp_roc_panel.R")

# 7. SUPP multivariate classifier (writes panel_A_*.csv)
source("04_Figures/F07/a_script/_supp_multivariate.R")

# 8. SUPP LOSO sensitivity (writes loso_auc/loso_auc_summary.csv)
source("04_Figures/F07/a_script/_supp_loso_sensitivity.R")

# 9. SUPP LOSO WGCNA refit (~13-15 min; writes loso_auc/loso_wgcna_refit_*.csv)
source("04_Figures/F07/a_script/_supp_loso_wgcna_refit.R")

# Build main composite (pA/pB via patchwork + cowplot tags)

# Restore paths clobbered by sourced helper scripts
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
DAT     <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

layout_cfg <- list(
  # Canvas dimensions (mm)
  comp_w  = 210,
  a_h     = 75,   # Panel A height (matches ggsave W=210, H=75)
  b_h     = 90,   # Panel B height (matches ggsave W=178, H=90)
  tag_a   = 10,   # spacer above panel A (mm)
  tag_b   = 12,   # spacer above panel B (mm)
  # Text size adjustments relative to composite_text_sizes() base
  tag_bump = 3,
  ttl_bump = 3,
  sub_bump = 2.5,
  # Label x positions (normalized)
  tag_x   = 4,    # mm from left edge (divided by comp_w in use)
  x_ttl   = 0.040,
  # Tag y nudge: small correction applied to tag vs title y
  tag_dy  = -0.002,
  # Subtitle drops below title y by this amount
  sub_offset_base = 0.016  # added to (1 / comp_h) in use
)

# Derived dimensions
COMP_W <- layout_cfg$comp_w
A_H    <- layout_cfg$a_h
B_H    <- layout_cfg$b_h
TAG_A  <- layout_cfg$tag_a
TAG_B  <- layout_cfg$tag_b
COMP_H <- TAG_A + A_H + TAG_B + B_H  # 187 mm

composite <- (plot_spacer() / pA / plot_spacer() / pB) +
  plot_layout(heights = c(TAG_A, A_H, TAG_B, B_H))

txt    <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag   + layout_cfg$tag_bump
TTL_SZ <- txt$title + layout_cfg$ttl_bump
SUB_SZ <- txt$subtitle + layout_cfg$sub_bump

TAG_X      <- layout_cfg$tag_x / COMP_W
X_TTL      <- layout_cfg$x_ttl
TAG_DY     <- layout_cfg$tag_dy
SUB_OFFSET <- layout_cfg$sub_offset_base + (1 / COMP_H)
DY5        <- 5 / COMP_H

TAG_Y_A <- 0.995 - DY5 + (3 / COMP_H)
# TAG_Y_B broken into named intermediates for clarity:
#   base    = fraction of canvas occupied by panel B + its spacer
#   spacer_frac = the spacer's proportion of canvas, scaled to ~85% of its height
#   fine_adj = small corrections totalling +4 - 5 - 10 mm converted to fractions
B_BASE_FRAC  <- (B_H + TAG_B) / COMP_H
SPACER_FRAC  <- (TAG_B / COMP_H) * 0.85
FINE_ADJ     <- (-DY5 * 0.4) - (5 / COMP_H) - (10 / COMP_H) + (4 / COMP_H)
TAG_Y_B <- B_BASE_FRAC + SPACER_FRAC + FINE_ADJ

composite <- composite & theme(plot.margin = margin(0, 4, 0, 5))

composite_final <- ggdraw(composite) +
  # Panel A tag + title + subtitle
  draw_label("A", x = TAG_X, y = TAG_Y_A - TAG_DY, size = TAG_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label(pA_title, x = TAG_X + X_TTL, y = TAG_Y_A,
             size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = TAG_X + X_TTL, y = TAG_Y_A - SUB_OFFSET,
             size = SUB_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  # Panel B tag + title + subtitle
  draw_label("B", x = TAG_X, y = TAG_Y_B - TAG_DY, size = TAG_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label(pB_title, x = TAG_X + X_TTL, y = TAG_Y_B,
             size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = TAG_X + X_TTL, y = TAG_Y_B - SUB_OFFSET,
             size = SUB_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

graphics.off()
pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "MAIN_F07_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F07_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)
message("F07 composite saved: [A: AUC bars] / [B: 2x3 hero grid]")

# Build supplementary xlsx (12 sheets)

f07 <- function(p) file.path(DAT, p)

message("=== F07 supplementary workbook ===")
f07_specs <- list(
  list(name="module_grid_summary",             path=f07("module_grid/module_grid_summary.csv")),
  list(name="module_grid_curves",              path=f07("module_grid/module_grid_curves.csv")),
  list(name="panel_A_classifier_auc",          path=f07("panel_A_multi_classifier_auc.csv")),
  list(name="panel_A_feature_stability",       path=f07("panel_A_feature_stability.csv")),
  list(name="panel_A_permutation",             path=f07("panel_A_permutation.csv")),
  list(name="panel_A_roc_curves",              path=f07("panel_A_roc_curves.csv")),
  list(name="classifier_pilot_summary",               path=f07("classifier_pilot_summary.csv")),
  list(name="classifier_pilot_curves",                path=f07("classifier_pilot_curves.csv")),
  list(name="panel_B_full_screen",             path=f07("panel_B_full_screen_bh.csv")),
  list(name="loso_auc_summary",                path=f07("loso_auc/loso_auc_summary.csv")),
  list(name="loso_wgcna_refit_summary",        path=f07("loso_auc/loso_wgcna_refit_summary.csv")),
  list(name="loso_wgcna_refit_mod_stability",  path=f07("loso_auc/loso_wgcna_refit_module_stability.csv"))
)
build_workbook(
  f07("F07_supplementary.xlsx"),
  title = "F07 \u2014 Figure 7 source data",
  description = "Phenotype-prediction outputs: univariate module-outcome ROCs, multivariate classifiers, LOSO cross-validation, age-stratified module-phenotype coupling, and per-module Jaccard stability.",
  overview_df = data.frame(
    Sheet = c(
      "module_grid_summary",
      "module_grid_curves",
      "panel_A_classifier_auc",
      "panel_A_feature_stability",
      "panel_A_permutation",
      "panel_A_roc_curves",
      "classifier_pilot_summary",
      "classifier_pilot_curves",
      "panel_B_full_screen",
      "loso_auc_summary",
      "loso_wgcna_refit_summary",
      "loso_wgcna_refit_mod_stability"),
    Description = c(
      "Per-module univariate ROC summary: AUC, permutation p, BH q for each module-outcome pair",
      "Per-module ROC curves: TPR/FPR coordinates for plotting",
      "Multivariate classifier AUC: raw-protein (k=10), phenotype-only, ME-stack, and delta-ME classifiers",
      "Feature stability: top-10 raw-protein classifier features across LOSO folds",
      "Permutation null: observed vs permuted AUC distributions for each classifier",
      "Multivariate ROC curves: TPR/FPR coordinates for each classifier",
      "Classifier pilot summary: module eigengene ROCs for age discrimination (top modules)",
      "Classifier pilot curves: TPR/FPR coordinates for classifier pilot ROCs",
      "Panel B full screen: 180-test module-phenotype correlations (9 modules \u00d7 2 sources \u00d7 5 outcomes \u00d7 2 strata) with BH correction",
      "LOSO fixed-module: leave-one-subject-out AUCs using fixed full-sample module definitions",
      "LOSO + WGCNA refit: full-refit LOSO \u2014 top-12 pair AUCs with network refit on n\u22121 subjects per fold, training modules matched to full-sample by Jaccard",
      "LOSO module stability: per-full-sample-module mean/min Jaccard of training-fold vs full-sample assignments, + count of folds where best Jaccard fell below 0.5"),
    stringsAsFactors = FALSE),
  sheet_specs = f07_specs
)
cleanup_after_workbook(f07_specs)
