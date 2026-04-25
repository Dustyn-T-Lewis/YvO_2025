# Figure 7 — Composite: Module Eigengene Predictive Power
#
# Capstone final structure (2026-04-14): 2 panels, both full-width.
#   A (top) : Baseline Module Eigengenes Classify Age
#             — AUC bars for 4 classifiers + per-bar ROC sparklines
#   B (bot) : Module – Phenotype Coupling Is Age-Dependent
#             — 2x3 hero scatter grid of top age-flipped/significant signals
#
# PCA + biplot demoted to SUPP (SUPP_pca_biplot.png) as redundant with A.
#
# Generates: MAIN_F07_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(patchwork)
library(cowplot)

# Source ALL data-producing panel scripts FIRST so every CSV needed by the
# xlsx build is present before cleanup_after_workbook() runs. The SUPP
# stitcher then becomes a pure composite-assembler (reads pre-rendered PNGs).
# panel_A.R hard-depends on c_data/module_grid/*.csv (written by
# SUPP_module_grid.R), so SUPP_module_grid.R must source before panel_A.R.
# Panel scripts each reassign `RPT` to panels/ subdir, so the stitcher's
# paths must be defined AFTER all sources.
source("04_Figures/F07/a_script/supp/panels/SUPP_module_grid.R")           # c_data/module_grid/*.csv + SUPP_F07_module_grid.{pdf,png}
source("04_Figures/F07/a_script/main/panels/panel_A.R")                    # MAIN_panel_A_auc.{png,pdf} + panel_A_*.csv
source("04_Figures/F07/a_script/main/panels/panel_B.R")                    # MAIN_panel_B_hero_grid.png + panel_B_full_screen_bh.csv
source("04_Figures/F07/a_script/supp/panels/SUPP_panel_B_grid.R")          # SUPP_F07_panel_B_grid.{pdf,png}
source("04_Figures/F07/a_script/supp/panels/prepare_roc_data.R")            # roc_pilot_{summary,curves}.csv
source("04_Figures/F07/a_script/supp/panels/SUPP_roc_panel.R")             # SUPP_F07_roc_panel.{pdf,png} + SUPP_*.csv
source("04_Figures/F07/a_script/supp/panels/SUPP_multivariate_classifier.R")  # SUPP_F07_multivariate_classifier.{pdf,png}
source("04_Figures/F07/a_script/supp/panels/SUPP_loso_sensitivity.R")          # loso_auc/loso_auc_summary.csv + SUPP_F07_loso_sensitivity.{pdf,png}
source("04_Figures/F07/a_script/supp/panels/SUPP_loso_wgcna_refit.R")          # loso_auc/loso_wgcna_refit_{summary,module_stability}.csv + SUPP_F07_loso_wgcna_refit.{pdf,png}

RPT_PNG <- "04_Figures/F07/b_reports/main/png"
RPT_PDF <- "04_Figures/F07/b_reports/main/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# --- Object-based composition (uses ggplot objects exported by panel scripts) ---
# Panel A ggsave dimensions: W=178, H=75  → aspect 178/75 = 2.373
# Panel B ggsave dimensions: W=178, H=90  → aspect 178/90 = 1.978
COMP_W  <- 210
A_H     <- 75                        # from panel_A.R ggsave height
B_H     <- 90                        # from panel_B.R ggsave height
TAG_A   <- 10                        # title strip above panel A
TAG_B   <- 12                        # title strip above panel B (room for title+subtitle+gap)
COMP_H  <- TAG_A + A_H + TAG_B + B_H  # 187 mm

composite <- (plot_spacer() / pA / plot_spacer() / pB) +
  plot_layout(heights = c(TAG_A, A_H, TAG_B, B_H))

txt    <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag + 3
TTL_SZ <- txt$title + 3
SUB_SZ <- txt$subtitle + 2.5

# Tag positions sit inside the dedicated top strips above each panel.
TAG_X   <- 4 / COMP_W
X_TTL   <- 0.040          # title starts right of tag (F02 convention)
TAG_DY  <- -0.002         # baseline-align tag with smaller title
SUB_OFFSET <- 0.016 + (1 / COMP_H)  # title-to-subtitle gap + 1mm down
DY5     <- 5 / COMP_H     # 5mm downward shift

TAG_Y_A <- 0.995 - DY5 + (3 / COMP_H)
TAG_Y_B <- (B_H + TAG_B) / COMP_H + (TAG_B / COMP_H) * 0.85 - DY5 * 0.4 - (5 / COMP_H) - (10 / COMP_H) + (4 / COMP_H)

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

graphics.off()  # close stale devices from sourced panels; cairo probe needs clean state
pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "MAIN_F07_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F07_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)
message("F07 composite saved: [A: AUC bars] / [B: 2x3 hero grid]")

# --- Supplementary Excel: one workbook, sheets keyed to figure panels ---
source("04_Figures/shared/figure_supplement_helpers.R")

f07 <- function(p) file.path("04_Figures/F07/c_data", p)

cat("=== F07 supplementary workbook ===\n")
# Spec list matches the current set of panel scripts sourced above. Stale entries
# for the 2026-04-15-superseded 4-panel SUPP and pre-redesign panel_C/D were
# removed to stop noisy SKIP lines.
f07_specs <- list(
  list(name="panel_A_classifier_auc",    path=f07("panel_A_multi_classifier_auc.csv")),
  list(name="panel_A_feature_stability", path=f07("panel_A_feature_stability.csv")),
  list(name="panel_A_permutation",       path=f07("panel_A_permutation.csv")),
  list(name="panel_A_roc_curves",        path=f07("panel_A_roc_curves.csv")),
  list(name="panel_B_full_screen_bh",    path=f07("panel_B_full_screen_bh.csv")),
  list(name="module_grid_summary",       path=f07("module_grid/module_grid_summary.csv")),
  list(name="module_grid_curves",        path=f07("module_grid/module_grid_curves.csv")),
  list(name="roc_pilot_summary",         path=f07("roc_pilot_summary.csv")),
  list(name="roc_pilot_curves",          path=f07("roc_pilot_curves.csv")),
  list(name="loso_auc_summary",          path=f07("loso_auc/loso_auc_summary.csv")),
  list(name="loso_wgcna_refit_summary",  path=f07("loso_auc/loso_wgcna_refit_summary.csv")),
  list(name="loso_wgcna_refit_mod_stability", path=f07("loso_auc/loso_wgcna_refit_module_stability.csv"))
)
build_workbook(
  f07("F07_supplementary.xlsx"),
  title = "F07 \u2014 Figure 7 source data",
  description = "Phenotype-prediction outputs: classifier AUC, feature stability, ROC curves, full module\u2013phenotype coupling screen, LOSO sensitivity.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_classifier_auc", "panel_A_feature_stability", "panel_A_permutation",
      "panel_A_roc_curves", "panel_B_full_screen_bh",
      "module_grid_summary", "module_grid_curves",
      "roc_pilot_summary", "roc_pilot_curves",
      "loso_auc_summary",
      "loso_wgcna_refit_summary",
      "loso_wgcna_refit_mod_stability"),
    Description = c(
      "Panel A: AUC + 95% CI per classifier (multivariate, baseline MEs)",
      "Panel A: module selection frequency across LOOCV folds",
      "Panel A: permutation p-values per classifier",
      "Panel A: ROC curve coordinates (sensitivity, specificity, threshold)",
      "Panel B: full screen Pearson r for module-trait pairs (BH-corrected)",
      "SUPP module grid: per-module classifier AUC summary",
      "SUPP module grid: per-module ROC curves",
      "ROC pilot: summary stats for 7 pilot classifiers (AUC + CI + perm_p)",
      "ROC pilot: per-classifier ROC curve coordinates",
      "LOSO sensitivity: in-sample vs LOSO eigengene-projection AUC for top-12 module-outcome pairs (module defs held fixed; addresses projection optimism, not module-definition circularity)",
      "LOSO + WGCNA refit: full-refit LOSO — top-12 pair AUCs with network refit on n-1 subjects per fold, training modules matched to full-sample by Jaccard; addresses module-definition circularity",
      "LOSO module stability: per-full-sample-module mean/min Jaccard of training-fold vs full-sample assignments, + count of folds where best Jaccard fell below 0.5"),
    stringsAsFactors = FALSE),
  sheet_specs = f07_specs
)
cleanup_after_workbook(f07_specs,
  extra_subdirs = c("04_Figures/F07/c_data/module_grid"))
