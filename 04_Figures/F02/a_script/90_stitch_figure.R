# F02 — Main Figure Composite (6-panel 3×2 grid)
# Sources 6 panel scripts and composes into a 3-column x 2-row layout.
#   Top:    A (PCA)    | B (Effect size) | C (DEPs per contrast)
#   Bottom: D (UpSet)  | E (Pathways)    | F (DEP rank location)
#
# Uniform column widths across rows guarantee:
#   - C and F align vertically (same right column)
#   - E right edge = B right edge (same middle column right edge)
#   - D left edge = A left edge (same left column left edge)
#
# pA, pB, pC, pF are raw ggplots → patchwork aligns their plot regions.
# pD (ggdraw/cowplot) and pE (patchwork with inset_element) → wrap_elements().
# Manual cowplot tag placement for consistent letters across raw + wrapped.
# Outputs: F02_composite_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

# Source panels — each creates a named plot object and saves standalone PDFs
source("04_Figures/F02/a_script/panels/panel_A.R")  # -> pA (ggplot, PCA)
source("04_Figures/F02/a_script/panels/panel_B.R")  # -> pB (ggplot, effect size)
source("04_Figures/F02/a_script/panels/panel_C.R")  # -> pC (ggplot, DEPs per contrast)
source("04_Figures/F02/a_script/panels/panel_D.R")  # -> pD (ggdraw, UpSet)
source("04_Figures/F02/a_script/panels/panel_E.R")  # -> pE (patchwork, Pathways)
source("04_Figures/F02/a_script/panels/panel_F.R")  # -> pF (ggplot, DEP rank location)

RPT_DIR <- "04_Figures/F02/b_reports"
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F02"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript_2026-04-16/02_Figures/F02_proteome_DEP"
dir.create(WRITING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BOX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BOX_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# Strip individual tags from raw ggplots — cowplot draw_label adds them uniformly
pA_clean <- pA + labs(tag = NULL)
pB_clean <- pB + labs(tag = NULL)
pC_clean <- pC + labs(tag = NULL)
pF_clean <- pF + labs(tag = NULL)

# 3x2 design with tiny spacer row between the two panel rows:
#   ABC   (row 1)
#   ###   (spacer)
#   DEF   (row 2)
layout <- "ABC\n###\nDEF"

ROW_TOP <- 0.485
SPACER  <- 0.03

composite <- pA_clean + pB_clean + pC_clean +
             wrap_elements(full = pD) +
             wrap_elements(full = pE) +
             pF_clean +
  plot_layout(
    design  = layout,
    # Middle column widened from 95 → 105 to expand E's plot panel and
    # give B/E more breathing room; A/C shrunk from 165 → 160 to keep
    # composite width constant at 425 mm of content.
    widths  = c(160, 105, 160),
    heights = c(ROW_TOP, SPACER, 1 - ROW_TOP - SPACER)
  ) &
  # F01-proportional title/subtitle sizes (see composite_text_sizes in shared/style.R)
  theme(plot.title    = element_text(size = composite_text_sizes(260)$title),
        plot.subtitle = element_text(size = composite_text_sizes(260)$subtitle))

# --- Manual tag placement via cowplot ---
# Column left-edge positions (normalized to composite width 440 mm, outer margin ~7.5 mm each side):
#   A cell left ≈ 7.5/440  ≈ 0.017        → X_LEFT  = 0.020 (1.3 mm past A cell left)
#   B cell left ≈ 167.5/440 ≈ 0.381       → X_MID   = 0.384 (1.3 mm past B cell left)
#   C cell left ≈ 272.5/440 ≈ 0.619       → X_RIGHT = 0.622 (1.3 mm past C cell left)
# Row top positions (normalized to total composite height):
#   Top row (A, B, C)    ≈ 0.995
#   Bottom row (D, E, F) ≈ 1 - ROW_TOP - SPACER + 0.005 ≈ 0.490
TAG_SZ  <- composite_text_sizes(260)$tag   # F01-proportional (15pt @ 155mm H -> 25pt @ 260mm H)
TOP_Y   <- 0.995 - 4/260                   # -4 mm tag shift (user-tuned)
BOT_Y   <- 1 - ROW_TOP - SPACER + 0.005 - 4/260
X_LEFT  <- 0.020
X_MID   <- 0.384
# C/F letters consistent ~1.3 mm past C cell left, matching A/D offset.
X_RIGHT <- 0.622

composite <- ggdraw(composite) +
  draw_label("A", x = X_LEFT,  y = TOP_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = X_MID,   y = TOP_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = X_RIGHT, y = TOP_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = X_LEFT,  y = BOT_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("E", x = X_MID,   y = BOT_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("F", x = X_RIGHT, y = BOT_Y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

# Composite size:
#   W: 160 (A/D) + 105 (B/E) + 160 (C/F) + margin ≈ 440 mm
#   H: 120 (top row) + 120 (bottom row) + spacer + margin ≈ 260 mm
COMP_W <- 440
COMP_H <- 260

ggsave(file.path(RPT_DIR, "MAIN_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "MAIN_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
file.copy(file.path(RPT_DIR, "MAIN_F02_composite.pdf"), file.path(WRITING_DIR, "01_main_composite_f02.pdf"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F02_composite.png"), file.path(WRITING_DIR, "01_main_composite_f02.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F02_composite.png"), file.path(BOX_DIR, "F02_main.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F02_composite.pdf"), file.path(BOX_DIR, "pdf", "F02_main.pdf"), overwrite = TRUE)

message("F02 composite (3x2, sourced panels, aligned) saved")

# --- Supplementary Excel: one sheet per panel (panel scripts already wrote CSVs) ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F02 supplementary workbook ===\n")
f02_specs <- list(
  list(name = "panel_A_pca_variance",   path = "04_Figures/F02/c_data/audit_panel_A_pca_variance_ci.csv"),
  list(name = "panel_A_betadisper",     path = "04_Figures/F02/c_data/audit_panel_A_betadisper.csv"),
  list(name = "panel_B_logfc_long",     path = "04_Figures/F02/c_data/panel_B_logfc_long.csv"),
  list(name = "panel_B_stats",          path = "04_Figures/F02/c_data/panel_B_stats.csv"),
  list(name = "panel_C_dep_fraction",   path = "04_Figures/F02/c_data/audit_panel_C_dep_fraction_ci.csv"),
  list(name = "panel_D_upset",          path = "04_Figures/F02/c_data/panel_D_upset_overlap_enrichment.csv"),
  list(name = "panel_E_fgsea_sig",      path = "04_Figures/F02/c_data/panel_E_fgsea_sig.csv"),
  list(name = "panel_E_fgsea_counts",   path = "04_Figures/F02/c_data/panel_E_fgsea_counts.csv"),
  list(name = "panel_F_barcode",        path = "04_Figures/F02/c_data/panel_F_barcode_enrichment.csv"),
  list(name = "SUPP_panel_C_cv_scatter", path = "04_Figures/F02/c_data/audit_panel_C_cv_scatter.csv"),
  list(name = "SUPP_panel_C_cv_delta",   path = "04_Figures/F02/c_data/audit_panel_C_cv_delta.csv"),
  list(name = "SUPP_panel_D_median_cv",  path = "04_Figures/F02/c_data/audit_panel_D_median_cv_ci.csv"),
  list(name = "SUPP_panel_D_wilcoxon",   path = "04_Figures/F02/c_data/audit_panel_D_wilcoxon_effects.csv"),
  list(name = "SUPP_panel_E_imputed",    path = "04_Figures/F02/c_data/audit_panel_E_imputed.csv"),
  list(name = "SUPP_panel_E_wilcoxon",   path = "04_Figures/F02/c_data/audit_panel_E_wilcoxon.csv")
)
build_workbook(
  "04_Figures/F02/c_data/F02_supplementary.xlsx",
  title = "F02 \u2014 Figure 2 source data",
  description = "Proteome + DEP overview: CV distributions, logFC density, PCA + PERMANOVA, DEP counts, UpSet intersections, fGSEA enrichment, rank barcode.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_pca_variance", "panel_A_betadisper",
      "panel_B_logfc_long", "panel_B_stats",
      "panel_C_dep_fraction",
      "panel_D_upset",
      "panel_E_fgsea_sig", "panel_E_fgsea_counts",
      "panel_F_barcode",
      "SUPP_panel_C_cv_scatter", "SUPP_panel_C_cv_delta",
      "SUPP_panel_D_median_cv",  "SUPP_panel_D_wilcoxon",
      "SUPP_panel_E_imputed",    "SUPP_panel_E_wilcoxon"),
    Description = c(
      "Panel A: PC1/PC2 variance explained with 95% bootstrap CI",
      "Panel A: dispersion homogeneity test (betadisper) per factor",
      "Panel B: long-format gene x contrast x logFC (all 4 contrasts)",
      "Panel B: per-contrast median |logFC| with bootstrap 95% CI + n_above_05",
      "Panel C: per-contrast DEP % of proteome with binomial CI",
      "Panel D: pairwise Fisher exact on Pi-score DEP set overlaps (BH)",
      "Panel E: significant fGSEA pathways (padj < 0.05) per contrast x database",
      "Panel E: significant fGSEA counts per contrast x database x direction (Up/Down)",
      "Panel F: DEP counts per contrast x direction (barcode)",
      "SUPP Panel C: per-protein CV% Pre vs Post + dCV Young vs Old",
      "SUPP Panel C: dCV per-protein Young vs Old long-format",
      "SUPP Panel D: median CV% by age x time + bootstrap CI",
      "SUPP Panel D: Wilcoxon effects on inter-individual CV%",
      "SUPP Panel E: per-subject intra-individual variability (imputed) + imputation fraction",
      "SUPP Panel E: Wilcoxon effects on intra-individual variability"),
    stringsAsFactors = FALSE),
  sheet_specs = f02_specs
)
file.copy("04_Figures/F02/c_data/F02_supplementary.xlsx", file.path(WRITING_DIR, "04_supp_data_f02.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F02/c_data/F02_supplementary.xlsx", file.path(BOX_DIR, "F02_source_data.xlsx"), overwrite = TRUE)
cleanup_after_workbook(f02_specs)
