# Figure 6 — Main Composite: WGCNA Module-Trait Associations
#
# Panel A (heatmap with preservation) — left 65%
# Panel B (stacked NES scatters: concordance + reversal) — right 35%
# Panel B title/subtitle — drawn by stitcher, aligned with Panel A title
# Panel B legend — positioned at same y as Panel A color legends
# Panel letters — slightly above-left of their respective titles
#
# Sources panel scripts which save standalone PNGs, then composites via cowplot.
# Outputs: MAIN_F06_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(patchwork)
library(cowplot)
library(png)
library(grid)

RPT_PDF  <- "04_Figures/F06/b_reports/main/pdf"
RPT_PNG  <- "04_Figures/F06/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# Source panels — each saves standalone PNGs to b_reports/main/png/panels/
source("04_Figures/F06/a_script/main/panels/panel_A.R")
source("04_Figures/F06/a_script/main/panels/panel_B.R")

# Re-define output dirs AFTER sourcing panels (which overwrite RPT_PNG/RPT_PDF)
RPT_PDF  <- "04_Figures/F06/b_reports/main/pdf"
RPT_PNG  <- "04_Figures/F06/b_reports/main/png"

# Read panel PNGs for raster compositing
PANEL_DIR <- "04_Figures/F06/b_reports/main/png/panels"
read_panel <- function(file, dir = PANEL_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob     <- read_panel("MAIN_panel_A_heatmap.png")
pB_grob     <- read_panel("MAIN_panel_B_scatters.png")
pB_leg_grob <- read_panel("MAIN_panel_B_legend.png")

COMP_W <- 470
COMP_H <- 300
TAG_SZ <- 16
TITLE_SZ <- 14
SUBTITLE_SZ <- 10

# --- Panel A structure (proportions within the PNG) ---
# Title area:    y = 0.92 - 1.00
# Brackets:      y = 0.82 - 0.92
# Heatmap grid:  y = 0.22 - 0.82
# X-axis labels: y = 0.14 - 0.22
# Color legends:  y = 0.06 - 0.14

A_TITLE_Y  <- 0.835   # title baseline — just above the panel tops
A_LEGEND_Y <- 0.14    # raised — bring protein key closer to plots
GRID_TOP   <- 0.80    # top of heatmap grid
GRID_BOT   <- 0.22    # bottom of heatmap grid

# Panel B scatters: align with heatmap grid, moved left closer to A
B_X        <- 0.43    # closer to panel A edge
B_W        <- 0.55    # fill remaining space
B_TITLE_X  <- 0.636   # title/letter sit above the actual scatter plot area
B_GRID_H   <- GRID_TOP - GRID_BOT  # 0.58

# --- Independently tunable label positions (mm from composite origin) ---
mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

# Panel A labels (mm) — letter baseline-aligned with title
A_LET_X  <- 15;   A_LET_Y  <- 250
A_TTL_X  <- 58;   A_TTL_Y  <- 250
A_SUB_X  <- 58;   A_SUB_Y  <- 244

# Panel B labels (mm) — letter baseline-aligned with title
B_LET_X  <- 285;  B_LET_Y  <- 250
B_TTL_X  <- 294;  B_TTL_Y  <- 250
B_SUB_X  <- 294;  B_SUB_Y  <- 244

# Crop canvas to content bounds (removes white space margins)
CROP_L <-  0.01;  CROP_R <- 0.80
CROP_B <-  0.17;  CROP_T <- 0.85
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B drawn FIRST (behind) — its white left margin hides behind Panel A
  draw_grob(pB_grob, x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg_grob, x = 0.70, y = 0.195,
            width = 0.24, height = 0.04, hjust = 0.5, vjust = 0.5) +
  # Panel A drawn SECOND (on top) — covers Panel B's white left margin
  draw_grob(pA_grob, x = 0.01, y = 0, width = 0.60, height = 0.96,
            hjust = 0, vjust = 0) +
  # Panel A: letter, title, subtitle
  draw_label("A", x = mm2x(A_LET_X), y = mm2y(A_LET_Y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("WGCNA Module\u2013Trait Associations",
             x = mm2x(A_TTL_X), y = mm2y(A_TTL_Y),
             size = TITLE_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label("10 modules | LMM (BH) | Stratified r (BH per-trait)",
             x = mm2x(A_SUB_X), y = mm2y(A_SUB_Y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  # Panel B: letter, title, subtitle
  draw_label("B", x = mm2x(B_LET_X), y = mm2y(B_LET_Y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Module\u2013Level NES Scatters",
             x = mm2x(B_TTL_X), y = mm2y(B_TTL_Y),
             size = TITLE_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label("fGSEA on module-member t-stat ranks",
             x = mm2x(B_SUB_X), y = mm2y(B_SUB_Y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "MAIN_F06_composite.pdf"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F06_composite.png"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       dpi = 300, limitsize = FALSE)
message("F06 composite saved: MAIN_F06_composite.{pdf,png}")

# --- Supplementary Excel: WGCNA module-trait analysis + per-panel data ---
source("04_Figures/shared/figure_supplement_helpers.R")

f06 <- function(p) file.path("04_Figures/F06/c_data", p)

cat("=== F06 supplementary workbook ===\n")

# Convert F07-consumed matrix RDS files to Excel-embeddable data frames (2026-04-15).
.shared       <- readRDS(f06("shared_objects.rds"))
.MEs          <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))),      "sample_id")
.me_pre       <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))),   "subject_key")
.me_post      <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))),  "subject_key")
.delta_me     <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")
.common_subj  <- data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE)

f06_specs <- list(
  # --- Panel/supp outputs written by active panel + supp scripts ---
  list(name="panel_A_heatmap",             path=f06("01_panel_A_heatmap_data.csv")),
  list(name="panel_A_correlation_CIs",     path=f06("01_panel_A_correlation_CIs.csv")),
  list(name="panel_B_module_fgsea",        path=f06("panel_B_module_fgsea.csv")),
  list(name="panel_B_triptych_zscores",    path=f06("03_panel_B_heatmap_zscores.csv")),
  list(name="panel_B_triptych_eigengene",  path=f06("03_panel_B_eigengene_data.csv")),
  list(name="panel_B_triptych_enrichment", path=f06("03_panel_B_triptych_enrichment.csv")),
  list(name="panel_D_hub_network",         path=f06("04_panel_D_hub_network.csv")),
  list(name="panel_D_hub_CIs",             path=f06("04_panel_D_hub_CIs.csv")),
  list(name="panel_E_preservation",        path=f06("05_panel_E_preservation.csv")),
  list(name="SUPP_a01_dendrogram",         path=f06("asupp_B_QC_dendrogram_SUPP_data.csv")),
  list(name="SUPP_a02_bicor_sensitivity",  path=f06("supp/a02_bicor_sensitivity.csv")),
  list(name="SUPP_a03_compartment_enrichment", path=f06("supp/a03_compartment_enrichment.csv")),
  list(name="SUPP_a05_sft_fit_indices",    path=f06("supp/a05_sft_fit_indices.csv")),
  list(name="SUPP_e02_dep_module_counts",  path=f06("supp/e02_dep_module_counts.csv")),
  list(name="SUPP_e02_dep_module_proportions", path=f06("supp/e02_dep_module_proportions.csv")),
  list(name="SUPP_lmm_stratified",         path=f06("supp_F_lmm_stratified_audit.csv")),
  # --- WGCNA reference outputs from YvO_WGCNA_run.R (wgcna/ subdir) ---
  list(name="WGCNA_module_assignments",    path=f06("wgcna/wgcna_module_assignments.csv")),
  list(name="WGCNA_mod_bio_labels",        path=f06("mod_bio_labels.csv")),
  list(name="WGCNA_hub_proteins",          path=f06("wgcna/wgcna_hub_proteins.csv")),
  list(name="WGCNA_module_enrichment",     path=f06("wgcna/wgcna_module_enrichment.csv")),
  list(name="WGCNA_module_GO_enrichment",  path=f06("wgcna/wgcna_module_GO_enrichment.csv")),
  list(name="WGCNA_module_trait_cor",      path=f06("wgcna/wgcna_module_trait_correlations.csv")),
  list(name="WGCNA_module_trait_pval_bh",  path=f06("wgcna/wgcna_module_trait_pvalues_bh.csv")),
  list(name="WGCNA_baseline_trait_cor",         path=f06("wgcna/wgcna_baseline_trait_correlations.csv")),
  list(name="WGCNA_baseline_trait_cor_young",   path=f06("wgcna/wgcna_baseline_trait_correlations_young.csv")),
  list(name="WGCNA_baseline_trait_cor_old",     path=f06("wgcna/wgcna_baseline_trait_correlations_old.csv")),
  list(name="WGCNA_baseline_pval_bh",         path=f06("wgcna/wgcna_baseline_trait_pvalues_bh.csv")),
  list(name="WGCNA_baseline_pval_bh_young",   path=f06("wgcna/wgcna_baseline_trait_pvalues_bh_young.csv")),
  list(name="WGCNA_baseline_pval_bh_old",     path=f06("wgcna/wgcna_baseline_trait_pvalues_bh_old.csv")),
  list(name="WGCNA_baseline_pval_raw_young",  path=f06("wgcna/wgcna_baseline_trait_pvalues_raw_young.csv")),
  list(name="WGCNA_baseline_pval_raw_old",    path=f06("wgcna/wgcna_baseline_trait_pvalues_raw_old.csv")),
  list(name="WGCNA_change_trait_cor",         path=f06("wgcna/wgcna_change_trait_correlations.csv")),
  list(name="WGCNA_change_trait_cor_young",   path=f06("wgcna/wgcna_change_trait_correlations_young.csv")),
  list(name="WGCNA_change_trait_cor_old",     path=f06("wgcna/wgcna_change_trait_correlations_old.csv")),
  list(name="WGCNA_change_pval_bh",           path=f06("wgcna/wgcna_change_trait_pvalues_bh.csv")),
  list(name="WGCNA_change_pval_bh_young",     path=f06("wgcna/wgcna_change_trait_pvalues_bh_young.csv")),
  list(name="WGCNA_change_pval_bh_old",       path=f06("wgcna/wgcna_change_trait_pvalues_bh_old.csv")),
  list(name="WGCNA_change_pval_raw_young",    path=f06("wgcna/wgcna_change_trait_pvalues_raw_young.csv")),
  list(name="WGCNA_change_pval_raw_old",      path=f06("wgcna/wgcna_change_trait_pvalues_raw_old.csv")),
  list(name="WGCNA_sft_summary",           path=f06("wgcna/wgcna_sft_summary.csv")),
  list(name="WGCNA_lmm_contrast_audit",    path=f06("wgcna/wgcna_lmm_contrast_audit.csv")),
  list(name="WGCNA_lmm_stratified_audit",  path=f06("wgcna/wgcna_lmm_stratified_audit.csv")),
  list(name="WGCNA_gs_phenotype_choices",  path=f06("wgcna/gs_phenotype_choices.csv")),
  # --- Metadata ---
  list(name="metadata_samples",            path=f06("meta.csv")),
  list(name="metadata_subj_age",           path=f06("subj_age.csv")),
  list(name="metadata_pheno_wide",         path=f06("pheno_wide.csv")),
  list(name="metadata_imp_annotations",    path=f06("imp_annotations.csv")),
  # --- ME matrices (from RDS -> data.frame for Excel) ---
  list(name="MEs",         df=.MEs),
  list(name="me_pre",      df=.me_pre),
  list(name="me_post",     df=.me_post),
  list(name="delta_me",    df=.delta_me),
  list(name="common_subj", df=.common_subj)
)

build_workbook(
  f06("F06_supplementary.xlsx"),
  title = "F06 \u2014 Figure 6 source data",
  description = "WGCNA module assignments, hub proteins, pathway enrichment, trait correlations, eigengene metadata (MEs, me_pre, me_post, \u0394me).",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_heatmap", "panel_A_correlation_CIs",
      "panel_B_module_fgsea",
      "panel_B_triptych_zscores", "panel_B_triptych_eigengene",
      "panel_B_triptych_enrichment",
      "panel_D_hub_network", "panel_D_hub_CIs",
      "panel_E_preservation",
      "SUPP_a01_dendrogram", "SUPP_a02_bicor_sensitivity",
      "SUPP_a03_compartment_enrichment", "SUPP_a05_sft_fit_indices",
      "SUPP_e02_dep_module_counts", "SUPP_e02_dep_module_proportions",
      "SUPP_lmm_stratified",
      "WGCNA_module_assignments", "WGCNA_mod_bio_labels", "WGCNA_hub_proteins",
      "WGCNA_module_enrichment", "WGCNA_module_GO_enrichment",
      "WGCNA_module_trait_cor", "WGCNA_module_trait_pval_bh",
      "WGCNA_baseline_trait_cor", "WGCNA_baseline_trait_cor_young", "WGCNA_baseline_trait_cor_old",
      "WGCNA_baseline_pval_bh", "WGCNA_baseline_pval_bh_young", "WGCNA_baseline_pval_bh_old",
      "WGCNA_baseline_pval_raw_young", "WGCNA_baseline_pval_raw_old",
      "WGCNA_change_trait_cor", "WGCNA_change_trait_cor_young", "WGCNA_change_trait_cor_old",
      "WGCNA_change_pval_bh", "WGCNA_change_pval_bh_young", "WGCNA_change_pval_bh_old",
      "WGCNA_change_pval_raw_young", "WGCNA_change_pval_raw_old",
      "WGCNA_sft_summary", "WGCNA_lmm_contrast_audit", "WGCNA_lmm_stratified_audit",
      "WGCNA_gs_phenotype_choices",
      "metadata_samples", "metadata_subj_age",
      "metadata_pheno_wide", "metadata_imp_annotations",
      "MEs", "me_pre", "me_post", "delta_me", "common_subj"),
    Description = c(
      "Panel A: 18-column heatmap data (modules x traits)",
      "Panel A: bootstrap CIs for heatmap correlations",
      "Panel B: per-module fGSEA collapsed pathway list (NES scatter input)",
      "SUPP triptych: heatmap z-scores for key modules",
      "SUPP triptych: per-module ME eigengene long-format (Pre vs Post by age)",
      "SUPP triptych: per-module enrichment results",
      "SUPP D: hub protein network edges per key module",
      "SUPP D: bootstrap CIs for hub kME",
      "SUPP E: module preservation Zsummary scores (a04 output)",
      "SUPP a01: sample/module dendrogram for QC composite",
      "SUPP a02: bicor vs Pearson sensitivity heatmap data",
      "SUPP a03: subcellular compartment enrichment per module",
      "SUPP a05: SFT fit indices per power",
      "SUPP e02: DEP-module counts per contrast",
      "SUPP e02: DEP-module proportions",
      "SUPP F: LMM stratified audit (legacy preservation)",
      "WGCNA: protein -> module color + gene symbol",
      "WGCNA: module ID, color, biological label, n_proteins",
      "WGCNA: top-10 hub proteins per module by kME",
      "WGCNA: multi-database ORA per module (H+KEGG+Reactome+GO:BP)",
      "WGCNA: GO:BP enrichment per module",
      "WGCNA: combined Pearson r (modules x traits)",
      "WGCNA: BH-corrected p-values (modules x traits)",
      "WGCNA: baseline (Pre) trait correlations - combined cohort",
      "WGCNA: baseline (Pre) trait correlations - Young only",
      "WGCNA: baseline (Pre) trait correlations - Old only",
      "WGCNA: baseline BH-corrected p-values - combined",
      "WGCNA: baseline BH-corrected p-values - Young only",
      "WGCNA: baseline BH-corrected p-values - Old only",
      "WGCNA: baseline raw p-values - Young only",
      "WGCNA: baseline raw p-values - Old only",
      "WGCNA: change (Post-Pre) trait correlations - combined cohort",
      "WGCNA: change trait correlations - Young only",
      "WGCNA: change trait correlations - Old only",
      "WGCNA: change BH-corrected p-values - combined",
      "WGCNA: change BH-corrected p-values - Young only",
      "WGCNA: change BH-corrected p-values - Old only",
      "WGCNA: change raw p-values - Young only",
      "WGCNA: change raw p-values - Old only",
      "WGCNA: soft-threshold power selection (R^2, connectivity)",
      "WGCNA: LMM contrast audit (lmer + emmeans + KR df, BH)",
      "WGCNA: LMM age-stratified audit",
      "WGCNA: gene-significance phenotype choices (panel_D_hub config)",
      "Metadata: sample-level (Col_ID, Group, Timepoint, Subject_ID)",
      "Metadata: subject -> age group mapping",
      "Metadata: per-subject phenotype data (baseline + change for all traits)",
      "Metadata: imputed-matrix protein annotations (uniprot -> gene + description)",
      "Module eigengenes per sample (rows = sample_id, cols = MEturquoise ... MEgrey)",
      "Pre-timepoint module eigengenes per subject (for F07 readers)",
      "Post-timepoint module eigengenes per subject (for F07 readers)",
      "Change in module eigengenes (Post - Pre) per subject (for F07 readers)",
      "Common subject keys across Pre/Post for subject-paired analyses"),
    stringsAsFactors = FALSE),
  sheet_specs = f06_specs
)
cleanup_after_workbook(f06_specs,
  extra_subdirs = character(0),
  preserve_patterns = c(
    "^00_input/", "^01_normalization/", "^02_Imputation/", "^03_DEP/",
    "^04_Figures/shared/",
    "^04_Figures/F06/c_data/wgcna/",
    "^04_Figures/F06/c_data/.*\\.rds$",
    "^04_Figures/F06/c_data/key_modules\\.txt$",
    "^04_Figures/F06/c_data/mod_bio_labels\\.csv$",
    "^04_Figures/F06/c_data/meta\\.csv$",
    "^04_Figures/F06/c_data/imp_annotations\\.csv$",
    "^04_Figures/F06/c_data/subj_age\\.csv$",
    "^04_Figures/F06/c_data/pheno_wide\\.csv$"))
