# Figure 6 — Main Composite: WGCNA Module-Trait Associations
#
# Panel A (heatmap with preservation) — left 65%
# Panel B (stacked NES scatters: concordance + reversal) — right 35%
# Panel B title/subtitle — drawn by stitcher, aligned with Panel A title
# Panel B legend — positioned at same y as Panel A color legends
# Panel letters — slightly above-left of their respective titles
#
# Sources panel scripts which save standalone PNGs, then composites via cowplot.
# Outputs: MAIN_F06_composite.{pdf,png} + F06_supplementary.xlsx

library(here)
source(here::here("04_Figures", "shared", "style.R"))

library(patchwork)
library(cowplot)
library(png)
library(grid)

BASE <- here::here("04_Figures", "F06")

RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# Source panels — each saves standalone PNGs to b_reports/main/png/panels/
source(here::here("04_Figures", "F06", "a_script", "_panel_A.R"))
source(here::here("04_Figures", "F06", "a_script", "_panel_B.R"))

# Re-define output dirs AFTER sourcing panels (which overwrite RPT_PNG/RPT_PDF)
RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")

# Read panel PNGs for raster compositing
PANEL_DIR <- file.path(BASE, "b_reports", "main", "png", "panels")
read_panel <- function(file, dir = PANEL_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob     <- read_panel("MAIN_panel_A_heatmap.png")
pB_grob     <- read_panel("MAIN_panel_B_scatters.png")
pB_leg_grob <- read_panel("MAIN_panel_B_legend.png")

# ══════════════════════════════════════════════════════════════
# LAYOUT CONFIG — all positioning constants in one place
# ══════════════════════════════════════════════════════════════

layout_cfg <- list(
  # Canvas dimensions (mm) — internal drawing surface before crop
  comp_w = 470, comp_h = 300,
  # Text sizes (pt)
  tag_sz      = 14,
  title_sz    = 12,
  subtitle_sz = 9,
  # Panel A structure (proportions within the PNG)
  a_title_y   = 0.835,
  a_legend_y  = 0.14,
  grid_top    = 0.80,
  grid_bot    = 0.22,
  # Panel B scatters: position + width in normalized coords
  b_x         = 0.43,
  b_w         = 0.55,
  b_title_x   = 0.636,
  # Panel B legend
  b_leg_x = 0.70, b_leg_y = 0.195, b_leg_w = 0.24, b_leg_h = 0.04,
  # Panel A grob position
  a_grob_x = 0.01, a_grob_y = 0, a_grob_w = 0.60, a_grob_h = 0.96,
  # Label positions (mm from composite origin)
  a_let_x = 15,  a_let_y = 248,
  a_ttl_x = 58,  a_ttl_y = 248,
  a_sub_x = 58,  a_sub_y = 243,
  b_let_x = 285, b_let_y = 248,
  b_ttl_x = 294, b_ttl_y = 248,
  b_sub_x = 294, b_sub_y = 243,
  # Crop canvas to content bounds (removes white space margins)
  crop_l = 0.01, crop_r = 0.80,
  crop_b = 0.17, crop_t = 0.85
)

COMP_W <- layout_cfg$comp_w
COMP_H <- layout_cfg$comp_h
TAG_SZ      <- layout_cfg$tag_sz
TITLE_SZ    <- layout_cfg$title_sz
SUBTITLE_SZ <- layout_cfg$subtitle_sz

# --- Panel A structure (proportions within the PNG) ---
A_TITLE_Y  <- layout_cfg$a_title_y
A_LEGEND_Y <- layout_cfg$a_legend_y
GRID_TOP   <- layout_cfg$grid_top
GRID_BOT   <- layout_cfg$grid_bot

# Panel B scatters: align with heatmap grid, moved left closer to A
B_X      <- layout_cfg$b_x
B_W      <- layout_cfg$b_w
B_TITLE_X <- layout_cfg$b_title_x
B_GRID_H <- GRID_TOP - GRID_BOT  # 0.58

# --- Independently tunable label positions (mm from composite origin) ---
mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

# Panel A labels (mm)
A_LET_X <- layout_cfg$a_let_x; A_LET_Y <- layout_cfg$a_let_y
A_TTL_X <- layout_cfg$a_ttl_x; A_TTL_Y <- layout_cfg$a_ttl_y
A_SUB_X <- layout_cfg$a_sub_x; A_SUB_Y <- layout_cfg$a_sub_y

# Panel B labels (mm)
B_LET_X <- layout_cfg$b_let_x; B_LET_Y <- layout_cfg$b_let_y
B_TTL_X <- layout_cfg$b_ttl_x; B_TTL_Y <- layout_cfg$b_ttl_y
B_SUB_X <- layout_cfg$b_sub_x; B_SUB_Y <- layout_cfg$b_sub_y

# Crop canvas to content bounds (removes white space margins)
CROP_L <- layout_cfg$crop_l; CROP_R <- layout_cfg$crop_r
CROP_B <- layout_cfg$crop_b; CROP_T <- layout_cfg$crop_t
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B drawn FIRST (behind) — its white left margin hides behind Panel A
  draw_grob(pB_grob, x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg_grob, x = layout_cfg$b_leg_x, y = layout_cfg$b_leg_y,
            width = layout_cfg$b_leg_w, height = layout_cfg$b_leg_h,
            hjust = 0.5, vjust = 0.5) +
  # Panel A drawn SECOND (on top) — covers Panel B's white left margin
  draw_grob(pA_grob, x = layout_cfg$a_grob_x, y = layout_cfg$a_grob_y,
            width = layout_cfg$a_grob_w, height = layout_cfg$a_grob_h,
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
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

DAT <- file.path(BASE, "c_data")
f06 <- function(p) file.path(DAT, p)

message("=== F06 supplementary workbook ===")

# Convert F07-consumed matrix RDS files to Excel-embeddable data frames (2026-04-15).
.shared       <- readRDS(f06("shared_objects.rds"))
.MEs          <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))),      "sample_id")
.me_pre       <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))),   "subject_key")
.me_post      <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))),  "subject_key")
.delta_me     <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")
.common_subj  <- data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE)

f06_specs <- list(
  # --- Panel outputs (heatmap + module fGSEA) ---
  list(name="panel_A_heatmap",             path=f06("01_panel_A_heatmap_data.csv")),
  list(name="panel_B_module_fgsea",        path=f06("panel_B_module_fgsea.csv")),
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
  list(name="WGCNA_lmm_contrast_check",    path=f06("wgcna/wgcna_lmm_contrast_check.csv")),
  list(name="WGCNA_lmm_stratified_check",  path=f06("wgcna/wgcna_lmm_stratified_check.csv")),
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
      "panel_A_heatmap",
      "panel_B_module_fgsea",
      "WGCNA_module_assignments", "WGCNA_mod_bio_labels", "WGCNA_hub_proteins",
      "WGCNA_module_enrichment", "WGCNA_module_GO_enrichment",
      "WGCNA_module_trait_cor", "WGCNA_module_trait_pval_bh",
      "WGCNA_baseline_trait_cor", "WGCNA_baseline_trait_cor_young", "WGCNA_baseline_trait_cor_old",
      "WGCNA_baseline_pval_bh", "WGCNA_baseline_pval_bh_young", "WGCNA_baseline_pval_bh_old",
      "WGCNA_baseline_pval_raw_young", "WGCNA_baseline_pval_raw_old",
      "WGCNA_change_trait_cor", "WGCNA_change_trait_cor_young", "WGCNA_change_trait_cor_old",
      "WGCNA_change_pval_bh", "WGCNA_change_pval_bh_young", "WGCNA_change_pval_bh_old",
      "WGCNA_change_pval_raw_young", "WGCNA_change_pval_raw_old",
      "WGCNA_sft_summary", "WGCNA_lmm_contrast_check", "WGCNA_lmm_stratified_check",
      "WGCNA_gs_phenotype_choices",
      "metadata_samples", "metadata_subj_age",
      "metadata_pheno_wide", "metadata_imp_annotations",
      "MEs", "me_pre", "me_post", "delta_me", "common_subj"),
    Description = c(
      "Panel A: 18-column heatmap data (modules x traits)",
      "Panels B & C: per-module fGSEA pathway NES for Training-Young, Training-Old, and Aging contrasts (source data for both NES scatter panels)",
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
      "WGCNA: LMM contrast check (lmer + emmeans + KR df, BH)",
      "WGCNA: LMM age-stratified check",
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

message("=== F06 main panels + xlsx complete ===")
