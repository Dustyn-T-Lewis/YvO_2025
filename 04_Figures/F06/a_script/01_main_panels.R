# Figure 6 — WGCNA module-trait associations composite.

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures/shared/style.R")

library(patchwork)
library(cowplot)
library(png)
library(grid)
library(dplyr)
library(tidyr)

BASE <- "04_Figures/F06"

RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

source("04_Figures/F06/a_script/_panel_A_module_heatmap.R")
source("04_Figures/F06/a_script/_panel_B_nes_scatters.R")

# Re-define after panel scripts overwrite these vars
RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")

PANEL_DIR <- file.path(BASE, "b_reports", "main", "png", "panels")
read_panel <- function(file, dir = PANEL_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob     <- read_panel("MAIN_panel_A_heatmap.png")
pB_grob     <- read_panel("MAIN_panel_B_scatters.png")
pB_leg_grob <- read_panel("MAIN_panel_B_legend.png")

layout_cfg <- list(
  comp_w = 470, comp_h = 300,
  tag_sz      = 14,
  title_sz    = 12,
  subtitle_sz = 9,
  a_title_y   = 0.835,
  a_legend_y  = 0.14,
  grid_top    = 0.80,
  grid_bot    = 0.22,
  b_x         = 0.43,
  b_w         = 0.55,
  b_title_x   = 0.636,
  b_leg_x = 0.70, b_leg_y = 0.195, b_leg_w = 0.24, b_leg_h = 0.04,
  a_grob_x = 0.01, a_grob_y = 0, a_grob_w = 0.60, a_grob_h = 0.96,
  a_let_x = 15,  a_let_y = 248,
  a_ttl_x = 58,  a_ttl_y = 248,
  a_sub_x = 58,  a_sub_y = 243,
  b_let_x = 285, b_let_y = 248,
  b_ttl_x = 294, b_ttl_y = 248,
  b_sub_x = 294, b_sub_y = 243,
  crop_l = 0.01, crop_r = 0.80,
  crop_b = 0.17, crop_t = 0.85
)

COMP_W <- layout_cfg$comp_w
COMP_H <- layout_cfg$comp_h
TAG_SZ      <- layout_cfg$tag_sz
TITLE_SZ    <- layout_cfg$title_sz
SUBTITLE_SZ <- layout_cfg$subtitle_sz

A_TITLE_Y  <- layout_cfg$a_title_y
A_LEGEND_Y <- layout_cfg$a_legend_y
GRID_TOP   <- layout_cfg$grid_top
GRID_BOT   <- layout_cfg$grid_bot

B_X       <- layout_cfg$b_x
B_W       <- layout_cfg$b_w
B_TITLE_X <- layout_cfg$b_title_x
B_GRID_H  <- GRID_TOP - GRID_BOT

mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

A_LET_X <- layout_cfg$a_let_x; A_LET_Y <- layout_cfg$a_let_y
A_TTL_X <- layout_cfg$a_ttl_x; A_TTL_Y <- layout_cfg$a_ttl_y
A_SUB_X <- layout_cfg$a_sub_x; A_SUB_Y <- layout_cfg$a_sub_y

B_LET_X <- layout_cfg$b_let_x; B_LET_Y <- layout_cfg$b_let_y
B_TTL_X <- layout_cfg$b_ttl_x; B_TTL_Y <- layout_cfg$b_ttl_y
B_SUB_X <- layout_cfg$b_sub_x; B_SUB_Y <- layout_cfg$b_sub_y

CROP_L <- layout_cfg$crop_l; CROP_R <- layout_cfg$crop_r
CROP_B <- layout_cfg$crop_b; CROP_T <- layout_cfg$crop_t
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B first (behind) so its white left margin is hidden under Panel A
  draw_grob(pB_grob, x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg_grob, x = layout_cfg$b_leg_x, y = layout_cfg$b_leg_y,
            width = layout_cfg$b_leg_w, height = layout_cfg$b_leg_h,
            hjust = 0.5, vjust = 0.5) +
  draw_grob(pA_grob, x = layout_cfg$a_grob_x, y = layout_cfg$a_grob_y,
            width = layout_cfg$a_grob_w, height = layout_cfg$a_grob_h,
            hjust = 0, vjust = 0) +
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
       device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F06_composite.png"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       dpi = 300)
message("F06 composite saved: MAIN_F06_composite.{pdf,png}")

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- file.path(BASE, "c_data")
f06 <- function(p) file.path(DAT, p)

# RDS matrix files converted to data frames for Excel embedding
.shared       <- readRDS(f06("shared_objects.rds"))
.MEs          <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))),      "sample_id")
.me_pre       <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))),   "subject_key")
.me_post      <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))),  "subject_key")
.delta_me     <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")

# Consolidate redundant cohort/metric/mode/check splits into long sheets
# Folds 22 sheets into 5 long-format sheets, while preserving sheet names that
# downstream consumers (notably F07 panels and supp scripts) read directly.
.read_long <- function(path, cohort, metric) {
  df <- safe_read(path)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  long <- tidyr::pivot_longer(df, -module, names_to = "trait", values_to = "value")
  long$cohort <- cohort
  long$metric <- metric
  long[, c("module", "trait", "cohort", "metric", "value")]
}
.collect <- function(stage) {
  variants <- list(
    list(c="combined", m="cor",     p=sprintf("wgcna/wgcna_%s_trait_correlations.csv", stage)),
    list(c="young",    m="cor",     p=sprintf("wgcna/wgcna_%s_trait_correlations_young.csv", stage)),
    list(c="old",      m="cor",     p=sprintf("wgcna/wgcna_%s_trait_correlations_old.csv", stage)),
    list(c="combined", m="pval_bh", p=sprintf("wgcna/wgcna_%s_trait_pvalues_bh.csv", stage)),
    list(c="young",    m="pval_bh", p=sprintf("wgcna/wgcna_%s_trait_pvalues_bh_young.csv", stage)),
    list(c="old",      m="pval_bh", p=sprintf("wgcna/wgcna_%s_trait_pvalues_bh_old.csv", stage)),
    list(c="young",    m="pval_raw",p=sprintf("wgcna/wgcna_%s_trait_pvalues_raw_young.csv", stage)),
    list(c="old",      m="pval_raw",p=sprintf("wgcna/wgcna_%s_trait_pvalues_raw_old.csv", stage)))
  do.call(rbind, lapply(variants, function(v) .read_long(f06(v$p), v$c, v$m)))
}
.baseline_trait_assoc <- .collect("baseline")
.change_trait_assoc   <- .collect("change")

# Module-trait associations (main combined design): combine cor + pval_bh
.mt_cor <- safe_read(f06("wgcna/wgcna_module_trait_correlations.csv"))
.mt_bh  <- safe_read(f06("wgcna/wgcna_module_trait_pvalues_bh.csv"))
.module_trait_assoc <- rbind(
  cbind(tidyr::pivot_longer(.mt_cor, -module, names_to = "trait", values_to = "value"),
        metric = "cor"),
  cbind(tidyr::pivot_longer(.mt_bh,  -module, names_to = "trait", values_to = "value"),
        metric = "pval_bh"))
.module_trait_assoc <- .module_trait_assoc[, c("module", "trait", "metric", "value")]

# Module enrichment: strict + relaxed (+ mode)
.enr_strict  <- safe_read(f06("wgcna/wgcna_module_enrichment.csv"))
.enr_relaxed <- safe_read(f06("03_panel_B_triptych_enrichment.csv"))
if (!is.null(.enr_strict))  .enr_strict$mode  <- "strict"
if (!is.null(.enr_relaxed)) .enr_relaxed$mode <- "relaxed"
.module_enrichment <- dplyr::bind_rows(.enr_strict, .enr_relaxed)

# LMM diagnostics: contrast + stratified (+ check)
.lmm_c <- safe_read(f06("wgcna/wgcna_lmm_contrast_check.csv"))
.lmm_s <- safe_read(f06("wgcna/wgcna_lmm_stratified_check.csv"))
if (!is.null(.lmm_c)) { .lmm_c$check <- "contrast";   .lmm_c$age_group <- NA_character_ }
if (!is.null(.lmm_s)) { .lmm_s$check <- "stratified"; .lmm_s$contrast  <- NA_character_ }
.lmm_diagnostics <- dplyr::bind_rows(.lmm_c, .lmm_s)

f06_specs <- list(
  list(name="panel_A_heatmap",                  path=f06("01_panel_A_heatmap_data.csv")),
  list(name="panel_B_module_fgsea",             path=f06("panel_B_module_fgsea.csv")),
  list(name="WGCNA_module_assignments",         path=f06("wgcna/wgcna_module_assignments.csv")),
  list(name="WGCNA_mod_bio_labels",             path=f06("mod_bio_labels.csv")),
  list(name="WGCNA_hub_proteins",               path=f06("wgcna/wgcna_hub_proteins.csv")),
  list(name="WGCNA_hub_network_edges",          path=f06("04_panel_D_hub_network.csv")),
  list(name="WGCNA_protein_zscores_by_group",   path=f06("03_panel_B_heatmap_zscores.csv")),
  list(name="WGCNA_module_enrichment",          df=.module_enrichment),
  list(name="WGCNA_module_trait_assoc",         df=.module_trait_assoc),
  list(name="WGCNA_baseline_trait_assoc",       df=.baseline_trait_assoc),
  list(name="WGCNA_change_trait_assoc",         df=.change_trait_assoc),
  list(name="WGCNA_sft_summary",                path=f06("wgcna/wgcna_sft_summary.csv")),
  list(name="WGCNA_lmm_diagnostics",            df=.lmm_diagnostics),
  list(name="WGCNA_module_preservation",        path=f06("05_panel_E_preservation.csv")),
  list(name="WGCNA_gs_phenotype_choices",       path=f06("wgcna/gs_phenotype_choices.csv")),
  list(name="metadata_subj_age",                path=f06("subj_age.csv")),
  list(name="metadata_pheno_wide",              path=f06("pheno_wide.csv")),
  list(name="MEs",         df=.MEs),
  list(name="me_pre",      df=.me_pre),
  list(name="me_post",     df=.me_post),
  list(name="delta_me",    df=.delta_me),
  list(name="common_subj", df=data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE))
)

build_workbook(
  f06("F06_supplementary.xlsx"),
  sheet_specs = f06_specs
)
cleanup_after_workbook(f06_specs,
  extra_subdirs = character(0),
  preserve_patterns = c(
    "^00_input/", "^01_normalization/", "^02_imputation/", "^03_DEP/",
    "^04_Figures/shared/",
    "^04_Figures/F06/c_data/wgcna/",
    "^04_Figures/F06/c_data/.*\\.rds$",
    "^04_Figures/F06/c_data/key_modules\\.txt$",
    "^04_Figures/F06/c_data/mod_bio_labels\\.csv$",
    "^04_Figures/F06/c_data/meta\\.csv$",
    "^04_Figures/F06/c_data/imp_annotations\\.csv$",
    "^04_Figures/F06/c_data/subj_age\\.csv$",
    "^04_Figures/F06/c_data/pheno_wide\\.csv$"))

message("F06 main panels + xlsx complete")
