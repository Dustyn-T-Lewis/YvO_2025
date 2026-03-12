################################################################################
#   YvO Figure 5 — Shared Setup
#   Sources: libraries, palettes, themes, helpers, and all shared data objects
#   Sourced by each panel_X.R script
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# DATA SCOPE
#   62 samples (31 subjects × 2 timepoints), ~2100 proteins, 14 non-grey
#   WGCNA modules.  Langfelder recommends n > 15 samples; we have 62.
#
# MODULE-TRAIT CORRELATIONS (computed in YvO_WGCNA_run.R, consumed here)
#   Pearson r via WGCNA::cor(); p-values via corPvalueStudent().
#   BH correction is applied GLOBALLY across all module × trait tests
#   (15 modules × 9 traits = 135 raw p-values flattened then adjusted).
#   See YvO_WGCNA_run.R lines 228–237.
#
#   IMPORTANT: corPvalueStudent() assumes n = nrow(datExpr) = 62 for all
#   traits.  However, phenotype traits (VL_thick_cm, Type_I_fCSA, etc.)
#   have missing data (n = 44–62).  cor(..., use = "pairwise.complete.obs")
#   adapts the correlation itself, but p-values use the full n, making them
#   LIBERAL for incomplete traits.  Panel B documents this caveat.
#
# SHARED OBJECTS exported to panel scripts:
#   meta, module_df, MEs, kME_all, imp_mat, datExpr, trait_cor, pval_bh, go_df
# ─────────────────────────────────────────────────────────────────────────────

# 0. LIBRARIES =================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(WGCNA)
  library(igraph)
  library(ggraph)
  library(cowplot)
  library(grid)
  library(png)
})

allowWGCNAThreads()
set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

RPT_DIR <- "04_Figures/F5/b_reports"
DAT_DIR <- "04_Figures/F5/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. CANONICAL PALETTES (shared across all figures) ============================
source("04_Figures/shared/palettes.R")
# KEY_TEXT, KEY_TITLE, KEY_ITEM etc. now centralized in palettes.R

# 2. HELPER FUNCTIONS ==========================================================
# clean_pathway_name(), darken_color(), sig_stars() are loaded from palettes.R above

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

# 3. LOAD DATA =================================================================

# MANUSCRIPT NOTE: 62 samples for WGCNA on ~2,100 proteins is within acceptable range
# (Langfelder recommends >15). O'Leary et al. 2025 used 15 participants successfully.

net <- readRDS("04_Figures/F5/c_data/wgcna/wgcna_network.rds")
module_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_assignments.csv", show_col_types = FALSE)
hub_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_hub_proteins.csv", show_col_types = FALSE)
trait_cor <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_trait_correlations.csv", show_col_types = FALSE)
pval_bh <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_trait_pvalues_bh.csv", show_col_types = FALSE)
go_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_GO_enrichment.csv", show_col_types = FALSE)

# Imputed data + metadata
imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$uniprot_id

# WGCNA convention: samples as rows, proteins as columns
datExpr <- t(imp_mat)

# Use DAList metadata for accurate group assignment
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta <- as.data.frame(dal$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID),
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time
)
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

# Phenotype columns from metadata
pheno_cols <- c("VL_thick_cm", "DXA_LBM_kg", "BMI",
                "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg")
for (pc in pheno_cols) {
  if (pc %in% names(dal_meta)) {
    vals <- dal_meta[[pc]]
    if (!is.numeric(vals)) vals <- as.numeric(as.character(vals))
    meta[[pc]] <- vals[match(meta$sample_id, dal_meta$Col_ID)]
  }
}

module_colors <- module_df$module_color
n_modules <- length(unique(module_colors[module_colors != "grey"]))
cat(sprintf("Loaded: %d modules, %d proteins, %d samples\n",
            n_modules, ncol(datExpr), nrow(datExpr)))

# 4. SHARED COMPUTATIONS =======================================================

# Module eigengenes — used by panels B, C, E, F
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes

# kME values for all proteins — used by panels D, E, and F6
kME_all <- signedKME(datExpr, MEs)

# 5. BIOLOGICAL MODULE LABELS =================================================
# Maps WGCNA color names to biology-derived labels for publication figures.
# Based on GO enrichment and top hub protein functions.
mod_bio_labels <- c(
  blue        = "Catabolic Metabolism",
  brown       = "Cytoskeletal Remodeling",
  turquoise   = "Protein Modification",
  green       = "Muscle Structure",
  black       = "Translation Machinery",
  pink        = "Immune/Complement",
  yellow      = "Ubiquitin-Proteasome",
  red         = "Red",
  magenta     = "Magenta"
)

# Helper: produce "Biology (color)" labels for axis display
mod_display_label <- function(color) {

  lbl <- mod_bio_labels[color]
  lbl[is.na(lbl)] <- str_to_title(color[is.na(lbl)])
  paste0(lbl, " (", str_to_title(color), ")")
}

cat("F5 setup loaded\n")
