################################################################################
#   YvO Figure 6 — Shared Setup
#   Loads all libraries, palettes, themes, helpers, data, and computed objects
#   that are needed by multiple panel scripts.
#
#   Source this file at the top of each panel_*.R script.
################################################################################
#
# STAT AUDIT — F6 Setup (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Data provenance: metadata derived from DAList .rds, not regex — PASS.
# 2. Subject pairing: uses subject_key = sub("_(Pre|Post)$", "", Col_ID),
#    with inner_join ensuring only paired subjects enter delta calculations.
# 3. Gene significance (gs_vl, gs_lbm): Pearson r of baseline expression
#    vs delta phenotype — no CI computed here; downstream panels should add.
# 4. kME computed on full datExpr (all timepoints). This is standard WGCNA
#    practice but the eigengene space mixes Pre+Post samples.
# 5. Added fisher_z_ci() helper: computes 95% CI for Pearson/partial r via
#    Fisher z-transform (used by Panels B, C, E).
# ---------------------------------------------------------------------------
#

# 0. LIBRARIES =================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)

  library(ggrepel)
  library(WGCNA)
  library(pROC)
  library(UpSetR)
  library(grid)
  library(gridExtra)
  library(png)
})

allowWGCNAThreads()
set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

# --- Output directories ---
RPT_DIR <- "04_Figures/F6/b_reports"
DAT_DIR <- "04_Figures/F6/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Canonical palettes (shared) ---
source("04_Figures/shared/palettes.R")
# KEY_TEXT, KEY_TITLE, KEY_ITEM etc. now centralized in palettes.R

mod_bio_labels <- c(
  brown       = "Brown",
  turquoise   = "Turquoise",
  yellow      = "Yellow",
  black       = "Black",
  blue        = "Blue",
  green       = "Green",
  red         = "Red",
  pink        = "Pink",
  magenta     = "Magenta",
  purple      = "Purple",
  cyan        = "Cyan",
  tan         = "Tan",
  salmon      = "Salmon",
  greenyellow = "Greenyellow"
)

# --- Outcome labels (used by panels B and C) ---
outcome_labels <- c(delta_VL  = "Delta VL (cm)",
                    delta_LBM = "Delta LBM (kg)")

# 1. LOAD DATA ================================================================

module_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_assignments.csv",
                       show_col_types = FALSE)
hub_df    <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_hub_proteins.csv",
                       show_col_types = FALSE)
trait_cor <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_trait_correlations.csv",
                       show_col_types = FALSE)

imp_df    <- read_csv("02_Imputation/c_data/01_imputed.csv",
                       show_col_types = FALSE)
ann_cols  <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
imp_mat   <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$uniprot_id
datExpr   <- t(imp_mat)                       # samples x proteins

# Use DAList metadata
dal      <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta <- as.data.frame(dal$metadata)

# Build metadata with a UNIQUE subject key derived from Col_ID
# (Subject_ID alone is ambiguous across age groups)
meta <- tibble(
  sample_id   = dal_meta$Col_ID,
  Subject_ID  = dal_meta$Subject_ID,
  age         = dal_meta$Group,
  time        = dal_meta$Timepoint,
  group       = dal_meta$Group_Time,
  subject_key = sub("_(Pre|Post)$", "", dal_meta$Col_ID)
)

# Attach phenotype columns
for (pc in c("VL_thick_cm", "DXA_LBM_kg")) {
  if (pc %in% names(dal_meta)) {
    meta[[pc]] <- as.numeric(as.character(dal_meta[[pc]]))
  }
}

# 2. COMPUTE SHARED OBJECTS ===================================================

# Module eigengenes (rows = samples, cols = ME*)
module_colors <- module_df$module_color
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes

# ---- Per-subject splits using subject_key ----
pre_meta  <- meta %>% filter(time == "Pre")
post_meta <- meta %>% filter(time == "Post")

# MEs use datExpr rownames (= sample_id)
me_pre_raw  <- MEs[pre_meta$sample_id, , drop = FALSE]
me_post_raw <- MEs[post_meta$sample_id, , drop = FALSE]
rownames(me_pre_raw)  <- pre_meta$subject_key
rownames(me_post_raw) <- post_meta$subject_key

common_subj <- intersect(rownames(me_pre_raw), rownames(me_post_raw))
me_pre  <- me_pre_raw[common_subj, , drop = FALSE]
me_post <- me_post_raw[common_subj, , drop = FALSE]
delta_me <- me_post - me_pre

# ---- Delta phenotype per subject_key ----
pheno_pre <- meta %>%
  filter(time == "Pre") %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg) %>%
  rename(VL_Pre = VL_thick_cm, LBM_Pre = DXA_LBM_kg)

pheno_post <- meta %>%
  filter(time == "Post") %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg) %>%
  rename(VL_Post = VL_thick_cm, LBM_Post = DXA_LBM_kg)

pheno_wide <- inner_join(pheno_pre, pheno_post, by = "subject_key") %>%
  mutate(delta_VL  = VL_Post  - VL_Pre,
         delta_LBM = LBM_Post - LBM_Pre) %>%
  filter(subject_key %in% common_subj)

# Subject-level age lookup (from Pre timepoint)
subj_age <- meta %>%
  filter(time == "Pre") %>%
  dplyr::select(subject_key, age) %>%
  distinct()

# ---- Gene significance: baseline protein vs delta phenotype ----
pre_expr <- datExpr[pre_meta$sample_id, , drop = FALSE]
rownames(pre_expr) <- pre_meta$subject_key
pre_expr <- pre_expr[common_subj, , drop = FALSE]

delta_vl_vec  <- pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)]
delta_lbm_vec <- pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)]

gs_vl <- cor(pre_expr, delta_vl_vec, use = "pairwise.complete.obs")
colnames(gs_vl) <- "GS_deltaVL"

gs_lbm <- cor(pre_expr, delta_lbm_vec, use = "pairwise.complete.obs")
colnames(gs_lbm) <- "GS_deltaLBM"

# ---- kME for all proteins (used by panels D, E) ----
kME_all <- signedKME(datExpr, MEs)

# ---- Baseline prediction correlations (used by panels B, D, E) ----
pred_cor <- tibble(module = colnames(me_pre)) %>%
  rowwise() %>%
  mutate(
    r_vl  = cor(me_pre[common_subj, module],
                pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_lbm = cor(me_pre[common_subj, module],
                pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    max_r = max(abs(r_vl), abs(r_lbm), na.rm = TRUE)
  ) %>% ungroup()

top3 <- pred_cor %>% arrange(desc(max_r)) %>% head(3) %>% pull(module)

# ---- Fisher z-transform CI helper (used by panels B, C, E) ----
# Returns 95% CI for Pearson or partial r using Fisher's z-transform.
# For partial r, effective n = n - k (where k = number of covariates).
# Reference: Bonett & Wright 2000, Psychometrika 65:23-28
fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z <- atanh(r)
  se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

cat(sprintf("Prediction data: %d subjects, %d modules, %d proteins\n",
            length(common_subj), ncol(MEs), ncol(datExpr)))

cat("F6 setup loaded\n")
