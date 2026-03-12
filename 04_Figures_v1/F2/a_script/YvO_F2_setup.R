################################################################################
#   Figure 2 — Shared Setup
#   Common packages, style constants, helpers, and data loading used by all
#   per-panel scripts (panel_AB.R, panel_C.R, etc.)
#
#   This script is idempotent: sourcing it multiple times has no side effects.
################################################################################

# === 1. PACKAGES ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(RedRibbon)
  library(msigdbr)
  library(fgsea)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggnewscale)
  library(ComplexHeatmap)
  library(circlize)
  # ggalluvial removed — manual Sankey construction below
})

# === 2. SETUP =================================================================

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

RPT_DIR <- "04_Figures/F2/b_reports"
DAT_DIR <- "04_Figures/F2/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
# Per-panel subdirectories for organized data output
for (pnl in c("panel_A", "panel_B", "panel_C", "panel_D", "panel_E", "panel_F", "panel_G", "shared"))
  dir.create(file.path(DAT_DIR, pnl), recursive = TRUE, showWarnings = FALSE)

# === 3. CANONICAL STYLE (from figure-style guide) ============================

CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
# KEY_TEXT, KEY_TITLE, KEY_ITEM etc. now centralized in palettes.R

SIG_COLORS <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both"       = "#2E7D32",
  "Sig Young only" = "#E05A4E",
  "Sig Old only"   = "#5DA5DA",
  "NS"             = "grey70"
)

SIG_LABEL_FILL <- c(
  "Interaction"    = alpha("#7B5EA7", 0.75),
  "Sig Both"       = alpha("#2E7D32", 0.75),
  "Sig Young only" = alpha("#E05A4E", 0.75),
  "Sig Old only"   = alpha("#5DA5DA", 0.75),
  "NS"             = alpha("grey70",  0.75)
)
SIG_LABEL_TEXT <- c(
  "Interaction"    = "white",
  "Sig Both"       = "white",
  "Sig Young only" = "white",
  "Sig Old only"   = "white",
  "NS"             = "white"
)

# THEME_PUB, THEME_FIG, TXT_* constants loaded from shared/palettes.R

# === 4. HELPERS ===============================================================
# clean_pathway_name(), darken_color(), sig_stars() are loaded from palettes.R
# THEME_PUB, THEME_FIG, TXT_* are loaded from palettes.R
# assign_go_slim_super(), SLIM_SUPER, bp_slim, etc. from go_slim_categories.R

source("04_Figures/shared/palettes.R")
source("04_Figures/shared/go_slim_categories.R")

classify_proteins <- function(pi_A, pi_B, pi_int, threshold = 0.05) {
  case_when(
    pi_int < threshold              ~ "Interaction",
    pi_A < threshold & pi_B < threshold ~ "Sig Both",
    pi_A < threshold                ~ "Sig Young only",
    pi_B < threshold                ~ "Sig Old only",
    TRUE                            ~ "NS"
  ) %>%
    factor(levels = c("Interaction", "Sig Both",
                       "Sig Young only", "Sig Old only", "NS"))
}

# --- Expanded 7-category classification (direction-split) -------------------
# Splits existing sig groups by logFC direction for richer heatmap display.
# Interaction proteins keep their single category (direction is implicit).

EXPANDED_ORDER <- c("Interaction",
                    "Sig Both Up", "Sig Both Down",
                    "Sig Young Up", "Sig Young Down",
                    "Sig Old Up", "Sig Old Down")

EXPANDED_COLORS <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both Up"    = "#2E7D32",
  "Sig Both Down"  = "#1B5E20",
  "Sig Young Up"   = "#E05A4E",
  "Sig Young Down" = "#B71C1C",
  "Sig Old Up"     = "#5DA5DA",
  "Sig Old Down"   = "#1565C0"
)

# ORA flank colors for Panel D RRHO quadrant enrichment
ORA_CONCORDANT_COL  <- "#FFB74D"
ORA_DISCORDANT_COL  <- "#64B5F6"

# ORA bar colors — one per RRHO quadrant
ORA_QUAD_COLORS <- c(
  "Concordant Up"               = "#E57373",   # warm red (both up)
  "Concordant Down"             = "#64B5F6",   # cool blue (both down)
  "Discordant (Y Up / O Down)"  = "#FFB74D",   # orange
  "Discordant (Y Down / O Up)"  = "#81C784"    # green
)

classify_expanded <- function(pi_Y, pi_O, pi_int, logFC_Y, logFC_O,
                               threshold = 0.05) {
  base <- classify_proteins(pi_Y, pi_O, pi_int, threshold)
  case_when(
    base == "Interaction"    ~ "Interaction",
    base == "Sig Both"   & logFC_Y > 0 ~ "Sig Both Up",
    base == "Sig Both"   & logFC_Y <= 0 ~ "Sig Both Down",
    base == "Sig Young only" & logFC_Y > 0 ~ "Sig Young Up",
    base == "Sig Young only" & logFC_Y <= 0 ~ "Sig Young Down",
    base == "Sig Old only"   & logFC_O > 0 ~ "Sig Old Up",
    base == "Sig Old only"   & logFC_O <= 0 ~ "Sig Old Down",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = EXPANDED_ORDER)
}

# === 5. DATA LOADING ==========================================================

message("Loading data...")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)

# Imputation status (MAR/MNAR/Complete per protein)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                          show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")
message(sprintf("  %d proteins with imputation status (%d imputed)",
                nrow(imputation_df), sum(imputation_df$imputed)))

# NOTE: fGSEA is run independently per contrast without cross-contrast correction.
# This is standard practice (each contrast = independent question), but increases
# the chance of spurious enrichment in at least one contrast. Acknowledge in methods.
fgsea_cache <- file.path(DAT_DIR, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) stop("fGSEA cache not found at ", fgsea_cache)
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))
