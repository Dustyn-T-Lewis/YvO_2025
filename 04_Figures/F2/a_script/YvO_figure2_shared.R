################################################################################
#   Figure 2 — Shared Setup
#   Common packages, style constants, helpers, and data loading used by all
#   per-panel scripts (YvO_panel_AB.R, YvO_panel_C.R, etc.)
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
for (pnl in c("panel_A", "panel_B", "panel_C", "panel_D", "panel_E", "panel_F", "shared"))
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

THEME_PUB <- theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(face = "bold", size = 9),
    plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 6.5),
    legend.key.size  = unit(3, "mm")
  )

# === 4. HELPERS ===============================================================
# clean_pathway_name(), darken_color(), sig_stars() are loaded from palettes.R
# (sourced implicitly via panel scripts or available from shared/palettes.R)

source("04_Figures/shared/palettes.R")

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

# === 5. DATA LOADING ==========================================================

message("Loading data...")
dep_df <- read_csv("03_DEP/c_data/combined_results.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)

# Imputation status (MAR/MNAR/Complete per protein)
imputation_df <- read_csv("02_Imputation/c_data/mar_mnar_classification.csv",
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
