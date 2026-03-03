################################################################################
#   Figure 2 — Shared Setup
#   Common packages, style constants, helpers, and data loading used by all
#   per-panel scripts (panel_A.R through panel_F.R).
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
  library(msigdbr)
  library(fgsea)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggnewscale)
  library(ggforce)
  library(openxlsx)
})

# === 2. SETUP =================================================================

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

source("05_Figures/shared/style.R")

RPT_DIR <- "05_Figures/F2/b_reports"
DAT_DIR <- "05_Figures/F2/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
for (pnl in c("panel_A", "panel_B", "panel_C", "panel_D",
              "panel_E", "panel_F", "shared"))
  dir.create(file.path(DAT_DIR, pnl), recursive = TRUE, showWarnings = FALSE)

# === 3. F2-SPECIFIC PALETTES =================================================

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

ORA_CONCORDANT_COL  <- "#FFB74D"
ORA_DISCORDANT_COL  <- "#64B5F6"

# === 4. HELPERS ===============================================================

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

imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                          show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")
message(sprintf("  %d proteins with imputation status (%d imputed)",
                nrow(imputation_df), sum(imputation_df$imputed)))

# fGSEA cache (lives in 04_Figures where it was originally computed)
fgsea_cache <- "04_Figures/F2/c_data/shared/fgsea_tstat_all_v2.csv"
if (!file.exists(fgsea_cache)) stop("fGSEA cache not found at ", fgsea_cache)
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# --- Ensure Interaction contrast is present in fGSEA results -----------------
if (!"Interaction" %in% unique(fgsea_all$contrast)) {
  message("  Computing fGSEA for Interaction contrast...")
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    dplyr::select(gs_name, gene_symbol)

  hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
  gobp_list     <- split(gobp$gene_symbol, gobp$gs_name)

  t_int <- dep_df %>%
    dplyr::select(gene, t_Interaction) %>%
    filter(!is.na(t_Interaction)) %>%
    deframe()

  res_h <- fgseaMultilevel(hallmark_list, t_int,
                           minSize = 15, maxSize = 200, eps = 0) %>%
    as_tibble() %>%
    mutate(contrast = "Interaction", database = "Hallmark",
           leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  res_bp <- fgseaMultilevel(gobp_list, t_int,
                            minSize = 15, maxSize = 200, eps = 0) %>%
    as_tibble() %>%
    mutate(contrast = "Interaction", database = "GO:BP",
           leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  fgsea_all <- bind_rows(fgsea_all, res_h, res_bp)
  message(sprintf("  Added %d Interaction fGSEA results", nrow(res_h) + nrow(res_bp)))
}

message(sprintf("Loaded %d proteins, %d fGSEA results (%s contrasts)",
                nrow(dep_df), nrow(fgsea_all),
                paste(unique(fgsea_all$contrast), collapse = ", ")))
