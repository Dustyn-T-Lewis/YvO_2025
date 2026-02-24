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
  library(RRHO2)
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
setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")

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
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

SIG_COLORS <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both"       = "#2E7D32",
  "Sig Young only" = "#E05A4E",
  "Sig Old only"   = "#5DA5DA",
  "NS"             = "grey55"
)

SIG_LABEL_FILL <- c(
  "Interaction"    = alpha("#7B5EA7", 0.75),
  "Sig Both"       = alpha("#2E7D32", 0.75),
  "Sig Young only" = alpha("#E05A4E", 0.75),
  "Sig Old only"   = alpha("#5DA5DA", 0.75),
  "NS"             = alpha("grey55",  0.75)
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

clean_pathway_name <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_remove("^REACTOME_") %>%
    str_remove("^KEGG_MEDICUS_") %>%
    str_remove("^KEGG_") %>%
    str_remove("^WP_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Mtorc1", "mTORC1") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f ", "E2F ") %>%
    str_replace("Dna ", "DNA ") %>%
    str_replace("Rna ", "RNA ") %>%
    str_replace("Tnfa ", "TNFa ") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("G2m ", "G2M ") %>%
    str_replace("Il6 ", "IL6 ") %>%
    str_replace("Il2 ", "IL2 ") %>%
    str_replace("Kras ", "KRAS ") %>%
    str_replace("P53 ", "p53 ") %>%
    str_replace("Tgf ", "TGF ") %>%
    str_replace("Nf Kb", "NF-kB") %>%
    str_replace("Atp ", "ATP ") %>%
    str_replace("Nadh ", "NADH ") %>%
    str_replace("Oxidative Phosphorylation", "OXPHOS") %>%
    str_replace("External Encapsulating Structure Or.*",
                "Extracellular Matrix Organization") %>%
    str_replace("Enzyme Linked Receptor Protein Signaling.*",
                "Receptor Protein Signaling") %>%
    str_trunc(45, ellipsis = "...")
}

sig_stars <- function(padj) {
  case_when(padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ "")
}

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

darken_color <- function(col, factor = 0.5) {
  rgb_vals <- col2rgb(col) / 255
  sapply(seq_along(col), function(i)
    rgb(rgb_vals[1, i] * factor, rgb_vals[2, i] * factor, rgb_vals[3, i] * factor))
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
