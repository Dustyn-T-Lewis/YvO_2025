#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# volcano_ring.R
# Self-contained utility for creating a circular volcano-in-ring composite plot.
# Architecture: standard Cartesian ggplot with ggforce::geom_arc_bar() arcs;
# NO coord_polar(). Everything lives in one coord_fixed() space.
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Packages ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
  library(scales)
})

# ─── Style constants ─────────────────────────────────────────────────────────
CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
NS_COLOR   <- "grey70"
KEY_TEXT    <- 2.2
KEY_TITLE   <- 2.3

# ─── clean_pathway_name (local copy, self-contained) ─────────────────────────
clean_pathway_name <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
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
    str_replace("Ifn", "IFN") %>%
    str_replace("Pi3k", "PI3K") %>%
    str_replace("Akt", "AKT") %>%
    str_replace("Mtor", "mTOR") %>%
    str_replace("Oxidative Phosphorylation", "OXPHOS") %>%
    str_trunc(45, ellipsis = "...")
}

# ═══════════════════════════════════════════════════════════════════════════════
# make_volcano_ring() — main entry point (stub)
# ═══════════════════════════════════════════════════════════════════════════════
make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title           = NULL,
                              n_terms         = 12,
                              gap_degrees     = 3,
                              start_offset    = 90,
                              databases       = c("Hallmark", "GO:BP", "GO:CC"),
                              volcano_radius  = 3.5,
                              tick_r0         = 4.0,
                              tick_r1         = 5.2,
                              arc_r0          = 5.4,
                              arc_r1          = 6.0,
                              label_r         = 6.5,
                              fc_thresh       = log2(1.5),
                              p_thresh        = 0.05,
                              up_color        = DIR_COLORS["Up"],
                              down_color      = DIR_COLORS["Down"],
                              ns_color        = NS_COLOR,
                              point_size      = 0.6,
                              point_alpha     = 0.5,
                              label_size      = 2.0,
                              ring_data_override = NULL) {

  # Stub: returns empty canvas

  ggplot() + theme_void()
}
