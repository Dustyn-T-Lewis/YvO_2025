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
# prepare_ring_data() — select GO terms and compute arc geometry
# ═══════════════════════════════════════════════════════════════════════════════
#' @param go_df       fGSEA results with columns: pathway, padj, NES,
#'                    leadingEdge (semicolon-separated), contrast, database
#' @param contrast    Character, which contrast to filter to
#' @param n_terms     Number of top terms to display
#' @param gap_degrees Angular gap (degrees) between adjacent arcs
#' @param start_offset Offset in degrees so the first arc starts at 12 o'clock
#'                    (90 = top; standard math convention is CCW from +x axis)
#' @param databases   Character vector of database names to include
#' @return tibble with one row per term, including arc geometry columns
prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms      = 12,
                              gap_degrees  = 3,
                              start_offset = 90,
                              databases    = c("Hallmark", "GO:BP", "GO:CC")) {

  # Filter to target contrast, databases, and significance

  ring <- go_df %>%
    filter(contrast == !!contrast,
           database %in% databases,
           padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = n_terms)

  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms found for contrast '",
            contrast, "' — returning empty tibble")
    return(tibble())
  }

  # Compute equal-width arc angles
  arc_width_deg <- (360 - n * gap_degrees) / n

  ring <- ring %>%
    mutate(
      term_idx    = row_number(),
      # Each arc: start at offset, then index * (width + gap)
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      # Convert to radians (ggforce uses radians)
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg   * pi / 180,
      mid_rad     = mid_deg   * pi / 180,
      # Clean labels for display
      clean_label = clean_pathway_name(pathway),
      # Parse leading edge genes into a list-column
      gene_list   = str_split(leadingEdge, ";")
    )

  ring
}

# ═══════════════════════════════════════════════════════════════════════════════
# build_tick_data() — map DE genes onto arcs as radial ticks
# ═══════════════════════════════════════════════════════════════════════════════
#' @param ring_data   Output of prepare_ring_data()
#' @param de_df       DE results with columns: gene, logFC_<contrast>
#' @param contrast    Contrast name (used to select logFC column)
#' @param tick_r0     Inner radius of tick ring
#' @param tick_r1     Outer radius of tick ring
#' @return tibble with one row per tick: gene, logFC, direction, angle_rad,
#'         x0, y0, x1, y1, term_idx, pathway
build_tick_data <- function(ring_data,
                            de_df,
                            contrast,
                            tick_r0 = 4.0,
                            tick_r1 = 5.2) {

  if (nrow(ring_data) == 0) return(tibble())

  logfc_col <- paste0("logFC_", contrast)

  # Build a simple lookup: gene -> logFC
  gene_lfc <- de_df %>%
    select(gene, lfc = all_of(logfc_col)) %>%
    filter(!is.na(lfc)) %>%
    distinct(gene, .keep_all = TRUE)

  # Angular padding (radians) inside each arc boundary
  pad_rad <- 0.5 * pi / 180   # 0.5 degrees

  tick_rows <- map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) return(tibble())

    # Evenly space ticks within the arc (with padding)
    arc_start <- row$start_rad + pad_rad
    arc_end   <- row$end_rad   - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad  # safety

    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)

    matched <- gene_lfc %>% filter(gene %in% genes_in_arc)
    # Align order to genes_in_arc
    matched <- matched[match(genes_in_arc, matched$gene), ] %>%
      filter(!is.na(gene))

    # Re-trim in case of NAs from match
    n_final <- nrow(matched)
    if (n_final == 0) return(tibble())
    tick_angles <- tick_angles[seq_len(n_final)]

    tibble(
      gene      = matched$gene,
      logFC     = matched$lfc,
      direction = ifelse(matched$lfc > 0, "Up", "Down"),
      angle_rad = tick_angles,
      x0        = tick_r0 * cos(tick_angles),
      y0        = tick_r0 * sin(tick_angles),
      x1        = tick_r1 * cos(tick_angles),
      y1        = tick_r1 * sin(tick_angles),
      term_idx  = row$term_idx,
      pathway   = row$pathway
    )
  })

  tick_rows
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
