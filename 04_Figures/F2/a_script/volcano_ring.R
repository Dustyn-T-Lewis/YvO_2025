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
# build_volcano_layers() — center volcano plot as ggplot layer list
# ═══════════════════════════════════════════════════════════════════════════════
#' @param de_df         DE results data frame
#' @param contrast      Contrast name
#' @param volcano_radius Radius of the circle the volcano is inscribed in
#' @param fc_thresh     logFC threshold for significance lines
#' @param p_thresh      P-value threshold for significance line
#' @param up_color      Color for upregulated points
#' @param down_color    Color for downregulated points
#' @param ns_color      Color for non-significant points
#' @param point_size    Point size
#' @param point_alpha   Point alpha
#' @return list of ggplot layers; attributes: x_max, y_max, n_up, n_down,
#'         scale_x (function mapping raw logFC -> plot x),
#'         scale_y (function mapping raw -log10p -> plot y)
build_volcano_layers <- function(de_df,
                                 contrast,
                                 volcano_radius = 3.5,
                                 fc_thresh      = log2(1.5),
                                 p_thresh       = 0.05,
                                 up_color       = DIR_COLORS["Up"],
                                 down_color     = DIR_COLORS["Down"],
                                 ns_color       = NS_COLOR,
                                 point_size     = 0.6,
                                 point_alpha    = 0.5) {

  logfc_col <- paste0("logFC_", contrast)
  pval_col  <- paste0("P.Value_", contrast)
  pi_col    <- paste0("pi_score_", contrast)

  # Prepare volcano data
  vdf <- de_df %>%
    transmute(
      gene     = gene,
      logFC    = .data[[logfc_col]],
      pvalue   = .data[[pval_col]],
      pi_score = .data[[pi_col]],
      neg_log10p = -log10(pvalue)
    ) %>%
    filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p))

  # Classify significance using pi-score

  vdf <- vdf %>%
    mutate(
      direction = case_when(
        pi_score < 0.05 & logFC > 0 ~ "Up",
        pi_score < 0.05 & logFC < 0 ~ "Down",
        TRUE                         ~ "NS"
      )
    )

  n_up   <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")

  # Compute data ranges (symmetric x, 0-based y)
  x_data_max <- max(abs(vdf$logFC), na.rm = TRUE)
  y_data_max <- max(vdf$neg_log10p, na.rm = TRUE)

  # Use a small margin so points don't sit on the boundary
  margin <- 0.92
  vr <- volcano_radius * margin

  # Linear scaling functions: data space -> plot Cartesian space
  scale_x_fn <- function(val) val / x_data_max * vr
  scale_y_fn <- function(val) val / y_data_max * vr

  # Apply scaling
  vdf <- vdf %>%
    mutate(
      x_plot = scale_x_fn(logFC),
      y_plot = scale_y_fn(neg_log10p) - vr  # shift so y=0 maps to -vr (bottom)
    )
  # Now y ranges from -vr (p=1) to 0 (max -log10p)... actually we want
  # the volcano centered. Let's center it: y goes from -vr to +vr
  # with 0 at the middle of the y range.
  # Better approach: map [0, y_data_max] -> [-vr, +vr]
  vdf <- vdf %>%
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  # Separate NS and significant for layering
  vdf_ns  <- vdf %>% filter(direction == "NS")
  vdf_sig <- vdf %>% filter(direction != "NS")

  layers <- list()

  # --- Circle boundary (optional subtle guide) ---
  circle_df <- tibble(
    angle = seq(0, 2 * pi, length.out = 200),
    cx    = volcano_radius * cos(angle),
    cy    = volcano_radius * sin(angle)
  )
  layers$circle <- geom_path(
    data = circle_df, aes(x = cx, y = cy),
    color = "grey85", linewidth = 0.3, linetype = "dotted",
    inherit.aes = FALSE
  )

  # --- NS points (background, smaller, more transparent) ---
  layers$ns_points <- geom_point(
    data = vdf_ns, aes(x = x_plot, y = y_plot),
    color = ns_color, size = point_size * 0.7, alpha = point_alpha * 0.5,
    inherit.aes = FALSE
  )

  # --- Significant points (on top) ---
  layers$sig_points <- geom_point(
    data = vdf_sig, aes(x = x_plot, y = y_plot, color = direction),
    size = point_size, alpha = point_alpha,
    inherit.aes = FALSE
  )

  # --- Color scale for direction ---
  layers$color_scale <- scale_color_manual(
    values = c(Up = unname(up_color), Down = unname(down_color)),
    guide  = "none"
  )

  # --- Axis lines ---
  # x-axis (horizontal, at y corresponding to -log10(1)=0 -> y_plot = -vr)
  # Actually, let's put the x-axis at the y=0 line in data space
  # Data y=0 -> y_plot = -vr. That's the bottom. Better: draw at the data center.
  # For a classic volcano, x-axis is at bottom (p=1) and y increases upward.
  # In our scaled space, y_plot for neg_log10p=0 is -vr.
  y_axis_pos <- -vr  # where -log10(p)=0

  layers$x_axis <- geom_segment(
    data = tibble(x = -vr, xend = vr, y = y_axis_pos, yend = y_axis_pos),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey40", linewidth = 0.3,
    inherit.aes = FALSE
  )

  # y-axis (vertical, at logFC=0 -> x_plot = 0)
  layers$y_axis <- geom_segment(
    data = tibble(x = 0, xend = 0, y = -vr, yend = vr),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey40", linewidth = 0.3,
    inherit.aes = FALSE
  )

  # --- Dashed threshold lines ---
  # Vertical lines at +/- fc_thresh
  fc_x_pos <- scale_x_fn(fc_thresh)
  fc_x_neg <- scale_x_fn(-fc_thresh)

  layers$fc_thresh_pos <- geom_segment(
    data = tibble(x = fc_x_pos, xend = fc_x_pos, y = -vr, yend = vr),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50", linewidth = 0.2, linetype = "dashed",
    inherit.aes = FALSE
  )
  layers$fc_thresh_neg <- geom_segment(
    data = tibble(x = fc_x_neg, xend = fc_x_neg, y = -vr, yend = vr),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50", linewidth = 0.2, linetype = "dashed",
    inherit.aes = FALSE
  )

  # Horizontal line at P = p_thresh
  p_y <- (-log10(p_thresh) / y_data_max) * 2 * vr - vr
  layers$p_thresh <- geom_segment(
    data = tibble(x = -vr, xend = vr, y = p_y, yend = p_y),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50", linewidth = 0.2, linetype = "dashed",
    inherit.aes = FALSE
  )

  # --- Axis tick marks and labels ---
  # X-axis ticks: nice intervals for logFC
  x_breaks_data <- pretty(c(-x_data_max, x_data_max), n = 5)
  x_breaks_data <- x_breaks_data[abs(x_breaks_data) <= x_data_max * 1.05]
  x_breaks_plot <- scale_x_fn(x_breaks_data)
  tick_len <- vr * 0.03

  x_tick_df <- tibble(
    x = x_breaks_plot,
    y    = y_axis_pos,
    yend = y_axis_pos - tick_len,
    label = as.character(round(x_breaks_data, 1))
  )

  layers$x_ticks <- geom_segment(
    data = x_tick_df, aes(x = x, y = y, xend = x, yend = yend),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )
  layers$x_tick_labels <- geom_text(
    data = x_tick_df, aes(x = x, y = yend - tick_len * 0.8, label = label),
    size = 1.6, color = "grey30",
    inherit.aes = FALSE
  )

  # Y-axis ticks: nice intervals for -log10(p)
  y_breaks_data <- pretty(c(0, y_data_max), n = 5)
  y_breaks_data <- y_breaks_data[y_breaks_data >= 0 & y_breaks_data <= y_data_max * 1.05]
  y_breaks_plot <- (y_breaks_data / y_data_max) * 2 * vr - vr

  y_tick_df <- tibble(
    y = y_breaks_plot,
    x    = 0,
    xend = -tick_len,
    label = as.character(round(y_breaks_data, 0))
  )

  layers$y_ticks <- geom_segment(
    data = y_tick_df, aes(x = x, y = y, xend = xend, yend = y),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )
  layers$y_tick_labels <- geom_text(
    data = y_tick_df, aes(x = xend - tick_len * 0.8, y = y, label = label),
    size = 1.6, color = "grey30",
    inherit.aes = FALSE
  )

  # --- Axis title labels ---
  layers$x_title <- annotate(
    "text", x = 0, y = -vr - tick_len * 5,
    label = expression(log[2]~"Fold Change"),
    size = 2.0, color = "grey20"
  )
  layers$y_title <- annotate(
    "text", x = -vr * 0.45, y = vr * 0.15,
    label = expression(-log[10]~"(P)"),
    size = 2.0, color = "grey20", angle = 90
  )

  # --- P threshold label ---
  layers$p_label <- annotate(
    "text", x = vr * 0.85, y = p_y + tick_len * 1.5,
    label = paste0("P = ", p_thresh),
    size = 1.5, color = "grey40", fontface = "italic"
  )

  # --- DEP count annotations ---
  layers$n_up_label <- annotate(
    "text", x = vr * 0.7, y = vr * 0.85,
    label = paste0(n_up, " Up"),
    size = 1.8, color = up_color, fontface = "bold"
  )
  layers$n_down_label <- annotate(
    "text", x = -vr * 0.7, y = vr * 0.85,
    label = paste0(n_down, " Down"),
    size = 1.8, color = down_color, fontface = "bold"
  )

  # Attach metadata as attributes
  attr(layers, "x_data_max") <- x_data_max
  attr(layers, "y_data_max") <- y_data_max
  attr(layers, "n_up")       <- n_up
  attr(layers, "n_down")     <- n_down
  attr(layers, "scale_x")    <- scale_x_fn
  attr(layers, "scale_y")    <- scale_y_fn

  layers
}

# ═══════════════════════════════════════════════════════════════════════════════
# build_ring_layers() — tick ring + enrichment arc ring
# ═══════════════════════════════════════════════════════════════════════════════
#' @param ring_data Output of prepare_ring_data()
#' @param tick_data Output of build_tick_data()
#' @param tick_r0   Inner radius of tick ring
#' @param tick_r1   Outer radius of tick ring
#' @param arc_r0    Inner radius of enrichment arc ring
#' @param arc_r1    Outer radius of enrichment arc ring
#' @param up_color  Color for upregulated ticks
#' @param down_color Color for downregulated ticks
#' @return list of ggplot layers
build_ring_layers <- function(ring_data,
                              tick_data,
                              tick_r0    = 4.0,
                              tick_r1    = 5.2,
                              arc_r0     = 5.4,
                              arc_r1     = 6.0,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"]) {

  if (nrow(ring_data) == 0) return(list())

  layers <- list()

  # ── 1. Gray background arcs for tick ring ──────────────────────────────────
  # geom_arc_bar needs: x0, y0, r0 (inner), r (outer), start, end (radians)
  layers$tick_bg <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
        start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15,
    inherit.aes = FALSE
  )

  # ── 2. Guide lines within tick ring ────────────────────────────────────────
  # 3 evenly spaced radii between tick_r0 and tick_r1
  guide_radii <- seq(tick_r0, tick_r1, length.out = 5)[2:4]  # inner 3

  for (gr in guide_radii) {
    layers[[paste0("guide_r_", round(gr, 2))]] <- geom_arc(
      data = ring_data,
      aes(x0 = 0, y0 = 0, r = gr, start = start_rad, end = end_rad),
      color = "grey82", linewidth = 0.1,
      inherit.aes = FALSE
    )
  }

  # ── 3. Protein ticks ──────────────────────────────────────────────────────
  if (nrow(tick_data) > 0) {
    layers$ticks <- geom_segment(
      data = tick_data,
      aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  # ── 4. Enrichment arcs (outer ring, fill = NES) ───────────────────────────
  layers$enrich_arcs <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1,
        start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )

  # ── 5. Fill scale for NES ─────────────────────────────────────────────────
  nes_max <- max(abs(ring_data$NES), na.rm = TRUE) * 1.05
  layers$fill_scale <- scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0,
    limits = c(-nes_max, nes_max),
    name = "NES"
  )

  layers
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
