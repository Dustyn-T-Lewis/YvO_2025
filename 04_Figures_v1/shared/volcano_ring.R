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
  library(patchwork)
  library(scales)
})

# ─── Style constants ─────────────────────────────────────────────────────────
CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
NS_COLOR   <- "grey70"
# KEY_TEXT, KEY_TITLE etc. now centralized in palettes.R
if (!exists("KEY_TEXT")) source("04_Figures/shared/palettes.R")

# ─── clean_ring_label (ring-specific with smart abbreviations) ────────────────
clean_ring_label <- function(name) {
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
    # --- Smart abbreviations for ring labels ---
    str_replace("Mitochondrial", "Mito.") %>%
    str_replace("Ubiquinone", "UQ") %>%
    str_replace("Organization", "Org.") %>%
    str_replace("Cytoskeleton", "Cytoskel.") %>%
    str_replace("Microtubule", "MT") %>%
    str_replace("Respiratory", "Resp.") %>%
    str_replace("Electron Transport", "ETC") %>%
    str_replace("Synthesis Coupled", "Synth.-Coupled") %>%
    str_replace("Regulation Of", "Reg.") %>%
    str_replace("Signaling Pathway", "Signaling") %>%
    str_replace("Biosynthetic Process", "Biosynthesis") %>%
    str_replace("Catabolic Process", "Catabolism") %>%
    str_replace("Metabolic Process", "Metabolism") %>%
    str_replace("Based Process", "Process") %>%
    str_replace("Response To", "Resp. to") %>%
    str_replace("Extracellular Matrix", "ECM") %>%
    str_replace("Epithelial Mesenchymal Transition", "EMT") %>%
    str_trunc(40, ellipsis = "...")
}

# ═══════════════════════════════════════════════════════════════════════════════
# prepare_ring_data() — select GO terms and compute arc geometry
# ═══════════════════════════════════════════════════════════════════════════════
#' @param go_df       fGSEA results with columns: pathway, padj, NES,
#'                    leadingEdge (semicolon-separated), contrast, database
#' @param contrast    Character, which contrast to filter to
#' @param n_terms     Number of top terms to display
#' @param gap_degrees Angular gap (degrees) between adjacent arcs
#' @param start_offset Offset in degrees; 0 = 12 o'clock in ggforce convention
#'                    (CW from top; ggforce uses x=sin(θ), y=cos(θ))
#' @param databases   Character vector of database names to include
#' @return tibble with one row per term, including arc geometry columns
prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms      = 12,
                              gap_degrees  = 3,
                              start_offset = 0,
                              databases    = c("Hallmark", "GO:BP")) {

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
      clean_label = clean_ring_label(pathway),
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
    dplyr::select(gene, lfc = all_of(logfc_col)) %>%
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
      x0        = tick_r0 * sin(tick_angles),
      y0        = tick_r0 * cos(tick_angles),
      x1        = tick_r1 * sin(tick_angles),
      y1        = tick_r1 * cos(tick_angles),
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
  # STAT AUDIT: pi-score = -log10(p)*|logFC|; threshold 0.05 is a heuristic
  # combining statistical significance with biological effect size. See
  # panel_AB.R audit block for full justification.

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

  # Apply scaling: x maps [-x_data_max, +x_data_max] -> [-vr, +vr]
  #                 y maps [0, y_data_max]             -> [-vr, +vr]
  vdf <- vdf %>%
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  # Separate NS and significant for layering
  vdf_ns  <- vdf %>% filter(direction == "NS")
  vdf_sig <- vdf %>% filter(direction != "NS")

  layers <- list()

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

  # --- DEP count badges (colored fill, white bold text) ---
  layers$n_up_label <- annotate(
    "label",
    x = vr * 0.55, y = vr * 0.55,
    label = paste0(n_up, " Up"),
    size = 2.8,
    color = "white", fill = up_color,
    fontface = "bold",
    label.padding = unit(2, "pt"),
    label.r = unit(1.5, "pt"),
    linewidth = 0
  )
  layers$n_down_label <- annotate(
    "label",
    x = -vr * 0.55, y = vr * 0.55,
    label = paste0(n_down, " Down"),
    size = 2.8,
    color = "white", fill = down_color,
    fontface = "bold",
    label.padding = unit(2, "pt"),
    label.r = unit(1.5, "pt"),
    linewidth = 0
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
                              arc_r0     = 5.2,
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

  # ── 5. Fill scale for NES (clamped so both hemispheres are saturated) ─────
  nes_clamp <- 2.0
  layers$fill_scale <- scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0,
    limits = c(-nes_clamp, nes_clamp),
    oob = scales::squish,
    name = "NES"
  )

  layers
}

# ═══════════════════════════════════════════════════════════════════════════════
# build_label_layer() — radial term labels outside the outer ring
# ═══════════════════════════════════════════════════════════════════════════════
#' @param ring_data  Output of prepare_ring_data()
#' @param label_r    Radius at which to place label text
#' @param label_size Text size for labels
#' @return list with geom_text layers
build_label_layer <- function(ring_data,
                              label_r    = 7.0,
                              label_size = 2.8) {

  if (nrow(ring_data) == 0) return(list())

  # Compute label positions — text is always HORIZONTAL (angle = 0).
  # Labels are centered on the arc centroid (not angular midpoint) so they

  # sit at the visual balance-point of each curved bar.
  label_df <- ring_data %>%
    mutate(
      # Arc centroid direction (accounts for curvature)
      arc_span   = end_rad - start_rad,
      cx         = (cos(start_rad) - cos(end_rad)) / arc_span,
      cy         = (sin(end_rad) - sin(start_rad)) / arc_span,
      # Project to label radius
      centroid_r = sqrt(cx^2 + cy^2),
      x_label    = label_r * cx / centroid_r,
      y_label    = label_r * cy / centroid_r
    )

  layers <- list()

  layers$labels <- geom_label(
    data = label_df,
    aes(x = x_label, y = y_label, label = clean_label),
    hjust = 0.5, vjust = 0.5,
    size = label_size, color = "grey20", angle = 0,
    fill = "grey96", alpha = 0.85,
    linewidth = 0.15, label.padding = unit(1.5, "pt"),
    label.r = unit(0.1, "lines"),
    inherit.aes = FALSE
  )

  layers
}

# ═══════════════════════════════════════════════════════════════════════════════
# make_volcano_ring() — main entry point
# ═══════════════════════════════════════════════════════════════════════════════
make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title           = NULL,
                              n_terms         = 12,
                              gap_degrees     = 3,
                              start_offset    = 0,
                              databases       = c("Hallmark", "GO:BP"),
                              volcano_radius  = 3.5,
                              tick_r0         = 4.0,
                              tick_r1         = 5.2,
                              arc_r0          = 5.2,
                              arc_r1          = 6.0,
                              label_r         = 6.5,
                              fc_thresh       = log2(1.5),
                              p_thresh        = 0.05,
                              up_color        = DIR_COLORS["Up"],
                              down_color      = DIR_COLORS["Down"],
                              ns_color        = NS_COLOR,
                              point_size      = 0.6,
                              point_alpha     = 0.5,
                              label_size      = 2.8,
                              contrast_title     = NULL,
                              contrast_subtitle  = NULL,
                              title_size         = 12,
                              ring_data_override = NULL) {

  # ── 1. Prepare ring data ──────────────────────────────────────────────────
  if (!is.null(ring_data_override)) {
    ring_data <- ring_data_override
    message("make_volcano_ring: using ring_data_override (", nrow(ring_data), " terms)")
  } else {
    ring_data <- prepare_ring_data(
      go_df        = go_df,
      contrast     = contrast,
      n_terms      = n_terms,
      gap_degrees  = gap_degrees,
      start_offset = start_offset,
      databases    = databases
    )
  }

  # ── 2. Build tick data ────────────────────────────────────────────────────
  tick_data <- build_tick_data(
    ring_data = ring_data,
    de_df     = de_df,
    contrast  = contrast,
    tick_r0   = tick_r0,
    tick_r1   = tick_r1
  )

  # ── 3. Build layer groups ─────────────────────────────────────────────────
  volcano_layers <- build_volcano_layers(
    de_df          = de_df,
    contrast       = contrast,
    volcano_radius = volcano_radius,
    fc_thresh      = fc_thresh,
    p_thresh       = p_thresh,
    up_color       = up_color,
    down_color     = down_color,
    ns_color       = ns_color,
    point_size     = point_size,
    point_alpha    = point_alpha
  )

  ring_layers <- build_ring_layers(
    ring_data  = ring_data,
    tick_data  = tick_data,
    tick_r0    = tick_r0,
    tick_r1    = tick_r1,
    arc_r0     = arc_r0,
    arc_r1     = arc_r1,
    up_color   = up_color,
    down_color = down_color
  )

  label_layers <- build_label_layer(
    ring_data  = ring_data,
    label_r    = label_r,
    label_size = label_size
  )

  # ── 4. Compose into single ggplot ─────────────────────────────────────────
  # Layer order: ring bg -> volcano -> ring foreground -> labels
  # Note: volcano_layers includes its own scale_color_manual (for direction).
  # ring_layers includes scale_fill_gradient2 (for NES).
  # These use different aesthetics (color vs fill) so they coexist.

  p <- ggplot() +
    # Ring background arcs first (behind volcano)
    ring_layers$tick_bg +
    # Guide lines in tick ring
    ring_layers[grep("^guide_r_", names(ring_layers))] +
    # Volcano center (all layers)
    volcano_layers +
    # Protein ticks (on top of volcano edge)
    ring_layers$ticks +
    # Enrichment arcs (outermost ring)
    ring_layers$enrich_arcs +
    # NES fill scale
    ring_layers$fill_scale +
    # Labels
    label_layers +
    # Coordinate system and theme
    coord_fixed(
      xlim = c(-(label_r + 0.8), label_r + 0.8),
      ylim = c(-(label_r + 0.8), label_r + 2.0),
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(2, 2, 2, 2, "mm"),
      legend.position = "none"
    )

  # ── 5. Title ──────────────────────────────────────────────────────────────
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size, face = "bold",
                                      margin = margin(b = 5)))
  }

  # ── 5b. Contrast title/subtitle (rendered inside plot) ──────────────────
  if (!is.null(contrast_title)) {
    p <- p + annotate("text", x = 0, y = label_r + 1.8,
                       label = contrast_title,
                       size = title_size / .pt, fontface = "bold", hjust = 0.5)
    if (!is.null(contrast_subtitle)) {
      p <- p + annotate("text", x = 0, y = label_r + 1.2,
                         label = contrast_subtitle,
                         size = (title_size - 3) / .pt, fontface = "italic",
                         color = "grey30", hjust = 0.5)
    }
  }

  # ── 6. Attach metadata ────────────────────────────────────────────────────
  attr(p, "ring_data")  <- ring_data
  attr(p, "tick_data")  <- tick_data
  attr(p, "n_up")       <- attr(volcano_layers, "n_up")
  attr(p, "n_down")     <- attr(volcano_layers, "n_down")

  p
}

# ═══════════════════════════════════════════════════════════════════════════════
# build_panel_legend() — compact legend strip for individual panel exports
# ═══════════════════════════════════════════════════════════════════════════════
#' @return patchwork object: NES gradient bar | point legend
build_panel_legend <- function() {

  # Row 1: NES gradient bar (uses same clamp as ring arcs)
  nes_clamp <- 2.0
  gradient_df <- tibble(
    x   = seq(-nes_clamp, nes_clamp, length.out = 200),
    y   = 0,
    nes = seq(-nes_clamp, nes_clamp, length.out = 200)
  )

  p_gradient <- ggplot(gradient_df, aes(x = x, y = y, fill = nes)) +
    geom_tile(height = 0.6, width = diff(gradient_df$x[1:2]) * 1.05) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-nes_clamp, nes_clamp),
      guide = "none"
    ) +
    annotate("rect",
             xmin = -nes_clamp, xmax = nes_clamp,
             ymin = -0.3, ymax = 0.3,
             fill = NA, color = "grey40", linewidth = 0.3) +
    annotate("text", x = 0, y = 0.7,
             label = "Enrichment (NES)",
             size = KEY_TITLE, fontface = "bold", color = "grey20") +
    annotate("text", x = -nes_clamp * 0.85, y = -0.65,
             label = "Suppressed", size = KEY_TEXT, color = "#4393C3",
             fontface = "bold.italic") +
    annotate("text", x = nes_clamp * 0.85, y = -0.65,
             label = "Activated", size = KEY_TEXT, color = "#D6604D",
             fontface = "bold.italic") +
    coord_fixed(ratio = 1, clip = "off") +
    xlim(-nes_clamp * 1.3, nes_clamp * 1.3) +
    ylim(-1.2, 1.2) +
    theme_void()

  # Row 2: Horizontal point legend
  pt_df <- tibble(
    x     = c(0, 1.8, 3.6),
    y     = 0,
    color = c(DIR_COLORS["Up"], DIR_COLORS["Down"], NS_COLOR),
    label = c("Up (sig.)", "Down (sig.)", "Non-sig.")
  )

  p_points <- ggplot() +
    geom_point(data = pt_df, aes(x = x, y = y),
               color = pt_df$color, size = 2, alpha = 0.8) +
    geom_text(data = pt_df, aes(x = x, y = y, label = label),
              nudge_x = 0.25, hjust = 0, size = KEY_TEXT, color = "grey30") +
    xlim(-0.5, 6) +
    ylim(-0.8, 0.8) +
    theme_void()

  (p_gradient | p_points) + plot_layout(widths = c(1, 1.2))
}

# ═══════════════════════════════════════════════════════════════════════════════
# make_volcano_ring_pair() — two-panel assembly with shared legend
# ═══════════════════════════════════════════════════════════════════════════════
#' Build a side-by-side Young vs Old volcano-ring composite figure.
#'
#' Each contrast independently selects its top GO terms (balanced 5 up + 5 down
#' across Hallmark and GO:BP databases). Panels share the same arc geometry
#' layout but display different pathway sets. A shared legend strip is placed
#' below.
#'
#' @param de_df          Full combined DE results (wide format)
#' @param go_df          Full fGSEA results (all contrasts)
#' @param contrast_young Contrast name for the Young panel
#' @param contrast_old   Contrast name for the Old panel
#' @param n_terms        Number of top GO terms to display
#' @param title_young    Panel title for Young
#' @param title_old      Panel title for Old
#' @param output_dir     Base directory for figure outputs
#' @param save_outputs   Whether to save PDF/PNG and CSVs
#' @param ...            Additional arguments passed to make_volcano_ring()
#' @return Combined patchwork object (invisibly if save_outputs = TRUE)
make_volcano_ring_pair <- function(
    de_df,
    go_df,
    contrast_young = "Training_Young",
    contrast_old   = "Training_Old",
    n_terms        = 10,
    title_young    = "Training Effect (Young)",
    title_old      = "Training Effect (Old)",
    output_dir     = "04_Figures/F2",
    save_outputs   = TRUE,
    title_size           = 12,
    label_size           = 2.8,
    contrast_title_a     = NULL,
    contrast_subtitle_a  = NULL,
    contrast_title_b     = NULL,
    contrast_subtitle_b  = NULL,
    ...) {

  suppressPackageStartupMessages(library(cowplot))

  databases <- c("Hallmark", "GO:BP")

  # ── 1. Select per-panel GO terms (independent) ────────────────────────────
  # Each contrast independently picks its own top 5 up + 5 down terms from
  # Hallmark + GO:BP, balanced across databases. No shared/pooled selection.

  sig_young <- go_df %>%
    filter(contrast == contrast_young, database %in% databases, padj < 0.05)
  sig_old <- go_df %>%
    filter(contrast == contrast_old, database %in% databases, padj < 0.05)

  n_each <- ceiling(n_terms / 2)  # 5 up + 5 down

  # Helper: pick n terms from a significance-sorted pool, balanced across databases
  pick_direction <- function(sig_df, n_target) {
    pool <- sig_df %>% arrange(padj)
    interleave_db <- function(df, n) {
      hm <- df %>% filter(database == "Hallmark")
      bp <- df %>% filter(database == "GO:BP")
      n_hm <- min(nrow(hm), ceiling(n / 2))
      n_bp <- min(nrow(bp), n - n_hm)
      n_hm <- min(nrow(hm), n - n_bp)
      bind_rows(hm %>% slice_head(n = n_hm),
                bp %>% slice_head(n = n_bp))
    }
    interleave_db(pool, n_target) %>% slice_head(n = n_target)
  }

  # Young: top 5 up + 5 down
  top_terms_young <- bind_rows(
    pick_direction(sig_young %>% filter(NES > 0), n_each),
    pick_direction(sig_young %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_terms)

  # Old: top 5 up + 5 down
  top_terms_old <- bind_rows(
    pick_direction(sig_old %>% filter(NES > 0), n_each),
    pick_direction(sig_old %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_terms)

  message("make_volcano_ring_pair: Young — ", nrow(top_terms_young), " terms (",
          sum(top_terms_young$NES > 0), " up + ", sum(top_terms_young$NES < 0), " down)")
  message("make_volcano_ring_pair: Old — ", nrow(top_terms_old), " terms (",
          sum(top_terms_old$NES > 0), " up + ", sum(top_terms_old$NES < 0), " down)")

  if (nrow(top_terms_young) == 0) {
    stop("No significant terms found for contrast '", contrast_young,
         "' — cannot build paired panels")
  }
  if (nrow(top_terms_old) == 0) {
    stop("No significant terms found for contrast '", contrast_old,
         "' — cannot build paired panels")
  }

  # ── 2. Build independent ring geometry for each panel ────────────────────
  build_ring_for_contrast <- function(top_terms, contrast_name, go_df, databases) {
    real_rows <- go_df %>%
      filter(contrast == contrast_name, pathway %in% top_terms$pathway)
    go_subset <- real_rows %>%
      mutate(padj = match(pathway, top_terms$pathway) * 1e-10)
    ring <- prepare_ring_data(
      go_df        = go_subset,
      contrast     = contrast_name,
      n_terms      = nrow(top_terms),
      gap_degrees  = 3,
      start_offset = 0,
      databases    = databases
    )

    # Post-process: widen gaps at UP/DOWN hemisphere split points
    n <- nrow(ring)
    if (n >= 2) {
      gap_normal <- 3
      gap_split  <- 8
      start_offset <- 0

      # Gap pattern: split gaps at positions n_each (UP→DOWN) and n (last→first)
      # e.g. for 10 terms: [3, 3, 3, 3, 8, 3, 3, 3, 3, 8]
      gaps <- rep(gap_normal, n)
      gaps[n_each] <- gap_split   # between last UP and first DOWN
      gaps[n]      <- gap_split   # between last DOWN and first UP (wraps at 6 o'clock)

      total_gap <- sum(gaps)
      arc_width_deg <- (360 - total_gap) / n

      # Recalculate start/end for each arc
      cum_offset <- 0
      for (i in seq_len(n)) {
        gap_before <- if (i == 1) gaps[n] else gaps[i - 1]
        if (i == 1) {
          cum_offset <- 0
        } else {
          cum_offset <- cum_offset + arc_width_deg + gaps[i - 1]
        }
        ring$start_deg[i] <- start_offset + cum_offset
        ring$end_deg[i]   <- ring$start_deg[i] + arc_width_deg
        ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
        ring$start_rad[i] <- ring$start_deg[i] * pi / 180
        ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
        ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
      }
    }

    ring
  }

  ring_data_young <- build_ring_for_contrast(top_terms_young, contrast_young, go_df, databases)
  ring_data_old   <- build_ring_for_contrast(top_terms_old, contrast_old, go_df, databases)

  message("  Ring geometry: Young ", nrow(ring_data_young), " arcs, Old ", nrow(ring_data_old), " arcs")

  # ── 3. Create both panels ─────────────────────────────────────────────────
  message("  Building Young panel...")
  p_young <- make_volcano_ring(
    de_df              = de_df,
    go_df              = go_df,
    contrast           = contrast_young,
    title              = title_young,
    title_size         = title_size,
    label_size         = label_size,
    contrast_title     = contrast_title_a,
    contrast_subtitle  = contrast_subtitle_a,
    ring_data_override = ring_data_young,
    ...
  )

  message("  Building Old panel...")
  p_old <- make_volcano_ring(
    de_df              = de_df,
    go_df              = go_df,
    contrast           = contrast_old,
    title              = title_old,
    title_size         = title_size,
    label_size         = label_size,
    contrast_title     = contrast_title_b,
    contrast_subtitle  = contrast_subtitle_b,
    ring_data_override = ring_data_old,
    ...
  )

  # ── 4. Build shared legend ─────────────────────────────────────────────────
  # 4a. NES gradient bar (clamped to match ring arcs)
  nes_clamp <- 2.0

  gradient_df <- tibble(
    x   = seq(-nes_clamp, nes_clamp, length.out = 200),
    y   = 0,
    nes = seq(-nes_clamp, nes_clamp, length.out = 200)
  )

  p_gradient <- ggplot(gradient_df, aes(x = x, y = y, fill = nes)) +
    geom_tile(height = 0.6, width = diff(gradient_df$x[1:2]) * 1.05) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-nes_clamp, nes_clamp),
      guide = "none"
    ) +
    # Border around the gradient bar
    annotate("rect",
             xmin = -nes_clamp, xmax = nes_clamp,
             ymin = -0.3, ymax = 0.3,
             fill = NA, color = "grey40", linewidth = 0.3) +
    # Labels
    annotate("text", x = 0, y = 0.7,
             label = "Enrichment Direction (GSEA NES)",
             size = label_size + 0.7, fontface = "bold", color = "grey20") +
    annotate("text", x = -nes_clamp * 0.85, y = -0.65,
             label = "Suppressed", size = label_size, color = "#4393C3",
             fontface = "bold.italic") +
    annotate("text", x = 0, y = -0.65,
             label = "No Significant Change", size = label_size, color = "grey50",
             fontface = "bold.italic") +
    annotate("text", x = nes_clamp * 0.85, y = -0.65,
             label = "Activated", size = label_size, color = "#D6604D",
             fontface = "bold.italic") +
    coord_fixed(ratio = 1, clip = "off") +
    xlim(-nes_clamp * 1.3, nes_clamp * 1.3) +
    ylim(-1.2, 1.2) +
    theme_void()

  # 4b. Point color legend (stacked vertically for readability)
  point_legend_df <- tibble(
    x     = c(0, 0, 0),
    y     = c(0.5, 0, -0.5),
    color = c(DIR_COLORS["Up"], DIR_COLORS["Down"], NS_COLOR),
    label = c("Significantly Upregulated Protein",
              "Significantly Downregulated Protein",
              "Non-significant Protein")
  )

  p_points_legend <- ggplot() +
    # Draw circles
    geom_point(data = point_legend_df,
               aes(x = x, y = y),
               color = point_legend_df$color,
               size = 2.5, alpha = 0.8) +
    # Labels next to circles
    geom_text(data = point_legend_df,
              aes(x = x, y = y, label = label),
              nudge_x = 0.15, hjust = 0, size = label_size, color = "grey30") +
    xlim(-0.3, 4) +
    ylim(-1.0, 1.0) +
    theme_void()

  # Combine legend elements horizontally
  shared_legend <- (p_gradient | p_points_legend) +
    plot_layout(widths = c(1, 1.2))

  # ── 5. Assemble full composite ─────────────────────────────────────────────
  combined <- (p_young | p_old) / shared_legend +
    plot_layout(heights = c(8, 1.5))

  # ── 6. Export ring data CSVs ───────────────────────────────────────────────
  if (save_outputs) {
    data_dir   <- file.path(output_dir, "c_data", "panel_A")
    report_dir <- file.path(output_dir, "b_reports")

    dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

    # CSV export: remove list-columns (gene_list) for flat file
    ring_young_export <- ring_data_young %>%
      dplyr::select(-gene_list)
    ring_old_export <- ring_data_old %>%
      dplyr::select(-gene_list)

    write_csv(ring_young_export, file.path(data_dir, "ring_terms.csv"))
    write_csv(ring_old_export,   file.path(data_dir, "ring_terms_old.csv"))
    message("  Exported: ", file.path(data_dir, "ring_terms.csv"))
    message("  Exported: ", file.path(data_dir, "ring_terms_old.csv"))

    # ── 7. Save figures ─────────────────────────────────────────────────────

    # Use cairo_pdf if available, otherwise fall back to standard pdf
    pdf_device <- tryCatch(
      { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
      error = function(e) "pdf"
    )

    # 7a. Panel A (Young) and Panel B (Old) — individual ring composites + legend
    panel_w <- 160   # mm width
    panel_h <- 180   # mm height (extra 20mm for legend strip)

    # Build per-panel legends (uses clamped NES scale internally)
    legend_young <- build_panel_legend()
    legend_old   <- build_panel_legend()

    p_young_with_key <- p_young / legend_young + plot_layout(heights = c(6, 1))
    p_old_with_key   <- p_old   / legend_old   + plot_layout(heights = c(6, 1))

    ggsave(file.path(report_dir, "panel_A_volcano.pdf"), p_young_with_key,
           width = panel_w, height = panel_h, units = "mm", device = pdf_device)
    ggsave(file.path(report_dir, "panel_A_volcano.png"), p_young_with_key,
           width = panel_w, height = panel_h, units = "mm", dpi = 300)
    message("  Saved: ", file.path(report_dir, "panel_A_volcano.pdf"))
    message("  Saved: ", file.path(report_dir, "panel_A_volcano.png"))

    ggsave(file.path(report_dir, "panel_B_volcano.pdf"), p_old_with_key,
           width = panel_w, height = panel_h, units = "mm", device = pdf_device)
    ggsave(file.path(report_dir, "panel_B_volcano.png"), p_old_with_key,
           width = panel_w, height = panel_h, units = "mm", dpi = 300)
    message("  Saved: ", file.path(report_dir, "panel_B_volcano.pdf"))
    message("  Saved: ", file.path(report_dir, "panel_B_volcano.png"))
  }

  # ── Return ─────────────────────────────────────────────────────────────────
  attr(combined, "ring_data_young") <- ring_data_young
  attr(combined, "ring_data_old")   <- ring_data_old
  attr(combined, "p_young")         <- p_young
  attr(combined, "p_old")           <- p_old

  invisible(combined)
}
