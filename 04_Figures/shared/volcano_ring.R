#!/usr/bin/env Rscript
# volcano_ring.R — Circular volcano-in-ring composite plot utility
# Standard Cartesian ggplot with ggforce::geom_arc_bar(); NO coord_polar().
# Sources style.R for palettes and sizing constants.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
  library(patchwork)
  library(scales)
})

if (!exists("FIG_THEME")) source("04_Figures/shared/style.R")

clean_ring_label <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    # Capitalisation fixes
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
    # Smart abbreviations for ring labels
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
    str_wrap(width = 18)
}

prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms      = 12,
                              gap_degrees  = 3,
                              start_offset = 0,
                              databases    = c("Hallmark", "GO:BP")) {

  ring <- go_df %>%
    filter(contrast == !!contrast,
           database %in% databases,
           padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = n_terms)

  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms for '", contrast, "'")
    return(tibble())
  }

  arc_width_deg <- (360 - n * gap_degrees) / n

  ring %>%
    mutate(
      term_idx    = row_number(),
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg   * pi / 180,
      mid_rad     = mid_deg   * pi / 180,
      clean_label = clean_ring_label(pathway),
      gene_list   = str_split(leadingEdge, ";")
    )
}

build_tick_data <- function(ring_data,
                            de_df,
                            contrast,
                            tick_r0 = 4.0,
                            tick_r1 = 5.2) {

  if (nrow(ring_data) == 0) return(tibble())

  logfc_col <- paste0("logFC_", contrast)

  gene_lfc <- de_df %>%
    dplyr::select(gene, lfc = all_of(logfc_col)) %>%
    filter(!is.na(lfc)) %>%
    distinct(gene, .keep_all = TRUE)

  pad_rad <- 0.5 * pi / 180

  map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) return(tibble())

    arc_start <- row$start_rad + pad_rad
    arc_end   <- row$end_rad   - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad

    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)

    matched <- gene_lfc %>%
      filter(gene %in% genes_in_arc) %>%
      `[`(match(genes_in_arc, .$gene), ) %>%
      filter(!is.na(gene))

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
}

build_volcano_layers <- function(de_df,
                                 contrast,
                                 volcano_radius = 3.5,
                                 fc_thresh      = log2(1.5),
                                 p_thresh       = 0.05,
                                 up_color       = DIR_COLORS["Up"],
                                 down_color     = DIR_COLORS["Down"],
                                 ns_color       = DIR_COLORS["NS"],
                                 point_size     = 0.6,
                                 point_alpha    = 0.5,
                                 count_label_size = 2.8) {

  logfc_col <- paste0("logFC_", contrast)
  pval_col  <- paste0("P.Value_", contrast)
  pi_col    <- paste0("pi_score_", contrast)

  vdf <- de_df %>%
    transmute(
      gene       = gene,
      logFC      = .data[[logfc_col]],
      pvalue     = .data[[pval_col]],
      pi_score   = .data[[pi_col]],
      neg_log10p = -log10(pvalue)
    ) %>%
    filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p)) %>%
    mutate(
      direction = case_when(
        pi_score < 0.05 & logFC > 0 ~ "Up",
        pi_score < 0.05 & logFC < 0 ~ "Down",
        TRUE                         ~ "NS"
      )
    )

  n_up   <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")

  x_data_max <- max(abs(vdf$logFC), na.rm = TRUE)
  y_data_max <- max(vdf$neg_log10p, na.rm = TRUE)

  margin <- 0.92
  vr <- volcano_radius * margin

  scale_x_fn <- function(val) val / x_data_max * vr
  scale_y_fn <- function(val) val / y_data_max * vr

  vdf <- vdf %>%
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  vdf_ns  <- vdf %>% filter(direction == "NS")
  vdf_sig <- vdf %>% filter(direction != "NS")

  layers <- list(
    ns_points = geom_point(
      data = vdf_ns, aes(x = x_plot, y = y_plot),
      color = ns_color, size = point_size * 0.9, alpha = point_alpha * 0.6,
      inherit.aes = FALSE
    ),
    sig_points = geom_point(
      data = vdf_sig, aes(x = x_plot, y = y_plot, color = direction),
      size = point_size * 1.3, alpha = point_alpha * 1.3, stroke = 0.3,
      inherit.aes = FALSE
    ),
    color_scale = scale_color_manual(
      values = c(Up = unname(up_color), Down = unname(down_color)),
      guide  = "none"
    ),
    x_axis_line = annotate(
      "segment", x = -vr * 0.42, xend = vr * 0.42, y = -vr, yend = -vr,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "both", length = unit(1.2, "mm"), type = "closed")
    ),
    x_axis_up = annotate(
      "text", x = vr * 0.45, y = -vr,
      label = "up", size = count_label_size * 0.6,
      color = unname(up_color), fontface = "italic", hjust = 0
    ),
    x_axis_down = annotate(
      "text", x = -vr * 0.45, y = -vr,
      label = "down", size = count_label_size * 0.6,
      color = unname(down_color), fontface = "italic", hjust = 1
    ),
    x_axis_label = annotate(
      "text", x = 0, y = -vr - 0.35,
      label = "log2 FC", size = count_label_size * 0.6,
      color = "grey40", fontface = "italic"
    ),
    y_axis_line = annotate(
      "segment", x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "last", length = unit(1.2, "mm"), type = "closed")
    ),
    y_axis_label = annotate(
      "text", x = 0, y = vr * 1.04,
      label = expression(-log[10]~italic(p)), size = count_label_size * 0.6,
      color = "grey40"
    ),
    n_up_box = annotate(
      "label", x = vr * 0.5, y = vr * 0.85,
      label = n_up, size = count_label_size,
      color = "black", fill = alpha(up_color, 0.9), fontface = "bold",
      label.padding = unit(2.5, "pt"), label.r = unit(2, "pt"),
      linewidth = 0.4
    ),
    n_up_text = annotate(
      "text", x = vr * 0.5, y = vr * 0.85,
      label = n_up, size = count_label_size,
      color = "white", fontface = "bold"
    ),
    n_down_box = annotate(
      "label", x = -vr * 0.5, y = vr * 0.85,
      label = n_down, size = count_label_size,
      color = "black", fill = alpha(down_color, 0.9), fontface = "bold",
      label.padding = unit(2.5, "pt"), label.r = unit(2, "pt"),
      linewidth = 0.4
    ),
    n_down_text = annotate(
      "text", x = -vr * 0.5, y = vr * 0.85,
      label = n_down, size = count_label_size,
      color = "white", fontface = "bold"
    )
  )

  attr(layers, "x_data_max") <- x_data_max
  attr(layers, "y_data_max") <- y_data_max
  attr(layers, "n_up")       <- n_up
  attr(layers, "n_down")     <- n_down
  attr(layers, "scale_x")    <- scale_x_fn
  attr(layers, "scale_y")    <- scale_y_fn

  layers
}

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

  layers$tick_bg <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
        start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15,
    inherit.aes = FALSE
  )

  guide_radii <- seq(tick_r0, tick_r1, length.out = 5)[2:4]
  for (gr in guide_radii) {
    layers[[paste0("guide_r_", round(gr, 2))]] <- geom_arc(
      data = ring_data,
      aes(x0 = 0, y0 = 0, r = gr, start = start_rad, end = end_rad),
      color = "grey82", linewidth = 0.1,
      inherit.aes = FALSE
    )
  }

  if (nrow(tick_data) > 0) {
    layers$ticks <- geom_segment(
      data = tick_data,
      aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  layers$enrich_arcs <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1,
        start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )

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

build_label_layer <- function(ring_data,
                              label_r    = 7.0,
                              label_size = 2.8) {

  if (nrow(ring_data) == 0) return(list())

  label_df <- ring_data %>%
    mutate(
      arc_span   = end_rad - start_rad,
      cx         = (cos(start_rad) - cos(end_rad)) / arc_span,
      cy         = (sin(end_rad) - sin(start_rad)) / arc_span,
      centroid_r = sqrt(cx^2 + cy^2),
      x_label    = label_r * cx / centroid_r,
      y_label    = label_r * cy / centroid_r
    )

  list(
    labels = geom_label(
      data = label_df,
      aes(x = x_label, y = y_label, label = clean_label),
      hjust = 0.5, vjust = 0.5,
      size = label_size * 0.65, color = "grey20", fontface = "bold",
      angle = 0, fill = "grey96", alpha = 0.85, lineheight = 0.8,
      linewidth = 0.15, label.padding = unit(1.5, "pt"),
      label.r = unit(0.1, "lines"),
      inherit.aes = FALSE
    )
  )
}

# min_size excludes small gene sets prone to tissue-irrelevant GO artifacts
# (Reimand et al. 2019, Nat Protocols S3.4).
select_ring_terms <- function(go_df, contrast_name, n_each = 6,
                              databases = c("Hallmark", "GO:BP"),
                              min_size = 15) {
  sig <- go_df %>%
    filter(contrast == contrast_name, database %in% databases,
           padj < 0.05, size >= min_size)

  pick_direction <- function(sig_df, n_target) {
    pool <- sig_df %>% arrange(padj)
    hm <- pool %>% filter(database == "Hallmark")
    bp <- pool %>% filter(database == "GO:BP")
    n_hm <- min(nrow(hm), ceiling(n_target / 2))
    n_bp <- min(nrow(bp), n_target - n_hm)
    n_hm <- min(nrow(hm), n_target - n_bp)
    bind_rows(hm %>% slice_head(n = n_hm),
              bp %>% slice_head(n = n_bp)) %>%
      slice_head(n = n_target)
  }

  bind_rows(
    pick_direction(sig %>% filter(NES > 0), n_each),
    pick_direction(sig %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_each * 2)
}

center_ring_angles <- function(ring, n_up) {
  n <- nrow(ring)
  if (n < 2 || n_up < 1) return(ring)
  # Midpoint of the Up block = average of first arc start and last Up arc end
  up_mid <- (ring$start_deg[1] + ring$end_deg[min(n_up, n)]) / 2
  offset <- 90 - up_mid   # center Up block at 90° (right side)
  ring$start_deg <- ring$start_deg + offset
  ring$end_deg   <- ring$end_deg   + offset
  ring$mid_deg   <- ring$mid_deg   + offset
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring
}

build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = 6,
                                 databases = c("Hallmark", "GO:BP")) {
  real_rows <- go_df %>%
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  go_subset <- real_rows %>%
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3, start_offset = 0,
    databases = databases
  )

  n <- nrow(ring)
  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    gaps[min(n_each, n)] <- gap_split
    gaps[n]              <- gap_split
    arc_width_deg <- (360 - sum(gaps)) / n
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_width_deg + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_width_deg
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
    ring <- center_ring_angles(ring, min(n_each, n))
  }
  ring
}

make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title              = NULL,
                              contrast_title     = NULL,
                              contrast_subtitle  = NULL,
                              title_size         = 12,
                              n_terms            = 12,
                              gap_degrees        = 3,
                              start_offset       = 0,
                              databases          = c("Hallmark", "GO:BP"),
                              volcano_radius     = 3.5,
                              tick_r0            = 4.0,
                              tick_r1            = 5.2,
                              arc_r0             = 5.2,
                              arc_r1             = 6.0,
                              label_r            = 6.5,
                              fc_thresh          = log2(1.5),
                              p_thresh           = 0.05,
                              up_color           = DIR_COLORS["Up"],
                              down_color         = DIR_COLORS["Down"],
                              ns_color           = DIR_COLORS["NS"],
                              point_size         = 0.6,
                              point_alpha        = 0.5,
                              label_size         = 2.8,
                              count_label_size   = 2.8,
                              ring_data_override = NULL) {

  if (!is.null(ring_data_override)) {
    ring_data <- ring_data_override
  } else {
    ring_data <- prepare_ring_data(
      go_df = go_df, contrast = contrast, n_terms = n_terms,
      gap_degrees = gap_degrees, start_offset = start_offset,
      databases = databases
    )
  }

  tick_data <- build_tick_data(
    ring_data = ring_data, de_df = de_df, contrast = contrast,
    tick_r0 = tick_r0, tick_r1 = tick_r1
  )

  volcano_layers <- build_volcano_layers(
    de_df = de_df, contrast = contrast, volcano_radius = volcano_radius,
    fc_thresh = fc_thresh, p_thresh = p_thresh,
    up_color = up_color, down_color = down_color, ns_color = ns_color,
    point_size = point_size, point_alpha = point_alpha,
    count_label_size = count_label_size
  )

  ring_layers <- build_ring_layers(
    ring_data = ring_data, tick_data = tick_data,
    tick_r0 = tick_r0, tick_r1 = tick_r1,
    arc_r0 = arc_r0, arc_r1 = arc_r1,
    up_color = up_color, down_color = down_color
  )

  label_layers <- build_label_layer(
    ring_data = ring_data, label_r = label_r, label_size = label_size
  )

  p <- ggplot() +
    ring_layers$tick_bg +
    ring_layers[grep("^guide_r_", names(ring_layers))] +
    volcano_layers +
    ring_layers$ticks +
    ring_layers$enrich_arcs +
    ring_layers$fill_scale +
    label_layers +
    coord_fixed(
      xlim = c(-(label_r + 1.5), label_r + 1.5),
      ylim = c(-(label_r + 1.5), label_r + 1.8),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.margin = margin(1, 2, 2, 2, "mm"), legend.position = "none")

  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size, face = "bold",
                                      margin = margin(b = 1)))
  }

  # Contrast title/subtitle rendered inside the plot (above outermost ring)
  if (!is.null(contrast_title)) {
    p <- p + annotate("text", x = 0, y = label_r + 1.5,
                       label = contrast_title,
                       size = title_size / .pt, fontface = "bold", hjust = 0.5)
    if (!is.null(contrast_subtitle)) {
      p <- p + annotate("text", x = 0, y = label_r + 0.7,
                         label = contrast_subtitle,
                         size = (title_size - 3) / .pt, fontface = "italic",
                         color = "grey30", hjust = 0.5)
    }
  }

  attr(p, "ring_data")  <- ring_data
  attr(p, "tick_data")  <- tick_data
  attr(p, "n_up")       <- attr(volcano_layers, "n_up")
  attr(p, "n_down")     <- attr(volcano_layers, "n_down")

  p
}

make_volcano_ring_pair <- function(
    de_df,
    go_df,
    contrast_young       = "Training_Young",
    contrast_old         = "Training_Old",
    n_terms              = 10,
    title_young          = "Training Effect (Young)",
    title_old            = "Training Effect (Old)",
    title_size           = 12,
    label_size           = 2.8,
    count_label_size     = 2.8,
    contrast_title_a     = NULL,
    contrast_subtitle_a  = NULL,
    contrast_title_b     = NULL,
    contrast_subtitle_b  = NULL,
    output_dir           = "04_Figures/F2",
    save_outputs         = TRUE,
    ...) {

  databases <- c("Hallmark", "GO:BP")

  sig_young <- go_df %>%
    filter(contrast == contrast_young, database %in% databases,
           padj < 0.05, size >= 15)
  sig_old <- go_df %>%
    filter(contrast == contrast_old, database %in% databases,
           padj < 0.05, size >= 15)

  n_each <- ceiling(n_terms / 2)

  pick_direction <- function(sig_df, n_target) {
    pool <- sig_df %>% arrange(padj)
    hm <- pool %>% filter(database == "Hallmark")
    bp <- pool %>% filter(database == "GO:BP")
    n_hm <- min(nrow(hm), ceiling(n_target / 2))
    n_bp <- min(nrow(bp), n_target - n_hm)
    n_hm <- min(nrow(hm), n_target - n_bp)
    bind_rows(hm %>% slice_head(n = n_hm),
              bp %>% slice_head(n = n_bp)) %>%
      slice_head(n = n_target)
  }

  top_terms_young <- bind_rows(
    pick_direction(sig_young %>% filter(NES > 0), n_each),
    pick_direction(sig_young %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_terms)

  top_terms_old <- bind_rows(
    pick_direction(sig_old %>% filter(NES > 0), n_each),
    pick_direction(sig_old %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_terms)

  if (nrow(top_terms_young) == 0)
    stop("No significant terms for '", contrast_young, "'")
  if (nrow(top_terms_old) == 0)
    stop("No significant terms for '", contrast_old, "'")

  ring_data_young <- build_ring_with_gaps(top_terms_young, contrast_young, go_df, n_each, databases)
  ring_data_old   <- build_ring_with_gaps(top_terms_old, contrast_old, go_df, n_each, databases)

  p_young <- make_volcano_ring(
    de_df = de_df, go_df = go_df, contrast = contrast_young,
    title = title_young, title_size = title_size, label_size = label_size,
    count_label_size = count_label_size,
    contrast_title = contrast_title_a, contrast_subtitle = contrast_subtitle_a,
    ring_data_override = ring_data_young, ...
  )
  p_old <- make_volcano_ring(
    de_df = de_df, go_df = go_df, contrast = contrast_old,
    title = title_old, title_size = title_size, label_size = label_size,
    count_label_size = count_label_size,
    contrast_title = contrast_title_b, contrast_subtitle = contrast_subtitle_b,
    ring_data_override = ring_data_old, ...
  )

  if (save_outputs) {
    data_dir   <- file.path(output_dir, "c_data", "panel_A")
    report_dir <- file.path(output_dir, "b_reports")
    dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

    write_csv(ring_data_young %>% dplyr::select(-gene_list),
              file.path(data_dir, "ring_terms.csv"))
    write_csv(ring_data_old %>% dplyr::select(-gene_list),
              file.path(data_dir, "ring_terms_old.csv"))

    pdf_device <- get_pdf_device()
    panel_w <- 160
    panel_h <- 160   # no legend strip — panels only

    ggsave(file.path(report_dir, "panel_A_volcano.pdf"), p_young,
           width = panel_w, height = panel_h, units = "mm", device = pdf_device)
    ggsave(file.path(report_dir, "panel_A_volcano.png"), p_young,
           width = panel_w, height = panel_h, units = "mm", dpi = 300)
    ggsave(file.path(report_dir, "panel_B_volcano.pdf"), p_old,
           width = panel_w, height = panel_h, units = "mm", device = pdf_device)
    ggsave(file.path(report_dir, "panel_B_volcano.png"), p_old,
           width = panel_w, height = panel_h, units = "mm", dpi = 300)
  }

  result <- list(p_young = p_young, p_old = p_old,
                 ring_data_young = ring_data_young,
                 ring_data_old = ring_data_old)
  invisible(result)
}
