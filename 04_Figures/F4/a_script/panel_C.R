################################################################################
#   Figure 4 — Panel C: Per-cluster Sankey triptych (heatmap | Sankey | bars)
#
#   Requires from setup: core_proteins, group_z, GROUP_COLS, GROUP_LABS,
#                         CLUSTER_COLORS, DB_COLORS, enrich_top,
#                         protein_pathway_links, optimal_k, THEME_PUB,
#                         clean_pathway_name, make_sigmoid_ribbon,
#                         reorder_within, scale_y_reordered,
#                         RPT_DIR, FIG_W, FIG_H, COL_WIDTHS, row_heights
#   Outputs: panels_C (list of ggplot objects), build_panel_C (function)
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Enrichment bars:
#    - Bars show -log10(p.adjust) from ORA. BH correction is applied within
#      each enricher() call in setup (per database per cluster). This is
#      the standard approach: Hallmark, GO:BP, GO:CC are separate hypothesis
#      families. Cross-database global correction is not standard for ORA
#      and would be overly conservative.                                PASS
#
# 2. Gene count annotation:
#    - Numeric labels on bars show the number of cluster genes in each
#      pathway (Count). This is purely descriptive.                     PASS
#
# 3. Heatmap z-score cap:
#    - Z_CAP = 2 truncates extreme values for color visualization. This
#      is a display parameter, not a statistical operation.             PASS
#
# 4. Sankey (1:1 greedy assignment):
#    - Visualization heuristic only. Each protein assigned to its most
#      significant pathway (lowest p.adjust) containing it. Not a
#      statistical claim. Documented in setup audit.                    PASS
#
# 5. Top 3 per database per cluster:
#    - slice_min(p.adjust, n = 3) selects the 3 most significant terms
#      per database per cluster after BH correction and rrvgo reduction.
#      This is a display filter, not an additional statistical test.    PASS
#
# 6. Multiple comparison scope summary:
#    - Level 1: Within each enricher() call, BH across all terms in that
#      database for that cluster.
#    - Level 2: rrvgo removes redundant GO terms (semantic similarity
#      > 0.85), keeping parent representatives. This is a conservative
#      filtering step that reduces the displayed term count.
#    - Level 3: Top 3 per database is a display cutoff.
#    - No additional cross-cluster correction applied. This is standard:
#      each cluster is an independent biological question.              PASS
# ---------------------------------------------------------------------------

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building Panel C: per-cluster Sankey triptych ===\n")

Z_CAP <- 2

# --- Shared x-axis max for enrichment bars across all clusters ---
shared_x_max_C <- 30   # fixed cap; break indicator for bars exceeding this

cat(sprintf("  Shared x-axis max (enrichment bars): %.2f\n", shared_x_max_C))

# --- build_panel_C function ---
build_panel_C <- function(cl_id, show_xlab, shared_x_max) {

  cl_label <- paste0("C", cl_id)
  cl_col   <- CLUSTER_COLORS[cl_label]

  # -- Genes for this cluster, sorted by descending membership --
  cl_core <- core_proteins %>%
    filter(cluster == cl_label) %>%
    arrange(desc(membership))
  cl_genes <- cl_core$gene

  n_genes <- length(cl_genes)
  cat(sprintf("  Panel C — %s: %d core genes\n", cl_label, n_genes))

  # -- Links and enrichment for this cluster --
  cl_links <- protein_pathway_links %>%
    filter(cluster == cl_label)
  cl_enrich <- enrich_top %>%
    filter(cluster == cl_label) %>%
    arrange(p.adjust)

  # Mapped pathways (exclude Unmapped)
  mapped_links <- cl_links %>% filter(pathway != "Unmapped")
  mapped_pathways <- unique(mapped_links$pathway)
  n_pw <- length(mapped_pathways)

  cat(sprintf("    Mapped pathways: %d, Mapped proteins: %d, Unmapped: %d\n",
              n_pw, nrow(mapped_links), sum(cl_links$pathway == "Unmapped")))

  # ===========================================================================
  # Sub-panel C1: Heatmap Strip
  # ===========================================================================

  ht_genes <- intersect(cl_genes, rownames(group_z))
  ht_mat   <- group_z[ht_genes, , drop = FALSE]

  ht_mat[ht_mat >  Z_CAP] <-  Z_CAP
  ht_mat[ht_mat < -Z_CAP] <- -Z_CAP

  ht_long <- as.data.frame(ht_mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = all_of(GROUP_COLS),
                 names_to = "group", values_to = "z") %>%
    mutate(
      gene  = factor(gene, levels = rev(ht_genes)),
      group = factor(group, levels = GROUP_COLS)
    )

  ht_colors <- c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")

  p_ht <- ggplot(ht_long, aes(x = group, y = gene, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = ht_colors,
      limits  = c(-Z_CAP, Z_CAP),
      oob     = scales::squish,
      guide   = "none"
    ) +
    scale_x_discrete(
      labels = if (show_xlab) GROUP_LABS else NULL,
      expand = expansion(0)
    ) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = NULL, y = NULL) +
    THEME_PUB +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.x     = if (show_xlab) element_text(size = 5, lineheight = 0.85)
                          else element_blank(),
      axis.ticks.x    = if (show_xlab) element_line() else element_blank(),
      axis.ticks.length = unit(0, "pt"),
      panel.border    = element_rect(colour = cl_col, linewidth = 0.6, fill = NA),
      legend.position = "none",
      plot.margin     = margin(t = 1, r = -5, b = if (show_xlab) 4 else 1, l = 0)
    )

  # ===========================================================================
  # Sub-panel C2: Sigmoid Sankey
  # ===========================================================================

  if (n_pw == 0) {
    p_sankey <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No mapped\npathways",
               size = 2.5, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 1, l = 0))
  } else {

    Y_SPAN  <- n_genes
    X_GENE  <- 0.05    # gene bars flush at left edge
    X_PW    <- 1.55    # pathway bars
    BAR_W   <- 0.06

    # -- Gene bar positions (top to bottom) --
    gene_h   <- Y_SPAN / (n_genes * 1.15)
    gene_gap <- if (n_genes > 1) (Y_SPAN - n_genes * gene_h) / (n_genes - 1) else 0

    gene_bars <- tibble(
      gene  = cl_genes,
      idx   = seq_along(cl_genes),
      y_top = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
      y_bot = y_top - gene_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left  = X_GENE - BAR_W / 2,
      x_right = X_GENE + BAR_W / 2,
      fill  = cl_col
    )

    # -- Pathway bar positions (top to bottom) --
    pw_h   <- Y_SPAN / (n_pw * 1.4)
    pw_gap <- if (n_pw > 1) (Y_SPAN - n_pw * pw_h) / (n_pw - 1) else 0

    # Order pathways to match enrichment order (lowest p.adjust first at top)
    pw_order <- cl_enrich %>%
      filter(Description %in% mapped_pathways) %>%
      pull(Description)
    # Add any mapped_pathways not in cl_enrich at the end
    pw_order <- c(pw_order, setdiff(mapped_pathways, pw_order))

    pw_bars <- tibble(
      pathway = pw_order,
      idx     = seq_along(pw_order),
      y_top   = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
      y_bot   = y_top - pw_h,
      y_ctr   = (y_top + y_bot) / 2,
      x_left  = X_PW - BAR_W / 2,
      x_right = X_PW + BAR_W / 2
    ) %>%
      left_join(
        mapped_links %>% dplyr::select(pathway, database) %>% distinct(),
        by = "pathway"
      ) %>%
      mutate(fill = DB_COLORS[database])

    # -- Slot allocation within each pathway bar --
    # Each pathway bar is subdivided into slots (one per gene mapping to it)
    slot_data <- mapped_links %>%
      group_by(pathway) %>%
      mutate(slot_idx = row_number(), n_slots = n()) %>%
      ungroup() %>%
      left_join(pw_bars %>% dplyr::select(pathway, y_top, y_bot), by = "pathway") %>%
      mutate(
        slot_h   = (y_top - y_bot) / n_slots,
        slot_top = y_top - (slot_idx - 1) * slot_h,
        slot_bot = slot_top - slot_h,
        slot_ctr = (slot_top + slot_bot) / 2
      )

    # -- Build ribbon data --
    ribbon_input <- slot_data %>%
      left_join(gene_bars %>% dplyr::select(gene, y_top_gene = y_top,
                                             y_bot_gene = y_bot, y_ctr_gene = y_ctr),
                by = "gene") %>%
      mutate(
        ribbon_id = paste0(cl_label, "_", gene, "_", pathway),
        cross_dist = abs(y_ctr_gene - slot_ctr)
      ) %>%
      arrange(desc(cross_dist))  # far-crossing ribbons drawn first (behind)

    ribbon_polys <- pmap_dfr(ribbon_input, function(gene, pathway, database, cluster,
                                                     slot_idx, n_slots,
                                                     y_top, y_bot,
                                                     slot_h, slot_top, slot_bot, slot_ctr,
                                                     y_top_gene, y_bot_gene, y_ctr_gene,
                                                     ribbon_id, cross_dist) {
      make_sigmoid_ribbon(
        x0 = X_GENE + BAR_W / 2,
        x1 = X_PW   - BAR_W / 2,
        y0_top = y_top_gene,
        y0_bot = y_bot_gene,
        y1_top = slot_top,
        y1_bot = slot_bot,
        ribbon_id = ribbon_id
      ) %>%
        mutate(fill = DB_COLORS[database])
    })

    # -- Build the plot --
    # Gene bar rectangles as data frame for geom_rect
    gene_rect <- gene_bars %>%
      dplyr::select(x_left, x_right, y_bot, y_top, fill)
    pw_rect <- pw_bars %>%
      dplyr::select(x_left, x_right, y_bot, y_top, fill)

    p_sankey <- ggplot() +
      # Ribbons (drawn first = behind bars)
      geom_polygon(
        data = ribbon_polys,
        aes(x = x, y = y, group = ribbon_id, fill = fill),
        alpha = 0.25, color = NA
      ) +
      # Gene bars
      geom_rect(
        data = gene_rect,
        aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
        color = NA
      ) +
      # Pathway bars
      geom_rect(
        data = pw_rect,
        aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
        color = NA
      ) +
      # Pathway labels: LEFT of pathway bars
      geom_text(
        data = pw_bars,
        aes(x = x_left - 0.04, y = y_ctr,
            label = str_trunc(clean_pathway_name(pathway), 40, ellipsis = "...")),
        hjust = 1, size = 2.5, fontface = "bold", color = "grey20"
      ) +
      scale_fill_identity() +
      coord_cartesian(xlim = c(0.0, 1.65),
                      ylim = c(0, Y_SPAN),
                      expand = FALSE, clip = "off") +
      theme_void() +
      theme(
        axis.ticks.length = unit(0, "pt"),
        plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 1, l = 0)
      )
  }

  # ===========================================================================
  # Sub-panel C3: Enrichment Bars
  # ===========================================================================

  if (nrow(cl_enrich) == 0) {
    p_bars <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nenrichment",
               size = 2.5, color = "grey50") +
      theme_void() +
      theme(
        axis.ticks.length = unit(0, "pt"),
        plot.margin = margin(t = 1, r = 1, b = if (show_xlab) 4 else 1, l = 0)
      )
  } else {

    bar_data <- cl_enrich %>%
      mutate(
        neg_log10_p = -log10(p.adjust),
        clean_name  = clean_pathway_name(Description),
        gene_count  = as.numeric(sub("/.*", "", GeneRatio)),
        db_fill     = DB_COLORS[database]
      )

    # Order pathways by -log10(p.adjust) within cluster
    bar_data <- bar_data %>%
      mutate(clean_name = reorder_within(clean_name, neg_log10_p, cluster))

    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3) +
      geom_text(aes(label = gene_count),
                hjust = -0.2, size = KEY_TEXT, fontface = "bold") +
      # Break indicator for capped bars
      {if (any(bar_data$neg_log10_p > shared_x_max))
        geom_text(
          data = bar_data %>% filter(neg_log10_p > shared_x_max),
          aes(x = shared_x_max, y = clean_name),
          label = "//", size = 2.0, fontface = "bold", color = "grey40",
          hjust = 0.5
        )
      } +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      scale_fill_identity() +
      scale_x_continuous(
        limits = c(0, shared_x_max),
        oob    = scales::squish,
        expand = expansion(mult = c(0, 0.08)),
        breaks = seq(0, shared_x_max, by = 5),
        name   = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_reordered() +
      labs(y = NULL) +
      THEME_PUB +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = if (show_xlab) element_text(size = 5) else element_blank(),
        axis.ticks.x = if (show_xlab) element_line() else element_blank(),
        axis.title.x = if (show_xlab) element_text(size = 6) else element_blank(),
        axis.ticks.length = unit(0, "pt"),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
        plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 1, l = 0)
      )
  }

  # ===========================================================================
  # Assemble the triptych
  # ===========================================================================

  (p_ht | p_sankey | p_bars) + plot_layout(widths = c(0.18, 0.38, 0.44))
}

# --- Build all clusters ---
panels_C <- lapply(1:optimal_k, function(i) {
  build_panel_C(i, show_xlab = (i == optimal_k), shared_x_max = shared_x_max_C)
})

cat("  Panel C: built", length(panels_C), "triptych rows\n")

# --- Save individual panel PDF ---
col_C <- wrap_plots(panels_C, ncol = 1, heights = row_heights)
ggsave(file.path(RPT_DIR, "panel_C_triptych.pdf"), col_C,
       width = FIG_W * COL_WIDTHS[3], height = FIG_H * 0.90,
       units = "mm", device = pdf)

cat("  Panel C complete\n")
