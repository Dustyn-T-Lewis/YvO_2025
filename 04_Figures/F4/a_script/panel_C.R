# Figure 4 — Panel C: Per-cluster Sankey triptych (heatmap | Sankey | bars)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(scales)
})

RPT <- "04_Figures/F4/b_reports"
DAT <- "04_Figures/F4/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CORE_THRESH <- 0.5
core_proteins <- read_csv(file.path(DAT, "06_mfuzz_assignments.csv"),
                          show_col_types = FALSE) %>%
  filter(membership >= CORE_THRESH)

group_z <- readRDS(file.path(DAT, "group_z.rds"))

enrich_top <- read_csv(file.path(DAT, "03_panel_C_enrichment.csv"),
                       show_col_types = FALSE)

protein_pathway_links <- read_csv(file.path(DAT, "04_panel_C_sankey_links.csv"),
                                  show_col_types = FALSE)

GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
GROUP_LABS <- c("Young\nPre", "Young\nPost", "Old\nPre", "Old\nPost")
optimal_k  <- 4

row_heights <- core_proteins %>%
  count(cluster) %>%
  mutate(h = n / sum(n)) %>%
  pull(h)

pdf_device <- get_pdf_device()

message("Panel C: per-cluster Sankey triptych...")

PC_W <- 280  # panel width mm
PC_H <- 300  # stacked height mm
Z_CAP <- 2

shared_x_max_C <- 15

txt_annot <- scale_text(BASE_GENE, PC_W)
txt_tick  <- scale_text(BASE_GENE, PC_W)
txt_axis  <- scale_text(BASE_STAT, PC_W)

build_panel_C <- function(cl_id, show_xlab, shared_x_max) {

  cl_label <- paste0("C", cl_id)
  cl_col   <- CLUSTER_COLORS[cl_label]

  cl_core <- core_proteins %>%
    filter(cluster == cl_label) %>%
    arrange(desc(membership))
  cl_genes <- cl_core$gene
  n_genes <- length(cl_genes)

  cl_links <- protein_pathway_links %>% filter(cluster == cl_label)
  cl_enrich <- enrich_top %>% filter(cluster == cl_label) %>% arrange(p.adjust)

  mapped_links <- cl_links %>% filter(pathway != "Unmapped")
  mapped_pathways <- unique(mapped_links$pathway)
  n_pw <- length(mapped_pathways)

  # Heatmap strip
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
      colours = ht_colors, limits = c(-Z_CAP, Z_CAP),
      oob = scales::squish, guide = "none"
    ) +
    scale_x_discrete(
      labels = if (show_xlab) GROUP_LABS else NULL,
      expand = expansion(0)
    ) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.x     = if (show_xlab) element_text(size = txt_tick, face = "bold", lineheight = 0.85)
                           else element_blank(),
      axis.ticks.x    = if (show_xlab) element_line() else element_blank(),
      axis.ticks.length = unit(0, "pt"),
      panel.border    = element_rect(colour = cl_col, linewidth = 0.6, fill = NA),
      legend.position = "none",
      plot.margin     = margin(t = 1, r = 2, b = if (show_xlab) 4 else 1, l = 0)
    )

  # Sigmoid Sankey
  if (n_pw == 0) {
    p_sankey <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No mapped\npathways",
               size = txt_annot, fontface = "bold", color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 1, l = 0))
  } else {
    Y_SPAN  <- n_genes
    X_GENE  <- 0.05
    X_PW    <- 1.85
    BAR_W   <- 0.06

    gene_h   <- Y_SPAN / (n_genes * 1.15)
    gene_gap <- if (n_genes > 1) (Y_SPAN - n_genes * gene_h) / (n_genes - 1) else 0

    gene_bars <- tibble(
      gene = cl_genes, idx = seq_along(cl_genes),
      y_top = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
      y_bot = y_top - gene_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left = X_GENE - BAR_W / 2, x_right = X_GENE + BAR_W / 2,
      fill = cl_col
    )

    pw_h   <- Y_SPAN / (n_pw * 1.4)
    pw_gap <- if (n_pw > 1) (Y_SPAN - n_pw * pw_h) / (n_pw - 1) else 0

    pw_order <- cl_enrich %>%
      filter(Description %in% mapped_pathways) %>%
      pull(Description)
    pw_order <- c(pw_order, setdiff(mapped_pathways, pw_order))

    pw_bars <- tibble(
      pathway = pw_order, idx = seq_along(pw_order),
      y_top = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
      y_bot = y_top - pw_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left = X_PW - BAR_W / 2, x_right = X_PW + BAR_W / 2
    ) %>%
      left_join(mapped_links %>% dplyr::select(pathway, database) %>% distinct(),
                by = "pathway") %>%
      mutate(fill = DB_COLORS[database])

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

    ribbon_input <- slot_data %>%
      left_join(gene_bars %>% dplyr::select(gene, y_top_gene = y_top,
                                             y_bot_gene = y_bot, y_ctr_gene = y_ctr),
                by = "gene") %>%
      mutate(ribbon_id = paste0(cl_label, "_", gene, "_", pathway),
             cross_dist = abs(y_ctr_gene - slot_ctr)) %>%
      arrange(desc(cross_dist))

    ribbon_polys <- pmap_dfr(ribbon_input, function(gene, pathway, database, cluster,
                                                     slot_idx, n_slots,
                                                     y_top, y_bot,
                                                     slot_h, slot_top, slot_bot, slot_ctr,
                                                     y_top_gene, y_bot_gene, y_ctr_gene,
                                                     ribbon_id, cross_dist) {
      make_sigmoid_ribbon(
        x0 = X_GENE + BAR_W / 2, x1 = X_PW - BAR_W / 2,
        y0_top = y_top_gene, y0_bot = y_bot_gene,
        y1_top = slot_top, y1_bot = slot_bot,
        ribbon_id = ribbon_id
      ) %>% mutate(fill = DB_COLORS[database])
    })

    gene_rect <- gene_bars %>% dplyr::select(x_left, x_right, y_bot, y_top, fill)
    pw_rect   <- pw_bars %>% dplyr::select(x_left, x_right, y_bot, y_top, fill)

    p_sankey <- ggplot() +
      geom_polygon(data = ribbon_polys,
                   aes(x = x, y = y, group = ribbon_id, fill = fill),
                   alpha = 0.25, color = NA) +
      geom_rect(data = gene_rect,
                aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
                color = NA) +
      geom_rect(data = pw_rect,
                aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
                color = NA) +
      geom_text(data = pw_bars,
                aes(x = x_left - 0.04, y = y_ctr,
                    label = str_trunc(clean_pathway_name(pathway), 40, ellipsis = "...")),
                hjust = 1, size = txt_annot, fontface = "bold", color = "grey20") +
      scale_fill_identity() +
      coord_cartesian(xlim = c(0.0, 1.95), ylim = c(0, Y_SPAN),
                      expand = FALSE, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 1, l = 0))
  }

  # Enrichment bars
  if (nrow(cl_enrich) == 0) {
    p_bars <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nenrichment",
               size = txt_annot, fontface = "bold", color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 1, b = if (show_xlab) 4 else 1, l = 0))
  } else {
    bar_data <- cl_enrich %>%
      mutate(
        neg_log10_p = -log10(p.adjust),
        clean_name  = clean_pathway_name(Description),
        gene_count  = as.numeric(sub("/.*", "", GeneRatio)),
        db_fill     = DB_COLORS[database]
      ) %>%
      mutate(clean_name = reorder_within(clean_name, neg_log10_p, cluster))

    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3) +
      geom_text(aes(label = gene_count),
                hjust = -0.2, size = txt_annot, fontface = "bold") +
      {if (any(bar_data$neg_log10_p > shared_x_max))
        geom_text(data = bar_data %>% filter(neg_log10_p > shared_x_max),
                  aes(x = shared_x_max, y = clean_name),
                  label = "//", size = txt_annot, fontface = "bold", color = "grey40",
                  hjust = 0.5)
      } +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      scale_fill_identity() +
      scale_x_continuous(
        limits = c(0, shared_x_max), oob = scales::squish,
        expand = expansion(mult = c(0, 0.08)),
        breaks = seq(0, shared_x_max, by = 5),
        labels = function(x) ifelse(x %% 5 == 0, as.character(as.integer(x)), ""),
        name = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_reordered() +
      labs(y = NULL) +
      FIG_THEME +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = if (show_xlab) element_text(size = txt_tick, face = "bold") else element_blank(),
        axis.ticks.x = if (show_xlab) element_line() else element_blank(),
        axis.title.x = if (show_xlab) element_text(size = txt_axis, face = "bold") else element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
        legend.position = "none",
        plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 1, l = 0)
      )
  }

  (p_ht | p_sankey | p_bars) + plot_layout(widths = c(0.18, 0.38, 0.44))
}

panels_C <- lapply(1:optimal_k, function(i) {
  build_panel_C(i, show_xlab = (i == optimal_k), shared_x_max = shared_x_max_C)
})

col_C <- wrap_plots(panels_C, ncol = 1, heights = row_heights)
ggsave(file.path(RPT, "panel_C_triptych.pdf"), col_C,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_triptych.png"), col_C,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

message("  Panel C saved")
