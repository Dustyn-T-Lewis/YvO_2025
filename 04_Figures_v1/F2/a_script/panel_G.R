################################################################################
#   Figure 2 — Panel G: Per-Contrast ORA Bar Chart (4 Databases)
#
#   Runs over-representation analysis (hypergeometric test, BH correction)
#   on each F2 significance class gene set using 4 databases:
#   Hallmark, GO:BP, GO:MF, GO:CC. Gene symbols via enricher() throughout.
#
#   Top 3 terms per database per contrast, displayed as horizontal bars
#   faceted by contrast with database indicated by label prefix.
#
#   Gene lists (direction-collapsed, mutually exclusive from expanded classes):
#     Interaction: all "Interaction"
#     Sig Both:    "Sig Both Up" + "Sig Both Down"
#     Sig Young:   "Sig Young Up" + "Sig Young Down"
#     Sig Old:     "Sig Old Up" + "Sig Old Down"
#
#   Generates:
#     b_reports/panel_G_ora.pdf, panel_G_ora.png
#     c_data/panel_G/ora_per_contrast.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_F2_setup.R")

message("Panel G: per-contrast ORA bar chart (4 databases)...")

dir.create(file.path(DAT_DIR, "panel_G"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. BUILD GENE LISTS (direction-collapsed from expanded classification)
# ==============================================================================

if (!exists("expanded_df")) {
  expanded_df <- dep_df %>%
    filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
    mutate(
      expanded_cat = classify_expanded(
        pi_score_Training_Young, pi_score_Training_Old, pi_score_Interaction,
        logFC_Training_Young, logFC_Training_Old
      )
    ) %>%
    filter(!is.na(expanded_cat))
}

gene_lists <- list(
  "Interaction" = expanded_df %>%
    filter(expanded_cat == "Interaction") %>% pull(gene),
  "Sig Both"    = expanded_df %>%
    filter(expanded_cat %in% c("Sig Both Up", "Sig Both Down")) %>% pull(gene),
  "Sig Young"   = expanded_df %>%
    filter(expanded_cat %in% c("Sig Young Up", "Sig Young Down")) %>% pull(gene),
  "Sig Old"     = expanded_df %>%
    filter(expanded_cat %in% c("Sig Old Up", "Sig Old Down")) %>% pull(gene)
)

gene_lists <- gene_lists[sapply(gene_lists, length) > 0]

for (nm in names(gene_lists))
  message(sprintf("  %s: %d genes", nm, length(gene_lists[[nm]])))

# ==============================================================================
# 2. DATABASES: Hallmark + GO:BP + GO:MF + GO:CC (gene symbols via msigdbr)
# ==============================================================================

db_specs <- list(
  Hallmark = list(collection = "H",  subcollection = NA),
  `GO:BP`  = list(collection = "C5", subcollection = "GO:BP"),
  `GO:MF`  = list(collection = "C5", subcollection = "GO:MF"),
  `GO:CC`  = list(collection = "C5", subcollection = "GO:CC")
)

t2g_list <- lapply(names(db_specs), function(db_name) {
  spec <- db_specs[[db_name]]
  if (is.na(spec$subcollection)) {
    msigdbr(species = "Homo sapiens", collection = spec$collection)
  } else {
    msigdbr(species = "Homo sapiens", collection = spec$collection,
            subcollection = spec$subcollection)
  }
}) %>% setNames(names(db_specs))

# Convert to TERM2GENE format (term, gene)
t2g_list <- lapply(t2g_list, function(df) {
  df %>% select(term = gs_name, gene = gene_symbol) %>% distinct()
})

universe_symbols <- unique(dep_df$gene)
message(sprintf("  Universe: %d gene symbols", length(universe_symbols)))

# ==============================================================================
# 3. RUN ORA: 4 databases x N contrasts, top 3 per db per contrast
# ==============================================================================

TOP_N <- 3

run_db_ora <- function(gene_set, contrast_name, t2g, db_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    enricher(gene = gene_set, universe = universe_symbols, TERM2GENE = t2g,
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1,
             minGSSize = 10, maxGSSize = 500),
    error = function(e) NULL)
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    as.data.frame(res) %>%
      mutate(contrast = contrast_name, database = db_name,
             pathway_label = clean_pathway_name(ID)) %>%
      arrange(pvalue, nchar(ID)) %>%
      distinct(geneID, .keep_all = TRUE) %>%
      slice_head(n = TOP_N)
  } else {
    tibble()
  }
}

ora_results <- list()
for (contrast_name in names(gene_lists)) {
  genes <- gene_lists[[contrast_name]]
  for (db_name in names(t2g_list)) {
    ora_results[[length(ora_results) + 1]] <-
      run_db_ora(genes, contrast_name, t2g_list[[db_name]], db_name)
  }
}
ora_combined <- bind_rows(ora_results)

message(sprintf("  ORA complete: %d significant terms across %d contrast-database pairs",
                nrow(ora_combined),
                n_distinct(paste(ora_combined$contrast, ora_combined$database))))

# ==============================================================================
# 4. EXPORT CSV
# ==============================================================================

write_csv(ora_combined, file.path(DAT_DIR, "panel_G", "ora_per_contrast.csv"))
message("  Exported ora_per_contrast.csv")

# ==============================================================================
# 5. UNIFIED BAR CHART (no faceting, contrast shading bands, functional colors)
# ==============================================================================

if (nrow(ora_combined) > 0) {

  # Contrast band colors (F1-style shading)
  CONTRAST_BAND_COLORS <- c(
    "Interaction" = "#7B5EA7",
    "Sig Both"    = "#2E7D32",
    "Sig Young"   = "#E05A4E",
    "Sig Old"     = "#5DA5DA"
  )

  # Database display order and short prefix
  DB_PREFIX <- c(Hallmark = "[H]", `GO:BP` = "[BP]", `GO:MF` = "[MF]", `GO:CC` = "[CC]")

  contrast_order <- c("Interaction", "Sig Both", "Sig Young", "Sig Old")
  contrast_order <- contrast_order[contrast_order %in% unique(ora_combined$contrast)]
  db_order <- c("Hallmark", "GO:BP", "GO:MF", "GO:CC")

  # Gene counts per contrast for band labels
  n_genes_per_contrast <- setNames(
    sapply(gene_lists, length),
    names(gene_lists)
  )

  plot_df <- ora_combined %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      pathway_label = stringr::str_trunc(pathway_label, 38, ellipsis = "..."),
      bar_label = paste(DB_PREFIX[database], pathway_label),
      contrast = factor(contrast, levels = contrast_order),
      database = factor(database, levels = db_order)
    ) %>%
    group_by(contrast, database) %>%
    slice_head(n = TOP_N) %>%
    ungroup()

  # Functional classification via keyword classifier
  plot_df$func_category <- classify_pathway_func(plot_df$ID)
  plot_df$func_fill     <- CONSOLIDATED_COLORS[plot_df$func_category]
  plot_df$func_fill[is.na(plot_df$func_fill)] <- "#D0D0D0"

  # Order: group by contrast, then database, then significance (bottom to top)
  plot_df <- plot_df %>%
    arrange(contrast, database, neg_log10_padj) %>%
    mutate(
      pathway_uid = paste(bar_label, contrast, sep = "___"),
      pathway_uid = factor(pathway_uid, levels = unique(pathway_uid)),
      y_pos = as.numeric(pathway_uid)
    )

  # Compute contrast shading bands
  band_df <- plot_df %>%
    group_by(contrast) %>%
    summarise(
      ymin = min(y_pos) - 0.5,
      ymax = max(y_pos) + 0.5,
      y_mid = mean(y_pos),
      .groups = "drop"
    ) %>%
    mutate(
      n_genes = n_genes_per_contrast[as.character(contrast)],
      band_fill = CONTRAST_BAND_COLORS[as.character(contrast)]
    )

  # Compute database sub-bands within each contrast
  db_band_df <- plot_df %>%
    group_by(contrast, database) %>%
    summarise(
      ymin = min(y_pos) - 0.4,
      ymax = max(y_pos) + 0.4,
      y_mid = mean(y_pos),
      .groups = "drop"
    )

  pG <- ggplot(plot_df, aes(x = neg_log10_padj, y = pathway_uid)) +
    # Contrast shading bands
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = band_df$ymin, ymax = band_df$ymax,
             fill = band_df$band_fill,
             alpha = 0.15, color = "grey70", linewidth = 0.2) +
    # Database sub-band separators (dotted lines between database groups)
    annotate("segment",
             x = 0, xend = Inf,
             y = db_band_df$ymax,
             color = "grey60", linewidth = 0.15, linetype = "dotted") +
    # Bars colored by functional category
    geom_col(aes(fill = func_fill), color = "grey30", linewidth = 0.2, width = 0.6) +
    # Value labels at bar tips
    geom_text(aes(label = sprintf("%.1f", neg_log10_padj)),
              hjust = -0.15, size = TXT_ORA_BAR, fontface = "bold", color = "grey30") +
    # Contrast group labels inside bands (left side)
    annotate("text",
             x = -Inf, y = band_df$y_mid,
             label = sprintf("%s (n=%d)", band_df$contrast, band_df$n_genes),
             hjust = -0.05, fontface = "bold", size = TXT_ORA_BAR * 1.1,
             color = "grey20") +
    scale_fill_identity() +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.25))) +
    labs(
      x = expression(-log[10](p[adj])),
      y = NULL,
      title = "Per-Contrast Over-Representation Analysis",
      subtitle = sprintf("Top %d per database | Hallmark + GO:BP/MF/CC | BH correction | Bars colored by function",
                         TOP_N)
    ) +
    THEME_FIG +
    theme(
      axis.text.y  = element_text(size = rel(0.7)),
      axis.text.x  = element_text(size = rel(0.7)),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin = margin(2, 4, 2, 2, "mm")
    )

  # Dynamic height: uniform bars, no facet overhead
  n_bars <- nrow(plot_df)
  fig_h  <- max(120, n_bars * 10 + 60)

  ggsave(file.path(RPT_DIR, "panel_G_ora.pdf"), pG,
         width = 220, height = fig_h, units = "mm", device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(RPT_DIR, "panel_G_ora.png"), pG,
         width = 220, height = fig_h, units = "mm", dpi = 300, limitsize = FALSE)

  message(sprintf("  Panel G saved: %d terms across %d contrasts, 220 x %d mm",
                  n_bars, n_distinct(plot_df$contrast), fig_h))
} else {
  message("  Panel G skipped (no significant ORA results)")
  pG <- ggplot() + theme_void() + labs(title = "Panel G: No significant ORA results")
}

message("  Panel G complete")
