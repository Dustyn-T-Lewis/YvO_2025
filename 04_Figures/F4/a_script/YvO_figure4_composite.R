################################################################################
#   Figure 4 — Composite Assembly + Supplementary Excel Workbook
#
#   Sources all per-panel scripts, assembles with patchwork into a single
#   Figure_4.pdf/png, and exports F4_supplementary.xlsx.
#
#   Panels:
#     A — FCM cluster profiles (spaghetti + ribbon)
#     B — PCA highlight-on-grey scatter
#     C — Per-cluster Sankey triptych (heatmap | Sankey | enrichment bars)
#     D — Cluster synthesis Sankey + stacked bars
#     Supp — Dmin elbow plot
################################################################################

# === 0. Source all panel scripts (setup loaded automatically by panel_A) =====

source("04_Figures/F4/a_script/panel_A.R")          # -> panels_A
source("04_Figures/F4/a_script/panel_B.R")          # -> panels_B, pca_scores
source("04_Figures/F4/a_script/panel_C.R")          # -> panels_C
source("04_Figures/F4/a_script/panel_D.R")          # -> panel_D, cluster_theme_counts
source("04_Figures/F4/a_script/supp_dmin_elbow.R")  # -> p_dmin

message("All F4 panel scripts sourced - assembling composite...")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    FIGURE ASSEMBLY (Legend + Headers + Panels)              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("=== Assembling final Figure 4 ===\n")

# ============================================================================
# PER-PANEL KEYS — separate legends under each panel column
# ============================================================================

cat("  Building per-panel keys...\n")

box_w <- 0.5
box_h <- 0.35
item_gap <- 0.6

# --- Key A/B: Age legend (Young / Old) ---
{
  ab_items <- list()
  x_cursor <- 0
  ab_items <- bind_rows(ab_items, tibble(
    x = x_cursor, y = 0, label = "Age:", type = "header", fill = NA_character_
  ))
  x_cursor <- x_cursor + 1.0
  for (nm in c("Young", "Old")) {
    ab_items <- bind_rows(ab_items, tibble(
      x = x_cursor, y = 0, label = NA_character_, type = "swatch",
      fill = AGE_COLORS[nm]
    ))
    ab_items <- bind_rows(ab_items, tibble(
      x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
      fill = NA_character_
    ))
    x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
  }
  sw_ab <- ab_items %>% filter(type == "swatch")
  hd_ab <- ab_items %>% filter(type == "header")
  it_ab <- ab_items %>% filter(type == "item")

  pKey_AB <- ggplot() +
    geom_rect(data = sw_ab,
              aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
              fill = sw_ab$fill, color = KEY_BORDER, linewidth = KEY_LW) +
    geom_text(data = hd_ab, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
    geom_text(data = it_ab, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
    coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))
}

# --- Key C: Database + Z-score ---
{
  c_items <- list()
  x_cursor <- 0
  c_items <- bind_rows(c_items, tibble(
    x = x_cursor, y = 0, label = "Database:", type = "header", fill = NA_character_
  ))
  x_cursor <- x_cursor + 1.8
  db_entries <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B", "GO:CC" = "#7E57C2")
  for (nm in names(db_entries)) {
    c_items <- bind_rows(c_items, tibble(
      x = x_cursor, y = 0, label = NA_character_, type = "swatch",
      fill = db_entries[nm]
    ))
    c_items <- bind_rows(c_items, tibble(
      x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
      fill = NA_character_
    ))
    x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
  }
  x_cursor <- x_cursor + 0.8
  c_items <- bind_rows(c_items, tibble(
    x = x_cursor, y = 0, label = "Z-score:", type = "header", fill = NA_character_
  ))
  x_cursor <- x_cursor + 1.5
  z_cols <- c("#2166AC", "white", "#B2182B")
  z_labs <- c("-2", "0", "+2")
  for (zi in seq_along(z_cols)) {
    c_items <- bind_rows(c_items, tibble(
      x = x_cursor, y = 0, label = NA_character_, type = "swatch",
      fill = z_cols[zi]
    ))
    c_items <- bind_rows(c_items, tibble(
      x = x_cursor + box_w + 0.05, y = 0, label = z_labs[zi], type = "item",
      fill = NA_character_
    ))
    x_cursor <- x_cursor + box_w + 0.05 + nchar(z_labs[zi]) * 0.22 + item_gap * 0.6
  }
  sw_c <- c_items %>% filter(type == "swatch")
  hd_c <- c_items %>% filter(type == "header")
  it_c <- c_items %>% filter(type == "item")

  pKey_C <- ggplot() +
    geom_rect(data = sw_c,
              aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
              fill = sw_c$fill, color = KEY_BORDER, linewidth = KEY_LW) +
    geom_text(data = hd_c, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
    geom_text(data = it_c, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
    coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))
}

# --- Key D: Cluster legend (C1-C4) ---
{
  d_items <- list()
  x_cursor <- 0
  d_items <- bind_rows(d_items, tibble(
    x = x_cursor, y = 0, label = "Cluster:", type = "header", fill = NA_character_
  ))
  x_cursor <- x_cursor + 1.6
  active_cls <- paste0("C", seq_len(optimal_k))
  for (cl in active_cls) {
    d_items <- bind_rows(d_items, tibble(
      x = x_cursor, y = 0, label = NA_character_, type = "swatch",
      fill = CLUSTER_COLORS[cl]
    ))
    d_items <- bind_rows(d_items, tibble(
      x = x_cursor + box_w + 0.1, y = 0, label = cl, type = "item",
      fill = NA_character_
    ))
    x_cursor <- x_cursor + box_w + 0.1 + nchar(cl) * 0.22 + item_gap
  }
  sw_d <- d_items %>% filter(type == "swatch")
  hd_d <- d_items %>% filter(type == "header")
  it_d <- d_items %>% filter(type == "item")

  pKey_D <- ggplot() +
    geom_rect(data = sw_d,
              aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
              fill = sw_d$fill, color = KEY_BORDER, linewidth = KEY_LW) +
    geom_text(data = hd_d, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
    geom_text(data = it_d, aes(x = x, y = 0, label = label),
              hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
    coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))
}

cat("  Per-panel keys built\n")

# ============================================================================
# STACK COLUMNS: proportional row heights per cluster
# ============================================================================

cat("  Stacking columns with proportional row heights...\n")
cat(sprintf("  Row height proportions: %s\n",
            paste(sprintf("%.3f", row_heights), collapse = ", ")))

col_A_inner <- wrap_plots(panels_A, ncol = 1, heights = row_heights)
y_label_grob <- wrap_elements(
  grid::textGrob("Z-score (group means)", rot = 90,
                 gp = grid::gpar(fontsize = 7, fontface = "bold"))
)
col_A_composed <- (y_label_grob | col_A_inner) + plot_layout(widths = c(0.08, 0.92))
col_A <- wrap_elements(full = col_A_composed)
col_B <- wrap_plots(panels_B, ncol = 1, heights = row_heights)
col_C <- wrap_plots(panels_C, ncol = 1, heights = row_heights)
col_D <- panel_D  # already full-height (single ggplot)

# ============================================================================
# HEADER ROW: panel labels
# ============================================================================

cat("  Building header row...\n")

make_header <- function(label) {
  ggplot() +
    labs(title = label) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0))
}

header_A <- make_header("A  Cluster Profiles")
header_B <- make_header("B  Cluster Geometry")
header_C <- make_header("C  Protein\u2013Pathway Mapping")
header_D <- make_header("D  Cluster Synthesis")

header_row <- (header_A | header_B | header_C | header_D) +
  plot_layout(widths = COL_WIDTHS)

# ============================================================================
# ASSEMBLE FULL FIGURE
# ============================================================================

cat("  Assembling full figure...\n")

body_row <- (col_A | col_B | col_C | col_D) +
  plot_layout(widths = COL_WIDTHS)

key_row <- (pKey_AB | plot_spacer() | pKey_C | pKey_D) +
  plot_layout(widths = COL_WIDTHS)

fig4 <- header_row / body_row / key_row +
  plot_layout(heights = c(0.04, 0.90, 0.06)) +
  plot_annotation(
    title = "Proteomic Response Archetypes to Resistance Training",
    subtitle = sprintf("Mfuzz FCM (k = %d, m = %.2f) on per-subject delta matrix (%d core proteins, membership >= %.1f) | Bootstrap ARI = %.3f (95%% CI [%.3f, %.3f])",
                       optimal_k, m_est, nrow(core_proteins), CORE_THRESH,
                       mean(boot_ari), quantile(boot_ari, 0.025), quantile(boot_ari, 0.975)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

cat("  Figure 4 assembled\n")

# ============================================================================
# SAVE FINAL OUTPUTS
# ============================================================================

cat("  Saving final outputs...\n")

ggsave(file.path(RPT_DIR, "Figure_4.pdf"), fig4,
       width = FIG_W, height = FIG_H, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "Figure_4.png"), fig4,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300)

cat(sprintf("Figure 4 saved to %s\n", RPT_DIR))

cat("=== Figure 4 composite assembly complete ===\n")

# === 2. Supplementary Excel workbook =========================================

library(openxlsx)

# Sheet 1: A_cluster_profiles — centroid data from panel A
sheet_A <- centroid_export

# Sheet 2: B_pca_scores — PCA scores from panel B
sheet_B <- pca_scores %>%
  dplyr::select(gene, PC1, PC2, cluster, membership)

# Sheet 3: C_enrichment — enrichment results from panel C
sheet_C <- enrich_top %>%
  dplyr::select(cluster, database, ID, Description, GeneRatio,
                pvalue, p.adjust, qvalue, Count)

# Sheet 4: C_sankey_links — protein-pathway links from panel C
sheet_C2 <- protein_pathway_links %>%
  dplyr::select(gene, pathway, database, cluster)

# Sheet 5: D_cluster_themes — cluster-theme counts from panel D
sheet_D <- cluster_theme_counts

# Sheet 6: _metadata — figure description
sheet_meta <- tibble(
  item = c(
    "figure_title",
    "description",
    "clustering_method",
    "parameters",
    "enrichment_method",
    "sample_sizes",
    "significance_threshold",
    "software",
    "date_generated"
  ),
  value = c(
    "Figure 4: Proteomic Response Archetypes to Resistance Training",
    "Mfuzz FCM clustering on per-subject delta matrix, with ORA enrichment (Hallmark, GO:BP, GO:CC) and rrvgo redundancy reduction.",
    "Fuzzy c-means (Mfuzz) with multi-start optimization (50 starts)",
    sprintf("k = %d, m = %.3f, core threshold = %.1f, bootstrap ARI = %.3f (95%% CI [%.3f, %.3f])",
            optimal_k, m_est, CORE_THRESH, mean(boot_ari),
            quantile(boot_ari, 0.025), quantile(boot_ari, 0.975)),
    "ORA (enricher, BH-adjusted p < 0.05) + rrvgo GO reduction (threshold 0.85)",
    sprintf("%d proteins, %d core proteins (membership >= %.1f), %d subjects",
            nrow(cluster_assign), nrow(core_proteins), CORE_THRESH,
            ncol(delta_mat)),
    "p.adjust < 0.05 (enrichment), membership >= 0.5 (core proteins)",
    "R, Mfuzz, clusterProfiler, rrvgo, patchwork",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

wb <- createWorkbook()
addWorksheet(wb, "A_cluster_profiles");  writeData(wb, "A_cluster_profiles",  sheet_A)
addWorksheet(wb, "B_pca_scores");        writeData(wb, "B_pca_scores",        sheet_B)
addWorksheet(wb, "C_enrichment");        writeData(wb, "C_enrichment",        sheet_C)
addWorksheet(wb, "C_sankey_links");      writeData(wb, "C_sankey_links",      sheet_C2)
addWorksheet(wb, "D_cluster_themes");    writeData(wb, "D_cluster_themes",    sheet_D)
addWorksheet(wb, "_metadata");           writeData(wb, "_metadata",           sheet_meta)

# Bold headers on every sheet
hs <- createStyle(textDecoration = "bold")
for (sn in names(wb)) {
  addStyle(wb, sn, hs, rows = 1, cols = seq_len(20), gridExpand = TRUE)
  setColWidths(wb, sn, cols = seq_len(20), widths = "auto")
}

xlsx_path <- file.path(DAT_DIR, "F4_supplementary.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved %s (%d sheets)\n", xlsx_path, length(names(wb))))

cat("Figure 4 composite assembly complete.\n")
