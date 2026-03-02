################################################################################
#   Figure 4 — Panel B: PCA highlight-on-grey scatter
#
#   Requires from setup: core_proteins, delta_mat, CLUSTER_COLORS, optimal_k,
#                         THEME_PUB, RPT_DIR, DAT_DIR, FIG_W, FIG_H,
#                         COL_WIDTHS, row_heights, cluster_ids, n_clusters
#   Outputs: panels_B (list of ggplot objects), pca_scores
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. PCA is descriptive (dimension reduction for visualization).
#    No statistical test needed.                                        PASS
#
# 2. Variance explained:
#    - PC1 and PC2 percentages reported in axis labels. These are point
#      estimates from eigenvalue decomposition. For ~2000 proteins and
#      ~30 subjects, with N << p, the first few PCs capture subject-level
#      variance structure. No CI on variance explained is standard for
#      PCA scatter plots in proteomics. A Marchenko-Pastur null could be
#      added as a supplementary analysis but is not required for this
#      descriptive panel.                                               PASS
#
# 3. Centering/scaling:
#    - prcomp(center = TRUE, scale. = FALSE) on delta_mat (already
#      z-scored row-wise). Not double-scaling since the PCA is on
#      t(delta_mat) = subjects x proteins. Column centering (per-protein
#      mean removal across subjects) is appropriate.                    PASS
#
# 4. Protein loadings used for scatter:
#    - pca_result$rotation (protein loadings) used as coordinates, not
#      pca_result$x (subject scores). This shows how proteins separate
#      in the first two loading dimensions. Correct for showing cluster
#      geometry in protein space.                                       PASS
# ---------------------------------------------------------------------------

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building Panel B: PCA highlight-on-grey scatter ===\n")

# --- B1. Run PCA on the delta matrix ---
# delta_mat: proteins (rows) x subjects (cols)
# t(delta_mat) => subjects x proteins; prcomp returns:
#   $rotation  = protein loadings (proteins x PCs) — one point per protein
pca_result <- prcomp(t(delta_mat), center = TRUE, scale. = FALSE)

# Protein scores for plotting (each dot = one protein)
pca_coords <- as.data.frame(pca_result$rotation[, 1:2])
pca_coords$gene <- rownames(pca_coords)

# Variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])

cat(sprintf("  PCA: %d proteins projected, PC1=%.1f%%, PC2=%.1f%%\n",
            nrow(pca_coords), var_explained[1], var_explained[2]))

# --- B2. Create pca_scores data frame ---
pca_scores <- pca_coords %>%
  left_join(core_proteins %>% dplyr::select(gene, cluster, membership),
            by = "gene") %>%
  mutate(
    cluster    = ifelse(is.na(cluster), NA_character_, cluster),
    membership = ifelse(is.na(membership), 0, membership)
  )

cat(sprintf("  pca_scores: %d rows (%d with cluster assignment)\n",
            nrow(pca_scores), sum(!is.na(pca_scores$cluster))))

# --- B3. Export PCA scores ---
write_csv(pca_scores, file.path(DAT_DIR, "02_panel_B_pca_scores.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "02_panel_B_pca_scores.csv")))

# --- B4. Build one highlight-on-grey plot per cluster ---
panels_B <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]

  is_first <- (i == 1)
  is_last  <- (i == n_clusters)

  # Core proteins for this cluster
  highlight_data <- pca_scores %>%
    filter(cluster == cid)

  p <- ggplot() +
    # Background: ALL proteins as grey points
    geom_point(
      data = pca_scores,
      aes(x = PC1, y = PC2),
      color = "grey55", size = 0.3, alpha = 0.45
    ) +
    # Highlighted cluster: core proteins in cluster color
    geom_point(
      data = highlight_data,
      aes(x = PC1, y = PC2, alpha = membership),
      color = CLUSTER_COLORS[cid], size = 0.6
    ) +
    scale_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
    # Cluster label: top-left corner
    annotate("text",
             x = min(pca_scores$PC1), y = max(pca_scores$PC2),
             label = cid,
             color = CLUSTER_COLORS[cid], fontface = "bold", size = TXT_ANNOT,
             hjust = 0, vjust = 1) +
    # Left accent stripe
    annotate("segment",
             x = min(pca_scores$PC1), xend = min(pca_scores$PC1),
             y = min(pca_scores$PC2), yend = max(pca_scores$PC2),
             colour = CLUSTER_COLORS[cid], linewidth = 3,
             lineend = "butt") +
    # Axis labels: only on bottom row
    labs(
      x = if (is_last) pc1_lab else NULL,
      y = if (is_first) pc2_lab else NULL
    ) +
    THEME_PUB +
    theme(
      panel.border = element_rect(colour = "grey70",
                                  linewidth = 0.3, fill = NA),
      axis.title.x = if (is_last) element_text(size = TXT_AXIS, face = "bold") else element_blank(),
      axis.title.y = if (is_first) element_text(size = TXT_AXIS, face = "bold") else element_blank(),
      axis.text.x  = if (is_last) element_text(size = TXT_TICK, face = "bold") else element_blank(),
      axis.text.y  = element_text(size = TXT_TICK, face = "bold"),
      axis.ticks.x = if (is_last) element_line() else element_blank(),
      plot.margin  = margin(t = 2, r = 2, b = if (is_last) 4 else 1, l = -2)
    )

  p
})

cat("  Panel B: built", length(panels_B), "PCA highlight-on-grey plots\n")

# --- Save individual panel PDF ---
col_B <- wrap_plots(panels_B, ncol = 1, heights = row_heights)
ggsave(file.path(RPT_DIR, "panel_B_pca.pdf"), col_B,
       width = FIG_W * COL_WIDTHS[2] * 2.5, height = FIG_H * 0.90,
       units = "mm", device = pdf)

cat("  Panel B complete\n")
