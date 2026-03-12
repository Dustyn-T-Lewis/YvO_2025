################################################################################
#   Figure 4 — Panel A: Cluster profiles with subject spaghetti
#
#   Requires from setup: core_proteins, abund_mat, sample_meta, group_z,
#                         GROUP_COLS, GROUP_LABS, CLUSTER_COLORS, AGE_COLORS,
#                         optimal_k, THEME_PUB, RPT_DIR
#   Outputs: panels_A (list of ggplot objects), n_clusters, cluster_ids
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. What is plotted:
#    - Each protein's group-level z-score (mean across subjects within each
#      Age x Time group, then z-scored row-wise). The ribbon is mean +/- SE
#      across PROTEINS within the cluster, not across subjects.
#    - This is a descriptive visualization of the cluster centroid pattern,
#      not a formal statistical test.                                   PASS
#
# 2. SE ribbons:
#    - SE = sd / sqrt(n_proteins_in_cluster). This describes precision of
#      the cluster mean profile, which is appropriate for showing how
#      tightly the cluster members agree.                               PASS
#    - NOTE: The y-axis range uses 1.96 * SE for limits, giving approximate
#      95% CI bounds on the cluster mean. This is correct for the purpose
#      of setting axis limits.                                          PASS
#
# 3. Missing formal test:
#    - No Age x Time interaction test is performed per cluster. This is
#      intentional: the clusters are DEFINED by their delta pattern, so
#      testing whether groups differ in their response within a cluster
#      would be circular (clusters were identified from the same data).
#      A formal test would require an independent validation cohort.
#      The visualization correctly shows the pattern without overclaiming.
#                                                                       PASS
#
# 4. Effect sizes:
#    - Not applicable for a descriptive cluster profile visualization.  PASS
#
# 5. Sample size:
#    - Core proteins per cluster (membership >= 0.5) typically range from
#      ~200-800 per cluster. Large enough for stable mean estimation.   PASS
# ---------------------------------------------------------------------------

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building Panel A: cluster profiles with subject spaghetti ===\n")

# --- A1. Reshape to long format ---
panel_a_long <- as.data.frame(group_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(GROUP_COLS),
               names_to = "group",
               values_to = "z_score") %>%
  left_join(core_proteins %>% dplyr::select(gene, cluster, membership), by = "gene") %>%
  mutate(
    age      = ifelse(str_detect(group, "^Young"), "Young", "Old"),
    time     = ifelse(str_detect(group, "_Post$"), "Post", "Pre"),
    time_num = ifelse(time == "Pre", 1, 2)
  )

cat(sprintf("  Long format: %d rows\n", nrow(panel_a_long)))

# --- A2. Per-cluster summary: mean_z and se_z ---
panel_a_summary <- panel_a_long %>%
  group_by(cluster, age, time, time_num) %>%
  summarise(
    mean_z = mean(z_score, na.rm = TRUE),
    se_z   = sd(z_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- A3. Shared y-axis range across all clusters ---
y_range <- panel_a_summary %>%
  summarise(
    y_lo = min(mean_z - 1.96 * se_z, na.rm = TRUE),
    y_hi = max(mean_z + 1.96 * se_z, na.rm = TRUE)
  )
y_pad   <- (y_range$y_hi - y_range$y_lo) * 0.1
y_limits <- c(y_range$y_lo - y_pad, y_range$y_hi + y_pad)

cat(sprintf("  Shared y-axis range: [%.2f, %.2f]\n", y_limits[1], y_limits[2]))

# --- A4. Build one plot per cluster ---
cluster_ids <- paste0("C", seq_len(optimal_k))
n_clusters  <- length(cluster_ids)

panels_A <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]

  # Subset data for this cluster
  cl_data    <- panel_a_long %>% filter(cluster == cid)
  cl_summary <- panel_a_summary %>% filter(cluster == cid)

  n_total <- n_distinct(cl_data$gene)
  n_core  <- n_total  # already filtered to core_proteins

  # Biological label from top Hallmark enrichment
  bio_label <- top_hallmark$label[top_hallmark$cluster == cid]
  if (length(bio_label) == 0) bio_label <- ""

  # Determine if this is the first (top) or last (bottom) cluster
  is_first <- (i == 1)
  is_last  <- (i == n_clusters)

  # Build plot
  p <- ggplot() +
    # Subject spaghetti: ultra-thin lines
    geom_line(
      data = cl_data,
      aes(x = time_num, y = z_score, group = interaction(gene, age),
          colour = age, alpha = membership),
      linewidth = 0.15
    ) +
    # SE ribbon
    geom_ribbon(
      data = cl_summary,
      aes(x = time_num, ymin = mean_z - se_z, ymax = mean_z + se_z,
          fill = age, group = age),
      alpha = 0.15
    ) +
    # Centroid lines
    geom_line(
      data = cl_summary,
      aes(x = time_num, y = mean_z, colour = age, group = age),
      linewidth = 1.2
    ) +
    # Centroid points
    geom_point(
      data = cl_summary,
      aes(x = time_num, y = mean_z, colour = age),
      size = 2.5
    ) +
    # Left accent stripe (cluster identity)
    annotate("segment",
             x = 0.7, xend = 0.7,
             y = y_limits[1], yend = y_limits[2],
             colour = CLUSTER_COLORS[cid], linewidth = 3,
             lineend = "butt") +
    # Scales
    scale_colour_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    scale_alpha_continuous(range = c(0.02, 0.15), guide = "none") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = if (is_last) c("Pre", "Post") else NULL,
      limits = c(0.7, 2.3),
      expand = expansion(0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      name   = NULL
    ) +
    # Labels
    labs(
      title    = sprintf("C%d: %s", i, bio_label),
      subtitle = sprintf("(n = %d)", n_core),
      x        = if (is_last) "Time" else NULL
    ) +
    # Theme
    THEME_PUB +
    theme(
      plot.title       = element_text(colour = CLUSTER_COLORS[cid],
                                      face = "bold", size = TXT_TITLE, hjust = 0.5),
      plot.subtitle    = element_text(colour = "grey30", face = "bold.italic",
                                      size = TXT_SUBTITLE, hjust = 0.5),
      panel.border     = element_rect(colour = "grey70",
                                      linewidth = 0.3, fill = NA),
      axis.title.y     = element_blank(),
      axis.title.x     = if (is_last) element_text(size = TXT_AXIS, face = "bold") else element_blank(),
      axis.text.y      = element_text(size = TXT_TICK, face = "bold"),
      axis.text.x      = if (is_last) element_text(size = TXT_TICK, face = "bold") else element_blank(),
      axis.ticks.x     = if (is_last) element_line() else element_blank(),
      plot.margin      = margin(t = 2, r = -2, b = if (is_last) 4 else 1, l = 2)
    )

  p
})

cat("  Panel A: built", length(panels_A), "cluster profile plots\n")

# --- Save individual panel PDF ---
col_A_inner <- wrap_plots(panels_A, ncol = 1, heights = row_heights)
ggsave(file.path(RPT_DIR, "panel_A_profiles.pdf"), col_A_inner,
       width = FIG_W * COL_WIDTHS[1] * 2.5, height = FIG_H * 0.90,
       units = "mm", device = pdf)

cat("  Panel A complete\n")
