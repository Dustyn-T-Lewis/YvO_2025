################################################################################
#   Figure 4 — Panel D: Cluster synthesis Sankey + stacked bars
#
#   Requires from setup: core_proteins, protein_pathway_links, theme_links,
#                         CLUSTER_COLORS, THEME_COLORS, optimal_k, row_heights,
#                         make_sigmoid_ribbon, assign_theme, DAT_DIR, RPT_DIR,
#                         FIG_W, FIG_H, COL_WIDTHS
#   Outputs: panel_D (ggplot object)
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Theme assignment:
#    - Rule-based regex mapping of pathway names to biological themes
#      (assign_theme() in setup). This is a synthesis/summary heuristic
#      for narrative clarity, not a statistical test. No statistical
#      validation needed.                                               PASS
#
# 2. Stacked bars (protein counts per theme):
#    - Purely descriptive: counts of proteins assigned to each theme,
#      segmented by cluster. No statistical test (e.g., chi-squared for
#      uneven distribution) is needed because the theme categories are
#      derived from the data (via pathway enrichment), so testing their
#      uniformity would be circular.                                    PASS
#
# 3. Sankey ribbons:
#    - Proportional flow visualization from clusters to themes. Width
#      proportional to protein count. Descriptive.                      PASS
#
# 4. "Other" category exclusion:
#    - Proteins assigned to "Other" theme (no regex match) are excluded
#      from the Sankey but retained in exported CSV. This is a display
#      decision that does not affect statistical conclusions.           PASS
# ---------------------------------------------------------------------------

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building Panel D: Cluster synthesis Sankey ===\n")

# --- Step 1: Theme curation ---------------------------------------------------

cat(sprintf("  Theme assignment: %d mapped proteins across themes\n", nrow(theme_links)))
cat("  Theme distribution:\n")
print(table(theme_links$theme))

# Compute summaries
cluster_theme_counts <- theme_links %>%
  count(cluster, theme, name = "n_proteins") %>%
  arrange(cluster, desc(n_proteins))

theme_totals <- cluster_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))

theme_order <- theme_totals$theme

cat("  Theme totals (ordered):\n")
for (i in seq_len(nrow(theme_totals))) {
  cat(sprintf("    %s: %d proteins\n", theme_totals$theme[i], theme_totals$total[i]))
}

# Export
write_csv(cluster_theme_counts, file.path(DAT_DIR, "05_panel_D_cluster_themes.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "05_panel_D_cluster_themes.csv")))

# --- Step 2: Cluster-to-theme Sankey coordinate system ------------------------

# Exclude "Other" from the Sankey diagram — it adds no biological insight
# and compresses the meaningful themes.  Keep it in the exported CSV.
sankey_theme_counts <- cluster_theme_counts %>% filter(theme != "Other")
sankey_theme_totals <- sankey_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))
sankey_theme_order <- sankey_theme_totals$theme

cat(sprintf("  Sankey themes (excl. Other): %d themes, %d proteins\n",
            length(sankey_theme_order), sum(sankey_theme_totals$total)))

D_Y_SPAN <- 100
D_X_CL   <- 1.0
D_X_TH   <- 4.0
D_BAR_W  <- 0.20   # wider bars so cluster labels are legible

# Cluster mapped protein counts (only non-Other themed proteins)
cluster_mapped <- sankey_theme_counts %>%
  group_by(cluster) %>%
  summarise(n_mapped = sum(n_proteins), .groups = "drop") %>%
  arrange(cluster)

# Ensure all active clusters are present
active_clusters <- sort(unique(sankey_theme_counts$cluster))
n_active <- length(active_clusters)

cat(sprintf("  Active clusters with mapped proteins: %d\n", n_active))

# --- Cluster bars (left side) ---
# Center each cluster bar on its corresponding Panel C row (proportional to core protein count)
# row_heights are proportional to core protein count per cluster (computed in setup)
row_cum <- cumsum(row_heights)
row_tops <- D_Y_SPAN - c(0, head(row_cum, -1)) * D_Y_SPAN
row_bots <- D_Y_SPAN - row_cum * D_Y_SPAN
# Nudge centers downward to compensate for patchwork inter-panel gaps in A/B/C
pw_gap_nudge <- c(0, -0.015, -0.03, -0.045) * D_Y_SPAN
row_ctrs <- (row_tops + row_bots) / 2 + pw_gap_nudge

cl_bars <- cluster_mapped %>%
  mutate(
    cl_idx  = match(cluster, paste0("C", 1:optimal_k)),
    row_h   = (row_tops - row_bots)[cl_idx],
    row_ctr = row_ctrs[cl_idx],
    bar_h   = (n_mapped / max(n_mapped)) * row_h * 0.85,
    y_ctr   = row_ctr,
    y_top   = y_ctr + bar_h / 2,
    y_bot   = y_ctr - bar_h / 2,
    x_left  = D_X_CL - D_BAR_W / 2,
    x_right = D_X_CL + D_BAR_W / 2,
    fill    = CLUSTER_COLORS[cluster]
  )

# --- Theme bars (right side) ---
# One bar per theme (excl. Other), ordered by total protein count (descending)
n_themes <- length(sankey_theme_order)
th_gap_frac <- 0.03
th_usable   <- D_Y_SPAN * (1 - th_gap_frac * (n_themes - 1) / n_themes)
th_gap_size <- if (n_themes > 1) (D_Y_SPAN - th_usable) / (n_themes - 1) else 0
th_total    <- sum(sankey_theme_totals$total)

th_bars <- sankey_theme_totals %>%
  mutate(
    theme = factor(theme, levels = sankey_theme_order),
    bar_h   = (total / th_total) * th_usable,
    y_top   = D_Y_SPAN - cumsum(c(0, head(bar_h + th_gap_size, -1))),
    y_bot   = y_top - bar_h,
    y_ctr   = (y_top + y_bot) / 2,
    x_left  = D_X_TH - D_BAR_W / 2,
    x_right = D_X_TH + D_BAR_W / 2,
    fill    = THEME_COLORS[as.character(theme)]
  ) %>%
  arrange(theme)

# --- Ribbon construction: cluster side tracking ---
# For each cluster, track cumulative allocation of its bar space across themes
# For each theme, track cumulative allocation of its bar space across clusters

# Initialize cumulative trackers
cl_cum <- setNames(rep(0, n_active), active_clusters)   # fraction used per cluster
th_cum <- setNames(rep(0, n_themes), sankey_theme_order) # fraction used per theme

# Build ribbons by iterating themes (top to bottom) then clusters within each theme
ribbon_list <- list()
ribbon_idx  <- 0

for (th in sankey_theme_order) {
  # Get cluster contributions to this theme
  th_contribs <- sankey_theme_counts %>%
    filter(theme == th) %>%
    arrange(cluster)

  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next

    ribbon_idx <- ribbon_idx + 1

    # --- Cluster bar side: slice proportional to this cluster's mapped count ---
    cl_row   <- cl_bars %>% filter(cluster == cl)
    cl_n     <- cl_row$n_mapped
    cl_h     <- cl_row$bar_h
    cl_y_top <- cl_row$y_top

    # This ribbon occupies fraction n/cl_n of the cluster bar
    frac_cl  <- n / cl_n
    # Start from cumulative position (top of bar going down)
    y0_top <- cl_y_top - cl_cum[cl] * cl_h
    y0_bot <- y0_top - frac_cl * cl_h
    cl_cum[cl] <- cl_cum[cl] + frac_cl

    # --- Theme bar side: slice proportional to theme total ---
    th_row   <- th_bars %>% filter(theme == th)
    th_n     <- th_row$total
    th_h     <- th_row$bar_h
    th_y_top <- th_row$y_top

    frac_th  <- n / th_n
    y1_top <- th_y_top - th_cum[th] * th_h
    y1_bot <- y1_top - frac_th * th_h
    th_cum[th] <- th_cum[th] + frac_th

    # Build sigmoid ribbon polygon
    ribbon_poly <- make_sigmoid_ribbon(
      x0 = D_X_CL + D_BAR_W / 2,
      x1 = D_X_TH - D_BAR_W / 2,
      y0_top = y0_top,
      y0_bot = y0_bot,
      y1_top = y1_top,
      y1_bot = y1_bot,
      ribbon_id = paste0("D_", cl, "_", th)
    ) %>%
      mutate(
        cluster = cl,
        theme   = th,
        fill    = CLUSTER_COLORS[cl]
      )

    ribbon_list[[ribbon_idx]] <- ribbon_poly
  }
}

D_ribbons <- bind_rows(ribbon_list)
cat(sprintf("  Built %d ribbons for Panel D Sankey\n", ribbon_idx))

# --- Step 3: Build stacked bar geometry in the same coordinate system ---------
# Layout zones (x-coordinates):
#   D_X_CL (1.0)              — cluster bars
#   D_X_TH (4.0)              — theme bars (with labels to the left)
#   D_X_BAR_START (4.25)      — stacked bars begin (snug after theme bars)
#   D_X_BAR_END (~5.75)       — stacked bars end + total labels

D_X_BAR_START <- D_X_TH + D_BAR_W / 2 + 0.15
D_MAX_BAR_LEN <- 1.5
max_theme_count <- max(sankey_theme_totals$total)

# Stacked bar height: use a fixed bar height centered on each theme bar center
D_SBAR_H <- 2.2   # fixed height for all stacked bars

# Build stacked bar rectangles: for each theme, one rect per cluster contribution
stacked_rects <- list()
stacked_idx <- 0

for (th in sankey_theme_order) {
  th_row <- th_bars %>% filter(theme == th)
  th_contribs <- sankey_theme_counts %>%
    filter(theme == th) %>%
    arrange(cluster)

  # Center stacked bar on theme bar center
  bar_y_ctr <- th_row$y_ctr
  bar_y_top <- bar_y_ctr + D_SBAR_H / 2
  bar_y_bot <- bar_y_ctr - D_SBAR_H / 2

  # Stack cluster segments left to right
  x_cursor <- D_X_BAR_START

  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next

    stacked_idx <- stacked_idx + 1
    seg_w <- (n / max_theme_count) * D_MAX_BAR_LEN

    stacked_rects[[stacked_idx]] <- tibble(
      xmin = x_cursor,
      xmax = x_cursor + seg_w,
      ymin = bar_y_bot,
      ymax = bar_y_top,
      cluster = cl,
      theme   = th,
      fill    = CLUSTER_COLORS[cl]
    )

    x_cursor <- x_cursor + seg_w
  }
}

D_stacked <- bind_rows(stacked_rects)

# --- Reference gridlines for stacked bars ---
grid_intervals <- seq(50, max_theme_count, by = 50)
D_grid_x <- D_X_BAR_START + (grid_intervals / max_theme_count) * D_MAX_BAR_LEN
D_grid_df <- tibble(x = D_grid_x)

# Total labels at bar tips
D_bar_totals <- sankey_theme_totals %>%
  mutate(
    x_tip = D_X_BAR_START + (total / max_theme_count) * D_MAX_BAR_LEN
  ) %>%
  left_join(th_bars %>% dplyr::select(theme, y_ctr), by = "theme")

cat(sprintf("  Built %d stacked bar segments\n", stacked_idx))

# --- Build the combined Sankey + stacked bars plot ---
p_D_sankey <- ggplot() +
  # Reference gridlines (behind everything)
  geom_segment(
    data = D_grid_df,
    aes(x = x, xend = x),
    y = min(th_bars$y_bot) - D_SBAR_H,
    yend = max(th_bars$y_top) + D_SBAR_H,
    color = "grey85", linewidth = 0.2, linetype = "dotted"
  ) +
  # Ribbons (drawn first = behind bars)
  geom_polygon(
    data = D_ribbons,
    aes(x = x, y = y, group = ribbon_id, fill = fill),
    alpha = 0.30, color = NA
  ) +
  # Cluster bars (left)
  geom_rect(
    data = cl_bars,
    aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
    fill = cl_bars$fill, color = "black", linewidth = 0.3
  ) +
  # Cluster labels (inside bars, white text)
  geom_text(
    data = cl_bars,
    aes(x = (x_left + x_right) / 2, y = y_ctr, label = cluster),
    color = "white", fontface = "bold", size = 2.8
  ) +
  # Theme bars (right side of Sankey, colored)
  geom_rect(
    data = th_bars,
    aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
    fill = th_bars$fill, color = "black", linewidth = 0.3
  ) +
  # Theme labels (LEFT of colored bars)
  geom_text(
    data = th_bars,
    aes(x = x_left - 0.08, y = y_ctr, label = theme),
    hjust = 1, size = 2.5, fontface = "bold",
    color = "grey20"
  ) +
  # Stacked bars (aligned to theme bar y-centers, fixed height)
  geom_rect(
    data = D_stacked,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = D_stacked$fill, color = "black", linewidth = 0.2
  ) +
  # Total count at bar tips
  geom_text(
    data = D_bar_totals,
    aes(x = x_tip + 0.06, y = y_ctr, label = total),
    hjust = 0, size = KEY_TEXT, fontface = "bold", color = "grey30"
  ) +
  # X-axis label for stacked bars
  annotate("text",
           x = D_X_BAR_START + D_MAX_BAR_LEN / 2,
           y = -3,
           label = "Protein count", size = 2.2, color = "grey40") +
  scale_fill_identity() +
  coord_cartesian(
    xlim = c(0.5, D_X_BAR_START + D_MAX_BAR_LEN + 0.7),
    ylim = c(-5, D_Y_SPAN + 2),
    clip = "off",
    expand = FALSE
  ) +
  theme_void() +
  theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 0))

cat("  Panel D combined plot built\n")

# --- Step 4: Panel D is now a single unified plot ---
panel_D <- p_D_sankey

# --- Save individual panel PDF ---
ggsave(file.path(RPT_DIR, "panel_D_synthesis.pdf"), panel_D,
       width = FIG_W * COL_WIDTHS[4], height = FIG_H * 0.90,
       units = "mm", device = pdf)

cat("  Panel D composed\n")

cat("  Panel D complete\n")
