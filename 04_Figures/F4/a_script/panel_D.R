# Figure 4 — Panel D: Cluster synthesis Sankey + stacked bars

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(patchwork)
})

RPT <- "04_Figures/F4/b_reports"
DAT <- "04_Figures/F4/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CORE_THRESH <- 0.5
core_proteins <- read_csv(file.path(DAT, "06_mfuzz_assignments.csv"),
                          show_col_types = FALSE) %>%
  filter(membership >= CORE_THRESH)

protein_pathway_links <- read_csv(file.path(DAT, "04_panel_C_sankey_links.csv"),
                                  show_col_types = FALSE)

theme_links <- read_csv(file.path(DAT, "theme_links.csv"),
                        show_col_types = FALSE)

optimal_k <- 4

row_heights <- core_proteins %>%
  count(cluster) %>%
  mutate(h = n / sum(n)) %>%
  pull(h)

n_other_theme  <- sum(theme_links$theme == "Other")
n_unmapped_all <- sum(protein_pathway_links$pathway == "Unmapped")
n_excluded     <- n_other_theme + n_unmapped_all
pct_excluded   <- 100 * n_excluded / nrow(core_proteins)

pdf_device <- get_pdf_device()

message("Panel D: cluster synthesis Sankey + stacked bars...")

PD_W <- 220  # panel width mm
PD_H <- 300  # panel height mm

txt_annot <- scale_text(BASE_GENE, PD_W)

cluster_theme_counts <- theme_links %>%
  count(cluster, theme, name = "n_proteins") %>%
  arrange(cluster, desc(n_proteins))

theme_totals <- cluster_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))

write_csv(cluster_theme_counts, file.path(DAT, "05_panel_D_cluster_themes.csv"))

sankey_theme_counts <- cluster_theme_counts %>% filter(theme != "Other")
sankey_theme_totals <- sankey_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))
sankey_theme_order <- sankey_theme_totals$theme

D_Y_SPAN <- 100
D_X_CL   <- 0.4
D_X_TH   <- 3.4
D_BAR_W  <- 0.20

cluster_mapped <- sankey_theme_counts %>%
  group_by(cluster) %>%
  summarise(n_mapped = sum(n_proteins), .groups = "drop") %>%
  arrange(cluster)

active_clusters <- sort(unique(sankey_theme_counts$cluster))
n_active <- length(active_clusters)

row_cum <- cumsum(row_heights)
row_tops <- D_Y_SPAN - c(0, head(row_cum, -1)) * D_Y_SPAN
row_bots <- D_Y_SPAN - row_cum * D_Y_SPAN
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

cl_cum <- setNames(rep(0, n_active), active_clusters)
th_cum <- setNames(rep(0, n_themes), sankey_theme_order)

ribbon_list <- list()
ribbon_idx  <- 0

for (th in sankey_theme_order) {
  th_contribs <- sankey_theme_counts %>% filter(theme == th) %>% arrange(cluster)
  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next
    ribbon_idx <- ribbon_idx + 1

    cl_row   <- cl_bars %>% filter(cluster == cl)
    cl_n     <- cl_row$n_mapped
    cl_h     <- cl_row$bar_h
    cl_y_top <- cl_row$y_top

    frac_cl  <- n / cl_n
    y0_top <- cl_y_top - cl_cum[cl] * cl_h
    y0_bot <- y0_top - frac_cl * cl_h
    cl_cum[cl] <- cl_cum[cl] + frac_cl

    th_row   <- th_bars %>% filter(theme == th)
    th_n     <- th_row$total
    th_h     <- th_row$bar_h
    th_y_top <- th_row$y_top

    frac_th  <- n / th_n
    y1_top <- th_y_top - th_cum[th] * th_h
    y1_bot <- y1_top - frac_th * th_h
    th_cum[th] <- th_cum[th] + frac_th

    ribbon_list[[ribbon_idx]] <- make_sigmoid_ribbon(
      x0 = D_X_CL + D_BAR_W / 2, x1 = D_X_TH - D_BAR_W / 2,
      y0_top = y0_top, y0_bot = y0_bot,
      y1_top = y1_top, y1_bot = y1_bot,
      ribbon_id = paste0("D_", cl, "_", th)
    ) %>% mutate(cluster = cl, theme = th, fill = CLUSTER_COLORS[cl])
  }
}
D_ribbons <- bind_rows(ribbon_list)

D_X_BAR_START <- D_X_TH + D_BAR_W / 2 + 0.15
D_MAX_BAR_LEN <- 1.5
max_theme_count <- max(sankey_theme_totals$total)
D_SBAR_H <- 2.2

stacked_rects <- list()
stacked_idx <- 0

for (th in sankey_theme_order) {
  th_row <- th_bars %>% filter(theme == th)
  th_contribs <- sankey_theme_counts %>% filter(theme == th) %>% arrange(cluster)
  bar_y_ctr <- th_row$y_ctr
  bar_y_top <- bar_y_ctr + D_SBAR_H / 2
  bar_y_bot <- bar_y_ctr - D_SBAR_H / 2
  x_cursor <- D_X_BAR_START

  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next
    stacked_idx <- stacked_idx + 1
    seg_w <- (n / max_theme_count) * D_MAX_BAR_LEN
    stacked_rects[[stacked_idx]] <- tibble(
      xmin = x_cursor, xmax = x_cursor + seg_w,
      ymin = bar_y_bot, ymax = bar_y_top,
      cluster = cl, theme = th, fill = CLUSTER_COLORS[cl]
    )
    x_cursor <- x_cursor + seg_w
  }
}
D_stacked <- bind_rows(stacked_rects)

grid_intervals <- seq(50, max_theme_count, by = 50)
D_grid_x <- D_X_BAR_START + (grid_intervals / max_theme_count) * D_MAX_BAR_LEN
D_grid_df <- tibble(x = D_grid_x)

D_bar_totals <- sankey_theme_totals %>%
  mutate(x_tip = D_X_BAR_START + (total / max_theme_count) * D_MAX_BAR_LEN) %>%
  left_join(th_bars %>% dplyr::select(theme, y_ctr), by = "theme")

panel_D <- ggplot() +
  geom_segment(data = D_grid_df,
               aes(x = x, xend = x),
               y = min(th_bars$y_bot) - D_SBAR_H,
               yend = max(th_bars$y_top) + D_SBAR_H,
               color = "grey85", linewidth = 0.2, linetype = "dotted") +
  geom_polygon(data = D_ribbons,
               aes(x = x, y = y, group = ribbon_id, fill = fill),
               alpha = 0.30, color = NA) +
  geom_rect(data = cl_bars,
            aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
            fill = cl_bars$fill, color = "black", linewidth = 0.3) +
  geom_text(data = cl_bars,
            aes(x = (x_left + x_right) / 2, y = y_ctr, label = cluster),
            color = "white", fontface = "bold", size = txt_annot) +
  geom_rect(data = th_bars,
            aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
            fill = th_bars$fill, color = "black", linewidth = 0.3) +
  geom_text(data = th_bars,
            aes(x = x_left - 0.08, y = y_ctr, label = theme),
            hjust = 1, size = txt_annot, fontface = "bold", color = "grey20") +
  geom_rect(data = D_stacked,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = D_stacked$fill, color = "black", linewidth = 0.2) +
  geom_text(data = D_bar_totals,
            aes(x = x_tip + 0.06, y = y_ctr, label = total),
            hjust = 0, size = txt_annot, fontface = "bold", color = "grey30") +
  annotate("text",
           x = D_X_BAR_START + D_MAX_BAR_LEN / 2, y = -3,
           label = "Protein count", size = txt_annot, fontface = "bold", color = "grey40") +
  annotate("text",
           x = D_X_BAR_START + D_MAX_BAR_LEN / 2, y = -6,
           label = sprintf("%d proteins (%.0f%%) not assigned to a theme",
                           n_excluded, pct_excluded),
           size = txt_annot * 0.85, fontface = "italic", color = "grey50") +
  scale_fill_identity() +
  coord_cartesian(
    xlim = c(0.05, D_X_BAR_START + D_MAX_BAR_LEN + 0.7),
    ylim = c(-5, D_Y_SPAN + 2),
    clip = "off", expand = FALSE
  ) +
  theme_void() +
  theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 0))

ggsave(file.path(RPT, "panel_D_synthesis.pdf"), panel_D,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_synthesis.png"), panel_D,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

message("  Panel D saved")
