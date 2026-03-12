################################################################################
#   Figure 4 — Panel A: Cluster profiles with subject spaghetti
#   Self-contained: loads pre-computed data from c_data/
#   Generates: panel_A_profiles.pdf/png
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
})

# ── Paths & constants ────────────────────────────────────────────────────────
RPT <- "04_Figures_v2/F4/b_reports"
DAT <- "04_Figures_v2/F4/c_data"
GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
optimal_k  <- 4

# ── Load data ────────────────────────────────────────────────────────────────
core_proteins <- read_csv(file.path(DAT, "06_mfuzz_assignments.csv"),
                          show_col_types = FALSE) %>%
  filter(membership >= 0.5)

group_z      <- readRDS(file.path(DAT, "group_z.rds"))
top_hallmark <- read_csv(file.path(DAT, "top_hallmark.csv"), show_col_types = FALSE)

# ── Row heights (proportional to cluster size) ───────────────────────────────
row_heights <- core_proteins %>%
  count(cluster) %>%
  arrange(cluster) %>%
  pull(n)
row_heights <- row_heights / sum(row_heights)

pdf_device <- get_pdf_device()

message("Panel A: cluster profiles with subject spaghetti...")

PA_W <- 120  # panel width mm
PA_H <- 300  # stacked height mm

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

# --- A2. Per-cluster summary ---
panel_a_summary <- panel_a_long %>%
  group_by(cluster, age, time, time_num) %>%
  summarise(
    mean_z = mean(z_score, na.rm = TRUE),
    se_z   = sd(z_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- A3. Shared y-axis range ---
y_range <- panel_a_summary %>%
  summarise(y_lo = min(mean_z - 1.96 * se_z), y_hi = max(mean_z + 1.96 * se_z))
y_pad   <- (y_range$y_hi - y_range$y_lo) * 0.1
y_limits <- c(y_range$y_lo - y_pad, y_range$y_hi + y_pad)

# --- A4. Text sizes ---
txt_title <- scale_text(BASE_PATHWAY, PA_W) * 1.3
txt_sub   <- scale_text(BASE_STAT, PA_W)
txt_axis  <- scale_text(BASE_STAT, PA_W)
txt_tick  <- scale_text(BASE_GENE, PA_W)

# --- A5. Build one plot per cluster ---
cluster_ids <- paste0("C", seq_len(optimal_k))
n_clusters  <- length(cluster_ids)

panels_A <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]
  cl_data    <- panel_a_long %>% filter(cluster == cid)
  cl_summary <- panel_a_summary %>% filter(cluster == cid)
  n_core <- n_distinct(cl_data$gene)

  bio_label <- top_hallmark$label[top_hallmark$cluster == cid]
  if (length(bio_label) == 0) bio_label <- ""

  is_last <- (i == n_clusters)

  ggplot() +
    geom_line(data = cl_data,
              aes(x = time_num, y = z_score, group = interaction(gene, age),
                  colour = age, alpha = membership),
              linewidth = 0.15) +
    geom_ribbon(data = cl_summary,
                aes(x = time_num, ymin = mean_z - se_z, ymax = mean_z + se_z,
                    fill = age, group = age),
                alpha = 0.15) +
    geom_line(data = cl_summary,
              aes(x = time_num, y = mean_z, colour = age, group = age),
              linewidth = 1.2) +
    geom_point(data = cl_summary,
               aes(x = time_num, y = mean_z, colour = age),
               size = 2.5) +
    annotate("segment", x = 0.7, xend = 0.7,
             y = y_limits[1], yend = y_limits[2],
             colour = CLUSTER_COLORS[cid], linewidth = 3, lineend = "butt") +
    scale_colour_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    scale_alpha_continuous(range = c(0.02, 0.15), guide = "none") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = if (is_last) c("Pre", "Post") else NULL,
      limits = c(0.7, 2.3), expand = expansion(0)
    ) +
    scale_y_continuous(limits = y_limits, name = NULL) +
    labs(
      title    = sprintf("C%d: %s", i, bio_label),
      subtitle = sprintf("(n = %d)", n_core),
      x        = if (is_last) "Time" else NULL
    ) +
    FIG_THEME +
    theme(
      plot.title    = element_text(colour = CLUSTER_COLORS[cid], face = "bold",
                                   size = txt_title, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey30", face = "bold.italic",
                                   size = txt_sub, hjust = 0.5),
      panel.border  = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
      axis.title.y  = element_blank(),
      axis.title.x  = if (is_last) element_text(size = txt_axis, face = "bold") else element_blank(),
      axis.text.y   = element_text(size = txt_tick, face = "bold"),
      axis.text.x   = if (is_last) element_text(size = txt_tick, face = "bold") else element_blank(),
      axis.ticks.x  = if (is_last) element_line() else element_blank(),
      legend.position = "none",
      plot.margin   = margin(t = 2, r = -2, b = if (is_last) 4 else 1, l = 2)
    )
})

# --- Save stacked panel ---
col_A <- wrap_plots(panels_A, ncol = 1, heights = row_heights)
ggsave(file.path(RPT, "panel_A_profiles.pdf"), col_A,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_A_profiles.png"), col_A,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

message("  Panel A saved")
