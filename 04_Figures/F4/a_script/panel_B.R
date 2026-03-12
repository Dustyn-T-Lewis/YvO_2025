# Figure 4 — Panel B: PCA highlight-on-grey scatter

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(patchwork)
})

RPT <- "04_Figures/F4/b_reports"
DAT <- "04_Figures/F4/c_data"
optimal_k <- 4

core_proteins <- read_csv(file.path(DAT, "06_mfuzz_assignments.csv"),
                          show_col_types = FALSE) %>%
  filter(membership >= 0.5)

delta_mat <- readRDS(file.path(DAT, "delta_mat.rds"))

row_heights <- core_proteins %>%
  count(cluster) %>%
  arrange(cluster) %>%
  pull(n)
row_heights <- row_heights / sum(row_heights)

pdf_device <- get_pdf_device()

message("Panel B: PCA highlight-on-grey scatter...")

PB_W <- 120  # panel width mm
PB_H <- 300  # stacked height mm

pca_result <- prcomp(t(delta_mat), center = TRUE, scale. = FALSE)
pca_coords <- as.data.frame(pca_result$rotation[, 1:2])
pca_coords$gene <- rownames(pca_coords)

var_explained <- summary(pca_result)$importance[2, 1:2] * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])

pca_scores <- pca_coords %>%
  left_join(core_proteins %>% dplyr::select(gene, cluster, membership), by = "gene") %>%
  mutate(
    cluster    = ifelse(is.na(cluster), NA_character_, cluster),
    membership = ifelse(is.na(membership), 0, membership)
  )

write_csv(pca_scores, file.path(DAT, "02_panel_B_pca_scores.csv"))

txt_annot <- scale_text(BASE_GENE, PB_W)
txt_axis  <- scale_text(BASE_STAT, PB_W)
txt_tick  <- scale_text(BASE_GENE, PB_W)

cluster_ids <- paste0("C", seq_len(optimal_k))
n_clusters  <- length(cluster_ids)

panels_B <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]
  is_first <- (i == 1)
  is_last  <- (i == n_clusters)

  highlight_data <- pca_scores %>% filter(cluster == cid)

  ggplot() +
    geom_point(data = pca_scores,
               aes(x = PC1, y = PC2),
               color = "grey55", size = 0.3, alpha = 0.45) +
    geom_point(data = highlight_data,
               aes(x = PC1, y = PC2, alpha = membership),
               color = CLUSTER_COLORS[cid], size = 0.6) +
    scale_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
    annotate("text",
             x = min(pca_scores$PC1), y = max(pca_scores$PC2),
             label = cid,
             color = CLUSTER_COLORS[cid], fontface = "bold", size = txt_annot,
             hjust = 0, vjust = 1) +
    annotate("segment",
             x = min(pca_scores$PC1), xend = min(pca_scores$PC1),
             y = min(pca_scores$PC2), yend = max(pca_scores$PC2),
             colour = CLUSTER_COLORS[cid], linewidth = 3, lineend = "butt") +
    labs(
      x = if (is_last) pc1_lab else NULL,
      y = if (is_first) pc2_lab else NULL
    ) +
    FIG_THEME +
    theme(
      panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
      axis.title.x = if (is_last) element_text(size = txt_axis, face = "bold") else element_blank(),
      axis.title.y = if (is_first) element_text(size = txt_axis, face = "bold") else element_blank(),
      axis.text.x  = if (is_last) element_text(size = txt_tick, face = "bold") else element_blank(),
      axis.text.y  = element_text(size = txt_tick, face = "bold"),
      axis.ticks.x = if (is_last) element_line() else element_blank(),
      legend.position = "none",
      plot.margin  = margin(t = 2, r = 2, b = if (is_last) 4 else 1, l = -2)
    )
})

col_B <- wrap_plots(panels_B, ncol = 1, heights = row_heights)
ggsave(file.path(RPT, "panel_B_pca.pdf"), col_B,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_pca.png"), col_B,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

message("  Panel B saved")
