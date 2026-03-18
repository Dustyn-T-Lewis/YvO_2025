# F2 Panel E: RRHO Concordance Map + Per-Quadrant ORA
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(msigdbr)
  library(fgsea)
})

PE_W <- 220

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr_df <- dep_df %>%
  transmute(gene, t_young = t_Training_Young, t_old = t_Training_Old) %>%
  filter(!is.na(t_young) & !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)
rank_young <- rr_df$gene[order(-rr_df$t_young)]   # index 1 = most UP
rank_old   <- rr_df$gene[order(-rr_df$t_old)]    # index 1 = most UP

# ~200 grid points per axis
step <- max(1, floor(n_shared / 200))
indices <- seq(1, n_shared, by = step)
if (tail(indices, 1) != n_shared) indices <- c(indices, n_shared)
n_grid <- length(indices)

# Compute -log10(p) hypergeometric overlap matrix (two-sided)
hmat <- matrix(0, nrow = n_grid, ncol = n_grid)
for (ii in seq_along(indices)) {
  top_young <- rank_young[1:indices[ii]]
  for (jj in seq_along(indices)) {
    top_old <- rank_old[1:indices[jj]]
    overlap <- length(intersect(top_young, top_old))
    p_upper <- phyper(overlap - 1, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = FALSE)
    p_lower <- phyper(overlap, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = TRUE)
    p_val <- min(2 * min(p_upper, p_lower), 1)
    expected <- indices[ii] * indices[jj] / n_shared
    sign_factor <- ifelse(overlap >= expected, 1, -1)
    hmat[ii, jj] <- sign_factor * (-log10(max(p_val, .Machine$double.xmin)))
  }
}

nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

# Both sorted descending (UP first):
# Bottom-left (small x = Y UP, small y = O UP) = Concordant Up (UU)
# Top-right   (large x = Y DOWN, large y = O DOWN) = Concordant Down (DD)
# Top-left    (small x = Y UP,   large y = O DOWN) = Discordant (Y Up / O Down) (UD)
# Bottom-right(large x = Y DOWN, small y = O UP) = Discordant (Y Down / O Up) (DU)
max_UU <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_DD <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_UD <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_DU <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

mid_young <- rank_young[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
n_UU <- length(intersect(mid_young, mid_old))
bot_young <- rank_young[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_DD <- length(intersect(bot_young, bot_old))
n_UD <- length(intersect(mid_young, bot_old))
n_DU <- length(intersect(bot_young, mid_old))

txt_quad <- scale_text(BASE_QUADRANT, PE_W)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

txt_stat_h <- scale_text(BASE_STAT, PE_W) * 0.75

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = expression(-log[10](P)),
                        guide = guide_colorbar(
                          barwidth = unit(25, "mm"), barheight = unit(2.5, "mm"),
                          title.position = "left", title.hjust = 1,
                          title.theme = element_text(size = 7, face = "bold"))) +
  # Quadrant boundary lines
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.4, alpha = 0.7) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.4, alpha = 0.7) +
  # Quadrant names — x: Y UP→DOWN (left=UP), y: O UP→DOWN (bottom=UP)
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = "Concordant Up",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = "Concordant Down",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = "Discordant\nY Up / O Down",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = "Discordant\nY Down / O Up",
           color = "white", fontface = "bold", size = txt_quad) +
  # Quadrant stats (smaller, below names)
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5 - mid_c * 0.15,
           label = sprintf("max = %.1f | n = %d", max_UU, n_UU),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5,
           y = mid_c + (nc - mid_c) * 0.5 - (nc - mid_c) * 0.15,
           label = sprintf("max = %.1f | n = %d", max_DD, n_DD),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r * 0.5,
           y = mid_c + (nc - mid_c) * 0.5 - (nc - mid_c) * 0.18,
           label = sprintf("max = %.1f | n = %d", max_UD, n_UD),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5 - mid_c * 0.18,
           label = sprintf("max = %.1f | n = %d", max_DU, n_DU),
           color = "white", size = txt_stat_h) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Threshold-Free Concordance (RRHO)",
    subtitle = sprintf("Two-sided hypergeometric | %d shared genes | step = %d",
                        n_shared, step),
    x = expression("Training (Young) rank"~(Up %->% Down)),
    y = expression("Training (Old) rank"~(Up %->% Down))
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 2)),
    axis.title.y     = element_text(margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.margin    = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  coord_fixed(ratio = 1)

# Hotspot genes extracted at peak hypergeometric enrichment per quadrant
# (Cahill et al. 2018 RRHO2), not the arbitrary midpoint used for heatmap counts.

find_peak <- function(mat, rows, cols) {
  sub_mat <- mat[rows, cols, drop = FALSE]
  peak <- which(sub_mat == max(sub_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(i = rows[peak[1]], j = cols[peak[2]])
}

peak_UU <- find_peak(hmat, 1:mid_r, 1:mid_c)
peak_DD <- find_peak(hmat, (mid_r+1):nr, (mid_c+1):nc)
peak_UD <- find_peak(hmat, 1:mid_r, (mid_c+1):nc)
peak_DU <- find_peak(hmat, (mid_r+1):nr, 1:mid_c)

hotspot_genes <- list(
  UU = intersect(rank_young[1:indices[peak_UU$i]], rank_old[1:indices[peak_UU$j]]),
  DD = intersect(rank_young[indices[peak_DD$i]:n_shared], rank_old[indices[peak_DD$j]:n_shared]),
  UD = intersect(rank_young[1:indices[peak_UD$i]], rank_old[indices[peak_UD$j]:n_shared]),
  DU = intersect(rank_young[indices[peak_DU$i]:n_shared], rank_old[1:indices[peak_DU$j]])
)

hotspot_export <- bind_rows(
  tibble(quadrant = "UU", gene = hotspot_genes$UU),
  tibble(quadrant = "DD", gene = hotspot_genes$DD),
  tibble(quadrant = "UD", gene = hotspot_genes$UD),
  tibble(quadrant = "DU", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

pw_collection_E <- build_pathway_collection(min_size = 10, max_size = 500)
all_genes_E <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes_E,
      pathways       = pw_collection_E,
      jaccard_cutoff = 0.5,
      min_size       = 10,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) > 0) {
    res %>%
      mutate(quadrant = quadrant_name,
             pathway_label = clean_pathway_name(pathway),
             ID = pathway,
             p.adjust = padj,
             GeneRatio = paste0(overlap, "/", length(gene_set)),
             geneID = sapply(overlapGenes, paste, collapse = "/")) %>%
      arrange(padj, size)
  } else {
    tibble()
  }
}

ora_UU <- run_quadrant_ora(hotspot_genes$UU, "Concordant Up")
ora_DD <- run_quadrant_ora(hotspot_genes$DD, "Concordant Down")
ora_UD <- run_quadrant_ora(hotspot_genes$UD, "Discordant (Y Up / O Down)")
ora_DU <- run_quadrant_ora(hotspot_genes$DU, "Discordant (Y Down / O Up)")

ora_concordant  <- bind_rows(ora_UU, ora_DD)
ora_discordant  <- bind_rows(ora_UD, ora_DU)

write_csv(ora_concordant,  file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_discordant,  file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

ora_all <- bind_rows(ora_UU, ora_DD, ora_UD, ora_DU)
txt_ora <- scale_text(BASE_STAT, PE_W)

if (nrow(ora_all) == 0) {
  pE_ora <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No enrichment at padj < 0.05",
             size = txt_ora, color = "grey50", fontface = "italic") +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm"))
} else {
  all_quadrant_names <- c("Concordant Up", "Concordant Down",
                          "Discordant (Y Up / O Down)", "Discordant (Y Down / O Up)")

  LABEL_MIN_OVERLAP <- 20
  n_labeled <- sum(ora_all$overlap >= LABEL_MIN_OVERLAP)
  n_unlabeled <- nrow(ora_all) - n_labeled

  plot_df <- ora_all %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      gene_ratio     = overlap / size,
      pathway_label  = ifelse(overlap >= LABEL_MIN_OVERLAP,
                              str_wrap(clean_pathway_name(pathway), width = 28),
                              NA_character_),
      quadrant       = factor(quadrant, levels = all_quadrant_names)
    )

  pE_ora <- ggplot(plot_df, aes(x = gene_ratio, y = neg_log10_padj,
                                color = quadrant)) +
    geom_point(aes(size = overlap), alpha = 0.75) +
    geom_label_repel(
      aes(label = pathway_label, fill = quadrant),
      color = "white", fontface = "bold",
      size        = txt_ora * 0.8,
      lineheight  = 0.85,
      na.rm       = TRUE,
      max.overlaps = 30,
      box.padding  = 0.5,
      point.padding = 0.3,
      segment.size  = 0.2,
      segment.color = "grey50",
      min.segment.length = 0.2,
      label.padding = unit(1, "pt"),
      label.r       = unit(1, "pt"),
      label.size     = 0.3,
      show.legend  = FALSE,
      seed = 42
    ) +
    scale_color_manual(values = ORA_QUAD_COLORS_F2, name = "Quadrant") +
    scale_fill_manual(values = ORA_QUAD_COLORS_F2, guide = "none") +
    scale_size_continuous(range = c(2, 10), name = "Overlap") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = "Enriched Pathways by Concordance Quadrant",
         subtitle = if (n_unlabeled > 0)
           sprintf("Labeled: overlap >= %d (%d of %d terms)",
                   LABEL_MIN_OVERLAP, n_labeled, nrow(ora_all)) else NULL,
         x = "Gene Ratio (overlap / pathway size)",
         y = expression(-log[10](p[adj]))) +
    FIG_THEME +
    theme(
      plot.title     = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.text    = element_text(size = 8),
      legend.title   = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      plot.margin = margin(2, 4, 2, 2, "mm")
    ) +
    guides(
      color = guide_legend(title = "Quadrant", override.aes = list(size = 4), nrow = 2),
      size = guide_legend(title = "Overlap", nrow = 1)
    )
}

pE <- (pE_heat | pE_ora) + plot_layout(widths = c(1, 1))

PE_TOTAL_W <- 380
PE_TOTAL_H <- 210
ggsave(file.path(RPT, "panel_E_rrho2.pdf"), pE,
       width = PE_TOTAL_W, height = PE_TOTAL_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_rrho2.png"), pE,
       width = PE_TOTAL_W, height = PE_TOTAL_H, units = "mm", dpi = 300)

rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_overlap = c(n_UU, n_DD, n_UD, n_DU),
  n_hotspot_genes = c(length(hotspot_genes$UU), length(hotspot_genes$DD),
                      length(hotspot_genes$UD), length(hotspot_genes$DU)),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

message("F2 Panel E done")
