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

txt_quad <- scale_text(BASE_QUADRANT, PE_W)

# Hotspot genes extracted at peak hypergeometric enrichment per quadrant
# (Cahill et al. 2018 RRHO2) — used for heatmap annotation and ORA.
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

n_UU <- length(hotspot_genes$UU)
n_DD <- length(hotspot_genes$DD)
n_UD <- length(hotspot_genes$UD)
n_DU <- length(hotspot_genes$DU)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

txt_stat_h <- scale_text(BASE_STAT, PE_W) * 0.75

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = expression(sign %*% -log[10](P)),
    guide = guide_colorbar(
      barwidth = unit(25, "mm"), barheight = unit(2.5, "mm"),
      title.position = "left", title.hjust = 1,
      title.theme = element_text(size = 7, face = "bold"))) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.4, alpha = 0.7) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.4, alpha = 0.7) +
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
  # Quadrant stats: hotspot gene count at peak enrichment
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

  MAX_PER_QUAD <- 12
  bar_df <- ora_all %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      pathway_label  = str_trunc(clean_pathway_name(pathway), 40),
      quadrant       = factor(quadrant, levels = all_quadrant_names)
    ) %>%
    group_by(quadrant) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = MAX_PER_QUAD) %>%
    ungroup() %>%
    arrange(quadrant, neg_log10_padj) %>%
    mutate(uid = fct_inorder(paste0(pathway_label, "___", quadrant)))

  n_per_quad <- bar_df %>% count(quadrant) %>% deframe()
  n_shown    <- nrow(bar_df)
  n_total    <- nrow(ora_all)

  pE_ora <- ggplot(bar_df, aes(x = neg_log10_padj, y = uid, fill = quadrant)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = overlap), hjust = -0.3, size = txt_ora * 0.7,
              color = "grey30") +
    scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
    scale_fill_manual(values = ORA_QUAD_COLORS_F2, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(quadrant ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(quadrant = function(x) str_wrap(x, 18))) +
    labs(title = "Enriched Pathways by Concordance Quadrant",
         subtitle = if (n_shown < n_total)
           sprintf("Top %d per quadrant (%d terms total)", MAX_PER_QUAD, n_total)
         else sprintf("%d terms total", n_total),
         x = expression(-log[10](p[adj])),
         y = NULL) +
    FIG_THEME +
    theme(
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 8, hjust = 0.5),
      strip.text.y     = element_text(size = 7, face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y  = element_text(size = 7),
      plot.margin  = margin(2, 4, 2, 2, "mm")
    )
}

pE <- (pE_heat | pE_ora) + plot_layout(widths = c(1, 1.3))

PE_TOTAL_W <- 400
PE_TOTAL_H <- 220
ggsave(file.path(RPT, "panel_E_rrho2.pdf"), pE,
       width = PE_TOTAL_W, height = PE_TOTAL_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_rrho2.png"), pE,
       width = PE_TOTAL_W, height = PE_TOTAL_H, units = "mm", dpi = 300)

rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

message("F2 Panel E done")
