# F05 Panel D: RRHO2 Reversal Map + Per-Quadrant ORA
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

PW <- 220; PH <- 260
RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

pdf_device <- get_pdf_device()
PD_W <- 220

rr_df <- dep_df %>%
  transmute(gene, t_aging = t_Aging, t_old = t_Training_Old) %>%
  filter(!is.na(t_aging), !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared   <- nrow(rr_df)
message(sprintf("  %d shared genes for RRHO", n_shared))
rank_aging <- rr_df$gene[order(-rr_df$t_aging)]    # index 1 = most UP with aging
rank_old   <- rr_df$gene[order(rr_df$t_old)]       # index 1 = most DOWN with training


step    <- max(1L, floor(n_shared / 200))
indices <- seq(1L, n_shared, by = step)
if (tail(indices, 1) != n_shared) indices <- c(indices, n_shared)
n_grid  <- length(indices)
message(sprintf("  Grid: %d x %d (step = %d)", n_grid, n_grid, step))

hmat <- matrix(0, nrow = n_grid, ncol = n_grid)
for (ii in seq_along(indices)) {
  top_aging <- rank_aging[1:indices[ii]]
  for (jj in seq_along(indices)) {
    top_old  <- rank_old[1:indices[jj]]
    overlap  <- length(intersect(top_aging, top_old))
    p_upper  <- phyper(overlap - 1, indices[ii], n_shared - indices[ii],
                       indices[jj], lower.tail = FALSE)
    p_lower  <- phyper(overlap, indices[ii], n_shared - indices[ii],
                       indices[jj], lower.tail = TRUE)
    p_val    <- min(2 * min(p_upper, p_lower), 1)
    expected <- indices[ii] * indices[jj] / n_shared
    sign_fac <- ifelse(overlap >= expected, 1, -1)
    hmat[ii, jj] <- sign_fac * (-log10(max(p_val, .Machine$double.xmin)))
  }
}
nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

# Gene sets extracted at the point of maximum hypergeometric enrichment
# within each quadrant (Cahill et al. 2018 RRHO2 framework). This selects the
# statistically optimal rank threshold pair rather than an arbitrary midpoint.

# Aging sorted descending (UP first), Training sorted ascending (DOWN first)
# Bottom-left (small x = Aging UP, small y = Training DOWN) = Reversed (Aging Up / Training Down)
# Top-right   (large x = Aging DOWN, large y = Training UP) = Reversed (Aging Down / Training Up)
# Top-left    (small x = Aging UP,   large y = Training UP) = Exacerbated Up
# Bottom-right(large x = Aging DOWN, small y = Training DOWN) = Exacerbated Down
max_rev_AgUp_TrDn <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_rev_AgDn_TrUp <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_exac_Up       <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_exac_Down     <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

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
  UU = intersect(rank_aging[1:indices[peak_UU$i]], rank_old[1:indices[peak_UU$j]]),
  DD = intersect(rank_aging[indices[peak_DD$i]:n_shared],
                 rank_old[indices[peak_DD$j]:n_shared]),
  UD = intersect(rank_aging[1:indices[peak_UD$i]],
                 rank_old[indices[peak_UD$j]:n_shared]),
  DU = intersect(rank_aging[indices[peak_DU$i]:n_shared],
                 rank_old[1:indices[peak_DU$j]])
)

# Hotspot gene counts (at peak enrichment, not arbitrary midpoint)
n_rev_AgUp_TrDn <- length(hotspot_genes$UU)
n_rev_AgDn_TrUp <- length(hotspot_genes$DD)
n_exac_Up       <- length(hotspot_genes$UD)
n_exac_Down     <- length(hotspot_genes$DU)


pw_collection_D <- build_pathway_collection(min_size = 10, max_size = 500)
all_genes <- rr_df$gene

run_ora <- function(gene_set, label) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes,
      pathways       = pw_collection_D,
      jaccard_cutoff = 0.5,
      min_size       = 10,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) == 0) return(tibble())
  res %>%
    mutate(quadrant = label,
           pathway_label = clean_pathway_name(pathway),
           ID = pathway,
           p.adjust = padj,
           GeneRatio = paste0(overlap, "/", length(gene_set)),
           geneID = sapply(overlapGenes, paste, collapse = "/")) %>%
    arrange(padj, size)
}

ora_reversal <- bind_rows(
  run_ora(hotspot_genes$UU, "Reversed (Aging Up / Training Down)"),
  run_ora(hotspot_genes$DD, "Reversed (Aging Down / Training Up)"))

ora_exacerbation <- bind_rows(
  run_ora(hotspot_genes$UD, "Exacerbated Up"),
  run_ora(hotspot_genes$DU, "Exacerbated Down"))

write_csv(ora_reversal, file.path(DAT, "panel_D", "rrho2_ora_concordant.csv"))
write_csv(ora_exacerbation, file.path(DAT, "panel_D", "rrho2_ora_discordant.csv"))

# Exacerbation quadrant ORA: document empty result if no pathways survived BH
if (nrow(ora_exacerbation) == 0) {
  write_csv(tibble(note = "No pathways enriched in exacerbation quadrants (padj<0.05)"),
            file.path(DAT, "panel_D", "rrho2_ora_discordant_note.csv"))
}


txt_quad <- scale_text(BASE_QUADRANT, PD_W)
txt_stat <- scale_text(BASE_STAT, PD_W)
txt_ora  <- scale_text(BASE_STAT, PD_W)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

txt_stat_h <- scale_text(BASE_STAT, PD_W) * 0.75

pD_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
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
           label = "Reversed\nAging Up / Training Down",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5,
           y = mid_c + (nc - mid_c) * 0.5,
           label = "Reversed\nAging Down / Training Up",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = "Exacerbated Up",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = "Exacerbated Down",
           color = "white", fontface = "bold", size = txt_quad) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5 - mid_c * 0.18,
           label = sprintf("max = %.1f | n = %d", max_rev_AgUp_TrDn, n_rev_AgUp_TrDn),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5,
           y = mid_c + (nc - mid_c) * 0.5 - (nc - mid_c) * 0.18,
           label = sprintf("max = %.1f | n = %d", max_rev_AgDn_TrUp, n_rev_AgDn_TrUp),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r * 0.5,
           y = mid_c + (nc - mid_c) * 0.5 - (nc - mid_c) * 0.15,
           label = sprintf("max = %.1f | n = %d", max_exac_Up, n_exac_Up),
           color = "white", size = txt_stat_h) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5 - mid_c * 0.15,
           label = sprintf("max = %.1f | n = %d", max_exac_Down, n_exac_Down),
           color = "white", size = txt_stat_h) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Threshold-Free Reversal (RRHO)",
       subtitle = sprintf("Two-sided hypergeometric | %d shared genes | step = %d",
                          n_shared, step),
       x = expression("Aging rank"~(Up %->% Down)),
       y = expression("Training (Old) rank"~(Down %->% Up))) +
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

ora_all <- bind_rows(ora_reversal, ora_exacerbation)

if (nrow(ora_all) > 0) {

  all_quadrant_names <- names(ORA_QUAD_COLORS_F3)
  # Short labels for facet strips
  quad_short <- c(
    "Reversed (Aging Up / Training Down)"  = "Reversed\n(Up \u2192 Down)",
    "Reversed (Aging Down / Training Up)"  = "Reversed\n(Down \u2192 Up)",
    "Exacerbated Up"                       = "Exacerbated Up",
    "Exacerbated Down"                     = "Exacerbated Down"
  )

  MAX_PER_QUAD <- 12
  bar_df <- ora_all %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      pathway_label  = str_trunc(clean_pathway_name(pathway), 40),
      quadrant       = factor(quadrant, levels = all_quadrant_names)
    ) %>%
    filter(!is.na(quadrant)) %>%
    group_by(quadrant, .drop = FALSE) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = MAX_PER_QUAD) %>%
    ungroup() %>%
    filter(!is.na(neg_log10_padj)) %>%
    arrange(quadrant, neg_log10_padj) %>%
    mutate(uid = fct_inorder(paste0(pathway_label, "___", quadrant)))

  n_shown <- nrow(bar_df)
  n_total <- nrow(ora_all)

  pD_ora <- ggplot(bar_df, aes(x = neg_log10_padj, y = uid, fill = quadrant)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = overlap), hjust = -0.3, size = txt_ora * 0.7,
              color = "grey30") +
    scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
    scale_fill_manual(values = ORA_QUAD_COLORS_F3, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(quadrant ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(quadrant = quad_short)) +
    labs(title = "Enriched Pathways by Reversal Quadrant",
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

} else {
  pD_ora <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("No enrichment at padj < 0.05\n(%d shared genes)",
                             n_shared),
             size = txt_ora, color = "grey50", fontface = "italic") +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm"))
}

pD <- (pD_heat | pD_ora) + plot_layout(widths = c(1, 1.3))

PD_TOTAL_W <- 400
PD_TOTAL_H <- 220
ggsave(file.path(RPT, "panel_D_rrho2.pdf"), pD,
       width = PD_TOTAL_W, height = PD_TOTAL_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_rrho2.png"), pD,
       width = PD_TOTAL_W, height = PD_TOTAL_H, units = "mm", dpi = 300)

n_ora_rev_UpDn <- if (nrow(ora_reversal) > 0) sum(ora_reversal$quadrant == "Reversed (Aging Up / Training Down)") else 0L
n_ora_rev_DnUp <- if (nrow(ora_reversal) > 0) sum(ora_reversal$quadrant == "Reversed (Aging Down / Training Up)") else 0L
n_ora_exac_Up  <- if (nrow(ora_exacerbation) > 0) sum(ora_exacerbation$quadrant == "Exacerbated Up") else 0L
n_ora_exac_Dn  <- if (nrow(ora_exacerbation) > 0) sum(ora_exacerbation$quadrant == "Exacerbated Down") else 0L

tibble(
  quadrant       = c("Reversed_AgingUp_TrainingDown", "Reversed_AgingDown_TrainingUp",
                     "Exacerbated_Up", "Exacerbated_Down"),
  max_neg_log10p = round(c(max_rev_AgUp_TrDn, max_rev_AgDn_TrUp,
                           max_exac_Up, max_exac_Down), 2),
  n_overlap      = c(n_rev_AgUp_TrDn, n_rev_AgDn_TrUp, n_exac_Up, n_exac_Down),
  n_hotspot_genes = c(length(hotspot_genes$UU), length(hotspot_genes$DD),
                      length(hotspot_genes$UD), length(hotspot_genes$DU)),
  n_ora_pathways  = c(n_ora_rev_UpDn, n_ora_rev_DnUp, n_ora_exac_Up, n_ora_exac_Dn),
  n_shared_genes  = n_shared
) %>%
  write_csv(file.path(DAT, "panel_D", "rrho2_summary.csv"))

cat("F05 Panel D done\n")
