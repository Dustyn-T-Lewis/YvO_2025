# F3 Panel D: RRHO2 Reversal Map + Per-Quadrant ORA
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
RPT <- "04_Figures/F3/b_reports"
DAT <- "04_Figures/F3/c_data"
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

mid_aging <- rank_aging[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
bot_aging <- rank_aging[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_rev_AgUp_TrDn <- length(intersect(mid_aging, mid_old))
n_rev_AgDn_TrUp <- length(intersect(bot_aging, bot_old))
n_exac_Up       <- length(intersect(mid_aging, bot_old))
n_exac_Down     <- length(intersect(bot_aging, mid_old))

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


txt_quad <- scale_text(BASE_QUADRANT, PD_W)
txt_stat <- scale_text(BASE_STAT, PD_W)
txt_ora  <- scale_text(BASE_STAT, PD_W)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

txt_stat_h <- scale_text(BASE_STAT, PD_W) * 0.75

pD_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = expression(-log[10](P)),
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

  TOP_N_LABEL <- 10
  top_ids <- ora_all %>%
    arrange(desc(overlap)) %>%
    slice_head(n = TOP_N_LABEL) %>%
    pull(ID)

  ora_plot_df <- ora_all %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           gene_ratio     = overlap / size,
           pathway_label  = ifelse(ID %in% top_ids,
                                   str_wrap(clean_pathway_name(ID), width = 28),
                                   NA_character_),
           quadrant       = factor(quadrant, levels = names(ORA_QUAD_COLORS_F3)))

  # Three-zone piecewise: 0-4 ~65%, 4-8 ~25%, 8-20 ~10% of visual height
  B1 <- 4; B2 <- 8; YMAX <- 20
  # Target visual proportions: 0–B1 = 65%, B1–B2 = 25%, B2–YMAX = 10%
  V1 <- 6.5; V2 <- 2.5; V3 <- 1.0  # visual units allocated to each zone
  pw_fwd <- function(y) {
    ifelse(y <= B1, y / B1 * V1,
    ifelse(y <= B2, V1 + (y - B1) / (B2 - B1) * V2,
                    V1 + V2 + (y - B2) / (YMAX - B2) * V3))
  }
  pw_inv <- function(t) {
    ifelse(t <= V1, t / V1 * B1,
    ifelse(t <= V1 + V2, B1 + (t - V1) / V2 * (B2 - B1),
                         B2 + (t - V1 - V2) / V3 * (YMAX - B2)))
  }
  pw_trans <- scales::trans_new("piecewise", transform = pw_fwd, inverse = pw_inv)

  pD_ora <- ggplot(ora_plot_df,
                   aes(x = gene_ratio, y = neg_log10_padj, color = quadrant)) +
    geom_point(aes(size = overlap), alpha = 0.75) +
    geom_label_repel(
      aes(label = pathway_label, fill = quadrant),
      color = "white", fontface = "bold",
      size = txt_ora * 0.8, lineheight = 0.85, na.rm = TRUE,
      max.overlaps = 30, box.padding = 0.5, point.padding = 0.3,
      segment.size = 0.2, segment.color = "grey50",
      min.segment.length = 0.2,
      label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
      label.size = 0.3, show.legend = FALSE, seed = 42
    ) +
    scale_color_manual(
      values = ORA_QUAD_COLORS_F3, name = "Quadrant", drop = FALSE,
      labels = c("Reversed (Aging Up / Training Down)"  = "Reversed (Up->Down)",
                 "Reversed (Aging Down / Training Up)"  = "Reversed (Down->Up)",
                 "Exacerbated Up" = "Exacerbated Up",
                 "Exacerbated Down" = "Exacerbated Down")
    ) +
    scale_fill_manual(values = ORA_QUAD_COLORS_F3, guide = "none",
                      drop = FALSE) +
    scale_size_continuous(range = c(2, 10), name = "Overlap") +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    scale_y_continuous(trans = pw_trans,
                       breaks = c(0, 2, 4, 6, 8, 20),
                       expand = expansion(mult = c(0.04, 0.06))) +
    labs(title = "Enriched Pathways by Reversal Quadrant",
         subtitle = sprintf("Top %d labeled by overlap (%d terms total)",
                            TOP_N_LABEL, nrow(ora_all)),
         x = "Gene Ratio (overlap / pathway size)",
         y = expression(-log[10](p[adj]))) +
    FIG_THEME +
    theme(panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.box       = "horizontal",
          legend.text      = element_text(size = 8),
          legend.title     = element_text(size = 9, face = "bold"),
          plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.margin      = margin(2, 4, 2, 2, "mm")) +
    guides(
      color = guide_legend(title = "Quadrant", override.aes = list(size = 4),
                           ncol = 2),
      size  = guide_legend(title = "Overlap", nrow = 1)
    )

  y_range <- pw_fwd(YMAX) - pw_fwd(0)
  npc_8  <- 0.04 + (pw_fwd(B2)   - pw_fwd(0)) / y_range * 0.92
  npc_20 <- 0.04 + (pw_fwd(YMAX) - pw_fwd(0)) / y_range * 0.92
  npc_mid <- (npc_8 + npc_20) / 2   # midpoint of the compressed zone

  sl <- 0.008   # half-length of each slash
  gap <- 0.005  # gap between the two slashes

  slash_grob <- grobTree(
    # White rectangle to mask the axis line in the gap
    rectGrob(x = unit(0, "npc"), y = unit(npc_mid, "npc"),
             width = unit(0.025, "npc"), height = unit(gap * 3, "npc"),
             gp = gpar(fill = "white", col = NA)),
    # Lower slash
    segmentsGrob(x0 = unit(-0.01, "npc"), x1 = unit(0.01, "npc"),
                 y0 = unit(npc_mid - gap - sl, "npc"),
                 y1 = unit(npc_mid - gap + sl, "npc"),
                 gp = gpar(col = "grey30", lwd = 1.2)),
    # Upper slash
    segmentsGrob(x0 = unit(-0.01, "npc"), x1 = unit(0.01, "npc"),
                 y0 = unit(npc_mid + gap - sl, "npc"),
                 y1 = unit(npc_mid + gap + sl, "npc"),
                 gp = gpar(col = "grey30", lwd = 1.2))
  )

  pD_ora <- pD_ora +
    annotation_custom(slash_grob, xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf) +
    coord_cartesian(clip = "off")

} else {
  pD_ora <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("No enrichment at padj < 0.05\n(%d shared genes)",
                             n_shared),
             size = txt_ora, color = "grey50", fontface = "italic") +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm"))
}

pD <- (pD_heat | pD_ora) + plot_layout(widths = c(1, 1))

PD_TOTAL_W <- 380
PD_TOTAL_H <- 210
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

cat("F3 Panel D done\n")
