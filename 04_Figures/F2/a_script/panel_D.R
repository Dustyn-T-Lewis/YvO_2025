################################################################################
#   Figure 2 — Panel D: RRHO Concordance Map
#   Rank-Rank Hypergeometric Overlap (Plaisier et al. 2010, NAR 38:e169)
#   Pure-R implementation to avoid RedRibbon C stack overflow.
#
#   Generates: panel_D_rrho2.pdf
#              + c_data/panel_D/rrho2_summary.csv, c_data/panel_D/rrho2_matrix.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_F2_setup.R")

# ==============================================================================
# PANEL D — RRHO Concordance Map (pure-R hypergeometric)
# ==============================================================================

message("Panel D: RRHO concordance map...")

# Prepare ranked lists sorted from most upregulated to most downregulated
rr_df <- dep_df %>%
  transmute(gene, t_young = t_Training_Young, t_old = t_Training_Old) %>%
  filter(!is.na(t_young) & !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)
message(sprintf("  %d shared genes for RRHO", n_shared))

# Rank lists: index 1 = most upregulated (highest t-stat)
rank_young <- rr_df$gene[order(-rr_df$t_young)]
rank_old   <- rr_df$gene[order(-rr_df$t_old)]

# Step size: ~200 grid points per axis for reasonable resolution
step <- max(1, floor(n_shared / 200))
indices <- seq(1, n_shared, by = step)
if (tail(indices, 1) != n_shared) indices <- c(indices, n_shared)
n_grid <- length(indices)
message(sprintf("  Grid: %d x %d (step = %d)", n_grid, n_grid, step))

# Compute -log10(p) hypergeometric overlap matrix
# hmat[i,j] = -log10(phyper) for overlap of top-i in list1 and top-j in list2
hmat <- matrix(0, nrow = n_grid, ncol = n_grid)
for (ii in seq_along(indices)) {
  top_young <- rank_young[1:indices[ii]]
  for (jj in seq_along(indices)) {
    top_old <- rank_old[1:indices[jj]]
    overlap <- length(intersect(top_young, top_old))
    # Two-sided: take minimum of upper and lower tail
    p_upper <- phyper(overlap - 1, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = FALSE)
    p_lower <- phyper(overlap, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = TRUE)
    p_val <- 2 * min(p_upper, p_lower)
    p_val <- min(p_val, 1)
    # Sign convention: positive = enriched (more overlap than expected)
    expected <- indices[ii] * indices[jj] / n_shared
    sign_factor <- ifelse(overlap >= expected, 1, -1)
    hmat[ii, jj] <- sign_factor * (-log10(max(p_val, .Machine$double.xmin)))
  }
}
message("  Hypergeometric matrix computed")

nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

# Quadrant max values
max_UU <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_DD <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_UD <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_DU <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

# Quadrant counts (at the exact midpoint split)
mid_young <- rank_young[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
n_UU <- length(intersect(mid_young, mid_old))
bot_young <- rank_young[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_DD <- length(intersect(bot_young, bot_old))
n_UD <- length(intersect(mid_young, bot_old))
n_DU <- length(intersect(bot_young, mid_old))

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

pD <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = expression(-log[10](P)),
                        guide = guide_colorbar(barwidth = unit(3, "cm"),
                                               barheight = unit(0.3, "cm"),
                                               title.position = "left",
                                               title.theme = element_text(size = 5.5, vjust = 0.8))) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = sprintf("Concordant Up\nmax = %.1f\nn = %d", max_UU, n_UU),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Concordant Down\nmax = %.1f\nn = %d", max_DD, n_DD),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Discordant\nY Up / O Down\nmax = %.1f\nn = %d", max_UD, n_UD),
           color = "white", fontface = "bold", size = 2.0) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Discordant\nY Down / O Up\nmax = %.1f\nn = %d", max_DU, n_DU),
           color = "white", fontface = "bold", size = 2.0) +
  # X-axis (Training Young) direction labels
  annotate("text", x = 1, y = -nc * 0.04,
           label = "<- Most upregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = nr, y = -nc * 0.04,
           label = "Most downregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  # Y-axis (Training Old) direction labels
  annotate("text", x = -nr * 0.04, y = 1, angle = 90,
           label = "<- Most upregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -nr * 0.04, y = nc, angle = 90,
           label = "Most downregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Threshold-Free Concordance (RRHO)",
    subtitle = sprintf("Two-sided hypergeometric overlap | %d shared genes | step = %d",
                        n_shared, step),
    x = "Training (Young) rank",
    y = "Training (Old) rank"
  ) +
  THEME_PUB +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(RPT_DIR, "panel_D_rrho2.pdf"), pD,
       width = 180, height = 180, units = "mm", device = pdf)

# Summary CSV
rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_overlap = c(n_UU, n_DD, n_UD, n_DU),
  matrix_rows = nr,
  matrix_cols = nc,
  n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"))

# Export full matrix
rrho2_export <- as.data.frame(hmat)
colnames(rrho2_export) <- paste0("col_", 1:nc)
rrho2_export$row <- 1:nr
rrho2_export <- rrho2_export %>% dplyr::select(row, everything())
write_csv(rrho2_export, file.path(DAT_DIR, "panel_D", "rrho2_matrix.csv"))

message("  Panel D saved")
