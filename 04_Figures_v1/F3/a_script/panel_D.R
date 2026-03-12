################################################################################
#   Figure 3 — Panel D: RRHO Reversal Map
#   Threshold-free overlap: t_Aging vs t_Training_Old
#   Pure-R implementation (Plaisier et al. 2010, NAR 38:e169) to avoid
#   RedRibbon C stack overflow.
#
#   Generates: panel_D_rrho2.pdf
#              + c_data/panel_D/rrho2_summary.csv, c_data/panel_D/rrho2_matrix.csv
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Hypergeometric test sidedness:
#    - Two-sided test used: p = 2 * min(p_upper, p_lower), capped at 1.  PASS
#    - Consistent with F2 Panel D (concordance map).                     PASS
#    - Two-sided is appropriate here because both enrichment (more overlap
#      than expected) and depletion (less than expected) are biologically
#      meaningful — enrichment = reversal or exacerbation depending on
#      quadrant, depletion = random/no systematic relationship.          PASS
#
# 2. Sign convention:
#    - Signed -log10(p): positive when overlap >= expected (enriched),
#      negative when overlap < expected (depleted). Matches Plaisier
#      et al. original convention.                                       PASS
#
# 3. Rank orientation:
#    - Aging ranks: index 1 = most upregulated (highest t-stat).         PASS
#    - Training_Old ranks: index 1 = most DOWNregulated (lowest t-stat). PASS
#    - This arrangement places Aging Up / Training Down (= Reversed)
#      in the top-left quadrant, the standard RRHO reversal layout.      PASS
#
# 4. Grid resolution:
#    - step = floor(n_shared / 200) gives ~200 x 200 grid.              PASS
#    - Fine enough to capture local hotspots without excessive compute.  PASS
#
# 5. Multiple testing:
#    - RRHO is a visualization tool, not a formal hypothesis test.
#      Individual cell p-values are not corrected for multiple
#      comparisons; the matrix is interpreted qualitatively by
#      hotspot location and magnitude. This is standard practice
#      (Plaisier et al. 2010).                                           PASS
#
# 6. No CIs needed: RRHO is a descriptive heatmap, not a point estimate.
#    Quadrant max -log10(p) values are reported for comparison.          PASS
# ---------------------------------------------------------------------------

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

message("Panel D: RRHO reversal map...")

# Prepare ranked lists sorted from most upregulated to most downregulated
rr_df <- dep_df %>%
  transmute(gene, t_aging = t_Aging, t_old = t_Training_Old) %>%
  filter(!is.na(t_aging) & !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)
message(sprintf("  %d shared genes for RRHO", n_shared))

# Rank lists:
# Aging:  index 1 = most upregulated (highest t-stat)
# Training_Old: index 1 = most DOWNregulated (lowest t-stat)
#   → top-left = Aging Up ∩ Training Down = Reversed
rank_aging <- rr_df$gene[order(-rr_df$t_aging)]
rank_old   <- rr_df$gene[order(rr_df$t_old)]

# Step size: ~200 grid points per axis
step <- max(1, floor(n_shared / 200))
indices <- seq(1, n_shared, by = step)
if (tail(indices, 1) != n_shared) indices <- c(indices, n_shared)
n_grid <- length(indices)
message(sprintf("  Grid: %d x %d (step = %d)", n_grid, n_grid, step))

# Compute -log10(p) hypergeometric overlap matrix
hmat <- matrix(0, nrow = n_grid, ncol = n_grid)
for (ii in seq_along(indices)) {
  top_aging <- rank_aging[1:indices[ii]]
  for (jj in seq_along(indices)) {
    top_old <- rank_old[1:indices[jj]]
    overlap <- length(intersect(top_aging, top_old))
    p_upper <- phyper(overlap - 1, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = FALSE)
    p_lower <- phyper(overlap, indices[ii], n_shared - indices[ii],
                      indices[jj], lower.tail = TRUE)
    p_val <- 2 * min(p_upper, p_lower)
    p_val <- min(p_val, 1)
    expected <- indices[ii] * indices[jj] / n_shared
    sign_factor <- ifelse(overlap >= expected, 1, -1)
    hmat[ii, jj] <- sign_factor * (-log10(max(p_val, .Machine$double.xmin)))
  }
}
message("  Hypergeometric matrix computed")

nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

# Quadrant max values (Training_Old axis flipped: col 1 = most downregulated)
# In ggplot: x = row (Aging, left=Up), y = col (Training_Old, bottom=Down)
# Bottom-left (small x, small y): Aging Up ∩ Training Down = Reversed (AgUp/TrDn)
# Top-right (large x, large y): Aging Down ∩ Training Up = Reversed (AgDn/TrUp)
# Top-left (small x, large y): Aging Up ∩ Training Up = Exacerbated Up
# Bottom-right (large x, small y): Aging Down ∩ Training Down = Exacerbated Down
max_rev_AgUp_TrDn <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_rev_AgDn_TrUp <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_exac_Up       <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_exac_Down     <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

# Quadrant counts at midpoint
mid_aging <- rank_aging[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
n_rev_AgUp_TrDn <- length(intersect(mid_aging, mid_old))
bot_aging <- rank_aging[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_rev_AgDn_TrUp <- length(intersect(bot_aging, bot_old))
n_exac_Up       <- length(intersect(mid_aging, bot_old))
n_exac_Down     <- length(intersect(bot_aging, mid_old))

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
           label = sprintf("Reversed\nAging Up / Training Down\nmax = %.1f\nn = %d", max_rev_AgUp_TrDn, n_rev_AgUp_TrDn),
           color = "white", fontface = "bold", size = 2.0) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Reversed\nAging Down / Training Up\nmax = %.1f\nn = %d", max_rev_AgDn_TrUp, n_rev_AgDn_TrUp),
           color = "white", fontface = "bold", size = 2.0) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Exacerbated Up\nmax = %.1f\nn = %d", max_exac_Up, n_exac_Up),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Exacerbated Down\nmax = %.1f\nn = %d", max_exac_Down, n_exac_Down),
           color = "white", fontface = "bold", size = 2.5) +
  # X-axis (Aging) direction labels
  annotate("text", x = 1, y = -nc * 0.04,
           label = "<- Most upregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = nr, y = -nc * 0.04,
           label = "Most downregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  # Y-axis (Training Old) direction labels (flipped: bottom = most downregulated)
  annotate("text", x = -nr * 0.04, y = 1, angle = 90,
           label = "<- Most downregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -nr * 0.04, y = nc, angle = 90,
           label = "Most upregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Threshold-Free Reversal (RRHO)",
    subtitle = sprintf("Two-sided hypergeometric overlap | %d shared genes | step = %d",
                        n_shared, step),
    x = "Aging rank",
    y = "Training (Old) rank"
  ) +
  THEME_PUB +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(RPT_DIR, "panel_D_rrho2.pdf"), pD,
       width = 180, height = 180, units = "mm", device = pdf)

# Summary CSV
rrho2_meta <- tibble(
  quadrant = c("Reversed_AgingUp_TrainingDown", "Reversed_AgingDown_TrainingUp",
               "Exacerbated_Up", "Exacerbated_Down"),
  max_neg_log10_pvalue = round(c(max_rev_AgUp_TrDn, max_rev_AgDn_TrUp,
                                  max_exac_Up, max_exac_Down), 2),
  n_overlap = c(n_rev_AgUp_TrDn, n_rev_AgDn_TrUp, n_exac_Up, n_exac_Down),
  matrix_rows = nr,
  matrix_cols = nc,
  n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"))

# Full matrix export
rrho2_export <- as.data.frame(hmat)
colnames(rrho2_export) <- paste0("col_", 1:nc)
rrho2_export$row <- 1:nr
rrho2_export <- rrho2_export %>% dplyr::select(row, everything())
write_csv(rrho2_export, file.path(DAT_DIR, "panel_D", "rrho2_matrix.csv"))

message("  Panel D saved")
