################################################################################
#   Figure 3 — Panel D: RRHO2 Reversal Map + Flanking ORA Bars
#   Rank-Rank Hypergeometric Overlap (Plaisier et al. 2010, NAR 38:e169)
#   Pure-R implementation with per-quadrant ORA pathway enrichment.
#
#   Rank orientation (reversal layout):
#     Aging:        index 1 = highest t-stat  (most UP with aging)
#     Training_Old: index 1 = lowest  t-stat  (most DOWN with training)
#   Top-left = Aging Up / Training Down = REVERSED
#
#   Generates: panel_D_rrho2.pdf
#              + c_data/panel_D/rrho2_summary.csv
#              + c_data/panel_D/rrho2_ora_concordant.csv   (reversal pathways)
#              + c_data/panel_D/rrho2_ora_discordant.csv   (exacerbation pathways)
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

message("Panel D: RRHO reversal map + ORA bars...")

# 4-color ORA palette: one per RRHO quadrant (reversal framing)
# Local definition to ensure quadrant labels match exactly
ORA_QUAD_COLORS <- c(
  "Reversed (Aging Up / Training Down)"  = "#E57373",
  "Reversed (Aging Down / Training Up)"  = "#64B5F6",
  "Exacerbated Up"                       = "#FFB74D",
  "Exacerbated Down"                     = "#81C784"
)

# === 1. RANK GENES BY T-STAT FOR BOTH CONTRASTS ==============================

rr_df <- dep_df %>%
  transmute(gene, t_aging = t_Aging, t_old = t_Training_Old) %>%
  filter(!is.na(t_aging), !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared   <- nrow(rr_df)
message(sprintf("  %d shared genes for RRHO", n_shared))
rank_aging <- rr_df$gene[order(-rr_df$t_aging)]    # index 1 = most UP with aging
rank_old   <- rr_df$gene[order(rr_df$t_old)]        # index 1 = most DOWN with training

# === 2. BUILD HYPERGEOMETRIC GRID (~200x200) ==================================

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
message("  Hypergeometric matrix computed")

nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

# === 3. IDENTIFY QUADRANT PEAKS AND HOTSPOT GENES ============================
# UU (small i, small j) = Aging Up / Training Down = REVERSED
# DD (large i, large j) = Aging Down / Training Up = REVERSED
# UD (small i, large j) = Aging Up / Training Up   = EXACERBATED Up
# DU (large i, small j) = Aging Down / Training Down = EXACERBATED Down

max_rev_AgUp_TrDn <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_rev_AgDn_TrUp <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_exac_Up       <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_exac_Down     <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

# Midpoint overlap counts
mid_aging <- rank_aging[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
bot_aging <- rank_aging[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_rev_AgUp_TrDn <- length(intersect(mid_aging, mid_old))
n_rev_AgDn_TrUp <- length(intersect(bot_aging, bot_old))
n_exac_Up       <- length(intersect(mid_aging, bot_old))
n_exac_Down     <- length(intersect(bot_aging, mid_old))

# Peak-based hotspot gene extraction
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

# === 4. RUN ORA ON REVERSAL AND EXACERBATION GENE SETS =======================

hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
gobp_t2g <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
all_t2g  <- bind_rows(hallmark_t2g, gobp_t2g) %>% distinct()
all_genes <- rr_df$gene

run_ora <- function(gene_set, label) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    enricher(gene = gene_set, universe = all_genes, TERM2GENE = all_t2g,
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1,
             minGSSize = 10, maxGSSize = 500),
    error = function(e) NULL)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble())
  as.data.frame(res) %>%
    mutate(quadrant = label, pathway_label = clean_pathway_name(ID)) %>%
    # Deduplicate: GO terms sharing identical gene sets produce identical
    # statistics. Keep the shortest-named term per unique gene set.
    arrange(pvalue, nchar(ID)) %>%
    distinct(geneID, .keep_all = TRUE) %>%
    slice_head(n = 5)
}

# Reversal = off-diagonal: UU (Aging Up / Training Down) + DD (Aging Down / Training Up)
ora_reversal <- bind_rows(
  run_ora(hotspot_genes$UU, "Reversed (Aging Up / Training Down)"),
  run_ora(hotspot_genes$DD, "Reversed (Aging Down / Training Up)"))

# Exacerbation = diagonal: UD (Aging Up / Training Up) + DU (Aging Down / Training Down)
ora_exacerbation <- bind_rows(
  run_ora(hotspot_genes$UD, "Exacerbated Up"),
  run_ora(hotspot_genes$DU, "Exacerbated Down"))

dir.create(file.path(DAT_DIR, "panel_D"), recursive = TRUE, showWarnings = FALSE)
write_csv(ora_reversal, file.path(DAT_DIR, "panel_D", "rrho2_ora_concordant.csv"))
write_csv(ora_exacerbation, file.path(DAT_DIR, "panel_D", "rrho2_ora_discordant.csv"))

# === 5. BUILD RRHO HEATMAP + ORA BAR FLANKS ==================================

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

pD_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "viridis", name = expression(-log[10](P)),
    guide = guide_colorbar(barwidth = unit(4, "cm"), barheight = unit(0.4, "cm"),
                           title.position = "left",
                           title.theme = element_text(size = 11, face = "bold",
                                                      vjust = 0.8))) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.5) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = sprintf("Reversed\nAging Up / Training Down\nmax = %.1f\nn = %d",
                           max_rev_AgUp_TrDn, n_rev_AgUp_TrDn),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5,
           y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Reversed\nAging Down / Training Up\nmax = %.1f\nn = %d",
                           max_rev_AgDn_TrUp, n_rev_AgDn_TrUp),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Exacerbated Up\nmax = %.1f\nn = %d",
                           max_exac_Up, n_exac_Up),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Exacerbated Down\nmax = %.1f\nn = %d",
                           max_exac_Down, n_exac_Down),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  # Axis direction labels
  annotate("text", x = 1, y = -nc * 0.04,
           label = "<- Most upregulated", hjust = 0, size = TXT_STATS,
           color = "grey30") +
  annotate("text", x = nr, y = -nc * 0.04,
           label = "Most downregulated ->", hjust = 1, size = TXT_STATS,
           color = "grey30") +
  annotate("text", x = -nr * 0.04, y = 1, angle = 90,
           label = "<- Most downregulated", hjust = 0, size = TXT_STATS,
           color = "grey30") +
  annotate("text", x = -nr * 0.04, y = nc, angle = 90,
           label = "Most upregulated ->", hjust = 1, size = TXT_STATS,
           color = "grey30") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off", expand = FALSE) +
  labs(title = "Threshold-Free Reversal (RRHO)",
       subtitle = sprintf("Two-sided hypergeometric | %d shared genes | step = %d",
                          n_shared, step),
       x = "Aging rank",
       y = "Training (Old) rank") +
  THEME_FIG +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(8, 8, 8, 8, "mm"))

# --- Unified ORA bar chart (single panel, 4-color per-quadrant) ---
ora_all <- bind_rows(ora_reversal, ora_exacerbation)

if (nrow(ora_all) > 0) {
  ora_plot_df <- ora_all %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           pathway_label = str_trunc(clean_pathway_name(ID), 40, ellipsis = "...")) %>%
    group_by(quadrant) %>% slice_head(n = 3) %>% ungroup() %>%
    mutate(quadrant = factor(quadrant, levels = names(ORA_QUAD_COLORS)))

  pD_ora <- ggplot(ora_plot_df,
                   aes(x = neg_log10_padj, y = reorder(pathway_label, neg_log10_padj),
                       fill = quadrant)) +
    geom_col(color = "grey30", linewidth = 0.2, width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", neg_log10_padj)),
              hjust = -0.15, size = TXT_ORA_BAR, fontface = "bold", color = "grey30") +
    scale_fill_manual(values = ORA_QUAD_COLORS, name = "Quadrant") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    facet_wrap(~quadrant, ncol = 1, scales = "free_y") +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    THEME_FIG +
    theme(axis.text.y = element_text(size = rel(0.7)),
          axis.text.x = element_text(size = rel(0.7)),
          legend.position = "none",
          strip.text = element_text(size = rel(0.7), face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.spacing = unit(2, "mm"),
          plot.margin = margin(2, 4, 2, 2, "mm"))
} else {
  pD_ora <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("No enrichment at padj < 0.05\n(%d shared genes)",
                             n_shared),
             size = TXT_LABEL, color = "grey50", fontface = "italic") +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm"))
}

# === 6. ASSEMBLE VIA PATCHWORK ================================================

pD <- (pD_heat | pD_ora) +
  plot_layout(widths = c(3, 2)) &
  theme(legend.position = "bottom")

# === 7. SAVE PDF + CSVs ======================================================

ggsave(file.path(RPT_DIR, "panel_D_rrho2.pdf"), pD,
       width = 350, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_D_rrho2.png"), pD,
       width = 350, height = 200, units = "mm", dpi = 300)

# Count ORA pathways per quadrant
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
  write_csv(file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"))

message("  Panel D saved (RRHO + ORA flanks)")
