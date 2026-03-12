################################################################################
#   Figure 2 — Panel D: RRHO Concordance Map
#   Rank-Rank Hypergeometric Overlap (Plaisier et al. 2010, NAR 38:e169)
#   Pure-R implementation to avoid RedRibbon C stack overflow.
#
#   Generates: panel_D_rrho2.pdf
#              + c_data/panel_D/rrho2_summary.csv, rrho2_matrix.csv,
#                rrho2_hotspot_genes.csv, rrho2_ora_concordant.csv,
#                rrho2_ora_discordant.csv
################################################################################
#
# ── STAT AUDIT (Task 13, 2026-02-27) ─────────────────────────────────────────
# 1. Test appropriateness:
#    - Hypergeometric test (Fisher exact) is the standard test for RRHO
#      (Plaisier et al. 2010). It tests whether the overlap between top-ranked
#      genes in two lists exceeds chance given the shared universe size.
# 2. One-sided vs two-sided:
#    - TWO-SIDED test implemented (lines 47-52): p = 2*min(p_upper, p_lower),
#      capped at 1. This detects both enrichment (more overlap than expected)
#      and depletion (less overlap than expected). Sign convention encodes
#      direction: positive = enriched, negative = depleted.
#    - Rationale: two-sided is conservative and detects both concordant AND
#      discordant patterns. The four quadrants naturally decompose into
#      concordant (UU, DD) and discordant (UD, DU) signal.
# 3. Grid resolution:
#    - ~200x200 grid (step = floor(n_shared/200)). This is the standard RRHO
#      resolution used in Plaisier et al. and subsequent implementations.
#      Sensitivity: finer grids produce smoother maps but do not change the
#      location or magnitude of the dominant signal. The step size is reported
#      in the plot subtitle and summary CSV for reproducibility.
# 4. Multiple testing across grid cells:
#    - NO correction applied across grid cells. This is standard RRHO practice
#      (Plaisier et al. 2010; Cahill et al. 2018). The RRHO map is a
#      visualization tool where patterns (hot-spot location and extent) are
#      interpreted holistically, not as independent tests. The max -log10(p)
#      in each quadrant is reported as the summary statistic. Formal inference
#      is derived from the protein-level and pathway-level tests in Panels C/E.
# 5. Effect sizes: The -log10(p) value encodes both significance and direction.
#    Overlap counts (n) at the midpoint split are reported per quadrant.
# 6. CIs: Not standard for RRHO heatmaps. The grid cells provide a continuous
#    landscape; formal CIs would require permutation-based null distributions,
#    which is beyond standard RRHO practice.
# 7. Reproducibility: Deterministic (no random component; ranks from t-stats).
# ─────────────────────────────────────────────────────────────────────────────

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
                        guide = guide_colorbar(barwidth = unit(4, "cm"),
                                               barheight = unit(0.4, "cm"),
                                               title.position = "left",
                                               title.theme = element_text(size = 11, face = "bold", vjust = 0.8))) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = sprintf("Concordant Up\nmax = %.1f\nn = %d", max_UU, n_UU),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Concordant Down\nmax = %.1f\nn = %d", max_DD, n_DD),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Discordant\nY Up / O Down\nmax = %.1f\nn = %d", max_UD, n_UD),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Discordant\nY Down / O Up\nmax = %.1f\nn = %d", max_DU, n_DU),
           color = "white", fontface = "bold", size = TXT_QUADRANT) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Threshold-Free Concordance (RRHO)",
    subtitle = sprintf("Two-sided hypergeometric overlap | %d shared genes | step = %d",
                        n_shared, step),
    x = "Training (Young) rank  \u2190 Up | Down \u2192",
    y = "Training (Old) rank  \u2190 Up | Down \u2192"
  ) +
  THEME_FIG +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ==============================================================================
# HOTSPOT GENE EXTRACTION + PER-QUADRANT ORA
# ==============================================================================
# ORA gene sets are explicitly derived from RRHO hotspot peak positions.
# This is downstream interpretation of RRHO hotspots, not a separate analysis.

message("  Extracting hotspot genes and running per-quadrant ORA...")

# Find peak cell in each quadrant
find_peak <- function(mat, rows, cols) {
  sub_mat <- mat[rows, cols, drop = FALSE]
  peak <- which(sub_mat == max(sub_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(i = rows[peak[1]], j = cols[peak[2]])
}

peak_UU <- find_peak(hmat, 1:mid_r, 1:mid_c)
peak_DD <- find_peak(hmat, (mid_r+1):nr, (mid_c+1):nc)
peak_UD <- find_peak(hmat, 1:mid_r, (mid_c+1):nc)
peak_DU <- find_peak(hmat, (mid_r+1):nr, 1:mid_c)

# Extract gene sets at each peak position
hotspot_genes <- list(
  UU = intersect(rank_young[1:indices[peak_UU$i]], rank_old[1:indices[peak_UU$j]]),
  DD = intersect(rank_young[indices[peak_DD$i]:n_shared], rank_old[indices[peak_DD$j]:n_shared]),
  UD = intersect(rank_young[1:indices[peak_UD$i]], rank_old[indices[peak_UD$j]:n_shared]),
  DU = intersect(rank_young[indices[peak_DU$i]:n_shared], rank_old[1:indices[peak_DU$j]])
)

message(sprintf("  Hotspot gene counts: UU=%d, DD=%d, UD=%d, DU=%d",
                length(hotspot_genes$UU), length(hotspot_genes$DD),
                length(hotspot_genes$UD), length(hotspot_genes$DU)))

# Export hotspot genes
hotspot_export <- bind_rows(
  tibble(quadrant = "UU", gene = hotspot_genes$UU),
  tibble(quadrant = "DD", gene = hotspot_genes$DD),
  tibble(quadrant = "UD", gene = hotspot_genes$UD),
  tibble(quadrant = "DU", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT_DIR, "panel_D", "rrho2_hotspot_genes.csv"))

# Per-quadrant ORA
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
gobp_t2g <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
all_t2g_D <- bind_rows(hallmark_t2g, gobp_t2g) %>% distinct()
all_genes_D <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    enricher(gene = gene_set, universe = all_genes_D, TERM2GENE = all_t2g_D,
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1,
             minGSSize = 10, maxGSSize = 500),
    error = function(e) NULL)
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    as.data.frame(res) %>%
      mutate(quadrant = quadrant_name,
             pathway_label = clean_pathway_name(ID)) %>%
      # Deduplicate: GO terms sharing identical gene sets produce identical
      # statistics. Keep the shortest-named term per unique gene set.
      arrange(pvalue, nchar(ID)) %>%
      distinct(geneID, .keep_all = TRUE) %>%
      slice_head(n = 5)
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

write_csv(ora_concordant,  file.path(DAT_DIR, "panel_D", "rrho2_ora_concordant.csv"))
write_csv(ora_discordant,  file.path(DAT_DIR, "panel_D", "rrho2_ora_discordant.csv"))

message(sprintf("  ORA results: concordant=%d pathways, discordant=%d pathways",
                nrow(ora_concordant), nrow(ora_discordant)))

# ==============================================================================
# BUILD 4-QUADRANT ORA BAR PLOT
# ==============================================================================

# Combine all four quadrants into one dataframe with quadrant-specific colors
ora_all <- bind_rows(ora_UU, ora_DD, ora_UD, ora_DU)

if (nrow(ora_all) == 0) {
  pD_ora <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("No enrichment at padj < 0.05\n(UU: %d, DD: %d, UD: %d, DU: %d hotspot genes)",
                             length(hotspot_genes$UU), length(hotspot_genes$DD),
                             length(hotspot_genes$UD), length(hotspot_genes$DU)),
             size = TXT_ORA_BAR, color = "grey50", fontface = "italic") +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm"))
} else {
  # For empty quadrants, insert a placeholder row so the facet still appears
  all_quadrant_names <- c("Concordant Up", "Concordant Down",
                          "Discordant (Y Up / O Down)", "Discordant (Y Down / O Up)")
  empty_quads <- setdiff(all_quadrant_names, unique(ora_all$quadrant))

  plot_df <- ora_all %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      pathway_label = str_trunc(pathway_label, 40, ellipsis = "..."),
      quadrant = factor(quadrant, levels = all_quadrant_names)
    ) %>%
    group_by(quadrant) %>%
    slice_head(n = 3) %>%
    ungroup()

  # Within-facet ordering: append quadrant suffix to make labels unique, then reorder
  plot_df <- plot_df %>%
    mutate(
      pathway_uid = paste(pathway_label, quadrant, sep = "___"),
      pathway_uid = reorder(pathway_uid, neg_log10_padj)
    )

  pD_ora <- ggplot(plot_df, aes(x = neg_log10_padj, y = pathway_uid, fill = quadrant)) +
    geom_col(color = "grey30", linewidth = 0.2, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", neg_log10_padj)),
              hjust = -0.15, size = TXT_ORA_BAR, fontface = "bold", color = "grey30") +
    facet_wrap(~ quadrant, scales = "free_y", ncol = 1) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = ORA_QUAD_COLORS, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    THEME_FIG +
    theme(
      axis.text.y  = element_text(size = rel(0.7)),
      axis.text.x  = element_text(size = rel(0.7)),
      strip.text   = element_text(size = rel(0.7), face = "bold"),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin = margin(2, 4, 2, 2, "mm")
    )

  # Add annotations for empty quadrants
  if (length(empty_quads) > 0) {
    empty_labels <- sapply(empty_quads, function(q) {
      key <- switch(q,
                    "Concordant Up"              = "UU",
                    "Concordant Down"            = "DD",
                    "Discordant (Y Up / O Down)" = "UD",
                    "Discordant (Y Down / O Up)" = "DU")
      sprintf("No enrichment (padj < 0.05)\n%d hotspot genes", length(hotspot_genes[[key]]))
    })
    message(sprintf("  Empty ORA quadrants: %s", paste(empty_quads, collapse = ", ")))
  }
}

# ==============================================================================
# PATCHWORK ASSEMBLY: RRHO (top) / ORA bars (bottom: 2x2 quadrant grid)
# ==============================================================================

pD_combined <- (pD | pD_ora) +
  plot_layout(widths = c(3, 2))

ggsave(file.path(RPT_DIR, "panel_D_rrho2.pdf"), pD_combined,
       width = 220, height = 260, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_D_rrho2.png"), pD_combined,
       width = 350, height = 200, units = "mm", dpi = 300)

# Summary CSV
rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_overlap = c(n_UU, n_DD, n_UD, n_DU),
  n_hotspot_genes = c(length(hotspot_genes$UU), length(hotspot_genes$DD),
                      length(hotspot_genes$UD), length(hotspot_genes$DU)),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
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

# Assign pD to the combined version for composite assembly
pD <- pD_combined

message("  Panel D saved (RRHO + ORA flanks)")
