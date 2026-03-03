################################################################################
#   Figure 2 — Panel E: RRHO2 Concordance Map + Flanking ORA Bars
#   Rank-Rank Hypergeometric Overlap (Plaisier et al. 2010, NAR 38:e169)
#   Pure-R implementation with per-quadrant ORA pathway enrichment.
#
#   Generates: panel_E_rrho2.pdf
#              + c_data/panel_E/rrho2_summary.csv
#              + c_data/panel_E/rrho2_ora_concordant.csv
#              + c_data/panel_E/rrho2_ora_discordant.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")

message("Panel E: RRHO concordance map...")

# === 1. RANK GENES BY T-STAT FOR BOTH CONTRASTS ==============================

rr_df <- dep_df %>%
  transmute(gene, t_young = t_Training_Young, t_old = t_Training_Old) %>%
  filter(!is.na(t_young), !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared   <- nrow(rr_df)
rank_young <- rr_df$gene[order(-rr_df$t_young)]
rank_old   <- rr_df$gene[order(-rr_df$t_old)]

# === 2. BUILD HYPERGEOMETRIC GRID (~200x200) ==================================

step    <- max(1L, floor(n_shared / 200))
indices <- seq(1L, n_shared, by = step)
if (tail(indices, 1) != n_shared) indices <- c(indices, n_shared)
n_grid  <- length(indices)

hmat <- matrix(0, nrow = n_grid, ncol = n_grid)
for (ii in seq_along(indices)) {
  top_young <- rank_young[1:indices[ii]]
  for (jj in seq_along(indices)) {
    top_old  <- rank_old[1:indices[jj]]
    overlap  <- length(intersect(top_young, top_old))
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

# === 3. IDENTIFY QUADRANT PEAKS AND HOTSPOT GENES ============================

max_UU <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_DD <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_UD <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_DU <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

# Midpoint overlap counts
mid_young <- rank_young[1:indices[mid_r]]
mid_old   <- rank_old[1:indices[mid_c]]
bot_young <- rank_young[(indices[mid_r]+1):n_shared]
bot_old   <- rank_old[(indices[mid_c]+1):n_shared]
n_UU <- length(intersect(mid_young, mid_old))
n_DD <- length(intersect(bot_young, bot_old))
n_UD <- length(intersect(mid_young, bot_old))
n_DU <- length(intersect(bot_young, mid_old))

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
  DD = intersect(rank_young[indices[peak_DD$i]:n_shared],
                 rank_old[indices[peak_DD$j]:n_shared]),
  UD = intersect(rank_young[1:indices[peak_UD$i]],
                 rank_old[indices[peak_UD$j]:n_shared]),
  DU = intersect(rank_young[indices[peak_DU$i]:n_shared],
                 rank_old[1:indices[peak_DU$j]])
)

# === 4. RUN ORA ON CONCORDANT AND DISCORDANT GENE SETS =======================

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
             minGSSize = 3, maxGSSize = 500),
    error = function(e) NULL)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble())
  as.data.frame(res) %>%
    mutate(quadrant = label, pathway_label = clean_pathway_name(ID)) %>%
    arrange(pvalue) %>%
    slice_head(n = 5)
}

ora_concordant <- bind_rows(
  run_ora(hotspot_genes$UU, "Concordant Up"),
  run_ora(hotspot_genes$DD, "Concordant Down"))
ora_discordant <- bind_rows(
  run_ora(hotspot_genes$UD, "Discordant (Y Up / O Down)"),
  run_ora(hotspot_genes$DU, "Discordant (Y Down / O Up)"))

write_csv(ora_concordant, file.path(DAT_DIR, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_discordant, file.path(DAT_DIR, "panel_E", "rrho2_ora_discordant.csv"))

# === 5. BUILD RRHO HEATMAP + ORA BAR FLANKS ==================================

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "viridis", name = expression(-log[10](P)),
    guide = guide_colorbar(barwidth = unit(3, "cm"), barheight = unit(0.3, "cm"),
                           title.position = "left",
                           title.theme = element_text(size = 9, face = "bold",
                                                      vjust = 0.8))) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.5) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed",
             color = "white", linewidth = 0.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = sprintf("Concordant Up\nmax = %.1f\nn = %d", max_UU, n_UU),
           color = "white", fontface = "bold", size = 3.5) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5,
           y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Concordant Down\nmax = %.1f\nn = %d", max_DD, n_DD),
           color = "white", fontface = "bold", size = 3.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Discordant\nY Up / O Down\nmax = %.1f\nn = %d",
                           max_UD, n_UD),
           color = "white", fontface = "bold", size = 3.5) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Discordant\nY Down / O Up\nmax = %.1f\nn = %d",
                           max_DU, n_DU),
           color = "white", fontface = "bold", size = 3.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(title = "Threshold-Free Concordance (RRHO)",
       subtitle = sprintf("Two-sided hypergeometric | %d shared genes | step = %d",
                          n_shared, step),
       x = "Training (Young) rank  \u2190 Up | Down \u2192",
       y = "Training (Old) rank  \u2190 Up | Down \u2192") +
  FIG_THEME +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(2, 2, 2, 2, "mm"))

# --- ORA bar builder ---
build_ora_bars <- function(ora_df, bar_fill, empty_label) {
  if (nrow(ora_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = empty_label,
                 size = 3, color = "grey50", fontface = "italic") +
        theme_void() + theme(plot.margin = margin(2, 2, 2, 2, "mm")))
  }
  plot_df <- ora_df %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           pathway_label = str_trunc(pathway_label, 40, ellipsis = "..."),
           pathway_label = factor(pathway_label, levels = rev(unique(pathway_label)))) %>%
    group_by(quadrant) %>% slice_head(n = 3) %>% ungroup()
  ggplot(plot_df, aes(x = neg_log10_padj, y = pathway_label)) +
    geom_col(fill = bar_fill, color = "grey30", linewidth = 0.2, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", neg_log10_padj)),
              hjust = -0.15, size = 2.5, fontface = "bold", color = "grey30") +
    facet_wrap(~ quadrant, ncol = 1, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          strip.text  = element_text(size = 7, face = "bold"),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin = margin(2, 4, 2, 2, "mm"))
}

conc_empty <- sprintf("No enrichment at padj < 0.05\n(UU: %d, DD: %d hotspot genes)",
                      length(hotspot_genes$UU), length(hotspot_genes$DD))
disc_empty <- sprintf("No enrichment at padj < 0.05\n(UD: %d, DU: %d hotspot genes)",
                      length(hotspot_genes$UD), length(hotspot_genes$DU))
pE_conc <- build_ora_bars(ora_concordant, ORA_CONCORDANT_COL, conc_empty)
pE_disc <- build_ora_bars(ora_discordant, ORA_DISCORDANT_COL, disc_empty)

# === 6. ASSEMBLE VIA PATCHWORK ================================================

pE <- pE_heat / (pE_conc | pE_disc) +
  plot_layout(heights = c(4, 1))

# === 7. SAVE PDF + CSVs ======================================================

ggsave(file.path(RPT_DIR, "panel_E_rrho2.pdf"), pE,
       width = 220, height = 260, units = "mm", device = cairo_pdf)

# Count ORA pathways per quadrant from already-computed results
n_ora_UU <- sum(ora_concordant$quadrant == "Concordant Up")
n_ora_DD <- sum(ora_concordant$quadrant == "Concordant Down")
n_ora_UD <- sum(ora_discordant$quadrant == "Discordant (Y Up / O Down)")
n_ora_DU <- sum(ora_discordant$quadrant == "Discordant (Y Down / O Up)")

tibble(
  quadrant       = c("Concordant_Up", "Concordant_Down",
                     "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10p = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_overlap      = c(n_UU, n_DD, n_UD, n_DU),
  n_hotspot_genes = c(length(hotspot_genes$UU), length(hotspot_genes$DD),
                      length(hotspot_genes$UD), length(hotspot_genes$DU)),
  n_ora_pathways  = c(n_ora_UU, n_ora_DD, n_ora_UD, n_ora_DU),
  n_shared_genes  = n_shared
) %>%
  write_csv(file.path(DAT_DIR, "panel_E", "rrho2_summary.csv"))

message("  Panel E saved (RRHO + ORA flanks)")
