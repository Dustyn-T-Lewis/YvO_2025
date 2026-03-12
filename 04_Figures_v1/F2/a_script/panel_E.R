################################################################################
#   Figure 2 — Panel E: fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP)
#   Generates: panel_E_nes_bubble.pdf, panel_E_nes_bubble.png
#              + c_data/panel_E/nes_scatter.csv, nes_scatter_stats.csv
################################################################################
#
# ── STAT AUDIT (Task 13, 2026-02-27) ─────────────────────────────────────────
# 1. Test appropriateness:
#    - Pearson r on NES values tests linear concordance of pathway-level
#      enrichment between Young and Old training responses. NES values are
#      continuous and approximately normally distributed (Gaussian by CLT
#      since NES is a normalized sum over ranked genes).
# 2. Assumption checking:
#    - NES is mean-zero under the null, approximately Gaussian. Pearson r is
#      appropriate. With N=30-80 pathways, the correlation test has moderate
#      power; effect-size interpretation is primary.
# 3. Weighting by set size:
#    - Currently UNWEIGHTED. Weighting by set size would up-weight broad
#      pathways (100+ genes) which are more precisely estimated. However,
#      unweighted is standard in GSEA concordance analyses (e.g., Subramanian
#      et al. 2005) and avoids giving disproportionate influence to a few
#      large gene sets. A weighted Pearson r is computed and exported as a
#      sensitivity check but not used for the primary display.
# 4. Effect sizes: r IS the effect size. Now reported with 95% CI.
# 5. CIs: ADDED — 95% CI via cor.test() for unweighted Pearson r.
#    Weighted r CI via Fisher z-transformation.
# 6. Multiple comparisons: Single correlation test; no MTC needed.
# 7. Reproducibility: Deterministic (NES from fGSEA; correlation is exact).
# ─────────────────────────────────────────────────────────────────────────────

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_F2_setup.R")

# ==============================================================================
# PANEL E — fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP only)
# ==============================================================================

message("Panel E: NES scatter (Hallmark + GO:BP reduced)...")

# --- 1. Filter to Hallmark + GO:BP only ---
fgsea_hbp <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO:BP"),
         contrast %in% c("Training_Young", "Training_Old", "Interaction"))

# --- 2. Reduce GO:BP terms with rrvgo (threshold 0.5) ---
gobp_sig_names <- fgsea_hbp %>%
  filter(database == "GO:BP", padj < 0.05) %>%
  pull(pathway) %>% unique()

gobp_keep <- gobp_sig_names  # fallback: keep all

if (length(gobp_sig_names) > 5) {
  gobp_msigdb <- msigdbr(species = "Homo sapiens", collection = "C5",
                          subcollection = "GO:BP") %>%
    dplyr::select(gs_name, gs_exact_source) %>% distinct()
  name_to_goid <- setNames(gobp_msigdb$gs_exact_source, gobp_msigdb$gs_name)
  gobp_go_ids <- unique(na.omit(name_to_goid[gobp_sig_names]))

  if (length(gobp_go_ids) > 5) {
    tryCatch({
      hsGO <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")
      simMatrix <- calculateSimMatrix(gobp_go_ids, orgdb = "org.Hs.eg.db",
                                       semdata = hsGO,
                                       ont = "BP", method = "Rel")
      # Build scores vector (mean -log10 padj across contrasts for each GO term)
      gobp_scores_df <- fgsea_hbp %>%
        filter(database == "GO:BP", padj < 0.05) %>%
        dplyr::select(pathway, padj) %>%
        left_join(gobp_msigdb %>% rename(pathway = gs_name), by = "pathway") %>%
        filter(gs_exact_source %in% gobp_go_ids) %>%
        group_by(gs_exact_source) %>%
        summarise(score = mean(-log10(padj)), .groups = "drop")
      scores_vec <- setNames(gobp_scores_df$score, gobp_scores_df$gs_exact_source)
      reducedTerms <- reduceSimMatrix(simMatrix, scores = scores_vec,
                                       threshold = 0.5,
                                       orgdb = "org.Hs.eg.db")
      goid_to_name <- setNames(gobp_msigdb$gs_name, gobp_msigdb$gs_exact_source)
      # Keep only parent (representative) terms, not all children
      parent_go_ids <- unique(reducedTerms$parent)
      gobp_keep <- unique(na.omit(goid_to_name[parent_go_ids]))
      message(sprintf("  rrvgo reduced GO:BP from %d to %d terms",
                      length(gobp_sig_names), length(gobp_keep)))
    }, error = function(e) {
      message("  rrvgo failed: ", e$message, " — keeping all GO:BP")
    })
  }
}

# Keep Hallmark pathways + reduced GO:BP pathways
fgsea_filtered <- fgsea_hbp %>%
  filter(database == "Hallmark" | pathway %in% gobp_keep)

# --- 3. Pivot wide and classify ---
fgsea_wide <- fgsea_filtered %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Training_Young), !is.na(NES_Training_Old)) %>%
  mutate(set_size = coalesce(size_Training_Young, size_Training_Old))

fgsea_sig <- fgsea_wide %>%
  filter(
    (!is.na(padj_Training_Young) & padj_Training_Young < 0.05) |
    (!is.na(padj_Training_Old)   & padj_Training_Old < 0.05) |
    (!is.na(padj_Interaction)    & padj_Interaction < 0.05)
  ) %>%
  mutate(
    sig_Y = !is.na(padj_Training_Young) & padj_Training_Young < 0.05,
    sig_O = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    sig_I = !is.na(padj_Interaction) & padj_Interaction < 0.05,
    significance = case_when(
      sig_I         ~ "Interaction",
      sig_Y & sig_O ~ "Sig Both",
      sig_Y         ~ "Sig Young only",
      sig_O         ~ "Sig Old only",
      TRUE          ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS)),
    pathway_label = clean_pathway_name(pathway)
  )

message(sprintf("  %d pathways after filtering (Hallmark: %d, GO:BP: %d)",
                nrow(fgsea_sig),
                sum(fgsea_sig$database == "Hallmark"),
                sum(fgsea_sig$database == "GO:BP")))

nes_cor <- cor.test(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old,
                    conf.level = 0.95)
# Spearman rho for NES concordance (primary metric, matches manuscript)
nes_rho <- cor.test(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old,
                    method = "spearman", conf.level = 0.95)
# Spearman CI via Fisher z-transformation
n_pw_rho <- nrow(fgsea_sig)
rho_z_nes <- atanh(nes_rho$estimate)
rho_se_nes <- 1 / sqrt(n_pw_rho - 3)
rho_ci_nes <- tanh(rho_z_nes + c(-1, 1) * qnorm(0.975) * rho_se_nes)

nes_lim <- max(abs(c(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old))) * 1.15

# Weighted Pearson r (sensitivity check: weight by set_size)
if (requireNamespace("weights", quietly = TRUE)) {
  nes_r_wt <- weights::wtd.cor(fgsea_sig$NES_Training_Young,
                                fgsea_sig$NES_Training_Old,
                                weight = fgsea_sig$set_size)[1, "correlation"]
} else {
  # Manual weighted correlation
  w <- fgsea_sig$set_size / sum(fgsea_sig$set_size)
  mx <- sum(w * fgsea_sig$NES_Training_Young)
  my <- sum(w * fgsea_sig$NES_Training_Old)
  cov_xy <- sum(w * (fgsea_sig$NES_Training_Young - mx) * (fgsea_sig$NES_Training_Old - my))
  sd_x <- sqrt(sum(w * (fgsea_sig$NES_Training_Young - mx)^2))
  sd_y <- sqrt(sum(w * (fgsea_sig$NES_Training_Old - my)^2))
  nes_r_wt <- cov_xy / (sd_x * sd_y)
}
# Weighted r CI via Fisher z-transformation
n_pw <- nrow(fgsea_sig)
wt_z <- atanh(nes_r_wt)
wt_se <- 1 / sqrt(n_pw - 3)
wt_ci <- tanh(wt_z + c(-1, 1) * qnorm(0.975) * wt_se)

# Export NES scatter statistics
nes_stats <- tibble(
  metric = c("Spearman_rho", "Pearson_r_unweighted", "Pearson_r_weighted_by_set_size"),
  estimate = c(nes_rho$estimate, nes_cor$estimate, nes_r_wt),
  ci_lower = c(rho_ci_nes[1], nes_cor$conf.int[1], wt_ci[1]),
  ci_upper = c(rho_ci_nes[2], nes_cor$conf.int[2], wt_ci[2]),
  p_value  = c(nes_rho$p.value, nes_cor$p.value, NA_real_),
  n_pathways = c(n_pw, n_pw, n_pw),
  note = c("95% CI via Fisher z-transformation (primary metric)",
           "95% CI from cor.test()",
           "95% CI via Fisher z; weight = set_size (sensitivity check)")
)
write_csv(nes_stats, file.path(DAT_DIR, "panel_E", "nes_scatter_stats.csv"))
message(sprintf("  NES Spearman rho = %.3f [%.3f, %.3f], Pearson r = %.3f [%.3f, %.3f], weighted r = %.3f [%.3f, %.3f]",
                nes_rho$estimate, rho_ci_nes[1], rho_ci_nes[2],
                nes_cor$estimate, nes_cor$conf.int[1], nes_cor$conf.int[2],
                nes_r_wt, wt_ci[1], wt_ci[2]))

# Quadrant counts
nq1 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old > 0)
nq2 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old > 0)
nq3 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old < 0)
nq4 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old < 0)

# Label pathways with set_size >= 50
label_pw <- fgsea_sig %>%
  filter(set_size >= 50)

# Bubble border: Hallmark = solid black, GO:BP = transparent
fgsea_sig <- fgsea_sig %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "NS"             ~ 0.45,
      significance == "Interaction"    ~ 0.60,
      significance == "Sig Both"       ~ 0.75,
      significance == "Sig Old only"   ~ 0.85,
      significance == "Sig Young only" ~ 0.85,
      TRUE ~ 0.60
    )
  )

label_pw <- label_pw %>%
  mutate(
    label_fill = SIG_LABEL_FILL[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT[as.character(significance)]
  ) %>%
  # Condense verbose pathway labels for readability
  mutate(pathway_label = pathway_label %>%
    str_replace("Neg(ative)? Reg(ulation)? Of Programmed Cell Death", "Anti-Apoptosis") %>%
    str_replace("Pos(itive)? Reg(ulation)? Of Cellular Component.*", "Pos. Reg. Cell Comp. Org.") %>%
    str_replace("Positive Regulation Of ", "Pos. Reg. ") %>%
    str_replace("Negative Regulation Of ", "Neg. Reg. ") %>%
    str_replace("Regulation Of ", "Reg. ") %>%
    str_replace("Process Utilizing Autophagic Me.*", "Autophagy") %>%
    str_replace("Post Translational Protein Modi.*", "Protein PTMs") %>%
    str_replace("Proton Motive Force Driven.*", "PMF-Driven ATP Synthesis") %>%
    str_replace("ATP Synthesis Coupled Electron.*", "ETC / ATP Synthesis") %>%
    str_replace("Ribose Phosphate Biosynthetic.*", "Ribose-P Biosynthesis") %>%
    str_replace("Purine Containing Compound Bio.*", "Purine Biosynthesis") %>%
    str_replace("Modified Amino Acid Metabolic.*", "Modified AA Metabolism") %>%
    str_replace("Sulfur Compound Metabolic.*", "Sulfur Metabolism") %>%
    str_replace("Proteolysis Involved In Protei.*", "Protein Proteolysis") %>%
    str_replace("Muscle Structure Development", "Muscle Development") %>%
    str_replace("Mrna Metabolic.*", "mRNA Metabolism") %>%
    str_replace("Small Molecule Catabolic.*", "Small Molecule Catabolism") %>%
    str_replace("Lipid Catabolic Process", "Lipid Catabolism") %>%
    str_replace("Microtubule Based Process", "Microtubule Processes") %>%
    str_replace("Proton Transmembrane Transport", "H+ Transmembrane Transport") %>%
    str_replace(" Process$", "")
  ) %>%
  # Per-label nudge: bias same-colored labels vertically so they cluster,

  # while still using a single ggrepel layer for cross-group repulsion
  mutate(nudge_y = case_when(
    significance == "Interaction"    ~ -0.12,
    significance == "Sig Young only" ~  0.12,
    significance == "Sig Both"       ~  0.18,
    significance == "Sig Old only"   ~ -0.18,
    TRUE ~ 0
  )) %>%
  arrange(significance)

# Sort for plotting (less prominent drawn first)
plot_df <- fgsea_sig %>%
  mutate(draw_order = factor(significance,
    levels = c("NS", "Interaction", "Sig Both", "Sig Old only", "Sig Young only"))) %>%
  arrange(draw_order)

pE <- ggplot(plot_df, aes(x = NES_Training_Young, y = NES_Training_Old)) +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Bubbles: fill = significance, size = set_size, border = database
  geom_point(aes(fill = significance, size = set_size),
             shape = 21,
             color = plot_df$border_col,
             alpha = plot_df$bubble_alpha,
             stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.8, color = "black"))) +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200),
                        guide = guide_legend(
                          order = 2,
                          override.aes = list(fill = "grey60",
                                              color = "black", alpha = 0.7))) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title = "Pathway-Level Concordance (fGSEA)",
    subtitle = sprintf("Hallmark + GO:BP (rrvgo-reduced) | padj < 0.05 | %d pathways | \u03c1 = %.3f [%.3f, %.3f], p %s",
                       nrow(fgsea_sig), nes_rho$estimate,
                       rho_ci_nes[1], rho_ci_nes[2],
                       ifelse(nes_rho$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", nes_rho$p.value))),
    x = "NES (Training Young)",
    y = "NES (Training Old)"
  ) +
  THEME_FIG +
  theme(legend.position = "none")

# Pathway labels — single layer (all labels repel each other to prevent
# cross-color overlap) with per-significance nudge_y for same-color clustering
pE <- pE +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = TXT_GENE, fontface = "bold",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = c(-nes_lim * 0.9, nes_lim * 0.9),
                   ylim = c(-nes_lim * 0.9, nes_lim * 0.9))

# Quadrant count labels (on top of pathway labels)
pE <- pE +
  annotate("label", x = nes_lim, y = nes_lim,
           label = sprintf("Concordant Up\u2002n = %d", nq1),
           hjust = 1, vjust = 1, size = TXT_QUADRANT, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = -nes_lim,
           label = sprintf("Concordant Down\u2002n = %d", nq3),
           hjust = 0, vjust = 0, size = TXT_QUADRANT, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = nes_lim,
           label = sprintf("Discordant\u2002n = %d", nq2),
           hjust = 0, vjust = 1, size = TXT_QUADRANT, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = nes_lim, y = -nes_lim,
           label = sprintf("Discordant\u2002n = %d", nq4),
           hjust = 1, vjust = 0, size = TXT_QUADRANT, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt"))

ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.pdf"), pE,
       width = 200, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.png"), pE,
       width = 200, height = 200, units = "mm", dpi = 300)

# Clean CSV
fgsea_sig %>%
  transmute(
    pathway,
    pathway_label,
    database,
    NES_Training_Young = round(NES_Training_Young, 3),
    NES_Training_Old   = round(NES_Training_Old, 3),
    padj_Training_Young = signif(padj_Training_Young, 4),
    padj_Training_Old   = signif(padj_Training_Old, 4),
    padj_Interaction    = signif(padj_Interaction, 4),
    significance        = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Training_Young) + abs(NES_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_E", "nes_scatter.csv"))

message("  Panel E saved")
