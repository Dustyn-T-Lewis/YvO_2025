################################################################################
#   Figure 2 — Panel E: fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP)
#   Generates: panel_E_nes_bubble.pdf, panel_E_nes_bubble.png
#              + c_data/panel_E/nes_scatter.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_F2_setup.R")

# ==============================================================================
# PANEL E — fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP only)
# ==============================================================================

message("Panel E: NES scatter (Hallmark + GO:BP reduced)...")

# --- 1. Filter to Hallmark + GO:BP only ---
fgsea_hbp <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO:BP"),
         contrast %in% c("Training_Young", "Training_Old", "Interaction"))

# --- 2. Reduce GO:BP terms with rrvgo (threshold 0.85) ---
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

nes_cor <- cor.test(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old)
nes_lim <- max(abs(c(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old))) * 1.15

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
  coord_cartesian(xlim = c(-3.5, 2.5), ylim = c(-2.5, 3)) +
  labs(
    title = "Pathway-Level Concordance (fGSEA)",
    subtitle = sprintf("Hallmark + GO:BP (rrvgo-reduced) | padj < 0.05 | %d pathways | r = %.2f, p %s",
                       nrow(fgsea_sig), nes_cor$estimate,
                       ifelse(nes_cor$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", nes_cor$p.value))),
    x = "NES (Training Young)",
    y = "NES (Training Old)"
  ) +
  THEME_PUB +
  theme(legend.position = "none")

# Pathway labels — single layer (all labels repel each other to prevent
# cross-color overlap) with per-significance nudge_y for same-color clustering
pE <- pE +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = 2.2, fontface = "bold",
                   max.overlaps = 30,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = c(-3.2, 2.2),
                   ylim = c(-2.2, 2.7))

# Quadrant count labels (on top of pathway labels)
pE <- pE +
  annotate("label", x = 2.5, y = 3.0,
           label = sprintf("Concordant Up\u2002n = %d", nq1),
           hjust = 1, vjust = 1, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -3.5, y = -2.5,
           label = sprintf("Concordant Down\u2002n = %d", nq3),
           hjust = 0, vjust = 0, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -3.5, y = 3.0,
           label = sprintf("Discordant\u2002n = %d", nq2),
           hjust = 0, vjust = 1, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = 2.5, y = -2.5,
           label = sprintf("Discordant\u2002n = %d", nq4),
           hjust = 1, vjust = 0, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt"))

# --- Hand-built legend: three columns side-by-side, items stacked vertically ---

sig_levels_e <- c("Interaction", "Sig Both", "Sig Young only", "Sig Old only")
ks_e <- 0.15

# Column 1: Significance
sig_key_df <- tibble(
  x = 0, y = rev(seq_along(sig_levels_e)) * ks_e,
  label = sig_levels_e,
  fill  = unname(SIG_COLORS[sig_levels_e])
)
# Column 2: Set size
size_breaks_e <- c(20, 50, 100)
size_range_e  <- c(2, 8)
size_key_df <- tibble(
  x = 3.5, y = rev(seq_along(size_breaks_e)) * ks_e,
  label = as.character(size_breaks_e),
  pt_size = scales::rescale(size_breaks_e, to = size_range_e, from = c(20, 200))
)
# Column 3: Database
db_key_df <- tibble(
  x = 6.0, y = c(2, 1) * ks_e,
  label = c("Hallmark", "GO:BP"),
  border = c("black", "grey75"),
  stroke = c(0.8, 1.2)
)

title_y_e <- (max(length(sig_levels_e), length(size_breaks_e), 2) + 1) * ks_e

pE_key <- ggplot() +
  # Significance column
  annotate("text", x = 0, y = title_y_e, label = "Significance",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_point(data = sig_key_df, aes(x = x, y = y),
             shape = 21, size = 3.5, fill = sig_key_df$fill,
             color = "black", stroke = 0.8) +
  geom_text(data = sig_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = KEY_ITEM, color = KEY_ITEM_COL) +
  # Set size column
  annotate("text", x = 3.5, y = title_y_e, label = "Set size",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_point(data = size_key_df, aes(x = x, y = y),
             shape = 21, size = size_key_df$pt_size, fill = "grey60",
             color = "black", alpha = 0.7) +
  geom_text(data = size_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = KEY_ITEM, color = KEY_ITEM_COL) +
  # Database column
  annotate("text", x = 6.0, y = title_y_e, label = "Database",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_point(data = db_key_df, aes(x = x, y = y),
             shape = 21, size = 3.5, fill = "grey70",
             color = db_key_df$border, stroke = db_key_df$stroke) +
  geom_text(data = db_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = KEY_ITEM, color = KEY_ITEM_COL) +
  scale_x_continuous(limits = c(-0.3, 8.5)) +
  scale_y_continuous(limits = c(0, title_y_e + ks_e)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pE_combined <- pE / pE_key + plot_layout(heights = c(0.90, 0.10))

ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.pdf"), pE_combined,
       width = 200, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.png"), pE_combined,
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
