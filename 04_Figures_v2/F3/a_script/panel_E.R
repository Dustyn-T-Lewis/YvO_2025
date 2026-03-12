# ── F3 Panel E: fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP) ────────────
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(msigdbr)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
})

PW <- 200; PH <- 200
RPT <- "04_Figures_v2/F3/b_reports"
DAT <- "04_Figures_v2/F3/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)
  f2_cache <- "04_Figures_v2/F2/c_data/shared/fgsea_tstat_all_v2.csv"
  f1_cache <- "04_Figures_v2/F1/c_data/06_panel_F_fgsea_results.csv"
  if (file.exists(f2_cache)) {
    file.copy(f2_cache, fgsea_cache)
    message("  Copied fGSEA cache from F2")
  } else if (file.exists(f1_cache)) {
    file.copy(f1_cache, fgsea_cache)
    message("  Copied fGSEA cache from F1")
  } else {
    stop("fGSEA cache not found")
  }
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

pdf_device <- get_pdf_device()
message("Panel E: NES scatter (Hallmark + GO:BP reduced)...")

PE_W <- 200  # panel width mm

fgsea_hbp <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO:BP"),
         contrast %in% c("Aging", "Training_Old"))

gobp_sig_names <- fgsea_hbp %>%
  filter(database == "GO:BP", padj < 0.05) %>%
  pull(pathway) %>% unique()

gobp_keep <- gobp_sig_names

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
      parent_go_ids <- unique(reducedTerms$parent)
      gobp_keep <- unique(na.omit(goid_to_name[parent_go_ids]))
      message(sprintf("  rrvgo reduced GO:BP from %d to %d terms",
                      length(gobp_sig_names), length(gobp_keep)))
    }, error = function(e) {
      message("  rrvgo failed: ", e$message, " -- keeping all GO:BP")
    })
  }
}

fgsea_filtered <- fgsea_hbp %>%
  filter(database == "Hallmark" | pathway %in% gobp_keep)

fgsea_wide <- fgsea_filtered %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Aging), !is.na(NES_Training_Old)) %>%
  mutate(set_size = coalesce(size_Aging, size_Training_Old))

fgsea_sig <- fgsea_wide %>%
  filter(
    (!is.na(padj_Aging)        & padj_Aging < 0.05) |
    (!is.na(padj_Training_Old) & padj_Training_Old < 0.05)
  ) %>%
  mutate(
    sig_A = !is.na(padj_Aging) & padj_Aging < 0.05,
    sig_T = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    significance = case_when(
      sig_A & sig_T ~ "Sig Both",
      sig_A         ~ "Sig Aging only",
      sig_T         ~ "Sig Training only",
      TRUE          ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F3)),
    pathway_label = clean_pathway_name(pathway)
  )

message(sprintf("  %d pathways after filtering (Hallmark: %d, GO:BP: %d)",
                nrow(fgsea_sig),
                sum(fgsea_sig$database == "Hallmark"),
                sum(fgsea_sig$database == "GO:BP")))

nes_cor <- cor.test(fgsea_sig$NES_Aging, fgsea_sig$NES_Training_Old)
nes_lim <- max(abs(c(fgsea_sig$NES_Aging, fgsea_sig$NES_Training_Old))) * 1.15

n_rev_tl  <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old > 0)
n_rev_br  <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old < 0)
n_exac_tr <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old > 0)
n_exac_bl <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old < 0)

n_rev_pw <- n_rev_tl + n_rev_br
n_total_pw <- nrow(fgsea_sig)
pw_rev_frac <- n_rev_pw / n_total_pw
pw_rev_binom <- binom.test(n_rev_pw, n_total_pw)
pw_rev_ci <- pw_rev_binom$conf.int * 100

message(sprintf("  NES correlation: r = %.3f [%.3f, %.3f], p = %.2g",
                nes_cor$estimate, nes_cor$conf.int[1], nes_cor$conf.int[2],
                nes_cor$p.value))
message(sprintf("  Pathway reversal: %d/%d (%.1f%%) [%.1f, %.1f]",
                n_rev_pw, n_total_pw, pw_rev_frac * 100,
                pw_rev_ci[1], pw_rev_ci[2]))

txt_gene <- scale_text(BASE_GENE, PE_W)
txt_quad <- scale_text(BASE_QUADRANT, PE_W)

label_pw <- fgsea_sig %>%
  filter(set_size >= 50) %>%
  mutate(
    label_fill = SIG_LABEL_FILL_F3[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F3[as.character(significance)]
  ) %>%
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
    str_replace("Organic Acid Catabolic.*", "Organic Acid Catabolism") %>%
    str_replace("Fatty Acid Catabolic.*", "Fatty Acid Catabolism") %>%
    str_replace("Membraneless Organelle Assembly", "Membraneless Org. Assembly") %>%
    str_replace("Muscle Cell Development", "Muscle Cell Dev.") %>%
    str_replace(" Process$", "")
  ) %>%
  mutate(nudge_y = case_when(
    significance == "Sig Both"          ~  0.15,
    significance == "Sig Aging only"    ~ -0.15,
    significance == "Sig Training only" ~  0.10,
    TRUE ~ 0
  )) %>%
  arrange(significance)

fgsea_sig <- fgsea_sig %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "NS"                ~ 0.45,
      significance == "Sig Both"          ~ 0.75,
      significance == "Sig Aging only"    ~ 0.85,
      significance == "Sig Training only" ~ 0.85,
      TRUE ~ 0.60
    )
  )

plot_df <- fgsea_sig %>%
  mutate(draw_order = factor(significance,
    levels = c("NS", "Sig Training only", "Sig Aging only", "Sig Both", "Reversal"))) %>%
  arrange(draw_order)

pE <- ggplot(plot_df, aes(x = NES_Aging, y = NES_Training_Old)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(aes(fill = significance, size = set_size),
             shape = 21, color = plot_df$border_col,
             alpha = plot_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = txt_gene, fontface = "bold",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Reversed  n = %d", n_rev_br),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Reversed  n = %d", n_rev_tl),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Exacerbated  n = %d", n_exac_tr),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Exacerbated  n = %d", n_exac_bl),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title = "Pathway-Level Reversal (fGSEA)",
    subtitle = sprintf("Hallmark + GO:BP (rrvgo-reduced) | padj < 0.05 | %d pathways\nr = %.2f [%.2f, %.2f], p %s | %.0f%% reversed",
                       nrow(fgsea_sig), nes_cor$estimate,
                       nes_cor$conf.int[1], nes_cor$conf.int[2],
                       ifelse(nes_cor$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", nes_cor$p.value)),
                       pw_rev_frac * 100),
    x = "NES (Aging)",
    y = "NES (Training Old)"
  ) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_E_nes_scatter.pdf"), pE,
       width = PE_W, height = PE_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_nes_scatter.png"), pE,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)

fgsea_sig %>%
  transmute(
    pathway, pathway_label, database,
    NES_Aging        = round(NES_Aging, 3),
    NES_Training_Old = round(NES_Training_Old, 3),
    padj_Aging        = signif(padj_Aging, 4),
    padj_Training_Old = signif(padj_Training_Old, 4),
    significance      = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Aging) + abs(NES_Training_Old))) %>%
  write_csv(file.path(DAT, "panel_E", "nes_scatter.csv"))

cat("F3 Panel E done\n")
