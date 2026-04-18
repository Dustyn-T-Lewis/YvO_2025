# F05 Panel B: fGSEA NES Scatter (GO Slim + Hallmark, a priori) — Reversal
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)
library(ggrepel)

PW <- 200; PH <- 200
RPT_PNG <- "04_Figures/F05/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F05/b_reports/main/pdf/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

fgsea_cache <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
stopifnot("fGSEA cache missing: frozen cache — see 04_Figures/shared/README.md" =
            file.exists(fgsea_cache))
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# Filter to a priori collection: GO Slim + Hallmark
fgsea_hg <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO Slim"),
         contrast %in% c("Aging", "Training_Old"))

# Pivot to wide format
fgsea_wide <- fgsea_hg %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Aging), !is.na(NES_Training_Old)) %>%
  mutate(set_size = coalesce(size_Aging, size_Training_Old))

# Classify significance
fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_A = !is.na(padj_Aging) & padj_Aging < 0.05,
    sig_T = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    significance = case_when(
      sig_A & sig_T ~ "Sig Both",
      sig_A         ~ "Sig Aging only",
      sig_T         ~ "Sig Training only",
      TRUE          ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F3)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21)
  )

fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant",
                nrow(fgsea_wide),
                sum(fgsea_wide$database == "Hallmark"),
                sum(fgsea_wide$database == "GO Slim"),
                nrow(fgsea_sig)))

# Spearman correlation
nes_cor_all <- cor.test(fgsea_wide$NES_Aging, fgsea_wide$NES_Training_Old, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Aging, fgsea_sig$NES_Training_Old, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Aging, fgsea_wide$NES_Training_Old))) * 1.15

# Quadrant counts (sig only)
n_rev_tl  <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old > 0)
n_rev_br  <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old < 0)
n_exac_tr <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old > 0)
n_exac_bl <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old < 0)
n_rev_pw  <- n_rev_tl + n_rev_br
n_total_sig <- nrow(fgsea_sig)
pw_rev_frac <- if (n_total_sig > 0) n_rev_pw / n_total_sig else 0

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]))

# Sizes
txt_pw   <- scale_text(BASE_PATHWAY, PW) * 0.90  # pathway labels (reduced for composite)
txt_quad <- scale_text(BASE_QUADRANT, PW) * 0.85  # quadrant labels (harmonized with F06 module NES scatter)

# Conservative display-label shortenings (display-only; pathway IDs unchanged).
# Same map shared with F04 panel B for cross-figure consistency.
PANEL_B_DISPLAY_OVERRIDES <- c(
  "Unfolded Protein Response"         = "UPR",
  "Ribosome Biogenesis"               = "Ribo Bio",
  "Amino Acid Metabolism"             = "AA Metabolism",
  "Plasma Membrane Protein Loc."      = "PM Protein Loc.",
  "Cytoplasmic Translation"           = "Cytoplasmic Transl.",
  "Mitochondrial Organization"        = "Mito Organization",
  "Precursor Metabolites & Energy"    = "Precursor Metab. & Energy",
  "Mitochondrial Transport"           = "Mito Transport",
  "Mitochondrial Protein Import"      = "Mito Protein Import",
  "Mitochondrial Protein Degradation" = "Mito Protein Deg."
)

# Label data — apply display-label overrides for long pathway names.
# Preserves meaning; only shortens for figure readability.
label_pw <- fgsea_sig %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F3[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F3[as.character(significance)],
    pathway_label  = pathway_label %>%
      str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") %>%
      str_replace("Muscle System.*", "Muscle System") %>%
      str_replace("Ketone Metabolic.*", "Ketone Metabolism") %>%
      str_replace("^Trna Metabolic.*", "tRNA Metabolism") %>%
      str_replace("^Establishment Or Maintenance Of Cell Polarity$", "Cell Polarity") %>%
      str_replace("^Generation Of Precursor Metabolites And Energy$",
                  "Precursor Metabolites & Energy") %>%
      str_replace("^Protein Localization To Plasma Membrane$",
                  "Plasma Membrane Protein Loc.") %>%
      str_replace("^Extracellular Matrix Organization$", "ECM Organization") %>%
      str_replace("^Epithelial Mesenchymal Transition$", "EMT") %>%
      str_replace("^Microtubule-Based Movement$", "Microtubule Movement") %>%
      str_replace("^Mitochondrion Organization$", "Mitochondrial Organization") %>%
      str_replace("^UV Response Dn$", "UV Response (Down)") %>%
      dplyr::recode(!!!PANEL_B_DISPLAY_OVERRIDES)
  )

# Split for layered plotting
ns_df  <- fgsea_wide %>% filter(significance == "NS")
sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
  mutate(draw_order = factor(significance,
    levels = c("Sig Training only", "Sig Aging only", "Sig Both"))) %>%
  arrange(draw_order)

# Subtitle
rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", \u03c1(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "GO Slim + Hallmark | %d pathways (%d sig.) | \u03c1 = %.2f [%.2f, %.2f], %s%s | %.0f%% reversed | negative \u03c1 = training opposes aging effects",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", nes_cor_all$p.value)),
  rho_sig_str, pw_rev_frac * 100
)

pB <- ggplot(mapping = aes(x = NES_Aging, y = NES_Training_Old)) +
  # Blue = reversed (off-diagonal)
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  # Red = exacerbated (diagonal)
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points
  geom_point(data = ns_df, aes(shape = database),
             size = 2.5, fill = "grey70", color = "grey55", alpha = 0.40, stroke = 0.4) +
  # Sig points — larger
  geom_point(data = sig_df, aes(fill = significance, size = set_size, shape = database),
             color = ifelse(sig_df$database == "Hallmark", "black", "grey65"),
             alpha = 0.80, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
  scale_shape_manual(values = c("Hallmark" = 24, "GO Slim" = 21), name = "Database") +
  scale_size_continuous(range = c(4, 12), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  # Pathway labels — larger
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_pw, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.5, segment.color = "grey30",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.9, point.padding = 0.6,
                   force = 14, force_pull = 0.2,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant labels — mirrored from Panel A scatter style
  annotate("label", x = nes_lim, y = -nes_lim,
           label = sprintf("Reversed  n = %d", n_rev_br),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#4393C3", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = nes_lim,
           label = sprintf("Reversed  n = %d", n_rev_tl),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#4393C3", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = nes_lim, y = nes_lim,
           label = sprintf("Exacerbated  n = %d", n_exac_tr),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#D6604D", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = -nes_lim,
           label = sprintf("Exacerbated  n = %d", n_exac_bl),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#D6604D", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(title = "Pathway-Level Reversal (fGSEA)",
       subtitle = subtitle_str,
       x = "NES (Aging)",
       y = "NES (Training Old)") +
  FIG_THEME +
  theme(  # Harmonized with F06 module NES scatter — tighter axis/legend typography
    axis.text         = element_text(size = 9,  face = "bold", color = "grey30"),
    axis.title        = element_text(size = 10, face = "bold"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 13, face = "bold", color = "grey25"),
    legend.text       = element_text(size = 11, color = "grey20"),
    legend.key.size   = unit(6, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.spacing.x  = unit(6, "mm"),
    legend.box.margin = margin(t = -1),
    plot.margin       = margin(5, 5, 5, 5)
  ) +
  guides(fill  = "none",
         shape = guide_legend(nrow = 1, order = 1,
                               keyheight = unit(6, "mm"),
                               keywidth  = unit(6, "mm"),
                               override.aes = list(size = 5, fill = "grey50")),
         size  = guide_legend(nrow = 1, order = 2,
                               keyheight = unit(6, "mm"),
                               keywidth  = unit(6, "mm")))

ggsave(file.path(RPT_PNG, "MAIN_panel_B_nes_scatter.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_nes_scatter.pdf"), pB,
       width = PW, height = PH, units = "mm", device = pdf_device)

# Export all terms
fgsea_wide %>%
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
  write_csv(file.path(DAT, "panel_B", "nes_scatter.csv"))

cat("F05 Panel B done\n")
