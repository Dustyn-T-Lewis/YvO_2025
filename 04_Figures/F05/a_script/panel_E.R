# F3 Panel E: fGSEA NES Scatter (GO Slim + Hallmark, a priori)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

PW <- 200; PH <- 200
RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)
  f2_cache <- "04_Figures/F04/c_data/shared/fgsea_tstat_all_v2.csv"
  f1_cache <- "04_Figures/shared/fgsea_f1_panel_F.csv"
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
PE_W <- 200

# Filter to a priori collection: GO Slim + Hallmark
fgsea_hg <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO Slim"),
         contrast %in% c("Aging", "Training_Old"))

# Pivot ALL terms to wide format
fgsea_wide <- fgsea_hg %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Aging), !is.na(NES_Training_Old)) %>%
  mutate(set_size = coalesce(size_Aging, size_Training_Old))

# Classify significance for ALL terms
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
    pathway_label = clean_pathway_name(pathway)
  )

# Subset significant terms (for labeling and quadrant counts)
fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant",
                nrow(fgsea_wide),
                sum(fgsea_wide$database == "Hallmark"),
                sum(fgsea_wide$database == "GO Slim"),
                nrow(fgsea_sig)))

# Spearman correlation on ALL terms (primary) and sig-only (secondary)
# NES distributions violate normality (Shapiro-Wilk p < 0.0001); Spearman is appropriate
nes_cor_all <- cor.test(fgsea_wide$NES_Aging, fgsea_wide$NES_Training_Old, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Aging, fgsea_sig$NES_Training_Old, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Aging, fgsea_wide$NES_Training_Old))) * 1.15

# Quadrant counts on sig terms only
n_rev_tl  <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old > 0)
n_rev_br  <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old < 0)
n_exac_tr <- sum(fgsea_sig$NES_Aging > 0 & fgsea_sig$NES_Training_Old > 0)
n_exac_bl <- sum(fgsea_sig$NES_Aging < 0 & fgsea_sig$NES_Training_Old < 0)

n_rev_pw <- n_rev_tl + n_rev_br
n_total_sig <- nrow(fgsea_sig)
pw_rev_frac <- if (n_total_sig > 0) n_rev_pw / n_total_sig else 0
pw_rev_binom <- if (n_total_sig > 0) binom.test(n_rev_pw, n_total_sig) else NULL
pw_rev_ci <- if (!is.null(pw_rev_binom)) pw_rev_binom$conf.int * 100 else c(NA, NA)

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f], p = %.2g",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
                nes_cor_all$p.value))
if (!is.null(nes_cor_sig)) {
  nes_ci_sig <- fisher_z_ci(nes_cor_sig$estimate, nrow(fgsea_sig))
  message(sprintf("  NES Spearman (sig): rho = %.3f [%.3f, %.3f], p = %.2g",
                  nes_cor_sig$estimate, nes_ci_sig[1], nes_ci_sig[2],
                  nes_cor_sig$p.value))
}
message(sprintf("  Pathway reversal: %d/%d sig (%.1f%%) [%.1f, %.1f]",
                n_rev_pw, n_total_sig, pw_rev_frac * 100,
                pw_rev_ci[1], pw_rev_ci[2]))

txt_gene <- scale_text(BASE_GENE, PE_W)
txt_quad <- scale_text(BASE_QUADRANT, PE_W)

# Labels on all sig terms (collection is already curated: 69 GO Slim + Hallmark)
label_pw <- fgsea_sig %>%
  mutate(
    label_fill = SIG_LABEL_FILL_F3[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F3[as.character(significance)]
  ) %>%
  mutate(pathway_label = pathway_label %>%
    str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") %>%
    str_replace("Muscle System.*", "Muscle System") %>%
    str_replace("Ketone Metabolic.*", "Ketone Metabolism")
  ) %>%
  mutate(nudge_y = case_when(
    significance == "Sig Both"          ~  0.15,
    significance == "Sig Aging only"    ~ -0.15,
    significance == "Sig Training only" ~  0.10,
    TRUE ~ 0
  )) %>%
  arrange(significance)

# Split data for layered plotting: NS behind, sig on top
ns_df  <- fgsea_wide %>% filter(significance == "NS")
sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "Sig Both"          ~ 0.75,
      significance == "Sig Aging only"    ~ 0.85,
      significance == "Sig Training only" ~ 0.85,
      TRUE ~ 0.60
    ),
    draw_order = factor(significance,
      levels = c("Sig Training only", "Sig Aging only", "Sig Both", "Reversal"))
  ) %>%
  arrange(draw_order)

# Subtitle with Spearman rho (all terms primary, sig secondary)
rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", rho(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "GO Slim + Hallmark (a priori) | %d pathways (%d sig.)\nrho(all) = %.2f [%.2f, %.2f], p %s%s | %.0f%% reversed",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "< 0.001", sprintf("= %.3f", nes_cor_all$p.value)),
  rho_sig_str, pw_rev_frac * 100
)

pE <- ggplot(mapping = aes(x = NES_Aging, y = NES_Training_Old)) +
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
  geom_point(data = ns_df, size = 1.5, fill = "grey70",
             shape = 21, color = "grey55", alpha = 0.40, stroke = 0.4) +
  geom_point(data = sig_df, aes(fill = significance, size = set_size),
             shape = 21, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = txt_gene, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.35, force = 5, force_pull = 0.3,
                   label.padding = unit(1.2, "pt"),
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
    subtitle = subtitle_str,
    x = "NES (Aging)",
    y = "NES (Training Old)"
  ) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_E_nes_scatter.pdf"), pE,
       width = PE_W, height = PE_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_nes_scatter.png"), pE,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)

# Export ALL terms (not just significant)
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
  write_csv(file.path(DAT, "panel_E", "nes_scatter.csv"))

cat("F3 Panel E done\n")
