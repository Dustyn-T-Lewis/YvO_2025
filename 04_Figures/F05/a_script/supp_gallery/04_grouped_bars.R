# Supplementary Enrichment Gallery — Diverging Bars by Biological Theme (F05)
# All significant pathways included (no top-N selection).
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

RPT <- "04_Figures/F05/b_reports/supp"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

rev <- readRDS(file.path(DAT, "prep_reversal.rds"))

# All sig pathways (no top-N selection)
rev_sel <- rev %>%
  filter(sig_Aging | sig_TO)

rev_bars <- rev_sel %>%
  select(pathway, pathway_label, bio_theme, NES_Aging, NES_TO, padj_Aging,
         padj_TO, sig_Aging, sig_TO) %>%
  pivot_longer(
    cols      = c(NES_Aging, NES_TO),
    names_to  = "contrast",
    values_to = "NES"
  ) %>%
  mutate(
    contrast = recode(contrast,
                      NES_Aging = "Aging",
                      NES_TO    = "Training (Old)"),
    sig = ifelse(contrast == "Aging", sig_Aging, sig_TO),
    alpha_val = ifelse(sig, 1, 0.4),
    pathway_label = clean_pathway_name(pathway, max_chars = 40)
  )

n_rows <- n_distinct(rev_bars$pathway)
fig_h <- max(150, 6 * n_rows + 40)

p_rev <- ggplot(rev_bars,
                 aes(x = NES,
                     y = reorder_within(pathway_label, NES, bio_theme),
                     fill = contrast, alpha = alpha_val)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(bio_theme ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_fill_manual(values = c("Aging"           = unname(CONTRAST_COLORS["Aging"]),
                                "Training (Old)"  = unname(CONTRAST_COLORS["Training_Old"])),
                    name = "Contrast") +
  scale_alpha_identity() +
  labs(
    title    = "Aging vs Training Response by Biological Theme",
    subtitle = sprintf("All %d significant pathways; opposing bars indicate reversal", n_rows),
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y    = element_text(size = 7),
    strip.text.y   = element_text(angle = 0, size = 7, hjust = 0),
    legend.position = "bottom",
    panel.spacing  = unit(2, "mm"),
    panel.border   = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "g_grouped_bars.pdf"), p_rev,
       width = 200, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "g_grouped_bars.png"), p_rev,
       width = 200, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Grouped bar plot saved.\n")
