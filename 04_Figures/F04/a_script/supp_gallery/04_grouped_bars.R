# Supplementary Enrichment Gallery — Diverging Bars by Biological Theme (F04)
# All significant pathways included (no top-N selection).
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

RPT <- "04_Figures/F04/b_reports/supp"
DAT <- "04_Figures/F04/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

blunt <- readRDS(file.path(DAT, "prep_blunting.rds"))

# All sig pathways (no top-N selection)
blunt_sel <- blunt %>%
  filter(sig_TY | sig_TO)

# Reshape to long for paired bars
blunt_bars <- blunt_sel %>%
  select(pathway, pathway_label, bio_theme, NES_TY, NES_TO, padj_TY, padj_TO,
         sig_TY, sig_TO) %>%
  pivot_longer(
    cols      = c(NES_TY, NES_TO),
    names_to  = "contrast",
    values_to = "NES"
  ) %>%
  mutate(
    contrast = recode(contrast,
                      NES_TY = "Training (Young)",
                      NES_TO = "Training (Old)"),
    sig = ifelse(contrast == "Training (Young)", sig_TY, sig_TO),
    alpha_val = ifelse(sig, 1, 0.4),
    pathway_label = clean_pathway_name(pathway, max_chars = 40)
  )

n_rows <- n_distinct(blunt_bars$pathway)
fig_h <- max(150, 6 * n_rows + 40)

p_blunt <- ggplot(blunt_bars,
                   aes(x = NES,
                       y = reorder_within(pathway_label, NES, bio_theme),
                       fill = contrast, alpha = alpha_val)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(bio_theme ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_fill_manual(values = c("Training (Young)" = unname(CONTRAST_COLORS["Training_Young"]),
                                "Training (Old)"  = unname(CONTRAST_COLORS["Training_Old"])),
                    name = "Contrast") +
  scale_alpha_identity() +
  labs(
    title    = "Training Response by Biological Theme",
    subtitle = sprintf("All %d significant pathways; opacity reflects padj < 0.05", n_rows),
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

ggsave(file.path(RPT, "g_grouped_bars.pdf"), p_blunt,
       width = 200, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "g_grouped_bars.png"), p_blunt,
       width = 200, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Grouped bar plot saved.\n")
