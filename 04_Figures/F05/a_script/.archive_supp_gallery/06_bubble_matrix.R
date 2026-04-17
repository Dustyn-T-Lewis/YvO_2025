# Supplementary Enrichment Gallery — Bubble Matrix (F05: Reversal)
# All classified pathways: circle size = -log10(padj), fill = NES.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)

RPT <- "04_Figures/F05/b_reports/supp"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

NES_GRADIENT <- scale_fill_gradient2(
  low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
  name = "NES", limits = c(-3.5, 3.5), oob = scales::squish
)

rev <- readRDS(file.path(DAT, "prep_reversal.rds"))

# All classified pathways (no top-N selection)
rev_sel <- rev %>%
  filter(pattern != "Other") %>%
  mutate(pw_label = make.unique(clean_pathway_name(pathway, max_chars = 45), sep = " "))

rev_bub <- rev_sel %>%
  select(pathway, pw_label, pattern, database,
         NES_Aging, NES_TO, padj_Aging, padj_TO) %>%
  pivot_longer(
    cols      = matches("^(NES|padj)_"),
    names_to  = c(".value", "contrast"),
    names_pattern = "(NES|padj)_(.*)"
  ) %>%
  mutate(
    contrast = recode(contrast, Aging = "Aging", TO = "Training\n(Old)"),
    contrast = factor(contrast, levels = c("Aging", "Training\n(Old)")),
    neg_log_p = pmin(-log10(padj), 20),
    neg_log_p = ifelse(is.na(padj) | padj >= 1, 0, neg_log_p)
  )

pw_order <- rev_sel %>%
  arrange(pattern, desc(abs(NES_Aging))) %>%
  pull(pw_label)
rev_bub$pw_label <- factor(rev_bub$pw_label, levels = rev(pw_order))

n_pw <- n_distinct(rev_bub$pw_label)
fig_h <- max(120, 5 * n_pw + 40)

p_rev <- ggplot(rev_bub, aes(x = contrast, y = pw_label)) +
  geom_point(aes(size = neg_log_p, fill = NES),
             shape = 21, stroke = 0.4, color = "grey30") +
  NES_GRADIENT +
  scale_size_continuous(range = c(0.5, 6),
                        name = expression(-log[10](p[adj])),
                        breaks = c(2, 5, 10, 15)) +
  labs(
    title    = "Enrichment: Aging Reversal",
    subtitle = sprintf("All %d classified pathways", n_pw),
    x = NULL, y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y     = element_text(size = 6.5),
    axis.text.x     = element_text(size = 9),
    legend.position = "bottom",
    legend.box      = "horizontal",
    panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "e_bubble_matrix.pdf"), p_rev,
       width = 170, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "e_bubble_matrix.png"), p_rev,
       width = 170, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Bubble matrix plot saved.\n")
