# Supplementary Enrichment Gallery — Reversal Ratio Lollipop (F05)
# All pathways included; alpha encodes denominator significance.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

RPT <- "04_Figures/F05/b_reports/supp"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

PATTERN_COLORS <- c(
  Reversed            = "#7B5EA7",
  Exacerbated         = "#FF8F00",
  "Aging-specific"    = "#E05A4E",
  "Training-specific" = "#5DA5DA",
  Other               = "grey60"
)

# Reversal ratio (-NES_TO / NES_Aging) for all pathways
rev <- readRDS(file.path(DAT, "prep_reversal.rds"))

rev_lollipop <- rev %>%
  filter(abs(NES_Aging) > 0.01) %>%
  mutate(
    ratio = pmin(pmax(-NES_TO / NES_Aging, -1.5), 1.5),
    sig_denom = sig_Aging,
    pathway_label = clean_pathway_name(pathway, max_chars = 50),
    pathway_label = fct_reorder(pathway_label, ratio)
  )

n_pw <- nrow(rev_lollipop)
fig_h <- max(120, 5 * n_pw + 30)

p_rev <- ggplot(rev_lollipop, aes(x = ratio, y = pathway_label)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#FFCDD2", alpha = 0.08) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#FFE0B2", alpha = 0.08) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#E1BEE7", alpha = 0.08) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7B5EA7",
             linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = ratio, yend = pathway_label,
                   color = pattern, alpha = ifelse(sig_denom, 1, 0.4)),
               linewidth = 0.5) +
  geom_point(aes(color = pattern, alpha = ifelse(sig_denom, 1, 0.4)), size = 2) +
  annotate("text", x = 1, y = Inf, label = "Full reversal", hjust = -0.1,
           vjust = 1.5, size = 3.2, color = "#7B5EA7", fontface = "italic") +
  annotate("text", x = 0, y = Inf, label = "No effect", hjust = 1.1,
           vjust = 1.5, size = 3.2, color = "grey50", fontface = "italic") +
  scale_color_manual(values = PATTERN_COLORS, name = "Pattern") +
  scale_alpha_identity() +
  labs(
    title    = "Aging Reversal Ratio by Pathway",
    subtitle = sprintf("-NES(Training Old) / NES(Aging) for all %d pathways (ratio capped at \u00b11.5)", n_pw),
    x = "Reversal ratio (positive = reversed)",
    y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y = element_text(size = 6.5),
    legend.position = "bottom",
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "f_waterfall.pdf"), p_rev,
       width = 200, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "f_waterfall.png"), p_rev,
       width = 200, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Waterfall lollipop plot saved.\n")
