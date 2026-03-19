# Supplementary Enrichment Gallery — NES Quadrant Scatter (F05: Reversal)
# NES vs NES scatter with quadrant backgrounds and pattern coloring.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
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

# Aging vs Training Old
rev <- readRDS(file.path(DAT, "prep_reversal.rds"))

rev <- rev %>%
  mutate(
    label = case_when(
      pattern == "Reversed" ~ pathway_label,
      pattern == "Aging-specific" & rank(-abs(NES_Aging)) <= 4 ~ pathway_label,
      pattern == "Training-specific" & rank(-abs(NES_TO)) <= 3 ~ pathway_label,
      TRUE ~ NA_character_
    )
  )

r_val  <- cor(rev$NES_Aging, rev$NES_TO, use = "complete.obs")
n_rev  <- sum(rev$pattern == "Reversed", na.rm = TRUE)

p_rev <- ggplot(rev, aes(NES_Aging, NES_TO)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#F3E5F5", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#F3E5F5", alpha = 0.4) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FFF3E0", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFF3E0", alpha = 0.4) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_point(aes(color = pattern, size = set_size), alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 30,
                  segment.size = 0.2, segment.color = "grey50",
                  box.padding = 0.6, min.segment.length = 0.2,
                  force = 2, force_pull = 0.5, seed = 42) +
  scale_color_manual(values = PATTERN_COLORS, name = "Pattern") +
  scale_size_continuous(range = c(1, 5), name = "Set size", guide = "none") +
  labs(
    title    = "Pathway-Level Reversal: Aging vs Training (Old)",
    subtitle = sprintf("r = %.2f, %d reversed pathways (dashed = full reversal)",
                        r_val, n_rev),
    x = "NES \u2014 Aging",
    y = "NES \u2014 Training (Old)"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)) +
  coord_fixed()

ggsave(file.path(RPT, "b_nes_scatter.pdf"), p_rev,
       width = 200, height = 200, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "b_nes_scatter.png"), p_rev,
       width = 200, height = 200, units = "mm", dpi = 300)

cat("NES scatter plot saved.\n")
