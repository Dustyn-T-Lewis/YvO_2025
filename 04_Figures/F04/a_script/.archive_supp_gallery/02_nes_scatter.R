# Supplementary Enrichment Gallery — NES Quadrant Scatter (F04: Blunting)
# NES vs NES scatter with quadrant backgrounds and pattern coloring.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)
library(ggrepel)

RPT <- "04_Figures/F04/b_reports/supp"
DAT <- "04_Figures/F04/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

PATTERN_COLORS <- c(
  Blunted            = "#E05A4E",
  Concordant         = "#2E7D32",
  "Old-specific"     = "#5DA5DA",
  "Interaction-only" = "#7B5EA7",
  Discordant         = "#FF8F00",
  Other              = "grey60"
)

# Training Young vs Training Old
blunt <- readRDS(file.path(DAT, "prep_blunting.rds"))

blunt <- blunt %>%
  mutate(
    label = case_when(
      pattern == "Concordant" & rank(-abs(NES_TY)) <= 4 ~ pathway_label,
      pattern == "Blunted" & rank(-abs(NES_TY)) <= 5 ~ pathway_label,
      pattern == "Old-specific" & rank(-abs(NES_TO)) <= 3 ~ pathway_label,
      pattern == "Interaction-only" & rank(-abs(NES_Int)) <= 2 ~ pathway_label,
      TRUE ~ NA_character_
    )
  )

r_val  <- cor(blunt$NES_TY, blunt$NES_TO, use = "complete.obs")
n_conc <- sum(sign(blunt$NES_TY) == sign(blunt$NES_TO), na.rm = TRUE)
pct_conc <- round(100 * n_conc / nrow(blunt), 1)

p_blunt <- ggplot(blunt, aes(NES_TY, NES_TO)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E8F5E9", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#E3F2FD", alpha = 0.4) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#E3F2FD", alpha = 0.4) +
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
    title    = "Pathway-Level Training Response: Young vs Old",
    subtitle = sprintf("r = %.2f, %d/%d concordant (%.1f%%)",
                        r_val, n_conc, nrow(blunt), pct_conc),
    x = "NES \u2014 Training (Young)",
    y = "NES \u2014 Training (Old)"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)) +
  coord_fixed()

ggsave(file.path(RPT, "b_nes_scatter.pdf"), p_blunt,
       width = 200, height = 200, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "b_nes_scatter.png"), p_blunt,
       width = 200, height = 200, units = "mm", dpi = 300)

cat("NES scatter plot saved.\n")
