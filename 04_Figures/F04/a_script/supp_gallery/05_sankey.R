# Supplementary Enrichment Gallery — Pathway Fate Alluvial (F04: Blunting)
# Left = TY direction, Right = TO direction. Flows colored by pattern.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
})

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

dir_label <- function(nes, sig) {
  case_when(
    !sig           ~ "NS",
    nes > 0        ~ "Up",
    nes < 0        ~ "Down",
    TRUE           ~ "NS"
  )
}

blunt <- readRDS(file.path(DAT, "prep_blunting.rds"))

blunt_sankey <- blunt %>%
  mutate(
    TY_dir = dir_label(NES_TY, sig_TY) %>% factor(levels = c("Up", "Down", "NS")),
    TO_dir = dir_label(NES_TO, sig_TO) %>% factor(levels = c("Up", "Down", "NS"))
  ) %>%
  count(TY_dir, TO_dir, pattern) %>%
  filter(n > 0)

p_blunt <- ggplot(blunt_sankey,
                   aes(axis1 = TY_dir, axis2 = TO_dir, y = n)) +
  geom_alluvium(aes(fill = pattern), width = 1/6, alpha = 0.7,
                curve_type = "sigmoid") +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Training (Young)", "Training (Old)"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = PATTERN_COLORS, name = "Pattern") +
  labs(
    title    = "Pathway Fate: Training Young \u2192 Training Old",
    subtitle = "Blunting = flows from Up/Down to NS",
    y = "Number of pathways"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_blank()
  )

ggsave(file.path(RPT, "c_sankey.pdf"), p_blunt,
       width = 180, height = 140, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "c_sankey.png"), p_blunt,
       width = 180, height = 140, units = "mm", dpi = 300)

cat("Sankey alluvial plot saved.\n")
