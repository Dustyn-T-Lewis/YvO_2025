# Supplementary Enrichment Gallery — Pathway Fate Alluvial (F05: Reversal)
# Left = Aging direction, Right = TO direction. Flows colored by pattern.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
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

dir_label <- function(nes, sig) {
  case_when(
    !sig           ~ "NS",
    nes > 0        ~ "Up",
    nes < 0        ~ "Down",
    TRUE           ~ "NS"
  )
}

rev <- readRDS(file.path(DAT, "prep_reversal.rds"))

rev_sankey <- rev %>%
  mutate(
    Aging_dir = dir_label(NES_Aging, sig_Aging) %>% factor(levels = c("Up", "Down", "NS")),
    TO_dir    = dir_label(NES_TO, sig_TO)       %>% factor(levels = c("Up", "Down", "NS"))
  ) %>%
  count(Aging_dir, TO_dir, pattern) %>%
  filter(n > 0)

p_rev <- ggplot(rev_sankey,
                 aes(axis1 = Aging_dir, axis2 = TO_dir, y = n)) +
  geom_alluvium(aes(fill = pattern), width = 1/6, alpha = 0.7,
                curve_type = "sigmoid") +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Aging", "Training (Old)"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = PATTERN_COLORS, name = "Pattern") +
  labs(
    title    = "Pathway Fate: Aging \u2192 Training (Old)",
    subtitle = "Reversed = flows crossing middle (Up\u2192Down or Down\u2192Up)",
    y = "Number of pathways"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "c_sankey.pdf"), p_rev,
       width = 180, height = 140, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "c_sankey.png"), p_rev,
       width = 180, height = 140, units = "mm", dpi = 300)

cat("Sankey alluvial plot saved.\n")
