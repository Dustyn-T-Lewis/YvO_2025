# F02 — Supplementary Figure Composite (Panels A + B + C)
# Sources the three F02 supp panel scripts and composes a 2-row layout:
#   A — CV scatter triptych (was Panel C, full width)
#   B — CV violins (was Panel D)  |  C — Imputed boxplots (was Panel E)
# Titles, subtitles, and tags placed at composite level via cowplot.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

source("04_Figures/F02/a_script/supp/panels/panel_C.R")
source("04_Figures/F02/a_script/supp/panels/panel_D.R")
source("04_Figures/F02/a_script/supp/panels/panel_E.R")

RPT_PDF <- "04_Figures/F02/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/F02/b_reports/supp/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# Row 1: panel A (CV scatter triptych, full width) — stripped pC12 + pC3
row_A <- (pC12 | pC3) + plot_layout(widths = c(2, 1))

# Row 2: panel B (CV violins) | panel C (imputed boxplots)
bottom_row <- (pD | pE) + plot_layout(widths = c(110, 190))

# Stack vertically
composite <- (wrap_elements(row_A) / wrap_elements(bottom_row)) +
  plot_layout(heights = c(135, 110)) &
  theme(plot.margin = margin(0, 0, 0, 0))

COMP_W <- 178
COMP_H <- 150
TAG_SZ <- composite_text_sizes(COMP_H)$tag
TTL_SZ <- composite_text_sizes(COMP_H)$title
SUB_SZ <- composite_text_sizes(COMP_H)$subtitle

# Layout positions (normalised coordinates)
# Row 1 (panel A): top of figure
Y_A <- 0.995
# Row 2 starts at y = 1 - 135/(135+110) ≈ 0.449
Y_B <- 0.449
# Panel C (bottom-right) starts at x = 110/(110+190) ≈ 0.367
X_C_TAG <- 0.367

X_TAG <- 0.02
X_TTL <- 0.08
X_SUB <- 0.08
TTL_OFFSET <- 0.000
SUB_OFFSET <- 0.028

composite <- ggdraw(composite) +
  # Panel A (top row, full width)
  draw_label("A",           x = X_TAG,  y = Y_A,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title,     x = X_TTL,  y = Y_A - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_subtitle,  x = X_SUB,  y = Y_A - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel B (bottom-left)
  draw_label("B",           x = X_TAG,  y = Y_B,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title,     x = X_TTL,  y = Y_B - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_subtitle,  x = X_SUB,  y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel C (bottom-right)
  draw_label("C",           x = X_C_TAG,       y = Y_B,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title,     x = X_C_TAG + 0.06, y = Y_B - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_subtitle,  x = X_C_TAG + 0.06, y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "SUPP_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

cat("F02 supp composite done\n")
