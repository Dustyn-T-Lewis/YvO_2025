# F01 — Supplementary Figure Composite (Panels A + B + C)
# Stacks three phenotype panels vertically:
#   A = Deadlift 1RM
#   B = Type II Fiber CSA
#   C = Type I Fiber CSA
# Titles, subtitles, and tags placed at composite level via cowplot.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

source("04_Figures/F01/a_script/supp/panels/panel_A.R")
source("04_Figures/F01/a_script/supp/panels/panel_B.R")
source("04_Figures/F01/a_script/supp/panels/panel_C.R")

RPT_PNG <- "04_Figures/F01/b_reports/supp/png"
RPT_PDF <- "04_Figures/F01/b_reports/supp/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# Compose each row (clean — no titles/tags/legends)
row_A <- (pSA_left | pSA_right) + plot_layout(widths = c(0.65, 0.35))
row_B <- (pSB_left | pSB_right) + plot_layout(widths = c(0.65, 0.35))
row_C <- (pSC_left | pSC_right) + plot_layout(widths = c(0.65, 0.35))

# Stack vertically — wrap_elements needed for nested patchwork objects
composite <- (wrap_elements(row_A) / wrap_elements(row_B) / wrap_elements(row_C)) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(0, 0, 0, 0))

COMP_W <- 85    # J Physiol single-column
COMP_H <- 120
TAG_SZ <- composite_text_sizes(COMP_H)$tag
TTL_SZ <- composite_text_sizes(COMP_H)$title
SUB_SZ <- composite_text_sizes(COMP_H)$subtitle

# Row top positions (normalised to COMP_H): 3 equal rows
Y_A <- 0.995;  Y_B <- 0.665;  Y_C <- 0.332
X_TAG <- 0.02
X_TTL <- 0.08
X_SUB <- 0.08
TTL_OFFSET <- 0.000
SUB_OFFSET <- 0.028

composite <- ggdraw(composite) +
  # Panel A
  draw_label("A",           x = X_TAG, y = Y_A,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title,     x = X_TTL, y = Y_A - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_subtitle,  x = X_SUB, y = Y_A - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel B
  draw_label("B",           x = X_TAG, y = Y_B,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title,     x = X_TTL, y = Y_B - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_subtitle,  x = X_SUB, y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel C
  draw_label("C",           x = X_TAG, y = Y_C,              size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title,     x = X_TTL, y = Y_C - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_subtitle,  x = X_SUB, y = Y_C - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "SUPP_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

cat("F01 supp composite done\n")
