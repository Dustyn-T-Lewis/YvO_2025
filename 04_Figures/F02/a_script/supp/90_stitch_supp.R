# F02 — Supplementary Figure Composite (Panels A + B + C)
# Sources the three F02 supp panel scripts and composes a 2-row layout:
#   Row 1: A — CV scatter triptych (was Panel C, full width)
#   Row 2: B — CV violins (was Panel D)  |  C — Imputed boxplots (was Panel E)
# Tag font matches F01/F02 main composites (bold, size 15).

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages(library(patchwork))

source("04_Figures/F02/a_script/supp/panels/panel_C.R")
source("04_Figures/F02/a_script/supp/panels/panel_D.R")
source("04_Figures/F02/a_script/supp/panels/panel_E.R")

RPT_DIR <- "04_Figures/F02/b_reports/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

# Re-tag the lead sub-plot of pC (pC12) so the composite shows "A" rather
# than the panel script's standalone "C" tag. Re-build pC_a from the
# re-labelled pC12 + pC3 so the new tag is the only one rendered.
pC12_a <- pC12 + labs(tag = "A")
pC_a   <- (pC12_a | pC3) + plot_layout(widths = c(2, 1))
pD_b   <- pD + labs(tag = "B")
pE_c   <- pE + labs(tag = "C")

# 2-row layout: pC_a fills the top row; pD and pE sit side-by-side in the
# bottom row at proportions matching their natural plot widths.
bottom_row <- (pD_b | pE_c) + plot_layout(widths = c(110, 190))

composite <- (pC_a / bottom_row) +
  plot_layout(heights = c(135, 110)) &
  theme(plot.tag = element_text(face = "bold", size = 15))

COMP_W <- 300
COMP_H <- 250

ggsave(file.path(RPT_DIR, "SUPP_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_DIR, "SUPP_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

cat("F02 supp composite done\n")
