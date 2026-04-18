# F01 — Supplementary Figure Composite (Panels A + B + C)
# Stacks three phenotype panels vertically:
#   A = Deadlift 1RM
#   B = Type II Fiber CSA
#   C = Type I Fiber CSA

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages(library(patchwork))

source("04_Figures/F01/a_script/supp/panels/panel_A.R")
source("04_Figures/F01/a_script/supp/panels/panel_B.R")
source("04_Figures/F01/a_script/supp/panels/panel_C.R")

RPT_PNG <- "04_Figures/F01/b_reports/supp/png"
RPT_PDF <- "04_Figures/F01/b_reports/supp/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# Shared tweaks: tag nudged right+down, tighter title→subtitle gap
supp_tweak <- theme(
  plot.tag          = element_text(face = "bold", size = 15,
                                   margin = margin(t = 3, r = 0, b = 0, l = 4)),
  plot.title        = element_text(margin = margin(b = 1)),
  plot.subtitle     = element_text(margin = margin(t = 0))
)

# Trim margins to reduce inter-row gap (tighter than before: negative bottom/top)
pSA_left  <- pSA_left  + supp_tweak + theme(plot.margin = margin(3, 5, -2, 5))
pSA_right <- pSA_right + supp_tweak + theme(plot.margin = margin(3, 5, -2, 5))
pSB_left  <- pSB_left  + supp_tweak + theme(plot.margin = margin(-2, 5, -2, 5))
pSB_right <- pSB_right + supp_tweak + theme(plot.margin = margin(-2, 5, -2, 5))
pSC_left  <- pSC_left  + supp_tweak + theme(plot.margin = margin(-2, 5, 3, 5))
pSC_right <- pSC_right + supp_tweak + theme(plot.margin = margin(-2, 5, 3, 5))

# Compose each row (same widths as standalone panels)
row_A <- (pSA_left | pSA_right) + plot_layout(widths = c(0.65, 0.35))
row_B <- (pSB_left | pSB_right) + plot_layout(widths = c(0.65, 0.35))
row_C <- (pSC_left | pSC_right) + plot_layout(widths = c(0.65, 0.35))

# Stack vertically — wrap_elements needed for nested patchwork objects
composite <- wrap_elements(row_A) / wrap_elements(row_B) / wrap_elements(row_C) +
  plot_layout(heights = c(1, 1, 1))

COMP_W <- 170   # matches individual panel width
COMP_H <- 185   # tighter vertical spacing (was 195)

ggsave(file.path(RPT_PDF, "SUPP_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

cat("F01 supp composite done\n")
