#!/usr/bin/env Rscript
# F07 Supplementary Composite Stitch — 3 output files
#
# Reads pre-rendered PNGs from b_reports/supp/png/panels/ and composites them.
#
# Output 1: SUPP_F07_composite_main.pdf  — Page 1, panels A+B side by side
# Output 2: SUPP_F07_composite_loso.pdf  — Page 2, panels C-E (LOSO + classifier)
# Output 3: SUPP_F07_composite.png       — PNG of page 1
#
# Panel sources (all in b_reports/supp/png/panels/):
#   A = Per-module ROC grid          (SUPP_F07_module_grid.png)
#   B = Panel B full sweep            (SUPP_F07_panel_B_grid.png)
#   C = LOSO fixed-module sensitivity (SUPP_F07_loso_sensitivity.png)
#   D = LOSO WGCNA-refit              (SUPP_F07_loso_wgcna_refit.png)
#   E = Classifier decomposition      (SUPP_F07_multivariate_classifier.png)

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork); library(cowplot); library(png); library(grid)
})

BASE       <- "04_Figures/F07"
RPT_PNG    <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF    <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PANELS <- file.path(RPT_PNG, "panels")

read_panel <- function(file) {
  path <- file.path(RPT_PANELS, file)
  if (!file.exists(path)) stop("Missing: ", path)
  img <- readPNG(path)
  list(grob = rasterGrob(img, interpolate = TRUE),
       aspect = dim(img)[2] / dim(img)[1])   # width/height
}

TAG_SZ <- 8

# Page 1 panels (A + B, aspect-matched)
pA <- read_panel("SUPP_F07_module_grid.png")
pB <- read_panel("SUPP_F07_panel_B_grid.png")

COMP_H <- 166   # mm tall
wA_mm  <- COMP_H * pA$aspect
wB_mm  <- COMP_H * pB$aspect
COMP_W <- wA_mm + wB_mm

page1 <- (wrap_elements(full = pA$grob) |
           wrap_elements(full = pB$grob)) +
  plot_layout(widths = c(wA_mm, wB_mm)) &
  theme(plot.margin = margin(0, 0, 0, 0))

TAG_X_A <- 5 / COMP_W
TAG_X_B <- (wA_mm + 5) / COMP_W
TAG_Y   <- 0.992

page1_final <- ggdraw(page1) +
  draw_label("A", x = TAG_X_A, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = TAG_X_B, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1)

# Page 2 panels (C + D top row, E centered below)
pC <- read_panel("SUPP_F07_loso_sensitivity.png")
pD <- read_panel("SUPP_F07_loso_wgcna_refit.png")
pE <- read_panel("SUPP_F07_multivariate_classifier.png")

P2_W <- 178  # mm
P2_H <- 166  # mm

page2 <- ((wrap_elements(full = pC$grob) | wrap_elements(full = pD$grob)) /
            (plot_spacer() | wrap_elements(full = pE$grob) | plot_spacer())) +
  plot_layout(heights = c(1, 1.14),
              widths  = c(1, 2, 1)) &
  theme(plot.margin = margin(2, 2, 2, 2))

page2_final <- ggdraw(page2) +
  draw_label("C", x = 0.01, y = 0.99, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.51, y = 0.99, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("E", x = 0.26, y = 0.46, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1)

# Write outputs
graphics.off()
pdf_device <- get_pdf_device()

# PDF 1: page 1 (main panels A+B)
ggsave(file.path(RPT_PDF, "SUPP_F07_composite_main.pdf"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)

# PDF 2: page 2 (LOSO panels C-E)
ggsave(file.path(RPT_PDF, "SUPP_F07_composite_loso.pdf"), page2_final,
       width = P2_W, height = P2_H, units = "mm",
       device = pdf_device, limitsize = FALSE)

# PNG = page 1 only
ggsave(file.path(RPT_PNG, "SUPP_F07_composite.png"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message(sprintf("Wrote SUPP_F07_composite_main.pdf (%.0f x %.0f mm)", COMP_W, COMP_H))
message(sprintf("Wrote SUPP_F07_composite_loso.pdf (%.0f x %.0f mm)", P2_W, P2_H))
message(sprintf("Wrote SUPP_F07_composite.png (%.0f x %.0f mm)", COMP_W, COMP_H))
