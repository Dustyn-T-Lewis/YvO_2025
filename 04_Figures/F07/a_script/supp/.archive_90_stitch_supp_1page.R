# F07 Supplementary Composite Stitch
#
# Layout (side-by-side, height-matched, tight):
#   A (left)  = Per-module ROC grid  (SUPP_F07_module_grid.png)
#   B (right) = Panel B full sweep   (SUPP_F07_panel_B_grid.png)
#
# Output: 04_Figures/F07/b_reports/supp/SUPP_F07_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork); library(cowplot); library(png); library(grid)
})

# Data-producing panel scripts are sourced by MAIN (90_stitch_figure.R) so
# the MAIN xlsx workbook has all CSVs before cleanup. This stitcher runs
# after MAIN and composites the pre-rendered PNGs at b_reports/supp/panels/.

RPT <- "04_Figures/F07/b_reports/supp"
RPT_PANELS <- file.path(RPT, "panels")

read_panel <- function(file) {
  path <- file.path(RPT_PANELS, file)
  if (!file.exists(path)) stop("Missing: ", path)
  img <- readPNG(path)
  list(grob = rasterGrob(img, interpolate = TRUE),
       aspect = dim(img)[2] / dim(img)[1])   # width/height
}

pA <- read_panel("SUPP_F07_module_grid.png")
pB <- read_panel("SUPP_F07_panel_B_grid.png")

# Both panels rendered at same composite height. Compute widths from
# their native aspect ratios so nothing gets stretched.
COMP_H <- 420   # mm tall (~16.5 in)
wA_mm  <- COMP_H * pA$aspect
wB_mm  <- COMP_H * pB$aspect
COMP_W <- wA_mm + wB_mm

composite <- (wrap_elements(full = pA$grob) |
              wrap_elements(full = pB$grob)) +
  plot_layout(widths = c(wA_mm, wB_mm)) &
  theme(plot.margin = margin(0, 0, 0, 0))

TAG_SZ <- composite_text_sizes(COMP_H)$tag

TAG_X_A <- 5 / COMP_W
TAG_X_B <- (wA_mm + 5) / COMP_W
TAG_Y   <- 0.992

composite_final <- ggdraw(composite) +
  draw_label("A", x = TAG_X_A, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = TAG_X_B, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1)

graphics.off()  # close any stale devices left open by sourced panels
pdf_device <- get_pdf_device()

ggsave(file.path(RPT, "SUPP_F07_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "SUPP_F07_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message(sprintf("Wrote SUPP_F07_composite (%.0f x %.0f mm)", COMP_W, COMP_H))
