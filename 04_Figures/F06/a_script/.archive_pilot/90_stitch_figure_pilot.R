# Figure 6 — PILOT v5: Precise alignment of all elements
#
# Panel A (heatmap, no preservation) — left 65%
# Panel B (stacked scatters, no title/legend) — right 35%, aligned to heatmap grid
# Panel B title/subtitle — drawn by stitcher, aligned with Panel A title
# Panel B legend — positioned at same y as Panel A color legends
# Panel letters — slightly above-left of their respective titles

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(patchwork)
library(cowplot)
library(png)
library(grid)

PILOT_DIR <- "04_Figures/F06/b_reports/pilot"

read_panel <- function(file, dir = PILOT_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA     <- read_panel("PILOT_panel_A_heatmap.png")
pB     <- read_panel("PILOT_panel_B_scatters.png")
pB_leg <- read_panel("PILOT_panel_B_legend.png")

COMP_W <- 470
COMP_H <- 300
TAG_SZ <- 16
TITLE_SZ <- 14
SUBTITLE_SZ <- 9.5

# --- Panel A structure (proportions within the PNG) ---
# Title area:    y = 0.92 - 1.00
# Brackets:      y = 0.82 - 0.92
# Heatmap grid:  y = 0.22 - 0.82
# X-axis labels: y = 0.14 - 0.22
# Color legends:  y = 0.06 - 0.14

A_TITLE_Y  <- 0.835   # title baseline — just above the panel tops
A_LEGEND_Y <- 0.14    # raised — bring protein key closer to plots
GRID_TOP   <- 0.80    # top of heatmap grid
GRID_BOT   <- 0.22    # bottom of heatmap grid

# Panel B scatters: align with heatmap grid, shifted left
B_X        <- 0.42    # shifted left — y-axis labels tuck near Panel A edge
B_W        <- 0.56    # widened to fill
B_TITLE_X  <- 0.626   # title/letter sit above the actual scatter plot area
B_GRID_H   <- GRID_TOP - GRID_BOT  # 0.58

# --- Independently tunable label positions (mm from composite origin) ---
# Convert mm → normalized: x_norm = mm / COMP_W,  y_norm = mm / COMP_H
mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

# Panel A labels (mm)
A_LET_X  <- 15;   A_LET_Y  <- 254   # letter "A"
A_TTL_X  <- 44;   A_TTL_Y  <- 250   # title
A_SUB_X  <- 44;   A_SUB_Y  <- 243   # subtitle

# Panel B labels (mm)
B_LET_X  <- 283;  B_LET_Y  <- 254   # letter "B"
B_TTL_X  <- 293;  B_TTL_Y  <- 250   # title
B_SUB_X  <- 293;  B_SUB_Y  <- 243   # subtitle

# Crop canvas to content bounds (removes white space margins)
CROP_L <-  0.01;  CROP_R <- 0.80   # x range
CROP_B <-  0.17;  CROP_T <- 0.85   # y range
SAVE_W <- COMP_W * (CROP_R - CROP_L)  # scaled output width
SAVE_H <- COMP_H * (CROP_T - CROP_B)  # scaled output height

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B drawn FIRST (behind) — its white left margin hides behind Panel A
  draw_grob(pB, x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg, x = 0.70, y = 0.195,
            width = 0.24, height = 0.04, hjust = 0.5, vjust = 0.5) +
  # Panel A drawn SECOND (on top) — covers Panel B's white left margin
  draw_grob(pA, x = 0, y = 0, width = 0.60, height = 0.96,
            hjust = 0, vjust = 0) +
  # Panel A: letter, title, subtitle
  draw_label("A", x = mm2x(A_LET_X), y = mm2y(A_LET_Y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("WGCNA Module\u2013Trait Associations",
             x = mm2x(A_TTL_X), y = mm2y(A_TTL_Y),
             size = TITLE_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label("10 modules | LMM (BH) | Stratified r (BH per-trait)",
             x = mm2x(A_SUB_X), y = mm2y(A_SUB_Y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  # Panel B: letter, title, subtitle
  draw_label("B", x = mm2x(B_LET_X), y = mm2y(B_LET_Y),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Module\u2013Level NES Scatters",
             x = mm2x(B_TTL_X), y = mm2y(B_TTL_Y),
             size = TITLE_SZ, fontface = "bold",
             hjust = 0, vjust = 1) +
  draw_label("fGSEA on module-member t-stat ranks",
             x = mm2x(B_SUB_X), y = mm2y(B_SUB_Y),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

pdf_device <- get_pdf_device()

ggsave(file.path(PILOT_DIR, "PILOT_F06_composite.pdf"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(PILOT_DIR, "PILOT_F06_composite.png"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("PILOT v5 F06 composite saved")
