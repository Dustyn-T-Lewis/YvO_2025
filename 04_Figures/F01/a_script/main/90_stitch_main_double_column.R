# F01 — Double-Column Composite (A left | B/C stacked right)
# J Physiol double-column (178 mm): Panel A occupies the left column at
# full height; Panels B and C stack vertically on the right.
# Title-group sizes match the single-column MAIN_F01_composite (TAG=8,
# TTL=7, SUB=4).  B/C plot regions extend vertically to align with A.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(ggplot2)
})

source("04_Figures/F01/a_script/main/panels/panel_A.R")
source("04_Figures/F01/a_script/main/panels/panel_B.R")
source("04_Figures/F01/a_script/main/panels/panel_C.R")

RPT_PNG <- "04_Figures/F01/b_reports/main/png"
RPT_PDF <- "04_Figures/F01/b_reports/main/pdf"

# --- Compose B and C sub-panels (left pre/post | right delta) ---
pB_comp <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

# --- Double-column layout: A spans both rows left, B top-right, C bot-right ---
# Panel A bottom lowered ~1mm from prior (13→10pt) to increase height.
# Panel B gets minimal top margin; Panel C gets minimal bottom margin.
composite <- wrap_elements(full = pA + theme(plot.margin = margin(8, 2, 10, 2))) +
  wrap_elements(full = pB_comp & theme(plot.margin = margin(4, 2, 2, 2))) +
  wrap_elements(full = pC_comp & theme(plot.margin = margin(4, 2, 4, 2))) +
  plot_layout(
    design = "AB\nAC",
    widths  = c(0.35, 0.65),
    heights = c(1, 1)
  )

COMP_W <- 178   # J Physiol double-column
COMP_H <- 75    # taller than 65mm so B/C plots extend ±5mm toward A edges

# Text sizes: single-col reference minus 1 (TAG=7, TTL=6, SUB=4).
TAG_SZ <- 7
TTL_SZ <- 6
SUB_SZ <- 4

# --- Tag + title + subtitle placement via cowplot ---
X_A   <- 0.010
X_BC  <- 0.360
X_TTL <- 0.030   # title offset right of tag
SUB_OFFSET <- 0.028   # subtitle up 1.5mm from 0.035 (tighter to title)

# Row vertical positions — title up 1mm, tag up 2mm from prior.
Y_TOP <- 0.984     # title up 1mm from 0.971
Y_MID <- 0.516     # title up 1mm from 0.503
TAG_DY <- -0.004   # tag up 2mm total → sits slightly above title baseline

composite <- ggdraw(composite) +
  # Panel A (left, full height)
  draw_label("A", x = X_A,         y = Y_TOP - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = X_A + X_TTL, y = Y_TOP,     size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET,
           label = pA_subtitle, parse = TRUE, hjust = 0, vjust = 1,
           size = SUB_SZ / .pt, colour = "grey30") +
  # Panel B (top-right)
  draw_label("B", x = X_BC,         y = Y_TOP - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = X_BC + X_TTL, y = Y_TOP,     size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_BC + X_TTL, y = Y_TOP - SUB_OFFSET,
           label = pB_subtitle, parse = TRUE, hjust = 0, vjust = 1,
           size = SUB_SZ / .pt, colour = "grey30") +
  # Panel C (bottom-right)
  draw_label("C", x = X_BC,         y = Y_MID - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = X_BC + X_TTL, y = Y_MID,     size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_BC + X_TTL, y = Y_MID - SUB_OFFSET,
           label = pC_subtitle, parse = TRUE, hjust = 0, vjust = 1,
           size = SUB_SZ / .pt, colour = "grey30")

ggsave(file.path(RPT_PDF, "MAIN_F01_composite_double_col.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "MAIN_F01_composite_double_col.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
cat("F01 double-column composite done\n")
