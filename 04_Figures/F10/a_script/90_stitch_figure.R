# F7 — Composite Figure Assembly
# Layout: A on top, B | C in middle, D on bottom

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
})

# Source panel scripts (each creates its panel object)
source("04_Figures/F10/a_script/panel_A.R")
source("04_Figures/F10/a_script/panel_B.R")
source("04_Figures/F10/a_script/panel_C.R")
source("04_Figures/F10/a_script/panel_D.R")

RPT <- "04_Figures/F10/b_reports"

composite <- (
  panel_A /
  (panel_B | panel_C) /
  panel_D
) +
  plot_layout(heights = c(0.25, 0.45, 0.30)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

dev <- get_pdf_device()
ggsave(file.path(RPT, "F7_composite.pdf"), composite,
       width = 380, height = 340, units = "mm", device = dev)
ggsave(file.path(RPT, "F7_composite.png"), composite,
       width = 380, height = 340, units = "mm", dpi = 300, bg = "white")

message("F7 composite saved")
