# F06 Supplementary — Analysis Composite (Panels A–C within 02_analysis)
# Layout: vertical stack (each panel full width)
#   A (age-stratified LMM heatmap)
#   B (exemplar strip boxplots)
#   C (UpSet DEP-modules)
#
# Generates: 02_analysis/SUPP_analysis_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F06/b_reports/supp/02_analysis"

read_panel <- function(file, dir = RPT) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel("SUPP_age_stratified_lmm.png")
pB <- read_panel("SUPP_exemplar_strip.png")
pC <- read_panel("SUPP_upset_dep_modules.png")

layout <- "A\n#\nB\n#\nC"

composite <- wrap_elements(full = pA) +
  wrap_elements(full = pB) + wrap_elements(full = pC) +
  plot_layout(
    design = layout,
    heights = c(0.42, 0.012, 0.22, 0.012, 0.334),
    widths  = c(1)
  ) +
  plot_annotation(
    title    = "Supplementary Figure S-Analysis: Module\u2013Trait Companions",
    subtitle = "A: age-stratified LMM | B: eigengene exemplar strip | C: DEP\u2013module overlap",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 8, face = "bold.italic", color = "grey40", hjust = 0),
      plot.margin   = margin(2, 4, 2, 4)
    )
  )

# --- Tag placement via cowplot ---
TAG_SZ <- 12

# Vertical stack: each panel's top edge
# Heights: A=0.42, spacer=0.012, B=0.22, spacer=0.012, C=0.334
Y_A <- 0.975
Y_B <- 0.555
Y_C <- 0.325

composite <- ggdraw(composite) +
  draw_label("A", x = 0.02, y = Y_A, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.02, y = Y_B, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.02, y = Y_C, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

COMP_W <- 450
COMP_H <- 580

pdf_device <- get_pdf_device()

ggsave(file.path(RPT, "SUPP_analysis_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "SUPP_analysis_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 supp analysis composite saved: ",
        file.path(RPT, "SUPP_analysis_composite.{pdf,png}"))
