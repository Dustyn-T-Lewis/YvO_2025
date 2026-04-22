# F06 Supplementary — QC Composite (Panels A–D)
# Layout (3 rows, 5-column grid):
#   Row 1: A (soft threshold, full width)
#   Row 2: B (dendrogram, 3 cols)       |  C (compartment enrichment, 2 cols)
#   Row 3: D (bicor sensitivity, full width)
#
# Module-trait heatmap dropped 2026-04-15: duplicates MAIN F06 Panel A blocks.
# The standalone SUPP_module_trait_heatmap.{pdf,png} is retained for reference.
#
# Generates: 01_QC/SUPP_QC_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

RPT_SRC <- "04_Figures/F06/b_reports/supp/01_QC/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/supp/01_QC/pdf"
RPT_PNG <- "04_Figures/F06/b_reports/supp/01_QC/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel("SUPP_soft_threshold.png")
pB <- read_panel("SUPP_dendrogram.png")
pC <- read_panel("SUPP_compartment_enrichment.png")
pD <- read_panel("SUPP_bicor_sensitivity.png")

# Area-based layout: 5-column x 5-row grid (rows 2,4 are spacers)
design <- c(
  area(1, 1, 1, 5),   # A: row 1, full width
  area(3, 1, 3, 3),   # B: row 3, cols 1-3 (dendrogram, 60%)
  area(3, 4, 3, 5),   # C: row 3, cols 4-5 (compartment enrichment, 40%)
  area(5, 1, 5, 5)    # D: row 5, full width (bicor sensitivity)
)

SPACER <- 0.012

composite <- wrap_elements(full = pA) +
  wrap_elements(full = pB) +
  wrap_elements(full = pC) +
  wrap_elements(full = pD) +
  plot_layout(
    design  = design,
    heights = c(0.28, SPACER, 0.35, SPACER, 0.35),
    widths  = c(1, 1, 1, 1, 1)
  ) +
  plot_annotation(
    title    = "Supplementary Figure S-QC: WGCNA Quality Control",
    subtitle = "A: scale-free topology | B: sample dendrogram | C: compartment enrichment | D: bicor sensitivity",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 7, hjust = 0),
      plot.subtitle = element_text(size = 4, face = "bold.italic", color = "grey40", hjust = 0),
      plot.margin   = margin(2, 4, 2, 4)
    )
  )

# --- Tag placement via cowplot ---
TAG_SZ <- 8

# X positions: left edge of each column group
X_LEFT  <- 0.020   # cols 1-3 left edge (A, B, D)
X_RIGHT <- 0.615   # cols 4-5 left edge (C, E)

# Y positions: top of each content row
Y_R1 <- 0.965    # row 1 (A)
Y_R2 <- 0.675    # row 2 (B, C)
Y_R3 <- 0.330    # row 3 (D, E)

composite <- ggdraw(composite) +
  draw_label("a", x = X_LEFT,  y = Y_R1, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("b", x = X_LEFT,  y = Y_R2, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("c", x = X_RIGHT, y = Y_R2, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("d", x = X_LEFT,  y = Y_R3, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

COMP_W <- 178
COMP_H <- 155

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "SUPP_QC_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_QC_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 supp QC composite saved: SUPP_QC_composite.{pdf,png}")
