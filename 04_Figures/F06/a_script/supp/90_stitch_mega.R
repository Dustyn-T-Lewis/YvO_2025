#!/usr/bin/env Rscript
# F06 SUPP — Mega composite for manuscript submission
# Tiles 6 essential SUPP panels into two single-page PDFs + a PNG preview.
# Output:
#   b_reports/supp/SUPP_F06_composite_qc.pdf       (Page 1: A-D, 2x2 grid)
#   b_reports/supp/SUPP_F06_composite_analysis.pdf  (Page 2: E-F, side-by-side)
#   b_reports/supp/SUPP_F06_composite.png           (PNG of page 1 for preview)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

SUPP     <- "04_Figures/F06/b_reports/supp"
SUPP_PDF <- file.path(SUPP, "pdf")
SUPP_PNG <- file.path(SUPP, "png")
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

read_panel <- function(path, tag = NULL) {
  if (!file.exists(path)) { message("  SKIP: ", path); return(NULL) }
  img <- readPNG(path)
  p <- ggplot() +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

# === Page 1: QC (2x2 grid, panels A-D) ===
# Native ratios: A 2.0:1, B 2.0:1, C 1.38:1, D 1.29:1 → avg ~1.67:1
# 2×2 at 300mm wide → each cell 150mm × ~90mm → page 300×180mm
p1_panels <- list(
  read_panel(file.path(SUPP, "01_QC/png/panels/SUPP_soft_threshold.png"),        "A"),
  read_panel(file.path(SUPP, "01_QC/png/panels/SUPP_dendrogram.png"),            "B"),
  read_panel(file.path(SUPP, "01_QC/png/panels/SUPP_compartment_enrichment.png"),"C"),
  read_panel(file.path(SUPP, "01_QC/png/panels/SUPP_bicor_sensitivity.png"),     "D")
)
p1_panels <- Filter(Negate(is.null), p1_panels)
page1 <- if (length(p1_panels) > 0) {
  wrap_plots(p1_panels, ncol = 2) &
    theme(plot.tag = element_text(face = "bold", size = 14))
} else NULL

# === Page 2: Analysis (side-by-side, panels E-F) ===
# Native ratios: E 1.40:1, F 1.18:1 → avg ~1.29:1
# 1×2 at 300mm wide → each cell 150mm × ~116mm → page 300×120mm
p2_panels <- list(
  read_panel(file.path(SUPP, "02_analysis/png/panels/SUPP_age_stratified_lmm.png"),  "E"),
  read_panel(file.path(SUPP, "02_analysis/png/panels/SUPP_upset_dep_modules.png"),   "F")
)
p2_panels <- Filter(Negate(is.null), p2_panels)
page2 <- if (length(p2_panels) > 0) {
  wrap_plots(p2_panels, ncol = 2) &
    theme(plot.tag = element_text(face = "bold", size = 14))
} else NULL

# === Save Page 1: QC composite ===
PW1 <- 300; PH1 <- 180
out_qc <- file.path(SUPP_PDF, "SUPP_F06_composite_qc.pdf")
if (!is.null(page1)) {
  ggsave(out_qc, page1, width = PW1, height = PH1, units = "mm",
         device = pdf_device)
  message("  -> ", out_qc)
}

# === Save Page 2: Analysis composite ===
PW2 <- 300; PH2 <- 120
out_analysis <- file.path(SUPP_PDF, "SUPP_F06_composite_analysis.pdf")
if (!is.null(page2)) {
  ggsave(out_analysis, page2, width = PW2, height = PH2, units = "mm",
         device = pdf_device)
  message("  -> ", out_analysis)
}

# === PNG preview (page 1) ===
out_png <- file.path(SUPP_PNG, "SUPP_F06_composite.png")
if (!is.null(page1)) {
  ggsave(out_png, page1, width = PW1, height = PH1, units = "mm", dpi = 300)
  message("  -> ", out_png)
}

n_out <- sum(!is.null(page1), !is.null(page2))
message(sprintf("F06 SUPP mega: %d composites written", n_out))
