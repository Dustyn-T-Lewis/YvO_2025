# Figure 6 — Supplementary Panels
#
# Generates:
#   - QC panels: soft threshold, dendrogram, compartment enrichment, bicor sensitivity
#   - QC composite: SUPP_F06_composite.{pdf,png}
#   - Module panels: per-module triptychs, hub networks + composite, hub chord
#   - Preservation CSV (consumed by Panel A heatmap)
#
# Does NOT source YvO_WGCNA_run.R (run separately).

library(here)
source(here::here("04_Figures", "shared", "style.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))

BASE <- here::here("04_Figures", "F06")

# =============================================================================
# 1. QC panels (individual)
# =============================================================================

message("=== F06 SUPP: Running QC panel scripts ===")
source(here::here("04_Figures", "F06", "a_script", "_supp_qc_soft_threshold.R"))
source(here::here("04_Figures", "F06", "a_script", "_supp_qc_dendrogram.R"))
source(here::here("04_Figures", "F06", "a_script", "_supp_qc_compartment.R"))
source(here::here("04_Figures", "F06", "a_script", "_supp_qc_bicor.R"))

# =============================================================================
# 2. QC composite (from 90_stitch_QC.R logic)
# =============================================================================

message("=== F06 SUPP: Building QC composite ===")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

RPT_SRC <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
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

# Stacked layout: A full width, B full width, C+D side by side
bottom_row <- wrap_elements(full = pC) + wrap_elements(full = pD) +
  plot_layout(widths = c(1, 1))

composite <- wrap_elements(full = pA) /
             wrap_elements(full = pB) /
             bottom_row +
  plot_layout(heights = c(0.30, 0.35, 0.35)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(4, 6, 4, 6)
    )
  )

# Tag placement via cowplot
TAG_SZ <- 16

composite <- ggdraw(composite) +
  draw_label("A", x = 0.02, y = 0.960, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.02, y = 0.660, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.02, y = 0.320, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.52, y = 0.320, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

COMP_W <- 250
COMP_H <- 330

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "SUPP_F06_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_F06_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 supp QC saved: SUPP_F06_composite.{pdf,png}")

# =============================================================================
# 3. Module panels: triptychs + hub networks + chord
# =============================================================================

message("=== F06 SUPP: Running module triptychs (all modules) ===")
source(here::here("04_Figures", "F06", "a_script", "_supp_mod_triptych.R"))

message("=== F06 SUPP: Running module hub networks (all modules) ===")
source(here::here("04_Figures", "F06", "a_script", "_supp_mod_hub.R"))

message("=== F06 SUPP: Running hub chord diagram ===")
source(here::here("04_Figures", "F06", "a_script", "_supp_mod_chord.R"))

# =============================================================================
# 4. Preservation (long-running — 200 permutations)
# =============================================================================

message("=== F06 SUPP: Running module preservation ===")
source(here::here("04_Figures", "F06", "a_script", "_supp_preservation.R"))

message("=== F06 supp panels complete ===")
