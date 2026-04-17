#!/usr/bin/env Rscript
# F06 SUPP — Mega composite for manuscript submission
# Tiles all SUPP sub-composites + standalone panels into one multi-page PDF.
# Output: b_reports/supp/SUPP_F06_composite.{pdf,png} (page 1 = first-page PNG)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

SUPP <- "04_Figures/F06/b_reports/supp"
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

make_page <- function(panels, title, subtitle, w = 380, h = 500) {
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) return(NULL)
  wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(face = "italic", size = 9, color = "grey30"))
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))
}

# === Page 1: QC ===
p1_panels <- list(
  A = read_panel(file.path(SUPP, "01_QC/SUPP_soft_threshold.png"), "A"),
  B = read_panel(file.path(SUPP, "01_QC/SUPP_dendrogram.png"), "B"),
  C = read_panel(file.path(SUPP, "01_QC/SUPP_compartment_enrichment.png"), "C"),
  D = read_panel(file.path(SUPP, "01_QC/SUPP_module_trait_heatmap.png"), "D"),
  E = read_panel(file.path(SUPP, "01_QC/SUPP_bicor_sensitivity.png"), "E")
)
page1 <- make_page(p1_panels,
  "S7a — WGCNA Quality Control",
  "A: Soft threshold | B: Dendrogram | C: Compartment enrichment | D: Module-trait heatmap | E: Bicor sensitivity")

# === Page 2: Analysis ===
p2_panels <- list(
  A = read_panel(file.path(SUPP, "02_analysis/SUPP_age_stratified_lmm.png"), "A"),
  B = read_panel(file.path(SUPP, "02_analysis/SUPP_exemplar_strip.png"), "B"),
  C = read_panel(file.path(SUPP, "02_analysis/SUPP_upset_dep_modules.png"), "C"),
  D = read_panel(file.path(SUPP, "02_analysis/SUPP_me_trajectories.png"), "D"),
  E = read_panel(file.path(SUPP, "02_analysis/SUPP_preservation_vs_effect.png"), "E")
)
page2 <- make_page(p2_panels,
  "S7b — Module-Trait Analysis",
  "A: Age-stratified LMM | B: Exemplar strip | C: DEP-module UpSet | D: ME trajectories | E: Preservation vs effect")

# === Page 3: Module triptychs ===
p3_panels <- list(
  A = read_panel(file.path(SUPP, "03_module/SUPP_triptych_turquoise_lipid_catabolism.png"), "A"),
  B = read_panel(file.path(SUPP, "03_module/SUPP_triptych_blue_cell_cycle_proteostasis.png"), "B"),
  C = read_panel(file.path(SUPP, "03_module/SUPP_triptych_green_oxidative_phosphorylation.png"), "C"),
  D = read_panel(file.path(SUPP, "03_module/SUPP_triptych_brown_redox_glycolysis.png"), "D"),
  E = read_panel(file.path(SUPP, "03_module/SUPP_triptych_magenta_translation_initiation.png"), "E")
)
page3 <- make_page(p3_panels,
  "S7c — Module Triptychs (Key Modules)",
  "A: Turquoise (Lipid) | B: Blue (Cell Cycle) | C: Green (OxPhos) | D: Brown (Redox) | E: Magenta (Translation)")

# === Page 4: Hub networks ===
p4_panels <- list(
  A = read_panel(file.path(SUPP, "03_module/SUPP_hub_turquoise_lipid_catabolism.png"), "A"),
  B = read_panel(file.path(SUPP, "03_module/SUPP_hub_blue_cell_cycle_proteostasis.png"), "B"),
  C = read_panel(file.path(SUPP, "03_module/SUPP_hub_green_oxidative_phosphorylation.png"), "C"),
  D = read_panel(file.path(SUPP, "03_module/SUPP_hub_brown_redox_glycolysis.png"), "D"),
  E = read_panel(file.path(SUPP, "03_module/SUPP_hub_magenta_translation_initiation.png"), "E"),
  F = read_panel(file.path(SUPP, "03_module/SUPP_hub_chord.png"), "F")
)
page4 <- make_page(p4_panels,
  "S7d — Hub Protein Networks (Key Modules)",
  "A: Turquoise | B: Blue | C: Green | D: Brown | E: Magenta | F: Hub chord diagram")

# === Assemble multi-page PDF ===
pages <- Filter(Negate(is.null), list(page1, page2, page3, page4))
PW <- 380; PH <- 500

outpdf <- file.path(SUPP, "SUPP_F06_composite.pdf")
outpng <- file.path(SUPP, "SUPP_F06_composite.png")

pdf(outpdf, width = PW / 25.4, height = PH / 25.4)
for (p in pages) print(p)
dev.off()

# PNG = page 1 only (for writing dir preview)
ggsave(outpng, pages[[1]], width = PW, height = PH, units = "mm", dpi = 300)

message(sprintf("F06 SUPP mega-composite: %d pages -> %s", length(pages), outpdf))
