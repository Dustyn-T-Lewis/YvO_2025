#!/usr/bin/env Rscript
# F04 SUPP — Composite stitch
# Reads pre-rendered SUPP PNGs and tiles into a composite.
# Output: b_reports/supp/SUPP_F04_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT_DIR    <- "04_Figures/F04/b_reports/supp"
PANELS_DIR <- file.path(RPT_DIR, "panels")
pdf_device <- get_pdf_device()

# -- Shared legend strip (rendered once per composite) -----------------------

CONTRAST_COL <- c(
  Aging          = "#D6604D",
  Training_Young = "#4393C3",
  Training_Old   = "#F57C00",
  Interaction    = "#7B1FA2"
)
MULTI_COL <- "#388E3C"

render_legend_strip <- function(path, width_mm = 460, height_mm = 55) {
  png(path, width = width_mm, height = height_mm, units = "mm", res = 300,
      bg = "white")
  par(mar = c(0.3, 0.3, 0.3, 0.3))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  n_grad <- 60
  nes_cols <- colorRampPalette(c("#4393C3", "white", "#D6604D"))(n_grad)
  lfc_cols <- colorRampPalette(c("#4393C3", "white", "#D6604D"))(n_grad)

  # Block 1: Primary Contrast swatches
  text(0.01, 0.90, "Primary Contrast (protein)", adj = c(0, 1),
       cex = 0.95, font = 2)
  ctr_labs <- c(names(CONTRAST_COL), "Multi (3+)")
  ctr_cols <- c(CONTRAST_COL, Multi = MULTI_COL)
  x_sw <- seq(0.015, 0.30, length.out = length(ctr_labs))
  sw_w <- 0.03
  for (i in seq_along(ctr_labs)) {
    rect(x_sw[i], 0.48, x_sw[i] + sw_w, 0.68,
         col = ctr_cols[i], border = "grey40", lwd = 0.4)
    text(x_sw[i] + sw_w / 2, 0.40, gsub("_", " ", ctr_labs[i]),
         adj = c(0.5, 1), cex = 0.58)
  }

  # Block 2: Protein logFC gradient
  text(0.38, 0.90, "Protein logFC (inner bar)", adj = c(0, 1),
       cex = 0.95, font = 2)
  bx <- seq(0.38, 0.62, length.out = n_grad + 1)
  for (j in seq_len(n_grad)) {
    rect(bx[j], 0.48, bx[j + 1], 0.68, col = lfc_cols[j], border = NA)
  }
  rect(bx[1], 0.48, bx[n_grad + 1], 0.68, col = NA, border = "grey40", lwd = 0.4)
  text(bx[1], 0.38, "-2", adj = c(0, 1), cex = 0.65)
  text(mean(range(bx)), 0.38, "0", adj = c(0.5, 1), cex = 0.65)
  text(bx[n_grad + 1], 0.38, "+2", adj = c(1, 1), cex = 0.65)

  # Block 3: Pathway NES gradient
  text(0.68, 0.90, "Pathway NES (outer arc, panels C\u2013E)",
       adj = c(0, 1), cex = 0.95, font = 2)
  bx2 <- seq(0.68, 0.97, length.out = n_grad + 1)
  for (j in seq_len(n_grad)) {
    rect(bx2[j], 0.48, bx2[j + 1], 0.68, col = nes_cols[j], border = NA)
  }
  rect(bx2[1], 0.48, bx2[n_grad + 1], 0.68, col = NA, border = "grey40", lwd = 0.4)
  text(bx2[1], 0.38, "-3", adj = c(0, 1), cex = 0.65)
  text(mean(range(bx2)), 0.38, "0", adj = c(0.5, 1), cex = 0.65)
  text(bx2[n_grad + 1], 0.38, "+3", adj = c(1, 1), cex = 0.65)

  # Footnote
  text(0.50, 0.10,
       paste0("Panel B uses its own Pathway Direction (mean logFC) and ",
              "Pathway Identity keys, baked into the chord \u2014 see panel title."),
       adj = c(0.5, 0.5), cex = 0.55, font = 3, col = "grey40")

  dev.off()
}

shared_legend_path <- file.path(PANELS_DIR, "SUPP_shared_legend_chord.png")
render_legend_strip(shared_legend_path)

# -- Build composite ---------------------------------------------------------

read_panel <- function(path, tag = NULL) {
  if (!file.exists(path)) {
    message("  SKIP: ", basename(path))
    return(NULL)
  }
  img <- readPNG(path)
  p <- ggplot() +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                       interpolate = TRUE) +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1))
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

panels <- list(
  A = read_panel(file.path(PANELS_DIR, "SUPP_enrichment_heatmap.png"), "A"),
  B = read_panel(file.path(PANELS_DIR, "SUPP_ora_chord_sq.png"), "B"),
  C = read_panel(file.path(PANELS_DIR, "SUPP_panel_F_fgsea_TY_sq.png"), "C"),
  D = read_panel(file.path(PANELS_DIR, "SUPP_panel_F_fgsea_TO_sq.png"), "D"),
  E = read_panel(file.path(PANELS_DIR, "SUPP_panel_F_fgsea_Inter_sq.png"), "E"),
  L = read_panel(shared_legend_path, NULL)
)
panels <- Filter(Negate(is.null), panels)

layout <- "
AA
BC
DE
LL
"

composite <- wrap_plots(panels, design = layout,
                         widths  = c(1, 1),
                         heights = c(184, 230, 230, 55)) +
  plot_annotation(
    title = "Supplementary Figure S3 \u2014 Training Concordance Diagnostics",
    subtitle = paste0(
      "A: Pathway-level enrichment heatmap (direction + pattern). ",
      "B: ORA\u2013DEP chord across Training DEPs. ",
      "C\u2013E: fGSEA leading-edge chords (top 3 up + 3 down per contrast). ",
      "Shared colour keys below."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(face = "italic", size = 10,
                                    color = "grey30", lineheight = 1.15),
      plot.margin   = margin(6, 6, 6, 6))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 16),
        plot.margin = margin(1, 1, 1, 1))

COMP_W <- 460
COMP_H <- 720

ggsave(file.path(RPT_DIR, "SUPP_F04_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "SUPP_F04_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F04 SUPP composite saved -> %s (%d panels)",
                RPT_DIR, length(panels)))
