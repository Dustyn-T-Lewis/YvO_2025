# F01 — Single-Column Composite (Panels A + B + C stacked vertically)
# For inline embedding in a two-column LaTeX manuscript (~3.3 in / 84 mm wide).
# Sources the same panel scripts as 90_stitch_figure.R; does NOT modify them.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(magick))

source("04_Figures/F01/a_script/main/panels/panel_A.R")
source("04_Figures/F01/a_script/main/panels/panel_B.R")
source("04_Figures/F01/a_script/main/panels/panel_C.R")

RPT_PNG     <- "04_Figures/F01/b_reports/main/png"
RPT_PDF     <- "04_Figures/F01/b_reports/main/pdf"

# --- Helpers ----
# Shrink geom_signif p-value text for compact format
shrink_signif <- function(p, new_size) {
  for (i in seq_along(p$layers)) {
    if ("textsize" %in% names(p$layers[[i]]$aes_params))
      p$layers[[i]]$aes_params$textsize <- new_size
  }
  p
}

# Crop whitespace from rendered image (actual pixel cropping, not margin compression)
crop_whitespace <- function(path, fuzz = 5) {
  img <- image_read(path)
  img <- image_trim(img, fuzz = fuzz)
  image_write(img, path)
  invisible(path)
}

# --- Build compact panel versions for single-column width ----
# Reduce text sizes for the narrower format.
# plot.tag.location = "panel" places the tag inside the panel area
# rather than in a dedicated margin column, eliminating tag-caused whitespace.
compact_theme <- theme(
  plot.title         = element_text(size = 6, face = "bold", margin = margin(b = 0.5)),
  plot.subtitle      = element_text(size = 4.5, margin = margin(t = 0)),
  plot.tag           = element_text(size = 8, face = "bold",
                                     margin = margin(t = -2, l = 2)),
  plot.tag.position  = "topleft",
  plot.tag.location  = "panel",
  axis.title         = element_text(size = 5),
  axis.text          = element_text(size = 4.5),
  legend.text        = element_text(size = 4.5),
  legend.title       = element_text(size = 5),
  panel.border       = element_rect(colour = "grey30", fill = NA, linewidth = 0.3),
  plot.margin        = margin(2, 2, 0, 2)
)

# Panel A: standalone bar chart — shrink text + p-value labels.
# Use plot.tag.location = "plot" so tag sits above the title (not inside panel).
# B/C use "panel" (tag inside panel border), so align them horizontally by
# shifting A slightly and B/C leftward via negative tag margin.
pA_sc <- shrink_signif(pA, 1.5) + compact_theme +
  theme(plot.tag.location = "plot",
        plot.margin = margin(7, 2, -1, 2))

# Panels B and C: each has left (pre/post) + right (delta) sub-panel.
# Y-axis title: fully specified element_text to avoid partial-merge issues
# with the parent axis.title from compact_theme / FIG_THEME.
bc_y_close <- theme(axis.title.y = element_text(face = "bold", size = 5,
                                                  margin = margin(r = -15)))
pB_left_sc  <- shrink_signif(pB_left, 1.5)  + compact_theme + bc_y_close +
  theme(plot.tag.location = "plot",
        plot.margin = margin(3, 2, -1, 2))
pB_right_sc <- shrink_signif(pB_right, 1.5) + compact_theme +
  theme(plot.margin = margin(3, 2, -1, 2))
pC_left_sc  <- shrink_signif(pC_left, 1.5)  + compact_theme + bc_y_close +
  theme(plot.tag.location = "plot",
        plot.margin = margin(3, 2, 1, 2))
pC_right_sc <- shrink_signif(pC_right, 1.5) + compact_theme +
  theme(plot.margin = margin(3, 2, 1, 2))

pB_sc <- (pB_left_sc | pB_right_sc) + plot_layout(widths = c(0.65, 0.35))
pC_sc <- (pC_left_sc | pC_right_sc) + plot_layout(widths = c(0.65, 0.35))

# Stack A / B / C vertically
composite_sc <- (pA_sc / pB_sc / pC_sc) +
  plot_layout(heights = c(1.1, 0.8, 0.8))

# --- Output dimensions ---
# Render slightly wider than target; magick trims excess whitespace after.
SC_W <- 86    # mm (slightly oversized; trimmed below)
SC_H <- 125   # mm

png_path <- file.path(RPT_PNG, "MAIN_F01_composite_single_col.png")
pdf_path <- file.path(RPT_PDF, "MAIN_F01_composite_single_col.pdf")

ggsave(pdf_path, composite_sc,
       width = SC_W, height = SC_H, units = "mm", device = get_pdf_device())
ggsave(png_path, composite_sc,
       width = SC_W, height = SC_H, units = "mm", dpi = 300)

# Actual pixel cropping — trims uniform-color borders from all edges
crop_whitespace(png_path)

cat("F01 single-column composite done\n")
