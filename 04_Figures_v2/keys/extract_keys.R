################################################################################
#   Legend Extraction Utilities
#   Shared helpers for extracting legends from panel plot objects
#   and assembling them into per-figure key composites.
#
#   Usage: source this file, then call extract_key() on any ggplot object.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
})

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

# Output directory
KEYS_RPT <- "04_Figures_v2/keys/b_reports"
dir.create(KEYS_RPT, recursive = TRUE, showWarnings = FALSE)

#' Extract legend from a ggplot object
#'
#' Temporarily restores the legend (which was stripped with legend.position="none"
#' in panel scripts), extracts it as a standalone grob, and returns it.
#'
#' @param plot_obj  A ggplot object
#' @param position  Legend position to use for extraction (default "bottom")
#' @param ...       Additional theme arguments passed to the legend
#' @return A legend grob (or NULL if no legend is present)
extract_key <- function(plot_obj, position = "bottom", ...) {
  p_with_legend <- plot_obj +
    theme(legend.position = position,
          legend.text     = element_text(size = 8),
          legend.title    = element_text(size = 9, face = "bold"),
          legend.key.size = unit(4, "mm"),
          ...)

  tryCatch(
    cowplot::get_legend(p_with_legend),
    error = function(e) {
      message("  Warning: could not extract legend: ", conditionMessage(e))
      NULL
    }
  )
}

#' Wrap a legend grob as a ggplot-compatible object for patchwork
#'
#' @param legend_grob  A grob from extract_key() or get_legend()
#' @return A wrap_elements() object usable in patchwork layouts
wrap_key <- function(legend_grob) {
  if (is.null(legend_grob)) return(plot_spacer())
  wrap_elements(legend_grob)
}

#' Save a key composite
#'
#' @param plot_obj  The assembled key composite (patchwork or ggplot)
#' @param filename  Output filename (without path), e.g. "F1_keys"
#' @param width     Width in mm
#' @param height    Height in mm
save_key <- function(plot_obj, filename, width = 280, height = 100) {
  pdf_device <- get_pdf_device()

  ggsave(file.path(KEYS_RPT, paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, units = "mm",
         device = pdf_device, limitsize = FALSE)
  ggsave(file.path(KEYS_RPT, paste0(filename, ".png")), plot_obj,
         width = width, height = height, units = "mm",
         dpi = 300, limitsize = FALSE)

  message(sprintf("  Key saved: %s.pdf/png (%d x %d mm)", filename, width, height))
}

cat("Key extraction utilities loaded\n")
