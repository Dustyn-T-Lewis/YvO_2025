# F06 Supplementary — Per-Module Panels (individual, not composite)
# Each module: triptych (top) + hub network (bottom) saved as a standalone PDF+PNG
# 10 modules ordered by protein count (largest first)
#
# Reads per-module PNGs from panel_B_triptych.R and panel_D_hub.R
# Generates: supp/panels/modules/{pdf,png}/SUPP_F06_module_{color}_{slug}.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

RPT_SRC <- "04_Figures/F06/b_reports/supp/png/modules"
RPT_PDF <- "04_Figures/F06/b_reports/supp/pdf/modules"
RPT_PNG <- "04_Figures/F06/b_reports/supp/png/modules"
DAT     <- "04_Figures/F06/c_data"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- All modules ordered by size (largest first) ---
mod_assign <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
ALL_MODULES <- mod_assign %>%
  filter(module_color != "grey") %>%
  count(module_color, sort = TRUE) %>%
  pull(module_color)

mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
pathway_slug <- setNames(
  gsub("[/ ]+", "_", tolower(mod_bio$bio_label)),
  mod_bio$module_color
)
display_labels <- setNames(mod_bio$display_label, mod_bio$module_color)

mod_counts <- mod_assign %>%
  filter(module_color != "grey") %>%
  count(module_color) %>%
  tibble::deframe()

# --- Read panel PNGs ---
read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) {
    message("  WARNING: Missing ", path)
    return(NULL)
  }
  rasterGrob(readPNG(path), interpolate = TRUE)
}

# --- Build and save one module panel ---
MODULE_W <- 280
MODULE_H <- 260

message("Building individual module panels for ", length(ALL_MODULES), " modules...")

for (mod in ALL_MODULES) {
  slug <- pathway_slug[mod]
  if (is.na(slug)) slug <- mod

  trip_fname <- sprintf("SUPP_triptych_%s_%s.png", mod, slug)
  hub_fname  <- sprintf("SUPP_hub_%s_%s.png", mod, slug)

  trip_grob <- read_panel(trip_fname)
  hub_grob  <- read_panel(hub_fname)

  if (is.null(trip_grob) || is.null(hub_grob)) {
    message("  SKIP ", mod, ": missing panel PNG(s)")
    next
  }

  # Module header
  n_prot <- mod_counts[mod]
  label_text <- sprintf("%s  (%d proteins)", display_labels[mod], n_prot)

  header <- ggplot() +
    annotate("rect", xmin = 0, xmax = 0.06, ymin = 0.15, ymax = 0.85,
             fill = mod, color = "black", linewidth = 0.4) +
    annotate("text", x = 0.09, y = 0.5, label = label_text,
             hjust = 0, size = 4.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 4, 0, 4))

  # Stack: header (5%) + triptych (50%) + hub (45%)
  panel <- wrap_elements(full = header) /
           wrap_elements(full = trip_grob) /
           wrap_elements(full = hub_grob) +
    plot_layout(heights = c(0.05, 0.50, 0.45))

  out_name <- sprintf("SUPP_F06_module_%s_%s", mod, slug)

  ggsave(file.path(RPT_PDF, paste0(out_name, ".pdf")), panel,
         width = MODULE_W, height = MODULE_H, units = "mm",
         device = pdf_device, limitsize = FALSE)
  ggsave(file.path(RPT_PNG, paste0(out_name, ".png")), panel,
         width = MODULE_W, height = MODULE_H, units = "mm",
         dpi = 300, limitsize = FALSE)

  message(sprintf("  Saved %s", out_name))
}

message("F06 supp module panels complete (",
        length(ALL_MODULES), " individual module panels)")
