# F03 SUPP — Composite stitch
# Sources 4 supplementary panels (B, C, D, E) and composes 2x2 grid.
# Output: b_reports/supp/SUPP_F03_composite.{pdf,png} (178 x 155 mm)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(cowplot)
})

source("04_Figures/F03/a_script/supp/panels/panel_B_phist.R")       # -> pB (stripped)
source("04_Figures/F03/a_script/supp/panels/panel_C_ma.R")          # -> pC (stripped)
source("04_Figures/F03/a_script/supp/panels/panel_D_sensitivity.R") # -> pD (stripped)
source("04_Figures/F03/a_script/supp/panels/panel_E_outlier.R")     # -> pE (stripped)

RPT_PDF <- "04_Figures/F03/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/F03/b_reports/supp/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_W <- 178
COMP_H <- 155
CTS    <- composite_text_sizes(COMP_H)

# Compose 2x2 grid — panels arrive pre-stripped
grid <- (pB | pC) / (pD | pE)

composite <- ggdraw(grid) +
  # --- Supra-title ---
  draw_label("F03 Supplementary \u2014 DEP diagnostics",
             x = 0.01, y = 0.998, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$title) +
  # --- Panel a (top-left) ---
  draw_label("A", x = 0.02, y = 0.965, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$tag) +
  draw_label(pB_title, x = 0.06, y = 0.965, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$title) +
  draw_label(pB_subtitle, x = 0.06, y = 0.937, hjust = 0, vjust = 1,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30") +
  # --- Panel b (top-right) ---
  draw_label("B", x = 0.52, y = 0.965, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$tag) +
  draw_label(pC_title, x = 0.56, y = 0.965, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$title) +
  draw_label(pC_subtitle, x = 0.56, y = 0.937, hjust = 0, vjust = 1,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30") +
  # --- Panel c (bottom-left) ---
  draw_label("C", x = 0.02, y = 0.50, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$tag) +
  draw_label(pD_title, x = 0.06, y = 0.50, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$title) +
  draw_label(pD_subtitle, x = 0.06, y = 0.472, hjust = 0, vjust = 1,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30") +
  # --- Panel d (bottom-right) ---
  draw_label("D", x = 0.52, y = 0.50, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$tag) +
  draw_label(pE_title, x = 0.56, y = 0.50, hjust = 0, vjust = 1,
             fontface = "bold", size = CTS$title) +
  draw_label(pE_subtitle, x = 0.56, y = 0.472, hjust = 0, vjust = 1,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30")

ggsave(file.path(RPT_PDF, "SUPP_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F03 SUPP composite saved -> %s / %s", RPT_PDF, RPT_PNG))
