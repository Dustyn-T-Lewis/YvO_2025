# F03 SUPP — Composite stitch
# Sources 4 supplementary panels (B, C, D, E) and composes 2x2 grid.
# Output: b_reports/supp/SUPP_F03_composite.{pdf,png} (300 x 260 mm)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

source("04_Figures/F03/a_script/supp/panels/panel_B_phist.R")       # -> pB
source("04_Figures/F03/a_script/supp/panels/panel_C_ma.R")          # -> pC
source("04_Figures/F03/a_script/supp/panels/panel_D_sensitivity.R") # -> pD
source("04_Figures/F03/a_script/supp/panels/panel_E_outlier.R")     # -> pE

RPT_PDF <- "04_Figures/F03/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/F03/b_reports/supp/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

composite <- (pB | pC) / (pD | pE) +
  plot_annotation(
    title = "F03 Supplementary — DEP diagnostics",
    subtitle = sprintf(
      "p-value distribution | MA plots | imputation rank stability | outlier sensitivity | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(face = "italic", size = 10,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 15))

COMP_W <- 300
COMP_H <- 260

ggsave(file.path(RPT_PDF, "SUPP_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F03 SUPP composite saved -> %s / %s", RPT_PDF, RPT_PNG))
