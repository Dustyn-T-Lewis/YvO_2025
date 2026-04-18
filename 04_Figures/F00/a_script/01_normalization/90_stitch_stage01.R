# F00 Stage 01 — Normalization detail composite
# 2x2 grid, 200x160 mm. Panels: A filter cascade, B missingness, C eta-squared, D outlier consensus.
# Output: 04_Figures/F00/b_reports/01_normalization/{pdf,png}/F00_stage01_normalization.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

# Load shared data once for all panels
int01 <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

# Source panels — each creates a named plot object
source("04_Figures/F00/a_script/01_normalization/panels/panel_A.R")  # -> pA_s01
source("04_Figures/F00/a_script/01_normalization/panels/panel_B.R")  # -> pB_s01
source("04_Figures/F00/a_script/01_normalization/panels/panel_C.R")  # -> pC_s01
source("04_Figures/F00/a_script/01_normalization/panels/panel_D.R")  # -> pD_s01

RPT_PNG <- "04_Figures/F00/b_reports/01_normalization/png"
RPT_PDF <- "04_Figures/F00/b_reports/01_normalization/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

composite_s01 <- (pA_s01 | pB_s01) / (pC_s01 | pD_s01) +
  plot_annotation(
    title = "YvO Proteomics Pipeline \u2014 Stage 01: Normalization detail",
    subtitle = sprintf(
      "Filter cascade | per-protein missingness | biological signal | outlier consensus | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(face = "italic", size = 9,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

COMP_W <- 200; COMP_H <- 160

ggsave(file.path(RPT_PDF, "F00_stage01_normalization.pdf"),
       composite_s01, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_PNG, "F00_stage01_normalization.png"),
       composite_s01, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

# --- Stage excel ---
source("04_Figures/shared/figure_supplement_helpers.R")

s01_specs <- list(
  list(name = "panel_A", path = "04_Figures/F00/c_data/01_normalization/panel_A_filter_cascade.csv"),
  list(name = "panel_B", path = "04_Figures/F00/c_data/01_normalization/panel_B_missingness.csv"),
  list(name = "panel_C", path = "04_Figures/F00/c_data/01_normalization/panel_C_eta_squared.csv"),
  list(name = "panel_D", path = "04_Figures/F00/c_data/01_normalization/panel_D_outlier_consensus.csv")
)
build_workbook(
  "04_Figures/F00/c_data/01_normalization/F00_stage01_normalization.xlsx",
  title = "F00 Stage 01 \u2014 Normalization QC source data",
  description = "Panels A\u2013D source data for normalization detail composite.",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C", "panel_D"),
    Description = c(
      "Panel A: Protein filter cascade (retained vs removed per step)",
      "Panel B: Per-protein missingness distribution (post-filter)",
      "Panel C: Eta-squared distribution (biological signal)",
      "Panel D: Outlier consensus diagnostic (Mahalanobis + QC flags)"),
    stringsAsFactors = FALSE),
  sheet_specs = s01_specs
)
cleanup_after_workbook(s01_specs)

message("Stage 01 composite + excel done")
