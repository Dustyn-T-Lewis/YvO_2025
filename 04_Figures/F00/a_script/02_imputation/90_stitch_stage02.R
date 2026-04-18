# F00 Stage 02 — Imputation detail composite
# 2x2 grid, 200x160 mm. Panels: A miss-class scatter, B benchmark, C MNAR audit, D sample intensity.
# Output: 04_Figures/F00/b_reports/02_imputation/{pdf,png}/F00_stage02_imputation.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

# Load shared data once for all panels
int02 <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

# Source panels
source("04_Figures/F00/a_script/02_imputation/panels/panel_A.R")  # -> pA_s02
source("04_Figures/F00/a_script/02_imputation/panels/panel_B.R")  # -> pB_s02
source("04_Figures/F00/a_script/02_imputation/panels/panel_C.R")  # -> pC_s02
source("04_Figures/F00/a_script/02_imputation/panels/panel_D.R")  # -> pD_s02

RPT_PNG <- "04_Figures/F00/b_reports/02_imputation/png"
RPT_PDF <- "04_Figures/F00/b_reports/02_imputation/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

composite_s02 <- (pA_s02 | pB_s02) / (pC_s02 | pD_s02) +
  plot_annotation(
    title = "YvO Proteomics Pipeline \u2014 Stage 02: Imputation detail",
    subtitle = sprintf(
      "Missingness x class | method benchmark | MNAR audit | sample integrity | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(face = "italic", size = 9,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

COMP_W <- 200; COMP_H <- 160

ggsave(file.path(RPT_PDF, "F00_stage02_imputation.pdf"),
       composite_s02, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_PNG, "F00_stage02_imputation.png"),
       composite_s02, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

# --- Stage excel ---
source("04_Figures/shared/figure_supplement_helpers.R")

s02_specs <- list(
  list(name = "panel_A", path = "04_Figures/F00/c_data/02_imputation/panel_A_miss_class_scatter.csv"),
  list(name = "panel_B", path = "04_Figures/F00/c_data/02_imputation/panel_B_benchmark_top5.csv"),
  list(name = "panel_C", path = "04_Figures/F00/c_data/02_imputation/panel_C_mnar_audit.csv"),
  list(name = "panel_D", path = "04_Figures/F00/c_data/02_imputation/panel_D_sample_intensity.csv")
)
build_workbook(
  "04_Figures/F00/c_data/02_imputation/F00_stage02_imputation.xlsx",
  title = "F00 Stage 02 \u2014 Imputation QC source data",
  description = "Panels A\u2013D source data for imputation detail composite.",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C", "panel_D"),
    Description = c(
      "Panel A: Per-protein missingness x MAR/MNAR classification scatter",
      "Panel B: 23-method benchmark composite scores (top 5)",
      "Panel C: MNAR audit imputation shift distribution",
      "Panel D: Sample-level mean intensity pre vs post imputation"),
    stringsAsFactors = FALSE),
  sheet_specs = s02_specs
)
cleanup_after_workbook(s02_specs)

message("Stage 02 composite + excel done")
