# F00 Stage 03 — DEP detail composite
# 2x2 grid, 200x160 mm. Panels: A DEP counts, B effect-size, C power, D blunting.
# Output: 04_Figures/F00/b_reports/03_DEP/{pdf,png}/F00_stage03_DEP.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(patchwork)
  library(ggplot2)
})

# Load shared data once for all panels
da_summ <- read_csv("03_DEP/c_data/02_DA_summary.csv", show_col_types = FALSE)
eff     <- read_csv("03_DEP/c_data/07_effect_size_bootstrap.csv", show_col_types = FALSE)
pwr_dat <- read_csv("03_DEP/c_data/08_power_analysis.csv", show_col_types = FALSE)
blunt   <- read_csv("03_DEP/c_data/06_blunting_diagnostics.csv", show_col_types = FALSE)
combo   <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

# Source panels
source("04_Figures/F00/a_script/03_DEP/panels/panel_A.R")  # -> pA_s03
source("04_Figures/F00/a_script/03_DEP/panels/panel_B.R")  # -> pB_s03
source("04_Figures/F00/a_script/03_DEP/panels/panel_C.R")  # -> pC_s03
source("04_Figures/F00/a_script/03_DEP/panels/panel_D.R")  # -> pD_s03

RPT_PNG <- "04_Figures/F00/b_reports/03_DEP/png"
RPT_PDF <- "04_Figures/F00/b_reports/03_DEP/pdf"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Strip individual panels for composite ---
pA_s03 <- strip_for_composite(pA_s03)
pB_s03 <- strip_for_composite(pB_s03)
pC_s03 <- strip_for_composite(pC_s03)
pD_s03 <- strip_for_composite(pD_s03)

COMP_W <- 200; COMP_H <- 160
txt <- composite_text_sizes(COMP_H)

composite_s03 <- (pA_s03 | pB_s03) / (pC_s03 | pD_s03) +
  plot_annotation(
    title = "YvO Proteomics Pipeline \u2014 Stage 03: DEP detail",
    subtitle = sprintf(
      "DEP counts | effect size | power | training blunting | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = txt$title),
                  plot.subtitle = element_text(face = "italic", size = txt$subtitle,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = txt$tag))

ggsave(file.path(RPT_PDF, "F00_stage03_DEP.pdf"),
       composite_s03, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_PNG, "F00_stage03_DEP.png"),
       composite_s03, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

# --- Stage excel ---
source("04_Figures/shared/figure_supplement_helpers.R")

s03_specs <- list(
  list(name = "panel_A", path = "04_Figures/F00/c_data/03_DEP/panel_A_dep_counts.csv"),
  list(name = "panel_B", path = "04_Figures/F00/c_data/03_DEP/panel_B_effect_size.csv"),
  list(name = "panel_C", path = "04_Figures/F00/c_data/03_DEP/panel_C_power_analysis.csv"),
  list(name = "panel_D", path = "04_Figures/F00/c_data/03_DEP/panel_D_blunting.csv")
)
build_workbook(
  "04_Figures/F00/c_data/03_DEP/F00_stage03_DEP.xlsx",
  title = "F00 Stage 03 \u2014 DEP QC source data",
  description = "Panels A\u2013D source data for DEP detail composite.",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C", "panel_D"),
    Description = c(
      "Panel A: DEP counts per contrast x threshold (heatmap data)",
      "Panel B: Effect-size bootstrap median |logFC| with 95% CI",
      "Panel C: Power analysis per contrast (min detectable Cohen d)",
      "Panel D: Training blunting scatter (|logFC| Young vs Old)"),
    stringsAsFactors = FALSE),
  sheet_specs = s03_specs
)
cleanup_after_workbook(s03_specs)

message("Stage 03 composite + excel done")
