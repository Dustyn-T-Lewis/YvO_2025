#!/usr/bin/env Rscript
# F01 — Phenotype Figure: Master Orchestrator
# Supp panels first (CSVs needed for xlsx), then main (builds xlsx), then supp composite

source(here::here("04_Figures_v2", "shared_functions", "shared_figure_supplement_helpers.R"))

DAT <- here::here("04_Figures_v2", "F01", "c_data")

# supp CSVs first so the workbook can include them
source(here::here("04_Figures_v2", "F01", "a_script", "02_supp_panels.R"))
source(here::here("04_Figures_v2", "F01", "a_script", "01_main_panels.R"))

f01_specs <- list(
  list(name = "panel_A", path = file.path(DAT, "panel_A_training_volume.csv")),
  list(name = "panel_B", path = file.path(DAT, "panel_B_dxa_lbm.csv")),
  list(name = "panel_C", path = file.path(DAT, "panel_C_vl_thickness.csv")),
  list(name = "SUPP_panel_A", path = file.path(DAT, "supp", "panel_A_deadlift_1rm.csv")),
  list(name = "SUPP_panel_B", path = file.path(DAT, "supp", "panel_B_type_II_fcsa.csv")),
  list(name = "SUPP_panel_C", path = file.path(DAT, "supp", "panel_C_type_I_fcsa.csv"))
)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  title = "F01 \u2014 Phenotype source data",
  description = "Main panels A\u2013C and supplementary panels A\u2013C.",
  sheet_specs = f01_specs
)
cleanup_after_workbook(f01_specs,
  extra_subdirs = file.path(DAT, "supp")
)

message("F01 complete: composites in 04_Figures_v2/F01/b_reports, table in ", DAT)
