# F01 Supplementary Panel A: Deadlift 1RM
# Thin wrapper — all logic lives in shared_prepost_panel.R

cfg <- list(
  dv_col        = "deadlift_1rm_kg",
  y_label       = "Deadlift 1RM (kg)",
  delta_label   = expression(bold(Delta ~ "1RM (kg)")),
  title         = "Deadlift 1RM",
  tag           = "a",
  output_prefix = "pSA",
  file_tag      = "panel_A_deadlift_1rm",
  audit_file    = "panel_A_deadlift_1rm.csv",
  file_prefix   = "SUPP",
  rpt_png       = "04_Figures/F01/b_reports/supp/png/panels",
  rpt_pdf       = "04_Figures/F01/b_reports/supp/pdf/panels",
  dat           = "04_Figures/F01/c_data/supp",
  coerce_cols     = TRUE,
  filter_complete = TRUE
)

source("04_Figures/F01/a_script/shared_prepost_panel.R")
