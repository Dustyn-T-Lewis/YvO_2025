# F01 Panel B: DXA Lean Body Mass
# Thin wrapper — all logic lives in shared_prepost_panel.R

cfg <- list(
  dv_col       = "DXA_LBM_kg",
  y_label      = "DXA LBM (kg)",
  delta_label  = expression(bold(Delta ~ "DXA LBM (kg)")),
  title        = "DXA Lean Body Mass",
  tag          = "b",
  output_prefix = "pB",
  file_tag     = "panel_B_dxa_lbm",
  audit_file   = "panel_B_dxa_lbm.csv",
  file_prefix  = "MAIN",
  rpt_png      = "04_Figures/F01/b_reports/main/png/panels",
  rpt_pdf      = "04_Figures/F01/b_reports/main/pdf/panels",
  dat          = "04_Figures/F01/c_data",
  use_plotmath_subtitle = TRUE
)

source("04_Figures/F01/a_script/shared_prepost_panel.R")
