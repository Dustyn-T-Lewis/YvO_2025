# F01 Panel C: Vastus Lateralis Thickness
# Thin wrapper — all logic lives in shared_prepost_panel.R

cfg <- list(
  dv_col       = "VL_thick_cm",
  y_label      = "VL thickness (cm)",
  delta_label  = expression(bold(Delta ~ "VL thickness (cm)")),
  title        = "VL Thickness",
  tag          = "c",
  output_prefix = "pC",
  file_tag     = "panel_C_vl_thickness",
  audit_file   = "panel_C_vl_thickness.csv",
  file_prefix  = "MAIN",
  rpt_png      = "04_Figures/F01/b_reports/main/png/panels",
  rpt_pdf      = "04_Figures/F01/b_reports/main/pdf/panels",
  dat          = "04_Figures/F01/c_data",
  use_plotmath_subtitle = TRUE,
  y_breaks     = c(0, 0.5, 1.0),
  y_labels     = c("0", ".5", "1")
)

source("04_Figures/F01/a_script/shared_prepost_panel.R")
