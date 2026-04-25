# F01 Supplementary Panel B: Type II Fiber CSA
# Thin wrapper — all logic lives in shared_prepost_panel.R

cfg <- list(
  dv_col        = "Type_II_fCSA",
  y_label       = expression(bold("Type II fCSA (" * mu * m^2 * ")")),
  delta_label   = expression(bold(Delta ~ "Type II fCSA (" * mu * m^2 * ")")),
  title         = "Type II Fiber CSA",
  tag           = "b",
  output_prefix = "pSB",
  file_tag      = "panel_B_type_II_fcsa",
  audit_file    = "panel_B_type_II_fcsa.csv",
  file_prefix   = "SUPP",
  rpt_png       = "04_Figures/F01/b_reports/supp/png/panels",
  rpt_pdf       = "04_Figures/F01/b_reports/supp/pdf/panels",
  dat           = "04_Figures/F01/c_data/supp",
  coerce_cols     = TRUE,
  filter_complete = TRUE
)

source("04_Figures/F01/a_script/shared_prepost_panel.R")
