# F01 Supplementary Panel C: Type I Fiber CSA
# Thin wrapper — all logic lives in shared_prepost_panel.R

cfg <- list(
  dv_col        = "Type_I_fCSA",
  y_label       = expression(bold("Type I fCSA (" * mu * m^2 * ")")),
  delta_label   = expression(bold(Delta ~ "Type I fCSA (" * mu * m^2 * ")")),
  title         = "Type I Fiber CSA",
  tag           = "c",
  output_prefix = "pSC",
  file_tag      = "panel_C_type_I_fcsa",
  audit_file    = "panel_C_type_I_fcsa.csv",
  file_prefix   = "SUPP",
  rpt_png       = "04_Figures/F01/b_reports/supp/png/panels",
  rpt_pdf       = "04_Figures/F01/b_reports/supp/pdf/panels",
  dat           = "04_Figures/F01/c_data/supp",
  coerce_cols     = TRUE,
  filter_complete = TRUE
)

source("04_Figures/F01/a_script/shared_prepost_panel.R")
