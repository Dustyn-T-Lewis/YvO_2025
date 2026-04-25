# F04 Panel E: RRHO2 Concordance Map + Per-Quadrant ORA
# Config wrapper for shared/comparison_panels/panel_E_rrho2.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id     = "F04",
  t_col_1    = "t_Training_Young",
  t_col_2    = "t_Training_Old",
  rrho_labels = c("Training (Young)", "Training (Old)"),

  title       = "Threshold-Free Concordance (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm corners = concordant gene regulation | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Training (Young) rank"~(Up %->% Down)),
  axis_label_2 = expression("Training (Old) rank"~(Up %->% Down)),

  quadrant_labels = list(
    UU = "Concordant Up",
    DD = "Concordant Down",
    UD = "Discordant (Y\u2191 O\u2193)",
    DU = "Discordant (Y\u2193 O\u2191)"
  ),

  hotspot_export_names = list(
    UU = "UU", DD = "DD", UD = "UD", DU = "DU"
  ),

  ora_min_size = 15,
  ora_quadrant_names = list(
    UU = "Concordant Up",
    DD = "Concordant Down",
    UD = "Discordant (Y Up / O Down)",
    DU = "Discordant (Y Down / O Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UU", "ora_DD"),
    file_2_quads = c("ora_UD", "ora_DU")
  ),
  ora_colors = ORA_QUAD_COLORS_F2,

  summary_quadrant_names = list(
    UU = "Concordant_Up",       UU_slug = "concordant_up",
    DD = "Concordant_Down",     DD_slug = "concordant_down",
    UD = "Discordant_YoungUp_OldDown",  UD_slug = "discordant_y_up",
    DU = "Discordant_YoungDown_OldUp",  DU_slug = "discordant_y_down"
  ),

  rpt_png = "04_Figures/F04/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F04/b_reports/main/pdf/panels",
  dat     = "04_Figures/F04/c_data",

  supp = NULL   # F04 has no supplementary ORA bar chart
)

source("04_Figures/shared/comparison_panels/panel_E_rrho2.R")
