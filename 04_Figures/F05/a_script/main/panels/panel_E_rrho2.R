# F05 Panel E: RRHO2 Reversal Map + Per-Quadrant ORA
# Config wrapper for shared/comparison_panels/panel_E_rrho2.R

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id     = "F05",
  t_col_1    = "t_Aging",
  t_col_2    = "t_Training_Old",
  rrho_labels = c("Aging", "Training (Old)"),

  title       = "Threshold-Free Reversal (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm off-diagonal = training reverses aging | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Aging rank"~(Up %->% Down)),
  axis_label_2 = expression("Training (Old) rank"~(Up %->% Down)),

  quadrant_labels = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (Age\u2191 Tr\u2193)",
    DU = "Reversed (Age\u2193 Tr\u2191)"
  ),

  hotspot_export_names = list(
    UU = "Exacerbated_Up",
    DD = "Exacerbated_Down",
    UD = "Reversed_AgingUp_TrainingDown",
    DU = "Reversed_AgingDown_TrainingUp"
  ),

  ora_min_size = 10,
  ora_quadrant_names = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (Aging Up / Training Down)",
    DU = "Reversed (Aging Down / Training Up)"
  ),
  ora_grouped = list(
    file_1_quads     = c("ora_UD", "ora_DU"),    # reversal -> concordant file
    file_2_quads     = c("ora_UU", "ora_DD"),    # exacerbation -> discordant file
    note_if_empty_2  = "No pathways enriched in exacerbation quadrants (padj<0.05)"
  ),
  ora_colors = ORA_QUAD_COLORS_F3,

  summary_quadrant_names = list(
    UU = "Exacerbated_Up",                  UU_slug = "exacerbated_up",
    DD = "Exacerbated_Down",                DD_slug = "exacerbated_down",
    UD = "Reversed_AgingUp_TrainingDown",   UD_slug = "reversed_aging_up",
    DU = "Reversed_AgingDown_TrainingUp",   DU_slug = "reversed_aging_down"
  ),

  rpt_png = "04_Figures/F05/b_reports/main/png/panels",
  rpt_pdf = "04_Figures/F05/b_reports/main/pdf/panels",
  dat     = "04_Figures/F05/c_data",

  supp = list(
    rpt_png        = "04_Figures/F05/b_reports/supp/png/panels",
    rpt_pdf        = "04_Figures/F05/b_reports/supp/pdf/panels",
    ora_bar_title  = "Enriched Pathways by Reversal Quadrant",
    ora_quad_order = c("Reversed (Aging Up / Training Down)",
                       "Reversed (Aging Down / Training Up)",
                       "Exacerbated Up", "Exacerbated Down"),
    ora_quad_short = c(
      "Reversed (Aging Up / Training Down)"  = "Reversed\n(Up \u2192 Down)",
      "Reversed (Aging Down / Training Up)"  = "Reversed\n(Down \u2192 Up)",
      "Exacerbated Up"                       = "Exacerbated Up",
      "Exacerbated Down"                     = "Exacerbated Down"
    )
  )
)

source("04_Figures/shared/comparison_panels/panel_E_rrho2.R")
