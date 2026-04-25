# F05 Panel D: fry Rotation Test -- Aging Reversal
# Config wrapper for shared/comparison_panels/panel_D_fry.R
#
# Tests whether Aging-significant protein sets collectively reverse with
# training in Old subjects.
#
# Circularity caveat: Aging (Old_Pre - Young_Pre) and Training_Old (Old_Post -
# Old_Pre) share Old_Pre with opposite signs -> structural negative correlation.
# Pi < 0.05 threshold reduces noise-driven bias vs nominal P < 0.05.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id           = "F05",
  contrast_source  = "Aging",
  contrast_test    = "Training_Old",
  set_prefix       = "aging",

  # Reversal: opposite direction expected
  expected_up      = "Down",
  expected_down    = "Up",
  driving_up_sign  = "neg",   # aging-up set reversal = t_TO < 0
  driving_dn_sign  = "pos",   # aging-down set reversal = t_TO > 0

  has_circularity  = TRUE,

  # Barcode titles
  up_title_fmt = "Aging-Up DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t",
  dn_title_fmt = "Aging-Down DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t%s",
  fig_color    = unname(CONTRAST_COLORS["Aging"]),  # #4CAF50
  stat_corner_up = "bottomleft",
  stat_corner_dn = "topright",

  # ORA labels
  ora_flank_up_label = "Reversed (Up\u2192Down)",
  ora_flank_dn_label = "Reversed (Down\u2192Up)",
  ora_supp_up_label  = "Aging-Up (Reversed)",
  ora_supp_dn_label  = "Aging-Down (Reversed)",
  ora_supp_title     = "Leading-Edge ORA: fry Driving Proteins (Reversal)",
  ora_supp_subtitle  = "Hypergeometric ORA on reversal-driving proteins | top 3 per set",

  # Label abbreviation map
  label_map = c(
    "Amino Acid Catabolic Process"                           = "AA Catabol.",
    "Fatty Acid Catabolic Process"                           = "FA Catabol.",
    "Amino Acid Metabolic Process"                           = "AA Metabolic Process",
    "Establishment Or Maintenance Of Cell Polarity"          = "Cell Polarity",
    "Generation Of Precursor Metabolites And Energy"         = "Precursor Metab. & Energy",
    "Protein Localization To Plasma Membrane"                = "Plasma Membrane Protein Loc.",
    "Organic Acid Catabolic Process"                         = "Organic Acid Catabolism",
    "Membraneless Organelle Assembly"                        = "Membraneless Org. Assembly"
  ),
  force_inside_labels = c("AA Catabol.", "FA Catabol."),
  long_label_mode     = "truncate",

  # Titles
  title        = "fry Gene-Set Rotation Test: Aging Reversal",
  subtitle_fmt = "Rotation-based set test (exact GSEA analogue) | Circularity r = %.3f | dupCor = %.3f | n = %d proteins",

  # Output directories
  panel_w     = 178,
  rpt_png     = "04_Figures/F05/b_reports/main/png/panels",
  rpt_pdf     = "04_Figures/F05/b_reports/main/pdf/panels",
  rpt_sup_png = "04_Figures/F05/b_reports/supp/png/panels",
  rpt_sup_pdf = "04_Figures/F05/b_reports/supp/pdf/panels",
  dat         = "04_Figures/F05/c_data"
)

source("04_Figures/shared/comparison_panels/panel_D_fry.R")
