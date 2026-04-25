# F04 Panel D: fry Rotation Test -- Training Concordance
# Config wrapper for shared/comparison_panels/panel_D_fry.R
#
# Tests whether Training_Young-significant protein sets respond concordantly
# in Old subjects.  No circularity: Training_Young (Young_Post - Young_Pre)
# and Training_Old (Old_Post - Old_Pre) share no contrast terms.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

cfg <- list(
  fig_id           = "F04",
  contrast_source  = "Training_Young",
  contrast_test    = "Training_Old",
  set_prefix       = "ty",

  # Concordance: same direction expected

  expected_up      = "Up",
  expected_down    = "Down",
  driving_up_sign  = "pos",   # up-set driving proteins have t_TO > 0
  driving_dn_sign  = "neg",   # down-set driving proteins have t_TO < 0

  has_circularity  = FALSE,

  # Barcode titles
  up_title_fmt = "Tr.(Y)-Up DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t",
  dn_title_fmt = "Tr.(Y)-Down DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t%s",
  fig_color    = unname(DIR_COLORS["Up"]),       # #D6604D
  stat_corner_up = "topright",
  stat_corner_dn = "bottomleft",

  # ORA labels
  ora_flank_up_label = "Concordant (Up)",
  ora_flank_dn_label = "Concordant (Down)",
  ora_supp_up_label  = "Up DEPs",
  ora_supp_dn_label  = "Down DEPs",
  ora_supp_title     = "Leading-Edge ORA: fry Driving Proteins",
  ora_supp_subtitle  = "Hypergeometric ORA on concordant driving proteins | top 3 per set",

  # Label abbreviation map
  label_map = c(
    "Dynein Recruitment To The Kinetochore"                  = "Dynein-Kinetochore Rec.",
    "Microtubule-Based Movement"                             = "MT-Based Movement",
    "Regulation Of Mitotic Cell Cycle Phase Transition"       = "Reg. Mitotic Cell Cycle Trans.",
    "Negative Regulation Of Nuclear Division"                 = "Neg. Reg. Nuclear Division",
    "Integrin Cell Surface Interactions"                      = "Integrin Cell Surface Int.",
    "Scavenging By Class A Receptors"                        = "Class A Rec. Scav.",
    "Degradation Of The Extracellular Matrix"                = "ECM Degradation",
    "Collagen Chain Trimerization"                           = "Collagen Trimerization",
    "Respiratory Electron Transport"                          = "Respiratory ETC",
    "ATP Synthesis Coupled Electron Transport"                = "ATP Synth. (ETC)",
    "Non Integrin Membrane Ecm Interactions"                  = "Non-integrin ECM",
    "Formation Of The Dystrophin Glycoprotein Complex Dgc"   = "Dystrophin Complex",
    "Cargo Recognition For Clathrin Mediated Endocytosis"    = "Clathrin Cargo Rec.",
    "Class A Receptor Scavenging"                            = "Class A Receptor Scav.",
    "Collagen Trimerization"                                  = "Collagen Trimerization",
    "Hcmv Early Events"                                       = "HCMV Early Events"
  ),
  force_inside_labels = NULL,
  long_label_mode     = "wrap",

  # Titles
  title        = "fry Gene-Set Rotation Test: Training Concordance",
  subtitle_fmt = "Rotation-based set test (exact GSEA analogue, accounts for inter-gene correlation) | dupCor = %.3f | n = %d proteins",

  # Output directories
  panel_w     = 178,
  rpt_png     = "04_Figures/F04/b_reports/main/png/panels",
  rpt_pdf     = "04_Figures/F04/b_reports/main/pdf/panels",
  rpt_sup_png = "04_Figures/F04/b_reports/supp/png/panels",
  rpt_sup_pdf = "04_Figures/F04/b_reports/supp/pdf/panels",
  dat         = "04_Figures/F04/c_data"
)

source("04_Figures/shared/comparison_panels/panel_D_fry.R")
