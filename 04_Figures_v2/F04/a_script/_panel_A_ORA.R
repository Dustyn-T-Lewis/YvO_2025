#!/usr/bin/env Rscript
# F04 Panel A ORA wrapper — sets cfg, sources the shared engine.
# Training Concordance: scatter (Training_Young x Training_Old) + quadrant ORA.
# Sourced by 01_main_panels.R — expects style.R already loaded.
# Exports (via engine): composite, scatter_df, n_sig, n_enrich, r_spear.

cfg <- list(
  fig_id = "F04",
  base = here::here("04_Figures_v2", "F04"),
  x_logfc_col = "logFC_Training_Young",
  pi_x_col = "pi_score_Training_Young",
  pi_int_col = "pi_score_Interaction",
  classify_fn = function(pi_x, pi_to, pi_int) {
    classify_proteins_f2(pi_x, pi_to, pi_int)
  },
  ref_slope = 1,
  repel_point_padding = 0.2,
  x_axis_label_expr = expression(log[2] * FC ~ "(Training Young)"),
  sig_colors = SIG_COLORS_F2,
  sig_label_fill = SIG_LABEL_FILL_F2,
  sig_label_text = SIG_LABEL_TEXT_F2,
  quad_names = list(
    TR = "Concordant Up",
    BL = "Concordant Down",
    BR = "Discordant (Y Up / O Down)",
    TL = "Discordant (Y Down / O Up)"
  ),
  ora_bind_order = c("TR", "TL", "BL", "BR"),
  quad_disp = list(
    TR = "Concordant Up",
    BL = "Concordant Down",
    TL = "Discordant (Y\u2193 O\u2191)",
    BR = "Discordant (Y\u2191 O\u2193)"
  ),
  key_lvls = c("Sig Both", "Interaction", "Sig Young only", "Sig Old only"),
  key_display = c("Sig Both", "Interaction", "Sig Young", "Sig Old"),
  key_x = c(1.25, 1.80, 2.45, 3.05),
  key_xlim = c(0.2, 4.5),
  bar_label_inside_frac = 0.50,
  force_outside = character(0),
  force_outside_right = character(0),
  bar_axis_title_mode = "axis_text",
  title = "Training Concordance: Quadrant ORA",
  display_labels = c(
    "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin Endocytosis",
    "The Role Of Gtse1 In G2 M Progression After G2 Checkpoint" = "GTSE1 G2/M Progression",
    "Reference Mitochondrial Complex Ucp1 In Thermogenesis" = "Mito Complex UCP1",
    "Reference Electron Transfer In Complex I" = "e- Transfer\nin Complex I",
    "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex",
    "Striated Muscle Contraction" = "Striated Muscle\nContraction",
    "Metabolism Of Vitamins And Cofactors" = "Vitamins/Cofactors Metab.",
    "Asparagine N Linked Glycosylation" = "Asparagine N-Linked\nGlycosylation",
    "Regulation Of Expression Of Slits And Robos" = "Slit/Robo Expression\nRegulation",
    "Ribosome Associated Quality Control" = "Ribosome QC",
    "Cellular Response To Starvation" = "Starvation Response",
    "ATP Synthesis Coupled Electron Transport" = "ATP Synthesis e- Transport",
    "Respiratory Electron Transport" = "Respiratory ETC",
    "Complex I Biogenesis" = "Complex I Biogenesis",
    "Reference Translation Initiation" = "Translation Initiation",
    "Non Integrin Membrane Ecm Interactions" = "Non-Integrin ECM",
    "Polyol Metabolic Process" = "Polyol Metabolism",
    "Rac1 Gtpase Cycle" = "RAC1 GTPase Cycle",
    "Rac3 Gtpase Cycle" = "RAC3 GTPase Cycle"
  )
)

source(here::here("04_Figures_v2", "shared_functions", "F04-F05_comparison_panels", "panel_A_ora.R"))
