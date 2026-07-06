#!/usr/bin/env Rscript
# F05 Panel A ORA wrapper — sets cfg, sources the shared engine.
# Aging Reversal: scatter (Aging x Training_Old) + quadrant ORA.
# Sourced by 01_main_panels.R — expects style.R already loaded.
# Exports (via engine): composite, scatter_df, n_sig, n_enrich, r_spear.

cfg <- list(
  fig_id = "F05",
  base = "04_Figures_v2/F05",
  x_logfc_col = "logFC_Aging",
  pi_x_col = "pi_score_Aging",
  pi_int_col = NULL,
  classify_fn = function(pi_x, pi_to, pi_int) {
    classify_proteins_f3(pi_x, pi_to)
  },
  ref_slope = -1,
  repel_point_padding = 0.3,
  x_axis_label_expr = expression(log[2] * FC ~ "(Aging)"),
  sig_colors = SIG_COLORS_F3,
  sig_label_fill = SIG_LABEL_FILL_F3,
  sig_label_text = SIG_LABEL_TEXT_F3,
  quad_names = list(
    TR = "Exacerbated Up",
    BL = "Exacerbated Down",
    BR = "Reversed (Aging Up / Training Down)",
    TL = "Reversed (Aging Down / Training Up)"
  ),
  ora_bind_order = c("TL", "TR", "BL", "BR"),
  quad_disp = list(
    TR = "Exacerbated Up",
    BL = "Exacerbated Down",
    TL = "Reversed (Age\u2193 Tr\u2191)",
    BR = "Reversed (Age\u2191 Tr\u2193)"
  ),
  key_lvls = c("Sig Both", "Sig Aging only", "Sig Training only"),
  key_display = c("Sig Both", "Sig Aging", "Sig Training"),
  key_x = c(1.25, 1.90, 2.55),
  key_xlim = c(0.2, 4.0),
  bar_label_inside_frac = 0.10,
  force_outside = c("ECM Organization", "ECM Proteoglycans", "NCAM Signaling"),
  force_outside_right = c("Protein Localization to Mito."),
  bar_axis_title_mode = "raw5",
  title = "Aging Reversal: Quadrant ORA",
  display_labels = c(
    "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex",
    "Endosomal Sorting Complex Required For Transport Escrt" = "ESCRT Complex",
    "Adaptive Immune Response Based On Somatic Recombination Of Immune Receptors Built From Immunoglobulin Superfamily Domains" = "Adaptive Immune (Ig)",
    "Modified Amino Acid Metabolic Process" = "Modified AA Metabolism",
    "Cell Surface Interactions At The Vascular Wall" = "Vascular Wall Int.",
    "Fcgamma Receptor Fcgr Dependent Phagocytosis" = "Fc\u03b3R Phagocytosis",
    "Monocarboxylic Acid Catabolic Process" = "Monocarb.\nAcid Catab.",
    "Extracellular Matrix Organization" = "ECM Organization",
    "Ecm Proteoglycans" = "ECM Proteoglycans",
    "Eukaryotic Translation Elongation" = "Translation Elongation",
    "Reference Translation Initiation" = "Translation Initiation",
    "Mitochondrial Protein Import" = "Mito Protein Import",
    "Mitochondrial Protein Degradation" = "Mito Protein\nDegradation",
    "Protein Localization" = "Protein\nLocalization",
    "Protein Localization To Mitochondrion" = "Protein Localization to Mito.",
    "Small Molecule Catabolic Process" = "Small Mol. Catabolism",
    "Lipid Catabolic Process" = "Lipid Catabolism",
    "Lipid Modification" = "Lipid Modificat.",
    "Fatty Acid Beta Oxidation" = "FA Beta-Oxidation",
    "Ncam Signaling For Neurite Out Growth" = "NCAM Signaling",
    "Xenobiotic Metabolism" = "Xenobiotic\nMetab.",
    "Mitochondrial Transport" = "Mito Transport",
    "Non Integrin Membrane Ecm Interactions" = "Non-Integrin ECM",
    "Rrna Processing" = "rRNA Processing"
  )
)

source("04_Figures_v2/shared/comparison_panels/panel_A_ora.R")
