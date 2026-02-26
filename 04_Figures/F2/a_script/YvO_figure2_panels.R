################################################################################
#   Figure 2 — Panel Export Pipeline (Orchestrator)
#   Sources shared setup once, then runs each per-panel script in order.
#
#   Individual panel scripts can also be run standalone (they source shared
#   setup themselves if dep_df is not already loaded).
#
#   Panels:
#     A — Volcano: Training Effect (Young)
#     B — Volcano: Training Effect (Old)
#     C — Concordance scatter (logFC x logFC)
#     D — RRHO2 threshold-free concordance map
#     E — fGSEA NES scatter (Hallmark + rrvgo-reduced GO:BP)
#     F — Interaction DEPs: Multi-Contrast Response & Pathway Enrichment
#         (merged dumbbell | sankey | bar panel)
#
#   Outputs per panel:
#     b_reports/panel_A_volcano.pdf         + c_data/panel_A/volcano_young.csv
#     b_reports/panel_B_volcano.pdf         + c_data/panel_B/volcano_old.csv
#     b_reports/panel_C_concordance.pdf     + c_data/panel_C/concordance.csv
#     b_reports/panel_D_rrho2.pdf           + c_data/panel_D/rrho2_*.csv
#     b_reports/panel_E_nes_bubble.pdf      + c_data/panel_E/nes_scatter.csv
#     b_reports/panel_F_interaction.pdf     + c_data/panel_F/*.csv
################################################################################

# === Shared setup (packages, style, helpers, data) ===========================
source("04_Figures/F2/a_script/YvO_figure2_shared.R")

# === Panel scripts (each checks for dep_df before re-sourcing shared) ========
source("04_Figures/F2/a_script/YvO_panel_AB.R")
source("04_Figures/F2/a_script/YvO_panel_C.R")
source("04_Figures/F2/a_script/YvO_panel_D.R")
source("04_Figures/F2/a_script/YvO_panel_E.R")
source("04_Figures/F2/a_script/YvO_panel_F.R")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n", strrep("=", 61), "\nFigure 2 Panel Export Complete\n", strrep("=", 61), "\n")
cat("\nPDFs:\n")
for (f in c("panel_A_volcano", "panel_B_volcano", "panel_C_concordance",
            "panel_D_rrho2", "panel_E_nes_bubble", "panel_F_interaction"))
  cat(sprintf("  %s/%s.pdf\n", RPT_DIR, f))
cat("\nData (organized by panel):\n")
csv_map <- list(
  panel_A = "volcano_young.csv",
  panel_B = "volcano_old.csv",
  panel_C = "concordance.csv",
  panel_D = c("rrho2_summary.csv", "rrho2_matrix.csv"),
  panel_E = "nes_scatter.csv",
  panel_F = c("interaction_classification.csv", "interaction_dot.csv",
              "interaction_dot_long.csv", "sankey_dot.csv", "sankey_links.csv",
              "interaction_patterns.csv", "srplot_input.csv"),
  shared  = "fgsea_tstat_all_v2.csv"
)
for (pnl in names(csv_map))
  for (f in csv_map[[pnl]])
    cat(sprintf("  %s/%s/%s\n", DAT_DIR, pnl, f))
