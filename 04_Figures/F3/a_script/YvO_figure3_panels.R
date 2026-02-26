################################################################################
#   Figure 3 — Panel Export Pipeline (Orchestrator)
#   Sources shared setup once, then runs each per-panel script in order.
#
#   Individual panel scripts can also be run standalone.
#
#   Panels:
#     A — Volcano Ring: Aging Effect
#     B — Volcano Ring: Reversal Effect
#     C — Reversal scatter (logFC Aging vs logFC Training Old)
#     D — RRHO2 reversal map
#     E — NES scatter (pathway-level reversal)
#     F — Reversal classification (dumbbell + sankey + enrichment bars)
#
#   Outputs per panel:
#     b_reports/panel_A_volcano.pdf         + c_data/panel_A/*.csv
#     b_reports/panel_B_volcano.pdf         + c_data/panel_B/*.csv
#     b_reports/panel_C_reversal_scatter.pdf + c_data/panel_C/reversal_scatter.csv
#     b_reports/panel_D_rrho2.pdf           + c_data/panel_D/rrho2_*.csv
#     b_reports/panel_E_nes_scatter.pdf     + c_data/panel_E/nes_scatter.csv
#     b_reports/panel_F_classification.pdf  + c_data/panel_F/*.csv
################################################################################

# === Shared setup (packages, style, helpers, data) ===========================
source("04_Figures/F3/a_script/YvO_figure3_shared.R")

# === Panel scripts ============================================================
source("04_Figures/F3/a_script/YvO_panel_AB.R")
source("04_Figures/F3/a_script/YvO_panel_C.R")
source("04_Figures/F3/a_script/YvO_panel_D.R")
source("04_Figures/F3/a_script/YvO_panel_E.R")
source("04_Figures/F3/a_script/YvO_panel_F.R")

# === Summary ==================================================================
cat("\n", strrep("=", 61), "\nFigure 3 Panel Export Complete\n", strrep("=", 61), "\n")
cat("\nPDFs:\n")
for (f in c("panel_A_volcano", "panel_B_volcano", "panel_C_reversal_scatter",
            "panel_D_rrho2", "panel_E_nes_scatter", "panel_F_classification"))
  cat(sprintf("  %s/%s.pdf\n", RPT_DIR, f))
cat("\nData (organized by panel):\n")
csv_map <- list(
  panel_A = c("ring_terms.csv", "volcano_aging.csv"),
  panel_B = c("ring_terms.csv", "volcano_reversal.csv"),
  panel_C = "reversal_scatter.csv",
  panel_D = c("rrho2_summary.csv", "rrho2_matrix.csv"),
  panel_E = "nes_scatter.csv",
  panel_F = c("classification.csv", "sankey_links.csv", "enrichment_bars.csv")
)
for (pnl in names(csv_map))
  for (f in csv_map[[pnl]])
    cat(sprintf("  %s/%s/%s\n", DAT_DIR, pnl, f))
