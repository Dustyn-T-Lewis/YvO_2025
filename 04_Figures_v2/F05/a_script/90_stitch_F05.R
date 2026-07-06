#!/usr/bin/env Rscript
# F05 — Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main stitcher (composite + xlsx + cleanup), then supp stitcher (composite).
#
# Run order:
#   1. _supp_enrichment_heatmap.R  -> supp panel CSVs + PNGs (pre-generate)
#   2. 02_supp_panels.R            -> supp composite from pre-rendered PNGs
#   3. 01_main_panels.R            -> 5-panel main composite
#   4. Supplementary xlsx          -> build workbook + cleanup
#   5. Final cleanup

withr::local_dir(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures_v2/shared/style.R")

BASE <- "04_Figures_v2/F05"

message("=== F05: Pre-generating supp panel data + PNGs ===")
source("04_Figures_v2/F05/a_script/_supp_enrichment_heatmap.R")

message("=== F05: Running supp diagnostics composite ===")
source("04_Figures_v2/F05/a_script/02_supp_panels.R")

message("=== F05: Running main composite ===")
source("04_Figures_v2/F05/a_script/01_main_panels.R")

# Supplementary Excel: one workbook, sheets keyed to figure panels
source("04_Figures_v2/shared/figure_supplement_helpers.R")

message("=== F05 supplementary workbook ===")
enrichment_reversal_df <- read.csv(file.path(BASE, "c_data", "panel_supp", "enrichment_reversal.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
reversal_pathway_stats_df <- read.csv(file.path(BASE, "c_data", "panel_supp", "reversal_pathway_stats.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
f05_specs <- list(
  list(name = "panel_A_ora_quadrant", path = file.path(BASE, "c_data", "panel_A", "ora_quadrant.csv")),
  list(name = "panel_B_pattern_class", path = file.path(BASE, "c_data", "panel_B_heatmap", "pattern_classification.csv")),
  list(name = "panel_B_sankey", path = file.path(BASE, "c_data", "panel_B_heatmap", "sankey_links.csv")),
  list(name = "panel_B_bar", path = file.path(BASE, "c_data", "panel_B_heatmap", "bar_data.csv")),
  list(name = "panel_C_fry_results", path = file.path(BASE, "c_data", "panel_C_fry", "fry_results_all.csv")),
  list(name = "panel_C_fry_driving", path = file.path(BASE, "c_data", "panel_C_fry", "driving_proteins.csv")),
  list(name = "panel_D_nes_scatter", path = file.path(BASE, "c_data", "panel_D", "nes_scatter.csv")),
  list(name = "panel_E_rrho2_summary", path = file.path(BASE, "c_data", "panel_E", "rrho2_summary.csv")),
  list(name = "panel_E_rrho2_hotspot", path = file.path(BASE, "c_data", "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "panel_E_rrho2_ora_concord", path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_concordant.csv")),
  list(name = "panel_E_rrho2_ora_discord", path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_discordant.csv")),
  list(name = "panel_E_rrho2_discord_note", path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_discordant_note.csv")),
  list(name = "SUPP_enrichment_reversal", df = enrichment_reversal_df),
  list(name = "SUPP_reversal_pathway_stats", df = reversal_pathway_stats_df),
  list(name = "SUPP_ora_dedup", path = file.path(BASE, "c_data", "panel_supp", "SUPP_ora_dedup_sensitivity.csv")),
  list(name = "SUPP_r_bootstrap", path = file.path(BASE, "c_data", "panel_supp", "SUPP_r_bootstrap.csv")),
  list(name = "SUPP_reversal_threshold", path = file.path(BASE, "c_data", "panel_supp", "SUPP_reversal_threshold.csv")),
  list(name = "SUPP_goslim_bars", path = file.path(BASE, "c_data", "panel_supp", "SUPP_goslim_distribution.csv")),
  list(name = "SUPP_fry_leading", path = file.path(BASE, "c_data", "panel_supp", "SUPP_fry_leading_edge.csv")),
  list(name = "SUPP_fry_circularity", path = file.path(BASE, "c_data", "panel_supp", "SUPP_fry_circularity.csv"))
)
build_workbook(
  file.path(BASE, "c_data", "F05_supplementary.xlsx"),
  title = "F05 \u2014 Figure 5 source data",
  description = "Aging-reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  sheet_specs = f05_specs
)
cleanup_after_workbook(f05_specs,
  extra_subdirs = c(
    file.path(BASE, "c_data", "panel_A"),
    file.path(BASE, "c_data", "panel_B_heatmap"),
    file.path(BASE, "c_data", "panel_C_fry"),
    file.path(BASE, "c_data", "panel_D"),
    file.path(BASE, "c_data", "panel_E"),
    file.path(BASE, "c_data", "panel_F"),
    file.path(BASE, "c_data", "panel_F_fry"),
    file.path(BASE, "c_data", "panel_G"),
    file.path(BASE, "c_data", "reversal_tests"),
    file.path(BASE, "c_data", "panel_supp"),
    file.path(BASE, "c_data", "supp")
  )
)

remaining <- list.files(file.path(BASE, "c_data"),
  pattern = "\\.csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("F05 complete: composites in 04_Figures_v2/F05/b_reports, table in ", file.path(BASE, "c_data"))
