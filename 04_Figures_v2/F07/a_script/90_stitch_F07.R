#!/usr/bin/env Rscript
# F07 — Phenotype Prediction: Master Orchestrator
# Sources main panels (all data-generating scripts + composite + xlsx),
# then supp panels (reads pre-rendered PNGs, builds composites).

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures_v2/shared/figure_supplement_helpers.R")

DAT <- "04_Figures_v2/F07/c_data"

message("=== F07: Running main panels + composite + xlsx ===")
source("04_Figures_v2/F07/a_script/01_main_panels.R")

message("=== F07: Running supp composite ===")
source("04_Figures_v2/F07/a_script/02_supp_panels.R")

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT,
  pattern = "\\.csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("F07 complete: composites in 04_Figures_v2/F07/b_reports, table in ", DAT)
