#!/usr/bin/env Rscript
# F02 — Proteome + DEP Overview: Master Orchestrator

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures_v2/shared/figure_supplement_helpers.R")

DAT <- "04_Figures_v2/F02/c_data"

# Supp panels first (CSVs needed for xlsx)
source("04_Figures_v2/F02/a_script/02_supp_panels.R")

source("04_Figures_v2/F02/a_script/01_main_panels.R")

audit_csvs <- list.files(DAT, pattern = "^(audit_|panel_|SUPP_panel_)", full.names = TRUE)
specs <- lapply(audit_csvs, \(p) list(
  name = tools::file_path_sans_ext(basename(p)),
  path = p
))
build_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  title = "F02 \u2014 Proteome + DEP overview source data",
  description = "Main and supplementary panel source data.",
  sheet_specs = specs
)
cleanup_after_workbook(specs)

message("F02 complete: composites in 04_Figures_v2/F02/b_reports, table in ", DAT)
