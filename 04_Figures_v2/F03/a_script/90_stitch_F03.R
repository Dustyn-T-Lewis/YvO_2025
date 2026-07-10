#!/usr/bin/env Rscript
# F03 — Volcano Rings: Master Orchestrator


source(here::here("04_Figures_v2", "shared_functions", "shared_figure_supplement_helpers.R"))

DAT <- here::here("04_Figures_v2", "F03", "c_data")

source(here::here("04_Figures_v2", "F03", "a_script", "02_supp_panels.R"))
source(here::here("04_Figures_v2", "F03", "a_script", "01_main_panels.R"))

# Build xlsx from per-contrast xlsx sheets + ring terms + supp CSVs
DEP_XLSX <- here::here("03_DEP", "c_data", "03_DEP_results.xlsx")
CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

f03_specs <- lapply(CTRS, \(ctr) {
  list(
    name = paste0("MAIN_", ctr),
    df = as.data.frame(readxl::read_excel(DEP_XLSX, sheet = ctr))
  )
})
for (tag in c("A", "B", "C", "D")) {
  ring_path <- file.path(DAT, paste0("panel_", tag), "ring_terms.csv")
  if (file.exists(ring_path)) {
    f03_specs <- c(f03_specs, list(list(name = paste0("RING_", tag), path = ring_path)))
  }
}
supp_csvs <- list.files(file.path(DAT, "supp"), pattern = "\\.csv$", full.names = TRUE)
for (p in supp_csvs) {
  f03_specs <- c(f03_specs, list(list(
    name = paste0("SUPP_", tools::file_path_sans_ext(basename(p))),
    path = p
  )))
}

build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  title = "F03 \u2014 Volcano ring source data",
  description = "Per-contrast DEP tables, ring terms, and supp diagnostics.",
  sheet_specs = f03_specs
)
cleanup_after_workbook(f03_specs,
  extra_subdirs = c(
    file.path(DAT, "panel_A"), file.path(DAT, "panel_B"),
    file.path(DAT, "panel_C"), file.path(DAT, "panel_D"),
    file.path(DAT, "supp")
  )
)

message("F03 complete: composites in 04_Figures_v2/F03/b_reports, table in ", DAT)
