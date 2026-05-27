#!/usr/bin/env Rscript
# F03 — Volcano Rings: Master Orchestrator

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F03/c_data"

source("04_Figures/F03/a_script/02_supp_panels.R")
source("04_Figures/F03/a_script/01_main_panels.R")

# Build xlsx from per-contrast xlsx sheets + ring terms + supp CSVs
DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

f03_specs <- lapply(CTRS, \(ctr) {
  list(name = paste0("MAIN_", ctr),
       df = as.data.frame(readxl::read_excel(DEP_XLSX, sheet = ctr)))
})
for (tag in c("A", "B", "C", "D")) {
  ring_path <- file.path(DAT, paste0("panel_", tag), "ring_terms.csv")
  if (file.exists(ring_path))
    f03_specs <- c(f03_specs, list(list(name = paste0("RING_", tag), path = ring_path)))
}
supp_csvs <- list.files(file.path(DAT, "supp"), pattern = "\\.csv$", full.names = TRUE)
for (p in supp_csvs)
  f03_specs <- c(f03_specs, list(list(name = paste0("SUPP_", tools::file_path_sans_ext(basename(p))),
                                       path = p)))

build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  title = "F03 \u2014 Volcano ring source data",
  description = "Per-contrast DEP tables, ring terms, and supp diagnostics.",
  overview_df = data.frame(Sheet = sapply(f03_specs, `[[`, "name"),
                           Description = sapply(f03_specs, `[[`, "name")),
  sheet_specs = f03_specs)
cleanup_after_workbook(f03_specs,
  extra_subdirs = c(file.path(DAT, "panel_A"), file.path(DAT, "panel_B"),
                    file.path(DAT, "panel_C"), file.path(DAT, "panel_D"),
                    file.path(DAT, "supp")))

BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
if (nzchar(BOX) && dir.exists(BOX)) {
  RPT <- "04_Figures/F03/b_reports"
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F03_composite.pdf"),
            file.path(box_pdf, "MAIN_F03_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F03_composite.png"),
            file.path(box_png, "MAIN_F03_composite.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F03_composite.pdf"),
            file.path(box_fig_pdf, "S05_Figure_F03.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F03_composite.png"),
            file.path(box_fig_png, "S05_Figure_F03.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F03_supplementary.xlsx"),
            file.path(box_tbl, "S07_Table_F03.xlsx"), overwrite = TRUE)
  message("Copied F03 outputs to Box")
}

message("F03 complete")
