#!/usr/bin/env Rscript
# F02 — Proteome + DEP Overview: Master Orchestrator

source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

DAT <- here::here("04_Figures", "F02", "c_data")

# Supp panels first (CSVs needed for xlsx)
source(here::here("04_Figures", "F02", "a_script", "02_supp_panels.R"))

# Main panels + composite
source(here::here("04_Figures", "F02", "a_script", "01_main_panels.R"))

# Build xlsx
audit_csvs <- list.files(DAT, pattern = "^(audit_|panel_)", full.names = TRUE)
specs <- lapply(audit_csvs, \(p) list(name = tools::file_path_sans_ext(basename(p)),
                                       path = p))
build_workbook(
  file.path(DAT, "F02_supplementary.xlsx"),
  title = "F02 \u2014 Proteome + DEP overview source data",
  description = "Main and supplementary panel source data.",
  overview_df = data.frame(Sheet = sapply(specs, `[[`, "name"),
                           Description = sapply(specs, `[[`, "name")),
  sheet_specs = specs)
cleanup_after_workbook(specs)

# Box copy
BOX <- file.path("/Users/dtl0018/Library/CloudStorage/Box-Box",
                 "YvO_proteomics_manuscript")
if (dir.exists(BOX)) {
  RPT <- here::here("04_Figures", "F02", "b_reports")
  box_f02 <- file.path(BOX, "02_Figures", "F02_proteome_DEP")
  for (d in file.path(box_f02, c("main", "supp")))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F02_composite.pdf"),
            file.path(box_f02, "main/MAIN_F02_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F02_composite.png"),
            file.path(box_f02, "main/MAIN_F02_composite.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F02_composite.pdf"),
            file.path(box_f02, "supp/SUPP_F02_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F02_composite.png"),
            file.path(box_f02, "supp/SUPP_F02_composite.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F02_supplementary.xlsx"),
            file.path(BOX, "03_Supplementary_Tables/S08_F02_proteome_DEP.xlsx"),
            overwrite = TRUE)
  message("Copied F02 outputs to Box")
}

message("F02 complete")
