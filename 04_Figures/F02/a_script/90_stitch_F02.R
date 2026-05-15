#!/usr/bin/env Rscript
# F02 — Proteome + DEP Overview: Master Orchestrator

setwd(rprojroot::find_rstudio_root_file())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F02/c_data"

# Supp panels first (CSVs needed for xlsx)
source("04_Figures/F02/a_script/02_supp_panels.R")

# Main panels + composite
source("04_Figures/F02/a_script/01_main_panels.R")

# Build xlsx
audit_csvs <- list.files(DAT, pattern = "^(audit_|panel_|SUPP_panel_)", full.names = TRUE)
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
BOX <- Sys.getenv("YVO_BOX_DIR", unset = file.path(
  "/Users/dtl0018/Library/CloudStorage/Box-Box",
  "YvO_proteomics_manuscript"))
if (dir.exists(BOX)) {
  RPT <- "04_Figures/F02/b_reports"
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F02_composite.pdf"),
            file.path(box_pdf, "MAIN_F02_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F02_composite.png"),
            file.path(box_png, "MAIN_F02_composite.png"), overwrite = TRUE)
  # S4 Figure
  file.copy(file.path(RPT, "supp/pdf/SUPP_F02_composite.pdf"),
            file.path(box_fig_pdf, "S04_Figure_F02.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F02_composite.png"),
            file.path(box_fig_png, "S04_Figure_F02.png"), overwrite = TRUE)
  # S6 Table
  file.copy(file.path(DAT, "F02_supplementary.xlsx"),
            file.path(box_tbl, "S06_Table_F02.xlsx"), overwrite = TRUE)
  message("Copied F02 outputs to Box")
}

message("F02 complete")
