#!/usr/bin/env Rscript
# F01 — Phenotype Figure: Master Orchestrator
# Supp panels first (CSVs needed for xlsx), then main (builds xlsx), then supp composite

source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

DAT <- here::here("04_Figures", "F01", "c_data")

# Generate supp CSVs first so xlsx includes them
source(here::here("04_Figures", "F01", "a_script", "02_supp_panels.R"))

# Main panels + both composites + xlsx
source(here::here("04_Figures", "F01", "a_script", "01_main_panels.R"))

# Build xlsx with all panel data
f01_specs <- list(
  list(name = "panel_A",      path = file.path(DAT, "panel_A_training_volume.csv")),
  list(name = "panel_B",      path = file.path(DAT, "panel_B_dxa_lbm.csv")),
  list(name = "panel_C",      path = file.path(DAT, "panel_C_vl_thickness.csv")),
  list(name = "SUPP_panel_A", path = file.path(DAT, "supp", "panel_A_deadlift_1rm.csv")),
  list(name = "SUPP_panel_B", path = file.path(DAT, "supp", "panel_B_type_II_fcsa.csv")),
  list(name = "SUPP_panel_C", path = file.path(DAT, "supp", "panel_C_type_I_fcsa.csv")))

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  title = "F01 \u2014 Phenotype source data",
  description = "Main panels A\u2013C and supplementary panels A\u2013C.",
  overview_df = data.frame(
    Sheet = sapply(f01_specs, `[[`, "name"),
    Description = c("Training volume", "DXA LBM", "VL thickness",
                    "Deadlift 1RM", "Type II fCSA", "Type I fCSA")),
  sheet_specs = f01_specs)
cleanup_after_workbook(f01_specs,
                       extra_subdirs = file.path(DAT, "supp"))

# Box copy
BOX <- file.path("/Users/dtl0018/Library/CloudStorage/Box-Box",
                 "YvO_proteomics_manuscript")
if (dir.exists(BOX)) {
  RPT <- here::here("04_Figures", "F01", "b_reports")
  box_f01 <- file.path(BOX, "02_Figures", "F01_phenotype")
  for (d in file.path(box_f01, c("main", "supp")))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F01_composite.pdf"),
            file.path(box_f01, "main/MAIN_F01_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F01_composite.png"),
            file.path(box_f01, "main/MAIN_F01_composite.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F01_composite.pdf"),
            file.path(box_f01, "supp/SUPP_F01_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F01_composite.png"),
            file.path(box_f01, "supp/SUPP_F01_composite.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F01_supplementary.xlsx"),
            file.path(BOX, "Supplementary/S05_F01_phenotype.xlsx"), overwrite = TRUE)
  message("Copied F01 outputs to Box")
}

message("F01 complete")
