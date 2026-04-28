#!/usr/bin/env Rscript
# F07 — Phenotype Prediction: Master Orchestrator
# Sources main panels (all data-generating scripts + composite + xlsx),
# then supp panels (reads pre-rendered PNGs, builds composites).

source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

DAT <- here::here("04_Figures", "F07", "c_data")

message("=== F07: Running main panels + composite + xlsx ===")
source(here::here("04_Figures", "F07", "a_script", "01_main_panels.R"))

message("=== F07: Running supp composite ===")
source(here::here("04_Figures", "F07", "a_script", "02_supp_panels.R"))

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

# ── Copy to Box manuscript directory ─────────────────────────────────────────
BOX <- file.path("/Users/dtl0018/Library/CloudStorage/Box-Box",
                 "YvO_proteomics_manuscript")
if (dir.exists(BOX)) {
  RPT <- here::here("04_Figures", "F07", "b_reports")
  box_f07 <- file.path(BOX, "02_Figures", "F07_phenotype_prediction")
  for (d in file.path(box_f07, c("main", "supp")))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"),
            file.path(box_f07, "main/MAIN_F07_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"),
            file.path(box_f07, "main/MAIN_F07_composite.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F07_composite_main.pdf"),
            file.path(box_f07, "supp/SUPP_F07_composite_main.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F07_composite_loso.pdf"),
            file.path(box_f07, "supp/SUPP_F07_composite_loso.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F07_composite.png"),
            file.path(box_f07, "supp/SUPP_F07_composite.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F07_supplementary.xlsx"),
            file.path(BOX, "02_Figures/F07_phenotype_prediction/F07_source_data.xlsx"), overwrite = TRUE)
  file.copy(file.path(DAT, "F07_supplementary.xlsx"),
            file.path(BOX, "03_Supplementary_Tables/S13_F07_phenotype_pred.xlsx"), overwrite = TRUE)
  message("Copied F07 outputs to Box")
}

message("=== F07 complete ===")
