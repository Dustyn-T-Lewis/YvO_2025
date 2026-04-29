#!/usr/bin/env Rscript
# F00 — Pipeline QC: Orchestrator
# Sources the supp composite script, copies outputs to Box
# F00 has no main figure — produces 2 supp composites:
#   SUPP_F00_normalization (panels A–G)
#   SUPP_F00_imputation    (panels H–N)

source(here::here("04_Figures", "F00", "a_script", "01_main_panels.R"))

BOX <- file.path("/Users/dtl0018/Library/CloudStorage/Box-Box",
                 "YvO_proteomics_manuscript")
if (dir.exists(BOX)) {
  RPT <- here::here("04_Figures", "F00", "b_reports")
  DAT <- here::here("04_Figures", "F00", "c_data")
  box_qc <- file.path(BOX, "02_Figures", "F00_pipeline_QC")
  dir.create(file.path(box_qc, "supp"), recursive = TRUE, showWarnings = FALSE)
  # Supp composites (2 figures)
  file.copy(file.path(RPT, "supp", "pdf", "SUPP_F00_normalization.pdf"),
            file.path(box_qc, "supp", "S1a_Figure.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "png", "SUPP_F00_normalization.png"),
            file.path(box_qc, "supp", "S1a_Figure.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "pdf", "SUPP_F00_imputation.pdf"),
            file.path(box_qc, "supp", "S1b_Figure.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "png", "SUPP_F00_imputation.png"),
            file.path(box_qc, "supp", "S1b_Figure.png"), overwrite = TRUE)
  # Source data
  file.copy(file.path(DAT, "F00_supplementary.xlsx"),
            file.path(box_qc, "F00_source_data.xlsx"), overwrite = TRUE)
  file.copy(file.path(DAT, "F00_supplementary.xlsx"),
            file.path(BOX, "Supplementary/S04_pipeline_QC.xlsx"), overwrite = TRUE)
  message("Copied F00 outputs to Box")
}

message("F00 complete")
