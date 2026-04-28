#!/usr/bin/env Rscript
# F00 — Pipeline QC: Orchestrator
# Sources the main composite script, copies outputs to Box

source(here::here("04_Figures", "F00", "a_script", "01_main_panels.R"))

BOX <- file.path("/Users/dtl0018/Library/CloudStorage/Box-Box",
                 "YvO_proteomics_manuscript")
if (dir.exists(BOX)) {
  RPT <- here::here("04_Figures", "F00", "b_reports")
  DAT <- here::here("04_Figures", "F00", "c_data")
  box_qc <- file.path(BOX, "02_Figures", "F00_pipeline_QC")
  dir.create(file.path(box_qc, "main"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main", "pdf", "MAIN_F00_QC_composite.pdf"),
            file.path(box_qc, "main", "MAIN_F00_QC_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main", "png", "MAIN_F00_QC_page1.png"),
            file.path(box_qc, "main", "MAIN_F00_QC_page1.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "main", "png", "MAIN_F00_QC_page2.png"),
            file.path(box_qc, "main", "MAIN_F00_QC_page2.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F00_supplementary.xlsx"),
            file.path(box_qc, "F00_source_data.xlsx"), overwrite = TRUE)
  message("Copied F00 outputs to Box")
}

message("F00 complete")
