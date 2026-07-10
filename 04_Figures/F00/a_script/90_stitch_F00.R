#!/usr/bin/env Rscript
# F00 — Pipeline QC: Orchestrator
# Sources the supp composite script, copies outputs to Box
# F00 has no main figure — produces 2 supp composites:
#   SUPP_F00_normalization (panels A–G)
#   SUPP_F00_imputation    (panels H–N)

setwd(here::here())

source("04_Figures/F00/a_script/01_supp_panels.R")

BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
if (nzchar(BOX) && dir.exists(BOX)) {
  RPT <- "04_Figures/F00/b_reports"
  DAT <- "04_Figures/F00/c_data"
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "supp", "pdf", "SUPP_F00_normalization.pdf"),
            file.path(box_fig_pdf, "S01_Figure_Normalization.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "png", "SUPP_F00_normalization.png"),
            file.path(box_fig_png, "S01_Figure_Normalization.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "pdf", "SUPP_F00_imputation.pdf"),
            file.path(box_fig_pdf, "S02_Figure_Imputation.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp", "png", "SUPP_F00_imputation.png"),
            file.path(box_fig_png, "S02_Figure_Imputation.png"), overwrite = TRUE)
  file.copy(file.path(DAT, "F00_supplementary.xlsx"),
            file.path(box_tbl, "S04_Table_F00.xlsx"), overwrite = TRUE)
  message("Copied F00 outputs to Box")
}

message("F00 complete")
