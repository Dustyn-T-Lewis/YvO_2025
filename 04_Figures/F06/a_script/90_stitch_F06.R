# F06 master orchestrator: run YvO_WGCNA_run.R separately first.


setwd(here::here())

source("04_Figures/F06/a_script/01_main_panels.R")
source("04_Figures/F06/a_script/02_supp_panels.R")

BASE <- "04_Figures/F06"
BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
RPT <- file.path(BASE, "b_reports")

if (nzchar(BOX) && dir.exists(BOX)) {
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F06_composite.pdf"),
            file.path(box_pdf, "MAIN_F06_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F06_composite.png"),
            file.path(box_png, "MAIN_F06_composite.png"), overwrite = TRUE)
  # S8 Figure (network diagnostics + module triptychs)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F06_composite.pdf"),
            file.path(box_fig_pdf, "S08_Figure_F06.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F06_composite.png"),
            file.path(box_fig_png, "S08_Figure_F06.png"), overwrite = TRUE)
  # S8 Figure continued: module triptychs (individual pages)
  triptych_pdfs <- list.files(file.path(RPT, "supp/pdf/modules"),
                              pattern = "^SUPP_triptych_", full.names = TRUE)
  for (tp in triptych_pdfs)
    file.copy(tp, file.path(box_fig_pdf, basename(tp)), overwrite = TRUE)
  # S10 Table
  file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
            file.path(box_tbl, "S10_Table_F06.xlsx"), overwrite = TRUE)
  message("Copied F06 outputs to Box")
}

message("F06 complete")
