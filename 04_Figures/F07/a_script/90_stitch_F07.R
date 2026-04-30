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
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"),
            file.path(box_pdf, "MAIN_F07_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"),
            file.path(box_png, "MAIN_F07_composite.png"), overwrite = TRUE)
  # S9 Figure (composite_main; LOSO sensitivity is reference-only)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F07_composite_main.pdf"),
            file.path(box_fig_pdf, "S09_Figure_F07.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F07_composite.png"),
            file.path(box_fig_png, "S09_Figure_F07.png"), overwrite = TRUE)
  # S11 Table
  file.copy(file.path(DAT, "F07_supplementary.xlsx"),
            file.path(box_tbl, "S11_Table_F07.xlsx"), overwrite = TRUE)
  message("Copied F07 outputs to Box")
}

message("=== F07 complete ===")
