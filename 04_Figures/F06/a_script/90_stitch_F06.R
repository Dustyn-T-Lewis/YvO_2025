# F06 — Master Orchestrator
#
# NOTE: YvO_WGCNA_run.R is NOT sourced here — run it separately first.
#
# Run order:
#   1. 02_supp_panels.R  -> QC composite + module triptychs/hubs/chord + preservation
#   2. 01_main_panels.R  -> 2-panel main composite + F06_supplementary.xlsx
#
# Supp runs before main because preservation CSV (from _supp_preservation.R)
# is consumed by Panel A heatmap. Main runs last to build the xlsx after all
# panel CSVs are generated.

library(here)

message("=== F06: Running supp panels (QC + modules + preservation) ===")
source(here::here("04_Figures", "F06", "a_script", "02_supp_panels.R"))

message("=== F06: Running main composite + xlsx ===")
source(here::here("04_Figures", "F06", "a_script", "01_main_panels.R"))

# --- Copy to Box manuscript directory ---
BASE <- here::here("04_Figures", "F06")
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- file.path(BASE, "b_reports")

file.copy(file.path(RPT, "main/pdf/MAIN_F06_composite.pdf"),
          file.path(BOX, "02_Figures/F06_WGCNA/main/pdf/MAIN_F06_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F06_composite.png"),
          file.path(BOX, "02_Figures/F06_WGCNA/main/png/MAIN_F06_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F06_composite.pdf"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/pdf/SUPP_F06_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F06_composite.png"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/png/SUPP_F06_composite.png"), overwrite = TRUE)
# Copy individual module panels
mod_pdfs <- list.files(file.path(RPT, "supp/pdf/modules"),
                       pattern = "^SUPP_F06_module_.*\\.pdf$", full.names = TRUE)
mod_pngs <- list.files(file.path(RPT, "supp/png/modules"),
                       pattern = "^SUPP_F06_module_.*\\.png$", full.names = TRUE)
box_supp_pdf <- file.path(BOX, "02_Figures/F06_WGCNA/supp/pdf")
box_supp_png <- file.path(BOX, "02_Figures/F06_WGCNA/supp/png")
dir.create(box_supp_pdf, recursive = TRUE, showWarnings = FALSE)
dir.create(box_supp_png, recursive = TRUE, showWarnings = FALSE)
for (f in mod_pdfs) file.copy(f, file.path(box_supp_pdf, basename(f)), overwrite = TRUE)
for (f in mod_pngs) file.copy(f, file.path(box_supp_png, basename(f)), overwrite = TRUE)
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "02_Figures/F06_WGCNA/F06_source_data.xlsx"), overwrite = TRUE)
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "03_Supplementary_Tables/S12_F06_WGCNA.xlsx"), overwrite = TRUE)

message("=== F06 complete ===")
