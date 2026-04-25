# F06 — Master Orchestrator
# Sources WGCNA runner is NOT included (run separately).
# Sources main panels + supp panels, then stitchers.
#
# Run order:
#   1. main/90_stitch_main.R                    -> 2-panel main composite + xlsx
#   2. supp/01_QC panel scripts                 -> individual QC panels
#   3. supp/03_module/panel_B_triptych.R        -> per-module triptych PNGs (all 10)
#   4. supp/03_module/panel_D_hub.R             -> per-module hub network PNGs (all 10)
#   5. supp/90_stitch_mega.R                    -> QC composite + module composite

setwd(rprojroot::find_rstudio_root_file())

cat("=== F06: Running main composite + xlsx ===\n")
source("04_Figures/F06/a_script/main/90_stitch_main.R")

cat("=== F06: Running supp QC panel scripts ===\n")
source("04_Figures/F06/a_script/supp/01_QC/a05_soft_threshold_ggplot.R")
source("04_Figures/F06/a_script/supp/01_QC/a01_dendrogram.R")
source("04_Figures/F06/a_script/supp/01_QC/a03_compartment_enrichment.R")
source("04_Figures/F06/a_script/supp/01_QC/a02_bicor_sensitivity.R")

cat("=== F06: Running supp module triptychs (all 10 modules) ===\n")
source("04_Figures/F06/a_script/supp/03_module/panel_B_triptych.R")

cat("=== F06: Running supp module hub networks (all 10 modules) ===\n")
source("04_Figures/F06/a_script/supp/03_module/panel_D_hub.R")

cat("=== F06: Running supp composites ===\n")
source("04_Figures/F06/a_script/supp/90_stitch_mega.R")

# --- Copy to Box manuscript directory ---
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- "04_Figures/F06/b_reports"
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
file.copy("04_Figures/F06/c_data/F06_supplementary.xlsx",
          file.path(BOX, "02_Figures/F06_WGCNA/F06_source_data.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F06/c_data/F06_supplementary.xlsx",
          file.path(BOX, "03_Supplementary_Tables/S12_F06_WGCNA.xlsx"), overwrite = TRUE)

cat("=== F06 complete ===\n")
