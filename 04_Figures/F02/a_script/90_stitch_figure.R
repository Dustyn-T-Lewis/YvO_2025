# F02 — Master Orchestrator
# Sources supp panels first (so CSVs exist for xlsx), then main + supp stitchers.
#
# Run order:
#   1. supp/panels/*.R               -> supp panel CSVs (pre-generate for xlsx)
#   2. main/90_stitch_main.R         -> 6-panel main composite + xlsx + cleanup
#   3. supp/90_stitch_supp.R         -> 3-panel supp composite
#   4. Final cleanup                 -> remove any leftover CSVs from supp re-source

setwd(rprojroot::find_rstudio_root_file())

cat("=== F02: Pre-generating supp panel data ===\n")
source("04_Figures/F02/a_script/supp/panels/panel_C.R")
source("04_Figures/F02/a_script/supp/panels/panel_D.R")
source("04_Figures/F02/a_script/supp/panels/panel_E.R")

cat("=== F02: Running main composite + xlsx ===\n")
source("04_Figures/F02/a_script/main/90_stitch_main.R")

cat("=== F02: Running supp composite ===\n")
source("04_Figures/F02/a_script/supp/90_stitch_supp.R")

# --- Copy to Box manuscript directory ---
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- "04_Figures/F02/b_reports"
file.copy(file.path(RPT, "main/pdf/MAIN_F02_composite.pdf"),
          file.path(BOX, "02_Figures/F02_proteome_DEP/main/pdf/MAIN_F02_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F02_composite.png"),
          file.path(BOX, "02_Figures/F02_proteome_DEP/main/png/MAIN_F02_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F02_composite.pdf"),
          file.path(BOX, "02_Figures/F02_proteome_DEP/supp/pdf/SUPP_F02_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F02_composite.png"),
          file.path(BOX, "02_Figures/F02_proteome_DEP/supp/png/SUPP_F02_composite.png"), overwrite = TRUE)
file.copy("04_Figures/F02/c_data/F02_supplementary.xlsx",
          file.path(BOX, "02_Figures/F02_proteome_DEP/F02_source_data.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F02/c_data/F02_supplementary.xlsx",
          file.path(BOX, "03_Supplementary_Tables/S08_F02_proteome_DEP.xlsx"), overwrite = TRUE)

# Final cleanup: supp stitcher re-sources panels, which re-writes CSVs
remaining <- list.files("04_Figures/F02/c_data", pattern = "\\.csv$", full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}

cat("=== F02 complete ===\n")
