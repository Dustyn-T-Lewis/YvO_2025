# F07 — Master Orchestrator
# Sources main + supp stitchers.
#
# Run order:
#   1. main/90_stitch_main.R         -> 2-panel main composite + xlsx + cleanup
#   2. supp/90_stitch_supp.R         -> supp composite

setwd(rprojroot::find_rstudio_root_file())

cat("=== F07: Running main composite + xlsx ===\n")
source("04_Figures/F07/a_script/main/90_stitch_main.R")

cat("=== F07: Running supp composite ===\n")
source("04_Figures/F07/a_script/supp/90_stitch_supp.R")

# Final cleanup: remove any leftover CSVs
remaining <- list.files("04_Figures/F07/c_data", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}

# --- Copy to Box manuscript directory ---
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- "04_Figures/F07/b_reports"
file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"),
          file.path(BOX, "02_Figures/F07_phenotype_prediction/main/pdf/MAIN_F07_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"),
          file.path(BOX, "02_Figures/F07_phenotype_prediction/main/png/MAIN_F07_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F07_composite_main.pdf"),
          file.path(BOX, "02_Figures/F07_phenotype_prediction/supp/pdf/SUPP_F07_composite_main.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F07_composite_loso.pdf"),
          file.path(BOX, "02_Figures/F07_phenotype_prediction/supp/pdf/SUPP_F07_composite_loso.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F07_composite.png"),
          file.path(BOX, "02_Figures/F07_phenotype_prediction/supp/png/SUPP_F07_composite.png"), overwrite = TRUE)
file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx",
          file.path(BOX, "02_Figures/F07_phenotype_prediction/F07_source_data.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx",
          file.path(BOX, "03_Supplementary_Tables/S13_F07_phenotype_pred.xlsx"), overwrite = TRUE)

cat("=== F07 complete ===\n")
