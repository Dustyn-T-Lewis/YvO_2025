# F01 — Phenotype Figure: Master Orchestrator
# Sources supp panels first (so CSVs exist for xlsx), then main + supp stitchers.
#
# Run order:
#   1. supp/panels/*.R                          -> supp panel CSVs (pre-generate for xlsx)
#   2. main/90_stitch_main_single_column.R      -> single-column composite + xlsx + cleanup
#   3. main/90_stitch_main_double_column.R      -> double-column composite
#   4. supp/90_stitch_supp.R                    -> 3-panel supp composite
#   5. Final cleanup                            -> remove any leftover CSVs

setwd(rprojroot::find_rstudio_root_file())

cat("=== F01: Pre-generating supp panel data ===\n")
source("04_Figures/F01/a_script/supp/panels/panel_A.R")
source("04_Figures/F01/a_script/supp/panels/panel_B.R")
source("04_Figures/F01/a_script/supp/panels/panel_C.R")

cat("=== F01: Running main single-column composite + xlsx ===\n")
source("04_Figures/F01/a_script/main/90_stitch_main_single_column.R")

cat("=== F01: Running main double-column composite ===\n")
source("04_Figures/F01/a_script/main/90_stitch_main_double_column.R")

cat("=== F01: Running supp composite ===\n")
source("04_Figures/F01/a_script/supp/90_stitch_supp.R")

# --- Copy to Box manuscript directory ---
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- "04_Figures/F01/b_reports"
file.copy(file.path(RPT, "main/pdf/MAIN_F01_composite.pdf"),
          file.path(BOX, "02_Figures/F01_phenotype/main/pdf/MAIN_F01_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F01_composite.png"),
          file.path(BOX, "02_Figures/F01_phenotype/main/png/MAIN_F01_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "main/pdf/MAIN_F01_composite_double_col.pdf"),
          file.path(BOX, "02_Figures/F01_phenotype/main/pdf/MAIN_F01_composite_double_col.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F01_composite_double_col.png"),
          file.path(BOX, "02_Figures/F01_phenotype/main/png/MAIN_F01_composite_double_col.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F01_composite.pdf"),
          file.path(BOX, "02_Figures/F01_phenotype/supp/pdf/SUPP_F01_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F01_composite.png"),
          file.path(BOX, "02_Figures/F01_phenotype/supp/png/SUPP_F01_composite.png"), overwrite = TRUE)
file.copy("04_Figures/F01/c_data/F01_supplementary.xlsx",
          file.path(BOX, "02_Figures/F01_phenotype/F01_source_data.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F01/c_data/F01_supplementary.xlsx",
          file.path(BOX, "03_Supplementary_Tables/S07_F01_phenotype.xlsx"), overwrite = TRUE)

# Final cleanup: supp stitcher re-sources panels, which re-writes CSVs
remaining <- list.files("04_Figures/F01/c_data", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}
supp_dir <- "04_Figures/F01/c_data/supp"
if (dir.exists(supp_dir) && length(list.files(supp_dir)) == 0)
  unlink(supp_dir, recursive = TRUE)

cat("=== F01 complete ===\n")
