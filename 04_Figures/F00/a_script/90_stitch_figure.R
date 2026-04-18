# F00 — Pipeline Quality Control: Master Orchestrator
# Sources all 4 sub-stitchers (main overview + 3 stage detail composites),
# then builds the overall combined composite at b_reports/ root.
#
# Run order:
#   1. main/90_stitch_main.R         -> 6-panel overview (A-F)
#   2. 01_normalization/90_stitch_stage01.R -> 4-panel normalization detail
#   3. 02_imputation/90_stitch_stage02.R    -> 4-panel imputation detail
#   4. 03_DEP/90_stitch_stage03.R           -> 4-panel DEP detail
#   5. Overall composite (this file)        -> combined at b_reports/ root
#
# Output: 04_Figures/F00/b_reports/main/{pdf,png}/MAIN_F00_QC_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())

# --- Run all sub-stitchers ---
cat("=== F00: Running main overview ===\n")
source("04_Figures/F00/a_script/main/90_stitch_main.R")

cat("=== F00: Running stage 01 (normalization) ===\n")
source("04_Figures/F00/a_script/01_normalization/90_stitch_stage01.R")

cat("=== F00: Running stage 02 (imputation) ===\n")
source("04_Figures/F00/a_script/02_imputation/90_stitch_stage02.R")

cat("=== F00: Running stage 03 (DEP) ===\n")
source("04_Figures/F00/a_script/03_DEP/90_stitch_stage03.R")

cat("=== F00 complete: all composites + stage excels built ===\n")
