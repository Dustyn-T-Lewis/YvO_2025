# F03 — Master Orchestrator
# Sources supp panels first (so CSVs exist for xlsx), then main + supp stitchers.
#
# Run order:
#   1. supp/panels/*.R               -> supp panel CSVs (pre-generate for xlsx)
#   2. main/90_stitch_main.R         -> 4-panel main composite + xlsx + cleanup
#   3. main/90_stitch_single_col.R   -> single-column variant
#   4. supp/90_stitch_supp.R         -> 4-panel supp composite
#   5. Final cleanup                 -> remove any leftover CSVs + empty subdirs

setwd(rprojroot::find_rstudio_root_file())

cat("=== F03: Pre-generating supp panel data ===\n")
source("04_Figures/F03/a_script/supp/panels/panel_B_phist.R")
source("04_Figures/F03/a_script/supp/panels/panel_C_ma.R")
source("04_Figures/F03/a_script/supp/panels/panel_D_sensitivity.R")
source("04_Figures/F03/a_script/supp/panels/panel_E_outlier.R")

cat("=== F03: Running main composite + xlsx ===\n")
source("04_Figures/F03/a_script/main/90_stitch_main.R")

cat("=== F03: Running single-column variant ===\n")
source("04_Figures/F03/a_script/main/90_stitch_single_col.R")

cat("=== F03: Running supp composite ===\n")
source("04_Figures/F03/a_script/supp/90_stitch_supp.R")

# --- Sync to manuscript directories ---
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F03"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F03_volcanoes"
for (d in c(file.path(WRITING_DIR, "main/pdf"), file.path(WRITING_DIR, "main/png"),
            file.path(WRITING_DIR, "supp/pdf"), file.path(WRITING_DIR, "supp/png"),
            file.path(BOX_DIR, "main/pdf"), file.path(BOX_DIR, "main/png"),
            file.path(BOX_DIR, "supp/pdf"), file.path(BOX_DIR, "supp/png")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

RPT <- "04_Figures/F03/b_reports"
# Main composites
file.copy(file.path(RPT, "main/pdf/MAIN_F03_composite.pdf"), file.path(WRITING_DIR, "main/pdf/MAIN_F03_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F03_composite.png"), file.path(WRITING_DIR, "main/png/MAIN_F03_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "main/pdf/MAIN_F03_composite.pdf"), file.path(BOX_DIR, "main/pdf/MAIN_F03_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F03_composite.png"), file.path(BOX_DIR, "main/png/MAIN_F03_composite.png"), overwrite = TRUE)
# Single-col variant
file.copy(file.path(RPT, "main/pdf/MAIN_F03_composite_single_col.pdf"), file.path(WRITING_DIR, "main/pdf/MAIN_F03_composite_single_col.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F03_composite_single_col.png"), file.path(WRITING_DIR, "main/png/MAIN_F03_composite_single_col.png"), overwrite = TRUE)
# Supp composites
file.copy(file.path(RPT, "supp/pdf/SUPP_F03_composite.pdf"), file.path(WRITING_DIR, "supp/pdf/SUPP_F03_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F03_composite.png"), file.path(WRITING_DIR, "supp/png/SUPP_F03_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F03_composite.pdf"), file.path(BOX_DIR, "supp/pdf/SUPP_F03_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F03_composite.png"), file.path(BOX_DIR, "supp/png/SUPP_F03_composite.png"), overwrite = TRUE)
# Supplementary data
file.copy("04_Figures/F03/c_data/F03_supplementary.xlsx", file.path(WRITING_DIR, "F03_supplementary.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F03/c_data/F03_supplementary.xlsx", file.path(BOX_DIR, "F03_source_data.xlsx"), overwrite = TRUE)

cat("  Synced to WRITING_DIR + BOX_DIR\n")

# Final cleanup: supp stitcher re-sources panels, which re-writes CSVs
remaining <- list.files("04_Figures/F03/c_data", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}
supp_dir <- "04_Figures/F03/c_data/supp"
if (dir.exists(supp_dir) && length(list.files(supp_dir)) == 0)
  unlink(supp_dir, recursive = TRUE)

cat("=== F03 complete ===\n")
