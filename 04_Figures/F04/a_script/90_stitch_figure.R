# F04 — Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main stitcher (composite + xlsx + cleanup), then supp stitcher (composite).
#
# Run order:
#   1. supp/panels/*.R               -> supp panel CSVs + PNGs (pre-generate)
#   2. main/90_stitch_main.R         -> 5-panel main composite + xlsx + cleanup
#   3. supp/90_stitch_supp.R         -> supp composite from pre-rendered PNGs
#   4. Final cleanup                 -> remove any leftover CSVs

setwd(rprojroot::find_rstudio_root_file())

cat("=== F04: Pre-generating supp panel data + PNGs ===\n")
source("04_Figures/F04/a_script/supp/panels/enrichment_heatmap.R")
source("04_Figures/F04/a_script/supp/panels/fgsea_chord.R")
source("04_Figures/F04/a_script/supp/panels/ora_chord.R")
source("04_Figures/F04/a_script/supp/panels/quadrant_chord.R")
source("04_Figures/F04/a_script/supp/panels/pattern_ora.R")

cat("=== F04: Running main composite + xlsx ===\n")
source("04_Figures/F04/a_script/main/90_stitch_main.R")

cat("=== F04: Running supp composite ===\n")
source("04_Figures/F04/a_script/supp/90_stitch_supp.R")

# --- Sync to manuscript directories ---
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F04"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F04_training_concordance"
for (d in c(file.path(WRITING_DIR, "main/pdf"), file.path(WRITING_DIR, "main/png"),
            file.path(WRITING_DIR, "supp/pdf"), file.path(WRITING_DIR, "supp/png"),
            file.path(BOX_DIR, "main/pdf"), file.path(BOX_DIR, "main/png"),
            file.path(BOX_DIR, "supp/pdf"), file.path(BOX_DIR, "supp/png")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

RPT <- "04_Figures/F04/b_reports"
# Main
file.copy(file.path(RPT, "main/pdf/MAIN_F04_composite.pdf"), file.path(WRITING_DIR, "main/pdf/MAIN_F04_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F04_composite.png"), file.path(WRITING_DIR, "main/png/MAIN_F04_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "main/pdf/MAIN_F04_composite.pdf"), file.path(BOX_DIR, "main/pdf/MAIN_F04_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F04_composite.png"), file.path(BOX_DIR, "main/png/MAIN_F04_composite.png"), overwrite = TRUE)
# Supp
file.copy(file.path(RPT, "supp/pdf/SUPP_F04_composite.pdf"), file.path(WRITING_DIR, "supp/pdf/SUPP_F04_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F04_composite.png"), file.path(WRITING_DIR, "supp/png/SUPP_F04_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F04_composite.pdf"), file.path(BOX_DIR, "supp/pdf/SUPP_F04_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F04_composite.png"), file.path(BOX_DIR, "supp/png/SUPP_F04_composite.png"), overwrite = TRUE)
# Data
file.copy("04_Figures/F04/c_data/F04_supplementary.xlsx", file.path(WRITING_DIR, "F04_supplementary.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F04/c_data/F04_supplementary.xlsx", file.path(BOX_DIR, "F04_source_data.xlsx"), overwrite = TRUE)
cat("  Synced to WRITING_DIR + BOX_DIR\n")

# Final cleanup: remove any leftover CSVs
remaining <- list.files("04_Figures/F04/c_data", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}

cat("=== F04 complete ===\n")
