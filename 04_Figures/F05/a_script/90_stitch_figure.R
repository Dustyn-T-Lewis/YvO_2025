# F05 — Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main stitcher (composite + xlsx + cleanup), then supp stitcher (composite).
#
# Run order:
#   1. supp/panels/*.R               -> supp panel CSVs + PNGs (pre-generate)
#   2. main/90_stitch_main.R         -> 5-panel main composite + xlsx + cleanup
#   3. supp/90_stitch_supp.R         -> supp composite from pre-rendered PNGs
#   4. Final cleanup                 -> remove any leftover CSVs

setwd(rprojroot::find_rstudio_root_file())

cat("=== F05: Pre-generating supp panel data + PNGs ===\n")
source("04_Figures/F05/a_script/supp/panels/enrichment_heatmap.R")
# Chord / pattern-ORA scripts archived (not in manuscript supplementary):
# source("04_Figures/F05/a_script/supp/panels/fgsea_chord.R")
# source("04_Figures/F05/a_script/supp/panels/ora_chord.R")
# source("04_Figures/F05/a_script/supp/panels/quadrant_chord.R")
# source("04_Figures/F05/a_script/supp/panels/pattern_ora.R")

cat("=== F05: Running supp diagnostics composite ===\n")
source("04_Figures/F05/a_script/supp/90_stitch_supp.R")

cat("=== F05: Running main composite + xlsx ===\n")
source("04_Figures/F05/a_script/main/90_stitch_main.R")

# --- Copy to Box manuscript directory ---
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- "04_Figures/F05/b_reports"
file.copy(file.path(RPT, "main/pdf/MAIN_F05_composite.pdf"),
          file.path(BOX, "02_Figures/F05_aging_reversal/main/pdf/MAIN_F05_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F05_composite.png"),
          file.path(BOX, "02_Figures/F05_aging_reversal/main/png/MAIN_F05_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F05_diagnostics.pdf"),
          file.path(BOX, "02_Figures/F05_aging_reversal/supp/pdf/SUPP_F05_diagnostics.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F05_diagnostics.png"),
          file.path(BOX, "02_Figures/F05_aging_reversal/supp/png/SUPP_F05_diagnostics.png"), overwrite = TRUE)
file.copy("04_Figures/F05/c_data/F05_supplementary.xlsx",
          file.path(BOX, "02_Figures/F05_aging_reversal/F05_source_data.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F05/c_data/F05_supplementary.xlsx",
          file.path(BOX, "03_Supplementary_Tables/S11_F05_aging_reversal.xlsx"), overwrite = TRUE)

# Final cleanup: remove any leftover CSVs
remaining <- list.files("04_Figures/F05/c_data", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}

cat("=== F05 complete ===\n")
