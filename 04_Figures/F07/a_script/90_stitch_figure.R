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

# --- Sync to manuscript directories ---
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F07"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F07_phenotype_prediction"
for (d in c(file.path(WRITING_DIR, "main/pdf"), file.path(WRITING_DIR, "main/png"),
            file.path(WRITING_DIR, "supp/pdf"), file.path(WRITING_DIR, "supp/png"),
            file.path(BOX_DIR, "main/pdf"), file.path(BOX_DIR, "main/png"),
            file.path(BOX_DIR, "supp/pdf"), file.path(BOX_DIR, "supp/png")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

RPT <- "04_Figures/F07/b_reports"
# Main
file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"), file.path(WRITING_DIR, "main/pdf/MAIN_F07_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"), file.path(WRITING_DIR, "main/png/MAIN_F07_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"), file.path(BOX_DIR, "main/pdf/MAIN_F07_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"), file.path(BOX_DIR, "main/png/MAIN_F07_composite.png"), overwrite = TRUE)
# Supp composites
for (f in list.files(file.path(RPT, "supp/pdf"), pattern = "^SUPP_F07_composite", full.names = TRUE))
  file.copy(f, file.path(WRITING_DIR, "supp/pdf", basename(f)), overwrite = TRUE)
for (f in list.files(file.path(RPT, "supp/png"), pattern = "^SUPP_F07_composite", full.names = TRUE))
  file.copy(f, file.path(WRITING_DIR, "supp/png", basename(f)), overwrite = TRUE)
for (f in list.files(file.path(RPT, "supp/pdf"), pattern = "^SUPP_F07_composite", full.names = TRUE))
  file.copy(f, file.path(BOX_DIR, "supp/pdf", basename(f)), overwrite = TRUE)
for (f in list.files(file.path(RPT, "supp/png"), pattern = "^SUPP_F07_composite", full.names = TRUE))
  file.copy(f, file.path(BOX_DIR, "supp/png", basename(f)), overwrite = TRUE)
# LOSO panels (standalone supp outputs)
for (f in c("SUPP_F07_loso_sensitivity", "SUPP_F07_loso_wgcna_refit")) {
  file.copy(file.path(RPT, "supp/pdf/panels", paste0(f, ".pdf")), file.path(WRITING_DIR, "supp/pdf", paste0(f, ".pdf")), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/panels", paste0(f, ".png")), file.path(WRITING_DIR, "supp/png", paste0(f, ".png")), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/panels", paste0(f, ".pdf")), file.path(BOX_DIR, "supp/pdf", paste0(f, ".pdf")), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/panels", paste0(f, ".png")), file.path(BOX_DIR, "supp/png", paste0(f, ".png")), overwrite = TRUE)
}
# Data
file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx", file.path(WRITING_DIR, "F07_supplementary.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx", file.path(BOX_DIR, "F07_source_data.xlsx"), overwrite = TRUE)
cat("  Synced to WRITING_DIR + BOX_DIR\n")

cat("=== F07 complete ===\n")
