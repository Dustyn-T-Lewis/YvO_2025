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

# --- Sync to manuscript directories (skip silently if paths don't exist) ---
WRITING_DIR <- Sys.getenv("YVO_WRITING_DIR",
  unset = "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F07")
BOX_DIR <- Sys.getenv("YVO_BOX_DIR",
  unset = "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F07_phenotype_prediction")

SYNC_TARGETS <- c(WRITING_DIR, BOX_DIR)
SYNC_TARGETS <- SYNC_TARGETS[dir.exists(dirname(SYNC_TARGETS))]

if (length(SYNC_TARGETS) > 0) {
  for (tgt in SYNC_TARGETS) {
    for (d in c(file.path(tgt, "main/pdf"), file.path(tgt, "main/png"),
                file.path(tgt, "supp/pdf"), file.path(tgt, "supp/png")))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  RPT <- "04_Figures/F07/b_reports"
  for (tgt in SYNC_TARGETS) {
    # Main
    file.copy(file.path(RPT, "main/pdf/MAIN_F07_composite.pdf"), file.path(tgt, "main/pdf/MAIN_F07_composite.pdf"), overwrite = TRUE)
    file.copy(file.path(RPT, "main/png/MAIN_F07_composite.png"), file.path(tgt, "main/png/MAIN_F07_composite.png"), overwrite = TRUE)
    # Supp composites
    for (f in list.files(file.path(RPT, "supp/pdf"), pattern = "^SUPP_F07_composite", full.names = TRUE))
      file.copy(f, file.path(tgt, "supp/pdf", basename(f)), overwrite = TRUE)
    for (f in list.files(file.path(RPT, "supp/png"), pattern = "^SUPP_F07_composite", full.names = TRUE))
      file.copy(f, file.path(tgt, "supp/png", basename(f)), overwrite = TRUE)
    # LOSO panels (standalone supp outputs)
    for (f in c("SUPP_F07_loso_sensitivity", "SUPP_F07_loso_wgcna_refit")) {
      file.copy(file.path(RPT, "supp/pdf/panels", paste0(f, ".pdf")), file.path(tgt, "supp/pdf", paste0(f, ".pdf")), overwrite = TRUE)
      file.copy(file.path(RPT, "supp/png/panels", paste0(f, ".png")), file.path(tgt, "supp/png", paste0(f, ".png")), overwrite = TRUE)
    }
    # Data
    file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx", file.path(tgt, "F07_supplementary.xlsx"), overwrite = TRUE)
  }
  # Box gets a different xlsx filename
  if (dir.exists(dirname(BOX_DIR)) && BOX_DIR %in% SYNC_TARGETS) {
    file.copy("04_Figures/F07/c_data/F07_supplementary.xlsx", file.path(BOX_DIR, "F07_source_data.xlsx"), overwrite = TRUE)
  }
  cat("  Synced to", length(SYNC_TARGETS), "manuscript dir(s)\n")
} else {
  cat("  Sync skipped: no manuscript directories found on this machine\n")
}

cat("=== F07 complete ===\n")
