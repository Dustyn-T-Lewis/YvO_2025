# F06 — Master Orchestrator
# Sources WGCNA runner is NOT included (run separately).
# Sources main panels + supp panels, then stitchers.
#
# Run order:
#   1. main/90_stitch_main.R                    -> 3-panel main composite + xlsx
#   2. supp/01_QC/90_stitch_QC.R                -> QC composite
#   3. supp/02_analysis/90_stitch_analysis.R     -> analysis composite
#   4. supp/03_module/90_stitch_module_hub.R     -> hub composite
#   5. supp/03_module/90_stitch_module_triptych.R -> triptych composite
#   6. supp/90_stitch_mega.R                     -> mega SUPP composite

setwd(rprojroot::find_rstudio_root_file())

cat("=== F06: Running main composite + xlsx ===\n")
source("04_Figures/F06/a_script/main/90_stitch_main.R")

cat("=== F06: Running supp QC composite ===\n")
source("04_Figures/F06/a_script/supp/01_QC/90_stitch_QC.R")

cat("=== F06: Running supp analysis composite ===\n")
source("04_Figures/F06/a_script/supp/02_analysis/90_stitch_analysis.R")

cat("=== F06: Running supp module hub composite ===\n")
source("04_Figures/F06/a_script/supp/03_module/90_stitch_module_hub.R")

cat("=== F06: Running supp module triptych composite ===\n")
source("04_Figures/F06/a_script/supp/03_module/90_stitch_module_triptych.R")

cat("=== F06: Running supp mega composite ===\n")
source("04_Figures/F06/a_script/supp/90_stitch_mega.R")

# --- Sync to manuscript directories (skipped silently if dirs not present) ---
WRITING_DIR <- Sys.getenv("YVO_WRITING_DIR",
                          unset = "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F06")
BOX_DIR     <- Sys.getenv("YVO_BOX_DIR",
                          unset = "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F06_WGCNA")

SYNC_DIRS <- c(WRITING_DIR, BOX_DIR)
SYNC_DIRS <- SYNC_DIRS[dir.exists(dirname(SYNC_DIRS))]  # keep only dirs whose parent exists

if (length(SYNC_DIRS) > 0) {
  RPT <- "04_Figures/F06/b_reports"
  for (sync_dir in SYNC_DIRS) {
    for (d in c(file.path(sync_dir, "main/pdf"), file.path(sync_dir, "main/png"),
                file.path(sync_dir, "supp/pdf"), file.path(sync_dir, "supp/png")))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    # Main
    file.copy(file.path(RPT, "main/pdf/MAIN_F06_composite.pdf"),
              file.path(sync_dir, "main/pdf/MAIN_F06_composite.pdf"), overwrite = TRUE)
    file.copy(file.path(RPT, "main/png/MAIN_F06_composite.png"),
              file.path(sync_dir, "main/png/MAIN_F06_composite.png"), overwrite = TRUE)
    # Supp (multiple composites)
    for (f in list.files(file.path(RPT, "supp/pdf"), pattern = "^SUPP_F06_composite", full.names = TRUE))
      file.copy(f, file.path(sync_dir, "supp/pdf", basename(f)), overwrite = TRUE)
    for (f in list.files(file.path(RPT, "supp/png"), pattern = "^SUPP_F06_composite", full.names = TRUE))
      file.copy(f, file.path(sync_dir, "supp/png", basename(f)), overwrite = TRUE)
  }
  # Data (use original var names for distinct filenames)
  if (dir.exists(dirname(WRITING_DIR))) {
    file.copy("04_Figures/F06/c_data/F06_supplementary.xlsx",
              file.path(WRITING_DIR, "F06_supplementary.xlsx"), overwrite = TRUE)
  }
  if (dir.exists(dirname(BOX_DIR))) {
    file.copy("04_Figures/F06/c_data/F06_supplementary.xlsx",
              file.path(BOX_DIR, "F06_source_data.xlsx"), overwrite = TRUE)
  }
  cat("  Synced to available manuscript directories\n")
} else {
  cat("  Manuscript sync skipped (no sync directories found)\n")
}

cat("=== F06 complete ===\n")
