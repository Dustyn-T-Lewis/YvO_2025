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

# --- Sync to manuscript directories (conditional on directory existence) ---
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F02"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript/02_Figures/F02_proteome_DEP"
RPT <- "04_Figures/F02/b_reports"

sync_targets <- c(WRITING_DIR, BOX_DIR)
names(sync_targets) <- c("WRITING_DIR", "BOX_DIR")

for (tgt_name in names(sync_targets)) {
  tgt <- sync_targets[[tgt_name]]
  # Only sync if the parent cloud directory is mounted
  if (!dir.exists(dirname(tgt))) {
    cat(sprintf("  %s not mounted — skipping sync\n", tgt_name))
    next
  }
  for (d in c(file.path(tgt, "main/pdf"), file.path(tgt, "main/png"),
              file.path(tgt, "supp/pdf"), file.path(tgt, "supp/png")))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  file.copy(file.path(RPT, "main/pdf/MAIN_F02_composite.pdf"), file.path(tgt, "main/pdf/MAIN_F02_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F02_composite.png"), file.path(tgt, "main/png/MAIN_F02_composite.png"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F02_composite.pdf"), file.path(tgt, "supp/pdf/SUPP_F02_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F02_composite.png"), file.path(tgt, "supp/png/SUPP_F02_composite.png"), overwrite = TRUE)

  xlsx_name <- if (tgt_name == "BOX_DIR") "F02_source_data.xlsx" else "F02_supplementary.xlsx"
  file.copy("04_Figures/F02/c_data/F02_supplementary.xlsx", file.path(tgt, xlsx_name), overwrite = TRUE)

  cat(sprintf("  Synced to %s\n", tgt_name))
}

# Final cleanup: supp stitcher re-sources panels, which re-writes CSVs
remaining <- list.files("04_Figures/F02/c_data", pattern = "\\.csv$", full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  cat(sprintf("  final cleanup: removed %d leftover CSV(s)\n", length(remaining)))
}

cat("=== F02 complete ===\n")
