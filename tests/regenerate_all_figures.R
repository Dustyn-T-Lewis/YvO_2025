#!/usr/bin/env Rscript
# Regenerate all figures from upstream data.
# Run from project root: Rscript tests/regenerate_all_figures.R
# Useful for: verifying reproducibility, checking after pipeline changes.

setwd(rprojroot::find_rstudio_root_file())

figures <- list(
  "04_Figures/F01/a_script/YvO_F01_composite.R",
  "04_Figures/F02/a_script/YvO_F02_composite.R",
  "04_Figures/F03/a_script/YvO_F03_composite.R",
  "04_Figures/F04/a_script/YvO_F04_composite.R",
  "04_Figures/F05/a_script/YvO_F05_composite.R",
  "04_Figures/F06/a_script/YvO_F06_composite.R",
  "04_Figures/WGCNA_F07/a_script/YvO_WGCNA_F07_composite.R",
  "04_Figures/WGCNA_F08/a_script/YvO_WGCNA_F08_composite.R",
  "04_Figures/PLIER_F09/a_script/90_stitch_figure.R",
  "04_Figures/F10/a_script/90_stitch_figure.R"
)

cat("=== Regenerating all figures ===\n")
for (fig in figures) {
  if (!file.exists(fig)) {
    cat("  SKIP:", fig, "(not found)\n")
    next
  }
  cat("  Running:", fig, "...\n")
  result <- tryCatch(
    source(fig, local = new.env(parent = globalenv())),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    cat("  FAIL:", fig, "\n    ", conditionMessage(result), "\n")
  } else {
    cat("  OK:", fig, "\n")
  }
}
cat("=== Done ===\n")
