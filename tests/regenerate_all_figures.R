#!/usr/bin/env Rscript
# Regenerate all figures from upstream data.
# Run from project root: Rscript tests/regenerate_all_figures.R
# Useful for: verifying reproducibility, checking after pipeline changes.

setwd(rprojroot::find_rstudio_root_file())

figures <- list(
  "04_Figures/F0/a_script/YvO_F0_composite.R",
  "04_Figures/F1/a_script/YvO_figure1_composite.R",
  "04_Figures/F2/a_script/YvO_figure2_composite.R",
  "04_Figures/F3/a_script/YvO_figure3_composite.R",
  "04_Figures/F4/a_script/YvO_figure4_composite.R",
  "04_Figures/F5/a_script/YvO_figure5_composite.R",
  "04_Figures/F6/a_script/YvO_figure6_composite.R"
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
