# F03 Panel B: Training Young volcano ring (Young_Post - Young_Pre)
# Thin wrapper — defines spec and delegates to build_volcano_panel.R.
# Exposes pB for the composite stitcher.

spec <- list(contrast = "Training_Young", title = "Training Response (Young)",
             subtitle = "Young_Post \u2212 Young_Pre", tag = "B")

source("04_Figures/F03/a_script/main/panels/build_volcano_panel.R")
