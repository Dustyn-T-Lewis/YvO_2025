# F03 Panel C: Training Old volcano ring (Old_Post - Old_Pre)
# Thin wrapper — defines spec and delegates to build_volcano_panel.R.
# Exposes pC for the composite stitcher.

spec <- list(contrast = "Training_Old", title = "Training Response (Old)",
             subtitle = "Old_Post \u2212 Old_Pre", tag = "C")

source("04_Figures/F03/a_script/main/panels/build_volcano_panel.R")
