# F03 Panel D: Age x Training interaction volcano ring (Training_Old - Training_Young)
# Thin wrapper — defines spec and delegates to build_volcano_panel.R.
# Exposes pD for the composite stitcher.

spec <- list(contrast = "Interaction", title = "Age \u00d7 Training Interaction",
             subtitle = "Training_Old \u2212 Training_Young", tag = "D")

source("04_Figures/F03/a_script/main/panels/build_volcano_panel.R")
