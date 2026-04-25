# F03 Panel A: Aging volcano ring (Old_Pre - Young_Pre)
# Thin wrapper — defines spec and delegates to build_volcano_panel.R.
# Exposes pA for the composite stitcher.

spec <- list(contrast = "Aging", title = "Aging Effect",
             subtitle = "Old_Pre \u2212 Young_Pre", tag = "A")

source("04_Figures/F03/a_script/main/panels/build_volcano_panel.R")
