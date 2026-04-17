# F02 — Style overrides
# Sources canonical palettes/theme/helpers from shared/style.R then adds
# only the constants that differ for F02 (wider panels than default).
source("04_Figures/shared/style.R")

# Diverging heatmap anchors (Panel B ΔCV scatter)
HEATMAP_LO <- "#2166AC"
HEATMAP_HI <- "#B2182B"

# F02-specific annotation sizing (panels are wider than the 180mm reference)
BASE_COUNT <- 4.0
BASE_GENE  <- 3.8
BASE_STAT  <- 4.0
